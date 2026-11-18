r"""Build one catalogue and run every row of its load past the schema — with no database, and throw it away.

⛔ **The pangenome stage of a large cohort is a day on a whole node, and it is one transaction.** A
refusal in its last loader, or one value Postgres will not take, discards all of it. Both of the ways it
can fail are properties of the cohort's DATA meeting this code, and neither needs a database to find:

* a **refusal** the frame build or a loader raises — a similarity artifact of another shape, a `.loci`
  that disagrees, a parser that drops a value the frame counted;
* a **row the schema will not hold** — a string longer than its column, an integer past its type, a NULL
  where the column forbids one, text where bytes belong.

So this builds the frames exactly as the ingest does, derives the locus columns, and then drives the
ingest's OWN loaders with ids invented on the spot and a `copy_rows` that sends nothing: each row is
coerced as the real loader coerces it and checked against the SQLAlchemy column it was headed for. The
`.loci` oracle runs for real, against the files.

⚠ **What it cannot see, said here so a clean run is not read as more than it is:** constraints the database
evaluates across rows (uniqueness, foreign keys, the CHECKs that are not a length or a NaN), the genome
layer already in the database (gene counts are taken from the `.loci` and `_meta` files instead), and the
intergenic-gap FEATURE rows, which need the gap ids the database assigns.

    python -m bacatlas_backend.instruments.catalogue_dry_run --data-root … --set-key … --species-key … \\
        --model-label … --run-id … [--gene-universe …] [--assignment …] [--nseq-dir …] [--go-fold fold.json]

Exit 0 when nothing was refused and every row fits; 1 otherwise, naming each.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from sqlalchemy import ARRAY, BigInteger, Boolean, Enum, Float, Integer, LargeBinary, Numeric, SmallInteger, String

#: The integer ranges Postgres enforces, by SQLAlchemy type. ⚠ Checked most-specific first: `SmallInteger`
#: and `BigInteger` are both `Integer` subclasses.
INTEGER_RANGES = (
    (SmallInteger, 2**15),
    (BigInteger, 2**63),
    (Integer, 2**31),
)

#: How many offending values are quoted per (table, column, kind).
EXAMPLES = 3

#: ⛔ A btree index entry holds at most 2,704 bytes (Postgres 16, 8 kB pages). A `Text` column has no length
#: of its own, so the ONLY limit on it is this one, where it sits in a btree — and COPY refuses the row with
#: *"index row size … exceeds btree version 4 maximum 2704"*. The margin below it is for the entry's other
#: columns and header; a value over this is named as a risk, because a compressible one may still fit.
BTREE_ENTRY_BYTES = 2704
BTREE_TEXT_BYTES = 2600


class NeedsDatabase(RuntimeError):
    """A loader asked the session a question only a database can answer — that part is not dry-runnable."""


class DrySession:
    """What a loader may do to a session in a dry run: add an object, flush. Anything else is named."""

    def __init__(self) -> None:
        self.added: Counter = Counter()

    def add(self, row) -> None:
        """Count the object by its model; nothing is stored."""
        self.added[type(row).__name__] += 1

    def flush(self) -> None:
        """Nothing to flush."""
        return None

    def execute(self, *_args, **_kwargs):
        """A question for a database: refused by name, so the caller reports that part as not run."""
        raise NeedsDatabase("the loader queried the database")

    query = execute


class SchemaFit:
    """Every row a loader would COPY, checked against the column it is headed for; nothing is sent."""

    def __init__(self) -> None:
        self.rows: Counter = Counter()
        self.violations: Counter = Counter()
        self.examples: dict = {}
        self._btree: dict = {}

    def _note(self, table: str, column: str, kind: str, value) -> None:
        key = (table, column, kind)
        self.violations[key] += 1
        shown = self.examples.setdefault(key, [])
        if len(shown) < EXAMPLES:
            shown.append(repr(value)[:160])

    def _btree_limited(self, table, column) -> bool:
        """Is this column part of a btree — an index, the primary key or a unique constraint?"""
        key = (table.name, column.name)
        if key not in self._btree:
            in_index = any(
                column.name in index.columns and (index.dialect_options["postgresql"]["using"] or "btree") == "btree"
                for index in table.indexes
            )
            in_constraint = any(
                column.name in constraint.columns
                for constraint in table.constraints
                if constraint.__class__.__name__ in ("PrimaryKeyConstraint", "UniqueConstraint")
            )
            self._btree[key] = in_index or in_constraint
        return self._btree[key]

    def _check(self, table, column, value) -> None:
        kind = column.type
        if value is None:
            if not column.nullable and column.default is None and column.server_default is None:
                self._note(table.name, column.name, "NULL in a NOT NULL column", value)
            return
        if isinstance(kind, ARRAY):
            if not isinstance(value, list):
                self._note(table.name, column.name, "not a list in an array column", value)
            elif any(item is None for item in value):
                self._note(table.name, column.name, "NULL element in an array", value[:5])
            elif isinstance(kind.item_type, Integer):
                bound = next(limit for kind_, limit in INTEGER_RANGES if isinstance(kind.item_type, kind_))
                if any(isinstance(item, bool) or not isinstance(item, int) for item in value):
                    self._note(table.name, column.name, "array element that is not an integer", value[:5])
                elif any(not -bound <= item < bound for item in value):
                    self._note(table.name, column.name, "array element out of range", value[:5])
            elif isinstance(kind.item_type, String) and any(not isinstance(item, str) for item in value):
                self._note(table.name, column.name, "array element that is not text", value[:5])
            return
        if isinstance(kind, LargeBinary):
            if not isinstance(value, (bytes, bytearray, memoryview)):
                self._note(table.name, column.name, "not bytes in a bytea column", value)
            return
        if isinstance(kind, Enum):
            names = set(kind.enums)
            if str(value) not in names:
                self._note(table.name, column.name, f"not one of the enum's {len(names)} names", value)
            return
        if isinstance(kind, String):
            if not isinstance(value, str):
                self._note(table.name, column.name, "not text in a text column", value)
            elif kind.length is not None and len(value) > kind.length:
                shown = f"{len(value)} chars: {value}"
                self._note(table.name, column.name, f"longer than String({kind.length})", shown)
            elif "\x00" in value:
                self._note(table.name, column.name, "a NUL byte, which text cannot hold", value)
            elif len(value) > BTREE_TEXT_BYTES // 4 and self._btree_limited(table, column):
                size = len(value.encode("utf-8", "surrogatepass"))
                if size > BTREE_TEXT_BYTES:
                    shown = f"{size:,} bytes: {value[:80]}…"
                    self._note(
                        table.name, column.name,
                        f"may not fit a btree index entry ({BTREE_ENTRY_BYTES:,} bytes)", shown,
                    )  # fmt: skip
            return
        if isinstance(kind, Boolean):
            if not isinstance(value, bool):
                self._note(table.name, column.name, "not a boolean in a boolean column", value)
            return
        if isinstance(kind, Integer):
            bound = next(limit for kind_, limit in INTEGER_RANGES if isinstance(kind, kind_))
            if isinstance(value, bool) or not isinstance(value, int):
                self._note(table.name, column.name, "not an integer in an integer column", value)
            elif not -bound <= value < bound:
                self._note(table.name, column.name, f"outside ±{bound:,}", value)
            return
        if isinstance(kind, (Float, Numeric)):
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                self._note(table.name, column.name, "not a number in a numeric column", value)
                return
        if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
            self._note(table.name, column.name, "NaN or infinity", value)

    def copy_rows(self, _session, table, columns, rows, *, batch_rows: int = 0) -> int:
        """`staging_table_loader.copy_rows`' contract — same refusals, same coercion — sending nothing."""
        from bacatlas_backend.ingest.staging_table_loader import _refuse_text_in_a_binary_column, coerce_value

        unknown = [name for name in columns if name not in table.c]
        if unknown:
            raise KeyError(f"{table.name} has no column(s) {unknown!r}")
        targets = [table.c[name] for name in columns]
        binary = tuple(i for i, column in enumerate(targets) if isinstance(column.type, LargeBinary))
        written = 0
        for row in rows:
            if binary and not written:
                _refuse_text_in_a_binary_column(table, columns, binary, row)
            if len(row) != len(targets):
                self._note(table.name, "(row)", f"{len(row)} values for {len(targets)} columns", row[:4])
            for column, value in zip(targets, row, strict=False):
                self._check(table, column, coerce_value(value))
            written += 1
        self.rows[table.name] += written
        return written


def _fold(path: Path | None):
    """The GO fold the ingest reads from its reference tables, from a JSON export of them — or `None`."""
    if path is None:
        return None
    from bacatlas_backend.instruments.gene_ontology_fold import GeneOntologyFold

    stored = json.loads(path.read_text())
    fold = GeneOntologyFold(
        slim_of={key: frozenset(value) for key, value in stored["slim_of"].items()},
        ancestors={key: frozenset(value) for key, value in stored["ancestors"].items()},
        namespace_of={key: int(value) for key, value in stored["namespace_of"].items()},
        name_of=stored.get("name_of", {}),
    )
    if stored.get("version") and stored["version"] != fold.version():
        raise SystemExit(f"{path} says it is fold {stored['version']} and rebuilds as {fold.version()}")
    return fold


def _gene_counts(artifacts, samples: list[str], genome_id_by_sample: dict, workers: int):
    """Each genome's gene count without a database: its `.nseq` gene table's, or its `_meta` row count.

    ⛔ **The `.nseq`, never the `.loci`.** The genome stage stores `genome_genes.gene_count` from the
    `.nseq` header, and the oracle then refuses a `.loci` whose length is not `gene_count × 4`. Taking
    the count from the `.loci` instead made that length check compare the file with itself — a `.loci`
    of the wrong length passed the dry run and would have refused the load.
    """
    import pyarrow.parquet as parquet

    from bacatlas_backend.gene_blocks.nseq_gene_table import read_nseq_gene_table
    from bacatlas_backend.ingest.ingest_locus_catalogue import GeneCounts

    present = artifacts.nseq_present() if artifacts.nseq_root is not None else frozenset()

    def count(sample: str) -> tuple[str, int, bool]:
        if sample in present:
            table = read_nseq_gene_table(artifacts.per_genome_nseq(sample), expect_sample=sample)
            return sample, table.gene_count, True
        return sample, parquet.read_metadata(artifacts.per_genome(sample, "meta")).num_rows, False

    with ThreadPoolExecutor(max(1, workers)) as pool:
        counted = list(pool.map(count, samples))
    by_genome = {genome_id_by_sample[sample]: genes for sample, genes, _ in counted}
    from_meta = [sample for sample, _, from_nseq in counted if not from_nseq]
    return GeneCounts(by_genome=by_genome, omitted_reasons={}), from_meta


def main(argv: list[str] | None = None) -> int:
    """Run the dry run; 0 clean, 1 a refusal or a row that does not fit."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--set-key", required=True)
    parser.add_argument("--species-key", default=None)
    parser.add_argument("--model-label", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--gene-universe", type=Path, default=None)
    parser.add_argument("--assignment", type=Path, default=None)
    parser.add_argument("--nseq-dir", type=Path, default=None)
    parser.add_argument("--go-fold", type=Path, default=None, help="a JSON export of the reference GO fold")
    parser.add_argument("--allow-missing-strand", action="store_true")
    parser.add_argument("--read-workers", type=int, default=16)
    args = parser.parse_args(argv)

    from bacatlas_backend.ingest import ingest_locus_catalogue as loader
    from bacatlas_backend.ingest.artifact_locator import CatalogueArtifacts
    from bacatlas_backend.ingest.catalogue_frames import build_catalogue_frames, load_pfam_name_table
    from bacatlas_backend.ingest.phase_report import phase
    from bacatlas_backend.models.locus import Locus

    artifacts = CatalogueArtifacts(
        data_root=args.data_root, set_key=args.set_key, model_label=args.model_label, run_id=args.run_id,
        species_key=args.species_key, gene_universe_path=args.gene_universe, assignment_path=args.assignment,
        nseq_root=args.nseq_dir,
    )  # fmt: skip
    absent = [role for role, path, _ in artifacts.required() if not path.exists()]
    print(f"required artifacts absent: {absent or 'none'}  (not refused here; the build will, if it reads one)")
    fold = _fold(args.go_fold)
    fold_said = fold.version() if fold is not None else "none given — GO-folded verdicts take the no-reference path"
    print(f"annotation fold: {fold_said}")

    frames = build_catalogue_frames(
        artifacts, annotation_fold=fold, allow_missing_strand=args.allow_missing_strand
    )
    print(
        f"frames: {len(frames.samples):,} genomes · {len(frames.loci):,} loci · "
        f"{len(frames.gene_memberships):,} gene memberships · {len(frames.arrangements):,} arrangements · "
        f"{len(frames.offset_occupants):,} offset occupants · {0 if frames.gaps is None else len(frames.gaps):,} gaps"
    )
    for section, reason in sorted(frames.omitted.items()):
        print(f"  omitted {section}: {str(reason)[:200]}")

    refusals: list[str] = []
    oracle: tuple[int, int] | None = None
    fit = SchemaFit()
    session = DrySession()
    derived = loader._derive_locus_columns(
        frames, load_pfam_name_table(), loader._policy(), loader._reference_names()["cog_orthogroup"]
    )
    phase("dry run: derived locus columns")

    genome_id_by_sample = {sample: index + 1 for index, sample in enumerate(frames.samples)}
    locus_id_by_label = {str(label): index + 1 for index, label in enumerate(frames.loci["node"])}
    locus_id_by_ordinal = [None] * len(frames.loci)
    for label, ordinal in zip(frames.loci["node"], frames.loci["catalogue_ordinal"], strict=True):
        locus_id_by_ordinal[int(ordinal)] = locus_id_by_label[str(label)]
    genome_id_by_ordinal = [genome_id_by_sample[sample] for sample in frames.samples]
    gene_counts, counted_from_meta = _gene_counts(
        artifacts, frames.samples, genome_id_by_sample, args.read_workers
    )
    print(
        f"gene counts: {len(frames.samples) - len(counted_from_meta):,} genomes from their .nseq gene table, "
        f"{len(counted_from_meta):,} from their _meta row count {counted_from_meta[:5]}"
    )
    phase("dry run: gene counts")

    real_copy = loader.copy_rows
    loader.copy_rows = fit.copy_rows
    steps = (
        ("locus rows", lambda: fit.copy_rows(
            session, Locus.__table__, loader.LOCUS_COLUMNS,
            loader._locus_rows(frames, derived, 1, 1, genome_id_by_sample))),
        ("annotation entries", lambda: loader._load_annotation_entries(session, frames, locus_id_by_label)),
        ("uniref crosstab", lambda: loader._load_uniref_crosstab(session, frames, locus_id_by_label)),
        ("arrangements", lambda: loader._load_arrangements(
            session, frames, 1, locus_id_by_label, genome_id_by_ordinal)),
        ("offset occupants", lambda: loader._load_offset_occupants(
            session, frames, 1, locus_id_by_label, locus_id_by_ordinal)),
        ("genome locus assignments + the .loci oracle", lambda: loader._load_gene_memberships(
            session, frames, 1, genome_id_by_sample,
            loci_root=None if args.nseq_dir is None else Path(args.nseq_dir).expanduser(),
            gene_counts=gene_counts, read_workers=args.read_workers)),
        ("intergenic gaps", lambda: loader._load_intergenic_gaps(session, frames, 1, locus_id_by_label)),
        ("similarity", lambda: loader._load_cluster_similarity(session, frames, 1, locus_id_by_label, derived)),
    )  # fmt: skip
    try:
        for name, step in steps:
            try:
                result = step()
                print(f"  ran  {name}: {result}")
                if hasattr(result, "loci_checked"):
                    oracle = (result.loci_checked, result.loci_absent)
            except NeedsDatabase as needs:
                print(f"  PART {name}: {needs} — the rows before that question were checked, the rest were not")
            except Exception as error:  # noqa: BLE001 — each refusal is reported and the run goes on to the next
                refusals.append(f"{name}: {type(error).__name__}: {str(error)[:600]}")
                print(f"  REFUSED {name}: {type(error).__name__}: {str(error)[:600]}")
            phase(f"dry run: {name[:22]}")
    finally:
        loader.copy_rows = real_copy

    print("rows a load would write (not counting what a refusal cut short):")
    for table, count in sorted(fit.rows.items()):
        print(f"  {table:34s} {count:>13,}")
    for model, count in sorted(session.added.items()):
        print(f"  {model:34s} {count:>13,}  (session.add)")
    if fit.violations:
        total = sum(fit.violations.values())
        print(f"ROWS THE SCHEMA WOULD NOT HOLD — {total:,} values in {len(fit.violations)} column(s):")
        for (table, column, kind), count in sorted(fit.violations.items()):
            shown = "; ".join(fit.examples[(table, column, kind)])
            print(f"  {table}.{column}: {kind} — {count:,} value(s), e.g. {shown}")
    verdict = "CLEAN" if not refusals and not fit.violations else "NOT CLEAN"
    # ⚠ Said beside the verdict and not folded into it: CLEAN with no `.loci` compared is a clean run of a
    # weaker check, and a reader must be able to tell the two apart from this line alone.
    compared = (
        "the .loci oracle did not run"
        if oracle is None
        else f"the .loci oracle compared {oracle[0]:,} genomes byte for byte; {oracle[1]:,} have no .loci"
    )
    print(
        f"DRY RUN {verdict}: {len(refusals)} refusal(s), "
        f"{len(fit.violations)} column(s) with values the schema would not hold — {compared}"
    )
    print(f"peak memory of this process: {_peak_gib():.2f} GiB (pid {os.getpid()})")
    return 0 if verdict == "CLEAN" else 1


def _peak_gib() -> float:
    from bacatlas_backend.ingest.phase_report import _peak_bytes

    return _peak_bytes() / 2**30


if __name__ == "__main__":
    raise SystemExit(main())
