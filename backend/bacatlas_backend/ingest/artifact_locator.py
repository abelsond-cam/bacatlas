"""Where every input of one catalogue lives — resolved once, checked all at once.

⛔ **This is derived from the naming conventions, NOT from `meta.audit.sources`.** That key
under-reports: the published *E. coli* payload names two artifacts (`_pfam_concordance.tsv` and
`_audit_summary.json`) while the export read **three** — the homology waterfall CSV supplied the
collapse tier, the resolved threshold and the ESM distance pair, and is never listed. A locator
built from it would silently omit the file three columns come from.

⛔ **Every missing artifact is reported together, and named.** One `FileNotFoundError` per run means
finding a second gap costs a second run; at eight artifact families and two representations that is
how a load takes an afternoon. `MissingArtifacts` lists all of them, with what each would have
supplied, so one look is enough.

⚠ **Three filename traps, all of them measured on the real store:**

1. **The output filename doubles the set token.** `locus_browser_{set}_{label}.json`, and `{label}`
   itself begins with the set — `locus_browser_ecoli_ecoli_nuna4_g2_0.98_…json`. A locator that
   assumes one `ecoli` finds nothing.
2. **A run's siblings are DERIVED, not passed.** `export_payload._nearest_sibling` finds the nearest
   -loci CSV by string-replacing `_cluster_similarity_` in the similarity CSV's own path, so it is
   addressed relative to it rather than rebuilt from the label. Reproduced here for the same reason:
   one run's medians can then never be paired with another run's neighbour list, which is a
   mismatch nothing downstream could detect — both would be in range and plausible. (The retired
   neighbourhood map addressed its two siblings the same way.)
3. **The payload JSON is not under `analysis/`.** The exporter writes it beside the maps on CSD3,
   but the local pull lands it in `data/browser/`. It is an ORACLE here, never an input.
"""

from __future__ import annotations

import os
from collections import deque
from collections.abc import Callable, Iterable, Iterator
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from itertools import islice
from pathlib import Path
from typing import TypeVar

#: The two representations, in the payload's tab order — Bacformer first, because it is the axis
#: the method is about.
REPRESENTATIONS = ("bacformer", "esm")

#: ⭐ How many reader threads each per-genome I/O pass runs, and how far ahead of the serial consumer
#: they may get (per thread). On Lustre a cold `stat` or `open` is ~28–48 ms of pure latency and the
#: GIL is released for all of it, so a few dozen threads turn ~4 h of serial waiting into minutes —
#: while the window keeps what is held in memory to a few dozen genomes' inputs, never the cohort's.
DEFAULT_READ_WORKERS = 16
READ_AHEAD_PER_WORKER = 4

_Item = TypeVar("_Item")
_Read = TypeVar("_Read")


def listed_files(directory: Path, suffix: str) -> frozenset[str]:
    """The stems of every regular file in ``directory`` whose name ends in ``suffix`` — ONE scan.

    ⛔ **The replacement for a per-genome `is_file()`, and it answers the same question.** A cold
    `stat` on Lustre is ~28–48 ms, and the 71,435-genome cohort asked it once per genome per tree:
    33–57 min each time, for an answer one `readdir` gives in about a second. ⚠ An entry counts only
    where `DirEntry.is_file()` says so — which follows a symlink exactly as `Path.is_file()` does — so
    a directory or a dangling link of that name is absent here, as it was there. An absent directory
    is an empty set, as `is_file()` was `False` for every path under it.
    """
    try:
        with os.scandir(directory) as entries:
            return frozenset(
                entry.name[: -len(suffix)]
                for entry in entries
                if entry.name.endswith(suffix) and entry.is_file()
            )
    except (FileNotFoundError, NotADirectoryError):
        return frozenset()


def read_ahead_in_order(
    items: Iterable[_Item], read: Callable[[_Item], _Read], workers: int
) -> Iterator[tuple[_Item, _Read]]:
    """`(item, read(item))` for every item, IN THE ORDER GIVEN, with the reads run ahead on a thread pool.

    ⭐ The consumer stays serial — it is the SQL loop, or a COPY generator — and sees exactly the
    sequence it would have computed itself; only the waiting moves off its thread. ⛔ The order is
    the input's, never completion order: a genome's rows and every id the database assigns follow it.
    ``read`` must not raise for a per-item failure — it captures it, so that failure surfaces when the
    consumer reaches THAT item, as it did when the consumer read the file itself. ``workers <= 1``
    runs every read inline, in order, on the caller's thread: the same values with no pool.
    """
    if workers <= 1:
        for item in items:
            yield item, read(item)
        return
    iterator = iter(items)
    with ThreadPoolExecutor(max_workers=workers, thread_name_prefix="bacatlas-read") as pool:
        pending = deque((item, pool.submit(read, item)) for item in islice(iterator, workers * READ_AHEAD_PER_WORKER))
        while pending:
            item, future = pending.popleft()
            for following in islice(iterator, 1):
                pending.append((following, pool.submit(read, following)))
            yield item, future.result()


class MissingArtifacts(FileNotFoundError):
    """Raised with EVERY absent required artifact, not the first one."""

    def __init__(self, missing: list[tuple[str, Path, str]]) -> None:
        self.missing = missing
        listing = "\n  ".join(
            f"{role:28s} {path}\n{'':30s} -> supplies {supplies}" for role, path, supplies in missing
        )
        super().__init__(f"{len(missing)} required artifact(s) absent:\n  {listing}")


@dataclass(frozen=True)
class CatalogueArtifacts:
    """Every input for one (species, model) catalogue, addressed from one data root.

    ``set_key`` is the dataset token in filenames (`ecoli`, `kp`) and ``species_key`` is the browser
    key. They are the same for the two published species and are kept separate anyway, because
    ``meta.species`` in the parquets is a THIRD vocabulary (`ecoli` | `kpneumoniae`) and collapsing
    any two of the three is how a load silently drops a whole species.
    """

    data_root: Path
    set_key: str
    model_label: str
    run_id: str
    species_key: str | None = None
    #: ⭐ The two inputs whose DERIVED names are a convention a big cohort does not keep, named
    #: outright instead. ⛔ They are overrides and not a second naming rule: a cohort that follows
    #: the convention passes neither, and the derived paths below are unchanged for it.
    #:
    #: ⛔ **And they are why a symlink is the wrong fix.** `__post_init__` refuses unless
    #: `model_label` starts with `{set_key}_`, so choosing a `set_key` to make a symlink resolve
    #: also pins the model label, the Pfam path, the payload path and `species_key`. The set key is
    #: load-bearing in six places; a filename is one.
    gene_universe_path: Path | None = None
    assignment_path: Path | None = None
    #: ⭐ Where nuna published `{BS}.nseq` / `{BS}.loci` — the gene table the ingest COPIES rather than
    #: re-encodes. ⛔ Not derived from `data_root`: the probe's 200 live in the bacatlas repo under
    #: `data/seq/`, the 71,435-genome cohort's under the model's own tree on CSD3, and neither is under
    #: `proc/`. `None` means no tree was named, and then nothing is read from one — the probe's loads
    #: encode the block from the GFF exactly as before.
    nseq_root: Path | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "data_root", Path(self.data_root).expanduser())
        for field_name in ("gene_universe_path", "assignment_path", "nseq_root"):
            value = getattr(self, field_name)
            if value is not None:
                object.__setattr__(self, field_name, Path(value).expanduser())
        if self.species_key is None:
            object.__setattr__(self, "species_key", self.set_key)
        if not self.model_label.startswith(f"{self.set_key}_"):
            raise ValueError(
                f"model_label {self.model_label!r} does not begin with the set token "
                f"{self.set_key + '_'!r}. Every artifact below is addressed by one or the other, "
                "and a mismatch resolves to paths that exist for a DIFFERENT catalogue."
            )

    # ── roots ─────────────────────────────────────────────────────────────────────────────────
    @property
    def processed_root(self) -> Path:
        """Everything the pipeline produced, as opposed to what it was given."""
        return self.data_root / "proc"

    @property
    def analysis_root(self) -> Path:
        """The audit, the assignments, the run manifests and the browser artifacts."""
        return self.processed_root / "analysis"

    @property
    def per_genome_root(self) -> Path:
        """The `{BS}_{meta,product,strand,dbxref,noncoding}.parquet` store."""
        return self.processed_root / "embeddings" / "meta"

    @property
    def annotation_root(self) -> Path:
        """`gff/<datasetID>/<BS>/<BS>.bakta.gff3.gz` — ⭐ the SEQUENCE source, and the only one."""
        return self.data_root / "raw" / "gff"

    # ── the run ───────────────────────────────────────────────────────────────────────────────
    @property
    def gene_universe(self) -> Path:
        """One row per gene in the set — the universe, and the genome vocabulary's own source.

        ⚠ The derived name is a CONVENTION, not a law — nuna's own `gene_universe(proc, dset, *,
        close_seq=...)` says so, because the 13,573-genome Kp cohort's file is
        `kleb_close_seq_isosrc_id0.98.parquet`, which no `{dset}` spells. `gene_universe_path`
        names it where the convention does not reach.
        """
        if self.gene_universe_path is not None:
            return self.gene_universe_path
        return self.processed_root / "mmseqs" / f"{self.set_key}_close_seq.parquet"

    @property
    def assignment(self) -> Path:
        """One row per gene, exactly once — the loader raises otherwise.

        ⚠ The derived path is where `run_nuna_model` writes a probe run. A big run writes under
        `comparison/nuna/<model>/<run>/step4/` instead, which no rule here reaches, so
        `assignment_path` names it.

        ⭐ **The parquet sidecar resolves only when the TSV is ABSENT.** `write_assignments` emits
        both, the sidecar is 16.0× smaller (measured on kp71k: 1,549,393,569 B against
        24,753,251,725 B), and `load_assignments` reads either — but preferring it where both exist
        would change which file every published catalogue is built from, silently, as a side effect
        of a path change. That swap is an input-thinning decision with its own frame-identity gate,
        not a locator default.
        """
        if self.assignment_path is not None:
            return self.assignment_path
        derived = self.analysis_root / "step4_global_merge" / "assignments" / f"{self.run_id}.tsv"
        sidecar = derived.with_suffix(".parquet")
        return sidecar if sidecar.exists() and not derived.exists() else derived

    @property
    def run_manifest(self) -> Path:
        """⚠ Its `precluster_assign` is an ABSOLUTE CSD3 path, so the chain does not walk locally."""
        return self.analysis_root / "step4_global_merge" / "runs" / f"{self.run_id}.json"

    @property
    def pfam_architectures(self) -> Path:
        """One row per protein WITH a Pfam hit — absence is missing coverage, not a lack of domains."""
        return self.analysis_root / "pfam" / f"pfam_arch_{self.set_key}.parquet"

    # ── the audit ─────────────────────────────────────────────────────────────────────────────
    def _audit(self, suffix: str) -> Path:
        return self.analysis_root / "accessory_audit" / f"{self.model_label}_{suffix}"

    @property
    def homology_waterfall(self) -> Path:
        """⚠ **Not in `meta.audit.sources`**, and it supplied three of the payload's columns."""
        return self._audit("homology_waterfall.csv")

    @property
    def pfam_concordance(self) -> Path:
        """One row per locus, full coverage — the `class_clan` verdict is READ from here, never re-derived."""
        return self._audit("pfam_concordance.tsv")

    @property
    def audit_summary(self) -> Path:
        """The graded headline, and the guard that this audit describes THIS model."""
        return self._audit("audit_summary.json")

    @property
    def cluster_table(self) -> Path:
        """The six columns nothing else carries, plus the Bacformer distance pair.

        ⚠ **The published export did not read this** — it did not exist on CSD3 on 2026-08-25, so
        the page's `bac_d_*` were recomputed from the full embedding matrix instead. Measured
        against the shipped catalogue afterwards, the two derivations agree to 5.0e-5 on all 12,104
        measurable loci, which is the payload's own 4-significant-figure rounding. So this file is a
        sound source for them — but it is a *different artifact* than the page used, and the ingest
        records which one each column came from rather than leaving the reader to assume.
        """
        return self._audit("cluster_table.parquet")

    # ── set-to-set similarity ─────────────────────────────────────────────────────────────────
    #
    # ⛔ What these replaced on 2026-09-24: `catalogue_map` / `map_sibling` / `null_baseline`, the
    # neighbourhood map and its random-MEDOID-pair baseline. They went with the medoid geometry they
    # were a picture of. `catalogue_map.py` in nuna still writes its three files; nothing reads them.
    def cluster_similarity(self, representation: str) -> Path:
        """One row per locus: the set-to-set numbers, from `build_cluster_similarity`."""
        return (
            self.analysis_root
            / "locus_browser"
            / f"{self.model_label}_cluster_similarity_{representation}.csv"
        )

    def cluster_similarity_metadata(self, representation: str) -> Path:
        """⚠ The `.meta` sidecar is the ONLY record of `rep`, `forms`, `k` and the floor — not the CSV."""
        return self.cluster_similarity(representation).with_suffix(".meta")

    def cluster_similarity_audit(self, representation: str) -> Path:
        """The run's own JSON, carrying the floor's QUARTILES, which the `.meta` reduces to a median.

        ⚠ Optional by design: a median alone still gives the card a baseline, just not its spread.
        """
        return self.cluster_similarity(representation).with_suffix(".json")

    def cluster_nearest(self, representation: str) -> Path:
        """The five nearest other loci, addressed RELATIVE to the similarity CSV.

        ⛔ Derived, not named independently — `export_payload._nearest_sibling`'s own rule. The two
        files are written together by one run and are only meaningful together, so a second flag
        would be a second chance to pair one run's medians with another run's neighbour list, and
        nothing downstream could detect it: both would be in range and plausible.
        """
        path = self.cluster_similarity(representation)
        return path.with_name(path.name.replace("_cluster_similarity_", "_cluster_nearest_"))

    # ── the projection: genomes the model never clustered ─────────────────────────────────────
    @property
    def projection_root(self) -> Path:
        """Where nuna wrote the placements for genomes OUTSIDE this pangenome's 100.

        ⚠ nuna names this directory `{set}_{model_label}` with its own `--model-label`, which is the
        part after the set token — so the directory name is exactly THIS class's `model_label`, and
        the two conventions meet without a rule that has to be kept in step.

        ⛔ Not in `required()`. A catalogue is complete without any projected genome, and making the
        whole load depend on an optional feature's artifacts would be the wrong dependency.
        """
        return self.analysis_root / "locus_projection" / self.model_label

    def projected_genome_table(self, sample_id: str) -> Path:
        """One row per placed gene, with its ±5 window and the arrangement it matched."""
        return self.projection_root / f"{sample_id}_placed_windows.tsv"

    @property
    def placement_summary(self) -> Path:
        """The rule, k, n, and the assignment sha256 the placements were computed against."""
        return self.projection_root / "placement_summary.json"

    @property
    def projection_manifest(self) -> Path:
        """The search's own manifest — the assignment it searched and that file's digest."""
        return self.projection_root / "manifest.json"

    # ── the oracle, which is never an input ───────────────────────────────────────────────────
    @property
    def published_payload(self) -> Path:
        """⛔ **The acceptance ORACLE.** Reading it into the database would make the test circular.

        ⚠ Named from `set_key`, which is the SPECIES' data set. nuna's exporter names a payload
        `locus_browser_<dset>_<label>.json`, where `--dset` is the STATIC publication key — so the two
        spellings diverge for any catalogue whose page key is not just the species: this would resolve
        `locus_browser_ecoli_ecoli_nuna5_…` where such an export writes `locus_browser_ecoli-nuna5_…`.

        Correct as it stands, because only STATICALLY published catalogues have an oracle and those are
        `nuna4` only (see `published_catalogues`), where `set_key == dset`. The property is optional and
        never raises, so a catalogue with no shipped payload simply reports `False`. The day a nuna5
        payload is published statically, this takes that row's `dset` — not a glob, which is what
        nuna's own `publish_site.find_payload` uses and why its keys must be prefix-free.
        """
        return self.data_root / "browser" / f"locus_browser_{self.set_key}_{self.model_label}.json"

    # ── per-genome ────────────────────────────────────────────────────────────────────────────
    def per_genome(self, sample_id: str, kind: str) -> Path:
        """`kind` ∈ `meta | product | strand | dbxref | noncoding`."""
        return self.per_genome_root / f"{sample_id}_{kind}.parquet"

    def per_genome_present(self, kind: str) -> frozenset[str]:
        """Which genomes have a `{BS}_{kind}.parquet` (a regular file) — ONE scan of `per_genome_root`.

        ⛔ Never a per-genome `is_file()`: see `listed_files`. The directory holds five kinds per genome,
        ~357k entries at 71,435 genomes, and the scan reads their names and `d_type`s, not their inodes.
        """
        return listed_files(self.per_genome_root, f"_{kind}.parquet")

    def genome_annotation(self, sample_id: str, bakrep_dataset_id: str) -> Path:
        """The gzipped Bakta GFF — annotation AND, in its `##FASTA` block, the bases."""
        return self.annotation_root / bakrep_dataset_id / sample_id / f"{sample_id}.bakta.gff3.gz"

    # ── the published sequence files ──────────────────────────────────────────────────────────
    def _nseq_root_or_refuse(self) -> Path:
        if self.nseq_root is None:
            raise ValueError(
                "no .nseq tree was named (`nseq_root` / `--nseq-dir`), so there is no per-genome "
                ".nseq or .loci path to give. A caller that wants one must have been told where they "
                "are; a guessed path would resolve to nothing and read as a cohort with no sequences."
            )
        return self.nseq_root

    def per_genome_nseq(self, sample_id: str) -> Path:
        """`<nseq_root>/<BS>.nseq` — the genome: contig names, the gene table, the 2-bit bases."""
        return self._nseq_root_or_refuse() / f"{sample_id}.nseq"

    def per_genome_loci(self, sample_id: str) -> Path:
        """`<nseq_root>/<BS>.loci` — the model: which catalogue locus each gene landed in."""
        return self._nseq_root_or_refuse() / f"{sample_id}.loci"

    def nseq_present(self) -> frozenset[str]:
        """Which genomes have a `.nseq` under `nseq_root`, from ONE directory listing.

        ⛔ Never a per-genome `exists()`. A cold `stat` on Lustre costs ~28–48 ms, and the cohort this
        exists for has 71,435 genomes — an hour of metadata latency to learn what one `readdir` says in
        a second. Callers hold the set and look up in it.

        ⚠ `None` root → an empty set, not an error: the probe names no tree and loads as it always has.
        A NAMED root that is not a directory IS refused — that is a misconfiguration, and letting it
        read as "no genome has a .nseq" is the silent success this locator exists to prevent.
        """
        if self.nseq_root is None:
            return frozenset()
        if not self.nseq_root.is_dir():
            raise FileNotFoundError(
                f"the .nseq tree was named as {self.nseq_root} and it is not a directory. Refused "
                "rather than read as a cohort with no sequences."
            )
        suffix = ".nseq"
        return frozenset(name[: -len(suffix)] for name in os.listdir(self.nseq_root) if name.endswith(suffix))

    # ── checking ──────────────────────────────────────────────────────────────────────────────
    def required(self) -> list[tuple[str, Path, str]]:
        """`(role, path, what it supplies)` for everything that must exist to ingest at all."""
        items = [
            ("gene_universe", self.gene_universe, "the gene universe and the genome vocabulary"),
            ("assignment", self.assignment, "gene → locus, and the locus set"),
            ("run_manifest", self.run_manifest, "pangenome_step and the lineage edges"),
            ("pfam_architectures", self.pfam_architectures, "locus.pfam_annotated_member_count"),
            ("homology_waterfall", self.homology_waterfall, "collapse_tier, resolved_threshold, the ESM pair"),
            ("pfam_concordance", self.pfam_concordance, "pfam_architecture_count, pfam_concordance_class"),
            ("audit_summary", self.audit_summary, "pangenome_evaluation and the model-identity guard"),
        ]
        for rep in REPRESENTATIONS:
            items += [
                (f"cluster_similarity[{rep}]", self.cluster_similarity(rep), "locus_similarity"),
                (f"similarity_meta[{rep}]", self.cluster_similarity_metadata(rep),
                 "the form, the gene-pair floor and the measurable denominator"),
                (f"cluster_nearest[{rep}]", self.cluster_nearest(rep), "locus_nearest_locus"),
            ]
        return items

    def optional(self) -> list[tuple[str, Path, str]]:
        """Present or absent is a FACT about the run, recorded rather than treated as a failure."""
        return [
            ("cluster_table", self.cluster_table, "u50 impurity/coverage and the medoid member"),
            *[
                (f"similarity_audit[{rep}]", self.cluster_similarity_audit(rep),
                 "the gene-pair floor's QUARTILES; without it the card draws a bare tick")
                for rep in REPRESENTATIONS
            ],
            ("published_payload", self.published_payload, "the acceptance oracle — never an input"),
        ]

    def verify(self, *, roles: frozenset[str] | None = None) -> dict[str, bool]:
        """Raise naming EVERY absent required artifact; return which optional ones are present.

        ``roles`` narrows the requirement to the artifacts one stage actually reads. ⚠ It is named by
        the caller, never inferred: the genome layer is model-independent and reads only the gene
        universe, so demanding the catalogue's artifacts of it made a 71,435-genome load wait a day
        for a producer it does not depend on. A role that `required()` does not list is refused —
        a typo here would otherwise verify nothing and say everything was present.
        """
        required = self.required()
        if roles is not None:
            unknown = roles - {role for role, _, _ in required}
            if unknown:
                raise ValueError(f"not required artifacts: {sorted(unknown)}")
            required = [item for item in required if item[0] in roles]
        missing = [item for item in required if not item[1].exists()]
        if missing:
            raise MissingArtifacts(missing)
        return {role: path.exists() for role, path, _ in self.optional()}
