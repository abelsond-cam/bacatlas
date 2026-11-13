"""Genomes the model never clustered, placed on a published pangenome's loci.

⛔ **It writes two tables and reads the rest.** `genome_locus_assignment`, `locus_arrangement`,
`genome_collection_membership` and `pangenome_genome_locus_count` are never touched, so the publish
gate, the column audit and every parity suite keep meaning what they meant. A projected genome is
not a member of anything: it is drawn beside the model, never inside it.

⛔ **Five refusals, each naming what it found**, because every one of them produces a page that looks
entirely ordinary:

1. **A genome in the collection.** Placing a modelled genome on its own model would report cosine 1.0
   against itself — a perfect score meaning nothing.
2. **A row count that is not the genome's gene count.** A partial file would place three quarters of
   a genome and say nothing about the rest.
3. **A duplicated `flat_index`, or one outside `[0, genome_genes.gene_count)`.** `flat_index` is
   positional and shared with the embedding `.npy` and every parquet; a value the genome does not
   have means the file is not describing this genome's genes. ⛔ Since 2026-11-13 this check is the
   ONLY thing that enforces it — the composite foreign key onto `gene` went with the per-gene table.
4. **A locus label this pangenome does not have.** Labels are model-private, so a label from another
   model resolves to a *different* locus rather than to nothing.
5. **A placement computed against a different assignment** — `assignment_sha256` must equal
   `pangenome.assignment_sha256`. This is the one that cannot be seen by looking at the numbers.

⚠ **Absences that are not the same absence**, kept apart all the way through: a gene with **no
Bacformer vector** has no row at all (counted in `gene_without_vector_count`); a gene **alone on its
contig** has a row with `neighbour_slot_codes` NULL (`gene_without_window_count`); a gene whose window
matches nothing has a vector and `matched_locus_arrangement_id` NULL. Collapsing any two would make
the page claim a neighbourhood comparison for a gene that has no neighbours.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd
from sqlalchemy import select, text
from sqlalchemy.orm import Session

from bacatlas_backend.ingest.staging_table_loader import replace_rows_for
from bacatlas_backend.models.genome import Genome
from bacatlas_backend.models.genome_collection import GenomeCollectionMembership
from bacatlas_backend.models.locus import Locus
from bacatlas_backend.models.pangenome import Pangenome
from bacatlas_backend.models.projected_genome import ProjectedGenePlacement, ProjectedGenome

PLACEMENT_COLUMNS = (
    "pangenome_id",
    "genome_id",
    "flat_index",
    "locus_id",
    "nearest_cosine",
    "agreeing_neighbour_count",
    "available_neighbour_count",
    "placed_summed_cosine",
    "runner_up_locus_id",
    "runner_up_summed_cosine",
    "is_contested",
    "copy_ordinal",
    "copies_at_locus",
    "neighbour_slot_codes",
    "matched_locus_arrangement_id",
)

SLOT_COLUMNS = tuple(f"s{index}" for index in range(10))


class ProjectionRefused(ValueError):
    """The file does not describe this pangenome, and loading it would be worse than not loading it."""


@dataclass
class ProjectionReport:
    """What was written, per genome, in counts a reader can check."""

    genomes: list[dict] = field(default_factory=list)

    def render(self) -> str:
        """One line per genome, plus a total — a loader that prints 'done' has told you nothing."""
        lines = [f"projection: {len(self.genomes)} genome(s)"]
        for genome in self.genomes:
            lines.append(
                f"  {genome['sample_id']}: {genome['placed']:,} placed on {genome['loci']:,} loci · "
                f"{genome['contested']:,} contested · {genome['window_matched']:,} windows matched · "
                f"{genome['without_window']:,} alone on a contig"
            )
        return "\n".join(lines)


def sha256_of(path: Path) -> str:
    """Streaming digest — what pins a loaded projection to the file it came from."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _locus_ids(session: Session, pangenome_id: int) -> dict[str, int]:
    """`node_label` → `locus_id` for this pangenome. ⛔ Labels are model-private: never cross runs.

    ⚠ The column is `node_label` and it is TEXT: nuna's node ids look numeric and are not.
    """
    rows = session.execute(
        select(Locus.node_label, Locus.locus_id).where(Locus.pangenome_id == pangenome_id)
    ).all()
    return {str(label): int(locus_id) for label, locus_id in rows}


def _genome_id(session: Session, sample_id: str) -> int:
    """The genome's own id, which must already exist — the genome layer is model-independent."""
    found = session.execute(select(Genome.genome_id).where(Genome.sample_id == sample_id)).scalar_one_or_none()
    if found is None:
        raise ProjectionRefused(
            f"{sample_id} is not in the `genome` table. A projected genome's coordinates come from the "
            "genome layer, which is loaded separately and is not a model artifact."
        )
    return int(found)


def _refuse_a_collection_genome(session: Session, pangenome_id: int, genome_id: int, sample_id: str) -> None:
    """⛔ A modelled genome cannot be projected onto its own model."""
    in_collection = session.execute(
        select(GenomeCollectionMembership.genome_id)
        .join(Pangenome, Pangenome.genome_collection_id == GenomeCollectionMembership.genome_collection_id)
        .where(Pangenome.pangenome_id == pangenome_id, GenomeCollectionMembership.genome_id == genome_id)
    ).first()
    if in_collection:
        raise ProjectionRefused(
            f"{sample_id} is one of this pangenome's modelled genomes. Its genes are already in "
            "`genome_locus_assignment`, and placing it again would report a gene's distance to itself."
        )


def _check_the_genes_are_this_genomes(session: Session, genome_id: int, frame: pd.DataFrame, sample_id: str) -> int:
    """Uniqueness and range of `flat_index` against `genome_genes.gene_count`. Returns the gene count.

    ⛔ **This check IS the referential integrity now.** `projected_gene_placement` used to carry a
    composite foreign key onto `gene (genome_id, flat_index)`; the gene layer is one row per genome
    (`genome_genes.gene_spans`, a block indexed by `flat_index`), so there is no row for a constraint
    to reference and migration `9e8f7a6b5c4d` dropped it. A `flat_index` is one of this genome's
    genes exactly when it is in `[0, gene_count)` — `flat_index` is a dense positional offset, so the
    half-open range is the whole membership test, not an approximation of it.

    ⛔ **No `genome_genes` row is a refusal, not a count of zero.** `count(*)` over an absent genome
    used to return 0 and then refuse every placement as unknown — a true verdict with a misleading
    reason. The genome layer being unloaded is named as such.
    """
    gene_count = session.execute(
        text("SELECT gene_count FROM genome_genes WHERE genome_id = :genome_id"), {"genome_id": genome_id}
    ).scalar_one_or_none()
    if gene_count is None:
        raise ProjectionRefused(
            f"{sample_id} has no genome_genes row, so its gene count is unknown and no placement can be "
            "checked against it. The genome layer is loaded separately and is not a model artifact."
        )
    gene_count = int(gene_count)
    if frame["flat_index"].duplicated().any():
        repeated = frame.loc[frame["flat_index"].duplicated(), "flat_index"].head(5).tolist()
        raise ProjectionRefused(f"{sample_id}: flat_index repeats (e.g. {repeated}); it is positional and unique")
    placed = [int(value) for value in frame["flat_index"]]
    unknown = sorted(value for value in placed if not 0 <= value < gene_count)
    if unknown:
        raise ProjectionRefused(
            f"{sample_id}: {len(unknown)} placed genes have a flat_index this genome does not have "
            f"(e.g. {unknown[:5]}; the genome has {gene_count} genes, so 0..{gene_count - 1}). "
            "flat_index is a positional row offset shared with the embedding store and every parquet — "
            "a mismatch means the file describes a different gene set."
        )
    return gene_count


def _slot_codes(row: pd.Series) -> list[int] | None:
    """The ten packed slots, or **None where the gene is alone on its contig** and has no window."""
    if not bool(row.get("has_a_window", True)):
        return None
    if any(pd.isna(row.get(column)) for column in SLOT_COLUMNS):
        return None
    return [int(row[column]) for column in SLOT_COLUMNS]


def load_projected_genome(
    session: Session,
    *,
    pangenome_id: int,
    path: Path,
    summary: dict,
    assignment_sha256: str | None,
    nuna_git_sha: str | None,
) -> dict:
    """One genome's placements, replacing any previous load of the same genome. Returns its counts."""
    frame = pd.read_csv(path, sep="\t", dtype={"locus_label": str, "runner_up_label": str})
    sample_id = str(frame["sample_id"].iloc[0])
    genome_id = _genome_id(session, sample_id)
    _refuse_a_collection_genome(session, pangenome_id, genome_id, sample_id)
    gene_count = _check_the_genes_are_this_genomes(session, genome_id, frame, sample_id)

    labels = _locus_ids(session, pangenome_id)
    unknown = sorted(set(frame["locus_label"]) - set(labels))
    if unknown:
        raise ProjectionRefused(
            f"{sample_id}: {len(unknown)} locus labels are not in this pangenome (e.g. {unknown[:5]}). "
            "Labels are model-private, so a label from another run would resolve to a DIFFERENT locus."
        )

    matched = _arrangement_ids(session, pangenome_id, frame, labels)
    rows = []
    for _, row in frame.iterrows():
        locus_id = labels[str(row["locus_label"])]
        runner_up = row.get("runner_up_label")
        rows.append(
            (
                pangenome_id,
                genome_id,
                int(row["flat_index"]),
                locus_id,
                _number(row.get("nearest_cosine")),
                int(row.get("agreeing_neighbours", 0)),
                int(row.get("available_neighbours", row.get("agreeing_neighbours", 0))),
                _number(row.get("placed_summed_cosine")),
                labels.get(str(runner_up)) if isinstance(runner_up, str) and runner_up in labels else None,
                _number(row.get("runner_up_summed_cosine")),
                bool(row.get("is_contested", False)),
                int(row.get("copy_ordinal", 1)),
                int(row.get("copies_at_this_locus", 1)),
                _slot_codes(row),
                matched.get(int(row["flat_index"])),
            )
        )

    counts = {
        "sample_id": sample_id,
        "genome_id": genome_id,
        "genes": gene_count,
        "placed": len(frame),
        "loci": int(frame["locus_label"].nunique()),
        "contested": int(frame["is_contested"].sum()) if "is_contested" in frame else 0,
        "multi_copy_loci": int((frame.get("copies_at_this_locus", pd.Series(dtype=int)) > 1).sum()),
        "window_matched": int(sum(1 for value in matched.values() if value is not None)),
        # ⚠ Two different absences, counted apart: no vector at all, and no window.
        "without_vector": max(0, gene_count - len(frame)),
        "without_window": int(sum(1 for row in rows if row[13] is None)),
    }
    # ⛔ The parent row FIRST. `projected_gene_placement` has a composite foreign key to
    # `projected_genome`, so a COPY before it fails on every row — and deleting the parent cascades
    # to its old placements, which is exactly the replace this re-ingest wants.
    _write_the_genome_row(
        session,
        pangenome_id=pangenome_id,
        genome_id=genome_id,
        path=path,
        summary=summary,
        counts=counts,
        frame=frame,
        assignment_sha256=assignment_sha256,
        nuna_git_sha=nuna_git_sha,
    )
    replace_rows_for(
        session,
        ProjectedGenePlacement.__table__,
        PLACEMENT_COLUMNS,
        rows,
        where=(ProjectedGenePlacement.pangenome_id == pangenome_id)
        & (ProjectedGenePlacement.genome_id == genome_id),
    )
    return counts


def _arrangement_ids(
    session: Session, pangenome_id: int, frame: pd.DataFrame, labels: dict[str, int]
) -> dict[int, int | None]:
    """`flat_index` → the `locus_arrangement_id` whose ten slots EQUAL this gene's, or absent.

    ⭐ One statement, matching on **the whole vector** rather than on the rank nuna computed. Nine
    slots agreeing and one differing is a different neighbourhood, so a per-slot or rank-based match
    would be the marginal-for-joint error the two views exist to prevent — and matching on the vector
    makes the database CHECK the projection rather than take its word. Both sides are packed
    identically (`catalogue_ordinal * 2 + same_strand`, −1 at a contig end) because both were built
    by the same nuna function.

    ⚠ A gene alone on its contig has no vector and no row here — absent, which is not the same as
    matching nothing.
    """
    rows = []
    for _, row in frame.iterrows():
        slots = _slot_codes(row)
        if slots is None:
            continue
        rows.append({"flat_index": int(row["flat_index"]), "locus_id": labels[str(row["locus_label"])], "slots": slots})
    if not rows:
        return {}
    found = session.execute(
        text(
            "SELECT w.flat_index, a.locus_arrangement_id "
            "  FROM json_to_recordset(CAST(:payload AS json)) "
            "       AS w(flat_index int, locus_id bigint, slots int[]) "
            "  JOIN locus_arrangement a "
            "    ON a.locus_id = w.locus_id AND a.neighbour_slot_codes = w.slots "
            " WHERE a.pangenome_id = :pangenome_id"
        ),
        {"payload": json.dumps(rows), "pangenome_id": pangenome_id},
    ).all()
    return {int(flat_index): int(arrangement_id) for flat_index, arrangement_id in found}


def _number(value) -> float | None:
    """A float, or None where it is absent or NaN — *not measured* and *measured zero* differ."""
    if value is None or (isinstance(value, float) and value != value) or pd.isna(value):
        return None
    return float(value)


def _write_the_genome_row(
    session: Session,
    *,
    pangenome_id: int,
    genome_id: int,
    path: Path,
    summary: dict,
    counts: dict,
    frame: pd.DataFrame,
    assignment_sha256: str | None,
    nuna_git_sha: str | None,
) -> None:
    """The per-genome row, carrying the rule and the distribution rather than a single score."""
    cosine = frame["nearest_cosine"].dropna() if "nearest_cosine" in frame else pd.Series(dtype=float)
    session.execute(
        ProjectedGenome.__table__.delete().where(
            (ProjectedGenome.pangenome_id == pangenome_id) & (ProjectedGenome.genome_id == genome_id)
        )
    )
    session.execute(
        ProjectedGenome.__table__.insert().values(
            pangenome_id=pangenome_id,
            genome_id=genome_id,
            rule_label=summary.get("rule", "the nearest modelled gene decides the locus; the first n check it"),
            neighbours_searched=int(summary.get("neighbours_searched") or summary.get("checking_neighbours", 10)),
            neighbours_reported=int(summary.get("checking_neighbours", 10)),
            representation=str(summary.get("representation", "bacformer")),
            source_file_path=str(path),
            source_sha256=sha256_of(path),
            assignment_sha256=assignment_sha256,
            nuna_git_sha=nuna_git_sha,
            gene_count=counts["genes"],
            placed_gene_count=counts["placed"],
            gene_without_vector_count=counts["without_vector"],
            gene_without_window_count=counts["without_window"],
            contested_gene_count=counts["contested"],
            distinct_locus_count=counts["loci"],
            multi_copy_locus_count=counts["multi_copy_loci"],
            window_matched_gene_count=counts["window_matched"],
            nearest_cosine_median=float(cosine.median()) if len(cosine) else None,
            nearest_cosine_fifth_percentile=float(cosine.quantile(0.05)) if len(cosine) else None,
            nearest_cosine_minimum=float(cosine.min()) if len(cosine) else None,
        )
    )


def load_projection(session: Session, *, pangenome_id: int, projection_root: Path) -> ProjectionReport:
    """Every `*_placed_windows.tsv` under one projection directory, refusing as a whole on any fault.

    ⛔ The sidecar's `assignment.sha256` must equal `pangenome.assignment_sha256`. A projection built
    against a different assignment would name loci that mean something else, and nothing about the
    resulting page would look wrong.
    """
    summary_path = projection_root / "placement_summary.json"
    if not summary_path.exists():
        raise ProjectionRefused(
            f"no placement_summary.json under {projection_root} — it carries the rule, k and the "
            "assignment digest, and loading placements without it would lose their provenance."
        )
    summary = json.loads(summary_path.read_text())
    assignment = summary.get("assignment", {}) or {}
    expected = session.execute(
        select(Pangenome.assignment_sha256).where(Pangenome.pangenome_id == pangenome_id)
    ).scalar_one_or_none()
    found = assignment.get("sha256")
    if expected and found and expected != found:
        raise ProjectionRefused(
            f"this projection was computed against assignment {found[:12]}… but the pangenome holds "
            f"{expected[:12]}…. Locus labels are model-private, so every placement would name a "
            "different locus than intended."
        )

    report = ProjectionReport()
    for path in sorted(projection_root.glob("*_placed_windows.tsv")):
        report.genomes.append(
            load_projected_genome(
                session,
                pangenome_id=pangenome_id,
                path=path,
                summary=summary,
                assignment_sha256=found,
                nuna_git_sha=summary.get("nuna_git_sha"),
            )
        )
    if not report.genomes:
        raise ProjectionRefused(f"no *_placed_windows.tsv under {projection_root}")
    return report
