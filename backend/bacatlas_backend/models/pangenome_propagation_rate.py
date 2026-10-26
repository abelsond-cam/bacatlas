"""How often a node's own call is unanimous, per catalogue · vocabulary · prevalence band.

⭐ **Precomputed at ingest, never aggregated per request** — the same stance and the same wording as
`pangenome_genome_locus_count`, for a heavier reason. The equivalent at request time is
`within_node_propagation.compute_propagation`, which reads **every annotated gene row of the
catalogue**: measured at **0.22–1.94 s** per (catalogue, vocabulary) on the published catalogues, so
a cold *E. coli* Function tab paid ~5.8 s across the four before a process cache absorbed it. At the
80,000-genome design target the same scan is over ~412 M rows.

⛔ **These twenty rows per catalogue are REDUNDANT, and that is the point.** Every one of them is a
`GROUP BY` over `locus.propagation_verdict_*`, which is itself stored. The redundancy is not an
oversight to be optimised away later — it is the stage's own oracle: the table, the verdict columns
and the live gene-grain instrument are three representations of one measurement, and
`tests/test_locus_annotation_rollups.py` requires all three to agree. A single stored number has
nothing to be checked against; two have.

⚠ **It cannot be a view or a generated column** for exactly the reason it exists:
`tests/test_route_cost.py`'s `_catalogue_aggregations` flags any request-time statement matching
`GROUP BY … FROM locus`, and `test_ONLY_the_shell_aggregates_over_the_catalogue` asserts the locus
and keystroke routes issue none.
"""

from __future__ import annotations

from sqlalchemy import BigInteger, CheckConstraint, ForeignKey, Integer
from sqlalchemy.orm import Mapped, mapped_column

from bacatlas_backend.database import Base
from bacatlas_backend.models.enumerations import AnnotationKind, PrevalenceBand


class PangenomePropagationRate(Base):
    """One catalogue's unanimity tally for one vocabulary in one prevalence band.

    ⛔ **`checkable` and `one_gene` are DIFFERENT POPULATIONS and neither is a subset of the other.**
    `checkable` is what the rate is *measured on* — nodes with two or more annotated genes and a rung
    all of them state. `one_gene` is what the rate is *applied to* — nodes resting on a single
    annotated gene, unanimous by construction and therefore deliberately excluded from the
    measurement. They are stored side by side because the gap between them is the caveat, and a page
    that printed the rate without the population it is extrapolated to would be hiding it.

    ⛔⛔ **`checkable = 0` is a measured zero and NOT a rate of 0.0.** `PropagationRate.rate` returns
    `None` below `within_node_propagation.MIN_NODES`, and RARE sits there in every catalogue measured
    so far — with *hundreds* of `one_gene` nodes beside it. *Not measurable* and *measured and bad*
    are different statements, and the enum-less zero here is the first place that distinction can be
    lost.

    ⚠ **Nodes in neither column are not an error.** `NO_COVERAGE` (no member gene carries the
    vocabulary) and `NOT_COMPARABLE` (several do and no shared rung exists) contribute to neither,
    so `checkable + one_gene` is the count of nodes the rate is *about*, never the band's node count.
    `locus.propagation_verdict_*` is where the other two are recoverable.
    """

    __tablename__ = "pangenome_propagation_rate"
    __table_args__ = (
        # ⚠ No secondary index, and that is deliberate rather than omitted: the only FK is
        # `pangenome_id`, which LEADS the primary key, so a catalogue drop uses it. The omission that
        # cost an 8-minute stall on `gene_locus_membership` was a FK that did *not* lead its key.
        CheckConstraint("unanimous >= 0 AND unanimous <= checkable", name="unanimous_within_checkable"),
        CheckConstraint("checkable >= 0 AND one_gene >= 0", name="populations_are_not_negative"),
    )

    #: ⚠ BigInteger, because `pangenome.pangenome_id` is.
    pangenome_id: Mapped[int] = mapped_column(
        BigInteger, ForeignKey("pangenome.pangenome_id", ondelete="CASCADE"), primary_key=True
    )
    #: ⛔ One of the FOUR kinds `within_node_propagation` supports — COG, GO slim, EC, KEGG. The
    #: enum carries seven; `gene_symbol`, `protein_product` and `pfam_architecture` have no rung
    #: ladder and no propagation question, so they have no row rather than a zero row.
    annotation_kind: Mapped[AnnotationKind] = mapped_column(primary_key=True)
    prevalence_band: Mapped[PrevalenceBand] = mapped_column(primary_key=True)

    #: Nodes the rate is MEASURED on: ≥ 2 annotated genes and a rung every one of them states.
    checkable: Mapped[int] = mapped_column(Integer, nullable=False)
    #: Of those, the nodes whose genes agreed at that rung.
    unanimous: Mapped[int] = mapped_column(Integer, nullable=False)
    #: Nodes the rate is APPLIED to: exactly one annotated gene. Never in `checkable`.
    one_gene: Mapped[int] = mapped_column(Integer, nullable=False)

    def __repr__(self) -> str:
        return (
            f"<PangenomePropagationRate {self.annotation_kind} {self.prevalence_band} "
            f"{self.unanimous}/{self.checkable} +{self.one_gene}>"
        )
