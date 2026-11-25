"""What annotation-transfer PREDICTS for a locus — one row per (locus, claim), never mixed with what was observed.

⛔⛔ **Its own table, never `locus_annotation_entry`** (`nuna/docs/model_evaluation/annotation_transfer_gains.md`
§8). That table holds what a locus's own genes were OBSERVED to carry. A predicted annotation sharing it
is one join away from being read as evidence, and nothing downstream would look wrong.

⭐ **One row per (locus, claim), because that is the unit the producer writes and the unit the accuracy
is about.** nuna's `summarise_transfer_gains.py` writes `{catalogue}_served_annotations.tsv` with one row
per (node, claim): a SET of classes (`predicted_classes`, `|`-joined), one `p_serve`, one donor. The
calibrated accuracy is the probability that the predicted set INTERSECTS the locus's true set
(`nuna.tl.annotation_transfer.measure_calibration.agrees`) — a property of the set, not of any one class
in it. Exploding the set into one row per class would repeat that number beside every class, where it
reads as "this class is right with probability p", which is not what was measured; and a row count would
stop being a prediction count (62,344 class rows against 39,009 predictions on the probe).

⚠ **A GO prediction is a CLOSED LINEAGE** — the donor's voted slim class plus every slim class above it
(`derive_node_claims.go_claims`), so `term_values` holds e.g. *plasma membrane* AND *membrane*. The
accuracy is scored on that whole set, so a match at the broadest class alone counts as correct. The
ingest stores the set most-specific-first (`ingest_served_annotations.order_lineage`) so the page can
lead with the class the donor actually states and still show what the number covers.

⚠ **A served locus does NOT always carry nothing for the claim** (review of 9cfa481, measured on the
probes). The producer serves a node whose OWN claim it reads as empty, and it reads GO as the rank-0 slim
class with namespace roots stripped AFTER the vote. nuna's `functional.top_go_slim` breaks a vote tie by
ascending accession, the roots (`GO:0003674`, `GO:0008150`, `GO:0005575`) sort first, so a node whose
genes all carry a root and a real class elects the root, has it stripped, and becomes a recipient — 426
of the 35,588 served GO predictions on the two probes sit beside an observed non-root class in the same
namespace. Whether to fix the producer or filter here is David's decision; BacAtlas serves what the
producer wrote and says nothing on the page that this contradicts.

Columns the producer writes and this table does NOT keep, and why:

* `n_recipient_genes`, `n_donor_genes` — checked at ingest against `locus.member_gene_count` on both
  sides. Equal by then, so storing them would be a second copy of a stored column. ⚠ That check is a
  size-histogram check, not an identity check (labels are near size order); identity is the file's
  `catalogue_run_id`, which must equal the catalogue's `run_id`.
* `catalogue_run_id`, `annotation_model_version` — the first is checked and dropped; the second is kept.
* `donor_catalogue`, `direction` — derivable from the donor locus's own pangenome.
* `rank_esm`, `rank_bac`, `both_directions`, `directions_agree` — model-reading diagnostics (the donor's
  rank in each space; the cross-direction winner rule's price). Nothing on the page reads them, and they
  stay in the artifact.
* `n_predicted_classes` — `cardinality(term_values)`.
"""

from __future__ import annotations

from dataclasses import dataclass

from sqlalchemy import BigInteger, CheckConstraint, Float, ForeignKey, Index, SmallInteger, String, Text, text
from sqlalchemy.dialects.postgresql import ARRAY
from sqlalchemy.orm import Mapped, mapped_column

from bacatlas_backend.database import Base
from bacatlas_backend.models.enumerations import ServedAnnotationKind

#: ⛔ The tag pattern, enforced by a CHECK on both columns that carry a version. Both version ladders have
#: a "v0.1" (nuna's clustering ladder and this one), so a bare `v0.1` is refused by the database itself,
#: not only by the ingest (`nuna/docs/annotation_model_versions.md`).
ANNOTATION_MODEL_VERSION_PATTERN = r"^at-v[0-9]+\.[0-9]+$"

#: ⚠ The medoid cosines' bound, with float32 slack. The medoids are L2-normalised in float32 (norms
#: 0.9999998-1.0000002 on the probe), so a cosine can exceed 1 by ~1e-6: the probe's E. coli file holds
#: `esm_medoid_sim = 1.000001`. A bound of exactly 1 would refuse a correct row.
COSINE_BOUND = 1.001


@dataclass(frozen=True)
class AnnotationModelVersion:
    """One version of the annotation model, as BacAtlas needs to know it.

    ⛔ **A file is loaded only if its version is registered here**, so a new version is a deliberate code
    change (its cut and its claims stated beside it), never something a file can introduce by spelling a
    new tag.
    """

    #: `at-v0.1` — the tag in file names, columns and the database.
    tag: str
    #: The probability at or above which a prediction is SERVED. ⚠ nuna owns the number
    #: (`score_transfer_calls.SERVE_CUT`); it is restated here only so the ingest can refuse a row below it.
    serve_cut: float
    #: The producer's claim keys this version predicts.
    claims: frozenset[str]

    @property
    def label(self) -> str:
        """How a reader sees it: `annotation-transfer v0.1`, never a bare `v0.1`."""
        return f"annotation-transfer {self.tag.removeprefix('at-')}"


#: ⭐ annotation-transfer v0.1 (David, 2026-11-14): the medoid model, both representations, one donor,
#: served at p ≥ 0.80, COG functional categories and GO slim in all three namespaces.
KNOWN_ANNOTATION_MODEL_VERSIONS: dict[str, AnnotationModelVersion] = {
    "at-v0.1": AnnotationModelVersion(tag="at-v0.1", serve_cut=0.80, claims=frozenset({"cog", "mf", "bp", "cc"})),
}

#: ⛔ The producer's claim key → how BacAtlas addresses it. The GO namespace index is the
#: `locus_annotation_entry.gene_ontology_namespace` contract (0 molecular function, 1 biological process,
#: 2 cellular component — `enumerations.GENE_ONTOLOGY_NAMESPACE_NAMES`). ⚠ A transposed entry here would
#: file every cellular-component prediction under molecular function and the page would look normal, so
#: the ingest checks each GO class's own namespace in the reference against the one this map assigns.
CLAIM_ADDRESS: dict[str, tuple[ServedAnnotationKind, int | None]] = {
    "cog": (ServedAnnotationKind.COG_FUNCTIONAL_CATEGORY, None),
    "mf": (ServedAnnotationKind.GENE_ONTOLOGY_SLIM, 0),
    "bp": (ServedAnnotationKind.GENE_ONTOLOGY_SLIM, 1),
    "cc": (ServedAnnotationKind.GENE_ONTOLOGY_SLIM, 2),
}

def annotation_model_label(tag: str | None) -> str | None:
    """`at-v0.1` → `annotation-transfer v0.1`; a tag this build does not know keeps its own spelling."""
    if tag is None:
        return None
    known = KNOWN_ANNOTATION_MODEL_VERSIONS.get(tag)
    return known.label if known is not None else tag


class ServedLocusAnnotation(Base):
    """One served prediction: a locus, one claim, the predicted class set, its accuracy, and its donor."""

    __tablename__ = "served_locus_annotation"
    __table_args__ = (
        # ⛔ ONE prediction per (locus, claim), whatever the version — a catalogue carries one version at a
        # time and the ingest replaces it whole. ⚠ The namespace is NULL for COG and Postgres treats NULLs
        # as distinct in a unique index, so it is COALESCEd (the `locus_annotation_entry` pattern); the
        # CHECK below is what makes the coalesced value meaningful. It also serves the read path, which
        # selects by `locus_id`.
        Index(
            "uq_served_locus_annotation__locus_kind_namespace",
            "locus_id",
            "annotation_kind",
            text("COALESCE(gene_ontology_namespace, -1)"),
            unique=True,
        ),
        # ⚠ `donor_locus_id` is a SECOND FK to `locus` and has no index without this one: deleting a
        # catalogue would then scan this table once per locus (the 8-minute stall on
        # `gene_locus_membership` was exactly that omission).
        Index("ix_served_locus_annotation__donor_locus_id", "donor_locus_id"),
        CheckConstraint(
            "(annotation_kind = 'GENE_ONTOLOGY_SLIM') = (gene_ontology_namespace IS NOT NULL)",
            name="namespace_matches_kind",
        ),
        CheckConstraint(
            "gene_ontology_namespace IS NULL OR gene_ontology_namespace BETWEEN 0 AND 2",
            name="namespace_is_known",
        ),
        CheckConstraint(
            "cardinality(term_values) >= 1 AND array_position(term_values, NULL) IS NULL",
            name="terms_are_stated",
        ),
        CheckConstraint("estimated_accuracy >= 0 AND estimated_accuracy <= 1", name="accuracy_is_a_probability"),
        CheckConstraint(
            f"esm_medoid_similarity BETWEEN -{COSINE_BOUND} AND {COSINE_BOUND} "
            f"AND bacformer_medoid_similarity BETWEEN -{COSINE_BOUND} AND {COSINE_BOUND}",
            name="similarities_are_cosines",
        ),
        CheckConstraint(
            f"annotation_model_version ~ '{ANNOTATION_MODEL_VERSION_PATTERN}'", name="version_is_tagged"
        ),
        # ⚠ NO `nan_guards` here, and the omission is measured rather than forgotten: Postgres orders NaN
        # ABOVE every number, so `NaN <= 1` and `NaN BETWEEN … AND 1.001` are both FALSE and the two range
        # CHECKs above already refuse a NaN (`test_served_annotations` inserts one and expects the
        # refusal). The guards would also have needed names past Postgres's 63-byte identifier limit:
        # `ck_served_locus_annotation__bacformer_medoid_similarity_is_not_nan` is 66.
    )

    served_locus_annotation_id: Mapped[int] = mapped_column(BigInteger, primary_key=True, autoincrement=True)
    #: The RECIPIENT — the locus the prediction is about. Dropping its catalogue drops its predictions.
    locus_id: Mapped[int] = mapped_column(ForeignKey("locus.locus_id", ondelete="CASCADE"), nullable=False)
    annotation_kind: Mapped[ServedAnnotationKind] = mapped_column(nullable=False)
    #: 0/1/2 for GO (the `locus_annotation_entry` contract), NULL for COG categories.
    gene_ontology_namespace: Mapped[int | None] = mapped_column(SmallInteger, nullable=True)
    #: ⭐ The predicted class SET, most specific first. COG: category letters, alphabetical (peers — a
    #: multi-letter category states each). GO: the donor's slim class, then its slim ancestors.
    #: ⚠ Names are joined at SERVE time (GO through the reference fold, COG letters on the page), exactly
    #: as the observed GO and EC names are, so refreshing a reference never needs a catalogue re-ingest.
    term_values: Mapped[list[str]] = mapped_column(ARRAY(Text), nullable=False)
    #: ⭐ `p_serve` — the calibrated probability that this set intersects the locus's true set. **The page
    #: shows this number, per prediction** (David, 2026-10-28: no bands). It is a promise, licensed by the
    #: held-out delivery, never measured on the locus it is shown for.
    estimated_accuracy: Mapped[float] = mapped_column(Float, nullable=False)
    #: ⭐ The donor — a real locus whose medoid gene can be inspected, which is the reason the medoid design
    #: was chosen over an average. ⚠ It may sit in ANOTHER catalogue (E. coli → Kp is most of the probe's Kp
    #: predictions).
    #: ⛔⛔ **DEFERRABLE INITIALLY DEFERRED, checked at COMMIT.** Postgres fires referential-integrity
    #: triggers per deleted ROW, so a plain NO ACTION (or RESTRICT) donor key failed a catalogue reload
    #: whenever a donor locus was deleted before its recipient — before the recipient's CASCADE had removed
    #: the row (review of 9cfa481; the end-of-statement claim this comment used to make was wrong). Deferred,
    #: a reload holding both ends (Kp → Kp) passes, because by commit the cascade has removed every row.
    #: ⚠ A row whose donor is in ANOTHER catalogue still fails that catalogue's reload at commit — so the
    #: pangenome stage refuses BEFORE the frame build (`refuse_reload_while_donating`), naming the catalogue
    #: holding those rows and the command that clears them. No CASCADE: that would delete another
    #: catalogue's predictions silently.
    donor_locus_id: Mapped[int] = mapped_column(
        ForeignKey("locus.locus_id", deferrable=True, initially="DEFERRED"), nullable=False
    )
    #: The evidence behind the number: the cosine between the two loci's medoid genes in each space.
    esm_medoid_similarity: Mapped[float] = mapped_column(Float, nullable=False)
    bacformer_medoid_similarity: Mapped[float] = mapped_column(Float, nullable=False)
    #: `at-v0.1`. ⛔ Never a bare `v0.1` — the CHECK refuses it.
    annotation_model_version: Mapped[str] = mapped_column(String(32), nullable=False)

    def __repr__(self) -> str:
        return (
            f"<ServedLocusAnnotation locus={self.locus_id} {self.annotation_kind} "
            f"ns={self.gene_ontology_namespace} {self.term_values} p={self.estimated_accuracy}>"
        )
