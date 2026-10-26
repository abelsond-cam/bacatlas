"""The locus row — a syntelog locus in one pangenome, with everything the card shows.

Wide on purpose: these are scalars the page filters, sorts and prints, and a locus is fetched whole.
"""

from __future__ import annotations

from sqlalchemy import (
    BigInteger,
    ForeignKey,
    Index,
    Integer,
    String,
    Text,
    UniqueConstraint,
)
from sqlalchemy.dialects.postgresql import ARRAY
from sqlalchemy.orm import Mapped, mapped_column

from bacatlas_backend.database import Base
from bacatlas_backend.models.column_types import measurement, nan_guards
from bacatlas_backend.models.enumerations import (
    GeneOntologyAgreementVerdict,
    PrevalenceBand,
    PropagationVerdict,
)


class Locus(Base):
    """One node of one pangenome."""

    __tablename__ = "locus"
    __table_args__ = (
        UniqueConstraint("pangenome_id", "node_label"),
        UniqueConstraint("pangenome_id", "catalogue_ordinal"),
        Index("ix_locus__pangenome_id__prevalence_band", "pangenome_id", "prevalence_band"),
        Index("ix_locus__pangenome_id__member_genome_count", "pangenome_id", "member_genome_count"),
        Index("ix_locus__pangenome_id__interest_score", "pangenome_id", "interest_score"),
        # ⚠ The FK to `genome` needs this to delete a genome without scanning every locus. Rare, but
        # at 889,160 loci a genome retraction should not be an outage — see the note on
        # `gene_locus_membership`, where the same omission cost an 8-minute stall.
        Index("ix_locus__medoid_genome_id", "medoid_genome_id"),
        # ⭐⭐ **THE TRIGRAM INDEX IS WHAT MAKES SEARCH THE SAME QUESTION IT IS TODAY.**
        # `serving_at_scale.md` §6 could not choose between static prefix buckets — which LOSE
        # mid-word substring, so `ligase` stops finding *O-antigen ligase RfaL* — and trigram
        # postings. `gin_trgm_ops` accelerates `LIKE '%q%'` directly, so `app.js`'s
        # `HAY[i].indexOf(q)` is preserved rather than approximated. Without this index the same
        # query still returns the same rows, by sequential scan: correct, and unusable at 889k loci.
        Index(
            "ix_locus__search_text_trigram",
            "search_text",
            postgresql_using="gin",
            postgresql_ops={"search_text": "gin_trgm_ops"},
        ),
        *nan_guards(
            "syntenic_a5",
            "uniref50_impurity",
            "uniref50_coverage",
            "resolved_threshold",
            "seqid_coverage",
            "interest_score",
        ),
    )

    locus_id: Mapped[int] = mapped_column(BigInteger, primary_key=True, autoincrement=True)
    pangenome_id: Mapped[int] = mapped_column(
        ForeignKey("pangenome.pangenome_id", ondelete="CASCADE"), nullable=False
    )
    pathogen_species_id: Mapped[int] = mapped_column(
        ForeignKey("pathogen_species.pathogen_species_id"), nullable=False
    )

    #: ⛔ **TEXT, always.** Node labels look numeric and are not; read back as int64 they silently
    #: stop matching. This is also the URL hash and the durable identity across a re-export.
    node_label: Mapped[str] = mapped_column(String(64), nullable=False)

    #: `export_payload.node_order`'s index — the dense 0..n-1 position `arr.vec`, `ctx.nid`,
    #: `gaps.a/b`, `map_reps.near` and `.loci` all point at. ⚠ **Payload-local**: a re-export
    #: renumbers every one of them, which is why `node_label` and not this is the durable key.
    catalogue_ordinal: Mapped[int] = mapped_column(Integer, nullable=False)

    member_gene_count: Mapped[int] = mapped_column(Integer, nullable=False)
    member_genome_count: Mapped[int] = mapped_column(Integer, nullable=False)
    prevalence_band: Mapped[PrevalenceBand] = mapped_column(nullable=False)

    # ── naming ────────────────────────────────────────────────────────────────────────────────
    #: ⭐ **Materialised at ingest, NOT NULL.** `displayOf` falls through `inferredName` →
    #: `soleArch` → the locus's full Pfam list + the Pfam name table, then `modalProduct` → its full
    #: product list. Computing it per neighbour at request time is the fan-out problem; computing it
    #: once here turns it into a five-column lookup over ~20 ids.
    display_name: Mapped[str] = mapped_column(String(256), nullable=False)
    #: `bakta_symbol` | `pfam_architecture` | `product` | `label` — so the page can say WHY.
    display_name_source: Mapped[str] = mapped_column(String(32), nullable=False)
    display_name_source_accession: Mapped[str | None] = mapped_column(String(64), nullable=True)
    #: ⭐ Also materialised. `bestProduct` reads every product row with counts and a regex built
    #: from the modal symbol; it is not `products[0]`.
    best_product: Mapped[str | None] = mapped_column(Text, nullable=True)
    #: The modal REAL Bakta symbol, locus tags stripped. NULL = unnamed (payload `-1`).
    bakta_gene_symbol: Mapped[str | None] = mapped_column(String(128), nullable=True)
    named_member_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)

    # ── UniRef50 ──────────────────────────────────────────────────────────────────────────────
    uniref50_family_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    #: Families holding ≥10 % of labelled members — **the number the card claims with**, while the
    #: list below shows every family.
    uniref50_major_family_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    uniref50_labelled_member_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    #: ⚠ NULL below 5 labelled members — not measured, never 0.
    uniref50_impurity: Mapped[float | None] = measurement()
    uniref50_coverage: Mapped[float | None] = measurement()

    # ── Pfam ──────────────────────────────────────────────────────────────────────────────────
    #: ⛔ Members carrying ANY Pfam-A domain. `member_gene_count − this` is **missing coverage, not
    #: a competing architecture** — the architecture list drops nulls, so without this denominator
    #: the page cannot tell the two apart.
    pfam_annotated_member_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    #: The TRUE architecture count, which the capped top-N list cannot show.
    pfam_architecture_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    #: ⛔ The audit's own `class_clan` verdict, **read and never re-derived**. The page once
    #: re-derived it by a different rule and disagreed on 2 of 22,624 loci — *"a page that quotes
    #: the report must not be able to contradict it, ever."*
    pfam_concordance_class: Mapped[str | None] = mapped_column(String(32), nullable=True)

    # ── evidence ──────────────────────────────────────────────────────────────────────────────
    # ⛔⛔ **THE THREE COLUMNS BELOW ARE COMPUTED OVER AT MOST 50 MEMBERS, NOT OVER THE LOCUS.**
    # `accessory_audit.FAMILY_MEMBER_CAP = 50`, longest-first: the within-cluster all-vs-all is
    # quadratic in members, so each cluster is capped to bound the cost. Measured on the published
    # ecoli catalogue, **4,138 of the 12,104 loci that carry a tier sit at that cap** (34.2 %) —
    # their real member counts run to 143 — so for a third of the catalogue `resolved_threshold`
    # is *"the identity at which the 50 longest members first group"*, which is a different and
    # weaker claim than *"at which the members first group"*.
    # ⚠ `member_gene_count` above is NOT affected: it comes from the assignment via
    # `node_membership_profile`, never from the waterfall. Do not fill it from the audit.
    syntenic_a5: Mapped[float | None] = measurement()
    #: `mmseq@0.4`, `pfam_not_alignable`, `esm_homology`, `synteny_only`, … ⚠ `synteny only` names
    #: the EVIDENCE, not a mistake: non-homologous cargo in a genuinely conserved neighbourhood is
    #: what an identity threshold cannot see, and the reason the method exists.
    collapse_tier: Mapped[str | None] = mapped_column(String(64), nullable=True)
    collapse_bucket: Mapped[str | None] = mapped_column(String(64), nullable=True)
    #: ⚠ Over ≤50 longest members — see the block above.
    resolved_threshold: Mapped[float | None] = measurement()
    seqid_coverage: Mapped[float | None] = measurement()

    # ── the representative member ─────────────────────────────────────────────────────────────
    #: ⛔ The four `*_medoid_distance` columns, `embed_within_over_nearest` and the two
    #: `separation_percentile_*` went on 2026-09-24. Each measured a locus by reducing it to the
    #: member below — its MEDOID — and then measuring that one point; the replacement measures the
    #: whole set and lives in `locus_cluster_similarity`, at a grain (locus, representation) that
    #: four columns per quantity could not express.
    #:
    #: These two stay because they are not a geometry: they NAME the representative member, which is
    #: an ordinary fact about the locus and survives the construction that used it. ⚠ Nothing reads
    #: them today; they are kept for an exemplar view, not because anything is waiting for them.
    medoid_genome_id: Mapped[int | None] = mapped_column(ForeignKey("genome.genome_id"), nullable=True)
    medoid_flat_index: Mapped[int | None] = mapped_column(Integer, nullable=True)

    # ── function block ────────────────────────────────────────────────────────────────────────
    cog_annotated_member_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    #: ⛔ A COUNT, never a relation. COG ids are single-valued per gene, so every member set is a
    #: singleton and an agreement ladder could only ever say `single` or `disjoint`. Two COG ids in
    #: one locus is an ordinary consequence of grouping above family level.
    cog_distinct_id_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    #: ⛔⛔ **A SET OF SINGLE LETTERS, NOT ONE CATEGORY — the same defect as the gene layer's.**
    #: Bakta emits `cog_cat` as a CONCATENATED set (`EHJQ`), and this is the modal such set over the
    #: locus's annotated members. Measured on the published catalogues: **95 of 118 distinct values
    #: in *E. coli* and 93 of 117 in kp carry more than one letter**, to a maximum of four. Stored as
    #: `String(16)` it FITS, so nothing fails — and `WHERE modal_cog_category = 'C'` silently misses
    #: every multi-category locus, which is precisely how the gene-layer version of this went unseen.
    #: ⚠ The page prints them joined; the array is what makes `@> ARRAY['C']` mean what it says.
    modal_cog_categories: Mapped[list[str] | None] = mapped_column(ARRAY(String(1)), nullable=True)
    ec_annotated_member_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    kegg_annotated_member_count: Mapped[int | None] = mapped_column(Integer, nullable=True)

    #: Coverage BEFORE every verdict, against the locus size — a share against the annotated subset
    #: would read 100 % on a locus where one gene in forty carries a label.
    go_annotated_member_count_molecular_function: Mapped[int | None] = mapped_column(Integer, nullable=True)
    go_annotated_member_count_biological_process: Mapped[int | None] = mapped_column(Integer, nullable=True)
    go_annotated_member_count_cellular_component: Mapped[int | None] = mapped_column(Integer, nullable=True)
    go_verdict_molecular_function: Mapped[GeneOntologyAgreementVerdict | None] = mapped_column(nullable=True)
    go_verdict_biological_process: Mapped[GeneOntologyAgreementVerdict | None] = mapped_column(nullable=True)
    go_verdict_cellular_component: Mapped[GeneOntologyAgreementVerdict | None] = mapped_column(nullable=True)

    #: ⛔⛔ **Genes carrying ANY GO term — and it is NOT the sum of the three columns above.**
    #: `instruments/annotation_transfer.py:231-234` says those three cannot be summed; this is the
    #: measurement behind that sentence, read off `bacatlas_ingest` on 2026-10-26: **351,559 ecoli
    #: genes carry a GO term and the three namespace columns sum to 808,043** — 2.30× too many, an
    #: overstatement of 456,484 genes. kp: **220,638 against 481,685**, 2.18×. A gene annotated in
    #: two namespaces is counted twice and in all three, three times, and nothing about the inflated
    #: number looks wrong.
    #:
    #: ⭐ This column is what lets `function_inference_service._annotated_counts` stop reading gene
    #: rows. The other three vocabularies' counts were already exact — re-verified here against the
    #: gene-grain recount on **all four** catalogues, 0 mismatching loci of 17,531 / 15,670 / 18,300
    #: / 16,314 — so GO was the only reason that per-request join existed.
    go_annotated_member_count_any: Mapped[int | None] = mapped_column(Integer, nullable=True)

    #: ⛔ **Genes with a non-empty folded CLAIM in this namespace — a different population from
    #: `go_annotated_member_count_*` one block up, and always the SMALLER one.** That column counts a
    #: gene carrying a *term* in the namespace; this one counts a gene whose terms survive the slim
    #: fold and therefore state something comparable. A gene asserting only the namespace root — 452
    #: *E. coli* loci state nothing but `cellular_component` — is counted there and says nothing here.
    #:
    #: Measured on `bacatlas_ingest`, 2026-10-26, stored-annotated → claiming:
    #: ecoli **287,999 → 279,268** molecular function (−8,731), 276,052 → 268,562 biological process
    #: (−7,490), 243,992 → 238,278 cellular component (−5,714); kp −6,896 / −3,394 / −2,162.
    #: ⚠ `transfer/claims.py:102-107` says reading the annotated columns instead *"understated one
    #: baseline by 25,751 genes"*. That figure is from a different defect — summing `genes_carrying`
    #: over labelled nodes only (`scripts/measure_annotation_propagation_loop.py`) — and its sign is
    #: the opposite of the mechanism named in the same sentence: a root-only gene makes the annotated
    #: column the LARGER number. The three gaps above are this column's own, measured here.
    #:
    #: ⚠ Folded, so **meaningful only against `pangenome.annotation_fold_version`** — see there.
    go_claiming_member_count_molecular_function: Mapped[int | None] = mapped_column(Integer, nullable=True)
    go_claiming_member_count_biological_process: Mapped[int | None] = mapped_column(Integer, nullable=True)
    go_claiming_member_count_cellular_component: Mapped[int | None] = mapped_column(Integer, nullable=True)

    # ── propagation: does this node's own call agree with itself? ──────────────────────────────
    #: ⭐⭐ **`within_node_propagation`'s five-state decision, stored per locus per vocabulary.**
    #: One column each for COG, EC, KEGG and the GO slim — **four, not six**: the decision ANDs over
    #: GO's three facets (`unanimous_at_deepest_shared_rung`, the `verdicts.append` loop), so a node
    #: has ONE GO verdict and not one per namespace, and `compute_propagation` is keyed by
    #: `annotation_kind` alone. Three facet columns could not be aggregated back into the rate the
    #: page quotes without reassembling that AND, and nothing reads a per-facet propagation rate.
    #:
    #: ⛔ **NULL and `no_coverage` are different facts and the distinction is load-bearing.**
    #: `no_coverage` means *the ingest looked and no member gene carries this vocabulary*;
    #: NULL means *this ingest generation never computed it*. The catalogues loaded before this
    #: column existed hold NULL, and a reader treating NULL as `no_coverage` would conclude that
    #: 17,531 ecoli loci carry no COG when 4,993 of them do. This is the same shape of error as the
    #: six all-NULL columns recorded at `ingest/catalogue_frames.py:441-448`, where NULL was
    #: indistinguishable from a documented *"not measured"*.
    #:
    #: ⚠ The GO verdict is folded, so it is meaningful only against
    #: `pangenome.annotation_fold_version`. The other three use no fold.
    #:
    #: ⚠ **No page component reads these four, and that is recorded rather than implied.** What the
    #: Function tab reads per locus is the COUNTS — `frontend/src/lib/functionVocabulary.ts:205-213`
    #: derives its own "carried by N of M genes" from `annotated_gene_count` and `member_gene_count`,
    #: and even the `checkable` boolean the API already ships is read by nothing. These columns are
    #: here for two other reasons: they make `pangenome_propagation_rate` a `GROUP BY` that can be
    #: re-derived and re-checked in SQL after the gene rows are gone, and `discordant` is a
    #: queryable statement about the model — **22 of 17,531 ecoli loci and 7 of 15,670 kp loci** have
    #: COG calls that disagree with each other — which today costs a 1.3 M-row scan to find.
    #:
    #: ⭐ `not_comparable` is the state that was previously invisible: on the FOLDED path production
    #: runs, **81 ecoli and 33 kp loci** carry GO on two or more genes with no facet every one of
    #: them states. On the unfolded path it is 0 — which is why
    #: `tests/test_within_node_propagation.py`'s coverage identity passes while omitting the fold.
    propagation_verdict_cog: Mapped[PropagationVerdict | None] = mapped_column(nullable=True)
    propagation_verdict_ec: Mapped[PropagationVerdict | None] = mapped_column(nullable=True)
    propagation_verdict_kegg: Mapped[PropagationVerdict | None] = mapped_column(nullable=True)
    propagation_verdict_go: Mapped[PropagationVerdict | None] = mapped_column(nullable=True)

    # ── extent ────────────────────────────────────────────────────────────────────────────────
    #: ⚠ 1-based inclusive, hence `end - start + 1`. The block on the track is drawn
    #: `max(6, len_nt / 10)` px wide, and a NEIGHBOUR's value sets a neighbour block's width.
    median_gene_length_nt: Mapped[int | None] = mapped_column(Integer, nullable=True)
    gene_length_interquartile_range_nt: Mapped[int | None] = mapped_column(Integer, nullable=True)

    # ── the marginal denominator, kept as an array ────────────────────────────────────────────
    #: `ctx.obs` — member genes for which each of the ten offsets EXISTS, in `OFFSETS` order
    #: `[-5..-1, 1..5]`. Counted **before** the top-N cut, so "and N others" is honest and
    #: `obs < size` exposes contig-edge truncation. An array because a table whose only non-key
    #: column is one integer is a 40× overhead.
    context_observed_member_counts: Mapped[list[int]] = mapped_column(ARRAY(Integer), nullable=False)

    #: `arr.tot` — arrangements in TOTAL. ⛔ Never moved by any display cap, and never conflated
    #: with the number listed: the page keys sentences on this.
    total_arrangement_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)

    #: Member genes that sit in SOME arrangement — summed over ALL of them, capped or not.
    #:
    #: ⛔ **This exists so the API can keep two remainders apart that the published page never had
    #: to.** That payload is uncapped, so there `member_gene_count − Σ(listed)` really did mean
    #: *"member genes with no recorded neighbourhood — no coordinates for the gene, so no window"*.
    #: The API caps at 8, so the same subtraction silently swept in every member sitting in an
    #: arrangement past the cap: measured at **15,912 E. coli genes over 2,340 loci and 10,437 kp
    #: genes over 1,548**, each of which the page would have told a reader has no coordinates.
    #: Stored rather than aggregated per request so the locus view stays at one statement per table.
    arrangement_member_gene_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)

    #: DISTINCT genomes with at least one gene in SOME arrangement — the UNION, over all of them.
    #:
    #: ⛔ **Not derivable from `arrangement_member_gene_count`, and not a sum.** The page has two
    #: different sentences for an anchored genome that carries no arrangement here — *"has no gene at
    #: this locus"* and *"has no recorded neighbourhood at this locus"* — and only one of them is
    #: ever true (`app.js:1762-1771`). Which one is settled by whether every genome present reaches
    #: an arrangement, which is a GENOME question: a genome at ρ > 1 has two genes here, and one of
    #: them missing a window leaves it counted in the gene remainder while it is still perfectly
    #: present in an arrangement. `sum(member_genome_count)` is likewise wrong, and wrong upward, for
    #: the same reason — it double-counts exactly those genomes.
    #:
    #: Measured on the published payloads: **6.26 % of ecoli loci and 3.69 % of kp loci** have at
    #: least one genome present and in no arrangement (worst case 64 of 100), and at those loci
    #: *"this genome has no gene here"* is simply false.
    arrangement_member_genome_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)

    #: `render_page.ranking`'s score, precomputed. Drives the landing locus and the example chips,
    #: which were a render-time payload mutation and are now columns behind a verification gate.
    interest_score: Mapped[float | None] = measurement()

    #: The materialised search haystack: symbol + label + every product + every UniRef50 accession,
    #: lowercased. ⚠ Pfam/COG/GO/EC/KEGG are deliberately NOT in it — that is what `app.js` searches
    #: today, and adding them is a product change, not an implementation detail.
    search_text: Mapped[str | None] = mapped_column(Text, nullable=True)

    def __repr__(self) -> str:
        return f"<Locus {self.node_label} {self.display_name!r} n={self.member_gene_count}>"
