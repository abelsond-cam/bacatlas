"""The genomes a species' anchor control can offer — `meta.genomes` + `anchorSearch` + `GENOME_N`.

⭐ **An 80,000-element resident array becomes a query.** The published page shipped every accession
in `meta.genomes`, filtered it in the browser, and counted each genome's loci on first opening of the
box by walking its whole membership index. Here the filter is one statement and the counts are READ,
from `pangenome_genome_locus_count`, which ingest wrote — never aggregated per request, because the
aggregate is over ~412 M membership rows at the design target and this runs on every keystroke.

⛔ **The query means what `anchorSearch` meant** (`app.js::anchorSearch`): the query is trimmed and
upper-cased, a genome matches when its upper-cased accession CONTAINS it as a literal substring, and
an empty query matches every genome. Two details carry that exactly:

* **`strpos`, not `LIKE`.** `%` and `_` are wildcards to `LIKE` and plain characters to `indexOf`, so
  a `LIKE` needs escaping to mean the same thing (`locus_search_service` does it); `strpos` is a
  literal substring test with nothing to escape.
* **The query is upper-cased in Python, not in SQL.** Python's `str.upper` and JavaScript's
  `toUpperCase` both apply the full Unicode case mapping (`ß` → `SS`); a server-side `upper()` under a
  libc collation maps character by character and would not. Accessions are ASCII, where Postgres'
  `upper()` agrees with both — so that side stays in SQL.

⚠ **Ordered by `collection_genome_ordinal`**, the order `meta.genomes` published and the dropdown
listed. It is stored, never re-derived by sorting accessions.

⭐ **Genomes now carry a SOURCE GROUP, and it is DERIVED, not stored** — see `GenomeSourceGroup`.
The anchor control offers two lists because the cohort has two populations with different
properties: ~4,000 complete long-read NCBI genomes (short enough to browse) and ~75,000
AllTheBacteria short-read assemblies (which need a type-to-narrow search). Deriving the group from
`sample_id_kind` + `is_complete_assembly` rather than adding a column means no migration and no
re-ingest: both columns are already populated by `build_genome_set`.
"""

from __future__ import annotations

import enum
from dataclasses import dataclass, field

from sqlalchemy import case, func, over, select
from sqlalchemy.orm import Session

from bacatlas_backend.models.enumerations import SampleIdentifierKind
from bacatlas_backend.models.genome import Genome
from bacatlas_backend.models.genome_collection import GenomeCollectionMembership
from bacatlas_backend.models.pangenome import Pangenome
from bacatlas_backend.models.pangenome_genome_locus_count import PangenomeGenomeLocusCount


class GenomeSourceGroup(enum.Enum):
    """Which population a genome belongs to — **derived at query time, never a column.**

    ⛔ **Two axes, not one, and collapsing them loses a real distinction.** `sample_id_kind` says
    how the genome is keyed, which is also what it came from: a BioSample-keyed assembly is an
    AllTheBacteria one (`NEXT_anchor_to_genome.md`: all 281 probe manifest rows read `source = atb`
    and `ncbi_accession` is blank on every one), while an assembly stem (`GCF_…`/`GCA_…`) is from
    NCBI. `is_complete_assembly` then says whether an NCBI genome is closed. A `GCA_` draft
    long-read assembly is NCBI-keyed and NOT complete — `metadata_v2` holds 5,519 LRA rows against
    **4,017** `is_complete` — so it belongs in neither of the two lists the reader was promised, and
    `NCBI_DRAFT` exists so that it is *shown in a third group* rather than silently disappearing
    from a control that claims to list every genome.

    ⚠ `is_complete_assembly` is NULLABLE and NULL means *not applicable / not recorded*, never
    false. On the 280-genome probe cohort it is NULL on every row, so every probe genome is
    `ATB` — which is correct, and is why the NCBI group renders as a stated absence rather than as
    an empty box.
    """

    #: NCBI-keyed and `is_complete_assembly` — the "complete long-read genomes" list.
    NCBI_COMPLETE = "ncbi_complete"
    #: NCBI-keyed and not complete. ⚠ Not promised to the reader, but never hidden.
    NCBI_DRAFT = "ncbi_draft"
    #: BioSample-keyed: an AllTheBacteria short-read assembly.
    ATB = "atb"


#: ⛔ The derivation, as ONE SQL expression, so the group a row reports and the group it is counted
#: in cannot diverge — a Python-side second copy is exactly how a filter and its tally disagree.
SOURCE_GROUP_EXPRESSION = case(
    (Genome.sample_id_kind == SampleIdentifierKind.BIOSAMPLE, GenomeSourceGroup.ATB.value),
    (Genome.is_complete_assembly.is_(True), GenomeSourceGroup.NCBI_COMPLETE.value),
    else_=GenomeSourceGroup.NCBI_DRAFT.value,
)

#: What a request gets without asking. The page listed all 100; the design target is 80,000.
DEFAULT_GENOME_LIMIT = 100

#: The clamp. ⚠ A LIMIT, not a page: `truncated` says when there is more, and the reader narrows the
#: query rather than scrolling 80,000 accessions.
MAXIMUM_GENOME_LIMIT = 1000


@dataclass
class ListedGenome:
    """One genome the anchor control can offer, with both of its locus counts."""

    sample_id: str
    collection_genome_ordinal: int
    #: Distinct loci where the genome has ≥ 1 gene.
    locus_count: int
    #: Distinct loci where it sits in some arrangement — the published dropdown's "N loci".
    arrangement_locus_count: int
    #: Which of the anchor control's lists this genome belongs in. Derived — see `GenomeSourceGroup`.
    source_group: GenomeSourceGroup


@dataclass
class GenomeListing:
    """The matching genomes, and enough counts that a short list is never mistaken for a whole one."""

    query: str
    genome_count: int
    #: ⛔ Matches BEFORE the limit. Without it a list cut at 100 reads as "100 genomes match".
    matched_genome_count: int
    genomes: list[ListedGenome] = field(default_factory=list)
    #: ⭐ **How many genomes this CATALOGUE holds in each source group, ignoring `q` and the
    #: limit.** ⛔ Independent of the query on purpose: the control has to distinguish *"no complete
    #: NCBI genome matches what you typed"* from *"this catalogue contains no complete NCBI genomes
    #: at all"*, and a count that moved with `q` could say only the first. The second is the true
    #: answer on every catalogue loaded today, and a reader told the first would keep typing.
    group_counts: dict[GenomeSourceGroup, int] = field(default_factory=dict)
    #: Which group was asked for, or `None` for all of them.
    source_group: GenomeSourceGroup | None = None

    @property
    def truncated(self) -> bool:
        """More genomes matched than the limit let through."""
        return self.matched_genome_count > len(self.genomes)


def clamp_genome_limit(requested: int | None) -> int:
    """`limit` as asked, held to [1, `MAXIMUM_GENOME_LIMIT`]; absent means the default."""
    if requested is None:
        return DEFAULT_GENOME_LIMIT
    return max(1, min(MAXIMUM_GENOME_LIMIT, requested))


def count_genome_source_groups(
    session: Session, pangenome: Pangenome
) -> dict[GenomeSourceGroup, int]:
    """How many genomes of this catalogue sit in each source group — **one grouped statement.**

    ⛔ **This is CATALOGUE metadata and belongs on the shell, not on the picker's keystroke.** The
    picker is budgeted at one statement plus the species resolution and is asserted to contain no
    `GROUP BY` at all (`test_route_cost.py`), because at the design target it fires on every
    keypress. The shell is the one route permitted to aggregate over the catalogue, it is fetched
    once, and its answer is cacheable — so the counts ride there.

    ⚠ **Independent of any query, deliberately.** The control has to distinguish *"nothing matches
    what you typed"* from *"this catalogue holds none of these at all"*. Only the second is true of
    the complete-genome list on every catalogue loaded today, and a reader told the first would keep
    typing. A count that moved with `q` could not say it.

    ⭐ Groups with zero members are PRESENT with a zero, not absent: an absent key cannot be told
    from one the server forgot to compute.

    ⛔⛔ **It counts over the SAME population `list_genomes` can serve, and the join on
    `pangenome_genome_locus_count` is what makes that true.** Counting collection membership instead
    is the obvious version and it is wrong: the listing joins that table for its two per-genome
    counts, so a collection member without a count row is INVISIBLE to the list while still being
    counted by a membership-only census. The control would then say *"1 complete genome"* over a
    list that could never show it, and the reader would keep looking for it. In a published
    catalogue the publish gate forbids that state
    (`test_the_publish_gate_REFUSES_a_collection_genome_with_no_count_row`), so the two scopes agree
    there anyway — but agreeing *by construction* costs one join, and agreeing *by invariant* costs
    a bug the first time the invariant is relaxed. Caught by
    `test_the_row_and_the_FILTER_cannot_disagree`, which built exactly that state.
    """
    rows = session.execute(
        select(SOURCE_GROUP_EXPRESSION.label("source_group"), func.count())
        .select_from(PangenomeGenomeLocusCount)
        .join(Genome, Genome.genome_id == PangenomeGenomeLocusCount.genome_id)
        .join(
            GenomeCollectionMembership,
            (GenomeCollectionMembership.genome_id == PangenomeGenomeLocusCount.genome_id)
            & (GenomeCollectionMembership.genome_collection_id == pangenome.genome_collection_id),
        )
        .where(PangenomeGenomeLocusCount.pangenome_id == pangenome.pangenome_id)
        .group_by(SOURCE_GROUP_EXPRESSION)
    ).all()
    counts = {group: 0 for group in GenomeSourceGroup}
    for value, count in rows:
        counts[GenomeSourceGroup(value)] = count
    return counts


def list_genomes(
    session: Session,
    pangenome: Pangenome,
    *,
    query: str,
    limit: int,
    source_group: GenomeSourceGroup | None = None,
    group_counts: dict[GenomeSourceGroup, int] | None = None,
) -> GenomeListing:
    """Every genome of the pangenome's collection whose accession contains `query` — ONE statement.

    ⭐ `count(*) OVER ()` is evaluated before `LIMIT`, so the number of matches rides on every row the
    limit keeps and the statement count does not grow with the limit, the collection or the query.
    Zero rows back means zero matched.

    ⭐ `source_group` narrows to one of the anchor control's lists. It is a `WHERE` on the same
    expression the `source_group` column reports, so a row's reported group and the filter that
    admitted it cannot disagree — a Python-side second copy of the rule is exactly how a filter and
    its tally drift apart.

    ⚠ `group_counts` is passed IN rather than computed here, because computing it needs a `GROUP BY`
    this route is asserted not to contain. `count_genome_source_groups` is the shell's job.
    """
    trimmed = (query or "").strip()
    matched = over(func.count())
    statement = (
        select(
            Genome.sample_id,
            GenomeCollectionMembership.collection_genome_ordinal,
            PangenomeGenomeLocusCount.locus_count,
            PangenomeGenomeLocusCount.arrangement_locus_count,
            matched,
            SOURCE_GROUP_EXPRESSION.label("source_group"),
        )
        .join(Genome, Genome.genome_id == PangenomeGenomeLocusCount.genome_id)
        .join(
            GenomeCollectionMembership,
            (GenomeCollectionMembership.genome_id == PangenomeGenomeLocusCount.genome_id)
            & (GenomeCollectionMembership.genome_collection_id == pangenome.genome_collection_id),
        )
        .where(PangenomeGenomeLocusCount.pangenome_id == pangenome.pangenome_id)
        .order_by(GenomeCollectionMembership.collection_genome_ordinal)
        .limit(limit)
    )
    if trimmed:
        # ⛔ Literal substring, case-folded as `anchorSearch` folds it — see the module docstring.
        statement = statement.where(func.strpos(func.upper(Genome.sample_id), trimmed.upper()) > 0)
    if source_group is not None:
        statement = statement.where(SOURCE_GROUP_EXPRESSION == source_group.value)

    rows = session.execute(statement).all()
    return GenomeListing(
        query=trimmed,
        genome_count=pangenome.genome_count,
        matched_genome_count=rows[0][4] if rows else 0,
        source_group=source_group,
        group_counts=group_counts or {},
        genomes=[
            ListedGenome(
                sample_id=sample_id,
                collection_genome_ordinal=ordinal,
                locus_count=locus_count,
                arrangement_locus_count=arrangement_locus_count,
                source_group=GenomeSourceGroup(group_value),
            )
            for sample_id, ordinal, locus_count, arrangement_locus_count, _, group_value in rows
        ],
    )
