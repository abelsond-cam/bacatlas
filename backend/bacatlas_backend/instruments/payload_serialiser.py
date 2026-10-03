"""Rebuild the published locus-browser payload from the database.

⭐⭐ **This is an instrument, not an endpoint.** The API serves its own, better shape — nulls with a
named ``absence_reason``, a resolved 6×6 cosine matrix, four separately named remainders. This
module exists to prove a different claim: that the database is a **lossless superset** of the JSON
catalogue the static site ships, checked against the real 17,531-locus payload rather than asserted.
Once that holds, the static site can be retired knowing nothing on it was lost.

⛔⛔ **EMISSION ORDER IS THE WHOLE PROBLEM.** ``_Intern`` assigns a string its index *on first use*,
so the pools in ``strings`` are in whatever order ``build_payload`` walked the data. Emit the same
strings in a different order and the payload is **correct and not identical** — every ``idx`` column
differs and a column-wise diff reports the catalogue changed. The order is therefore not a detail of
this file, it is its contract:

1. ``nodes`` first, **column by column in payload order**, because its dict literal is what interns
   ``name`` → ``tier`` → ``pfclass`` → ``cog_cat`` → ``go_0/1/2``;
2. then the ``lists`` loop, in ``build_payload``'s own order: sym, prod, u50, pfam, cog, ec, kegg, go;
3. then the ``u50`` **second pass**, which adds to ``prod``, ``pfam`` and ``sym`` — three pools that
   already have entries.

⭐ **The order is checked, not trusted.** :func:`build_payload_from_database` runs
:func:`verify_intern_walk` over its own output before returning, so a wrong emission order raises
here — with the pool and the column named — instead of surfacing later as thousands of changed
indices in a diff nobody can read.

⚠ **Every value is READ, never recomputed.** The ingest wrote the float columns through the
payload's own ``_floats``/``_sigfigs`` helpers, so ``locus.syntenic_a5`` already *is* ``nodes.a5``.
Re-rounding a stored value would risk double rounding; re-deriving one would risk disagreeing with
the report, which is the mistake ``pfam_concordance_class`` exists to prevent.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from sqlalchemy import select
from sqlalchemy.orm import Session

from bacatlas_backend.instruments.payload_string_pool import StringPool
from bacatlas_backend.models.locus import Locus

#: The payload schema this database was INGESTED from, and therefore the one a rebuild reproduces.
#:
#: ⛔ Pinned, not imported. It used to read ``nuna.SCHEMA_VERSION`` live, which was right only while the
#: two happened to agree: nuna bumped 14 → 16 on 2026-09-24 and this rebuild began claiming to be a
#: schema-16 payload while carrying schema-14 rows, with the byte-identity test as the sole signal. A
#: rebuild reproduces a SPECIFIC published artifact, so the number belongs to the artifact. It moves in
#: the same commit as the re-ingest that changes what is in these tables — never before it.
#:
#: ⛔ **14 → 16 on 2026-10-03, and it was ALREADY WRONG when this was written.** The rule above is
#: right and was broken in the other direction: migration `858439e8fd7c` (2026-09-24) dropped the four
#: `*_d_*` medoid columns and the map/null tables and added `locus_similarity`, and the re-ingest that
#: filled them ran on **2026-09-25** (`pangenome.ingested_at` 08:35 ecoli / 08:42 kp). The pin did not
#: move with it. So from 2026-09-25 this rebuild emitted a **schema-16-SHAPED payload stamped
#: `"schema": 14`** — carrying `sim`, omitting `map_reps`, `null` and the four `*_d_*` columns — which
#: no schema-14 reader could consume.
#:
#: ⚠ **Why nothing caught it for eight days.** `test_every_other_block_is_BYTE_IDENTICAL…` compares
#: `schema` against the oracle, and the oracle was the frozen **schema-14** page: the stale pin and the
#: stale oracle were wrong in the SAME direction, so they agreed. Two wrongs agreeing is exactly what
#: that test's own comment claims it is proof against. Nothing else looks: nuna's `diff_payloads` skips
#: the key outright (`verify_payload_invariance.py:120-121`). Re-exporting the oracle at schema 16 is
#: what separated them, and it failed immediately — `assert 16 == 14`.
INGESTED_PAYLOAD_SCHEMA = 16


def _payload_constants():
    """The constants `build_payload` itself used, imported rather than copied — code, not a memory.

    ⚠ ``SCHEMA_VERSION`` is deliberately NOT among them — see ``INGESTED_PAYLOAD_SCHEMA`` above. The rest
    are genuinely shared vocabulary (band order, offsets, the policy block) and would be a second thing to
    keep in step if copied.
    """
    from nuna.tl.locus_browser.export_payload import BAND_ORDER, NS_CODE, OFFSETS, POLICY

    return INGESTED_PAYLOAD_SCHEMA, tuple(OFFSETS), tuple(BAND_ORDER), dict(NS_CODE), POLICY


@dataclass
class CataloguePools:
    """The eleven ``strings`` pools, shared across blocks — which is why order is global, not local.

    ⚠ ``cog`` is one pool holding **two vocabularies**: ``nodes.cog_cat`` writes Bakta category
    strings into it, and ``lists.cog.idx`` then writes COG *ids*. That is not a bug to tidy; it is
    how every published catalogue was interned.
    """

    sym: StringPool = field(default_factory=StringPool)
    prod: StringPool = field(default_factory=StringPool)
    u50: StringPool = field(default_factory=StringPool)
    pfam: StringPool = field(default_factory=StringPool)
    tier: StringPool = field(default_factory=StringPool)
    pfclass: StringPool = field(default_factory=StringPool)
    cog: StringPool = field(default_factory=StringPool)
    go: StringPool = field(default_factory=StringPool)
    ec: StringPool = field(default_factory=StringPool)
    kegg: StringPool = field(default_factory=StringPool)
    goclass: StringPool = field(default_factory=StringPool)

    def tables(self) -> dict[str, list[str]]:
        """The ``strings`` block — the pools as they stand after every block has been emitted."""
        return {name: getattr(self, name).values for name in self.__dataclass_fields__}


def _built_date(built_at) -> str | None:
    """``meta.built`` — a plain ``YYYY-MM-DD``, whatever shape the column turns out to hold.

    ⛔ `Pangenome.built_at` is `Mapped[str | None]` over `String(32)`, and this call site read
    `built_at.date().isoformat()`. That is an `AttributeError` on a `str`, and it has never fired
    only because the column is NULL on every catalogue loaded so far — nothing in `ingest` writes
    it. The crash would arrive with the data, which is the worst time for it.
    """
    if not built_at:
        return None
    if hasattr(built_at, "date"):
        return built_at.date().isoformat()
    # ⚠ The leading date of an ISO string, which is what the exporter writes (`2026-10-03`).
    return str(built_at).split("T")[0].split(" ")[0]


def _enum_value(value) -> str | None:
    """An enum's payload spelling, or ``None`` — the DB stores enums, the payload stores strings."""
    return None if value is None else getattr(value, "value", value)


def load_catalogue_loci(session: Session, pangenome_id: int) -> list[Locus]:
    """Every locus of one pangenome, **in catalogue order**.

    ⛔ ``catalogue_ordinal`` is a stored ordinal, not a sort of convenience: it is the payload's own
    array index, written once at ingest from ``node_order``. Ordering by anything else — the label,
    the id, the size — produces a different catalogue that is internally consistent and wrong,
    because every ``nid``, ``near``, ``a``/``b`` and ``failures`` entry addresses this order.
    """
    loci = list(
        session.execute(
            select(Locus)
            .where(Locus.pangenome_id == pangenome_id)
            .order_by(Locus.catalogue_ordinal)
        ).scalars()
    )
    expected = list(range(len(loci)))
    actual = [locus.catalogue_ordinal for locus in loci]
    if actual != expected:
        # A hole or a duplicate here shifts every downstream index by one and nothing else notices.
        raise ValueError(
            f"pangenome {pangenome_id} has {len(loci):,} loci whose catalogue_ordinal is not "
            f"0…{len(loci) - 1:,} — the payload's array index cannot be reconstructed"
        )
    return loci


def node_block(loci: list[Locus], pools: CataloguePools) -> dict[str, list]:
    """The ``nodes`` block — **built column by column, in the payload's own column order**.

    ⛔ The order of the statements below is load-bearing. ``build_payload``'s ``node_block`` is a
    dict literal, evaluated top to bottom, and each interning column writes into a shared pool as it
    is evaluated. Reordering these lines reorders ``strings.sym``, ``strings.tier``,
    ``strings.pfclass``, ``strings.cog`` and ``strings.goclass`` — and therefore every index that
    addresses them.

    ⚠ ``len_nt`` and ``len_iqr`` are APPENDED at the end, never inserted: the invariance check reads
    ``nodes`` by position, so a column slipped into the middle makes every prior payload look
    different for no reason.
    """
    _, _, band_order, namespace_codes, _ = _payload_constants()
    block: dict[str, list] = {
        "label": [locus.node_label for locus in loci],
        "size": [locus.member_gene_count for locus in loci],
        "genomes": [locus.member_genome_count for locus in loci],
        "band": [band_order.index(_enum_value(locus.prevalence_band)) for locus in loci],
        "name": pools.sym.indices(locus.bakta_gene_symbol for locus in loci),
        "n_named": [locus.named_member_count for locus in loci],
        "n_u50": [locus.uniref50_family_count for locus in loci],
        "n_u50_major": [locus.uniref50_major_family_count for locus in loci],
        "n_u50_labelled": [locus.uniref50_labelled_member_count for locus in loci],
        "n_pfam": [locus.pfam_annotated_member_count for locus in loci],
        "a5": [locus.syntenic_a5 for locus in loci],
        "tier": pools.tier.indices(locus.collapse_tier for locus in loci),
        "pfclass": pools.pfclass.indices(locus.pfam_concordance_class for locus in loci),
        # ⚠ -1, not null: `_Intern` never saw this column — `build_payload` wrote the sentinel itself.
        "n_arch": [-1 if locus.pfam_architecture_count is None else locus.pfam_architecture_count for locus in loci],
        "resolved": [locus.resolved_threshold for locus in loci],
        # ⛔ `esm_d_intra` / `esm_d_near` / `bac_d_intra` / `bac_d_near` went on 2026-09-24 with the
        # card that read them. They were cosine DISTANCES to a locus's medoid; `sim` carries medians
        # over the whole set, in a block of its own rather than four columns here.
        "n_cog": [locus.cog_annotated_member_count for locus in loci],
        "cog_n": [locus.cog_distinct_id_count for locus in loci],
        # ⛔ Bakta writes a SET of COG categories as a letter run (`MV`, `KG`, `DN`), which the
        # schema stores as an array so it can be queried by letter. The payload carries the run, so
        # the join must preserve the stored order — the array is not a set to be sorted.
        "cog_cat": pools.cog.indices(
            None if not locus.modal_cog_categories else "".join(locus.modal_cog_categories)
            for locus in loci
        ),
        "n_ec": [locus.ec_annotated_member_count for locus in loci],
        "n_kegg": [locus.kegg_annotated_member_count for locus in loci],
    }
    # The GO verdicts and their coverage counts, per namespace — the verdict NEVER without the
    # count, because a list of labels cannot say whether members disagreed or were never annotated.
    for namespace, code in namespace_codes.items():
        block[f"go_{code}"] = pools.goclass.indices(
            _enum_value(getattr(locus, f"go_verdict_{namespace}")) for locus in loci
        )
    for namespace, code in namespace_codes.items():
        block[f"n_go_{code}"] = [
            getattr(locus, f"go_annotated_member_count_{namespace}") for locus in loci
        ]
    block["len_nt"] = [locus.median_gene_length_nt for locus in loci]
    block["len_iqr"] = [locus.gene_length_interquartile_range_nt for locus in loci]
    return block


#: ⛔⛔ **The order ``build_payload`` interns the vocabularies in — NOT the order it emits them.**
#: Its `for f, col, pool in (...)` loop runs before the `lists` dict literal is assembled, so the
#: pools are filled in *this* sequence while the emitted block happens to share it. They are two
#: different facts that agree today; this constant is the one that matters, and `u50` sitting third
#: rather than last is exactly the kind of detail a "tidy" alphabetical rewrite would destroy.
LIST_INTERN_ORDER: tuple[tuple[str, str], ...] = (
    ("sym", "GENE_SYMBOL"),
    ("prod", "PROTEIN_PRODUCT"),
    ("u50", "__crosstab__"),
    ("pfam", "PFAM_ARCHITECTURE"),
    ("cog", "COG_ORTHOGROUP"),
    ("ec", "EC_NUMBER"),
    ("kegg", "KEGG_ORTHOLOGY"),
    ("go", "GENE_ONTOLOGY_SLIM"),
)


def _run_lengths(per_locus: dict[int, list], n_loci: int, columns: dict[str, list]) -> dict:
    """A per-locus variable-length list as ``{"n": [...], <col>: [...]}`` — run lengths, not offsets.

    ⛔ ``n`` is length ``n_loci`` and carries a **zero for every locus with no entries**. A block
    listing only the loci that have something is a different, shorter array that silently
    misaligns every prefix sum after the first gap.
    """
    return {"n": [len(per_locus.get(index, ())) for index in range(n_loci)], **columns}


def annotation_rows(session: Session, pangenome_id: int) -> dict[str, list]:
    """Every annotation entry, grouped by list key, **in the payload's own row order**.

    ⛔ Two different orders, and using one for the other breaks the run lengths silently:

    * six vocabularies rank **within a locus** → ``(catalogue_ordinal, rank)``;
    * GO ranks within **(locus, namespace)** → ``(catalogue_ordinal, namespace, rank)``.

    ``build_payload`` uses `_ordered` for the first and `_ordered_ns` for the second, and its own
    comment says why: *"`top_go` ranks within (locus, namespace), so it sorts on namespace too or
    the run lengths break."*
    """
    from bacatlas_backend.models.locus_annotation import LocusAnnotationEntry

    rows = session.execute(
        select(LocusAnnotationEntry, Locus.catalogue_ordinal)
        .join(Locus, Locus.locus_id == LocusAnnotationEntry.locus_id)
        .where(Locus.pangenome_id == pangenome_id)
        .order_by(
            Locus.catalogue_ordinal,
            LocusAnnotationEntry.gene_ontology_namespace.nulls_first(),
            LocusAnnotationEntry.rank_within_locus,
        )
    ).all()
    grouped: dict[str, list] = {}
    for entry, ordinal in rows:
        grouped.setdefault(entry.annotation_kind.name, []).append((ordinal, entry))
    return grouped


def uniref_crosstab_rows(session: Session, pangenome_id: int) -> list:
    """The UniRef50 cross-tab in ``(catalogue_ordinal, rank)`` order — the payload's ``lists.u50``."""
    from bacatlas_backend.models.locus_annotation import LocusUnirefFamilyCrosstab

    return session.execute(
        select(LocusUnirefFamilyCrosstab, Locus.catalogue_ordinal)
        .join(Locus, Locus.locus_id == LocusUnirefFamilyCrosstab.locus_id)
        .where(Locus.pangenome_id == pangenome_id)
        .order_by(Locus.catalogue_ordinal, LocusUnirefFamilyCrosstab.rank_within_locus)
    ).all()


def list_block(session: Session, pangenome_id: int, n_loci: int, pools: CataloguePools) -> dict:
    """The whole ``lists`` block, interned in :data:`LIST_INTERN_ORDER`.

    ⚠ ``u50``'s five extra columns split two ways, and the split is the absence rule:
    ``prod``/``arch``/``sym`` are **pool indices**, so a missing value is ``-1``; ``npf``/``nsym``
    are **counts**, so a missing value is ``0``. Emitting ``-1`` for a count would make "no Pfam
    coverage" sort below zero, and emitting ``0`` for an index would name the pool's first string.
    """
    grouped = annotation_rows(session, pangenome_id)
    crosstab = uniref_crosstab_rows(session, pangenome_id)
    per_locus_crosstab: dict[int, list] = {}
    for family, ordinal in crosstab:
        per_locus_crosstab.setdefault(ordinal, []).append(family)

    interned: dict[str, list[int]] = {}
    for key, kind in LIST_INTERN_ORDER:
        pool = getattr(pools, key)
        if kind == "__crosstab__":
            interned[key] = pool.indices(family.uniref50_accession for family, _ in crosstab)
        else:
            interned[key] = pool.indices(entry.term_value for _, entry in grouped.get(kind, ()))

    # ⛔ The SECOND pass over the cross-tab, and it must come after every list above: it adds to
    # `prod`, `pfam` and `sym`, three pools that already hold entries. Running it earlier changes
    # three pools and every index into them.
    crosstab_product = pools.prod.indices(family.modal_bakta_product for family, _ in crosstab)
    crosstab_architecture = pools.pfam.indices(family.modal_pfam_architecture for family, _ in crosstab)
    crosstab_symbol = pools.sym.indices(family.modal_bakta_gene_symbol for family, _ in crosstab)

    def _grouped_lengths(kind: str) -> dict[int, list]:
        per_locus: dict[int, list] = {}
        for ordinal, entry in grouped.get(kind, ()):
            per_locus.setdefault(ordinal, []).append(entry)
        return per_locus

    block: dict[str, dict] = {}
    for key, kind in LIST_INTERN_ORDER:
        if kind == "__crosstab__":
            continue
        entries = grouped.get(kind, ())
        columns = {
            "idx": interned[key],
            "cnt": [entry.member_gene_count for _, entry in entries],
        }
        if key == "go":
            # GO carries its namespace per entry, so one list serves all three groups.
            columns["ns"] = [entry.gene_ontology_namespace for _, entry in entries]
        block[key] = _run_lengths(_grouped_lengths(kind), n_loci, columns)

    block["u50"] = _run_lengths(
        per_locus_crosstab,
        n_loci,
        {
            "idx": interned["u50"],
            "cnt": [family.member_gene_count for family, _ in crosstab],
            "prod": crosstab_product,
            "arch": crosstab_architecture,
            "npf": [family.pfam_annotated_member_count or 0 for family, _ in crosstab],
            "sym": crosstab_symbol,
            "nsym": [family.distinct_real_symbol_count or 0 for family, _ in crosstab],
        },
    )
    # Emitted in the payload's own key order, which is NOT the interning order above.
    return {key: block[key] for key in ("sym", "prod", "u50", "pfam", "cog", "ec", "kegg", "go")}


def genome_ordinal_map(session: Session, pangenome_id: int) -> dict[int, int]:
    """``genome_id`` → ``collection_genome_ordinal`` — the index ``arr.gid`` actually holds.

    ⛔⛔ **The schema stores a foreign key; the payload stores a position.** ``member_genome_ids``
    holds real ``genome_id`` values, which is right — a column of positions could not be joined to
    anything, and the ordinal is a property of the *collection*, not of the genome. But ``arr.gid``
    indexes into ``meta.genomes``, so the two are different integers over the same objects and
    reading one as the other silently names the wrong genome throughout: at 100 genomes both are
    small integers, both are in range, and every page renders.
    """
    from bacatlas_backend.models.genome_collection import GenomeCollectionMembership
    from bacatlas_backend.models.pangenome import Pangenome

    collection_id = session.execute(
        select(Pangenome.genome_collection_id).where(Pangenome.pangenome_id == pangenome_id)
    ).scalar_one()
    rows = session.execute(
        select(
            GenomeCollectionMembership.genome_id,
            GenomeCollectionMembership.collection_genome_ordinal,
        ).where(GenomeCollectionMembership.genome_collection_id == collection_id)
    ).all()
    ordinals = {genome_id: ordinal for genome_id, ordinal in rows}
    if sorted(ordinals.values()) != list(range(len(ordinals))):
        raise ValueError(
            f"collection {collection_id} has {len(ordinals):,} members whose ordinals are not "
            f"0…{len(ordinals) - 1} — every arr.gid would name the wrong genome"
        )
    return ordinals


def arrangement_block(session: Session, pangenome_id: int, loci: list[Locus]) -> dict:
    """The ``arr`` block — the **joint** view: whole ±5 neighbourhoods, and who carries each.

    ⛔ ``tot`` is read from ``locus.total_arrangement_count``, **not counted from the rows**. It is
    the TRUE total and the listed count is whatever the export cap admitted; conflating them is
    what lets a capped payload claim a locus has four neighbourhoods when it has thirty-seven.
    (The published catalogues happen to be uncapped — ``meta.top_arrangements`` is 0 and
    ``n == tot`` on all 33,201 loci — which is exactly why reading `tot` from a count would pass
    here and fail on the first capped export.)

    ⚠ ``gid`` rides on ``gen``: run *j* is exactly ``gen[j]`` long, concatenated in the same row
    order. Asserted before returning, because a mismatch makes the page read every arrangement's
    membership from the wrong offset — and every genome name it then shows is a real genome.
    """
    from bacatlas_backend.models.locus_arrangement import LocusArrangement

    rows = session.execute(
        select(LocusArrangement, Locus.catalogue_ordinal)
        .join(Locus, Locus.locus_id == LocusArrangement.locus_id)
        .where(Locus.pangenome_id == pangenome_id)
        .order_by(Locus.catalogue_ordinal, LocusArrangement.rank_within_locus)
    ).all()

    ordinal_of_genome_id = genome_ordinal_map(session, pangenome_id)
    counts = [0] * len(loci)
    genes, genomes, flips, vectors, members = [], [], [], [], []
    for arrangement, ordinal in rows:
        counts[ordinal] += 1
        genes.append(arrangement.member_gene_count)
        genomes.append(arrangement.member_genome_count)
        # ⚠ The payload stores an int, not a bool — `flip` is emitted through `_ints`.
        flips.append(int(arrangement.is_recorded_reverse_complement))
        vectors.extend(arrangement.neighbour_slot_codes)
        # ⚠ Sorted by ORDINAL, not by genome_id. Both are ascending in the published catalogues
        # because ingest happened to assign ids in accession order, so mapping alone would pass
        # here — and would silently emit a differently-ordered run the first time a collection is
        # built from an existing `genome` table. All 70,519 published runs are ascending ordinals.
        members.extend(sorted(ordinal_of_genome_id[gid] for gid in arrangement.member_genome_ids))

    block = {
        "n": counts,
        "cnt": genes,
        "gen": genomes,
        "flip": flips,
        "tot": [locus.total_arrangement_count for locus in loci],
        "vec": vectors,
        "gid": members,
    }
    if len(members) != sum(genomes):
        raise ValueError(
            f"arr.gid is {len(members):,} entries but arr.gen sums to {sum(genomes):,} — the page "
            "would read every arrangement's membership from the wrong offset"
        )
    return block


def context_block(session: Session, pangenome_id: int, loci: list[Locus]) -> dict:
    """The ``ctx`` block — the **marginal** view: one candidate at one position, ten slots per locus.

    ⛔ ``n`` and ``obs`` are flat arrays of ``n_loci × 10``, addressed as
    ``catalogue_ordinal * 10 + slot`` where ``slot`` is the position of the signed offset in
    ``OFFSETS`` — **not the offset itself**, which runs −5…−1, 1…5 and has no zero.

    ⛔⛔ ``obs`` is the denominator counted **before the top-N cut**, so it cannot be recovered by
    summing the occupant rows: a locus with more than ``top_neighbours`` candidates at a slot has a
    larger ``obs`` than its rows account for, and that difference is one of the payload's four
    distinct remainders. It is read from ``locus.context_observed_member_counts``, which is where
    the ingest stored it for exactly this reason.

    ⚠ ``nid`` holds payload **ordinals**, not database ids — the page holds indices everywhere.
    """
    from bacatlas_backend.models.locus_offset_occupant import LocusOffsetOccupant

    _, offsets, _, _, _ = _payload_constants()
    slot_of = {offset: index for index, offset in enumerate(offsets)}
    n_slots = len(offsets)
    ordinal_of_locus_id = {locus.locus_id: locus.catalogue_ordinal for locus in loci}

    rows = session.execute(
        select(LocusOffsetOccupant, Locus.catalogue_ordinal)
        .join(Locus, Locus.locus_id == LocusOffsetOccupant.locus_id)
        .where(Locus.pangenome_id == pangenome_id)
        .order_by(
            Locus.catalogue_ordinal,
            LocusOffsetOccupant.signed_offset,
            LocusOffsetOccupant.rank_within_offset,
        )
    ).all()

    slot_counts = [0] * (len(loci) * n_slots)
    neighbours, gene_counts, same_strand = [], [], []
    for occupant, ordinal in rows:
        slot_counts[ordinal * n_slots + slot_of[occupant.signed_offset]] += 1
        neighbours.append(ordinal_of_locus_id[occupant.neighbour_locus_id])
        gene_counts.append(occupant.member_gene_count)
        same_strand.append(occupant.same_strand_member_count)

    observed: list[int] = []
    for locus in loci:
        observed.extend(locus.context_observed_member_counts)
    return {
        "n": slot_counts,
        "obs": observed,
        "nid": neighbours,
        "cnt": gene_counts,
        "same": same_strand,
    }


#: The payload stores the variance score ×1000 and **clipped at 2.0**, as an integer. The column
#: holds the score itself, so both steps are re-applied here rather than stored twice.
VARIANCE_SCALE = 1000
VARIANCE_CLIP = 2.0



def _rounded(value: float | None, decimals: int) -> float | None:
    """`round`, but NULL stays NULL — an unmeasured baseline is not a baseline of zero."""
    return None if value is None else round(float(value), decimals)


def gaps_block(session: Session, pangenome_id: int, loci: list[Locus]) -> dict:
    """The ``gaps`` block — adjacencies, their length spread, and what was seen inside them.

    ⛔⛔ **The sparse triple is where the meaning INVERTS, and it inverts in only one direction.**
    The payload ships ``vi``/``vd``/``vmd`` sparsely because a JSON array cannot be sparse, so
    *absent* there means **every genome agrees** — 86.1 % of ecoli's gaps. The column is dense, so
    ``0.0`` carries that agreement explicitly and ``NULL`` means *this run did not measure variance*
    — a state the sparse form could not express at all. Re-sparsifying is therefore
    ``score > 0``, under which a measured zero and an unmeasured NULL both drop out, which is
    exactly what the published file did. **The database is strictly more informative than the
    payload here, and going back is lossy on purpose.**

    ⚠ ``q1 == q3`` means the MIDDLE HALF agrees; only ``mn == mx`` certifies "every genome agrees".
    Both pairs are carried because an earlier page reported the first as the second.
    """
    from bacatlas_backend.models.intergenic_gap import IntergenicGap, IntergenicGapFeature

    ordinal_of_locus_id = {locus.locus_id: locus.catalogue_ordinal for locus in loci}
    gaps = list(
        session.execute(
            select(IntergenicGap)
            .where(IntergenicGap.pangenome_id == pangenome_id)
            .order_by(IntergenicGap.intergenic_gap_id)
        ).scalars()
    )
    if not gaps:
        # ⚠ `{}` and not an empty scaffold: a cohort without the non-coding pass has nothing to say
        # here, and the page special-cases an absent block rather than an empty one.
        return {}

    features_by_gap: dict[int, list] = {}
    for feature in session.execute(
        select(IntergenicGapFeature)
        .join(IntergenicGap, IntergenicGap.intergenic_gap_id == IntergenicGapFeature.intergenic_gap_id)
        .where(IntergenicGap.pangenome_id == pangenome_id)
        .order_by(IntergenicGapFeature.intergenic_gap_id, IntergenicGapFeature.rank_within_gap)
    ).scalars():
        features_by_gap.setdefault(feature.intergenic_gap_id, []).append(feature)

    # ⚠ Two pools LOCAL to this block, not the global `strings` table — `intergenic_block` builds
    # its own `_Intern` pair, so their indices address `gaps.labels`/`gaps.types` and nothing else.
    labels, types = StringPool(), StringPool()
    per_gap_counts, label_indices, type_indices, feature_counts = [], [], [], []
    for gap in gaps:
        features = features_by_gap.get(gap.intergenic_gap_id, ())
        per_gap_counts.append(len(features))
        for feature in features:
            label_indices.append(labels.index(feature.feature_label))
            type_indices.append(types.index(feature.feature_type))
            feature_counts.append(feature.observed_genome_count)

    varying = [
        index
        for index, gap in enumerate(gaps)
        if gap.length_variance_score is not None and gap.length_variance_score > 0
    ]
    return {
        "a": [ordinal_of_locus_id[gap.flanking_locus_id_a] for gap in gaps],
        "b": [ordinal_of_locus_id[gap.flanking_locus_id_b] for gap in gaps],
        "n": [gap.observed_genome_count for gap in gaps],
        "nt": [gap.median_signed_length_nt for gap in gaps],
        "q1": [gap.quartile1_signed_length_nt for gap in gaps],
        "q3": [gap.quartile3_signed_length_nt for gap in gaps],
        "mn": [gap.minimum_signed_length_nt for gap in gaps],
        "mx": [gap.maximum_signed_length_nt for gap in gaps],
        "vi": varying,
        "vd": [
            int(round(min(gaps[index].length_variance_score, VARIANCE_CLIP) * VARIANCE_SCALE))
            for index in varying
        ],
        "vmd": [gaps[index].modal_length_nt for index in varying],
        "n_feat": [gap.distinct_named_feature_count for gap in gaps],
        "fn": per_gap_counts,
        "flab": label_indices,
        "ftyp": type_indices,
        "fcnt": feature_counts,
        "labels": labels.values,
        "types": types.values,
    }


# ⛔ `map_representation_blocks` was DELETED on 2026-09-24 with the neighbourhood map it rebuilt: a UMAP
# over every locus's MEDOID plus that medoid's 6×6 local geometry. `similarity_block` below is what
# replaced it, and it is a different measurement rather than a renaming.
#
# ⚠ Two things that block knew, worth keeping if any packed array returns here. The base64 arrays went
# through nuna's own `_b64_i16`/`_b64_i32` rather than a second little-endian packer, because an
# endianness mistake produces a valid base64 string, a valid picture, and the wrong loci. And its `-1`
# dropped a SLOT rather than a rank, so the surviving slot indices were what addressed the triangle —
# reading by rank drew one locus's distances on another, and the picture still looked like a picture.


def _headline_value(metric_name: str, value):
    """An audit headline value with its JSON type restored — **by integrality, not by name**.

    ⛔ `pangenome_evaluation.numeric_value` is one numeric column, so *whether the audit wrote an
    int or a float is genuinely not stored*. This is the one place in the payload where that bites:
    `17531` must not come back as `17531.0`, and `0.032587` must not come back as `0`.

    ⚠ A suffix rule was tried first and was wrong twice in twenty-one keys —
    `split_gene_rate_excl_singletons` does not end in `_rate` (so it was truncated to `0`) and
    `over_merge_gene_rate_num` contains `_rate` but is a count. Integrality is a property of the
    *value* rather than of the spelling, and it is right on every key here.

    ⚠ Its one blind spot, stated rather than discovered: **a rate that is exactly 0.0 or 1.0** comes
    back as an int. Only `no_homology_gene_rate` is such a value today, and only because the tier it
    counts was retired — so it already sits inside a recorded difference.
    """
    if value is None:
        return None
    number = float(value)
    return int(number) if number == int(number) else number


def _arrangement_cap(loci: list[Locus], listed_counts: list[int]) -> int:
    """The `--top-arrangements` the export ran with, recovered from what it shipped.

    ⭐ **0 means every arrangement was kept**, and it is a different statement from "four were": the
    page tells the reader whether a rarer neighbourhood exists or merely was not shipped.

    A locus is evidence of a cap only when its `total_arrangement_count` exceeds what was listed,
    and then the cap is what it was cut to. Uncapped catalogues have no such locus, so the answer is
    0 — which is what both published catalogues shipped, on all 33,201 loci.
    """
    # ⛔ `listed > 0` is not a tidy-up. **847 loci have no arrangement row at all** — their genes
    # never reached a window, because `ac` is an inner join on coordinates — and a locus that listed
    # NONE reveals no cap. Without the guard such a locus reads as "cut to zero", which is a cap of
    # 0, which is the encoding for *uncapped*: the two ends of the scale collide on the case that
    # actually occurs. (It happens to be harmless here only because all 847 also have `tot == 0`.)
    capped = [
        listed
        for locus, listed in zip(loci, listed_counts, strict=True)
        if listed > 0 and locus.total_arrangement_count > listed
    ]
    if not capped:
        return 0
    cap = max(capped)
    if min(capped) != cap:
        raise ValueError(
            f"loci were cut to between {min(capped)} and {cap} arrangements — no single "
            "--top-arrangements produces that, so the cap cannot be stated"
        )
    return cap


def similarity_block(session: Session, pangenome_id: int, loci: list[Locus], *, model_label: str) -> dict:
    """The `sim` block — set-to-set similarity per representation, **keyed by rep, not a list**.

    ⛔ **This replaced `map_representation_blocks` and `null_baseline_block` together**, because in
    the payload they replaced one card and its axis. The old pair emitted a UMAP over every locus's
    MEDOID, that medoid's 6×6 local geometry, and a distribution of random MEDOID pairs.

    ⚠ The two base64 arrays are emitted through nuna's own `_b64_i16`/`_b64_i32`, imported rather
    than reimplemented — a second little-endian packer is a second thing to get wrong, and an
    endianness mistake produces a valid base64 string and the wrong loci.

    ⚠ `near_i` is `n_loci × k` **raveled** and **-1 means absent**, which here is a genuinely ragged
    list rather than the retired slot-dropping sentinel: a locus whose shortlist held fewer than `k`
    others simply has fewer rows in `locus_nearest_locus`, and the remaining slots pad.
    """
    from nuna.tl.locus_browser.export_payload import _b64_i16, _b64_i32

    from bacatlas_backend.models.locus_similarity import (
        NEAREST_LOCUS_COUNT,
        LocusNearestLocus,
        LocusSimilarity,
        PangenomeSimilarityBaseline,
    )

    baselines = list(
        session.execute(
            select(PangenomeSimilarityBaseline)
            .where(PangenomeSimilarityBaseline.pangenome_id == pangenome_id)
            .order_by(PangenomeSimilarityBaseline.pangenome_similarity_baseline_id)
        ).scalars()
    )
    if not baselines:
        return {}

    ordinal_of = {locus.locus_id: locus.catalogue_ordinal for locus in loci}
    out: dict[str, dict] = {}
    for baseline in baselines:
        representation = baseline.representation
        rep = representation.value
        rows = session.execute(
            select(LocusSimilarity, Locus.catalogue_ordinal)
            .join(Locus, Locus.locus_id == LocusSimilarity.locus_id)
            .where(Locus.pangenome_id == pangenome_id, LocusSimilarity.representation == representation)
            .order_by(Locus.catalogue_ordinal)
        ).all()
        if len(rows) != len(loci):
            raise ValueError(
                f"{rep}: {len(rows):,} similarity rows for {len(loci):,} loci — every array in this "
                "block is positional, so a missing row shifts the whole catalogue by one"
            )

        columns = {
            "within": [row.within_similarity for row, _ in rows],
            "near": [row.nearest_similarity for row, _ in rows],
            "weak_in": [row.weak_own_similarity for row, _ in rows],
            "weak_out": [row.weak_other_similarity for row, _ in rows],
            "own": [row.own_neighbour_fraction for row, _ in rows],
        }

        near_index = [-1] * (len(loci) * NEAREST_LOCUS_COUNT)
        near_value = [0] * (len(loci) * NEAREST_LOCUS_COUNT)
        for nearest, ordinal in session.execute(
            select(LocusNearestLocus, Locus.catalogue_ordinal)
            .join(Locus, Locus.locus_id == LocusNearestLocus.locus_id)
            .where(Locus.pangenome_id == pangenome_id, LocusNearestLocus.representation == representation)
        ).all():
            slot = ordinal * NEAREST_LOCUS_COUNT + (nearest.rank - 1)
            near_index[slot] = ordinal_of[nearest.neighbour_locus_id]
            near_value[slot] = _quantised_cosine(nearest.cross_similarity)

        out[rep] = {
            "form": baseline.similarity_form,
            "k": baseline.nearest_locus_count,
            "cos_scale": COSINE_SCALE_FACTOR,
            **{key: _floats(values) for key, values in columns.items()},
            "near_i": _b64_i32(near_index),
            "near_v": _b64_i16(near_value),
            "floor": _rounded(baseline.floor_median, SIMILARITY_FLOOR_DECIMALS),
            "floor_q": (
                None
                if baseline.floor_p25 is None
                else [
                    _rounded(baseline.floor_p25, SIMILARITY_FLOOR_DECIMALS),
                    _rounded(baseline.floor_p75, SIMILARITY_FLOOR_DECIMALS),
                    _rounded(baseline.floor_p99, SIMILARITY_FLOOR_DECIMALS),
                ]
            ),
            "n_measurable": baseline.measurable_locus_count,
            "knn_k": baseline.neighbour_knn_k,
            # ⚠ A fact about the EXPORT rather than about the measurement, so it is rebuilt from the
            # model label rather than stored — the same distinction `META_FIELD_SOURCES` applies.
            "source": f"{model_label}_cluster_similarity_{rep}.csv",
        }
    return out


#: The payload stores cosines as int16 hundredths-of-a-basis-point; nuna's `_COS`, restated so this
#: side does not import a private name that has already been deleted once.
COSINE_SCALE_FACTOR = 10_000
SIMILARITY_FLOOR_DECIMALS = 6


def _quantised_cosine(value: float | None) -> int:
    if value is None:
        return 0
    scaled = round(float(value) * COSINE_SCALE_FACTOR)
    return max(-COSINE_SCALE_FACTOR, min(COSINE_SCALE_FACTOR, scaled))


def _floats(values, decimals: int = 4) -> list[float | None]:
    """Round to `decimals` places; `None` stays `None` — the payload's own `_floats` convention."""
    return [None if value is None else round(float(value), decimals) for value in values]


def meta_block(
    session: Session,
    pangenome_id: int,
    loci: list[Locus],
    *,
    species_key: str,
    listed_arrangement_counts: list[int],
) -> dict:
    """The `meta` block — the page's header, footer and genome vocabulary.

    ⛔ `genomes` is the **vocabulary `arr.gid` indexes into**, ordered by
    `collection_genome_ordinal` and never re-sorted here. Re-deriving it by sorting accessions is
    the implicit ordering this project punishes: it agrees today and stops agreeing the moment a
    collection is built in any other order.
    """
    from nuna.tl.locus_browser.export_payload import (
        AUDIT_HEADLINE_KEYS,
        SEQ_FLANK,
        TOP_NEIGHBOURS,
    )

    from bacatlas_backend.models.genome import Genome
    from bacatlas_backend.models.genome_collection import GenomeCollectionMembership
    from bacatlas_backend.models.pangenome import Pangenome, PangenomeEvaluation
    from bacatlas_backend.models.pathogen_species import PathogenSpecies

    _, offsets, band_order, _, policy = _payload_constants()
    pangenome = session.get(Pangenome, pangenome_id)
    species = session.get(PathogenSpecies, pangenome.pathogen_species_id)

    vocabulary = [
        accession
        for accession, in session.execute(
            select(Genome.sample_id)
            .join(
                GenomeCollectionMembership,
                GenomeCollectionMembership.genome_id == Genome.genome_id,
            )
            .where(
                GenomeCollectionMembership.genome_collection_id == pangenome.genome_collection_id
            )
            .order_by(GenomeCollectionMembership.collection_genome_ordinal)
        ).all()
    ]

    evaluation_detail: dict[str, str | None] = {}
    evaluation_detail: dict[str, str | None] = {}
    headline_rows: dict[str, float | None] = {}
    for row in session.execute(
        select(PangenomeEvaluation).where(PangenomeEvaluation.pangenome_id == pangenome_id)
    ).scalars():
        headline_rows[row.metric_name] = row.numeric_value
        evaluation_detail[row.metric_name] = row.detail

    # ⭐ The audit's own label, as the audit wrote it — read from the database rather than from the
    # checked-in `published_catalogues` triple. The registry is right about these two catalogues and
    # is the correct authority for *addressing artifacts*; it is the wrong authority for what a
    # pangenome IS, and an instrument that asks it here would keep passing after the database
    # started disagreeing with it.
    # ⚠ `nuna_model.label` is NOT this string — it carries no species prefix.
    model_label = evaluation_detail.get("label")
    if not model_label:
        raise ValueError(
            f"pangenome {pangenome_id} has no `label` evaluation row — `meta.model_label` and the "
            "audit block's own label would have to be guessed"
        )
    headline = {
        key: _headline_value(key, headline_rows[key])
        for key in AUDIT_HEADLINE_KEYS
        if key in headline_rows
    }

    # ⚠ Derived from the tier the audit assigned, so the lists cannot drift from the counts in the
    # report — and kept as TWO lists, because a Pfam conflict is weaker evidence than an over-merge
    # and merging them would state a conflict as a grade.
    failure_tiers = set(policy["failure_tiers"])
    contested = policy["contested_pfclass"]
    return {
        "species": species.scientific_name,
        # ⛔ `dset` is the STATIC PUBLICATION key — the `published.tsv` row this payload is the page for,
        # which is what `render_page` matches to know which entry in the species picker is "me". It is
        # NOT `pangenome.catalogue_key`, and the two are different vocabularies on purpose:
        #
        #     published.tsv (static site)   ecoli          kp            ← nuna4 keeps the legacy urls
        #     pangenome.catalogue_key       ecoli-nuna4    kp-nuna4      ← the service's address
        #
        # The service cannot key nuna4 as `ecoli`, because its resolver reads the hyphen: a key with one
        # pins a catalogue, a bare species follows its default. Make them equal and that rule collapses.
        # The static site cannot key nuna4 as `ecoli-nuna4` either, because `ecoli.html` is a published
        # url and renaming it breaks every link that exists.
        #
        # ⚠ An earlier note here promised "a browser key on the PANGENOME" as the fix. It was WRONG, and
        # the parity suite caught it: emitting `catalogue_key` makes a rebuilt E. coli payload say
        # `ecoli-nuna4` against a frozen page saying `ecoli`. The divergence is real; it is just not
        # this field's to resolve. The day a nuna5 catalogue is published STATICALLY it gains a
        # `published.tsv` row (`ecoli-nuna5`) and that row's key is what belongs here.
        "dset": species.species_key,
        "model_id": pangenome.run_id,
        "model_label": model_label,
        # ⚠ `built_at` is `Mapped[str | None]` on a `String(32)` column, so `.date()` on it is an
        # AttributeError waiting for the day something writes it. It is NULL on all four loaded
        # catalogues today, which is the only reason this has never fired. Take the leading date of
        # whatever string is there; a datetime, if one is ever mapped, still answers `.date()`.
        "built": _built_date(pangenome.built_at),
        "git_sha": pangenome.git_sha,
        "provenance": pangenome.provenance_rows,
        "omitted": pangenome.omitted_sections or {},
        "n_genomes": pangenome.genome_count,
        "n_genes": pangenome.gene_count,
        "n_loci": len(loci),
        "offsets": list(offsets),
        "bands": list(band_order),
        "genomes": vocabulary,
        "top_neighbours": int(TOP_NEIGHBOURS),
        # ⛔ Read from the DATA, never from `TOP_ARRANGEMENTS`. The constant is nuna's *default* (4)
        # and the published exports overrode it to run uncapped, so taking the default emits `4` for
        # a catalogue that shipped `0` — and 0 vs 4 is the difference between "a rarer neighbourhood
        # exists but was not shipped" and "there are no others", which is the whole reason the field
        # is in the payload.
        "top_arrangements": _arrangement_cap(loci, listed_arrangement_counts),
        "seq": {"dir": "data/seq", "flank": SEQ_FLANK},
        "policy": policy,
        "audit": {
            "label": model_label,
            "sources": [
                f"{model_label}_pfam_concordance.tsv",
                f"{model_label}_audit_summary.json",
            ],
            "headline": headline,
            "failures": [
                locus.catalogue_ordinal for locus in loci if locus.collapse_tier in failure_tiers
            ],
            "contested": [
                locus.catalogue_ordinal
                for locus in loci
                if locus.pfam_concordance_class == contested
            ],
        },
    }


def build_payload_from_database(session: Session, species_key: str) -> dict:
    """The whole payload, rebuilt — **and its own interning order checked before it is returned**.

    ⭐ The self-check is the point: a serialiser that emitted the blocks in the wrong order would
    otherwise produce thousands of changed indices in a diff nobody can read. Here it raises, naming
    the pool and the column.
    """
    from bacatlas_backend.instruments.payload_reproduction import verify_intern_walk
    from bacatlas_backend.models.pathogen_species import PathogenSpecies

    schema_version, _, _, _, _ = _payload_constants()
    pangenome_id = session.execute(
        select(PathogenSpecies.default_pangenome_id).where(
            PathogenSpecies.species_key == species_key
        )
    ).scalar_one()
    loci = load_catalogue_loci(session, pangenome_id)
    pools = CataloguePools()

    # ⛔ ORDER: nodes, then lists, then everything that does not intern. See the module docstring.
    nodes = node_block(loci, pools)
    lists = list_block(session, pangenome_id, len(loci), pools)
    # ⚠ None of the blocks below intern, so their order is free — `arr` is built first only because
    # `meta.top_arrangements` is recovered from what it listed.
    arrangements = arrangement_block(session, pangenome_id, loci)
    gaps = gaps_block(session, pangenome_id, loci)

    payload = {
        "schema": schema_version,
        "meta": meta_block(
            session,
            pangenome_id,
            loci,
            species_key=species_key,
            listed_arrangement_counts=arrangements["n"],
        ),
        "strings": pools.tables(),
        "nodes": nodes,
        "lists": lists,
        "arr": arrangements,

        "ctx": context_block(session, pangenome_id, loci),
        **({"gaps": gaps} if gaps else {}),
    }
    similarity = similarity_block(
        session, pangenome_id, loci, model_label=payload["meta"]["model_label"]
    )
    if similarity:
        payload["sim"] = similarity
    failures = verify_intern_walk(payload)
    if failures:
        raise ValueError(
            "the rebuilt payload's own interning order is wrong — the blocks were emitted in a "
            "different order from `build_payload`:\n  " + "\n  ".join(failures)
        )
    return payload
