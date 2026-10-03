r"""Does ESM carry function ACROSS species? — reliable first, then worthwhile.

Reads the cross-species node map that `nuna.tl.probe.cross_species_neighbours` writes (every
*K. pneumoniae* locus's nearest *E. coli* loci, measured over every cross gene pair) and asks David's
two questions, in his order:

    "I would first study the result at ESM similarity thresholds to check that the numbers are
     trustworthy at the same thresholds as within Kp. Then I would look at the total number of inferred
     labels it achieves — this is the main question, is it worthwhile. Only last, after showing
     worthwhile, would I worry about checking the Uniref50 bridge."   (David, 2026-10-03)

⛔⛔ **Nothing here goes near the web page.** No payload change, no schema bump, no card, no ingest.
This is a measurement script; whether cross-species transfer gets built is David's decision from the
numbers.

⭐ **The ladder, the tiers, the level folding and the agreement arithmetic are IMPORTED from
`instruments/annotation_transfer.py`** — the same code the API serves the within-species card from. The
whole question is whether *the existing ladder* transfers, so re-deriving the bands here would make the
comparison unanswerable. ⛔ **Do not move the cuts.**

⭐ **The within-kp column is recomputed from the database on every run, never transcribed.** It is the
thing the cross-species column is read against, and a frozen baseline decays silently — nuna
`PROJECT_STATE.md` §6, 2026-10-03: a stale pin and a stale oracle were wrong in the same direction and
agreed, so the test that existed to catch it reported nothing.

⚠ **One statistical parameter changes, and it is a decision rather than an implementation detail.**
Within a species the chance baseline is *"a random **other** node of this catalogue"*. Across species
the comparator is *"a random annotated node of the **donor** catalogue"* — so the pool is the donor's
claims alone (never the union, which would dilute it with the recipient's own distribution) and the
self-exclusion is dropped, because the focal recipient node is not a donor. `tally_cells(donor_pool=)`
binds both consequences to one argument. A wrong null rescales every lift and leaves agreement
untouched, so nothing in the output would look wrong.

Run from `backend/` with `BACATLAS_DATABASE_URL` set and both `*-nuna4` catalogues loaded:

    python scripts/measure_cross_species_transfer.py \\
        ../../nuna/data/proc/analysis/cross_species/ecoli_to_kp_<model>_cross_species_neighbours.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

from sqlalchemy import create_engine, text
from sqlalchemy.orm import Session

from bacatlas_backend.instruments.annotation_transfer import (
    LEVEL_LABEL,
    MIN_PAIRS,
    NOT_CALLED,
    QUOTED_LEVEL,
    TIERS,
    Cell,
    quotable_level,
    tally_cells,
    tier_for,
)
from bacatlas_backend.models.enumerations import AnnotationKind

#: The three linkage rules the map carries for every pair. ⭐ All three are measured, because David
#: deferred the choice and the map already holds them: `cross` is the median over every cross gene
#: pair, `cross_max` the closest single gene pair (his anisotropy hypothesis — *"mapping by closest
#: member is probably the best gauge of community"*), `cross_p75` between them.
RULES = ("cross", "cross_p75", "cross_max")
#: Worked cases printed under each reach table. A table of counts is not something that can be
#: checked by eye; the within-species result became credible when it could be read on named nodes.
EXAMPLES = 8
KINDS = (
    AnnotationKind.COG_ORTHOGROUP,
    AnnotationKind.GENE_ONTOLOGY_SLIM,
    AnnotationKind.KEGG_ORTHOLOGY,
    AnnotationKind.EC_NUMBER,
)
#: Catalogue keys, keyed by the species token the map's column names use.
CATALOGUE = {"ecoli": "ecoli-nuna4", "kp": "kp-nuna4"}


# ================================================================================ the map, off disk
def read_map(path: str, *, donor: str, recipient: str) -> tuple[list[dict], dict]:
    """The cross-species map plus its provenance — `(rows, provenance)`.

    ⛔ **Read with `csv`, so every value arrives as `str` and stays one.** `locus.node_label` is TEXT
    and node labels *look* numeric (`"1098"`); a reader that inferred dtypes would turn them into
    int64, and the join against the database would then match nothing at all while looking perfectly
    reasonable. The counts and similarities are cast explicitly, one column at a time.

    ⛔⛔ **THE DIRECTION IS NOT IN THE HEADER, AND A TEST IS WHAT ESTABLISHED THAT.** The header is
    symmetric — a map in either direction carries both `*_node` and both `n_*_genes` columns — so a
    header check reads an `ecoli -> kp` map as a `kp -> ecoli` one without complaint, and every rank
    then means the opposite of what it says: `rank_cross == 1` is *"this donor is this recipient's
    best"*, which inverted is a different and wrong claim. The direction is therefore taken from the
    run's provenance JSON beside the TSV, which records it explicitly, and the header check is kept
    as a second, independent assertion that catches a truncated or older map.

    ⚠ **The sidecar is REQUIRED, not optional.** It carries the shortlist's recall audits at both
    linkage rules, and `~/.claude/CLAUDE.md` is explicit that a comparison asserts its own coverage
    before it reports a difference — a map whose recall was never measured is not something to report
    transfer rates from.
    """
    sidecar = Path(path).with_suffix(".json")
    if not sidecar.exists():
        raise SystemExit(
            f"FATAL: no provenance beside the map — expected {sidecar}.\n"
            "  It is what records the DIRECTION (the TSV header is symmetric and cannot), the cap the"
            " walk is bounded by, and the shortlist recall at both linkage rules."
        )
    provenance = json.loads(sidecar.read_text())
    stated = (provenance.get("donor"), provenance.get("recipient"))
    if stated != (donor, recipient):
        raise SystemExit(
            f"FATAL: {sidecar.name} says {stated[1]} asks and {stated[0]} answers; this run was told "
            f"{recipient} asks and {donor} answers.\n"
            "  Read the wrong way round every rank means the opposite of what it says, and every "
            "label would transfer in the wrong direction while the table looked entirely reasonable."
        )

    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []
        wanted = [
            f"{recipient}_node",
            f"{donor}_node",
            f"n_{recipient}_genes",
            *RULES,
            *(f"rank_{rule}" for rule in RULES),
        ]
        missing = [column for column in wanted if column not in header]
        if missing:
            raise SystemExit(
                f"FATAL: {path} does not look like a {recipient}->{donor} map — missing {missing}.\n"
                f"  Its header is: {header}\n"
                f"  Pass --donor/--recipient matching the run that wrote it; a reversed map read "
                "forwards would transfer every label the wrong way and look fine."
            )
        rows = []
        for row in reader:
            rows.append(
                {
                    "recipient_node": row[f"{recipient}_node"],
                    "donor_node": row[f"{donor}_node"],
                    "recipient_genes": int(row[f"n_{recipient}_genes"]),
                    **{rule: float(row[rule]) for rule in RULES},
                    **{f"rank_{rule}": int(row[f"rank_{rule}"]) for rule in RULES},
                }
            )
    if not rows:
        raise SystemExit(f"FATAL: {path} has a valid header and no rows")
    return rows, provenance


# ============================================================================= the database side
def locus_index(session: Session, catalogue: str) -> dict[str, tuple[int, int]]:
    """`node_label` → `(locus_id, member_gene_count)` for one catalogue. Labels are TEXT, both ways."""
    rows = session.execute(
        text("""
            select l.node_label, l.locus_id, l.member_gene_count
              from locus l join pangenome p using (pangenome_id)
             where p.catalogue_key = :catalogue
        """),
        {"catalogue": catalogue},
    ).all()
    if not rows:
        raise SystemExit(f"FATAL: catalogue '{catalogue}' holds no loci — is it loaded?")
    return {str(label): (locus_id, members) for label, locus_id, members in rows}


def claims(session: Session, catalogue: str, kind: AnnotationKind):
    """Per locus, its top-rank call folded onto every rung — via the instrument, not a second copy.

    ⭐ `_claims` is the single implementation of the fold, and it is what `compute_calibration`
    (hence the live card) uses. `measure_cog_function_inference.py` reimplements the same SQL inline,
    which is a drift risk this script deliberately does not inherit.
    """
    from bacatlas_backend.instruments.annotation_transfer import _claims

    pangenome_id = session.execute(
        text("select pangenome_id from pangenome where catalogue_key = :catalogue"),
        {"catalogue": catalogue},
    ).scalar_one()
    levels, _terms, locus_count = _claims(session, pangenome_id, kind)
    return levels, locus_count


def within_species_edges(session: Session, catalogue: str) -> list[tuple[int, int, float]]:
    """Every stored ESM neighbour edge of one catalogue — the within-species column, recomputed.

    ⭐ **Every rank, not rank 1 only**, matching `compute_calibration`: 65 % of assignments come from
    ranks 2-5, and at fixed cosine the rank effect was <= 0.3 pp in the strong tiers and flipped sign
    between species in the weak ones. Pooling here is what makes the two columns the same measurement.
    """
    return [
        (locus_id, neighbour_id, cosine)
        for locus_id, neighbour_id, cosine in session.execute(
            text("""
                select n.locus_id, n.neighbour_locus_id, n.cross_similarity
                  from locus_nearest_locus n
                  join locus l on l.locus_id = n.locus_id
                  join pangenome p on p.pangenome_id = l.pangenome_id
                 where p.catalogue_key = :catalogue and n.representation = 'ESM'
                   and n.cross_similarity is not null
            """),
            {"catalogue": catalogue},
        ).all()
    ]


# ================================================================= SECTION 1 — is it RELIABLE?
def rate(cell: Cell | None) -> str:
    """A rate with its n, or the n alone below the floor. ⛔ Never a rate off 7 pairs."""
    if cell is None:
        return "         —"
    if cell.agreement is None:
        return f"  n={cell.pairs:<6}"
    return f"{cell.agreement:>6.1%} n={cell.pairs:<5}"


def reliability(
    cross_cells: dict[tuple[str, int], Cell],
    within_cells: dict[tuple[str, int], Cell],
    kind: AnnotationKind,
    *,
    rule: str,
    donor: str,
    recipient: str,
) -> None:
    """The one table this whole exercise exists to produce: the same tiers, two columns.

    ⛔ **Lift, never an odds ratio** (David, 2026-10-03). Where chance is 0.2 % the two read 394x and
    3,432x, and the odds ratio is unbounded as chance → 0 — it cannot be said as *"x more likely than
    a random annotated node"*, which is the sentence the card actually prints. ⚠ And the two columns'
    lifts are against DIFFERENT nulls by construction (donor pool vs own catalogue), which is exactly
    why both chance figures are printed rather than just the ratio.
    """
    print(f"\n  {kind.value} · ranked by '{rule}' — cross-species beside within-{recipient}")
    print(
        f"     {'tier':<11}{'rung':<18}{recipient + '->' + donor:>20}{'chance':>9}{'lift':>8}"
        f"{'   within-' + recipient:>20}{'chance':>9}{'lift':>8}"
    )
    for tier, _low, _high in TIERS:
        for level in sorted(LEVEL_LABEL[kind], reverse=True):
            here, there = cross_cells.get((tier, level)), within_cells.get((tier, level))
            if here is None and there is None:
                continue
            quoted = " ⭐" if QUOTED_LEVEL[kind].get(tier) == level else "   "
            print(
                f"     {tier:<11}{LEVEL_LABEL[kind][level]:<18}{rate(here):>20}"
                f"{(f'{here.chance:.3%}' if here else '—'):>9}"
                f"{(f'{here.lift:.0f}x' if here and here.lift else '—'):>8}"
                f"{rate(there):>20}"
                f"{(f'{there.chance:.3%}' if there else '—'):>9}"
                f"{(f'{there.lift:.0f}x' if there and there.lift else '—'):>8}{quoted}"
            )
    print(f"     ⭐ = the rung this tier is allowed to quote. A cell below {MIN_PAIRS} pairs shows its n and NO rate.")


def transfer_gate(cross_cells, within_cells, kind: AnnotationKind, *, recipient: str) -> list[str]:
    """Does the ladder TRANSFER? One line per quoted tier, as a difference with both n's.

    ⛔ Returns the warnings rather than printing a verdict: *"Gather the data, then STOP"*
    (`~/.claude/CLAUDE.md`). Whether a gap of this size means the ladder does not transfer is David's
    reading, not this script's — what belongs here is the gap, stated plainly, with what it rests on.
    """
    notes = []
    for tier, _low, _high in TIERS:
        level = QUOTED_LEVEL[kind].get(tier)
        here, there = cross_cells.get((tier, level)), within_cells.get((tier, level))
        if here is None or there is None or here.agreement is None or there.agreement is None:
            notes.append(
                f"     {tier:<11}{LEVEL_LABEL[kind][level]:<18} not measurable on both sides"
                f" (cross n={here.pairs if here else 0}, within n={there.pairs if there else 0})"
            )
            continue
        # in percentage POINTS, not as a ratio: 97 % against 90 % is +7 pp, and a reader who sees
        # "+0.08" cannot tell a 8 pp gap from an 8 % relative one
        delta_pp = 100 * (here.agreement - there.agreement)
        notes.append(
            f"     {tier:<11}{LEVEL_LABEL[kind][level]:<18}{here.agreement:>7.1%} vs "
            f"{there.agreement:>7.1%} within-{recipient}  →  {delta_pp:+6.1f} pp"
            f"   (n={here.pairs:,} / {there.pairs:,})"
        )
    return notes


# ================================================================ SECTION 2 — is it WORTHWHILE?
def candidates_from_map(rows: list[dict], *, rule: str, recipient_ids, donor_ids) -> dict:
    """The cross-species map as `recipient_locus_id -> [(rank, donor_locus_id, cosine)]` under one rule."""
    out = defaultdict(list)
    for row in rows:
        recipient_id = recipient_ids[row["recipient_node"]][0]
        out[recipient_id].append((row[f"rank_{rule}"], donor_ids[row["donor_node"]][0], row[rule]))
    return out


def candidates_from_database(session: Session, catalogue: str) -> dict:
    """The stored WITHIN-species neighbour graph, in the same shape — the third column, computed.

    ⭐ **The within-species ladder's reach is measured by the same `walk`, not quoted from the
    write-up.** It is the thing the cross-species reach has to be read against, and a number lifted
    out of a document keeps reading plausibly long after the catalogue moves underneath it. Feeding
    both through one implementation also means the two columns cannot differ by a definition: same
    tiers, same `quotable_level` fallback, same four states — only the donor set changes.
    """
    out = defaultdict(list)
    for locus_id, rank, neighbour_id, cosine in session.execute(
        text("""
            select n.locus_id, n.rank, n.neighbour_locus_id, n.cross_similarity
              from locus_nearest_locus n
              join locus l on l.locus_id = n.locus_id
              join pangenome p on p.pangenome_id = l.pangenome_id
             where p.catalogue_key = :catalogue and n.representation = 'ESM'
               and n.cross_similarity is not null
        """),
        {"catalogue": catalogue},
    ).all():
        out[locus_id].append((rank, neighbour_id, cosine))
    return out


def walk(
    candidates: dict,
    *,
    recipient_ids,
    donor_levels,
    recipient_levels,
    kind: AnnotationKind,
    label_of: dict | None = None,
    max_rank: int | None = None,
) -> dict:
    """Walk each unlabelled recipient node's donors by rank until one carries the vocabulary.

    ⛔ **Four states, counted separately, and they must sum to the unlabelled node total.** A node
    with no annotated donor and a node whose donor agreed are both "no disagreement" in a naive
    tally, and that is the failure `nuna/CLAUDE.md` §2 exists to forbid. *No candidate at all* is
    kept apart from *no annotated donor among its candidates* because they mean different things:
    the first is a node the neighbour structure never reached, the second is the dark set.

    ⚠ `candidates` is keyed on **locus_id**, not on node label, so the cross-species map and the
    stored within-species graph feed the same function. `label_of` only decorates the worked
    examples; nothing is measured through it.

    ⛔⛔ **`max_rank` exists because the two donor lists are DIFFERENT LENGTHS, and that alone would
    decide the comparison.** `locus_nearest_locus` stores **5** neighbours per locus (checked: ranks
    1-5, 4.77 on average, 0 duplicates) because the page shows five; the cross-species map carries
    up to `cap`, ~23 on average. A walk that may descend 23 ranks reaches an annotated donor more
    often than one that may descend 5 **whatever the embedding is doing** — so the reach is reported
    at its own depth *and* truncated to the baseline's, and the depth is named in every header.
    Reporting only the untruncated number would have flattered cross-species transfer for free.

    ⛔ The ranks within a locus must be distinct, or `sorted` silently falls back to ordering by
    donor id and the walk takes the lowest-numbered locus rather than the nearest one.
    """
    unlabelled = [(locus_id, genes) for locus_id, genes in recipient_ids.values() if locus_id not in recipient_levels]
    state = {
        "no_candidate_at_all": [0, 0],
        "no_annotated_donor": [0, 0],
        "too_remote_to_call": [0, 0],
        "labelled": [0, 0],
    }
    #: ⭐ every call the walk actually took: `(recipient_id, donor_id, tier, level, genes)`. The
    #: by-tier table and the worked examples are DERIVED from it rather than tallied beside it, so
    #: they cannot disagree — and Stage 4 stratifies this list instead of walking a second time.
    taken: list[tuple[int, int, str, int, int]] = []
    cosine_of: dict[tuple[int, int], float] = {}

    for locus_id, genes in unlabelled:
        mine = candidates.get(locus_id)
        if mine and len({rank for rank, _donor, _cosine in mine}) != len(mine):
            raise SystemExit(
                f"FATAL: locus {locus_id} has {len(mine)} candidates sharing "
                f"{len({r for r, _d, _c in mine})} distinct ranks — `sorted` would then order them by "
                "donor id, and the walk would take the lowest-numbered donor rather than the nearest"
            )
        if max_rank is not None and mine:
            mine = [candidate for candidate in mine if candidate[0] <= max_rank]
        if not mine:
            state["no_candidate_at_all"][0] += 1
            state["no_candidate_at_all"][1] += genes
            continue
        found = None
        for _rank, donor_id, cosine in sorted(mine):
            if donor_id in donor_levels:
                found = (donor_id, cosine)
                break
        if found is None:
            state["no_annotated_donor"][0] += 1
            state["no_annotated_donor"][1] += genes
            continue
        donor_id, cosine = found
        cosine_of[(locus_id, donor_id)] = cosine
        tier = tier_for(cosine)
        if tier == NOT_CALLED:
            state["too_remote_to_call"][0] += 1
            state["too_remote_to_call"][1] += genes
            continue
        level = quotable_level(kind, tier, donor_levels[donor_id])
        if level is None:
            state["no_annotated_donor"][0] += 1
            state["no_annotated_donor"][1] += genes
            continue
        state["labelled"][0] += 1
        state["labelled"][1] += genes
        taken.append((locus_id, donor_id, tier, level, genes))

    total = sum(count for count, _genes in state.values())
    if total != len(unlabelled):
        raise SystemExit(
            f"FATAL: the four states sum to {total:,} but there are {len(unlabelled):,} unlabelled "
            f"{kind.value} nodes — a node has been counted twice or not at all"
        )
    by_tier: dict[tuple[str, int], list[int]] = defaultdict(lambda: [0, 0])
    for _r, _d, tier, level, genes in taken:
        by_tier[(tier, level)][0] += 1
        by_tier[(tier, level)][1] += genes
    # ⭐ a handful of worked cases, because a table of counts is not something that can be checked by
    # eye. The within-species result became credible read on gumC and fimA by name.
    names = label_of or {}
    examples = [
        (names.get(r, str(r)), names.get(d, str(d)), tier, LEVEL_LABEL[kind][level], cosine_of[(r, d)])
        for r, d, tier, level, _g in taken[:EXAMPLES]
    ]
    return {
        "unlabelled": len(unlabelled),
        "state": state,
        "by_tier": dict(by_tier),
        "examples": examples,
        "taken": taken,
        "max_rank": max_rank,
        "depth_p50": _median_depth(candidates, max_rank),
    }


def _median_depth(candidates: dict, max_rank: int | None) -> float:
    """How many donors the walk was actually allowed to see, per locus — the comparison's own scale."""
    depths = sorted(
        len(mine) if max_rank is None else len([c for c in mine if c[0] <= max_rank]) for mine in candidates.values()
    )
    return float(depths[len(depths) // 2]) if depths else float("nan")


def report_reach(
    reach: dict, matched: dict, baseline: dict, kind: AnnotationKind, *, rule: str, recipient: str, donor: str
) -> None:
    """What the transfer would add, in nodes AND genes, at its own depth and at the baseline's.

    ⛔ **Three columns, because the middle one is the only fair one.** `locus_nearest_locus` stores
    five neighbours per locus and the cross-species map carries ~23, so the untruncated reach is
    partly a longer-list effect rather than an embedding effect. The depth each column was allowed
    is printed in its header; the two that share a depth are the comparison.

    ⚠ **Additive, never a replacement** — the two donor sets differ and a node can be reached by
    both, so the columns are not to be subtracted from one another.
    """
    for column in (reach, matched, baseline):
        if column["unlabelled"] != reach["unlabelled"]:
            raise SystemExit(
                f"FATAL: the columns describe different populations ({column['unlabelled']:,} vs "
                f"{reach['unlabelled']:,} unlabelled nodes). They must be the same nodes asked "
                "different questions, or the comparison is between two things."
            )
    print(f"\n  {kind.value} · ranked by '{rule}' — what a {donor}->{recipient} transfer would reach")
    heads = [
        (f"{donor}->{recipient}", f"all {reach['depth_p50']:.0f}", reach),
        (f"{donor}->{recipient}", f"top {matched['max_rank']}", matched),
        (f"within-{recipient}", f"all {baseline['depth_p50']:.0f}", baseline),
    ]
    print(f"     {'state':<28}" + "".join(f"{name:>16}" for name, _depth, _c in heads))
    print(f"     {'donors the walk may see':<28}" + "".join(f"{depth:>16}" for _n, depth, _c in heads))
    for name in reach["state"]:
        print(
            f"     {name.replace('_', ' '):<28}"
            + "".join(f"{column['state'][name][0]:>9,}/{column['state'][name][1]:<6,}" for _n, _d, column in heads)
        )
    print(f"     {'TOTAL unlabelled nodes':<28}{reach['unlabelled']:>16,}")
    print("     cells are nodes/genes. ⚠ ADDITIVE to the within-species gain, never a replacement;")
    print("     ⛔ and the two columns that share a depth are the only ones comparable to each other.")
    if reach["by_tier"]:
        print("\n     labelled cross-species (full depth), by the tier that earned it:")
        print(f"     {'tier':<11}{'rung':<18}{'nodes':>9}{'genes':>11}")
        for tier, _low, _high in TIERS:
            for level in sorted(LEVEL_LABEL[kind], reverse=True):
                cell = reach["by_tier"].get((tier, level))
                if cell:
                    print(f"     {tier:<11}{LEVEL_LABEL[kind][level]:<18}{cell[0]:>9,}{cell[1]:>11,}")
    if reach["examples"]:
        print("\n     worked cases — read a few by hand before trusting the totals:")
        for node, donor_node, tier, rung, cosine in reach["examples"]:
            print(f"       {recipient} {node:<10} <- {donor} {donor_node:<10} {cosine:.4f}  {tier:<11}{rung}")


# ============================================= SECTION 3 — WHY? the UniRef50 bridge (Stage 4)
#: ⛔ **THREE states, never two.** *not bridged* and *not measurable* are different claims, and
#: conflating them is the mistake `nuna/CLAUDE.md` retired the label `no_homology` for: it meant
#: *not measured*, not *nothing found*. A pair where either node carries no UniRef50 accession at
#: all has not been shown to lack a bridge — it has not been asked.
BRIDGE_STATES = ("bridged", "not_bridged", "unmeasurable")


def uniref50_sets(session: Session, catalogue: str) -> dict[int, frozenset[str]]:
    """Per locus, the set of UniRef50 accessions **its genes** carry — the gene-level bridge.

    ⛔ Gene level, not `locus_uniref_family_crosstab`, which keeps only the top 8 families per locus
    (`rank_within_locus` maxes at 7). A bridge is an *intersection*, so the single family that links
    two nodes can be a rare one in either of them — exactly what a cap drops. `capped_uniref50_sets`
    exists beside this to put a number on what the cap would have hidden, rather than an argument.
    """
    out: dict[int, set[str]] = defaultdict(set)
    for locus_id, accession in session.execute(
        text("""
            select distinct m.locus_id, f.uniref50_accession
              from gene_locus_membership m
              join locus l on l.locus_id = m.locus_id
              join pangenome p on p.pangenome_id = l.pangenome_id
              join gene_functional_annotation f
                   on f.genome_id = m.genome_id and f.flat_index = m.flat_index
             where p.catalogue_key = :catalogue and f.uniref50_accession is not null
        """),
        {"catalogue": catalogue},
    ).all():
        out[locus_id].add(accession)
    return {locus_id: frozenset(values) for locus_id, values in out.items()}


def capped_uniref50_sets(session: Session, catalogue: str) -> dict[int, frozenset[str]]:
    """The same sets as the published crosstab carries them — top 8 families per locus."""
    out: dict[int, set[str]] = defaultdict(set)
    for locus_id, accession in session.execute(
        text("""
            select c.locus_id, c.uniref50_accession
              from locus_uniref_family_crosstab c
              join locus l on l.locus_id = c.locus_id
              join pangenome p on p.pangenome_id = l.pangenome_id
             where p.catalogue_key = :catalogue
        """),
        {"catalogue": catalogue},
    ).all():
        out[locus_id].add(accession)
    return {locus_id: frozenset(values) for locus_id, values in out.items()}


def bridge_state(recipient: frozenset[str] | None, donor: frozenset[str] | None) -> str:
    """`bridged` | `not_bridged` | `unmeasurable` — see `BRIDGE_STATES`."""
    if not recipient or not donor:
        return "unmeasurable"
    return "bridged" if (recipient & donor) else "not_bridged"


def report_cap_cost(rows, *, recipient_ids, donor_ids, gene_level, capped) -> None:
    """How many pairs the 8-capped crosstab would have called differently — a number, not a claim."""
    moved: Counter[tuple[str, str]] = Counter()
    for row in rows:
        r_id, d_id = recipient_ids[row["recipient_node"]][0], donor_ids[row["donor_node"]][0]
        fine = bridge_state(gene_level["recipient"].get(r_id), gene_level["donor"].get(d_id))
        coarse = bridge_state(capped["recipient"].get(r_id), capped["donor"].get(d_id))
        moved[(fine, coarse)] += 1
    disagreeing = sum(n for (fine, coarse), n in moved.items() if fine != coarse)
    total = sum(moved.values())
    print(f"\n  what the 8-capped crosstab would have hidden, over all {total:,} shortlisted pairs:")
    print(f"     {disagreeing:,} pairs ({disagreeing / total:.2%}) get a different bridge state")
    for (fine, coarse), n in sorted(moved.items(), key=lambda kv: -kv[1]):
        flag = "   <- the cap's error" if fine != coarse else ""
        print(f"     gene-level {fine:<13} capped {coarse:<13}{n:>9,}{flag}")


def bridge_reliability(
    rows,
    *,
    recipient_ids,
    donor_ids,
    recipient_sets,
    donor_sets,
    both,
    donor_levels,
    kind: AnnotationKind,
    rule: str,
    recipient: str,
    donor: str,
) -> None:
    """⭐ Agreement by tier, SPLIT by whether a UniRef50 family bridges the pair.

    **This is the question Stage 4 exists for.** If the ladder only holds where UniRef50 already
    bridges the two nodes, cross-species ESM transfer adds nothing a sequence search would not have
    found. If it holds in the **not bridged** stratum, that is Nuna's founding claim measured —
    divergent orthologues and HGT islands that an identity threshold cannot see.

    ⚠ The chance baseline stays the **whole donor catalogue** in both strata: the comparator is still
    *"a random annotated donor node"*, which is the alternative to using this donor, and it does not
    become a different question because this particular donor happens to share a family.
    """
    split: dict[str, list] = {state: [] for state in BRIDGE_STATES}
    for row in rows:
        r_id, d_id = recipient_ids[row["recipient_node"]][0], donor_ids[row["donor_node"]][0]
        split[bridge_state(recipient_sets.get(r_id), donor_sets.get(d_id))].append((r_id, d_id, row[rule]))
    cells = {state: tally_cells(both, edges, kind, donor_pool=donor_levels) for state, edges in split.items()}
    print(f"\n  {kind.value} · ranked by '{rule}' — agreement SPLIT by the UniRef50 bridge")
    print(f"     {'tier':<11}{'rung':<18}" + "".join(f"{state:>22}" for state in BRIDGE_STATES))
    for tier, _low, _high in TIERS:
        level = QUOTED_LEVEL[kind].get(tier)
        if level is None:
            continue
        line = f"     {tier:<11}{LEVEL_LABEL[kind][level]:<18}"
        for state in BRIDGE_STATES:
            line += f"{rate(cells[state].get((tier, level))):>22}"
        print(line)
    print("     ⭐ 'not_bridged' is the interesting column: ESM reaching where 50 % identity does not.")
    print("     ⛔ 'unmeasurable' = one side carries NO UniRef50 at all — not asked, not answered.")


def bridge_reach(
    reach: dict, *, recipient_sets, donor_sets, kind: AnnotationKind, label_of: dict, recipient: str, donor: str
) -> None:
    """⭐ The labelled nodes, split by bridge state — how much of the reach sequence identity misses.

    Stratifies the walk's OWN record of what it took, so this cannot drift from the §5 totals.
    """
    rows: dict[tuple[str, str], list[int]] = defaultdict(lambda: [0, 0])
    unbridged_examples = []
    for recipient_id, donor_id, tier, _level, genes in reach["taken"]:
        state = bridge_state(recipient_sets.get(recipient_id), donor_sets.get(donor_id))
        rows[(tier, state)][0] += 1
        rows[(tier, state)][1] += genes
        if state == "not_bridged" and tier in (">= 0.99", "0.98-0.99") and len(unbridged_examples) < EXAMPLES:
            unbridged_examples.append(
                (label_of.get(recipient_id, recipient_id), label_of.get(donor_id, donor_id), tier)
            )
    print(f"\n  {kind.value} — the reach, split by the UniRef50 bridge (nodes/genes)")
    print(f"     {'tier':<11}" + "".join(f"{state:>22}" for state in BRIDGE_STATES))
    totals: dict[str, list[int]] = {state: [0, 0] for state in BRIDGE_STATES}
    for tier, _low, _high in TIERS:
        line = f"     {tier:<11}"
        for state in BRIDGE_STATES:
            nodes, genes = rows.get((tier, state), [0, 0])
            totals[state][0] += nodes
            totals[state][1] += genes
            line += f"{nodes:>12,}/{genes:<9,}"
        print(line)
    print(f"     {'TOTAL':<11}" + "".join(f"{totals[s][0]:>12,}/{totals[s][1]:<9,}" for s in BRIDGE_STATES))
    high = {
        s: [sum(rows.get((t, s), [0, 0])[i] for t in (">= 0.99", "0.98-0.99")) for i in (0, 1)] for s in BRIDGE_STATES
    }
    print(
        f"     {'>= 0.98':<11}"
        + "".join(f"{high[s][0]:>12,}/{high[s][1]:<9,}" for s in BRIDGE_STATES)
        + "   <- the tiers the ladder stands behind"
    )
    if unbridged_examples:
        print("     worked cases with NO UniRef50 bridge, at >= 0.98 — the ones worth reading by hand:")
        for node, donor_node, tier in unbridged_examples:
            print(f"       {recipient} {node:<10} <- {donor} {donor_node:<10} {tier}")


def modal_symbols(session: Session, catalogue: str) -> dict[int, str]:
    """Per locus, its modal Bakta gene symbol — the independent axis the bridge control needs."""
    return {
        locus_id: symbol
        for locus_id, symbol in session.execute(
            text("""
                select l.locus_id, a.term_value
                  from locus l join pangenome p using (pangenome_id)
                  join locus_annotation_entry a on a.locus_id = l.locus_id
                 where p.catalogue_key = :catalogue
                   and a.annotation_kind = 'GENE_SYMBOL' and a.rank_within_locus = 0
            """),
            {"catalogue": catalogue},
        ).all()
    }


def bridge_symbol_control(
    rows,
    *,
    recipient_ids,
    donor_ids,
    recipient_sets,
    donor_sets,
    symbols,
    rule: str,
    floor: float = 0.98,
) -> None:
    """⭐ Is `not_bridged` remote homology, or UniRef50 ASSIGNMENT noise? The control that asks.

    ⛔ **The distinction decides the headline.** `yjeT`/`yjeT`, `priB`/`priB` and `purT`/`purT` all read
    as *not bridged* while carrying the same Bakta symbol and sitting in different UniRef50 clusters
    (ecoli `UniRef50_P39771` against kp `UniRef50_P46927` for `purT`). Either they are genuinely below
    50 % identity — Nuna's founding claim — or UniRef50 put two true orthologues in different clusters,
    which would make the unbridged count an artifact of the bridge rather than a property of ESM.

    Bakta `GENE_SYMBOL` is a **third** axis, independent of both UniRef50 and the vocabulary being
    transferred, and it is **one-directional evidence**: symbol *agreement* is evidence the two nodes
    are the same gene; symbol *disagreement* is NOT evidence against it, because the two species'
    symbol vocabularies differ — which is why `nuna/CLAUDE.md` calls Bakta symbols the wrong
    cross-species *currency*. It is used here only in the direction it supports.

    ⚠ **The denominator is pairs where BOTH nodes are named**, reported beside the rate, because the
    unbridged set is named far less often than the bridged one — so the rate rests on a subset that
    may be the easier part of it.
    """
    tally: dict[str, list[int]] = {state: [0, 0, 0] for state in BRIDGE_STATES}
    for row in rows:
        if int(row[f"rank_{rule}"]) != 1 or row[rule] < floor:
            continue
        r_id, d_id = recipient_ids[row["recipient_node"]][0], donor_ids[row["donor_node"]][0]
        state = bridge_state(recipient_sets.get(r_id), donor_sets.get(d_id))
        recipient_symbol, donor_symbol = symbols["recipient"].get(r_id), symbols["donor"].get(d_id)
        tally[state][0] += 1
        if recipient_symbol and donor_symbol:
            tally[state][1] += 1
            tally[state][2] += recipient_symbol == donor_symbol
    print(f"\n  CONTROL — do the pairs share a Bakta SYMBOL? (rank-1 pairs at '{rule}' >= {floor})")
    print(f"     {'bridge state':<16}{'pairs':>9}{'both named':>12}{'same symbol':>13}{'agreement':>12}")
    for state in BRIDGE_STATES:
        pairs, named, same = tally[state]
        print(f"     {state:<16}{pairs:>9,}{named:>12,}{same:>13,}{(f'{same / named:.1%}' if named else '—'):>12}")
    print("     ⛔ ONE-DIRECTIONAL: agreement is evidence of orthology; disagreement is not evidence")
    print("        against it, because the two species' symbol vocabularies differ.")
    print("     ⚠ A high unbridged rate means UniRef50 SEPARATED true orthologues — it does not say")
    print("        whether they are below 50 % identity. Only measured identity answers that.")


# ======================= SECTION 4 — UniRef50 FIRST, then ESM: which arm infers what (David's ask)
#: A shared family must cover at least this share of the recipient node's genes to count for the
#: strict variant. The permissive variant takes ANY shared family — which is the right default for
#: David's point that *"not all the syntelogues will match the same uniref 50"*, and the strict one is
#: the sensitivity check against a bridge resting on one spurious rare family.
STRICT_BRIDGE_SHARE = 0.5


def uniref50_gene_counts(session: Session, catalogue: str) -> dict[int, dict[str, int]]:
    """Per locus, how many of its genes carry each UniRef50 accession.

    ⭐ Gene **counts**, not just the set, because David's point needs them: a bridge can rest on a
    minority of a node's genes, and whether that is enough is a decision the numbers should inform
    rather than one buried in a set intersection.
    """
    out: dict[int, dict[str, int]] = defaultdict(dict)
    for locus_id, accession, genes in session.execute(
        text("""
            select m.locus_id, f.uniref50_accession, count(*)
              from gene_locus_membership m
              join locus l on l.locus_id = m.locus_id
              join pangenome p on p.pangenome_id = l.pangenome_id
              join gene_functional_annotation f
                   on f.genome_id = m.genome_id and f.flat_index = m.flat_index
             where p.catalogue_key = :catalogue and f.uniref50_accession is not null
             group by 1, 2
        """),
        {"catalogue": catalogue},
    ).all():
        out[locus_id][accession] = genes
    return dict(out)


def uniref_first_then_esm(
    *,
    recipient_ids,
    donor_ids,
    recipient_counts,
    donor_counts,
    donor_levels,
    recipient_levels,
    esm_candidates,
    kind: AnnotationKind,
    label_of: dict,
    recipient: str,
    donor: str,
    esm_floor: float = 0.98,
) -> None:
    """⭐ Two arms in series: the UniRef50 bridge first (NO embedding), then ESM for what it missed.

    **David, 2026-10-03:** *"if we are inferring by uniref first (at any esm) and then by esm, the
    question is how many genes get inferred by each."*

    ⛔ **Arm A uses no embedding at all and is not restricted to the ESM shortlist.** A UniRef50-first
    transfer is a catalogue-wide join on the accession — it does not need a neighbour map, so limiting
    it to the ~21 shortlisted candidates would hand ESM a reach that belongs to the bridge. This arm is
    therefore the **baseline cross-species ESM has to beat**, and it was not measured until now.

    ⭐ **Two gene counts per arm, because they answer different questions** (David: *"not all the
    syntelogues will match the same uniref 50. BUT given the syntelogue node itself is very strongly
    uniform … we can then annotate the whole node from it"*): the genes **actually carrying** a shared
    family, and the **whole node** that a within-node imputation then covers. The second is the one that
    moves the coverage dial, and it rests on the measured 99.5-99.9 % within-node unanimity
    (`function_inference.md` §2), not on an assumption.
    """
    # an inverted index over the donor catalogue: family -> the donor loci that carry it
    holders: dict[str, set[int]] = defaultdict(set)
    for locus_id, families in donor_counts.items():
        for accession in families:
            holders[accession].add(locus_id)

    unlabelled = [(locus_id, genes) for locus_id, genes in recipient_ids.values() if locus_id not in recipient_levels]
    arms = {
        "A: UniRef50 bridge (any ESM)": [0, 0, 0],
        f"A-strict: bridge >= {STRICT_BRIDGE_SHARE:.0%} of genes": [0, 0, 0],
        f"B: ESM >= {esm_floor}, no bridge available": [0, 0, 0],
        "neither": [0, 0, 0],
    }
    examples = []
    for locus_id, node_genes in unlabelled:
        mine = recipient_counts.get(locus_id, {})
        # every donor locus sharing a family with this node, and how many of MY genes that family covers
        reach: dict[int, int] = defaultdict(int)
        for accession, my_genes in mine.items():
            for donor_id in holders.get(accession, ()):
                reach[donor_id] = max(reach[donor_id], my_genes)
        annotated = {d: shared for d, shared in reach.items() if d in donor_levels}
        if annotated:
            best = max(annotated, key=lambda d: (annotated[d], -d))
            shared_genes = annotated[best]
            arms["A: UniRef50 bridge (any ESM)"][0] += 1
            arms["A: UniRef50 bridge (any ESM)"][1] += shared_genes
            arms["A: UniRef50 bridge (any ESM)"][2] += node_genes
            key = f"A-strict: bridge >= {STRICT_BRIDGE_SHARE:.0%} of genes"
            if node_genes and shared_genes / node_genes >= STRICT_BRIDGE_SHARE:
                arms[key][0] += 1
                arms[key][1] += shared_genes
                arms[key][2] += node_genes
            if len(examples) < 4:
                examples.append((label_of.get(locus_id, locus_id), label_of.get(best, best), shared_genes, node_genes))
            continue
        # Arm B: the bridge offered nothing, so ESM is the only route left
        best_esm = None
        for _rank, donor_id, cosine in sorted(esm_candidates.get(locus_id, [])):
            if donor_id in donor_levels and cosine >= esm_floor:
                best_esm = (donor_id, cosine)
                break
        if best_esm:
            arms[f"B: ESM >= {esm_floor}, no bridge available"][0] += 1
            arms[f"B: ESM >= {esm_floor}, no bridge available"][2] += node_genes
        else:
            arms["neither"][0] += 1
            arms["neither"][2] += node_genes

    # ⛔ A reach without its accuracy is exactly what this whole study refuses to report. Arm A's
    # donors are chosen catalogue-wide with NO similarity requirement, so the "bridged" agreement
    # measured over the ESM shortlist (97.8-99.7 %) does not transfer to it — the shortlist excludes
    # the remote pairs an unrestricted join admits. Measured here on its own terms: recipient nodes
    # that DO carry the vocabulary, against the call their best bridged donor would have given,
    # per rung and pooled over all similarities.
    arm_a: dict[int, list[int]] = defaultdict(lambda: [0, 0])
    for locus_id, mine_levels in recipient_levels.items():
        reach: dict[int, int] = defaultdict(int)
        for accession, my_genes in recipient_counts.get(locus_id, {}).items():
            for donor_id in holders.get(accession, ()):
                reach[donor_id] = max(reach[donor_id], my_genes)
        annotated = {d: shared for d, shared in reach.items() if d in donor_levels}
        if not annotated:
            continue
        best = max(annotated, key=lambda d: (annotated[d], -d))
        theirs = donor_levels[best]
        for level in sorted(LEVEL_LABEL[kind], reverse=True):
            if mine_levels[level] and theirs[level]:
                arm_a[level][0] += 1
                arm_a[level][1] += bool(mine_levels[level] & theirs[level])

    # ⛔ And Arm B on ITS OWN population, for the same reason. Arm B is *nodes with no bridged
    # annotated donor anywhere*, which is not the same set as *nodes whose ESM-chosen donor happened
    # to be unbridged* — so the stratified table's `not_bridged` column is close to it but is not it.
    # Quoting one arm at the other's rate is the error both of these exist to avoid.
    arm_b: dict[int, list[int]] = defaultdict(lambda: [0, 0])
    for locus_id, mine_levels in recipient_levels.items():
        reach = {
            donor_id for accession in recipient_counts.get(locus_id, {}) for donor_id in holders.get(accession, ())
        }
        if reach & set(donor_levels):
            continue  # a bridge was available, so this node belongs to Arm A
        for _rank, donor_id, cosine in sorted(esm_candidates.get(locus_id, [])):
            if donor_id in donor_levels and cosine >= esm_floor:
                theirs = donor_levels[donor_id]
                for level in sorted(LEVEL_LABEL[kind], reverse=True):
                    if mine_levels[level] and theirs[level]:
                        arm_b[level][0] += 1
                        arm_b[level][1] += bool(mine_levels[level] & theirs[level])
                break

    print(f"\n  {kind.value} — each arm's OWN accuracy, measured on its own population")
    print(f"     {'rung':<20}{'A: bridge, any sim':>22}{'B: ESM >= ' + str(esm_floor) + ', no bridge':>30}")
    for level in sorted(LEVEL_LABEL[kind], reverse=True):
        cells = []
        for arm in (arm_a, arm_b):
            pairs, agreeing = arm[level]
            cells.append(f"{agreeing / pairs:.1%} n={pairs:,}" if pairs >= MIN_PAIRS else f"n={pairs:,}")
        print(f"     {LEVEL_LABEL[kind][level]:<20}{cells[0]:>22}{cells[1]:>30}")
    print("     ⚠ Column A is NOT the 'bridged' column of the stratified table — that one is restricted")
    print("        to pairs the ESM shortlist surfaced, and this arm deliberately is not. Column B is")
    print("        NOT the 'not_bridged' column either: that is 'the chosen donor was unbridged', this")
    print("        is 'no bridged donor existed anywhere'. Each arm is quoted at its own rate.")

    print(f"\n  {kind.value} — UniRef50 FIRST, then ESM: what each arm infers")
    print(f"     {'arm':<42}{'nodes':>9}{'bridged genes':>15}{'whole-node genes':>18}")
    for name, (nodes, bridged_genes, node_genes) in arms.items():
        shown = f"{bridged_genes:,}" if bridged_genes else "—"
        print(f"     {name:<42}{nodes:>9,}{shown:>15}{node_genes:>18,}")
    total = sum(arm[0] for name, arm in arms.items() if not name.startswith("A-strict"))
    print(f"     {'TOTAL unlabelled nodes':<42}{total:>9,}")
    if total != len(unlabelled):
        raise SystemExit(f"FATAL: the arms sum to {total:,} but there are {len(unlabelled):,} unlabelled nodes")
    print("     ⛔ Arm A uses NO embedding and is NOT limited to the ESM shortlist — it is the baseline")
    print("        cross-species ESM has to beat. Arm B is ESM's MARGINAL contribution over it.")
    print("     ⭐ 'bridged genes' are the genes actually carrying a shared family; 'whole-node genes'")
    print("        is what a within-node imputation then covers, at the measured 99.5-99.9 % unanimity.")
    if examples:
        print("     arm-A cases (shared genes / node genes):")
        for node, donor_node, shared, node_genes in examples:
            print(f"       {recipient} {node:<8} <- {donor} {donor_node:<8} {shared}/{node_genes}")


def main() -> None:
    """Section 1 (reliable?) then Section 2 (worthwhile?), for every vocabulary and every rule."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("map", help="the TSV nuna.tl.probe.cross_species_neighbours wrote")
    parser.add_argument("--donor", default="ecoli", choices=sorted(CATALOGUE))
    parser.add_argument("--recipient", default="kp", choices=sorted(CATALOGUE))
    parser.add_argument("--rules", nargs="+", default=list(RULES), choices=list(RULES))
    arguments = parser.parse_args()
    if arguments.donor == arguments.recipient:
        sys.exit(f"--donor and --recipient are both '{arguments.donor}'; this measures nothing")

    url = os.environ.get("BACATLAS_DATABASE_URL")
    if not url:
        sys.exit("BACATLAS_DATABASE_URL is not set — see backend/README.md")

    rows, provenance = read_map(arguments.map, donor=arguments.donor, recipient=arguments.recipient)
    donor_catalogue, recipient_catalogue = CATALOGUE[arguments.donor], CATALOGUE[arguments.recipient]
    print(f"cross-species transfer: {arguments.recipient} asks, {arguments.donor} answers")
    print(f"  map: {arguments.map}")
    print(
        f"  {len(rows):,} shortlisted pairs | "
        f"{len({row['recipient_node'] for row in rows}):,} {arguments.recipient} nodes | "
        f"{len({row['donor_node'] for row in rows}):,} {arguments.donor} nodes"
    )
    print(
        f"  rep={provenance.get('rep')} k={provenance.get('k')} cap={provenance.get('cap')} "
        f"rank_by={provenance.get('rank_by')}"
    )
    # ⛔ The map's OWN coverage, printed before any transfer number is reported from it. A recall
    # figure belongs to the rule it was measured under, so each is named with its rule.
    for audit in provenance.get("recall") or []:
        print(
            f"  shortlist recall [{audit['rule']} q={audit['q']}]: top1={audit['top1']:.4f} "
            f"top1_in_list={audit['top1_in_topk']:.4f} list_recall={audit['topk_recall']:.4f} "
            f"(n={audit['n_probe_loci']})"
        )
    if not provenance.get("recall"):
        print("  ⚠ the provenance carries NO recall audit — the shortlist's coverage is unmeasured")
    print(f"  a cell below {MIN_PAIRS} pairs reports its n and no rate")

    with Session(create_engine(url)) as database:
        donor_ids = locus_index(database, donor_catalogue)
        recipient_ids = locus_index(database, recipient_catalogue)
        unknown = ({row["donor_node"] for row in rows} - set(donor_ids)) | (
            {row["recipient_node"] for row in rows} - set(recipient_ids)
        )
        if unknown:
            sys.exit(
                f"FATAL: {len(unknown)} node labels in the map are in neither catalogue "
                f"(e.g. {sorted(unknown)[:5]}) — the map and the database describe different models"
            )
        within = within_species_edges(database, recipient_catalogue)
        within_graph = candidates_from_database(database, recipient_catalogue)
        # ⛔ Merge on locus_id, NOT on node label. The two catalogues share the label namespace —
        # both run '0', '1', '10', … — so `{**donor_ids, **recipient_ids}` lets the recipient's
        # labels overwrite the donor's and every donor locus_id falls out of the map. The symptom was
        # the donor column of the worked cases printing raw locus ids (`ecoli 304133`) instead of node
        # labels, which looked like a formatting choice rather than a lost lookup. locus_id IS
        # globally unique (checked: ecoli 302532-320062, kp 320063-335732, 0 shared), so keying on it
        # is safe where keying on the label is not.
        label_of = {locus_id: label for label, (locus_id, _g) in donor_ids.items()}
        label_of.update({locus_id: label for label, (locus_id, _g) in recipient_ids.items()})
        if len(label_of) != len(donor_ids) + len(recipient_ids):
            sys.exit(
                f"FATAL: {len(donor_ids) + len(recipient_ids) - len(label_of)} locus_ids are shared "
                "between the two catalogues — the merged claims dict and this label map both collide"
            )
        # ⛔ measured, not assumed: the cross-species reach is truncated to this before it is
        # compared, because a longer candidate list reaches an annotated donor more often whatever
        # the embedding is doing. `locus_nearest_locus` stores five because the page shows five.
        baseline_depth = max((len(v) for v in within_graph.values()), default=0)
        # Stage 4 — the UniRef50 bridge, at GENE level. Loaded once; it does not vary by vocabulary.
        gene_level = {
            "donor": uniref50_sets(database, donor_catalogue),
            "recipient": uniref50_sets(database, recipient_catalogue),
        }
        capped = {
            "donor": capped_uniref50_sets(database, donor_catalogue),
            "recipient": capped_uniref50_sets(database, recipient_catalogue),
        }
        symbols = {
            "donor": modal_symbols(database, donor_catalogue),
            "recipient": modal_symbols(database, recipient_catalogue),
        }
        counts = {
            "donor": uniref50_gene_counts(database, donor_catalogue),
            "recipient": uniref50_gene_counts(database, recipient_catalogue),
        }
        print(
            f"  within-{arguments.recipient} baseline: {len(within):,} stored ESM neighbour edges over "
            f"{len(within_graph):,} loci, at most {baseline_depth} per locus — recomputed from the "
            "database on every run, and the depth the cross-species reach is truncated to"
        )

        for side in ("donor", "recipient"):
            catalogue = donor_catalogue if side == "donor" else recipient_catalogue
            with_any = len(gene_level[side])
            print(
                f"  UniRef50 ({side}, {catalogue}): {with_any:,} loci carry at least one accession at "
                f"gene level, {len(capped[side]):,} in the 8-capped crosstab"
            )
        report_cap_cost(rows, recipient_ids=recipient_ids, donor_ids=donor_ids, gene_level=gene_level, capped=capped)
        for rule in arguments.rules:
            bridge_symbol_control(
                rows,
                recipient_ids=recipient_ids,
                donor_ids=donor_ids,
                recipient_sets=gene_level["recipient"],
                donor_sets=gene_level["donor"],
                symbols=symbols,
                rule=rule,
            )

        for kind in KINDS:
            donor_levels, donor_loci = claims(database, donor_catalogue, kind)
            recipient_levels, recipient_loci = claims(database, recipient_catalogue, kind)
            both = {**donor_levels, **recipient_levels}
            print(
                f"\n{'=' * 108}\n{kind.value}: "
                f"{len(donor_levels):,}/{donor_loci:,} {arguments.donor} nodes carry it, "
                f"{len(recipient_levels):,}/{recipient_loci:,} {arguments.recipient} nodes do"
            )
            # ⭐ recomputed, never transcribed — the column the cross-species one is judged against
            within_cells = tally_cells(recipient_levels, within, kind)

            # the recipient's OWN graph, through the SAME walk — the column the reach is read against
            baseline = walk(
                within_graph,
                recipient_ids=recipient_ids,
                donor_levels=recipient_levels,
                recipient_levels=recipient_levels,
                kind=kind,
                label_of=label_of,
            )
            for rule in arguments.rules:
                edges = [
                    (recipient_ids[row["recipient_node"]][0], donor_ids[row["donor_node"]][0], row[rule])
                    for row in rows
                ]
                # ⚠ the donor's claims alone are the null; see the module docstring
                cross_cells = tally_cells(both, edges, kind, donor_pool=donor_levels)
                reliability(
                    cross_cells, within_cells, kind, rule=rule, donor=arguments.donor, recipient=arguments.recipient
                )
                print("\n     DOES THE LADDER TRANSFER? the quoted rung only, as a difference:")
                for note in transfer_gate(cross_cells, within_cells, kind, recipient=arguments.recipient):
                    print(note)
                mapped = candidates_from_map(rows, rule=rule, recipient_ids=recipient_ids, donor_ids=donor_ids)
                shared = dict(
                    recipient_ids=recipient_ids,
                    donor_levels=donor_levels,
                    recipient_levels=recipient_levels,
                    kind=kind,
                    label_of=label_of,
                )
                reach = walk(mapped, **shared)
                report_reach(
                    reach,
                    # ⛔ truncated to the depth the stored within-species graph actually has, or the
                    # comparison is a list-length effect wearing an embedding's clothes
                    walk(mapped, max_rank=baseline_depth, **shared),
                    baseline,
                    kind,
                    rule=rule,
                    recipient=arguments.recipient,
                    donor=arguments.donor,
                )
                # ⭐ Stage 4 — WHY. Runs last and is read last: it explains a reach that §5 has
                # already established, and its numbers were kept out of §3 and §5 for that reason.
                bridge_reliability(
                    rows,
                    recipient_ids=recipient_ids,
                    donor_ids=donor_ids,
                    recipient_sets=gene_level["recipient"],
                    donor_sets=gene_level["donor"],
                    both=both,
                    donor_levels=donor_levels,
                    kind=kind,
                    rule=rule,
                    recipient=arguments.recipient,
                    donor=arguments.donor,
                )
                bridge_reach(
                    reach,
                    recipient_sets=gene_level["recipient"],
                    donor_sets=gene_level["donor"],
                    kind=kind,
                    label_of=label_of,
                    recipient=arguments.recipient,
                    donor=arguments.donor,
                )
                uniref_first_then_esm(
                    recipient_ids=recipient_ids,
                    donor_ids=donor_ids,
                    recipient_counts=counts["recipient"],
                    donor_counts=counts["donor"],
                    donor_levels=donor_levels,
                    recipient_levels=recipient_levels,
                    esm_candidates=mapped,
                    kind=kind,
                    label_of=label_of,
                    recipient=arguments.recipient,
                    donor=arguments.donor,
                )


if __name__ == "__main__":
    main()
