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
from collections import defaultdict
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
    by_tier: dict[tuple[str, int], list[int]] = defaultdict(lambda: [0, 0])
    examples: list[tuple] = []

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
        by_tier[(tier, level)][0] += 1
        by_tier[(tier, level)][1] += genes
        # ⭐ a handful of worked cases, because a table of counts is not something that can be
        # checked by eye. The within-species result became credible read on gumC and fimA by name.
        if len(examples) < EXAMPLES:
            names = label_of or {}
            examples.append(
                (
                    names.get(locus_id, str(locus_id)),
                    names.get(donor_id, str(donor_id)),
                    tier,
                    LEVEL_LABEL[kind][level],
                    cosine,
                )
            )

    total = sum(count for count, _genes in state.values())
    if total != len(unlabelled):
        raise SystemExit(
            f"FATAL: the four states sum to {total:,} but there are {len(unlabelled):,} unlabelled "
            f"{kind.value} nodes — a node has been counted twice or not at all"
        )
    return {
        "unlabelled": len(unlabelled),
        "state": state,
        "by_tier": dict(by_tier),
        "examples": examples,
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
        label_of = {locus_id: label for label, (locus_id, _g) in {**donor_ids, **recipient_ids}.items()}
        # ⛔ measured, not assumed: the cross-species reach is truncated to this before it is
        # compared, because a longer candidate list reaches an annotated donor more often whatever
        # the embedding is doing. `locus_nearest_locus` stores five because the page shows five.
        baseline_depth = max((len(v) for v in within_graph.values()), default=0)
        print(
            f"  within-{arguments.recipient} baseline: {len(within):,} stored ESM neighbour edges over "
            f"{len(within_graph):,} loci, at most {baseline_depth} per locus — recomputed from the "
            "database on every run, and the depth the cross-species reach is truncated to"
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
                report_reach(
                    walk(mapped, **shared),
                    # ⛔ truncated to the depth the stored within-species graph actually has, or the
                    # comparison is a list-length effect wearing an embedding's clothes
                    walk(mapped, max_rank=baseline_depth, **shared),
                    baseline,
                    kind,
                    rule=rule,
                    recipient=arguments.recipient,
                    donor=arguments.donor,
                )


if __name__ == "__main__":
    main()
