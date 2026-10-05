"""Every number in `nuna/docs/model_evaluation/function_inference.md`, from the live database.

Run from `backend/` with `BACATLAS_DATABASE_URL` set and both `*-nuna4` catalogues loaded:

    python scripts/measure_cog_function_inference.py

⭐ **The ladder, the level folding and the agreement measurement are IMPORTED from
`instruments/annotation_transfer.py` — the same code the API serves them from.** A second copy here
is exactly how a page and the document describing it come to disagree; this script adds only the
sections the page does not need (within-node consistency, what the two mechanisms assign in nodes and
genes, the neighbour-graph assortativity, and the Pfam cross-tab).
"""

from __future__ import annotations

import os
import sys
from collections import Counter, defaultdict
from statistics import median

from sqlalchemy import create_engine, text
from sqlalchemy.orm import Session

from bacatlas_backend.instruments.annotation_transfer import (
    ANNOTATED_COUNT_PREDICATE,
    LEVEL_LABEL,
    MIN_PAIRS,
    NOT_CALLED,
    QUOTED_LEVEL,
    TIERS,
    cog_levels,
    ec_levels,
    go_rung,
    quotable_level,
    single_rung,
    tally_cells,
    tier_for,
)
from bacatlas_backend.instruments.within_node_propagation import MIN_NODES, compute_propagation
from bacatlas_backend.models.enumerations import AnnotationKind
from bacatlas_backend.services.reference_name_service import gene_ontology_fold

CATALOGUES = ("ecoli-nuna4", "kp-nuna4")
#: ⭐ All four vocabularies, not the two this file was named for. GO and KEGG were measured for the
#: first time on 2026-10-03 after David asked why they were missing; GO turns out to cover MORE genes
#: than COG in *E. coli* (351,559 against 313,617), so leaving it out understated the model's reach.
KINDS = (
    AnnotationKind.COG_ORTHOGROUP,
    AnnotationKind.GENE_ONTOLOGY_SLIM,
    AnnotationKind.KEGG_ORTHOLOGY,
    AnnotationKind.EC_NUMBER,
)
BARS = (0.95, 0.90, 0.80)
#: ⛔ The annotated-member count is COUNTED from the gene rows rather than read from the matching
#: `locus.*_annotated_member_count` column, because GO has no such column — it has three, one per
#: namespace, and a gene annotated in two of them would be counted twice by their sum. Verified
#: against the columns that do exist: 0 mismatching nodes out of 17,531 and 15,670 on COG, EC and
#: KEGG, so the two routes agree wherever both are available.
#: ⛔ **Imported, not restated.** This was a fifth hand-written copy of the same four predicates;
#: the instrument derives them from `GENE_VALUE_COLUMN`, so the column and the predicate are one
#: fact. A local copy is how a script and the API come to filter on different things.
GENE_PREDICATE = ANNOTATED_COUNT_PREDICATE


def load(session: Session, catalogue: str, kind: AnnotationKind):
    """Per node: member and annotated gene counts, prevalence band, folded top call, neighbours."""
    rows = session.execute(
        text("""
            select l.locus_id, l.member_gene_count, l.prevalence_band,
                   a.term_value, l.modal_cog_categories
              from locus l join pangenome p using (pangenome_id)
              left join locus_annotation_entry a
                     on a.locus_id = l.locus_id and a.annotation_kind = :kind
                    and a.rank_within_locus = 0
             where p.catalogue_key = :catalogue
        """),
        {"kind": kind.name, "catalogue": catalogue},
    ).all()
    # ⛔ GROUPED, because rank 0 is not unique per locus: `top_go_slim` ranks within
    # (locus, NAMESPACE), so a GO-bearing locus has three rank-0 rows and indexing by locus in the
    # loop would keep whichever namespace arrived last.
    stated = defaultdict(list)
    members, bands, categories_by_locus = {}, {}, {}
    for locus_id, n_members, band, term, categories in rows:
        members[locus_id] = (n_members, 0)
        bands[locus_id] = band
        categories_by_locus[locus_id] = categories
        if term is not None:
            stated[locus_id].append(term)
    levels = {}
    for locus_id in members:
        said = stated.get(locus_id, [])
        if kind is AnnotationKind.EC_NUMBER:
            folded = ec_levels(said[0]) if said else None
        elif kind is AnnotationKind.COG_ORTHOGROUP:
            folded = cog_levels(said[0] if said else None, categories_by_locus[locus_id])
        elif kind is AnnotationKind.GENE_ONTOLOGY_SLIM:
            # ⛔ The SAME closed claim the service measures. A script that compared GO's slim classes
            # unclosed would print numbers the page could not reproduce — and this script is where
            # `function_inference.md` §5 gets its figures, so the two would diverge in the document.
            folded = go_rung(said, gene_ontology_fold(session))
        else:
            folded = single_rung(said)
        if folded and any(folded.values()):
            levels[locus_id] = folded

    edges = session.execute(
        text("""
            select n.locus_id, n.rank, n.neighbour_locus_id, n.cross_similarity
              from locus_nearest_locus n
              join locus l on l.locus_id = n.locus_id
              join pangenome p on p.pangenome_id = l.pangenome_id
             where p.catalogue_key = :catalogue and n.representation = 'ESM'
               and n.cross_similarity is not null
             order by n.locus_id, n.rank
        """),
        {"catalogue": catalogue},
    ).all()
    neighbours = defaultdict(list)
    for locus_id, rank, neighbour_id, cosine in edges:
        neighbours[locus_id].append((rank, neighbour_id, cosine))

    gene_calls: dict[int, int] = dict(
        session.execute(
            text(f"""
            select m.locus_id, count(*)
              from gene_locus_membership m
              join locus l on l.locus_id = m.locus_id
              join pangenome p on p.pangenome_id = l.pangenome_id
              join gene_functional_annotation f
                   on f.genome_id = m.genome_id and f.flat_index = m.flat_index
             where p.catalogue_key = :catalogue and {GENE_PREDICATE[kind]}
             group by 1
        """),
            {"catalogue": catalogue},
        ).all()
    )
    #: now that the per-gene count is in hand, it IS the annotated-member count for every axis
    members = {lid: (n, gene_calls.get(lid, 0)) for lid, (n, _) in members.items()}
    return members, bands, levels, neighbours, gene_calls


def within_node(session: Session, pangenome_id: int, kind: AnnotationKind) -> None:
    """Section 1 — among a node's annotated genes, do they agree, and how often per band?

    ⛔⛔ **This section was WRONG for GO and KEGG until 2026-10-04, and silently.** It read its
    values from a statement that selected `f.cog_id` and `f.ec_numbers` only, so the GO and KEGG
    rows measured **COG** agreement over a GO- or KEGG-filtered subset of genes. The numbers looked
    entirely reasonable — around 99.5 %, which is what COG agreement is — which is why nothing
    caught it: a wrong measurement of the right shape reads as a right one.

    ⭐ It now calls `instruments/within_node_propagation.compute_propagation`, the same code the
    Function tab serves its rate from, so there is ONE implementation and it cannot drift from the
    page again. Two things changed with it, both deliberate:

    * **per prevalence band, not pooled.** Nodes where unanimity can be checked have a median of 100
      genes and nodes resting on one annotated gene a median of 1, so a single rate quotes the
      behaviour of core nodes at a singleton. `rare` has zero checkable nodes and reports no rate.
    * **at the deepest rung all the genes state, not one rung at a time.** That is the test
      `tally_cells` applies to a pair, generalised to a node's N genes; the per-rung ladder it
      replaced could call `1.1.1.1` and `1.1.1.2` agreement by dropping to L3.
    """
    print("\n  1. WITHIN-NODE CONSISTENCY — do a node's annotated genes agree, per prevalence band?")
    print(f"     {'band':<12}{'one gene':>10}{'checkable':>11}{'unanimous':>11}{'rate':>9}{'95% CI':>16}")
    rates = compute_propagation(
        session, pangenome_id=pangenome_id, annotation_kind=kind, fold=gene_ontology_fold(session)
    )
    for band in ("CORE", "SOFT_CORE", "SHELL", "CLOUD", "RARE"):
        measured = rates.get(band)
        if measured is None:
            continue
        interval = measured.interval
        # ⛔ A dash, never a 0.0: below the floor the band has a COUNT and no rate.
        print(
            f"     {band.lower().replace('_', ' '):<12}{measured.one_gene:>10,}{measured.checkable:>11,}"
            f"{measured.unanimous:>11,}"
            f"{(f'{100 * measured.rate:.2f}%' if measured.rate is not None else '—'):>9}"
            f"{(f'{100 * interval[0]:.2f}-{100 * interval[1]:.2f}%' if interval else f'< {MIN_NODES} nodes'):>16}"
        )


def assign(members, levels, neighbours, gene_calls, kind: AnnotationKind):
    """Internal inference over nodes that HAVE a call; the neighbour walk over those that do not.

    The two populations are disjoint and together they exhaust the unannotated genes.
    """
    internal = {"checkable": [0, 0], "one gene": [0, 0]}
    transfer: dict = defaultdict(lambda: [0, 0])
    supplier_rank: Counter[int] = Counter()
    sizes = defaultdict(list)
    uncalled = [0, 0]
    unreachable = [0, 0]
    assigned: dict[int, tuple[str, int]] = {}
    for locus_id, (n_members, n_annotated) in members.items():
        if locus_id in levels:
            bucket = "one gene" if gene_calls.get(locus_id, 0) <= 1 else "checkable"
            internal[bucket][0] += 1
            internal[bucket][1] += max(0, n_members - n_annotated)
            sizes["has its own call"].append(n_members)
            continue
        hit = next(((r, nb, c) for r, nb, c in neighbours.get(locus_id, []) if nb in levels), None)
        if hit is None:
            unreachable[0] += 1
            unreachable[1] += n_members
            sizes["no annotated neighbour in 5"].append(n_members)
            continue
        rank, neighbour_id, cosine = hit
        tier = tier_for(cosine)
        level = quotable_level(kind, tier, levels[neighbour_id]) if tier != NOT_CALLED else None
        if level is None:
            uncalled[0] += 1
            uncalled[1] += n_members
            sizes["neighbour too remote"].append(n_members)
            continue
        transfer[(tier, level)][0] += 1
        transfer[(tier, level)][1] += n_members
        supplier_rank[rank] += 1
        assigned[locus_id] = (tier, level)
        sizes["assigned from a neighbour"].append(n_members)
    return internal, transfer, supplier_rank, sizes, uncalled, unreachable, assigned


def report(session: Session, catalogue: str, kind: AnnotationKind) -> dict:
    """Print sections 1-7 for one catalogue and one vocabulary."""
    pangenome_id = session.execute(
        text("select pangenome_id from pangenome where catalogue_key = :catalogue"),
        {"catalogue": catalogue},
    ).scalar_one()
    members, bands, levels, neighbours, gene_calls = load(session, catalogue, kind)
    cells = tally_cells(levels, [(lid, nb, c) for lid, rows in neighbours.items() for _r, nb, c in rows], kind)
    internal, transfer, supplier_rank, sizes, uncalled, unreachable, assigned = assign(
        members, levels, neighbours, gene_calls, kind
    )

    print(f"\n{'=' * 100}\n{catalogue} · {kind.value}")
    within_node(session, pangenome_id, kind)

    print("\n  2. THE LADDER — nodes with no call of their own, labelled from a neighbour")
    print(
        f"     {'tier':<12}{'quoted':<20}{'nodes':>7}{'genes':>9}{'agree':>8}{'95% CI':>10}"
        f"{'n':>7}{'chance':>8}{'lift':>8}"
    )
    for name, _, _ in TIERS:
        for level in sorted(LEVEL_LABEL[kind], reverse=True):
            nodes, genes = transfer[(name, level)]
            if not nodes:
                continue
            cell = cells.get((name, level))
            fallback = " ←fallback" if level != QUOTED_LEVEL[kind][name] else ""
            label = LEVEL_LABEL[kind][level] + fallback
            if cell is None or cell.agreement is None:
                print(
                    f"     {name:<12}{label:<20}{nodes:>7}{genes:>9,}"
                    f"{'too few pairs':>25}{cell.pairs if cell else 0:>7}"
                )
                continue
            low, high = cell.interval
            print(
                f"     {name:<12}{label:<20}{nodes:>7}{genes:>9,}{100 * cell.agreement:7.1f}%"
                f"{f'{100 * low:.0f}-{100 * high:.0f}':>10}{cell.pairs:>7}"
                f"{100 * cell.chance:7.2f}%{cell.lift:7.1f}x"
            )
    print(f"     {NOT_CALLED + ' not called':<32}{uncalled[0]:>7}{uncalled[1]:>9,}")
    print(f"     {'no annotated neighbour in 5':<32}{unreachable[0]:>7}{unreachable[1]:>9,}")
    print("     supplying rank: " + " · ".join(f"{r}: {supplier_rank[r]:,}" for r in sorted(supplier_rank)))

    print("\n  3. INTERNAL NODE INFERENCE — the node's modal call applied to its unannotated genes")
    print(
        f"     {'>= 2 annotated genes (checkable)':<38}{internal['checkable'][0]:>7} nodes"
        f"{internal['checkable'][1]:>9,} genes"
    )
    print(
        f"     {'exactly 1 annotated gene (⛔ no check)':<38}{internal['one gene'][0]:>7} nodes"
        f"{internal['one gene'][1]:>9,} genes"
    )

    print("\n  4. HOW BIG ARE THE NODES? — why a gene-weighted percentage barely moves")
    print(f"     {'group':<30}{'nodes':>8}{'genes':>10}{'mean':>8}{'median':>8}")
    for name, values in sorted(sizes.items(), key=lambda item: -sum(item[1])):
        print(
            f"     {name:<30}{len(values):>8,}{sum(values):>10,}{sum(values) / len(values):>8.1f}{median(values):>8.0f}"
        )

    total_genes = sum(n for n, _ in members.values())
    own_genes = sum(a for _, a in members.values())
    shortfall = total_genes - own_genes
    print(
        f"\n  5. COVERAGE — {own_genes:,} of {total_genes:,} genes already named "
        f"({100 * own_genes / total_genes:.1f} %); {len(levels):,} of {len(members):,} nodes "
        f"({100 * len(levels) / len(members):.1f} %)"
    )
    print(f"     {'bar':<8}{'nodes':>8}{'genes':>9}{'genes named':>20}{'of shortfall':>14}{'nodes named':>20}")
    for bar in BARS:
        cleared = [
            lid
            for lid, key in assigned.items()
            if cells.get(key) and cells[key].agreement is not None and cells[key].agreement > bar
        ]
        genes = internal["checkable"][1] + sum(members[lid][0] for lid in cleared)
        gene_shift = f"{100 * own_genes / total_genes:.1f} → {100 * (own_genes + genes) / total_genes:.1f} %"
        node_shift = (
            f"{100 * len(levels) / len(members):.1f} → {100 * (len(levels) + len(cleared)) / len(members):.1f} %"
        )
        print(
            f"     {f'> {bar:.2f}':<8}{internal['checkable'][0] + len(cleared):>8,}{genes:>9,}"
            f"{gene_shift:>20}{100 * genes / shortfall:>13.1f}%{node_shift:>20}"
        )

    targets = [lid for lid in members if lid not in levels and neighbours.get(lid)]
    base = len(levels) / len(members)
    total_edges = sum(len(neighbours[lid]) for lid in targets)
    named_edges = sum(1 for lid in targets for _r, nb, _c in neighbours[lid] if nb in levels)
    zero = sum(1 for lid in targets if not any(nb in levels for _r, nb, _c in neighbours[lid]))
    print("\n  6. THE DARK SET — is the 5-neighbour cap the limit, or does the unknown neighbour itself?")
    print(
        f"     base rate {100 * base:.1f} % of nodes carry the axis, but only "
        f"{100 * named_edges / total_edges:.1f} % of unannotated nodes' neighbours do"
    )
    print(
        f"     {100 * zero / len(targets):.1f} % have ZERO annotated neighbours, against "
        f"{100 * (1 - base) ** 5:.1f} % if status were spread at random "
        f"→ {(zero / len(targets)) / (1 - base) ** 5:.2f}x concentrated"
    )
    return {
        "cells": cells,
        "members": members,
        "bands": bands,
        "levels": levels,
        "assigned": assigned,
        "internal": internal,
        "gene_calls": gene_calls,
    }


def prevalence(result: dict, bar: float = 0.90) -> None:
    """Section 7 — where the lift lands.

    ⚠ Internal inference counts ONLY where it is checkable, exactly as §5 does: a
    single-annotated-gene node has no measured confidence and must not inflate a band.
    """
    members, bands, levels = result["members"], result["bands"], result["levels"]
    gene_calls, cells = result["gene_calls"], result["cells"]
    print(f"\n  7. WHERE THE LIFT LANDS — by prevalence band, bar > {bar:.2f}")
    print(f"     {'band':<11}{'nodes':>8}{'named':>8}{'+new':>7}{'→ nodes':>9}{'genes':>10}{'named':>8}{'→ genes':>9}")
    per = defaultdict(lambda: [0, 0, 0, 0, 0, 0])
    for locus_id, (n_members, n_annotated) in members.items():
        row = per[bands[locus_id]]
        row[0] += 1
        row[4] += n_members
        row[5] += n_annotated
        if locus_id in levels:
            row[1] += 1
            if gene_calls.get(locus_id, 0) > 1:
                row[3] += max(0, n_members - n_annotated)
        elif locus_id in result["assigned"]:
            cell = cells.get(result["assigned"][locus_id])
            if cell and cell.agreement is not None and cell.agreement > bar:
                row[2] += 1
                row[3] += n_members
    for name in ("CORE", "SOFT_CORE", "SHELL", "RARE", "CLOUD"):
        if name not in per:
            continue
        nodes, named, new, gained, genes, before = per[name]
        print(
            f"     {name:<11}{nodes:>8,}{named:>8,}{new:>7,}"
            f"{100 * (named + new) / nodes:>8.1f}%{genes:>10,}"
            f"{100 * before / genes:>7.1f}%{100 * (before + gained) / genes:>8.1f}%"
        )


def pfam_crosstab(session: Session, catalogue: str) -> None:
    """Controlled vocabulary × Pfam — what is uncharacterised rather than merely un-COGged."""
    rows = session.execute(
        text("""
            select exists (select 1 from locus_annotation_entry a
                            where a.locus_id = l.locus_id
                              and a.annotation_kind in ('COG_ORTHOGROUP','EC_NUMBER',
                                                        'KEGG_ORTHOLOGY','GENE_ONTOLOGY_SLIM')),
                   coalesce(l.pfam_architecture_count, 0) > 0,
                   count(*), sum(l.member_gene_count)
              from locus l join pangenome p using (pangenome_id)
             where p.catalogue_key = :catalogue
             group by 1, 2 order by 1, 2
        """),
        {"catalogue": catalogue},
    ).all()
    print(f"\n  8. {catalogue} — what is actually dark? controlled vocabulary × Pfam")
    print(f"     {'controlled vocab':<18}{'Pfam':<7}{'nodes':>8}{'genes':>10}")
    for has_vocabulary, has_pfam, nodes, genes in rows:
        print(f"     {str(bool(has_vocabulary)):<18}{str(bool(has_pfam)):<7}{nodes:>8,}{genes:>10,}")


if __name__ == "__main__":
    url = os.environ.get("BACATLAS_DATABASE_URL")
    if not url:
        sys.exit("BACATLAS_DATABASE_URL is not set — see backend/README.md")
    print(f"min pairs before a cell reports a rate: {MIN_PAIRS}")
    with Session(create_engine(url)) as database:
        for catalogue_key in CATALOGUES:
            for annotation_kind in KINDS:
                prevalence(report(database, catalogue_key, annotation_kind))
            pfam_crosstab(database, catalogue_key)
