"""Every number in `nuna/docs/model_evaluation/COG_function_inference.md`, from the live database.

Run from `backend/` with `BACATLAS_DATABASE_URL` set and both `*-nuna4` catalogues loaded:

    python scripts/measure_cog_function_inference.py

⭐ **The calibration is COMPUTED, never hard-coded.** The agreement rate a panel quotes beside a
neighbour is the thing the reader judges the suggestion by, so it has to be re-derived from whatever
catalogue is loaded. A literal copied out of the document would keep reading correctly long after the
catalogue moved underneath it — the failure mode recorded in `PROJECT_STATE.md` §6, 2026-10-03
("a frozen oracle cannot detect a drift it shares").

⚠ **Three shapes in the data that a naive version gets wrong**, each measured before it was handled:
  * an EC term value is a comma-joined SET of codes (1,034 of 9,398), so `split_part(v,'.',4)` on
    `1.6.5.9,7.1.1.-` returns `9,7` — a code neither side holds;
  * 2,257 of 9,398 carry `-` placeholders, so a pair can be INCOMPARABLE at a level without
    disagreeing, and a per-level denominator makes the rate non-monotone as the level relaxes;
  * a node with exactly ONE annotated gene is unanimous by construction — the absence of evidence,
    not evidence — and is counted apart from the checkable ones throughout.
"""
from __future__ import annotations

import os
import random
import sys
from collections import Counter, defaultdict
from statistics import median

import psycopg

#: David, 2026-10-03: below 0.96 merged into one tier (0.95-0.96 and 0.90-0.95 were
#: indistinguishable), and nothing below 0.90 is called at all — "definitely too remote to call".
TIERS = [(">= 0.99", 0.99, 1.01), ("0.98-0.99", 0.98, 0.99), ("0.97-0.98", 0.97, 0.98),
         ("0.96-0.97", 0.96, 0.97), ("0.90-0.96", 0.90, 0.96)]
NOT_CALLED = "< 0.90 not called"
#: the level quoted per tier — a DECISION (David, 2026-10-03), not a measurement
QUOTED = {
    "EC_NUMBER": {">= 0.99": 4, "0.98-0.99": 3, "0.97-0.98": 1, "0.96-0.97": 1, "0.90-0.96": 1},
    "COG_ORTHOGROUP": {">= 0.99": 2, "0.98-0.99": 2, "0.97-0.98": 1, "0.96-0.97": 1, "0.90-0.96": 1},
}
LEVELS = {"EC_NUMBER": [4, 3, 2, 1], "COG_ORTHOGROUP": [2, 1]}
LEVEL_NAME = {
    "EC_NUMBER": {4: "L4 full code", 3: "L3 sub-subclass", 2: "L2 subclass", 1: "L1 class"},
    "COG_ORTHOGROUP": {2: "L2 orthogroup", 1: "L1 category"},
}
COUNT_COLUMN = {"EC_NUMBER": "ec_annotated_member_count",
                "COG_ORTHOGROUP": "cog_annotated_member_count"}
GENE_PREDICATE = {"EC_NUMBER": "f.ec_numbers is not null", "COG_ORTHOGROUP": "f.cog_id is not null"}
CATALOGUES = ("ecoli-nuna4", "kp-nuna4")
BARS = (0.95, 0.90, 0.80)
MIN_PAIRS = 30          # below this a cell reports its n and no rate
NULL_DRAWS = 200
SEED = 20261003


def dsn() -> str:
    """The psycopg DSN, from the SQLAlchemy URL the rest of the backend uses."""
    url = os.environ.get("BACATLAS_DATABASE_URL")
    if not url:
        sys.exit("BACATLAS_DATABASE_URL is not set — see backend/README.md")
    return url.replace("postgresql+psycopg://", "postgresql://", 1)


# ----------------------------------------------------------------- the two ladders
def ec_claims(raw: str) -> dict[int, frozenset[str]]:
    """Level -> the k-field prefixes this annotation RESOLVES to; empty means it does not state one."""
    out: dict[int, set[str]] = {k: set() for k in (1, 2, 3, 4)}
    for code in (c.strip() for c in raw.split(",")):
        fields = code.split(".")
        for k in (1, 2, 3, 4):
            if len(fields) >= k and all(f.isdigit() for f in fields[:k]):
                out[k].add(".".join(fields[:k]))
    return {k: frozenset(v) for k, v in out.items()}


def cog_claims(cog_id: str | None, categories: list[str] | None) -> dict[int, frozenset[str]]:
    """The two COG rungs — the orthogroup accession, then its functional-category letters.

    ⚠ L1 is `locus.modal_cog_categories`: the modal CONCATENATED category set over the node's
    annotated members, not the category of the quoted orthogroup, so agreement is set overlap.
    """
    return {2: frozenset([cog_id]) if cog_id else frozenset(),
            1: frozenset(categories or ())}


def wilson(hits: int, n: int) -> tuple[float, float]:
    """A 95 % Wilson interval — it stays inside [0, 1] at the 99 %+ rates the top tiers reach."""
    if n == 0:
        return (float("nan"), float("nan"))
    p, d = hits / n, 1 + 1.96**2 / n
    centre = (p + 1.96**2 / (2 * n)) / d
    half = 1.96 * ((p * (1 - p) / n + 1.96**2 / (4 * n * n)) ** 0.5) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def tier_of(cosine: float) -> str:
    """Which similarity tier a neighbour's cosine falls in, or `NOT_CALLED`."""
    for name, lo, hi in TIERS:
        if lo <= cosine < hi:
            return name
    return NOT_CALLED


def load(conn, catalogue: str, kind: str):
    """Per node: member and annotated gene counts, prevalence band, top-rank call, and neighbours."""
    with conn.cursor() as cur:
        cur.execute(
            f"""
            select l.locus_id, l.member_gene_count, l.{COUNT_COLUMN[kind]}, l.prevalence_band,
                   a.term_value, l.modal_cog_categories
              from locus l join pangenome p using (pangenome_id)
              left join locus_annotation_entry a
                     on a.locus_id = l.locus_id and a.annotation_kind = %s
                    and a.rank_within_locus = 0
             where p.catalogue_key = %s
            """,
            (kind, catalogue),
        )
        members, bands, claims = {}, {}, {}
        for locus_id, n_members, n_annotated, prevalence, term, cats in cur.fetchall():
            members[locus_id] = (n_members, n_annotated or 0)
            bands[locus_id] = prevalence
            found = ec_claims(term) if (kind == "EC_NUMBER" and term) else (
                cog_claims(term, cats) if kind != "EC_NUMBER" else None)
            if found and any(found.values()):
                claims[locus_id] = found

        cur.execute(
            """
            select n.locus_id, n.rank, n.neighbour_locus_id, n.cross_similarity
              from locus_nearest_locus n
              join locus l on l.locus_id = n.locus_id
              join pangenome p on p.pangenome_id = l.pangenome_id
             where p.catalogue_key = %s and n.representation = 'ESM'
               and n.cross_similarity is not null
             order by n.locus_id, n.rank
            """,
            (catalogue,),
        )
        neighbours = defaultdict(list)
        for locus_id, rank, nb, cosine in cur.fetchall():
            neighbours[locus_id].append((rank, nb, cosine))

        #: how many of each node's genes carry a call — the checkable / unverifiable split
        cur.execute(
            f"""
            select m.locus_id, count(*)
              from gene_locus_membership m
              join locus l on l.locus_id = m.locus_id
              join pangenome p on p.pangenome_id = l.pangenome_id
              join gene_functional_annotation f
                   on f.genome_id = m.genome_id and f.flat_index = m.flat_index
             where p.catalogue_key = %s and {GENE_PREDICATE[kind]}
             group by 1
            """,
            (catalogue,),
        )
        gene_calls = dict(cur.fetchall())
    return members, bands, claims, neighbours, gene_calls


def gene_level_claims(conn, catalogue: str, kind: str):
    """Per node, its annotated genes' values grouped with a count — for the within-node check."""
    with conn.cursor() as cur:
        cur.execute(
            f"""
            select m.locus_id, f.cog_id, f.cog_categories, f.ec_numbers, count(*)
              from gene_locus_membership m
              join locus l on l.locus_id = m.locus_id
              join pangenome p on p.pangenome_id = l.pangenome_id
              join gene_functional_annotation f
                   on f.genome_id = m.genome_id and f.flat_index = m.flat_index
             where p.catalogue_key = %s and {GENE_PREDICATE[kind]}
             group by 1, 2, 3, 4
            """,
            (catalogue,),
        )
        groups = defaultdict(list)
        for locus_id, cog_id, _cats, ecs, n in cur.fetchall():
            if kind == "EC_NUMBER":
                groups[locus_id].append((frozenset(ecs or ()), n))
            else:
                groups[locus_id].append((frozenset([cog_id]) if cog_id else frozenset(), n))
    return groups


def calibrate(claims, neighbours, kind):
    """Agreement and the matched-composition null per (tier, level), pooled over ALL ranks.

    ⭐ All ranks, not rank 1: 65 % of assignments are supplied by ranks 2-5, that is the population
    a panel serves, and the rank effect at fixed cosine was <= 0.3 pp in the >= 0.99 tiers and
    flipped SIGN between species in the weak ones — no consistent effect, so splitting would fit noise.
    """
    levels = LEVELS[kind]
    rng = random.Random(SEED)
    pools = {k: [lid for lid, c in claims.items() if c[k]] for k in levels}
    observed = defaultdict(lambda: [0, 0])
    chance = defaultdict(lambda: [0, 0])
    for locus_id, mine in claims.items():
        for _rank, nb, cosine in neighbours.get(locus_id, []):
            if nb not in claims:
                continue
            tier = tier_of(cosine)
            if tier == NOT_CALLED:
                continue
            theirs = claims[nb]
            for k in levels:
                if not (mine[k] and theirs[k]):
                    continue
                observed[(tier, k)][0] += 1
                observed[(tier, k)][1] += bool(mine[k] & theirs[k])
                for _ in range(NULL_DRAWS):
                    partner = claims[pools[k][rng.randrange(len(pools[k]))]]
                    if partner[k]:
                        chance[(tier, k)][0] += 1
                        chance[(tier, k)][1] += bool(mine[k] & partner[k])
    return observed, chance


def assign(members, claims, neighbours, gene_calls, kind):
    """Internal inference over nodes that HAVE a call; the neighbour walk over those that do not.

    The two populations are disjoint and together they exhaust the unannotated genes.
    """
    levels = LEVELS[kind]
    internal = {"checkable": [0, 0], "one gene": [0, 0]}
    transfer: dict = defaultdict(lambda: [0, 0])
    supplier_rank: Counter[int] = Counter()
    sizes = defaultdict(list)
    uncalled = [0, 0]
    unreachable = [0, 0]
    assigned_nodes: dict[int, tuple[str, int]] = {}

    for locus_id, (n_members, n_annotated) in members.items():
        if locus_id in claims:
            bucket = "one gene" if gene_calls.get(locus_id, 0) <= 1 else "checkable"
            internal[bucket][0] += 1
            internal[bucket][1] += max(0, n_members - n_annotated)
            sizes["has its own call"].append(n_members)
            continue
        hit = next(((r, nb, c) for r, nb, c in neighbours.get(locus_id, []) if nb in claims), None)
        if hit is None:
            unreachable[0] += 1
            unreachable[1] += n_members
            sizes["no annotated neighbour in 5"].append(n_members)
            continue
        rank, nb, cosine = hit
        tier = tier_of(cosine)
        permitted = QUOTED[kind].get(tier)
        quotable = (next((k for k in levels if k <= permitted and claims[nb][k]), None)
                    if permitted else None)
        if quotable is None:
            uncalled[0] += 1
            uncalled[1] += n_members
            sizes["neighbour too remote"].append(n_members)
            continue
        transfer[(tier, quotable)][0] += 1
        transfer[(tier, quotable)][1] += n_members
        supplier_rank[rank] += 1
        assigned_nodes[locus_id] = (tier, quotable)
        sizes["assigned from a neighbour"].append(n_members)
    return internal, transfer, supplier_rank, sizes, uncalled, unreachable, assigned_nodes


# ----------------------------------------------------------------- report
def report(conn, catalogue: str, kind: str) -> dict:
    """Print all eight sections for one catalogue and one axis; return what §7 needs."""
    members, bands, claims, neighbours, gene_calls = load(conn, catalogue, kind)
    observed, chance = calibrate(claims, neighbours, kind)
    internal, transfer, supplier_rank, sizes, uncalled, unreachable, assigned = assign(
        members, claims, neighbours, gene_calls, kind)
    levels = LEVELS[kind]

    print(f"\n{'=' * 100}\n{catalogue} · {kind}")

    # --- 1. within-node consistency
    groups = gene_level_claims(conn, catalogue, kind)
    print("\n  1. WITHIN-NODE CONSISTENCY — do a node's annotated genes agree?")
    print(f"     {'level':<17}{'nodes':>8}{'1 gene':>9}{'checkable':>11}{'unanimous':>11}{'support':>9}")
    for k in levels:
        single = multi = unanimous = 0
        support = []
        for g in groups.values():
            folded = [(frozenset(v for v in (prefix(c, k) if kind == "EC_NUMBER" else c
                                             for c in vals) if v), n) for vals, n in g] \
                if kind == "EC_NUMBER" else [(vals, n) for vals, n in g]
            folded = [(v, n) for v, n in folded if v]
            genes = sum(n for _, n in folded)
            if genes == 0:
                continue
            tally: Counter[str] = Counter()
            for values, n in folded:
                for v in values:
                    tally[v] += n
            if genes == 1:
                single += 1
                continue
            multi += 1
            top = tally.most_common(1)[0][1]
            support.append(top / genes)
            unanimous += top == genes
        print(f"     {LEVEL_NAME[kind][k]:<17}{single + multi:>8,}{single:>9,}{multi:>11,}"
              f"{(f'{100 * unanimous / multi:.1f}%' if multi else '—'):>11}"
              f"{(f'{sum(support) / len(support):.3f}' if support else '—'):>9}")

    # --- 2. the ladder
    print("\n  2. THE LADDER — nodes with no call of their own, labelled from a neighbour")
    print(f"     {'tier':<12}{'quoted':<20}{'nodes':>7}{'genes':>9}{'agree':>8}{'95% CI':>10}"
          f"{'n':>7}{'chance':>8}{'lift':>7}")
    for name, _, _ in TIERS:
        for k in levels:
            nodes, genes = transfer[(name, k)]
            if not nodes:
                continue
            n, hits = observed[(name, k)]
            fallback = " ←fallback" if k != QUOTED[kind][name] else ""
            if n < MIN_PAIRS:
                print(f"     {name:<12}{LEVEL_NAME[kind][k] + fallback:<20}{nodes:>7}{genes:>9,}"
                      f"{'n too small':>25} {n:>6}")
                continue
            lo, hi = wilson(hits, n)
            cn, ch = chance[(name, k)]
            q = ch / cn if cn else float("nan")
            print(f"     {name:<12}{LEVEL_NAME[kind][k] + fallback:<20}{nodes:>7}{genes:>9,}"
                  f"{100 * hits / n:7.1f}%{f'{100 * lo:.0f}-{100 * hi:.0f}':>10}{n:>7}"
                  f"{100 * q:7.1f}%{(hits / n) / q if q else float('inf'):6.1f}x")
    print(f"     {NOT_CALLED:<32}{uncalled[0]:>7}{uncalled[1]:>9,}")
    print(f"     {'no annotated neighbour in 5':<32}{unreachable[0]:>7}{unreachable[1]:>9,}")
    print("     supplying rank: " + " · ".join(f"{r}: {supplier_rank[r]:,}"
                                               for r in sorted(supplier_rank)))

    # --- 3. internal inference, and the node sizes that explain the gene totals
    print("\n  3. INTERNAL NODE INFERENCE — the node's modal call applied to its unannotated genes")
    print(f"     {'>= 2 annotated genes (checkable)':<34}{internal['checkable'][0]:>7} nodes"
          f"{internal['checkable'][1]:>9,} genes")
    print(f"     {'exactly 1 annotated gene (⛔ no check)':<34}{internal['one gene'][0]:>7} nodes"
          f"{internal['one gene'][1]:>9,} genes")
    print("\n  4. HOW BIG ARE THE NODES? — why a gene-weighted percentage barely moves")
    print(f"     {'group':<30}{'nodes':>8}{'genes':>10}{'mean':>8}{'median':>8}")
    for name, values in sorted(sizes.items(), key=lambda kv: -sum(kv[1])):
        print(f"     {name:<30}{len(values):>8,}{sum(values):>10,}"
              f"{sum(values) / len(values):>8.1f}{median(values):>8.0f}")

    # --- 5. the three denominators, at each bar
    total_genes = sum(n for n, _ in members.values())
    own_genes = sum(a for _, a in members.values())
    shortfall = total_genes - own_genes
    print(f"\n  5. COVERAGE — {own_genes:,} of {total_genes:,} genes already named "
          f"({100 * own_genes / total_genes:.1f} %); {len(claims):,} of {len(members):,} nodes "
          f"({100 * len(claims) / len(members):.1f} %)")
    print(f"     {'bar':<8}{'nodes':>8}{'genes':>9}{'genes named':>20}{'of shortfall':>14}"
          f"{'nodes named':>20}")
    base_nodes = 100 * len(claims) / len(members)
    base_genes = 100 * own_genes / total_genes
    for bar in BARS:
        from_neighbours = [lid for lid, (tier, k) in assigned.items()
                           if observed[(tier, k)][0] >= MIN_PAIRS
                           and observed[(tier, k)][1] / observed[(tier, k)][0] > bar]
        nodes = internal["checkable"][0] + len(from_neighbours)
        genes = internal["checkable"][1] + sum(members[lid][0] for lid in from_neighbours)
        after_genes = 100 * (own_genes + genes) / total_genes
        after_nodes = 100 * (len(claims) + len(from_neighbours)) / len(members)
        print(f"     {f'> {bar:.2f}':<8}{nodes:>8,}{genes:>9,}"
              f"{f'{base_genes:.1f} → {after_genes:.1f} %':>20}{100 * genes / shortfall:>13.1f}%"
              f"{f'{base_nodes:.1f} → {after_nodes:.1f} %':>20}")

    # --- 6. the dark set
    targets = [lid for lid in members if lid not in claims and neighbours.get(lid)]
    base = len(claims) / len(members)
    nb_total = sum(len(neighbours[lid]) for lid in targets)
    nb_named = sum(1 for lid in targets for _r, nb, _c in neighbours[lid] if nb in claims)
    zero = sum(1 for lid in targets if not any(nb in claims for _r, nb, _c in neighbours[lid]))
    print("\n  6. THE DARK SET — is the 5-neighbour cap the limit, or does the unknown neighbour itself?")
    print(f"     base rate {100 * base:.1f} % of nodes carry the axis, but only "
          f"{100 * nb_named / nb_total:.1f} % of unannotated nodes' neighbours do")
    print(f"     {100 * zero / len(targets):.1f} % have ZERO annotated neighbours, against "
          f"{100 * (1 - base) ** 5:.1f} % if status were spread at random "
          f"→ {(zero / len(targets)) / (1 - base) ** 5:.2f}x concentrated")
    return {"observed": observed, "members": members, "bands": bands, "claims": claims,
            "assigned": assigned, "internal": internal, "gene_calls": gene_calls}


def prefix(code: str, level: int) -> str | None:
    """An EC code's first `level` fields, or None where it does not resolve that deep."""
    fields = code.split(".")
    return ".".join(fields[:level]) if len(fields) >= level and all(
        f.isdigit() for f in fields[:level]) else None


def prevalence(result: dict, kind: str, bar: float = 0.90) -> None:
    """Where the lift lands, by prevalence band.

    ⚠ Internal inference counts ONLY where it is checkable, exactly as §5 does — a
    single-annotated-gene node has no measured confidence and must not inflate a band.
    """
    members, bands, claims = result["members"], result["bands"], result["claims"]
    gene_calls = result["gene_calls"]
    print(f"\n  7. WHERE THE LIFT LANDS — by prevalence band, bar > {bar:.2f}")
    print(f"     {'band':<11}{'nodes':>8}{'named':>8}{'+new':>7}{'→ nodes':>9}"
          f"{'genes':>10}{'named':>8}{'→ genes':>9}")
    per = defaultdict(lambda: [0, 0, 0, 0, 0, 0])
    for locus_id, (n_members, n_annotated) in members.items():
        row = per[bands[locus_id]]
        row[0] += 1
        row[4] += n_members
        row[5] += n_annotated
        if locus_id in claims:
            row[1] += 1
            if gene_calls.get(locus_id, 0) > 1:
                row[3] += max(0, n_members - n_annotated)
        elif locus_id in result["assigned"]:
            tier, k = result["assigned"][locus_id]
            n, hits = result["observed"][(tier, k)]
            if n >= MIN_PAIRS and hits / n > bar:
                row[2] += 1
                row[3] += n_members
    for name in ("CORE", "SOFT_CORE", "SHELL", "RARE", "CLOUD"):
        if name not in per:
            continue
        nodes, named, new, gained, genes, before = per[name]
        print(f"     {name:<11}{nodes:>8,}{named:>8,}{new:>7,}{100 * (named + new) / nodes:>8.1f}%"
              f"{genes:>10,}{100 * before / genes:>7.1f}%{100 * (before + gained) / genes:>8.1f}%")


def pfam_crosstab(conn, catalogue: str) -> None:
    """Controlled vocabulary × Pfam — what is actually uncharacterised rather than merely un-COGged."""
    with conn.cursor() as cur:
        cur.execute(
            """
            select exists (select 1 from locus_annotation_entry a
                            where a.locus_id = l.locus_id
                              and a.annotation_kind in ('COG_ORTHOGROUP','EC_NUMBER',
                                                        'KEGG_ORTHOLOGY','GENE_ONTOLOGY_SLIM')),
                   coalesce(l.pfam_architecture_count, 0) > 0,
                   count(*), sum(l.member_gene_count)
              from locus l join pangenome p using (pangenome_id)
             where p.catalogue_key = %s
             group by 1, 2 order by 1, 2
            """,
            (catalogue,),
        )
        print(f"\n  8. {catalogue} — what is actually dark? controlled vocabulary × Pfam")
        print(f"     {'controlled vocab':<18}{'Pfam':<7}{'nodes':>8}{'genes':>10}")
        for cv, pfam, nodes, genes in cur.fetchall():
            print(f"     {str(bool(cv)):<18}{str(bool(pfam)):<7}{nodes:>8,}{genes:>10,}")


if __name__ == "__main__":
    with psycopg.connect(dsn()) as connection:
        for catalogue_key in CATALOGUES:
            for annotation_kind in ("COG_ORTHOGROUP", "EC_NUMBER"):
                outcome = report(connection, catalogue_key, annotation_kind)
                prevalence(outcome, annotation_kind)
            pfam_crosstab(connection, catalogue_key)
