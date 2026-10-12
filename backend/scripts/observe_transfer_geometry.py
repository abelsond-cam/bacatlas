r"""S1 — look at the geometry before modelling it. ⛔ No model is fitted here, on purpose.

David, 2026-10-12: *"This is hard to do in one step. It requires observation of the values and checking
the model at each step and carefully extending it."* So this script only measures and prints. Three
questions, all on artifacts already on disk, none needing cluster time:

1. ⭐ **Max or median?** Paired, on the recipients where the two linkage rules disagree about which
   annotated donor is nearest — McNemar, plus what kind of donor each rule reaches (size, annotation
   depth). ⛔ No size-correction model: David, *"I'm really reticent to get into games of modeling max
   with size … an overcomplex model is a bad model."*
2. ⭐ **What does the candidate ladder look like?** How often do the two representations name the same
   nearest annotated node; how often do the top-N agree on a label; what the similarity drop-off looks
   like when they agree against when they do not.
3. ⭐ **Which representation finds the better donor, per claim?** Paired on the same recipients.
   ⛔ Entered with no expectation: one earlier session concluded Bacformer wins for biological process,
   a later check said ESM wins everywhere, both were measured on data this work replaces, and they
   contradict each other.

⚠ **No held-out mask is needed for any of this**, and that is worth stating because the rest of the
programme needs one. Here the thing being predicted is the recipient's own claim and the predictor is
a *different* node's claim, so nothing the search sees is derived from the answer. The mask exists to
stop a propagated label being scored against itself; nothing is propagated here.

Run from `backend/`::

    python scripts/observe_transfer_geometry.py \
        --cross-esm ../../nuna/data/proc/analysis/cross_species/ecoli_to_kp_..._cross_species_neighbours.tsv \
        --cross-bacformer ../../nuna/data/proc/analysis/cross_species/..._BACFORMER_...tsv
"""

from __future__ import annotations

import argparse
import os
from collections import Counter

from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from bacatlas_backend.transfer import claims as claim_module
from bacatlas_backend.transfer.candidates import read_cross_map, within_species

CATALOGUES = ("ecoli-nuna4", "kp-nuna4")


def mcnemar(both_right: int, max_only: int, median_only: int, both_wrong: int) -> tuple[float, str]:
    """Exact-ish McNemar on the DISCORDANT pairs only — the whole point of pairing.

    ⛔ The concordant cells carry no information about which rule is better; a test that used them
    would be answering *"are these donors any good"*, not *"is one rule better than the other"*.
    Returns the normal-approximation statistic and a plain-English verdict.
    """
    discordant = max_only + median_only
    if discordant < 10:
        return float("nan"), f"too few discordant pairs ({discordant}) to say"
    statistic = (abs(max_only - median_only) - 1) ** 2 / discordant
    winner = "max" if max_only > median_only else "median"
    # chi-square, 1 df: 3.84 at p=0.05, 6.63 at p=0.01
    if statistic < 3.84:
        return statistic, f"no measurable difference (chi2={statistic:.2f}, needs 3.84)"
    return statistic, f"⭐ {winner} is better (chi2={statistic:.2f}, p<{'0.01' if statistic > 6.63 else '0.05'})"


def first_annotated(ordered: list[dict], label_to_id: dict[str, int], labelled: dict[int, frozenset[str]]):
    """The first donor in this order that actually carries the claim — or `None`."""
    for row in ordered:
        donor = label_to_id.get(row["donor_label"])
        if donor is not None and donor in labelled:
            return donor, row
    return None


def q1_max_vs_median(recipient, donor_cat, cross_map, provenance) -> None:
    """⭐ The paired experiment: where the two rules disagree, which one is right more often?"""
    print(f"\n{'=' * 112}\nQ1 — MAX OR MEDIAN? paired, on the recipients where they disagree\n{'=' * 112}")
    print(f"  map: {provenance['donor']} → {provenance['recipient']}, cap {provenance['cap']}, "
          f"{provenance['n_pairs_written']:,} pairs")
    print("  ⛔ 'disagree' means the first ANNOTATED donor differs, which is the choice that matters.\n")
    label_to_id = donor_cat.id_of_label()
    print(f"  {'claim':<22}{'scoreable':>10}{'disagree':>10}{'max right':>11}{'med right':>11}"
          f"{'max only':>10}{'med only':>10}   verdict")
    detail: dict[str, dict] = {}
    for kind in claim_module.CLAIM_KINDS:
        truth = recipient.claims[kind.key]
        labelled = donor_cat.claims[kind.key]
        cells = Counter()
        sizes = {"max": [], "median": []}
        depth = {"max": [], "median": []}
        scoreable = 0
        for node, true_claim in truth.items():
            ordered = cross_map.get(recipient.node_label[node])
            if not ordered:
                continue
            by_max = first_annotated(sorted(ordered, key=lambda r: r["rank_cross_max"]), label_to_id, labelled)
            by_med = first_annotated(sorted(ordered, key=lambda r: r["rank_cross"]), label_to_id, labelled)
            if by_max is None or by_med is None:
                continue
            scoreable += 1
            if by_max[0] == by_med[0]:
                continue
            right_max = claim_module.agrees(labelled[by_max[0]], true_claim)
            right_med = claim_module.agrees(labelled[by_med[0]], true_claim)
            cells[(right_max, right_med)] += 1
            for name, (donor, row) in (("max", by_max), ("median", by_med)):
                sizes[name].append(row["n_donor_genes"])
                depth[name].append(donor_cat.genes_carrying[kind.key].get(donor, 0))
        disagree = sum(cells.values())
        if not disagree:
            print(f"  {kind.label:<22}{scoreable:>10,}{0:>10}   —")
            continue
        max_only, med_only = cells[(True, False)], cells[(False, True)]
        right_max = cells[(True, True)] + max_only
        right_med = cells[(True, True)] + med_only
        _stat, verdict = mcnemar(cells[(True, True)], max_only, med_only, cells[(False, False)])
        print(f"  {kind.label:<22}{scoreable:>10,}{disagree:>10,}{right_max / disagree:>11.1%}"
              f"{right_med / disagree:>11.1%}{max_only:>10,}{med_only:>10,}   {verdict}")
        detail[kind.key] = {"sizes": sizes, "depth": depth}

    print("\n  H1 — what kind of donor does each rule REACH? (mean over the disagreements)")
    print(f"  {'claim':<22}{'max: donor genes':>18}{'median: donor genes':>21}"
          f"{'max: annotated':>16}{'median: annotated':>19}")
    for kind in claim_module.CLAIM_KINDS:
        if kind.key not in detail:
            continue
        d = detail[kind.key]
        mean = lambda xs: sum(xs) / len(xs) if xs else float("nan")  # noqa: E731
        print(f"  {kind.label:<22}{mean(d['sizes']['max']):>18.1f}{mean(d['sizes']['median']):>21.1f}"
              f"{mean(d['depth']['max']):>16.1f}{mean(d['depth']['median']):>19.1f}")
    print("  ⚠ 'annotated' = member genes of the donor that state the claim themselves — H2.")


def _first_annotated_walk(walk, representation, labelled):
    """The first neighbour in this representation's stored walk that carries the claim."""
    for name, _rank, neighbour, similarity in sorted(walk, key=lambda step: step[1]):
        if name == representation and neighbour in labelled:
            return neighbour, similarity
    return None


def q2_candidate_ladder(catalogue, walks) -> None:
    """⭐ Look at the values before modelling them — David's explicit instruction."""
    print(f"\n{'=' * 112}\nQ2 — THE CANDIDATE LADDER, {catalogue.catalogue} (within species)\n{'=' * 112}")
    print("  ⛔ Observation only. How much do the two representations actually disagree, and does")
    print("     disagreement look different from agreement?\n")
    print(f"  {'claim':<22}{'both reach':>11}{'same node':>11}{'rate':>8}{'agree on label':>15}{'rate':>8}"
          f"{'ESM sim: agree':>15}{'differ':>9}")
    for kind in claim_module.CLAIM_KINDS:
        labelled = catalogue.claims[kind.key]
        both = same = agree = 0
        similarity = {"agree": [], "differ": []}
        for node, walk in walks.items():
            if node in labelled:
                continue  # already annotated; nothing to transfer
            esm = _first_annotated_walk(walk, "ESM", labelled)
            bacformer = _first_annotated_walk(walk, "BACFORMER", labelled)
            if esm is None or bacformer is None:
                continue
            both += 1
            if esm[0] == bacformer[0]:
                same += 1
            concurring = claim_module.agrees(labelled[esm[0]], labelled[bacformer[0]])
            agree += concurring
            similarity["agree" if concurring else "differ"].append(esm[1])
        if not both:
            print(f"  {kind.label:<22}{0:>11}")
            continue
        mean = lambda xs: sum(xs) / len(xs) if xs else float("nan")  # noqa: E731
        print(f"  {kind.label:<22}{both:>11,}{same:>11,}{same / both:>8.1%}{agree:>15,}{agree / both:>8.1%}"
              f"{mean(similarity['agree']):>15.4f}{mean(similarity['differ']):>9.4f}")
    print("  ⚠ 'same node' = the two representations name the SAME first annotated donor.")
    print("    'agree on label' = their two donors carry intersecting claims, whether or not the same node.")


def q3_which_representation(catalogue, walks) -> None:
    """⭐ Paired: on the same recipients, whose first annotated donor is right more often?"""
    print(f"\n{'=' * 112}\nQ3 — WHICH REPRESENTATION FINDS THE BETTER DONOR, {catalogue.catalogue}\n{'=' * 112}")
    print("  ⛔ Paired on the recipients where the two name DIFFERENT donors — the concordant ones")
    print("     carry no information about which representation is better.\n")
    print(f"  {'claim':<22}{'differ':>9}{'ESM right':>11}{'BAC right':>11}{'ESM only':>10}{'BAC only':>10}   verdict")
    for kind in claim_module.CLAIM_KINDS:
        labelled = catalogue.claims[kind.key]
        truth = catalogue.claims[kind.key]
        cells = Counter()
        for node, true_claim in truth.items():
            walk = walks.get(node)
            if not walk:
                continue
            # ⛔ the node's own label must not be its own donor
            visible = {other: value for other, value in labelled.items() if other != node}
            esm = _first_annotated_walk(walk, "ESM", visible)
            bacformer = _first_annotated_walk(walk, "BACFORMER", visible)
            if esm is None or bacformer is None or esm[0] == bacformer[0]:
                continue
            cells[(claim_module.agrees(visible[esm[0]], true_claim),
                   claim_module.agrees(visible[bacformer[0]], true_claim))] += 1
        differ = sum(cells.values())
        if not differ:
            print(f"  {kind.label:<22}{0:>9}   —")
            continue
        esm_only, bac_only = cells[(True, False)], cells[(False, True)]
        both_right = cells[(True, True)]
        _stat, verdict = mcnemar(both_right, esm_only, bac_only, cells[(False, False)])
        verdict = verdict.replace("max is better", "ESM is better").replace("median is better", "Bacformer is better")
        print(f"  {kind.label:<22}{differ:>9,}{(both_right + esm_only) / differ:>11.1%}"
              f"{(both_right + bac_only) / differ:>11.1%}{esm_only:>10,}{bac_only:>10,}   {verdict}")


def main() -> None:
    """Measure and print. ⛔ Nothing is fitted, nothing is recorded."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cross-esm", required=True)
    parser.add_argument("--cross-bacformer", required=False)
    arguments = parser.parse_args()

    url = os.environ.get("BACATLAS_DATABASE_URL")
    if not url:
        raise SystemExit("FATAL: set BACATLAS_DATABASE_URL")
    engine = create_engine(url, future=True)

    with Session(engine) as session:
        loaded = {name: claim_module.load(session, name) for name in CATALOGUES}
        walks = {name: within_species(session, loaded[name]) for name in CATALOGUES}

    print(f"\n{'=' * 112}")
    print("S1 — OBSERVING THE TRANSFER GEOMETRY, BEFORE ANY MODEL")
    print(f"{'=' * 112}")
    for name, catalogue in loaded.items():
        counts = " · ".join(
            f"{kind.key} {len(catalogue.claims[kind.key]):,}" for kind in claim_module.CLAIM_KINDS
        )
        print(f"  {name:<12}{catalogue.node_total:,} nodes · {catalogue.gene_total:,} genes · annotated: {counts}")

    cross_map, provenance = read_cross_map(arguments.cross_esm, donor="ecoli", recipient="kp")
    q1_max_vs_median(loaded["kp-nuna4"], loaded["ecoli-nuna4"], cross_map, provenance)
    for name in CATALOGUES:
        q2_candidate_ladder(loaded[name], walks[name])
    for name in CATALOGUES:
        q3_which_representation(loaded[name], walks[name])

    print(f"\n{'=' * 112}")
    print("⛔ NOTHING HERE IS RECORDED, and no model has been fitted. These numbers go to David first.")
    print(f"{'=' * 112}\n")


if __name__ == "__main__":
    main()
