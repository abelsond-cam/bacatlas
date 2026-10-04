r"""ESM against Bacformer on the SAME cross-species node pairs — which ranks function better.

⛔⛔ **WHY THIS SCRIPT EXISTS RATHER THAN A COLUMN IN THE TIER TABLE.** The decided ladder is an **ESM**
calibration: 96.74 % of stored within-kp ESM neighbour edges sit at or above 0.90 against **1.17 %** of
Bacformer's (877 of 75,263), and 95.5 % of kp loci have their best Bacformer donor below the 0.90 floor
with **zero** at ≥0.99. Running a Bacformer map through those cuts empties every cell — not because
Bacformer carries no signal but because the thresholds belong to a different geometry. Re-cutting the
tiers per representation is a statistical decision (David's, 2026-10-03: *AUC first, re-calibrate only if
it earns it*), so this asks the prior question instead: **does the representation rank agreeing donors
above disagreeing ones at all?** A representation that cannot rank cannot be rescued by better cuts.

⛔ **AND WHY THE PAIRS MUST BE SHARED.** Each map's AUC over its *own* shortlist is not comparable: the
shortlists differ (ESM 365,105 pairs, Bacformer 557,834) and Bacformer's carries **50,906** disagreeing
COG pairs against ESM's **23,457**, so a shortlist holding more easy negatives earns a higher AUC for
free. `nuna/CLAUDE.md` §"Five things a harness cannot catch" (2): *a comparison must assert its own
coverage BEFORE it reports a difference.* Here both representations score the **intersection** of the two
shortlists, with the same labels, so the only thing that differs is the similarity.

⚠ **The intersection is a biased subset and the output says so**: it is the pairs *both* methods
independently surfaced, so it over-represents pairs that are close in both geometries. It is the only
common ground two different shortlists have, and the coverage is reported beside every number.

Run from `backend/` with `BACATLAS_DATABASE_URL` set and both `*-nuna4` catalogues loaded:

    python scripts/compare_representations_cross_species.py ESM=<esm map>.tsv BACFORMER=<bacformer map>.tsv
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from pathlib import Path

from sqlalchemy import create_engine, text
from sqlalchemy.orm import Session


def _transfer_module():
    """The single implementation of the fold, the ladder and the AUC — imported, never re-written."""
    path = Path(__file__).resolve().parent / "measure_cross_species_transfer.py"
    spec = importlib.util.spec_from_file_location("measure_cross_species_transfer", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


xs = _transfer_module()


def read_scores(path: str, *, donor: str, recipient: str) -> tuple[dict[tuple[str, str], float], dict]:
    """``{(recipient_node, donor_node): cross}`` plus the map's provenance, through the shared reader."""
    rows, provenance = xs.read_map(path, donor=donor, recipient=recipient)
    return {(row["recipient_node"], row["donor_node"]): row["cross"] for row in rows}, provenance


def main() -> None:
    """One AUC per (vocabulary, representation) over the shared pairs, with the coverage beside it."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "maps", nargs="+", metavar="NAME=PATH", help="two or more labelled maps, e.g. ESM=a.tsv BACFORMER=b.tsv"
    )
    parser.add_argument("--donor", default="ecoli")
    parser.add_argument("--recipient", default="kp")
    arguments = parser.parse_args()
    if len(arguments.maps) < 2:
        sys.exit("at least two maps are needed — this script exists to compare them")

    url = os.environ.get("BACATLAS_DATABASE_URL")
    if not url:
        sys.exit("BACATLAS_DATABASE_URL is not set — see backend/README.md")

    scores: dict[str, dict] = {}
    for item in arguments.maps:
        if "=" not in item:
            sys.exit(f"expected NAME=PATH, got {item!r}")
        name, path = item.split("=", 1)
        scores[name], provenance = read_scores(path, donor=arguments.donor, recipient=arguments.recipient)
        stated = (provenance.get("rep") or "esm").upper()
        # ⛔ the label must match what the map says it is, or the two columns are mislabelled and the
        # whole comparison reads backwards while looking entirely reasonable
        if stated != name.upper():
            sys.exit(
                f"FATAL: {Path(path).name} was built from {stated} but is labelled {name!r}. "
                "A mislabelled column inverts the comparison and nothing in the output would say so."
            )
        print(f"{name}: {len(scores[name]):,} shortlisted pairs  ({Path(path).name})")

    shared = set.intersection(*(set(table) for table in scores.values()))
    if not shared:
        sys.exit("FATAL: the maps share no node pair at all — are they the same direction and model?")
    print(
        f"\nshared by every shortlist: {len(shared):,} pairs — "
        + ", ".join(f"{len(shared) / len(table):.1%} of {name}'s" for name, table in scores.items())
    )
    print(
        "⚠ the intersection over-represents pairs close in EVERY geometry; it is the only common"
        " ground two different shortlists have."
    )

    catalogues = {"donor": f"{arguments.donor}-nuna4", "recipient": f"{arguments.recipient}-nuna4"}
    with Session(create_engine(url)) as database:
        ids = {
            side: {
                str(label): locus_id
                for label, locus_id in database.execute(
                    text("""
                        select l.node_label, l.locus_id
                          from locus l join pangenome p using (pangenome_id)
                         where p.catalogue_key = :catalogue
                    """),
                    {"catalogue": catalogue},
                ).all()
            }
            for side, catalogue in catalogues.items()
        }
        print(f"\n{'axis':<22}{'pairs':>9}{'agreeing':>10}" + "".join(f"{name + ' AUC':>16}" for name in scores))
        for kind in xs.KINDS:
            donor_levels, _ = xs.claims(database, catalogues["donor"], kind)
            recipient_levels, _ = xs.claims(database, catalogues["recipient"], kind)
            both = {**donor_levels, **recipient_levels}
            ladder = sorted(xs.LEVEL_LABEL[kind], reverse=True)
            scored: dict[str, list] = {name: [] for name in scores}
            agreeing = 0
            for recipient_label, donor_label in shared:
                recipient_id = ids["recipient"][recipient_label]
                donor_id = ids["donor"][donor_label]
                mine, theirs = both.get(recipient_id), both.get(donor_id)
                if mine is None or theirs is None:
                    continue
                level = next((lv for lv in ladder if mine[lv] and theirs[lv]), None)
                if level is None:
                    continue
                agree = bool(mine[level] & theirs[level])
                agreeing += agree
                for name in scores:
                    scored[name].append((scores[name][(recipient_label, donor_label)], agree))
            total = len(next(iter(scored.values())))
            cells = ""
            for name in scores:
                auc = xs.discrimination(scored[name])[0]
                cells += f"{auc:>16.3f}" if auc == auc else f"{'—':>16}"
            print(f"{kind.value:<22}{total:>9,}{agreeing:>10,}{cells}")
        print("\nAUC = P(an agreeing donor scores above a disagreeing one); 0.5 is no information.")
        print("⛔ Same pairs, same labels, only the similarity differs — which a per-map AUC cannot claim.")
        print("⛔ It says whether a cut COULD exist, never where one belongs.")


if __name__ == "__main__":
    main()
