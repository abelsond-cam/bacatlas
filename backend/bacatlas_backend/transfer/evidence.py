"""Pair evidence — ⛔ `esm_max` and `bac_max` for the **same** node pair, and the join that proves it.

Design §1.3. ⛔⛔ **This is the defect the whole rebuild exists to fix.** The old model put one pair's
ESM cosine beside a *different* pair's Bacformer cosine: it took ESM's nearest donor and Bacformer's
nearest donor, which S1 measured as **the same node only 29-38 % of the time**, and fitted a surface
over the two numbers as though they described one relationship. They did not.

The fix is not subtle — it is a join on `(recipient_node, donor_node)` — and the discipline is all in
what the join is *allowed to do*:

⛔ **OUTER, never inner.** A pair in ESM's list and absent from Bacformer's is the single most
informative row in the file: it is the disagreement, which is the subject of this work and not noise.
An inner join would delete exactly the evidence that motivated the rebuild.

⛔ **A missing cosine stays NaN.** Zero, the column mean, or a copy of the other representation are all
a **fabricated covariate**, and the first is worst — it reads as a confident *low* similarity rather
than as an absent measurement, which is the opposite of the truth (the pair was near enough for the
other representation to shortlist it).

⛔ **The fittable set is named once.** A model needs both covariates for the pair in front of it, so
it can only be fitted on `in_both` — and how large that population is, and what was left out, is
reported rather than left for a reader to infer from a row count.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

#: The artifact columns this reads. `cross_max` is the settled node-to-node similarity (design §2.1).
NEIGHBOUR_COLUMNS = ("recipient_node", "donor_node", "cross_max", "rank_cross_max")

#: The evidence frame's columns, in order, so a caller never guesses.
EVIDENCE_COLUMNS = (
    "recipient_node",
    "donor_node",
    "esm_max",
    "bac_max",
    "rank_esm",
    "rank_bacformer",
    "in_both",
)


@dataclass
class PairEvidence:
    """One direction's pairs with both similarities, plus ⛔ the coverage the pairs were joined under."""

    pairs: pd.DataFrame
    report: dict = field(default_factory=dict)
    label: str = ""

    def fittable(self) -> pd.DataFrame:
        """⭐ The only rows a model may be fitted on — both covariates describe **this** pair."""
        return self.pairs[self.pairs["in_both"]].reset_index(drop=True)


def _checked(frame: pd.DataFrame, which: str) -> pd.DataFrame:
    missing = [c for c in NEIGHBOUR_COLUMNS if c not in frame.columns]
    if missing:
        raise ValueError(f"{which} table is missing {missing}; expected {list(NEIGHBOUR_COLUMNS)}")
    if not len(frame):
        raise ValueError(
            f"{which} table has no rows — a join against it would return an empty fittable set, and an "
            "empty fit reports a perfect calibration on nothing"
        )
    dup = frame.duplicated(subset=["recipient_node", "donor_node"]).sum()
    if dup:
        raise ValueError(
            f"{which} table has {dup} duplicate (recipient, donor) pairs. The artifacts are one row per "
            "pair by construction; a fan-out here silently reweights the fitted set"
        )
    return frame[list(NEIGHBOUR_COLUMNS)]


def pair_evidence(esm: pd.DataFrame, bacformer: pd.DataFrame, *, label: str = "") -> PairEvidence:
    """Join the two representations' neighbour tables on the **pair** — the one seam that must hold.

    Both inputs are a neighbour artifact's :data:`NEIGHBOUR_COLUMNS` for one direction and one screen.
    ⚠ They are **not** expected to agree on which pairs exist; that disagreement is the measurement.
    """
    left = _checked(esm, "esm").rename(columns={"cross_max": "esm_max", "rank_cross_max": "rank_esm"})
    right = _checked(bacformer, "bacformer").rename(
        columns={"cross_max": "bac_max", "rank_cross_max": "rank_bacformer"}
    )
    merged = left.merge(right, on=["recipient_node", "donor_node"], how="outer", indicator=True)
    merged["in_both"] = merged["_merge"] == "both"

    # ⭐ Rank-1 agreement, reported beside the counts because S1's finding is that the two
    # representations usually name DIFFERENT nearest donors — and an accuracy quoted without it has
    # left out the reason this join exists.
    r1_esm = left[left["rank_esm"] == 1].set_index("recipient_node")["donor_node"]
    r1_bac = right[right["rank_bacformer"] == 1].set_index("recipient_node")["donor_node"]
    shared = r1_esm.index.intersection(r1_bac.index)
    agree = int((r1_esm.reindex(shared) == r1_bac.reindex(shared)).sum())

    counts = merged["_merge"].value_counts()
    n = int(len(merged))
    report = {
        "label": label,
        "n_pairs": n,
        "n_both": int(counts.get("both", 0)),
        "n_esm_only": int(counts.get("left_only", 0)),
        "n_bacformer_only": int(counts.get("right_only", 0)),
        "frac_both": (int(counts.get("both", 0)) / n) if n else float("nan"),
        "n_recipients_both_rank1": int(len(shared)),
        "n_recipients_rank1_agrees": agree,
        "frac_rank1_agrees": (agree / len(shared)) if len(shared) else float("nan"),
    }
    return PairEvidence(pairs=merged[list(EVIDENCE_COLUMNS)].reset_index(drop=True), report=report, label=label)
