"""Pair evidence — ⛔ both covariates must describe the SAME node pair.

Design §1.3. The defect this strand exists to fix put one pair's ESM cosine beside a *different*
pair's Bacformer cosine, because the two representations name the same first annotated donor only
29-38 % of the time (S1). A join that fills a missing side, or that silently drops it, reproduces that
defect with a different mechanism — so the join is outer, the misses are counted, and nothing is imputed.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from bacatlas_backend.transfer.evidence import PairEvidence, pair_evidence


def table(rows) -> pd.DataFrame:
    """A neighbour artifact's four relevant columns."""
    return pd.DataFrame(rows, columns=["recipient_node", "donor_node", "cross_max", "rank_cross_max"])


ESM = table([("r1", "d1", 0.99, 1), ("r1", "d2", 0.95, 2), ("r2", "d3", 0.90, 1)])
BAC = table([("r1", "d1", 0.70, 1), ("r1", "d9", 0.65, 2), ("r2", "d3", 0.55, 1)])


def test_a_pair_in_both_representations_carries_both_cosines():
    ev = pair_evidence(ESM, BAC)
    both = ev.pairs[ev.pairs.in_both]
    z = both.set_index(["recipient_node", "donor_node"])
    assert z.loc[("r1", "d1"), "esm_max"] == pytest.approx(0.99)
    assert z.loc[("r1", "d1"), "bac_max"] == pytest.approx(0.70)
    assert z.loc[("r2", "d3"), "bac_max"] == pytest.approx(0.55)


def test_the_join_is_OUTER_so_a_one_sided_pair_survives_and_is_visible():
    """⛔ Dropping it would hide the disagreement; the 29-38 % concordance is the subject, not noise."""
    ev = pair_evidence(ESM, BAC)
    seen = set(zip(ev.pairs.recipient_node, ev.pairs.donor_node, strict=True))
    assert ("r1", "d2") in seen, "ESM-only pair must survive"
    assert ("r1", "d9") in seen, "Bacformer-only pair must survive"
    assert len(ev.pairs) == 4


def test_a_MISSING_cosine_is_NaN_and_is_never_imputed():
    """⛔ A zero, a mean or a copy of the other representation would all be a fabricated covariate —
    and would read as a confident low similarity rather than as an absent measurement."""
    ev = pair_evidence(ESM, BAC)
    z = ev.pairs.set_index(["recipient_node", "donor_node"])
    assert np.isnan(z.loc[("r1", "d2"), "bac_max"])
    assert np.isnan(z.loc[("r1", "d9"), "esm_max"])
    assert not z["esm_max"].eq(0).any() and not z["bac_max"].eq(0).any()


def test_the_fittable_set_is_the_pairs_present_in_BOTH():
    """⭐ The model can only be fitted where both covariates describe this pair, so that population is
    named once and counted, rather than each caller deciding."""
    ev = pair_evidence(ESM, BAC)
    f = ev.fittable()
    assert len(f) == 2
    assert f[["esm_max", "bac_max"]].notna().all().all()


def test_the_concordance_IS_REPORTED_not_left_to_be_rediscovered():
    """⛔ 'A comparison must assert its own coverage BEFORE it reports a difference' (CLAUDE.md)."""
    ev = pair_evidence(ESM, BAC)
    assert ev.report["n_pairs"] == 4
    assert ev.report["n_both"] == 2
    assert ev.report["n_esm_only"] == 1
    assert ev.report["n_bacformer_only"] == 1
    assert ev.report["frac_both"] == pytest.approx(0.5)


def test_rank_1_agreement_is_reported_per_recipient():
    """⭐ S1's finding in one number: do the two representations name the SAME nearest donor? Here r1
    disagrees (d1 vs d1 — same) and r2 agrees, so 2 of 2. The column exists so a run cannot quote an
    accuracy without this beside it."""
    ev = pair_evidence(ESM, BAC)
    assert ev.report["n_recipients_both_rank1"] == 2
    assert ev.report["n_recipients_rank1_agrees"] == 2
    bac2 = table([("r1", "d2", 0.70, 1), ("r2", "d3", 0.55, 1)])
    assert pair_evidence(ESM, bac2).report["n_recipients_rank1_agrees"] == 1


def test_an_empty_representation_is_refused_rather_than_returning_an_empty_fit():
    with pytest.raises(ValueError, match="no rows"):
        pair_evidence(ESM, table([]))


def test_a_duplicated_pair_is_refused_because_it_would_silently_double_a_row():
    """⚠ The artifacts are one row per pair by construction; if that ever breaks, the join fans out and
    the fitted set is quietly reweighted."""
    dup = pd.concat([ESM, ESM.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="duplicate"):
        pair_evidence(dup, BAC)


def test_the_evidence_is_a_frame_with_the_declared_columns_in_order():
    assert list(pair_evidence(ESM, BAC).pairs.columns) == [
        "recipient_node", "donor_node", "esm_max", "bac_max", "rank_esm", "rank_bacformer", "in_both",
    ]  # fmt: skip


def test_pair_evidence_is_a_dataclass_carrying_its_own_provenance():
    ev = pair_evidence(ESM, BAC, label="ecoli_to_kp")
    assert isinstance(ev, PairEvidence)
    assert ev.label == "ecoli_to_kp"
    assert ev.report["label"] == "ecoli_to_kp"
