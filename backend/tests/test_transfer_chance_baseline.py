"""The chance baseline behind the Function tab's "vs chance" column — pinned as pure functions.

⭐ **Every test here runs WITHOUT a database**, and they are the only ones that do: the DB-backed tests
in `test_function_inference.py` skip when `BACATLAS_DATABASE_URL` is unset (as it is in CI), and a skip
is not a pass. Moved here from `test_cross_species_transfer.py` when the script that file tested was
deleted (it read the per-gene tables the schema dropped); these three test the INSTRUMENT, which is live.

⛔ A wrong null rescales every lift and leaves agreement untouched, so nothing in the output would look
wrong — which is why the pool and the self-exclusion are asserted by arithmetic rather than by eye.
"""

from __future__ import annotations

import pytest

from bacatlas_backend.instruments.annotation_transfer import MIN_PAIRS, _chance_by_rung, tally_cells
from bacatlas_backend.models.enumerations import AnnotationKind

COG = AnnotationKind.COG_ORTHOGROUP


def _cog(accession: str, category: str) -> dict:
    # ⚠ A rung is `(facet, depth)` since 2026-10-06; COG has one axis, so its facet is None.
    return {(None, 2): frozenset({accession}), (None, 1): frozenset({category})}


@pytest.fixture
def two_catalogues():
    """Four donor loci (two hold `E`, two hold `C`) and three recipient loci, all holding `E`.

    ⚠ The shape is chosen so the right null and the wrong one give DIFFERENT answers. With a single
    recipient locus they coincide: adding the focal node to the pool raises numerator and denominator
    by one each, and self-exclusion removes one from each again. Several recipients break the tie.
    """
    donor = {1: _cog("COG0001", "E"), 2: _cog("COG0002", "E"), 3: _cog("COG0003", "C"), 4: _cog("COG0004", "C")}
    recipient = {101: _cog("COG0001", "E"), 102: _cog("COG0001", "E"), 103: _cog("COG0002", "E")}
    return donor, recipient


def test_the_served_null_self_excludes_and_a_donor_pool_does_not(two_catalogues):
    """⛔ The default — what `compute_calibration` serves — is "a random OTHER node of this catalogue".

    5 of the 7 loci hold `E`; the focal node is one of them, so chance is (5−1)/(7−1). A donor pool
    (2 of 4 donor loci hold `E`, the focal node not among them) is exactly 0.5.
    """
    donor, recipient = two_catalogues
    both = {**donor, **recipient}
    edges = [(101, 1, 0.995), (102, 2, 0.995), (103, 1, 0.995)]

    served = tally_cells(both, edges, COG)[(">= 0.99", (None, 1))]
    assert served.chance == pytest.approx((5 - 1) / (7 - 1))

    across = tally_cells(both, edges, COG, donor_pool=donor)[(">= 0.99", (None, 1))]
    assert across.chance == pytest.approx(2 / 4)
    assert across.chance != pytest.approx(served.chance), "the fixture must be able to tell them apart"
    # ⛔ agreement is IDENTICAL under both nulls — which is why a wrong null is invisible
    assert across.agreeing == served.agreeing and across.pairs == served.pairs


def test_self_exclusion_removes_one_from_both_sides_of_the_ratio(two_catalogues):
    """The two modes of `_chance_by_rung`, on one pool, so the arithmetic is visible."""
    donor, _recipient = two_catalogues
    excluded = _chance_by_rung(donor, (None, 1), self_excluded=True)
    included = _chance_by_rung(donor, (None, 1), self_excluded=False)
    assert included(frozenset({"E"})) == pytest.approx(2 / 4)
    assert excluded(frozenset({"E"})) == pytest.approx((2 - 1) / (4 - 1))
    # A claim nothing in the pool holds. Self-excluded it goes NEGATIVE, because the formula assumes the
    # focal node holds its own claim — which every served call does, so the case never reaches a page.
    assert included(frozenset({"Z"})) == 0.0
    assert excluded(frozenset({"Z"})) == pytest.approx(-1 / 3)


def test_a_claim_nobody_holds_has_a_RATE_and_still_no_lift(two_catalogues):
    """Chance 0 ⇒ `lift` is None, never an unbounded ratio — and the guard on chance is what says so.

    ⚠ It must be asserted on a cell with at least `MIN_PAIRS` pairs. Below that `agreement` is None, so
    `lift` is None whatever chance is, and the test would pass with the guard deleted (the version this
    replaces fed it two edges and did exactly that).
    """
    donor, _recipient = two_catalogues
    recipients = {200 + i: _cog("COG9999", "Z") for i in range(MIN_PAIRS)}
    levels = {**donor, **recipients}
    edges = [(locus, 1, 0.995) for locus in recipients]
    cell = tally_cells(levels, edges, COG, donor_pool=donor)[(">= 0.99", (None, 1))]
    assert cell.pairs == MIN_PAIRS
    assert cell.agreement == 0.0, "a measured rate, so the chance guard — not the pair floor — decides lift"
    assert cell.chance == 0.0
    assert cell.lift is None
