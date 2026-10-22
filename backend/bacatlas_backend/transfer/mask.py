"""The held-out mask — ⛔ `L_open` is the only label a transfer step may read.

Design: `nuna/docs/model_evaluation/annotation_transfer_design.md` §2.3b step 1.

    L_open(u) = u's claim   if u carries it and u ∉ M
              = ∅           otherwise  (unannotated, or masked)

`M` is a random **20 %** of the nodes **carrying the claim**, seed 0, drawn per claim.

⛔⛔ **This is the one error that would make every number in this strand meaningless, and it does not
fail — it produces a better number.** A step that reads a masked node's real claim scores the model
against something it was shown; the accuracy goes up and nothing anywhere raises. So the two labels are
separate methods with separate names, and the truth accessor **refuses** every node it has no business
answering for:

======================================= =========================================================
:meth:`HeldOutMask.open_label`          every step. ∅ for masked **and** for unannotated nodes
:meth:`HeldOutMask.held_out_truth`      ⛔ the scorer alone, and only for a held-out node
======================================= =========================================================

⭐ **A masked node is deliberately indistinguishable from an unannotated one.** That is the whole
construction: the held-out nodes stand in for the nodes we actually want to annotate, so any step able
to tell the two apart is measuring something easier than the real problem.

⚠ **`held_out_truth` also refuses a node that has been propagated onto** — the design's named rule.
Once a predicted label has been written to a node, its "truth" is a prediction, and scoring against it
compares the model with itself and cannot fail. Call :meth:`note_propagated` as labels are written; it
records that a prediction was made and does **not** make it readable.

⚠ **Takes a plain claim map, not a `Session`.** `claims.Catalogue.claims[key]` is already
`{locus_id: frozenset[str]}`, so this module needs no database and its tests need no fixture — which
is why they could be written first.
"""

from __future__ import annotations

import numpy as np

#: The design's settings. ⛔ Named constants because a mask fraction passed positionally at one of
#: thirteen call sites, differing from the others, is invisible in every output.
HELD_OUT_FRACTION = 0.2
SEED = 0

#: Below this many held-out nodes a claim is not scoreable and the mask refuses to be built. ⛔ Not a
#: calibrated threshold — the smallest number for which a rate is arithmetically defined at all.
MIN_HELD_OUT = 1


class LabelLeak(Exception):
    """⛔ Raised when truth was asked for about a node that has no business having truth read.

    Its own class, not `ValueError`, because this is never a bad argument to be handled — it is a
    measurement about to be wrong, and nothing should be catching it.
    """


class HeldOutMask:
    """One claim's open labels and its held-out truth, kept apart on purpose.

    Build with :meth:`build`. The instance is otherwise read-only except for
    :meth:`note_propagated`, which only ever narrows what truth will answer.
    """

    __slots__ = ("_claims", "_held_out", "_propagated", "fraction", "seed")

    def __init__(self, claims: dict[int, frozenset[str]], held_out: frozenset[int], *, fraction: float, seed: int):
        self._claims = claims
        self._held_out = held_out
        self._propagated: set[int] = set()
        self.fraction = fraction
        self.seed = seed

    # ------------------------------------------------------------------------------------ building
    @classmethod
    def build(
        cls,
        claims: dict[int, frozenset[str]],
        *,
        fraction: float = HELD_OUT_FRACTION,
        seed: int = SEED,
    ) -> HeldOutMask:
        """Draw the held-out set from the nodes that **carry** the claim.

        ⛔ **Sorted before the draw.** The claim map arrives from a SQL query, whose row order is not
        guaranteed, so a draw that followed insertion order would hold out different nodes on two runs
        with the same seed — and every accuracy figure would be unreproducible while looking fine.
        Sorting makes the mask a function of the node ids and the seed, and of nothing else.

        ⛔ **Drawn from the carrying nodes only.** Masking an unannotated node holds back nothing and
        puts a row with no truth into the held-out denominator.
        """
        if not 0.0 < fraction < 1.0:
            raise ValueError(f"fraction must lie in (0, 1), got {fraction}")
        carrying = np.array(sorted(claims), dtype=np.int64)
        n_hold = int(round(len(carrying) * fraction))
        if n_hold < MIN_HELD_OUT:
            raise ValueError(
                f"too few nodes carry this claim to hold any back: {len(carrying)} carrying x "
                f"{fraction} rounds to {n_hold}. A mask of nothing makes every rate 0/0, which prints "
                "as a dash rather than as 'this claim was not scoreable' — so it is refused here."
            )
        rng = np.random.default_rng(seed)
        chosen = rng.choice(carrying, size=n_hold, replace=False)
        return cls(dict(claims), frozenset(int(node) for node in chosen), fraction=fraction, seed=seed)

    # -------------------------------------------------------------------------------- what it holds
    @property
    def held_out(self) -> frozenset[int]:
        """`M` — the nodes whose claim is hidden. ⚠ Membership is readable; the labels are not."""
        return self._held_out

    @property
    def n_carrying(self) -> int:
        """Nodes stating this claim at all — the population the mask was drawn from."""
        return len(self._claims)

    @property
    def n_open(self) -> int:
        """Nodes whose claim a transfer step may read."""
        return len(self._claims) - len(self._held_out)

    # -------------------------------------------------------------------------- L_open, for everyone
    def open_label(self, node: int) -> frozenset[str]:
        """`L_open(node)` — ⭐ **the only label any transfer step may read.**

        ∅ for a masked node, an unannotated one, and one this catalogue has never heard of. Those
        three are one answer on purpose: a step that could tell them apart would be solving an easier
        problem than the one being measured.
        """
        if node in self._held_out:
            return frozenset()
        return self._claims.get(node, frozenset())

    def is_held_out(self, node: int) -> bool:
        """Whether `node` is in `M`. ⚠ For the scorer's bookkeeping, never as a feature."""
        return node in self._held_out

    # ------------------------------------------------------------------ L_true, for the scorer alone
    def held_out_truth(self, node: int) -> frozenset[str]:
        """`L_true(node)` — ⛔ **the scorer's accessor, and nobody else's.**

        Refuses three populations, each because scoring on it would raise the measured accuracy
        without raising anything else:

        * **not held out** — the model was shown this label, so scoring on it is leakage;
        * **carries no claim** — there is no truth, and the empty set would score a prediction as
          wrong against nothing;
        * **propagated** — a predicted label has been written here, so its "truth" is a prediction
          and the comparison cannot fail.
        """
        if node in self._propagated:
            raise LabelLeak(
                f"node {node} has been propagated onto, so its label is a prediction — scoring "
                "against it compares the model with itself and cannot fail"
            )
        if node not in self._held_out:
            carried = node in self._claims
            raise LabelLeak(
                f"node {node} is not held out "
                + (
                    "and carries no claim, so there is no truth to score against"
                    if not carried
                    else "— the model could read this label, so scoring on it is leakage"
                )
            )
        if node not in self._claims:  # pragma: no cover — build() draws only from carrying nodes
            raise LabelLeak(f"node {node} is held out but carries no claim; the mask is inconsistent")
        return self._claims[node]

    def note_propagated(self, node: int) -> None:
        """Record that a predicted label was written to `node`.

        ⚠ This narrows what :meth:`held_out_truth` will answer and changes nothing else — in
        particular it does **not** make the prediction readable through :meth:`open_label`, because a
        step that could read its own output back as evidence would compound its own errors.
        """
        self._propagated.add(int(node))

    @property
    def n_propagated(self) -> int:
        """How many nodes have had a predicted label written to them."""
        return len(self._propagated)
