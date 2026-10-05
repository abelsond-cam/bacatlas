"""What a gene's GO terms claim — the closure, and the two ways it could go wrong.

⭐ This is David's fourth report, as a test file. On *Klebsiella* `gumC` the page said the members'
cellular-component classes *differ*: 16 genes say `plasma membrane`, one says `membrane`. They do
not differ — one is the parent of the other, and the comparison could not see it, because
`goslim_metagenomics` holds both as separate classes and the comparator is plain set algebra.

The fix has two halves and both can fail silently, in opposite directions:

* **too little** — a closure that stops at the first slim ancestor leaves siblings disjoint;
* **too much** — a closure that keeps the namespace roots makes *everything* agree, which is
  `vendor_go`'s standing objection to full ancestor closure and the reason it was rejected.

Every test below pins one side or the other.
"""

from __future__ import annotations

from bacatlas_backend.instruments.gene_ontology_fold import NAMESPACE_ROOTS, GeneOntologyFold

MEMBRANE, PLASMA, OUTER = "GO:0016020", "GO:0005886", "GO:0019867"
CYTOPLASM, INTRACELLULAR = "GO:0005737", "GO:0005622"
KINASE = "GO:0016301"
CC_ROOT, MF_ROOT = "GO:0005575", "GO:0003674"

#: A toy slim with the production shapes: a parent and two children all in the slim, a second branch
#: that shares nothing with the first, a namespace root, and a raw term that folds onto a class.
FOLD = GeneOntologyFold(
    slim_of={
        MEMBRANE: frozenset({MEMBRANE}),
        PLASMA: frozenset({PLASMA}),
        OUTER: frozenset({OUTER}),
        CYTOPLASM: frozenset({CYTOPLASM}),
        INTRACELLULAR: frozenset({INTRACELLULAR}),
        KINASE: frozenset({KINASE}),
        CC_ROOT: frozenset({CC_ROOT}),
        MF_ROOT: frozenset({MF_ROOT}),
        # A RAW term, which is what within-node propagation hands in.
        "GO:0031226": frozenset({PLASMA}),
        # A raw term that reaches no slim class self-maps, per `vendor_go.slim_of`.
        "GO:0099999": frozenset({"GO:0099999"}),
    },
    ancestors={
        PLASMA: frozenset({MEMBRANE}),
        OUTER: frozenset({MEMBRANE}),
        CYTOPLASM: frozenset({INTRACELLULAR}),
    },
    namespace_of={
        MEMBRANE: 2, PLASMA: 2, OUTER: 2, CYTOPLASM: 2, INTRACELLULAR: 2,
        KINASE: 0, CC_ROOT: 2, MF_ROOT: 0, "GO:0099999": 1,
    },
)  # fmt: skip


# ── the half that was missing ────────────────────────────────────────────────────────────────────
def test_a_class_and_its_PARENT_now_meet():
    """⭐ `gumC` exactly: 16 genes say plasma membrane, one says membrane."""
    assert FOLD.claims([PLASMA]) & FOLD.claims([MEMBRANE]) == {MEMBRANE}


def test_two_SIBLINGS_meet_at_their_shared_parent():
    # Plain set algebra could never see this: neither side states `membrane` at all.
    assert FOLD.claims([PLASMA]) & FOLD.claims([OUTER]) == {MEMBRANE}


def test_a_RAW_accession_folds_before_it_closes():
    # ⛔ The two measurement paths hand in different things — the ladder slim classes, within-node
    # propagation raw accessions — and both called `single_rung`, which wraps whatever it is given.
    # So one compared slim classes and the other compared raw terms, and the page printed both as
    # "agree". This is the one door for either.
    assert FOLD.claims(["GO:0031226"]) == {PLASMA, MEMBRANE}
    assert FOLD.claims(["GO:0031226"]) & FOLD.claims([MEMBRANE]) == {MEMBRANE}


# ── the half that must NOT happen ────────────────────────────────────────────────────────────────
def test_the_closure_does_NOT_collapse_two_different_branches():
    """⛔⛔ `vendor_go`'s standing objection, as a gate.

    Under FULL ancestor closure every cellular-component pair shares `GO:0005575` and `disjoint`
    becomes unreachable. Restricting the closure to the slim and stripping the roots is what keeps
    this assertion possible — and a change that re-admits the roots fails here, not in production.
    """
    assert FOLD.claims([MEMBRANE]) & FOLD.claims([CYTOPLASM]) == frozenset()
    assert FOLD.claims([KINASE]) & FOLD.claims([MEMBRANE]) == frozenset()


def test_a_ROOT_is_never_claimed_even_when_it_is_the_only_term_stated():
    # 452 ecoli loci state nothing but `cellular_component`. Their claim set is EMPTY — "annotated,
    # and stating nothing comparable" — which is neither agreement nor disagreement.
    assert FOLD.claims([CC_ROOT]) == frozenset()
    assert FOLD.claims([CC_ROOT, MF_ROOT]) == frozenset()
    assert not FOLD.claims([PLASMA]) & NAMESPACE_ROOTS


def test_a_root_stated_BESIDE_a_real_class_does_not_leak_in():
    # ⛔ The ordering inside `claims`: roots are dropped BEFORE the closure. Closing first and
    # stripping after is a different operation that only looks the same.
    assert FOLD.claims([PLASMA, CC_ROOT]) == {PLASMA, MEMBRANE}


# ── the states a caller must be able to tell apart ───────────────────────────────────────────────
def test_an_unfoldable_term_is_KEPT_not_dropped():
    # `vendor_go.slim_of`'s rule. Folding it away would make an annotated gene read as unannotated,
    # which is a different and false finding from "annotated with nothing comparable".
    assert FOLD.claims(["GO:0099999"]) == {"GO:0099999"}


def test_empty_and_missing_terms_are_tolerated():
    assert FOLD.claims([]) == frozenset()
    assert FOLD.claims(["", None]) == frozenset()  # type: ignore[list-item]


def test_an_UNLOADED_fold_behaves_AS_BEFORE_rather_than_emptying_the_claim():
    """⛔ A deployment that has not run `ingest --stage reference` must degrade, not break.

    With no reference the claim is the terms as stated, minus the roots — which is what this code
    did before the reference existed. An empty fold that returned an empty claim would silently
    report every GO pair as incomparable and the agreement rate would vanish.
    """
    empty = GeneOntologyFold()
    assert empty.claims([PLASMA, MEMBRANE]) == {PLASMA, MEMBRANE}
    assert empty.claims([PLASMA]) & empty.claims([MEMBRANE]) == frozenset(), "the OLD, wrong answer"
    assert empty.claims([CC_ROOT]) == frozenset(), "the roots go even with no reference"


# ── the namespace split, which the pooled claim conflates ────────────────────────────────────────
def test_the_split_separates_two_answers_to_two_different_questions():
    """⚠ Pooled, two nodes "agree on GO" if they share any class in any namespace.

    A molecular function and a cellular component are not rival answers to one question. This is why
    GO's chance baseline reads far above COG's, and splitting is what the rest of Stage 2 does with
    it.
    """
    split = FOLD.claims_by_namespace([PLASMA, KINASE])
    assert split == {2: frozenset({PLASMA, MEMBRANE}), 0: frozenset({KINASE})}


def test_a_namespace_the_gene_says_nothing_in_is_ABSENT_not_empty():
    # Absent means "not comparable here"; a caller must drop such a pair rather than score it as
    # either agreeing or disagreeing. An empty frozenset would read as the second.
    assert set(FOLD.claims_by_namespace([KINASE])) == {0}
    assert FOLD.claims_by_namespace([CC_ROOT]) == {}
