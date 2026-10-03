"""The derived columns, against the published site catalogue.

⭐ **`meta.landing` and `meta.examples` are a real oracle for `interest_score`.** They are `ranking`'s
own output, computed by `render_page` and shipped inside the site catalogue — so reproducing them
exercises `_interest`, `_pfam_concordance`, the payload's rounding and the stability of the sort, all
at once, over every one of the 17,531 loci.

⚠ They are in `data/{species}.json` and **not** in `data/browser/…json`. `render` MUTATES the payload
and the site build is the only caller that serialises afterwards; the export has neither key. Reading
the wrong file would skip these tests rather than fail them, so the fixture names which one it wants.

`display_name` has no such oracle — it exists only inside `app.js` — so it is tested on the
behaviours that file documents, plus a coverage census over the whole catalogue that says how many
loci each fall-through level actually names.
"""

import pytest

from bacatlas_backend.ingest.derive_locus_display import (
    SOURCE_BAKTA_SYMBOL,
    SOURCE_LABEL,
    SOURCE_PFAM_ARCHITECTURE,
    SOURCE_PRODUCT,
    best_product,
    display_name,
    search_text,
    sole_architecture,
)
from bacatlas_backend.ingest.derive_locus_ranking import (
    BARRED_SENTINEL,
    build_interest_inputs,
    interest_score,
    landing_index,
    pfam_concordance,
    ranking,
    similarity_index,
)
from tests.conftest import PUBLISHED_SITE_CATALOGUE_DIR
from tests.payload_oracle import load_catalogue

pytest.importorskip("pandas", reason="the ingest extra")

#: A minimal stand-in for the vendored table: accession → [short, description, IPR, IPR name, clan, …].
PFAM_NAMES = {
    "PF04932": ["Wzy_C", "O-Antigen ligase", "IPR007016", "O-antigen ligase", "", ""],
    "PF01464": ["SLT", "Transglycosylase SLT domain", "", "", "CL0037", "Lysozyme"],
    "PF00593": ["TonB_dep_Rec", "TonB dependent receptor", "", "", "CL0193", ""],
}


# ── display_name ───────────────────────────────────────────────────────────────────────────────
def test_a_real_bakta_symbol_wins_outright_and_says_so():
    result = display_name("3510", "wzi", [("PF04932", 90)], [("Surface protein", 90)], PFAM_NAMES)
    assert (result.name, result.source, result.source_accession) == ("wzi", SOURCE_BAKTA_SYMBOL, None)


def test_a_single_accession_architecture_names_the_locus_from_its_pfam_SHORT_NAME():
    """`Caps_assemb_Wzi -> Wzi`: the last underscore token, if it has the shape of a gene symbol."""
    names = {**PFAM_NAMES, "PF10670": ["Caps_assemb_Wzi", "Capsule assembly", "", "", "", ""]}
    result = display_name("3510", None, [("PF10670", 90)], [], names)
    assert (result.name, result.source, result.source_accession) == ("Wzi", SOURCE_PFAM_ARCHITECTURE, "PF10670")


def test_a_TWO_accession_architecture_names_nothing_because_it_names_no_gene():
    """⚠ `soleArch` is about the ARCHITECTURE having one accession, not the locus having one architecture."""
    assert sole_architecture([("PF04932,PF01464", 90)], PFAM_NAMES) is None


def test_an_architecture_whose_accession_is_not_in_the_table_names_nothing():
    """⛔ This is the `pfam_names` trap: the table is attached at RENDER time and is not in the export."""
    assert sole_architecture([("PF99999", 90)], PFAM_NAMES) is None


def test_a_short_name_whose_last_token_is_not_gene_shaped_falls_through():
    """⚠ `TonB_dep_Rec` DOES pass — `Rec` is gene-shaped — so the negative case needs a real one."""
    names = {**PFAM_NAMES, "PF00005": ["ABC_tran", "ABC transporter", "", "", "CL0023", ""]}
    result = display_name("3510", None, [("PF00005", 5)], [], names)
    assert result.source == SOURCE_LABEL
    assert display_name("3510", None, [("PF00593", 5)], [], PFAM_NAMES).name == "Rec"


def test_the_product_level_takes_a_symbol_named_before_family_or_like():
    result = display_name("3510", None, [], [("LysR family transcriptional regulator", 40)], PFAM_NAMES)
    assert (result.name, result.source) == ("LysR", SOURCE_PRODUCT)
    assert display_name("77", None, [], [("MarR-like protein", 8)], PFAM_NAMES).name == "MarR"


def test_a_locus_with_nothing_to_go_on_is_named_by_its_LABEL_and_never_left_NULL():
    result = display_name("10515", None, [], [], PFAM_NAMES)
    assert (result.name, result.source) == ("locus 10515", SOURCE_LABEL)


# ── best_product ───────────────────────────────────────────────────────────────────────────────
def test_the_best_product_EXTENDS_the_modal_rather_than_being_the_commonest():
    """The `rfaL` case: modal *Ligase*, and the same locus holds *O-antigen ligase RfaL*."""
    entries = [("Ligase", 60), ("O-antigen ligase RfaL", 30), ("Ligase", 0)]
    assert best_product("rfaL", entries[:2]) == "O-antigen ligase RfaL"


def test_a_candidate_that_would_become_a_DIFFERENT_protein_cannot_win():
    """Plain longest-wins picked *L-arabinose isomerase* for `araB`. It neither extends nor names."""
    entries = [("Ribulokinase", 50), ("L-arabinose isomerase", 40)]
    assert best_product("araB", entries) == "Ribulokinase"


def test_a_candidate_below_a_tenth_of_the_LISTED_counts_is_ignored():
    entries = [("Ligase", 100), ("O-antigen ligase RfaL", 5)]
    assert best_product("rfaL", entries) == "Ligase"


def test_a_vague_candidate_may_never_win_but_a_vague_MODAL_may_be_replaced():
    """⛔ The asymmetry is the point: swapping a vague modal for something specific is why this exists."""
    assert best_product(None, [("Wzi", 50), ("hypothetical protein Wzi", 30)]) == "Wzi"
    assert (
        best_product(None, [("hypothetical protein", 50), ("hypothetical protein, capsule assembly", 30)])
        == "hypothetical protein, capsule assembly"
    )


def test_a_case_only_variant_is_churn_and_never_a_swap():
    assert best_product(None, [("Ligase", 50), ("ligase", 40)]) == "Ligase"


def test_a_locus_with_no_products_has_no_best_product_rather_than_an_empty_string():
    assert best_product("wzi", []) is None


# ── search_text ────────────────────────────────────────────────────────────────────────────────
def test_the_haystack_is_symbol_label_products_and_uniref_and_NOT_pfam():
    """⛔ Adding Pfam/COG/GO/EC/KEGG is a product change, not an implementation detail."""
    text = search_text("3510", "wzi", [("Surface assembly protein Wzi", 90)], ["P30979", "A0A1"])
    assert text == "wzi 3510 surface assembly protein wzi p30979 a0a1"
    assert "pf04932" not in text


def test_a_locus_with_no_symbol_still_carries_its_label_so_it_stays_findable():
    assert search_text("10515", None, [], []).split() == ["10515"]


# ── the published-catalogue oracle ─────────────────────────────────────────────────────────────
@pytest.fixture(scope="module")
def scored_ecoli():
    """`(catalogue, scores)` for the published *E. coli* site catalogue."""
    path = PUBLISHED_SITE_CATALOGUE_DIR / "ecoli.json"
    if not path.exists():
        pytest.skip(f"the published site catalogue is not present at {path}")
    catalogue = load_catalogue(path)
    n = catalogue.n_loci
    architecture_counts = [
        [row.gene_count for row in catalogue.annotation_rows("pfam", i)] for i in range(n)
    ]
    concordance = pfam_concordance(
        architecture_counts, catalogue.nodes["n_pfam"], catalogue.nodes["size"]
    )
    rows = build_interest_inputs(
        syntenic_a5=catalogue.nodes["a5"],
        uniref50_major_family_count=catalogue.nodes["n_u50_major"],
        collapse_tier=[catalogue.string("tier", t) for t in catalogue.nodes["tier"]],
        resolved_threshold=catalogue.nodes["resolved"],
        prevalence_band_index=catalogue.nodes["band"],
        member_gene_count=catalogue.nodes["size"],
        concordance=concordance,
        pfam_concordance_class=[catalogue.string("pfclass", c) for c in catalogue.nodes["pfclass"]],
    )
    return catalogue, [interest_score(row, catalogue.meta["policy"]) for row in rows]


def test_the_ranking_reproduces_the_published_example_chips_over_the_WHOLE_catalogue(scored_ecoli):
    catalogue, scores = scored_ecoli
    assert len(scores) == catalogue.n_loci == 17_531, "coverage, stated: every locus was scored"
    assert ranking(scores)[:4] == catalogue.meta["examples"]


def test_the_landing_locus_reproduces_the_published_one_and_is_NOT_the_first_chip(scored_ecoli):
    """`landing` FILTERS the ranking for `fimA`, and on this catalogue that preference fires."""
    catalogue, scores = scored_ecoli
    order = ranking(scores)
    symbols = [catalogue.string("sym", value) for value in catalogue.nodes["name"]]
    landing = landing_index(order, symbols)
    assert landing == catalogue.meta["landing"] == 2811
    assert landing != catalogue.meta["examples"][0]
    assert str(symbols[landing]).lower() == "fima"


def test_negative_does_NOT_mean_barred_and_both_counts_are_measured(scored_ecoli):
    """⛔ 51 loci carry the sentinel; 912 more are negative and fully rankable."""
    _, scores = scored_ecoli
    assert sum(1 for score in scores if score == BARRED_SENTINEL) == 51
    assert sum(1 for score in scores if score < 0 and score != BARRED_SENTINEL) == 912


def _measurable(catalogue, representation: str) -> tuple[list[float | None], int]:
    """The published catalogue's separation percentiles, ranked with the helper the ingest uses.

    ⛔ **The MEASURE changed on 2026-09-24; the measurable SET did not — and this helper had to be
    re-pointed on 2026-10-03 to keep saying so.** It used to read `nodes.{esm,bac}_d_{intra,near}`,
    the retired medoid geometry, and its own docstring called those "frozen inside the published site
    catalogue — the parity oracle, which does not move". The oracle DID move: both pages were
    re-exported at schema 16, those four columns no longer exist, and this raised `KeyError:
    'esm_d_intra'`. The numbers it asserts are unchanged — 12,104 / 5,427 / 5,458 / 31 all hold
    against the `sim` block, because a singleton has no within-locus pair under either construction.

    ⛔ **`raw["sim"]`, never `raw.get("sim", {})`.** A missing block must raise here. `payload_oracle`
    already degrades silently on the two it lost (`.get("map_reps", [])`, `.get("null", {})`), and an
    accessor that returns empty would make `unmeasurable` the empty set — which passes the `<=`
    subset assertion below trivially. *Not looked at* must not be spellable as *no difference*.
    """
    block = catalogue.raw["sim"][representation]
    index = similarity_index(
        block["within"], block["near"], block["weak_in"], block["weak_out"], block["own"]
    )
    return index.separation_percentile, int(block["n_measurable"])


def test_the_separation_denominator_is_the_number_the_card_prints(scored_ecoli):
    """*"p12 of 12,104 loci"* — and the 5,427 singletons are NULL, never 0.000."""
    catalogue, _ = scored_ecoli
    for representation in ("esm", "bacformer"):
        percentile, declared = _measurable(catalogue, representation)
        measurable = sum(1 for value in percentile if value is not None)
        assert measurable == 12_104, representation
        assert sum(1 for value in percentile if value is None) == catalogue.n_loci - 12_104
        assert all(0.0 <= value <= 1.0 for value in percentile if value is not None)
        # ⭐ The block DECLARES its own measurable count; the ranking MEASURES one. Checking the two
        # against each other is free and is the only thing that catches a payload whose declaration
        # and whose columns describe different catalogues.
        assert declared == measurable, representation


def test_the_measurable_set_is_the_MEASUREMENT_and_not_the_prevalence_band(scored_ecoli):
    """⚠ `RARE` covers 5,458 loci and only 5,427 are unmeasurable — the extra 31 are paralogues
    inside one genome, which DO have a within-locus pair. Gating on the band would blank all 31."""
    catalogue, _ = scored_ecoli
    percentile, _declared = _measurable(catalogue, "esm")
    rare = catalogue.meta["bands"].index("rare")
    unmeasurable = {i for i, value in enumerate(percentile) if value is None}
    rare_loci = {i for i, band in enumerate(catalogue.nodes["band"]) if band == rare}
    assert len(unmeasurable) == 5_427
    assert len(rare_loci) == 5_458
    assert len(rare_loci - unmeasurable) == 31
    assert unmeasurable <= rare_loci
