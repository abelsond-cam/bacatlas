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
    SOURCE_CORROBORATED_SYMBOL,
    SOURCE_LABEL,
    SOURCE_PFAM_ARCHITECTURE,
    SOURCE_PFAM_PRODUCT_AGREEMENT,
    SOURCE_PRODUCT,
    best_product,
    display_name,
    names_in_text,
    search_text,
    sole_architecture,
    spelled_names,
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

#: No symbol vote at all — what every `display_name` call below that is not ABOUT the vote passes.
NO_VOTE = {
    "symbol_entries": [], "member_gene_count": 0, "named_member_count": 0, "cog_entries": [],
    "cog_names": {},
}


# ── display_name ───────────────────────────────────────────────────────────────────────────────
def test_a_real_bakta_symbol_wins_outright_and_says_so():
    result = display_name(
        "3510", "wzi", [("PF04932", 90)], [("Surface protein", 90)], PFAM_NAMES,
        symbol_entries=[("wzi", 90)], member_gene_count=90, named_member_count=90, cog_entries=[],
        cog_names={},
    )
    assert (result.name, result.source, result.source_accession) == ("wzi", SOURCE_BAKTA_SYMBOL, None)


def test_a_single_accession_architecture_names_the_locus_from_its_pfam_SHORT_NAME():
    """`Caps_assemb_Wzi -> Wzi`: the last underscore token, if it has the shape of a gene symbol."""
    names = {**PFAM_NAMES, "PF10670": ["Caps_assemb_Wzi", "Capsule assembly", "", "", "", ""]}
    result = display_name("3510", None, [("PF10670", 90)], [], names, **NO_VOTE)
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
    result = display_name("3510", None, [("PF00005", 5)], [], names, **NO_VOTE)
    assert result.source == SOURCE_LABEL
    assert display_name("3510", None, [("PF00593", 5)], [], PFAM_NAMES, **NO_VOTE).name == "Rec"


def test_the_product_level_takes_a_symbol_named_before_family_or_like():
    result = display_name(
        "3510", None, [], [("LysR family transcriptional regulator", 40)], PFAM_NAMES, **NO_VOTE
    )
    assert (result.name, result.source) == ("LysR", SOURCE_PRODUCT)
    assert display_name("77", None, [], [("MarR-like protein", 8)], PFAM_NAMES, **NO_VOTE).name == "MarR"


def test_a_locus_with_nothing_to_go_on_is_named_by_its_LABEL_and_never_left_NULL():
    result = display_name("10515", None, [], [], PFAM_NAMES, **NO_VOTE)
    assert (result.name, result.source) == ("locus 10515", SOURCE_LABEL)


# ── a thin symbol vote, settled by the locus's other annotation ────────────────────────────────
#: kp-probe100 locus 722 as `bacatlas_ingest` holds it: Pfam names BOTH `Wzz` and `Wzc` on every gene.
CPS_PFAM = {
    **PFAM_NAMES,
    "PF02706": ["Wzz", "Chain length determinant protein", "IPR003856", "", "", ""],
    "PF23607": ["WZC_N", "Tyrosine-protein kinase Wzc, N-terminal domain", "", "", "CL0841", "Wzc_PBSD"],
    "PF13807": ["GNVR", "G-rich domain on putative tyrosine kinase", "", "", "", ""],
    "PF13614": ["AAA_31", "AAA domain", "", "", "CL0023", ""],
    "PF01656": ["CbiA", "CobQ/CobB/MinD/ParA nucleotide binding domain", "", "", "CL0023", ""],
}
CPS_ARCHITECTURES = [
    ("PF02706,PF23607,PF13807,PF13614,PF01656", 64),
    ("PF02706,PF23607,PF13807,PF13614", 35),
    ("PF02706,PF23607,PF13807", 1),
]
CPS_PRODUCTS = [("Tyrosine autokinase", 30), ("Polysaccharide biosynthesis tyrosine autokinase", 29)]
CPS_COG = {"COG3206": "Exopolysaccharide export protein/domain GumC/Wzc1"}


def _cps(symbols, size=100):
    return display_name(
        "722", symbols[0][0] if symbols else None, CPS_ARCHITECTURES, CPS_PRODUCTS, CPS_PFAM,
        symbol_entries=symbols, member_gene_count=size,
        named_member_count=sum(count for _, count in symbols), cog_entries=[("COG3206", 1)],
        cog_names=CPS_COG,
    )


def test_the_gumC_tie_goes_to_wzc_which_every_member_genes_Pfam_names():
    """⭐ kp 722: `gumC` ×1, `wzc` ×1 of 100, and the alphabet used to pick `gumC`.

    COG3206 names both, one gene each; PF23607 names `wzc` on all 100 — and the accession that did
    is what the page cites.
    """
    result = _cps([("gumC", 1), ("wzc", 1)])
    assert (result.name, result.source, result.source_accession) == (
        "wzc", SOURCE_CORROBORATED_SYMBOL, "PF23607"
    )


def test_a_name_no_member_gene_carries_is_NEVER_proposed_however_well_Pfam_supports_it():
    """⛔ `Wzz` and `Wzc` both sit on all 100 genes. With only `gumC` in the vote, `gumC` stands.

    Pfam names a family by its founder, so letting it PROPOSE a name was measured to rename `glpR`
    `deoR` and `etk` `wzc` — it may only choose between the names the genes carry. ⚠ The one
    exception is `pfam_product_agreement` (under a twentieth named, and the PRODUCTS spell the name
    too, on half the genes); these products (*"Tyrosine autokinase"*) spell neither.
    """
    result = _cps([("gumC", 1)])
    assert (result.name, result.source, result.source_accession) == ("gumC", SOURCE_BAKTA_SYMBOL, None)


def test_a_DECISIVE_vote_is_never_overridden_by_a_family_name():
    """`yhaJ` 95 / `lysR` 5 with every product saying *LysR family* — measured, and `yhaJ` must stand."""
    result = display_name(
        "26", "yhaJ", [], [("LysR family transcriptional regulator", 100)], PFAM_NAMES,
        symbol_entries=[("yhaJ", 95), ("lysR", 5)], member_gene_count=100, named_member_count=100,
        cog_entries=[], cog_names={},
    )
    assert (result.name, result.source) == ("yhaJ", SOURCE_BAKTA_SYMBOL)


def test_the_vote_is_thin_only_STRICTLY_below_a_fifth_of_the_member_genes():
    """2 of 10 is exactly the threshold and decides; 2 of 11 is below it and does not."""
    assert _cps([("gumC", 2), ("wzc", 2)], size=10).source == SOURCE_BAKTA_SYMBOL
    assert _cps([("gumC", 2), ("wzc", 2)], size=11).name == "wzc"


def test_a_thin_vote_nothing_else_names_keeps_its_commonest_symbol():
    """With nothing else naming either symbol, each scores its own count, so the commonest stands."""
    result = display_name(
        "9", "abcD", [], [("Hypothetical protein", 3)], PFAM_NAMES,
        symbol_entries=[("abcD", 2), ("xyzE", 1)], member_gene_count=50, named_member_count=3,
        cog_entries=[], cog_names={},
    )
    assert (result.name, result.source) == ("abcD", SOURCE_BAKTA_SYMBOL)


def test_the_modal_symbol_and_the_listed_vote_must_be_one_vote():
    """The gate and the candidates read the listed vote; the fallback reads the modal column."""
    with pytest.raises(ValueError, match="one vote"):
        display_name(
            "722", "gumC", [], [], PFAM_NAMES,
            symbol_entries=[("wzc", 1), ("gumC", 1)], member_gene_count=100, named_member_count=2,
            cog_entries=[], cog_names={},
        )


# ── almost no gene named: the name Pfam and the products agree on ──────────────────────────────
#: kp-71k v0.32 locus 4008 as the beta dump holds it (2026-11-24), Pfam names from the vendored table.
WZI_PFAM = {
    **PFAM_NAMES,
    "PF14052": ["Caps_assemb_Wzi", "Capsule assembly protein Wzi", "", "", "", ""],
    "PF06082": ["YjbH", "Exopolysaccharide biosynthesis protein YbjH", "", "", "", ""],
}
WZI_ARCHITECTURES = [("PF14052", 69825), ("PF06082", 10)]
WZI_PRODUCTS = [
    ("capsule assembly Wzi family protein", 63810),
    ("Capsule assembly Wzi family protein", 920),
    ("Surface assembly of capsule", 275),
    ("Uncharacterized protein in cps region", 218),
    ("Wzi", 14),
    ("Lipoprotein", 7),
]
WZI_SYMBOLS = [("yjbH", 2), ("gfcC", 1)]


def _wzi(named=3, size=69884, products=WZI_PRODUCTS, architectures=WZI_ARCHITECTURES, symbols=WZI_SYMBOLS):
    return display_name(
        "4008", symbols[0][0] if symbols else None, architectures, products, WZI_PFAM,
        symbol_entries=symbols, member_gene_count=size, named_member_count=named, cog_entries=[],
        cog_names={},
    )


def test_three_named_genes_of_69884_do_not_name_the_locus_PF14052_and_the_products_say_Wzi():
    """⭐ kp-71k 4008 was `yjbH` on two genes. Pfam PF14052 and the products agree on `Wzi`."""
    result = _wzi()
    assert (result.name, result.source, result.source_accession) == (
        "Wzi", SOURCE_PFAM_PRODUCT_AGREEMENT, "PF14052"
    )


def test_the_rule_fires_only_STRICTLY_below_a_twentieth_of_the_genes_named():
    """5 of 100 named is the threshold and the Bakta vote stands (thin, so corroborated); 4 is below."""
    symbols = [("yjbH", 3), ("gfcC", 2)]
    at = _wzi(named=5, size=100, symbols=symbols,
              architectures=[("PF14052", 100)], products=[("Capsule assembly protein Wzi", 100)])
    assert at.source != SOURCE_PFAM_PRODUCT_AGREEMENT
    below = _wzi(named=4, size=100, symbols=[("yjbH", 2), ("gfcC", 2)],
                 architectures=[("PF14052", 100)], products=[("Capsule assembly protein Wzi", 100)])
    assert (below.name, below.source) == ("Wzi", SOURCE_PFAM_PRODUCT_AGREEMENT)


def test_the_agreed_name_must_be_on_HALF_the_genes_in_EACH_source():
    """Products spelling `Wzi` on 49 of 100 genes are not agreement enough; on 50 they are."""
    architectures = [("PF14052", 100)]
    short = _wzi(size=100, architectures=architectures,
                 products=[("capsule assembly Wzi family protein", 49), ("Lipoprotein", 51)])
    assert (short.name, short.source) == ("yjbH", SOURCE_BAKTA_SYMBOL)
    enough = _wzi(size=100, architectures=architectures,
                  products=[("capsule assembly Wzi family protein", 50), ("Lipoprotein", 50)])
    assert enough.name == "Wzi"
    thin_pfam = _wzi(size=100, architectures=[("PF14052", 49)],
                     products=[("capsule assembly Wzi family protein", 100)])
    assert thin_pfam.name == "yjbH"


def test_a_name_only_ONE_source_spells_never_names_the_locus_this_way():
    """Pfam alone says `Wzi` on every gene; the products do not — the Bakta vote stands."""
    result = _wzi(products=[("Lipoprotein", 69884)])
    assert (result.name, result.source) == ("yjbH", SOURCE_BAKTA_SYMBOL)


def test_a_locus_with_NO_named_gene_takes_the_agreed_name_before_the_pfam_short_name():
    """`GerE` is PF00196's short name; Pfam's description and the products both say LuxR."""
    names = {**PFAM_NAMES, "PF00196": ["GerE", "Bacterial regulatory proteins, luxR family", "", "", "", ""]}
    result = display_name(
        "4279", None, [("PF00196", 90)], [("LuxR family transcriptional regulator", 80)], names,
        symbol_entries=[], member_gene_count=100, named_member_count=0, cog_entries=[], cog_names={},
    )
    assert (result.name, result.source, result.source_accession) == (
        "LuxR", SOURCE_PFAM_PRODUCT_AGREEMENT, "PF00196"
    )


def test_the_agreed_name_is_shown_in_PROTEIN_case_so_it_never_reads_as_a_bakta_symbol():
    """The products spell it `luxR`; the page says `LuxR`, like every name it infers."""
    names = {**PFAM_NAMES, "PF00196": ["GerE", "Bacterial regulatory proteins, luxR family", "", "", "", ""]}
    result = display_name(
        "4279", None, [("PF00196", 90)], [("transcriptional regulator, luxR family", 80)], names,
        symbol_entries=[], member_gene_count=100, named_member_count=0, cog_entries=[], cog_names={},
    )
    assert result.name == "LuxR"


def test_a_TonB_dependent_receptor_is_never_named_TonB():
    """⛔ *"TonB-dependent receptor"* needs TonB; it is not TonB. 51 kp-71k loci were, before this."""
    names = {
        **PFAM_NAMES,
        "PF00593": ["TonB_dep_Rec_b-barrel", "TonB dependent receptor-like, beta-barrel", "", "", "", ""],
    }
    result = display_name(
        "7675", "oprC", [("PF00593", 100)], [("TonB-dependent receptor", 95)], names,
        symbol_entries=[("oprC", 2)], member_gene_count=100, named_member_count=2, cog_entries=[],
        cog_names={},
    )
    assert (result.name, result.source) == ("oprC", SOURCE_BAKTA_SYMBOL)
    assert "tonb" not in spelled_names("TonB-dependent receptor")
    assert "tonb" in spelled_names("Gram-negative bacterial TonB protein C-terminal")


def test_spelled_names_skip_lowercase_english_but_keep_symbol_and_protein_case():
    """`box` (*DEAD/DEAH box helicase*) and `and` are gene-shaped and name nothing."""
    assert spelled_names("DEAD/DEAH box helicase and") == {}
    assert spelled_names("Capsule assembly protein Wzi") == {"wzi": "Wzi"}
    assert spelled_names("host specificity factor tipJ") == {"tipj": "tipJ"}
    assert spelled_names("Tle1 phospholipase") == {"tle1": "Tle1", "tle": "Tle"}


def test_a_word_both_sources_capitalise_but_that_names_no_gene_is_refused():
    """`Ion_trans` and *"Ion transport domain-containing protein"*: `Ion` is not a gene (census)."""
    names = {**PFAM_NAMES, "PF00520": ["Ion_trans", "Ion transport protein", "", "", "", ""]}
    result = display_name(
        "942", None, [("PF00520", 100)], [("Ion transport domain-containing protein", 100)], names,
        symbol_entries=[], member_gene_count=100, named_member_count=0, cog_entries=[], cog_names={},
    )
    assert result.source != SOURCE_PFAM_PRODUCT_AGREEMENT


def test_names_in_text_reads_gene_cased_words_and_not_all_caps_ones():
    assert names_in_text("Exopolysaccharide export protein/domain GumC/Wzc1") >= {"gumc", "wzc1", "wzc"}
    assert "abc" not in names_in_text("ABC transporter ATP-binding protein")
    assert names_in_text("WZC N") == frozenset()
    assert names_in_text(None) == frozenset()


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
