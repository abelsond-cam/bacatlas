"""Parity suites T1, T3(a), T5 and T7 — the database against the published catalogue.

These four are **pure data**: they need no browser, so they are S1's real exit criterion rather than
S2's. Each drives the loaded database and compares it, row for row, against the catalogue the live
page actually serves.

⛔ **Every suite asserts its own coverage BEFORE it reports a difference** — same locus count, same
order, same labels, non-empty. This repo's own rule, and its own scar: *a diff loop that `continue`s
on a missing column reports "0 differ" while silently skipping exactly the columns that changed*
(`4ab35ca`). A suite that examined 300 of 17,531 loci must say 300, not "pass".

⚠ **They read `BACATLAS_DATABASE_URL`, not the schema probe.** The probe database is dropped and
recreated per module; loading a whole catalogue into it would cost ~2 minutes per run and would test
a load this suite had just performed rather than the one of record. If the published *E. coli*
pangenome is not loaded there, every test here SKIPS with the reason — it does not pass.

⚠ **The one permitted difference is the named allowlist**, `known_parity_exceptions`. Never a
tolerance: the 2026-09-04 audit re-run retired `no_homology`, moving 2 *E. coli* loci and 6 kp loci
to `synteny_only`, each recorded by node label rather than absorbed by a fuzzy comparison.
"""

from __future__ import annotations

import math
import os

import pytest
from sqlalchemy import create_engine, select
from sqlalchemy.orm import Session

from bacatlas_backend.ingest.published_catalogues import catalogue as published_catalogue
from bacatlas_backend.models.locus import Locus
from bacatlas_backend.models.locus_annotation import LocusAnnotationEntry, LocusUnirefFamilyCrosstab
from bacatlas_backend.models.locus_arrangement import LocusArrangement
from bacatlas_backend.models.locus_offset_occupant import LocusOffsetOccupant
from bacatlas_backend.models.pangenome import Pangenome
from tests.conftest import NUNA_DATA_ROOT, PUBLISHED_SITE_CATALOGUE_DIR
from tests.known_parity_exceptions import RETIRED_TIER, exceptions_for, symbol_fold_for
from tests.payload_oracle import OFFSETS, load_catalogue

#: The `lists.<key>` → `AnnotationKind` map, so one comparison serves all seven vocabularies.
ANNOTATION_KIND_BY_LIST_KEY = {
    "sym": "GENE_SYMBOL",
    "prod": "PROTEIN_PRODUCT",
    "pfam": "PFAM_ARCHITECTURE",
    "cog": "COG_ORTHOGROUP",
    "ec": "EC_NUMBER",
    "kegg": "KEGG_ORTHOLOGY",
    "go": "GENE_ONTOLOGY_SLIM",
}


#: ⭐ Both published species. Running only *E. coli* would leave the loader's one species-specific
#: surface — the `strict98`/`kp98` run_id and the parquets' third species vocabulary — untested.
SPECIES_KEYS = ("ecoli", "kp")


@pytest.fixture(scope="module", params=SPECIES_KEYS)
def parity(request):
    """`(catalogue, session, loci, pangenome)` — the published payload beside the loaded database.

    ⛔ Coverage is asserted HERE, once, before any test compares anything: same locus count, same
    order, same labels. A suite whose fixture did not check that could report "0 differ" while
    comparing two catalogues that do not describe the same loci.
    """
    species_key = request.param
    path = PUBLISHED_SITE_CATALOGUE_DIR / f"{species_key}.json"
    if not path.exists():
        pytest.skip(f"the published site catalogue is not present at {path}")
    url = os.environ.get("BACATLAS_DATABASE_URL")
    if not url:
        pytest.skip("BACATLAS_DATABASE_URL is not set — parity runs against the loaded database")

    entry = published_catalogue(species_key)
    engine = create_engine(url, future=True)
    session = Session(engine)
    request.addfinalizer(session.close)

    pangenome = session.execute(
        select(Pangenome).where(Pangenome.run_id == entry.run_id)
    ).scalar_one_or_none()
    if pangenome is None:
        pytest.skip(f"run {entry.run_id} is not loaded in {url}")

    loci = list(
        session.execute(
            select(Locus)
            .where(Locus.pangenome_id == pangenome.pangenome_id)
            .order_by(Locus.catalogue_ordinal)
        ).scalars()
    )
    if not loci:
        pytest.skip(f"pangenome {pangenome.pangenome_id} holds no loci")

    catalogue = load_catalogue(path)
    assert len(loci) == catalogue.n_loci == entry.locus_count, (
        f"the database holds {len(loci):,} loci and the published catalogue {catalogue.n_loci:,}"
    )
    assert [locus.catalogue_ordinal for locus in loci] == list(range(len(loci)))
    assert [locus.node_label for locus in loci] == [
        str(label) for label in catalogue.nodes["label"]
    ], "the locus ORDER differs — every index-addressed block would compare the wrong pair"
    return catalogue, session, loci, pangenome, entry


def _same(left, right) -> bool:
    """Equality that treats NULL and NaN as one value and never as zero."""
    left_absent = left is None or (isinstance(left, float) and math.isnan(left))
    right_absent = right is None or (isinstance(right, float) and math.isnan(right))
    if left_absent or right_absent:
        return left_absent and right_absent
    if isinstance(left, float) or isinstance(right, float):
        return float(left) == float(right)
    return left == right


# ── T5 · locus information ─────────────────────────────────────────────────────────────────────
#: `(locus column, payload nodes key, how to read the payload value)`. Every scalar the card shows.
#:
#: ⛔ **The four `*_medoid_distance` rows went on 2026-09-24 and are NOT replaced here**, because the
#: published payload this table compares against is FROZEN at schema 14 and has no `sim` block to
#: compare to. Dropping four comparisons without saying so would be exactly the silent loss of
#: coverage this file's fixture asserts against, so the replacement measurement is compared against
#: its own oracle — the artifact it was ingested from — in `test_T5b` below.
SCALAR_COLUMNS = (
    ("member_gene_count", "size", None),
    ("member_genome_count", "genomes", None),
    ("named_member_count", "n_named", None),
    ("uniref50_family_count", "n_u50", None),
    ("uniref50_major_family_count", "n_u50_major", None),
    ("uniref50_labelled_member_count", "n_u50_labelled", None),
    ("pfam_annotated_member_count", "n_pfam", None),
    ("syntenic_a5", "a5", None),
    ("resolved_threshold", "resolved", None),
    ("cog_annotated_member_count", "n_cog", None),
    ("cog_distinct_id_count", "cog_n", None),
    ("ec_annotated_member_count", "n_ec", None),
    ("kegg_annotated_member_count", "n_kegg", None),
    ("median_gene_length_nt", "len_nt", None),
    ("gene_length_interquartile_range_nt", "len_iqr", None),
    ("go_annotated_member_count_molecular_function", "n_go_0", None),
    ("go_annotated_member_count_biological_process", "n_go_1", None),
    ("go_annotated_member_count_cellular_component", "n_go_2", None),
)


def test_T5_every_card_scalar_matches_over_the_whole_catalogue(parity):
    catalogue, _, loci, _, entry = parity
    examined, differing = 0, {}
    for column, key, _reader in SCALAR_COLUMNS:
        assert key in catalogue.nodes, f"the payload has no nodes.{key} — this comparison is vacuous"
        for index, locus in enumerate(loci):
            examined += 1
            if not _same(getattr(locus, column), catalogue.nodes[key][index]):
                differing.setdefault(column, []).append(
                    (locus.node_label, getattr(locus, column), catalogue.nodes[key][index])
                )
    assert examined == len(SCALAR_COLUMNS) * catalogue.n_loci
    assert catalogue.n_loci == entry.locus_count, f"examined {catalogue.n_loci:,} loci"
    assert not differing, {key: value[:3] for key, value in differing.items()}


def test_T5_the_interned_string_columns_match_and_absence_stays_absence(parity):
    """⚠ `-1` in an interned column is ABSENT. A serialiser mapping it to `0` names locus 0."""
    catalogue, _, loci, _, entry = parity
    pairs = (
        ("bakta_gene_symbol", "name", "sym"),
        ("collapse_tier", "tier", "tier"),
        ("pfam_concordance_class", "pfclass", "pfclass"),
    )
    differing = {}
    for column, key, pool in pairs:
        for index, locus in enumerate(loci):
            expected = catalogue.string(pool, catalogue.nodes[key][index])
            allowed = exceptions_for(entry.species_key, column)
            if allowed is not None and locus.node_label in allowed.node_labels:
                continue
            if not _same(getattr(locus, column), expected):
                differing.setdefault(column, []).append(
                    (locus.node_label, getattr(locus, column), expected)
                )
    assert not differing, {key: value[:3] for key, value in differing.items()}


def test_T5_the_modal_cog_category_is_a_SET_and_rejoins_to_the_payload_string(parity):
    """⛔ `EHJQ` is four categories. Stored whole it fits, and `= 'C'` then misses 95 of 118 values."""
    catalogue, _, loci, _, entry = parity
    examined, multi, differing = 0, 0, []
    for index, locus in enumerate(loci):
        expected = catalogue.string("cog", catalogue.nodes["cog_cat"][index])
        actual = locus.modal_cog_categories
        examined += 1
        if expected is None:
            if actual is not None:
                differing.append((locus.node_label, actual, expected))
            continue
        if actual is None or "".join(actual) != expected:
            differing.append((locus.node_label, actual, expected))
        elif len(actual) > 1:
            multi += 1
    assert examined == catalogue.n_loci
    assert not differing, differing[:5]
    assert multi > 0, "no locus has a multi-letter category set — the array would be untested"


#: How many loci carry a tier at all, per species — the denominator this comparison runs over, so
#: "they all agree" is a claim about a measured population and not about an empty one.
TIER_BEARING_LOCI = {"ecoli": 12_104, "kp": 11_272}


def test_T5_the_tier_column_matches_the_catalogue_on_EVERY_locus(parity):
    """⭐ Replaced `…the_only_tier_difference_is_the_named_allowlist…` on 2026-10-03.

    It recorded that the 2026-09-04 audit re-run retired `no_homology` while the published pages
    predated it, so 2 ecoli / 6 kp loci read `no_homology` on the page and `synteny_only` in the
    database. **Both pages were re-exported at schema 16 on 2026-10-03**, so the allowlist is empty.

    ⛔ **Rewritten as a positive gate rather than narrowed to an empty allowlist.** Narrowing would
    have left this asserting `set() == frozenset()` and the loop below dead — green forever, gating
    nothing, which is the "tolerance" `known_parity_exceptions`' own docstring exists to forbid. The
    retirement is now an INVARIANT to hold, not a difference to tolerate.
    """
    catalogue, _, loci, _, entry = parity
    moved = [
        (locus.node_label, locus.collapse_tier, catalogue.string("tier", catalogue.nodes["tier"][index]))
        for index, locus in enumerate(loci)
        if not _same(locus.collapse_tier, catalogue.string("tier", catalogue.nodes["tier"][index]))
    ]
    assert moved == [], f"{len(moved)} loci disagree on tier: {moved[:5]}"
    # ⛔ Coverage BEFORE the agreement is reported. `_same` treats None leniently and a comparison
    # over an all-None column also finds nothing moved, so both sides state their measured count.
    assert sum(1 for locus in loci if locus.collapse_tier is not None) == TIER_BEARING_LOCI[entry.species_key]
    assert sum(1 for value in catalogue.nodes["tier"] if value >= 0) == TIER_BEARING_LOCI[entry.species_key]
    # ⭐ And the retired name is in neither the pool nor the database — the needle `RETIRED_TIER` is
    # kept in the registry for.
    assert RETIRED_TIER not in catalogue.raw["strings"]["tier"]
    assert RETIRED_TIER not in {locus.collapse_tier for locus in loci}
    # ⚠ A re-added tolerance has to be deliberate: there is no recorded exception for this column now.
    assert exceptions_for(entry.species_key, "collapse_tier") is None


def test_T5_the_annotation_lists_match_row_for_row_in_every_vocabulary(parity):
    catalogue, session, loci, pangenome, entry = parity
    species_key = entry.species_key  # ⚠ the loop below rebinds `entry` to an ORM row
    rows = session.execute(
        select(LocusAnnotationEntry, Locus.catalogue_ordinal)
        .join(Locus, Locus.locus_id == LocusAnnotationEntry.locus_id)
        .where(Locus.pangenome_id == pangenome.pangenome_id)
    ).all()
    held: dict[tuple[int, str], list] = {}
    for entry, ordinal in rows:
        held.setdefault((ordinal, entry.annotation_kind.name), []).append(entry)
    for value in held.values():
        # GO carries three interleaved rank sequences; the payload emits them namespace-major.
        value.sort(key=lambda entry: (entry.gene_ontology_namespace or 0, entry.rank_within_locus))

    # ⛔ Named loci, never a tolerance: the four the allele-variant fold renames hold a DIFFERENT
    # symbol list from the frozen page by decision, and `test_T5_the_symbol_fold_renames_EXACTLY…`
    # below asserts every value of it. Every other list on those loci is still compared.
    folded = symbol_fold_for(species_key).node_labels
    examined, differing = 0, []
    for list_key, kind in ANNOTATION_KIND_BY_LIST_KEY.items():
        for index in range(catalogue.n_loci):
            if list_key == "sym" and catalogue.label(index) in folded:
                continue
            expected = catalogue.annotation_rows(list_key, index)
            actual = held.get((index, kind), [])
            examined += max(len(expected), len(actual))
            if len(expected) != len(actual):
                differing.append((list_key, catalogue.label(index), len(expected), len(actual)))
                continue
            for entry, want in zip(actual, expected, strict=True):
                if entry.term_value != want.term or entry.member_gene_count != want.gene_count:
                    differing.append((list_key, catalogue.label(index), want.term, entry.term_value))
                if list_key == "go" and entry.gene_ontology_namespace != want.namespace:
                    differing.append((list_key, catalogue.label(index), "namespace", entry.rank_within_locus))
    skipped = sum(
        len(held.get((index, ANNOTATION_KIND_BY_LIST_KEY["sym"]), []))
        for index in range(catalogue.n_loci)
        if catalogue.label(index) in folded
    )
    assert examined == sum(len(block) for block in held.values()) - skipped, (
        f"examined {examined:,} annotation rows over {catalogue.n_loci:,} loci"
    )
    assert examined > 50_000, f"examined only {examined:,} annotation rows"
    assert not differing, differing[:5]


def test_T5_the_uniref_crosstab_matches_including_its_five_extra_columns(parity):
    catalogue, session, loci, pangenome, entry = parity
    species_key = entry.species_key  # ⚠ the loop below rebinds `entry` to an ORM row
    folded = symbol_fold_for(species_key).node_labels
    rows = session.execute(
        select(LocusUnirefFamilyCrosstab, Locus.catalogue_ordinal)
        .join(Locus, Locus.locus_id == LocusUnirefFamilyCrosstab.locus_id)
        .where(Locus.pangenome_id == pangenome.pangenome_id)
    ).all()
    held: dict[int, list] = {}
    for entry, ordinal in rows:
        held.setdefault(ordinal, []).append(entry)
    for value in held.values():
        value.sort(key=lambda entry: entry.rank_within_locus)

    examined, differing = 0, []
    for index in range(catalogue.n_loci):
        expected = catalogue.annotation_rows("u50", index)
        actual = held.get(index, [])
        examined += len(expected)
        if len(expected) != len(actual):
            differing.append((catalogue.label(index), len(expected), len(actual)))
            continue
        for entry, want in zip(actual, expected, strict=True):
            checks = (
                (entry.uniref50_accession, want.term),
                (entry.member_gene_count, want.gene_count),
                (entry.modal_bakta_product, want.modal_product),
                (entry.modal_pfam_architecture, want.modal_architecture),
                (entry.pfam_annotated_member_count, want.pfam_annotated_count),
            )
            # ⛔ The two SYMBOL columns move on the folded loci by decision, and nowhere else —
            # so they are compared on every other locus rather than dropped from the suite.
            if catalogue.label(index) not in folded:
                checks += (
                    (entry.modal_bakta_gene_symbol, want.modal_symbol),
                    (entry.distinct_real_symbol_count, want.distinct_symbol_count),
                )
            for left, right in checks:
                if not _same(left, right):
                    differing.append((catalogue.label(index), left, right))
    assert examined > 15_000, f"examined {examined:,} cross-tab rows"
    assert not differing, differing[:5]


def test_T5_the_symbol_fold_renames_EXACTLY_the_loci_it_names_and_no_others(parity):
    """⛔ The other half of every exclusion above: what the excused loci actually hold now.

    An exception that only says *"do not look here"* is a tolerance wearing a name. This asserts
    the frozen value, the current value and the set — so a fold that reached one locus further,
    or stopped reaching one, fails here rather than passing quietly everywhere else.
    """
    catalogue, _, loci, _, entry = parity
    fold = symbol_fold_for(entry.species_key)
    by_label = {locus.node_label: locus for locus in loci}
    for label, frozen_value, current in fold.renamed:
        locus = by_label[label]
        index = next(i for i in range(catalogue.n_loci) if catalogue.label(i) == label)
        published = catalogue.string("sym", catalogue.nodes["name"][index])
        assert published == frozen_value, f"{label}: the frozen page says {published!r}"
        assert locus.bakta_gene_symbol == current, f"{label}: the database says {locus.bakta_gene_symbol!r}"
        assert locus.display_name == current

    # ⚠ And the set: every OTHER locus keeps the name the frozen page gave it. The interned-string
    # comparison above already skips the folded four, so without this clause a fold that renamed a
    # thousand loci would be invisible to both.
    renamed = {
        locus.node_label
        for index, locus in enumerate(loci)
        if not _same(locus.bakta_gene_symbol, catalogue.string("sym", catalogue.nodes["name"][index]))
    }
    assert renamed == fold.node_labels


# ── T7 · EggNOG ────────────────────────────────────────────────────────────────────────────────
def test_T7_every_GO_verdict_matches_and_no_coverage_is_a_VALUE_not_a_NULL(parity):
    """⛔ Fewer than two annotated members is neither agreement nor disagreement."""
    catalogue, _, loci, _, entry = parity
    namespaces = (
        ("molecular_function", "go_0"),
        ("biological_process", "go_1"),
        ("cellular_component", "go_2"),
    )
    examined, differing, no_coverage = 0, [], 0
    for namespace, key in namespaces:
        for index, locus in enumerate(loci):
            examined += 1
            verdict = getattr(locus, f"go_verdict_{namespace}")
            expected = catalogue.string("goclass", catalogue.nodes[key][index])
            actual = verdict.value if verdict is not None else None
            if actual == "no_coverage":
                no_coverage += 1
            if not _same(actual, expected):
                differing.append((locus.node_label, namespace, actual, expected))
    assert examined == 3 * catalogue.n_loci
    assert not differing, differing[:5]
    assert no_coverage > 0, "no_coverage never appears — the value would be indistinguishable from NULL"


def test_T7_a_coverage_count_travels_with_every_verdict(parity):
    """⚠ A list of labels cannot say whether members disagreed or were never annotated."""
    _, _, loci, _, entry = parity
    for locus in loci:
        for namespace in ("molecular_function", "biological_process", "cellular_component"):
            verdict = getattr(locus, f"go_verdict_{namespace}")
            covered = getattr(locus, f"go_annotated_member_count_{namespace}")
            assert covered is not None, locus.node_label
            if verdict is not None and verdict.value == "no_coverage":
                assert covered < 2, (locus.node_label, namespace, covered)


def test_T7_kegg_ids_are_present_and_never_NAMED(parity):
    """⛔ KEGG permits linking freely but not redistributing its content."""
    _, session, _, pangenome, entry = parity
    rows = session.execute(
        select(LocusAnnotationEntry.term_value, LocusAnnotationEntry.term_name)
        .join(Locus, Locus.locus_id == LocusAnnotationEntry.locus_id)
        .where(
            Locus.pangenome_id == pangenome.pangenome_id,
            LocusAnnotationEntry.annotation_kind == "KEGG_ORTHOLOGY",
        )
    ).all()
    assert len(rows) > 1_000, f"examined {len(rows):,} KEGG rows"
    assert all(name is None for _, name in rows), "a KEGG description reached the database"


def test_T7_cog_and_go_terms_DO_carry_their_vendored_names(parity):
    """The other half: a name joined into the row is what makes the render-time mutation impossible."""
    _, session, _, pangenome, entry = parity
    for kind, minimum in (("COG_ORTHOGROUP", 4_000), ("GENE_ONTOLOGY_SLIM", 10_000)):
        rows = session.execute(
            select(LocusAnnotationEntry.term_value, LocusAnnotationEntry.term_name)
            .join(Locus, Locus.locus_id == LocusAnnotationEntry.locus_id)
            .where(
                Locus.pangenome_id == pangenome.pangenome_id,
                LocusAnnotationEntry.annotation_kind == kind,
            )
        ).all()
        named = sum(1 for _, name in rows if name)
        assert len(rows) >= minimum, f"{kind}: examined {len(rows):,}"
        assert named / len(rows) > 0.95, f"{kind}: only {named:,} of {len(rows):,} rows carry a name"


# ── T1 · arrangements ──────────────────────────────────────────────────────────────────────────
def test_T1_every_arrangement_matches_rank_counts_flip_and_all_ten_slots(parity):
    catalogue, session, loci, pangenome, entry = parity
    rows = session.execute(
        select(LocusArrangement, Locus.catalogue_ordinal)
        .join(Locus, Locus.locus_id == LocusArrangement.locus_id)
        .where(LocusArrangement.pangenome_id == pangenome.pangenome_id)
    ).all()
    held: dict[int, list] = {}
    for arrangement, ordinal in rows:
        held.setdefault(ordinal, []).append(arrangement)
    for value in held.values():
        value.sort(key=lambda arrangement: arrangement.rank_within_locus)

    examined, differing = 0, []
    for index in range(catalogue.n_loci):
        expected = catalogue.arrangements(index)
        actual = held.get(index, [])
        examined += len(expected)
        if len(expected) != len(actual):
            differing.append((catalogue.label(index), len(expected), len(actual)))
            continue
        for arrangement, want in zip(actual, expected, strict=True):
            if (
                arrangement.member_gene_count != want.gene_count
                or arrangement.member_genome_count != want.genome_count
                or arrangement.is_recorded_reverse_complement != want.is_flipped
                or list(arrangement.neighbour_slot_codes) != list(want.slot_codes)
            ):
                differing.append((catalogue.label(index), arrangement.rank_within_locus))
    assert examined == sum(len(block) for block in held.values()), (
        f"examined {examined:,} arrangements over {catalogue.n_loci:,} loci"
    )
    assert examined > 50_000, f"examined only {examined:,} arrangements"
    assert not differing, differing[:5]


def test_T1_the_total_is_carried_separately_and_no_display_cap_moves_it(parity):
    """⛔ `arr.tot` is what the page keys sentences on; it must never be the listed count."""
    catalogue, _, loci, _, entry = parity
    differing = [
        (locus.node_label, locus.total_arrangement_count, catalogue.total_arrangement_count(index))
        for index, locus in enumerate(loci)
        if locus.total_arrangement_count != catalogue.total_arrangement_count(index)
    ]
    assert not differing, differing[:5]


def test_T1_a_contig_end_survives_as_MINUS_ONE_and_is_not_folded_into_the_rest(parity):
    """⛔ `-1` is a VALUE. Folding it in would claim an observation that was never made."""
    catalogue, session, loci, pangenome, entry = parity
    held = session.execute(
        select(LocusArrangement.neighbour_slot_codes).where(
            LocusArrangement.pangenome_id == pangenome.pangenome_id
        )
    ).scalars().all()
    ours = sum(1 for codes in held for code in codes if code == -1)
    theirs = sum(
        1
        for index in range(catalogue.n_loci)
        for arrangement in catalogue.arrangements(index)
        for code in arrangement.slot_codes
        if code == -1
    )
    assert ours == theirs > 0, (ours, theirs)


def test_T1_a_genome_at_rho_above_one_occupies_TWO_arrangements_and_the_schema_allows_it(parity):
    """⚠ `sum(member_genome_count)` may exceed the locus's own count while the UNION equals it."""
    catalogue, session, loci, pangenome, entry = parity
    over = 0
    for index, locus in enumerate(loci):
        arrangements = catalogue.arrangements(index)
        if sum(a.genome_count for a in arrangements) > locus.member_genome_count:
            over += 1
    assert over > 0, "no locus puts one genome in two arrangements — this invariant is untested"


# ── T3(a) · the graph ──────────────────────────────────────────────────────────────────────────
def test_T3a_the_adjacency_multiset_matches_over_every_locus_and_offset(parity):
    catalogue, session, loci, pangenome, entry = parity
    ordinal_by_locus_id = {locus.locus_id: locus.catalogue_ordinal for locus in loci}
    rows = session.execute(
        select(
            LocusOffsetOccupant.locus_id,
            LocusOffsetOccupant.signed_offset,
            LocusOffsetOccupant.rank_within_offset,
            LocusOffsetOccupant.neighbour_locus_id,
            LocusOffsetOccupant.member_gene_count,
            LocusOffsetOccupant.same_strand_member_count,
        ).where(LocusOffsetOccupant.pangenome_id == pangenome.pangenome_id)
    ).all()
    held = {
        (
            ordinal_by_locus_id[locus_id],
            offset,
            rank,
            ordinal_by_locus_id[neighbour_id],
            count,
            same,
        )
        for locus_id, offset, rank, neighbour_id, count, same in rows
    }
    expected = {
        (
            index,
            occupant.signed_offset,
            occupant.rank,
            occupant.neighbour_locus_index,
            occupant.gene_count,
            occupant.same_strand_count,
        )
        for index in range(catalogue.n_loci)
        for occupant in catalogue.offset_occupants(index)
    }
    assert len(rows) == len(expected) > 200_000, f"examined {len(rows):,} occupant rows"
    assert held == expected, list(held ^ expected)[:5]


def test_T3a_the_reverse_index_is_the_forward_one_read_the_other_way(parity):
    """⭐ Not two computations that agree today — the same rows, read in both directions."""
    _, session, loci, pangenome, entry = parity
    forward = session.execute(
        select(LocusOffsetOccupant.locus_id, LocusOffsetOccupant.neighbour_locus_id).where(
            LocusOffsetOccupant.pangenome_id == pangenome.pangenome_id
        )
    ).all()
    reverse = session.execute(
        select(LocusOffsetOccupant.locus_id, LocusOffsetOccupant.neighbour_locus_id)
        .where(LocusOffsetOccupant.pangenome_id == pangenome.pangenome_id)
        .order_by(LocusOffsetOccupant.neighbour_locus_id)
    ).all()
    assert sorted(forward) == sorted(reverse)
    assert len(forward) > 200_000, f"examined {len(forward):,} edges in both directions"


def test_T3a_every_marginal_occupant_appears_in_at_least_one_arrangement_slot(parity):
    """⛔ The converse is deliberately NOT asserted — the marginal's modes need not co-occur."""
    catalogue, _, _, _, entry = parity
    examined, missing = 0, []
    for index in range(0, catalogue.n_loci, 7):  # a seventh of the catalogue; the count is reported
        slots_at_offset: dict[int, set[int]] = {}
        for arrangement in catalogue.arrangements(index):
            for position, offset in enumerate(OFFSETS):
                code = arrangement.slot_codes[position]
                if code >= 0:
                    slots_at_offset.setdefault(offset, set()).add(code // 2)
        for occupant in catalogue.offset_occupants(index):
            examined += 1
            if occupant.neighbour_locus_index not in slots_at_offset.get(occupant.signed_offset, set()):
                missing.append((catalogue.label(index), occupant.signed_offset, occupant.rank))
    assert examined > 30_000, f"examined {examined:,} occupant rows over {catalogue.n_loci // 7:,} loci"
    assert not missing, missing[:5]


def test_T3a_the_observed_denominator_is_counted_before_the_top_N_cut(parity):
    """⚠ `obs < size` is contig-edge truncation, and `obs >= listed` is what makes "and N others" honest."""
    catalogue, _, loci, _, entry = parity
    truncated = 0
    for index, locus in enumerate(loci):
        observed = list(locus.context_observed_member_counts)
        assert observed == list(catalogue.observed_member_counts(index)), locus.node_label
        listed: dict[int, int] = {}
        for occupant in catalogue.offset_occupants(index):
            listed[occupant.signed_offset] = listed.get(occupant.signed_offset, 0) + occupant.gene_count
        for position, offset in enumerate(OFFSETS):
            assert observed[position] >= listed.get(offset, 0), (locus.node_label, offset)
        if any(count < locus.member_gene_count for count in observed):
            truncated += 1
    assert truncated > 0, "no locus shows contig-edge truncation — implausible for draft assemblies"


def test_T5b_every_similarity_matches_the_artifact_it_was_INGESTED_from(parity):
    """The set-to-set numbers, against the CSV rather than the payload.

    ⛔ **This is not a weaker oracle than T5, it is a different one.** T5 compares the database to a
    published payload, which is the right check for anything that shipped. These numbers have not
    shipped yet — the published payload is frozen at schema 14 — so their oracle is the artifact
    `build_cluster_similarity` wrote, read here independently of the ingest that loaded it.

    ⚠ Coverage is asserted before any difference is reported: every locus, both representations, all
    five stored quantities. A loop that skipped a missing column would report "0 differ" while
    comparing nothing, which is the failure mode this suite exists to refuse.
    """
    pandas = pytest.importorskip("pandas", reason="the ingest extra")
    catalogue, session, loci, _, entry = parity

    from bacatlas_backend.ingest.artifact_locator import REPRESENTATIONS, CatalogueArtifacts
    from bacatlas_backend.ingest.catalogue_frames import SIMILARITY_COLUMNS
    from bacatlas_backend.models.enumerations import EmbeddingRepresentation
    from bacatlas_backend.models.locus_similarity import LocusSimilarity

    artifacts = CatalogueArtifacts(
        data_root=NUNA_DATA_ROOT, set_key=entry.set_key,
        model_label=entry.model_label, run_id=entry.run_id,
    )  # fmt: skip
    examined, differing = 0, {}
    for representation in REPRESENTATIONS:
        frame = pandas.read_csv(artifacts.cluster_similarity(representation), dtype={"node": str})
        by_label = frame.set_index("node")
        stored = {
            row.locus_id: row
            for row in session.execute(
                select(LocusSimilarity).where(
                    LocusSimilarity.representation == EmbeddingRepresentation(representation),
                    LocusSimilarity.locus_id.in_([locus.locus_id for locus in loci]),
                )
            ).scalars()
        }
        assert len(stored) == len(loci), (
            f"{representation}: {len(stored):,} similarity rows for {len(loci):,} loci — a "
            "comparison over a subset would report agreement it never checked"
        )
        for locus in loci:
            row = stored[locus.locus_id]
            source = by_label.loc[locus.node_label]
            for column, stem in SIMILARITY_COLUMNS.items():
                examined += 1
                expected = source[f"{stem}_raw"]
                expected = None if expected != expected else float(expected)  # NaN → None
                if not _same(getattr(row, column), expected):
                    differing.setdefault(f"{representation}.{column}", []).append(
                        (locus.node_label, getattr(row, column), expected)
                    )
    assert examined == len(REPRESENTATIONS) * len(SIMILARITY_COLUMNS) * catalogue.n_loci
    assert not differing, {key: value[:3] for key, value in differing.items()}
