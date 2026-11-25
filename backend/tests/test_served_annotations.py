"""annotation-transfer's served predictions — the table, the ingest that fills it, and the tab that reads it.

Three databases, in three sections, each for what only it can show:

* **the probe schema** (`BACATLAS_PROBE_DATABASE_URL`, rebuilt per module) — the CHECKs on the real
  tables, every refusal of the ingest on small hand-built files whose one defect is the point, the five
  answers of the tab, and a catalogue reload that holds predictions;
* **a loaded database** (`BACATLAS_DATABASE_URL`, both probe catalogues published) + **the real probe
  predictions** (`BACATLAS_SERVED_ANNOTATIONS_DIR`, CSD3 job 40993415's `served_annotations_at-v0.1/`) — the
  real files load, value for value, and reconcile with the producer's own gain table
  (`BACATLAS_TRANSFER_GAINS_TABLE`), inside a transaction rolled back;
* **the same loaded database with the predictions persisted** (`ingest --stage served-annotations`) — the
  Function tab serves each one with its own accuracy and its donor.

⭐ Every fixture value that a mutation could round, swap or transpose is DISTINCT and written out to full
precision (review of 9cfa481: two-decimal accuracies hid a round-up, equal species hid the donor join).
"""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

import pandas as pd
import pytest
from sqlalchemy import create_engine, func, select, text
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session

from bacatlas_backend.ingest.ingest_command_line import main as ingest_main
from bacatlas_backend.ingest.ingest_locus_catalogue import _delete_pangenome_layer
from bacatlas_backend.ingest.ingest_served_annotations import (
    CATALOGUE_RUN_COLUMN,
    DONOR_RUN_COLUMN,
    VERSION_COLUMN,
    ServedAnnotationsRefused,
    ServedDonorsBlockReload,
    clear_served_annotations,
    gene_ontology_slim_classes,
    load_served_annotations,
    newest_pangenome,
    order_lineage,
    record_none_produced,
    refuse_reload_while_donating,
)
from bacatlas_backend.ingest.published_catalogues import catalogue as published_catalogue
from bacatlas_backend.instruments.gene_ontology_fold import GeneOntologyFold
from bacatlas_backend.models.enumerations import ServedAnnotationKind, ServedAnnotationState
from bacatlas_backend.models.genome_collection import GenomeCollection
from bacatlas_backend.models.pangenome import Pangenome
from bacatlas_backend.models.pathogen_species import PathogenSpecies
from bacatlas_backend.models.served_locus_annotation import ServedLocusAnnotation
from bacatlas_backend.services.reference_name_service import clear_gene_ontology_fold_cache
from bacatlas_backend.services.served_annotation_service import load_predictions
from tests.conftest import make_locus

#: A small GO slim with every shape the lineage rules meet. Cellular component (2): plasma membrane ⊂
#: membrane. Molecular function (0): catalytic activity. ⭐ Biological process (1): a FOUR-class lineage —
#: transcription ⊂ RNA metabolic process ⊂ metabolic process and transcription ⊂ biosynthetic process ⊂
#: metabolic process — the real probe's most common GO set (×231), so the order of ancestors is checked.
PLASMA_MEMBRANE, MEMBRANE, CATALYTIC = "GO:0005886", "GO:0016020", "GO:0003824"
TRANSCRIPTION, RNA_METABOLIC, BIOSYNTHETIC, METABOLIC = "GO:0006351", "GO:0016070", "GO:0009058", "GO:0008152"
#: A raw (non-slim) GO term and a namespace root: both must be refused as predicted classes.
RAW_TERM, CC_ROOT = "GO:0000109", "GO:0005575"
FOLD = GeneOntologyFold(
    ancestors={
        PLASMA_MEMBRANE: frozenset({MEMBRANE}),
        TRANSCRIPTION: frozenset({RNA_METABOLIC, BIOSYNTHETIC, METABOLIC}),
        RNA_METABOLIC: frozenset({METABOLIC}),
        BIOSYNTHETIC: frozenset({METABOLIC}),
    },
    namespace_of={
        PLASMA_MEMBRANE: 2, MEMBRANE: 2, CATALYTIC: 0, RAW_TERM: 2, CC_ROOT: 2,
        TRANSCRIPTION: 1, RNA_METABOLIC: 1, BIOSYNTHETIC: 1, METABOLIC: 1,
    },  # fmt: skip
    name_of={PLASMA_MEMBRANE: "plasma membrane", MEMBRANE: "membrane", CATALYTIC: "catalytic activity"},
)
SLIM = frozenset({PLASMA_MEMBRANE, MEMBRANE, CATALYTIC, TRANSCRIPTION, RNA_METABOLIC, BIOSYNTHETIC, METABOLIC})
#: The lineage as it must be stored: the stated class, then each class before its own ancestors.
TRANSCRIPTION_LINEAGE = [TRANSCRIPTION, BIOSYNTHETIC, RNA_METABOLIC, METABOLIC]

SEED_RUN = "probe_run_res0.1_seed0"
OTHER_RUN = "other_run_res0.1_seed0"


@pytest.fixture(autouse=True)
def _no_cached_fold():
    """⛔ The fold is cached per PROCESS; a fold read from one database must not answer for another."""
    clear_gene_ontology_fold_cache()
    yield
    clear_gene_ontology_fold_cache()


# ── the probe schema: hand-built loci, hand-built files ─────────────────────────────────────────
@pytest.fixture()
def loci(session, seeded):
    """Four loci in the seed catalogue and one in a SECOND catalogue of a DIFFERENT species.

    ⚠ The donor "20" is created FIRST, so its locus_id is LOWER than its recipients' — the order in which
    a non-deferred donor key failed a catalogue reload (review of 9cfa481).
    """
    pangenome = seeded["pangenome"]
    species = PathogenSpecies(species_key="probe2", scientific_name="Klebsiella pneumoniae")
    session.add(species)
    session.flush()
    collection = GenomeCollection(
        pathogen_species_id=species.pathogen_species_id, collection_key="probeother",
        collection_label="Kp probe", roster_lineage=seeded["collection"].roster_lineage, genome_count=100,
    )  # fmt: skip
    session.add(collection)
    session.flush()
    other = Pangenome(
        run_id=OTHER_RUN, catalogue_key="probe2-other", model_label="other_x", ingest_generation=1,
        pathogen_species_id=species.pathogen_species_id, genome_collection_id=collection.genome_collection_id,
        genome_count=100, gene_count=1, locus_count=1,
    )  # fmt: skip
    session.add(other)
    session.flush()
    made = {}
    for key, ordinal, label, size, owner in (
        ("20", 20, "20", 40, pangenome),
        ("10", 10, "10", 5, pangenome),
        ("11", 11, "11", 7, pangenome),
        ("12", 12, "12", 2, pangenome),
        ("other:30", 30, "30", 9, other),
    ):
        locus = make_locus(seeded, ordinal=ordinal, label=label, member_gene_count=size,
                           pangenome_id=owner.pangenome_id, pathogen_species_id=owner.pathogen_species_id)  # fmt: skip
        session.add(locus)
        session.flush()
        made[key] = locus
    assert made["20"].locus_id < made["10"].locus_id, "the fixture's donor must sort before its recipients"
    return {"pangenome": pangenome, "other": other, **made}


def _row(**overrides):
    row = {
        "recipient_node": "10",
        "claim": "cog",
        "predicted_classes": "P|E",
        "n_predicted_classes": 2,
        "p_serve": 0.8367346938775511,
        "donor_node": "20",
        "esm_medoid_sim": 0.971234,
        "bac_medoid_sim": 0.814567,
        "rank_esm": 1,
        "rank_bac": 1,
        "n_recipient_genes": 5,
        "n_donor_genes": 40,
        "donor_catalogue": "probe",
        "direction": "probe->probe",
        "both_directions": False,
        "directions_agree": True,
        VERSION_COLUMN: "at-v0.1",
        CATALOGUE_RUN_COLUMN: SEED_RUN,
        DONOR_RUN_COLUMN: SEED_RUN,
    }
    row.update(overrides)
    return row


def _from_other(**overrides):
    """A row whose donor is locus 30 of the OTHER catalogue (another species)."""
    fields = {"donor_node": "30", "n_donor_genes": 9, "donor_catalogue": "other", DONOR_RUN_COLUMN: OTHER_RUN}
    return _row(**{**fields, **overrides})


def _write(tmp_path: Path, rows, *, drop: tuple[str, ...] = (), name: str = "probe_served_annotations.tsv") -> Path:
    path = tmp_path / name
    pd.DataFrame(rows).drop(columns=list(drop)).to_csv(path, sep="\t", index=False)
    return path


def _load(session, loci, path, *, pangenome=None, set_key="probe"):
    return load_served_annotations(
        session, pangenome=pangenome or loci["pangenome"], set_key=set_key, path=path, fold=FOLD, slim_classes=SLIM
    )


GOOD_ROWS = [
    _row(),
    _row(recipient_node="11", n_recipient_genes=7, claim="cc", predicted_classes=f"{MEMBRANE}|{PLASMA_MEMBRANE}",
         p_serve=0.80, esm_medoid_sim=0.912345, bac_medoid_sim=0.705678),
    _row(recipient_node="11", n_recipient_genes=7, claim="mf", predicted_classes=CATALYTIC, n_predicted_classes=1,
         p_serve=1.0, esm_medoid_sim=0.998877, bac_medoid_sim=0.443322),
    _from_other(recipient_node="12", n_recipient_genes=2, claim="bp", n_predicted_classes=4,
                predicted_classes="|".join([METABOLIC, RNA_METABOLIC, TRANSCRIPTION, BIOSYNTHETIC]),
                p_serve=0.9206349206349206, esm_medoid_sim=0.951111, bac_medoid_sim=0.662222),
]  # fmt: skip


def _stored(session):
    return {
        (row.locus_id, row.annotation_kind, row.gene_ontology_namespace): row
        for row in session.execute(select(ServedLocusAnnotation)).scalars()
    }


def test_the_good_file_stores_ONE_row_per_claim_value_for_value(session, loci, tmp_path):
    """⭐ One row per (locus, claim), the set most specific first, every number exactly as the file wrote it."""
    path = _write(tmp_path, GOOD_ROWS)
    report = _load(session, loci, path)
    assert (report.rows_in_file, report.rows_written, report.version) == (4, 4, "at-v0.1")
    stored = _stored(session)
    assert len(stored) == 4, "one row per (locus, claim) — never one per class"

    cog = stored[(loci["10"].locus_id, ServedAnnotationKind.COG_FUNCTIONAL_CATEGORY, None)]
    assert cog.term_values == ["E", "P"], "COG letters are peers, stored alphabetically"
    # ⛔ Full precision, and the two cosines each in their own column — a round-up or a swap dies here.
    assert cog.estimated_accuracy == 0.8367346938775511
    assert (cog.esm_medoid_similarity, cog.bacformer_medoid_similarity) == (0.971234, 0.814567)
    assert cog.donor_locus_id == loci["20"].locus_id

    cc = stored[(loci["11"].locus_id, ServedAnnotationKind.GENE_ONTOLOGY_SLIM, 2)]
    assert cc.term_values == [PLASMA_MEMBRANE, MEMBRANE], "the stated class first, then what it implies"
    assert cc.estimated_accuracy == 0.80, "the cut is inclusive: p ≥ 0.80 is served"
    mf = stored[(loci["11"].locus_id, ServedAnnotationKind.GENE_ONTOLOGY_SLIM, 0)]
    assert mf.term_values == [CATALYTIC] and mf.estimated_accuracy == 1.0
    bp = stored[(loci["12"].locus_id, ServedAnnotationKind.GENE_ONTOLOGY_SLIM, 1)]
    assert bp.term_values == TRANSCRIPTION_LINEAGE, "each class before its own ancestors, ties by accession"
    assert bp.donor_locus_id == loci["other:30"].locus_id, "the cross-catalogue donor resolved through its run"

    pangenome = loci["pangenome"]
    assert pangenome.served_annotation_state is ServedAnnotationState.LOADED
    assert pangenome.annotation_model_version == "at-v0.1"
    assert pangenome.served_annotation_sha256 == hashlib.sha256(path.read_bytes()).hexdigest()
    assert report.by_donor_catalogue == {"probe": 3, "other": 1}

    # ⚠ A second load REPLACES, says how much it replaced, and re-fingerprints.
    again = _load(session, loci, _write(tmp_path, GOOD_ROWS[:1], name="again.tsv"))
    assert (again.rows_written, again.rows_replaced) == (1, 4)
    assert session.execute(select(func.count()).select_from(ServedLocusAnnotation)).scalar_one() == 1
    assert pangenome.served_annotation_sha256 != report.sha256


def test_none_produced_and_cleared_are_STATED_and_drop_the_rows(session, loci, tmp_path):
    """⭐ Two explicit instructions, two different states — and neither leaves a row or a version behind."""
    _load(session, loci, _write(tmp_path, GOOD_ROWS))
    none = record_none_produced(session, pangenome=loci["pangenome"])
    assert none.rows_replaced == 4 and "NONE PRODUCED" in none.render()
    assert loci["pangenome"].served_annotation_state is ServedAnnotationState.NONE_PRODUCED
    assert (loci["pangenome"].annotation_model_version, loci["pangenome"].served_annotation_sha256) == (None, None)
    _load(session, loci, _write(tmp_path, GOOD_ROWS))
    cleared = clear_served_annotations(session, pangenome=loci["pangenome"])
    assert cleared.rows_replaced == 4 and "CLEARED" in cleared.render()
    assert loci["pangenome"].served_annotation_state is ServedAnnotationState.NOT_LOADED
    assert session.execute(select(func.count()).select_from(ServedLocusAnnotation)).scalar_one() == 0


@pytest.mark.parametrize(
    "dropped,match",
    [
        ((VERSION_COLUMN,), "No `annotation_model_version` column"),
        ((CATALOGUE_RUN_COLUMN,), "No `catalogue_run_id` / `donor_run_id`"),
        ((DONOR_RUN_COLUMN,), "No `catalogue_run_id` / `donor_run_id`"),
    ],
)
def test_a_file_without_its_VERSION_or_its_RUNS_is_refused_and_says_why(session, loci, tmp_path, dropped, match):
    with pytest.raises(ServedAnnotationsRefused, match=match):
        _load(session, loci, _write(tmp_path, GOOD_ROWS, drop=dropped))


@pytest.mark.parametrize(
    "versions,match",
    [
        (["v0.1"], "not registered"),  # ⛔ the bare spelling — nuna's clustering ladder has a v0.1 too
        (["at-v9.9"], "not registered"),
        (["at-v0.1", "at-v0.2"], "names 2 versions"),
        (["at-v0.1", ""], "empty annotation_model_version"),
    ],
)
def test_an_unknown_bare_mixed_or_empty_version_is_refused(session, loci, tmp_path, versions, match):
    rows = [_row(recipient_node=label, n_recipient_genes=size, **{VERSION_COLUMN: version})
            for (label, size), version in zip((("10", 5), ("11", 7)), versions, strict=False)]  # fmt: skip
    with pytest.raises(ServedAnnotationsRefused, match=match):
        _load(session, loci, _write(tmp_path, rows))


def test_a_file_from_ANOTHER_RUN_is_refused_by_its_run_id(session, loci, tmp_path):
    """⛔⛔ The identity check. Every label and every size below would pass; only the run id says so."""
    with pytest.raises(ServedAnnotationsRefused, match="its node labels belong to another run"):
        _load(session, loci, _write(tmp_path, [_row(**{CATALOGUE_RUN_COLUMN: OTHER_RUN})]))


def test_a_donor_run_that_is_NOT_LOADED_or_CONTRADICTS_its_catalogue_is_refused(session, loci, tmp_path):
    with pytest.raises(ServedAnnotationsRefused) as refused:
        _load(session, loci, _write(tmp_path, [_row(), _from_other(**{DONOR_RUN_COLUMN: "absent_run_seed0"})]))
    assert "1 of 2 row(s) take their donor from catalogue ['other'] (run absent_run_seed0)" in str(refused.value)
    assert "NOT LOADED" in str(refused.value) and "cannot be dropped" in str(refused.value)
    with pytest.raises(ServedAnnotationsRefused, match="contradicts their donor_run_id"):
        _load(session, loci, _write(tmp_path, [_row(**{DONOR_RUN_COLUMN: OTHER_RUN})]))


def test_a_label_that_resolves_to_a_locus_of_another_SIZE_is_refused(session, loci, tmp_path):
    """The size check, on both ends — a histogram check, kept beside the run id because it is free."""
    with pytest.raises(ServedAnnotationsRefused, match="DIFFERENT size") as refused:
        _load(session, loci, _write(tmp_path, [_row(n_recipient_genes=6), _row(recipient_node="11",
              n_recipient_genes=7, n_donor_genes=41)]))  # fmt: skip
    assert "recipient" in str(refused.value) and "donor" in str(refused.value), "both ends are checked"
    assert session.execute(select(func.count()).select_from(ServedLocusAnnotation)).scalar_one() == 0


def test_every_problem_is_named_at_once(session, loci, tmp_path):
    """One refusal per run means a second gap costs a second run."""
    rows = [
        _row(p_serve=0.79),  # below the cut
        _row(recipient_node="999"),  # resolves to nothing
        _row(recipient_node="11", n_recipient_genes=7, claim="ec"),  # not a v0.1 claim
        _row(recipient_node="12", n_recipient_genes=2, predicted_classes="S", n_predicted_classes=1),
    ]
    with pytest.raises(ServedAnnotationsRefused) as refused:
        _load(session, loci, _write(tmp_path, rows))
    message = str(refused.value)
    for expected in ("outside [0.8, 1]", "not loci of probe-nuna4", "not predicted by at-v0.1", "function unknown"):
        assert expected in message, (expected, message)


def _line_with(message: str, needle: str) -> str:
    lines = [line for line in message.splitlines() if needle in line]
    assert len(lines) == 1, (needle, message)
    return lines[0]


def test_a_GO_class_in_the_WRONG_namespace_is_refused_ON_ITS_OWN_ROW(session, loci, tmp_path):
    """⚠ Each reason is tied to the row that earned it, so a transposed claim → namespace map cannot pass
    by producing the same two messages from the opposite rows (review of 9cfa481)."""
    rows = [
        _row(claim="mf", predicted_classes=PLASMA_MEMBRANE, n_predicted_classes=1),
        _row(recipient_node="11", n_recipient_genes=7, claim="cc", predicted_classes=PLASMA_MEMBRANE,
             n_predicted_classes=1),  # its ancestor `membrane` is missing: not closed
    ]  # fmt: skip
    with pytest.raises(ServedAnnotationsRefused) as refused:
        _load(session, loci, _write(tmp_path, rows))
    assert f"mf:{PLASMA_MEMBRANE}" in _line_with(str(refused.value), "does not place in namespace 0")
    assert f"cc:{PLASMA_MEMBRANE}" in _line_with(str(refused.value), "not one slim class closed upward")


def test_a_GO_root_or_a_RAW_term_is_refused_as_a_predicted_class(session, loci, tmp_path):
    """The producer strips roots and emits slim classes only; either here means a different producer."""
    rows = [
        _row(claim="cc", predicted_classes=CC_ROOT, n_predicted_classes=1),
        _row(recipient_node="11", n_recipient_genes=7, claim="cc", predicted_classes=RAW_TERM, n_predicted_classes=1),
    ]
    with pytest.raises(ServedAnnotationsRefused) as refused:
        _load(session, loci, _write(tmp_path, rows))
    line = _line_with(str(refused.value), "not a slim class")
    assert f"cc:{CC_ROOT}" in line and f"cc:{RAW_TERM}" in line


def test_order_lineage_finds_the_one_stated_class_and_orders_what_it_implies():
    assert order_lineage(frozenset({MEMBRANE, PLASMA_MEMBRANE}), FOLD) == [PLASMA_MEMBRANE, MEMBRANE]
    assert order_lineage(frozenset(TRANSCRIPTION_LINEAGE), FOLD) == TRANSCRIPTION_LINEAGE
    assert order_lineage(frozenset({MEMBRANE}), FOLD) == [MEMBRANE]
    assert order_lineage(frozenset({PLASMA_MEMBRANE}), FOLD) is None, "an open lineage"
    assert order_lineage(frozenset({MEMBRANE, CATALYTIC}), FOLD) is None, "two stated classes"
    assert order_lineage(frozenset({TRANSCRIPTION, METABOLIC}), FOLD) is None, "ancestors missing between"


def test_the_database_itself_refuses_a_bare_version_a_NaN_and_a_kind_namespace_mismatch(session, loci):
    """⚠ Postgres orders NaN ABOVE every number, so the range CHECKs refuse it — asserted, not assumed."""
    base = dict(
        locus_id=loci["10"].locus_id, annotation_kind=ServedAnnotationKind.COG_FUNCTIONAL_CATEGORY,
        gene_ontology_namespace=None, term_values=["E"], estimated_accuracy=0.9,
        donor_locus_id=loci["20"].locus_id, esm_medoid_similarity=0.9, bacformer_medoid_similarity=0.9,
        annotation_model_version="at-v0.1",
    )  # fmt: skip
    for bad in (
        {"annotation_model_version": "v0.1"},
        {"estimated_accuracy": float("nan")},
        {"esm_medoid_similarity": float("nan")},
        {"estimated_accuracy": 1.2},
        {"gene_ontology_namespace": 2},
        {"term_values": []},
    ):
        with pytest.raises(IntegrityError), session.begin_nested():
            session.add(ServedLocusAnnotation(**{**base, **bad}))
            session.flush()


@pytest.mark.parametrize(
    "state,version,sha",
    [
        (ServedAnnotationState.LOADED, "v0.1", "a" * 64),  # bare version
        (ServedAnnotationState.LOADED, None, "a" * 64),  # loaded, no version
        (ServedAnnotationState.LOADED, "at-v0.1", None),  # loaded, no fingerprint
        (ServedAnnotationState.NOT_LOADED, "at-v0.1", None),  # a version outliving its rows
        (ServedAnnotationState.NONE_PRODUCED, None, "a" * 64),
    ],
)
def test_the_pangenome_cannot_hold_a_state_that_contradicts_its_version_or_file(session, loci, state, version, sha):
    with pytest.raises(IntegrityError), session.begin_nested():
        pangenome = loci["pangenome"]
        pangenome.served_annotation_state, pangenome.annotation_model_version = state, version
        pangenome.served_annotation_sha256 = sha
        session.flush()


def test_the_tab_gives_FIVE_different_answers(session, loci, tmp_path):
    """⛔ Calls; declined for this locus; none produced; not loaded; a version this build cannot show."""
    _load(session, loci, _write(tmp_path, GOOD_ROWS))
    served = load_predictions(session, pangenome=loci["pangenome"], locus=loci["11"])
    assert served["state"] == "loaded"
    assert served["annotation_model_label"] == "annotation-transfer v0.1" and served["serve_cut"] == 0.80
    assert [(call["annotation_kind"], call["gene_ontology_namespace"]) for call in served["calls"]] == [
        ("gene_ontology_slim", "molecular_function"),
        ("gene_ontology_slim", "cellular_component"),
    ]
    cc = served["calls"][1]
    assert (cc["esm_medoid_similarity"], cc["bacformer_medoid_similarity"]) == (0.912345, 0.705678)
    assert cc["terms"] == [{"value": PLASMA_MEMBRANE, "name": None}, {"value": MEMBRANE, "name": None}]
    assert cc["donor"] == {
        "node_label": "20", "display_name": "wzi", "gene_count": 40, "catalogue_key": "probe-nuna4",
        "species_name": "Escherichia coli", "same_catalogue": True,
    }  # fmt: skip
    # ⭐ A donor in ANOTHER catalogue of ANOTHER species: its own key, its own species, never this one's.
    across = load_predictions(session, pangenome=loci["pangenome"], locus=loci["12"])["calls"]
    assert [call["donor"] for call in across] == [{
        "node_label": "30", "display_name": "wzi", "gene_count": 9, "catalogue_key": "probe2-other",
        "species_name": "Klebsiella pneumoniae", "same_catalogue": False,
    }]  # fmt: skip
    assert [term["value"] for term in across[0]["terms"]] == TRANSCRIPTION_LINEAGE

    declined = load_predictions(session, pangenome=loci["pangenome"], locus=loci["20"])
    assert (declined["state"], declined["calls"], declined["annotation_model_version"]) == ("loaded", [], "at-v0.1")

    record_none_produced(session, pangenome=loci["pangenome"])
    none = load_predictions(session, pangenome=loci["pangenome"], locus=loci["11"])
    assert (none["state"], none["calls"], none["annotation_model_version"]) == ("none_produced", [], None)

    clear_served_annotations(session, pangenome=loci["pangenome"])
    unloaded = load_predictions(session, pangenome=loci["pangenome"], locus=loci["11"])
    assert (unloaded["state"], unloaded["calls"]) == ("not_loaded", [])

    # ⛔ Fail closed: rows from a version this build does not register are not shown at all.
    _load(session, loci, _write(tmp_path, GOOD_ROWS))
    loci["pangenome"].annotation_model_version = "at-v9.9"
    session.flush()
    unknown = load_predictions(session, pangenome=loci["pangenome"], locus=loci["11"])
    assert (unknown["state"], unknown["calls"], unknown["serve_cut"]) == ("unknown_model", [], None)
    assert unknown["annotation_model_label"] == "at-v9.9"


def test_the_fingerprint_moves_with_the_LOADED_FILE_not_only_the_version(session, loci, tmp_path):
    """The ETag's data part: a corrected file under the same version is a different fingerprint."""
    pangenome = loci["pangenome"]
    assert pangenome.served_annotation_fingerprint() == "not_loaded"
    _load(session, loci, _write(tmp_path, GOOD_ROWS))
    first = pangenome.served_annotation_fingerprint()
    _load(session, loci, _write(tmp_path, GOOD_ROWS[:2], name="corrected.tsv"))
    assert first.startswith("at-v0.1-") and pangenome.served_annotation_fingerprint() != first
    record_none_produced(session, pangenome=pangenome)
    assert pangenome.served_annotation_fingerprint() == "none_produced"


# ── a catalogue reload while predictions are loaded ─────────────────────────────────────────────
def test_a_reload_with_a_SAME_catalogue_donor_ordered_FIRST_passes_and_resets_the_state(session, loci, tmp_path):
    """⛔⛔ The deferred donor key (review of 9cfa481). Donor 20's locus_id is LOWER than its recipients',
    so a per-row (non-deferred) check fails the moment it is deleted; deferred, the recipients' cascade has
    removed every row by the time the check runs. `SET CONSTRAINTS ALL IMMEDIATE` runs it without a commit."""
    _load(session, loci, _write(tmp_path, GOOD_ROWS[:3]))  # all three take donor 20, in this catalogue
    _delete_pangenome_layer(session, loci["pangenome"].pangenome_id)
    session.execute(text("SET CONSTRAINTS ALL IMMEDIATE"))
    assert session.execute(select(func.count()).select_from(ServedLocusAnnotation)).scalar_one() == 0
    session.refresh(loci["pangenome"])
    assert loci["pangenome"].served_annotation_state is ServedAnnotationState.NOT_LOADED, (
        "⛔ a reload leaves the catalogue NOT LOADED — never naming a version over no rows"
    )
    assert loci["pangenome"].annotation_model_version is None


def test_a_reload_that_would_delete_ANOTHER_catalogues_donors_is_refused_FIRST(session, loci, tmp_path):
    """The other catalogue's rows name a donor here: refused before any frame is built, naming who and how."""
    other_file = _write(tmp_path, [_row(recipient_node="30", n_recipient_genes=9, donor_node="10", n_donor_genes=5,
                                        **{CATALOGUE_RUN_COLUMN: OTHER_RUN, DONOR_RUN_COLUMN: SEED_RUN})],
                        name="other.tsv")  # fmt: skip
    _load(session, loci, other_file, pangenome=loci["other"], set_key="other")
    with pytest.raises(ServedDonorsBlockReload) as refused:
        refuse_reload_while_donating(session, loci["pangenome"].pangenome_id)
    message = str(refused.value)
    assert "probe2-other: 1 prediction(s)" in message
    assert f"--clear-served-annotations --set-key other --model-label other_x --run-id {OTHER_RUN}" in message
    refuse_reload_while_donating(session, loci["other"].pangenome_id)  # nothing donates TO the other one
    # ⚠ And the guard is load-bearing: without it, the deferred key fails that reload at commit.
    with pytest.raises(IntegrityError), session.begin_nested():
        _delete_pangenome_layer(session, loci["pangenome"].pangenome_id)
        session.execute(text("SET CONSTRAINTS ALL IMMEDIATE"))
    clear_served_annotations(session, pangenome=loci["other"])
    refuse_reload_while_donating(session, loci["pangenome"].pangenome_id)  # cleared → allowed


# ── the CLI: an instruction, always, and a refusal that writes nothing ─────────────────────────
def _served_state(engine, run_id):
    with Session(engine) as fresh:
        pangenome = newest_pangenome(fresh, run_id)
        rows = fresh.execute(select(func.count()).select_from(ServedLocusAnnotation)).scalar_one()
        return pangenome.served_annotation_state, pangenome.annotation_model_version, rows


@pytest.mark.parametrize(
    "extra",
    [
        [],  # ⛔ no instruction: never inferred from a path nobody gave
        ["--served-annotations", "/nonexistent/served.tsv"],  # named and absent is NOT "none produced"
        ["--served-annotations", "/nonexistent/served.tsv", "--no-served-annotations"],  # two instructions
    ],
)
def test_the_stage_refuses_a_misconfiguration_WITHOUT_WRITING(engine, seeded, extra):
    url = engine.url.render_as_string(hide_password=False)
    before = _served_state(engine, SEED_RUN)
    argv = ["--database-url", url, "--stage", "served-annotations", "--set-key", "probe",
            "--model-label", "probe_x", "--run-id", SEED_RUN, *extra]  # fmt: skip
    try:
        code = ingest_main(argv)
    except SystemExit as exit_:  # argparse refuses two mutually exclusive flags itself
        code = exit_.code
    assert code == 2
    assert _served_state(engine, SEED_RUN) == before, "a refused instruction wrote something"


def test_the_stage_records_NONE_PRODUCED_only_when_told(engine, seeded):
    url = engine.url.render_as_string(hide_password=False)
    argv = ["--database-url", url, "--stage", "served-annotations", "--set-key", "probe",
            "--model-label", "probe_x", "--run-id", SEED_RUN]  # fmt: skip
    assert ingest_main([*argv, "--no-served-annotations"]) == 0
    assert _served_state(engine, SEED_RUN)[0] is ServedAnnotationState.NONE_PRODUCED
    assert ingest_main([*argv, "--clear-served-annotations"]) == 0
    assert _served_state(engine, SEED_RUN)[0] is ServedAnnotationState.NOT_LOADED


# ── a loaded database + the REAL probe predictions (CSD3 job 40993415) ──────────────────────────
SERVED_DIR = Path(os.environ.get("BACATLAS_SERVED_ANNOTATIONS_DIR", "/nonexistent")).expanduser()
GAINS_TABLE = Path(os.environ.get("BACATLAS_TRANSFER_GAINS_TABLE", "/nonexistent")).expanduser()


@pytest.fixture(scope="module")
def loaded_url():
    url = os.environ.get("BACATLAS_DATABASE_URL")
    if not url:
        pytest.skip("BACATLAS_DATABASE_URL is not set — the real predictions need both catalogues loaded")
    with Session(create_engine(url, future=True)) as session:
        published = session.execute(
            select(func.count()).select_from(PathogenSpecies).where(PathogenSpecies.default_pangenome_id.is_not(None))
        ).scalar_one()
    if published < 2:
        pytest.skip(f"only {published} species are published in {url}")
    return url


@pytest.fixture(scope="module")
def real_files():
    if not (SERVED_DIR / "ecoli_served_annotations.tsv").is_file():
        pytest.skip(f"the probe's served predictions are not at {SERVED_DIR} (BACATLAS_SERVED_ANNOTATIONS_DIR)")
    return SERVED_DIR


def _real_load(session, catalogue: str, path: Path):
    from bacatlas_backend.services.reference_name_service import gene_ontology_fold

    pangenome = newest_pangenome(session, published_catalogue(catalogue).run_id)
    return pangenome, load_served_annotations(
        session, pangenome=pangenome, set_key=catalogue, path=path,
        fold=gene_ontology_fold(session), slim_classes=gene_ontology_slim_classes(session),
    )  # fmt: skip


def test_the_REAL_file_without_its_version_and_run_columns_is_refused(loaded_url, real_files, tmp_path):
    """⛔ Job 39455187's format — the same 16 columns, no version, no runs — must not load."""
    old_format = tmp_path / "kp_old_format.tsv"
    pd.read_csv(real_files / "kp_served_annotations.tsv", sep="\t", dtype=str, keep_default_na=False).drop(
        columns=[VERSION_COLUMN, CATALOGUE_RUN_COLUMN, DONOR_RUN_COLUMN]
    ).to_csv(old_format, sep="\t", index=False)
    with Session(create_engine(loaded_url, future=True)) as session:
        with pytest.raises(ServedAnnotationsRefused, match="No `annotation_model_version` column"):
            _real_load(session, "kp", old_format)
        session.rollback()


def test_the_REAL_kp_file_is_refused_on_the_ecoli_catalogue_by_its_run_id(loaded_url, real_files):
    with Session(create_engine(loaded_url, future=True)) as session:
        with pytest.raises(ServedAnnotationsRefused, match="its node labels belong to another run"):
            _real_load(session, "ecoli", real_files / "kp_served_annotations.tsv")
        session.rollback()


CLAIM_OF = {
    ("COG_FUNCTIONAL_CATEGORY", None): "cog",
    ("GENE_ONTOLOGY_SLIM", 0): "mf",
    ("GENE_ONTOLOGY_SLIM", 1): "bp",
    ("GENE_ONTOLOGY_SLIM", 2): "cc",
}


@pytest.mark.parametrize("catalogue", ["ecoli", "kp"])
def test_the_REAL_probe_predictions_load_VALUE_FOR_VALUE_and_reconcile_with_the_gain_table(
    loaded_url, real_files, catalogue
):
    """⭐ 13,534 + 25,475 = 39,009 predictions; every stored accuracy and both cosines EQUAL to the file's, row
    by row; and each claim's count equal to `transfer_gains_by_catalogue.tsv`, which the producer counted
    separately. Inside one transaction, rolled back: the persisted load is the CLI's job."""
    if not GAINS_TABLE.is_file():
        pytest.skip(f"the gain table is not at {GAINS_TABLE} (BACATLAS_TRANSFER_GAINS_TABLE)")
    path = real_files / f"{catalogue}_served_annotations.tsv"
    gains = pd.read_csv(GAINS_TABLE, sep="\t")
    gains = gains.loc[gains["recipient"] == catalogue].set_index("claim")
    expected = pd.read_csv(path, sep="\t", dtype={"recipient_node": str, "donor_node": str})
    with Session(create_engine(loaded_url, future=True)) as session:
        pangenome, report = _real_load(session, catalogue, path)
        assert report.rows_written == report.rows_in_file == int(gains["gained_nodes"].sum())
        assert report.rows_written == {"ecoli": 13_534, "kp": 25_475}[catalogue]
        assert report.by_claim == {claim: int(n) for claim, n in gains["gained_nodes"].items()}
        other = "kp" if catalogue == "ecoli" else "ecoli"
        assert report.by_donor_catalogue == {
            catalogue: int(gains["from_interpolation"].sum()),
            other: int(gains["from_extrapolation"].sum()),
        }
        stored = pd.DataFrame(
            session.execute(
                text(
                    "select l.node_label as recipient_node, s.annotation_kind::text as kind, "
                    "s.gene_ontology_namespace as ns, s.estimated_accuracy, s.esm_medoid_similarity, "
                    "s.bacformer_medoid_similarity, d.node_label as stored_donor "
                    "from served_locus_annotation s join locus l on l.locus_id = s.locus_id "
                    "join locus d on d.locus_id = s.donor_locus_id where l.pangenome_id = :p"
                ),
                {"p": pangenome.pangenome_id},
            ).all()
        )
        session.rollback()
    stored["claim"] = [
        CLAIM_OF[(kind, None if pd.isna(ns) else int(ns))] for kind, ns in zip(stored.kind, stored.ns, strict=True)
    ]
    joined = expected.merge(stored, on=["recipient_node", "claim"], how="outer", indicator=True, validate="one_to_one")
    # ⛔ Coverage BEFORE the comparison: every file row has exactly one stored row, and nothing else is stored.
    assert (joined["_merge"] == "both").all() and len(joined) == len(expected)
    assert (joined["estimated_accuracy"] == joined["p_serve"]).all(), "an accuracy was changed on the way in"
    assert (joined["esm_medoid_similarity"] == joined["esm_medoid_sim"]).all()
    assert (joined["bacformer_medoid_similarity"] == joined["bac_medoid_sim"]).all()
    assert (joined["stored_donor"] == joined["donor_node"]).all()


# ── the loaded database WITH the predictions persisted: the tab itself ──────────────────────────
@pytest.fixture(scope="module")
def client(loaded_url):
    from bacatlas_backend.application_factory import create_application
    from bacatlas_backend.configuration import Configuration

    with Session(create_engine(loaded_url, future=True)) as session:
        carrying = session.execute(
            select(func.count()).select_from(Pangenome).where(
                Pangenome.served_annotation_state == ServedAnnotationState.LOADED
            )
        ).scalar_one()
    if not carrying:
        pytest.skip(f"no catalogue in {loaded_url} carries predictions; run `ingest --stage served-annotations`")
    return create_application(Configuration(database_url=loaded_url)).test_client()


def _a_locus_with_cross_catalogue_go_and_cog(loaded_url, catalogue_key):
    with Session(create_engine(loaded_url, future=True)) as session:
        return session.execute(
            text(
                "select l.node_label from locus l join pangenome p on p.pangenome_id = l.pangenome_id "
                "join served_locus_annotation s on s.locus_id = l.locus_id "
                "join locus d on d.locus_id = s.donor_locus_id "
                "where p.catalogue_key = :key group by l.locus_id, l.node_label "
                "having count(*) >= 2 and bool_or(s.annotation_kind = 'COG_FUNCTIONAL_CATEGORY') "
                "and bool_or(d.pangenome_id <> l.pangenome_id) "
                "and bool_or(cardinality(s.term_values) >= 2 and s.annotation_kind = 'GENE_ONTOLOGY_SLIM') "
                "order by l.node_label limit 1"
            ),
            {"key": catalogue_key},
        ).scalar_one()


def test_the_tab_serves_each_prediction_with_its_OWN_accuracy_and_a_linkable_donor(client, loaded_url):
    with Session(create_engine(loaded_url, future=True)) as session:
        kp = newest_pangenome(session, published_catalogue("kp").run_id)
        ecoli_key = newest_pangenome(session, published_catalogue("ecoli").run_id).catalogue_key
        kp_key, fingerprint = kp.catalogue_key, kp.served_annotation_fingerprint()
    label = _a_locus_with_cross_catalogue_go_and_cog(loaded_url, kp_key)
    response = client.get(f"/api/v1/catalogues/{kp_key}/loci/{label}/function")
    predictions = response.get_json()["predictions"]
    assert (predictions["state"], predictions["annotation_model_version"]) == ("loaded", "at-v0.1")
    assert predictions["annotation_model_label"] == "annotation-transfer v0.1"
    assert response.headers["ETag"].endswith(f'-{fingerprint}"'), "the tag tracks the loaded FILE"
    # ⚠ The stage writes nothing the shell serves, so the shell keeps the pangenome's own tag.
    shell = client.get(f"/api/v1/catalogues/{kp_key}")
    assert "served_annotations" not in shell.get_json()["pangenome"]["omitted_sections"]

    with Session(create_engine(loaded_url, future=True)) as session:
        stored = session.execute(
            text(
                "select s.term_values, s.estimated_accuracy, s.esm_medoid_similarity, s.bacformer_medoid_similarity, "
                "d.node_label as donor from served_locus_annotation s join locus l on l.locus_id = s.locus_id "
                "join pangenome p on p.pangenome_id = l.pangenome_id join locus d on d.locus_id = s.donor_locus_id "
                "where p.catalogue_key = :key and l.node_label = :label"
            ),
            {"key": kp_key, "label": label},
        ).all()
    assert len(predictions["calls"]) == len(stored) >= 2
    by_values = {tuple(row.term_values): row for row in stored}
    for call in predictions["calls"]:
        row = by_values[tuple(term["value"] for term in call["terms"])]
        assert call["estimated_accuracy"] == row.estimated_accuracy, "the row's own number, never a band"
        assert (call["esm_medoid_similarity"], call["bacformer_medoid_similarity"]) == (
            row.esm_medoid_similarity, row.bacformer_medoid_similarity,
        )  # fmt: skip
        assert call["donor"]["node_label"] == row.donor
        if call["annotation_kind"] == "gene_ontology_slim":
            assert all(term["name"] for term in call["terms"]), "GO classes are named from the reference"
        else:
            assert all(term["name"] is None for term in call["terms"]), "COG letters are named on the page"
    across = [call["donor"] for call in predictions["calls"] if not call["donor"]["same_catalogue"]]
    assert across and all(donor["catalogue_key"] == ecoli_key for donor in across)
    assert all(donor["species_name"] == "Escherichia coli" for donor in across)


def test_the_ladder_is_GONE_from_the_function_response(client):
    payload = client.get("/api/v1/species/kp/loci/722/function").get_json()
    assert "calibration" not in payload
    assert set(payload["inference"]) == {"vocabularies"}
    assert all(set(entry) == {"annotation_kind", "own"} for entry in payload["inference"]["vocabularies"])
