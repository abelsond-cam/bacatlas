"""Two catalogues of one species, held at once and addressed apart.

⭐ **What this protects.** A species may now hold several catalogues — nuna4 beside nuna5, or a
sensitive model beside a less-sensitive one — and the *display* picks which to show (David,
2026-09-25). The storage layer always allowed it; what is new is that a reader can address one.

The failure modes these cover, each of which renders perfectly while being wrong:

**A key that steals another key's payloads.** nuna's payload lookup globs `locus_browser_<key>_*`
and takes the newest match, so `ecoli` would also match `locus_browser_ecoli_nuna5_*.json` and, being
the later export, win — the live page then serves a different model's catalogue in silence. A hyphen
cannot occur inside a species token, so `ecoli-nuna5` is safe where `ecoli_nuna5` is not. Same shape
as `-excl` being a prefix of `-exclLOGP`.

**A resolver that follows the default when it was asked to pin.** A link carrying `ecoli-nuna5` must
keep meaning that clustering after the default moves; that is the entire reason the model is in the
URL. If `resolve_catalogue` ever fell back to the species' default, every shared link would silently
start showing whatever is current, and two people discussing one locus would be looking at different
clusterings while seeing the same URL.

**A second catalogue reaching into the first's rows.** Node labels and ordinals are model-private, so
loading nuna5 must leave nuna4's loci exactly as they were.
"""

from __future__ import annotations

import pytest

from bacatlas_backend.ingest.ingest_pangenome_run import catalogue_key_for
from bacatlas_backend.models.enumerations import ExclusivityForm, ExclusivityFormSource
from bacatlas_backend.models.pangenome import Pangenome
from bacatlas_backend.services.species_catalogue_service import (
    CatalogueNotFound,
    list_catalogues,
    resolve_catalogue,
    resolve_published_pangenome,
)

# ------------------------------------------------------------------ the key, and the prefix trap


def test_the_key_joins_species_and_model_with_a_HYPHEN():
    assert catalogue_key_for("ecoli", "nuna5") == "ecoli-nuna5"


def test_an_UNDERSCORE_key_is_refused_because_the_species_key_would_match_it():
    """⛔ `ecoli` globs `locus_browser_ecoli_*`, which matches `locus_browser_ecoli_nuna5_….json`.

    The newest export wins that glob, so publishing a second model would repoint the FIRST model's
    page with nothing in any output saying so. The hyphen is what makes the keys prefix-free.
    """
    with pytest.raises(Exception, match="underscore"):
        catalogue_key_for("ecoli", "nuna5", "ecoli_nuna5")


def test_a_deliberate_key_that_is_not_the_default_is_allowed():
    """The point is a capability, not a naming scheme — `ecoli-sensitive` is a legitimate key."""
    assert catalogue_key_for("ecoli", "nuna5", "ecoli-sensitive") == "ecoli-sensitive"


# --------------------------------------------------------------------------- two, side by side


def _second_catalogue(session, seeded, *, key="probe-nuna5", run_id="probe_run_nuna5"):
    """A second pangenome for the SAME species, sharing its genomes and collection."""
    pangenome = Pangenome(
        run_id=run_id,
        catalogue_key=key,
        model_label="nuna5_g2_0.98_3b0.5rhoPAIRMAX_3c0.2rhoCEILexclLOGP_step4g0.5rhoCEIL",
        ingest_generation=1,
        pathogen_species_id=seeded["species"].pathogen_species_id,
        genome_collection_id=seeded["collection"].genome_collection_id,
        nuna_model_id=seeded["model"].nuna_model_id,
        exclusivity_form=ExclusivityForm.EXCLUSION,
        exclusivity_form_source=ExclusivityFormSource.RUN_ID_TOKEN,
        run_id_exclusivity_token="-exclLOGP",
        genome_count=100,
        gene_count=489146,
        locus_count=2,
        is_published=True,
    )
    session.add(pangenome)
    session.flush()
    return pangenome


def test_two_catalogues_of_one_species_coexist(session, seeded):
    """⛔ There is no unique constraint on `pangenome(pathogen_species_id)`, and there must not be."""
    second = _second_catalogue(session, seeded)
    assert second.pathogen_species_id == seeded["pangenome"].pathogen_species_id
    assert second.pangenome_id != seeded["pangenome"].pangenome_id


def test_a_catalogue_key_resolves_to_ITS_OWN_pangenome_not_the_species_default(session, seeded):
    """The pin. `probe`'s default is the seed; asking for `probe-nuna5` must not return it."""
    second = _second_catalogue(session, seeded)
    assert resolve_catalogue(session, second.catalogue_key).pangenome_id == second.pangenome_id
    assert resolve_catalogue(session, "probe-nuna4").pangenome_id == seeded["pangenome"].pangenome_id


def test_the_bare_species_follows_the_DEFAULT_and_the_two_answers_differ(session, seeded):
    """⚠ The whole distinction, in one assertion: a species key follows, a catalogue key pins."""
    seeded["species"].default_pangenome_id = seeded["pangenome"].pangenome_id
    second = _second_catalogue(session, seeded)
    session.flush()

    followed = resolve_published_pangenome(session, "probe")
    pinned = resolve_catalogue(session, second.catalogue_key)
    assert followed.pangenome_id == seeded["pangenome"].pangenome_id
    assert pinned.pangenome_id == second.pangenome_id
    assert followed.pangenome_id != pinned.pangenome_id, "a pin that follows the default is not a pin"


def test_an_unknown_key_is_reported_rather_than_falling_back(session, seeded):
    """⛔ A miss must NOT degrade to the species' default — that is a silent wrong-model answer."""
    with pytest.raises(CatalogueNotFound):
        resolve_catalogue(session, "probe-nosuchmodel")


def test_the_picker_lists_by_is_published_not_by_the_default_pointer(session, seeded):
    """`is_published` is visibility; `default_pangenome_id` is one per species.

    A catalogue that is loaded but not offered is reachable by key and absent from the menu — which
    is how a model is staged before anyone is shown it.
    """
    seeded["pangenome"].is_published = True
    hidden = _second_catalogue(session, seeded, key="probe-hidden", run_id="probe_run_hidden")
    hidden.is_published = False
    session.flush()

    keys = [pangenome.catalogue_key for pangenome, _species, _model in list_catalogues(session)]
    assert "probe-nuna4" in keys
    assert "probe-hidden" not in keys, "unpublished is staged, not shown"
    assert resolve_catalogue(session, "probe-hidden").pangenome_id == hidden.pangenome_id, (
        "and it is still addressable by key"
    )


# ----------------------------------------------- the one row two catalogues share, and its guard


def _roster_of(monkeypatch, tmp_path, samples):
    """Point the roster ingest at a controlled genome list, without needing real artifacts."""
    from bacatlas_backend.ingest import ingest_genome_collection_roster as roster_module

    universe = tmp_path / "gene_universe.parquet"
    universe.write_text("not read — the vocabulary is stubbed, only the sha256 touches this file")
    monkeypatch.setattr(roster_module, "read_genome_vocabulary", lambda _artifacts: (samples, 1))

    class _Artifacts:
        gene_universe = universe
        set_key = "probe_ecoli"

    return roster_module, _Artifacts()


def test_a_REORDERED_roster_is_refused_when_another_catalogue_depends_on_the_order(
    session, seeded, monkeypatch, tmp_path
):
    """⛔⛔ The single most dangerous thing about holding two catalogues of one species.

    `genome_collection` is unique on `(pathogen_species_id, collection_key)`, so both models share
    ONE collection row — and the membership is deleted and rewritten on every pangenome ingest, with
    `collection_genome_ordinal = enumerate(samples)`. That ordinal is what `meta.genomes` is ordered
    by, and `arr.gid` indexes into `meta.genomes`. So a second model whose roster lists the same
    genomes in a different ORDER renames every genome on every arrangement of the model already
    loaded — and every name it then shows is a real genome, so nothing downstream can detect it.

    The pre-existing check refuses HOLES. This refuses REORDERING, which only became reachable when
    a species was allowed to hold more than one catalogue.
    """
    from bacatlas_backend.models.enumerations import SampleIdentifierKind
    from bacatlas_backend.models.genome import Genome

    second = Genome(
        pathogen_species_id=seeded["species"].pathogen_species_id,
        sample_id="SAMEA0000002",
        sample_id_kind=SampleIdentifierKind.BIOSAMPLE,
        strand_is_observed=True,
    )
    session.add(second)
    session.flush()

    ordered = ["SAMEA0000001", "SAMEA0000002"]
    roster_module, artifacts = _roster_of(monkeypatch, tmp_path, ordered)
    roster_module.ingest_genome_collection(
        session, artifacts, pathogen_species_id=seeded["species"].pathogen_species_id
    )
    session.flush()

    # The seed pangenome already points at this collection, so the order is load-bearing.
    _roster_of(monkeypatch, tmp_path, list(reversed(ordered)))
    with pytest.raises(roster_module.RosterError) as raised:
        roster_module.ingest_genome_collection(
            session, artifacts, pathogen_species_id=seeded["species"].pathogen_species_id
        )

    # ⚠ Assert WHY it refused, not merely that it did. A test that accepts any RosterError would
    # still pass if the guard never ran and the roster failed for an unrelated reason — which is
    # how a guard comes to be inert while its own test stays green.
    message = str(raised.value)
    assert "DIFFERENT roster" in message
    assert "probe-nuna4" in message, "it must name the catalogue that depends on the order"
    assert "ordinal 0" in message, "and where the two rosters first disagree"
    assert "SAMEA0000001" in message and "SAMEA0000002" in message


def test_an_IDENTICAL_roster_is_still_idempotent(session, seeded, monkeypatch, tmp_path):
    """⚠ The guard must not break re-ingesting the same model — that is the ordinary path."""
    ordered = ["SAMEA0000001"]
    roster_module, artifacts = _roster_of(monkeypatch, tmp_path, ordered)
    first, _ = roster_module.ingest_genome_collection(
        session, artifacts, pathogen_species_id=seeded["species"].pathogen_species_id
    )
    session.flush()
    again, _ = roster_module.ingest_genome_collection(
        session, artifacts, pathogen_species_id=seeded["species"].pathogen_species_id
    )
    assert first == again


# ------------------------------------------------- what publishing must NOT do, and the key guard


def test_publishing_one_catalogue_does_NOT_hide_the_other(session, seeded):
    """⛔⛔ The defect that defeated the whole feature: publishing nuna5 deleted nuna4 from the picker.

    `publish_pangenome` used to clear `is_published` on whatever the species pointer previously
    named. That was harmless while the column had no reader and a species had one catalogue — it
    retired a superseded GENERATION. Once `is_published` became picker visibility and the pointer's
    previous value became a different MODEL, the same line meant "promoting nuna5 removes nuna4",
    with nothing in any output saying so.

    ⚠ Publishing ADDS. Hiding a catalogue is a separate, explicit act.
    """
    from bacatlas_backend.ingest.publish_pangenome import publish_pangenome

    seeded["pangenome"].is_published = True
    seeded["species"].default_pangenome_id = seeded["pangenome"].pangenome_id
    second = _second_catalogue(session, seeded)
    session.flush()
    assert {p.catalogue_key for p, _s, _m in list_catalogues(session)} == {"probe-nuna4", "probe-nuna5"}

    publish_pangenome(session, run_id=second.run_id, force=True)
    session.flush()
    session.refresh(seeded["species"])

    offered = {p.catalogue_key for p, _s, _m in list_catalogues(session)}
    assert "probe-nuna5" in offered, "the newly published catalogue is offered"
    assert "probe-nuna4" in offered, "and the comparator the picker exists to show is STILL offered"
    assert seeded["species"].default_pangenome_id == second.pangenome_id, "the pointer still moved"


def test_a_key_with_NO_HYPHEN_is_refused_because_the_resolver_reads_the_hyphen(session, seeded):
    """⛔ The API tells a catalogue key from a species key by the hyphen, so a hyphenless key is
    either unreachable or — worse — silently resolves to the species DEFAULT instead of itself.

    `--catalogue-key ecoli` was accepted by the old guard, which exempted a key equal to its species.
    """
    with pytest.raises(Exception, match="no hyphen"):
        catalogue_key_for("ecoli", "nuna5", "sensitive")
    with pytest.raises(Exception, match="no hyphen"):
        catalogue_key_for("ecoli", "nuna5", "ecoli")


def test_an_underscore_ANYWHERE_is_refused_not_just_the_species_prefix():
    """⛔ Prefix-freeness is STRUCTURAL, not a comparison against whatever keys exist today.

    The payload glob is `locus_browser_<key>_*.json`, so key A captures key B exactly when B starts
    with `A + "_"` — which requires an underscore inside B. Forbidding the character makes the
    property hold for every key that will ever be added, in any order.

    The earlier guard checked one key against one species key and so still admitted `ecoli-nuna5`
    beside `ecoli-nuna5_damped`: the shorter key's glob matches the longer key's payload and, being
    the newer export, wins.
    """
    with pytest.raises(Exception, match="underscore"):
        catalogue_key_for("ecoli", "nuna5", "ecoli-nuna5_damped")
    with pytest.raises(Exception, match="underscore"):
        catalogue_key_for("kp", "nuna5", "ecoli_nuna5")


def test_every_key_the_registry_can_produce_satisfies_both_rules():
    """⚠ The default must never trip its own guard — check it against nuna's real model keys."""
    model_keys = [
        "nuna4", "nuna5", "nuna5damped", "nuna5logp02", "nuna5k250",
        "nuna5contk250", "nuna5logp02cont", "nuna6", "nuna6d", "nuna5esm097", "nuna5esm096",
    ]
    for species in ("ecoli", "kp"):
        for model in model_keys:
            key = catalogue_key_for(species, model)
            assert "-" in key and "_" not in key, key


def test_offering_a_catalogue_does_not_make_it_the_default(session, seeded):
    """⭐ Three separate decisions: LOAD it, OFFER it, SERVE it by default.

    They were one act only while a species had one catalogue. Keeping them apart is what lets nuna4
    stay the default and the rollback while nuna5 is selectable beside it — which is the shape the
    user asked for. Publishing was the only verb available, and publishing moves the pointer.
    """
    from bacatlas_backend.ingest.publish_pangenome import offer_catalogue

    seeded["species"].default_pangenome_id = seeded["pangenome"].pangenome_id
    second = _second_catalogue(session, seeded)
    second.is_published = False
    session.flush()
    assert "probe-nuna5" not in {p.catalogue_key for p, _s, _m in list_catalogues(session)}

    offer_catalogue(session, catalogue_key="probe-nuna5")
    session.flush()
    session.refresh(seeded["species"])

    assert "probe-nuna5" in {p.catalogue_key for p, _s, _m in list_catalogues(session)}
    assert seeded["species"].default_pangenome_id == seeded["pangenome"].pangenome_id, (
        "offering must not promote — the default is a separate decision"
    )


def test_a_catalogue_can_be_taken_back_out_of_the_picker(session, seeded):
    """And hiding is the same one boolean, not a delete — the catalogue stays addressable by key."""
    from bacatlas_backend.ingest.publish_pangenome import offer_catalogue

    second = _second_catalogue(session, seeded)
    session.flush()
    offer_catalogue(session, catalogue_key="probe-nuna5", offered=False)
    session.flush()

    assert "probe-nuna5" not in {p.catalogue_key for p, _s, _m in list_catalogues(session)}
    assert resolve_catalogue(session, "probe-nuna5").pangenome_id == second.pangenome_id


# ── the cohort label is what the picker SHOWS, so it is gated at the picker ───────────────────


def test_a_cohort_still_wearing_its_KEY_is_REFUSED_from_the_picker(session, seeded):
    """⛔ The picker's options are cohorts now, so a missing label is a reader-visible defect.

    The label defaults to the collection key — a machine token with no space in it — so without
    this guard `--offer` would put `kp-71k` in a menu where "Kp" belongs, and nothing downstream
    could tell that from a label somebody meant. Loading and addressing by key stay allowed; only
    showing it is refused.
    """
    from bacatlas_backend.ingest.publish_pangenome import PublishRefused, offer_catalogue

    collection = seeded["collection"]
    collection.collection_label = collection.collection_key
    _second_catalogue(session, seeded).is_published = False
    session.flush()

    with pytest.raises(PublishRefused, match="still its own key"):
        offer_catalogue(session, catalogue_key="probe-nuna5")

    session.flush()
    assert "probe-nuna5" not in {p.catalogue_key for p, _s, _m in list_catalogues(session)}, (
        "the refusal must leave the picker unchanged, not raise after writing the boolean"
    )


def test_the_SAME_catalogue_is_offered_once_its_cohort_has_a_label(session, seeded):
    """The guard's other side — otherwise it could refuse everything and still pass above."""
    from bacatlas_backend.ingest.publish_pangenome import offer_catalogue

    collection = seeded["collection"]
    collection.collection_label = collection.collection_key
    second = _second_catalogue(session, seeded)
    second.is_published = False
    session.flush()

    collection.collection_label = "E. coli probe"
    session.flush()
    offer_catalogue(session, catalogue_key="probe-nuna5")
    session.flush()

    assert "probe-nuna5" in {p.catalogue_key for p, _s, _m in list_catalogues(session)}


def test_HIDING_a_catalogue_is_never_refused_for_want_of_a_label(session, seeded):
    """⚠ A guard on the way in must not become a trap on the way out.

    Taking something out of the picker is the remedy for a bad label, so refusing it when the label
    is bad would leave a catalogue stuck on display with no way back.
    """
    from bacatlas_backend.ingest.publish_pangenome import offer_catalogue

    collection = seeded["collection"]
    second = _second_catalogue(session, seeded)
    session.flush()
    assert "probe-nuna5" in {p.catalogue_key for p, _s, _m in list_catalogues(session)}

    collection.collection_label = collection.collection_key
    session.flush()
    offer_catalogue(session, catalogue_key="probe-nuna5", offered=False)
    session.flush()

    assert "probe-nuna5" not in {p.catalogue_key for p, _s, _m in list_catalogues(session)}
    assert second.is_published is False
