"""The roster: the ordinal `arr.gid` indexes into, and the two ways it can be silently wrong.

Both failures here are invisible in the data — an ordering that is one place out renames every
genome on the page, and a collection with a hole in it does the same from the first gap onward.
Neither raises anywhere downstream, so each is asserted at the point it is created.
"""

import pytest
from sqlalchemy import select

from bacatlas_backend.ingest.artifact_locator import CatalogueArtifacts
from bacatlas_backend.ingest.ingest_genome_collection_roster import (
    RosterError,
    genome_id_by_ordinal,
    ingest_genome_collection,
    read_genome_vocabulary,
)
from bacatlas_backend.ingest.ingest_pathogen_species import ingest_pathogen_species
from bacatlas_backend.models.enumerations import SampleIdentifierKind
from bacatlas_backend.models.genome import Genome
from bacatlas_backend.models.genome_collection import GenomeCollection, GenomeCollectionMembership

pandas = pytest.importorskip("pandas", reason="the ingest extra")


def _stage_genome_rows(session, species_id, sample_ids):
    """The genome layer, reduced to what the roster join needs — one row per accession."""
    session.add_all(
        [
            Genome(
                pathogen_species_id=species_id,
                sample_id=sample_id,
                sample_id_kind=SampleIdentifierKind.BIOSAMPLE,
            )
            for sample_id in sample_ids
        ]
    )
    session.flush()


def _fake_universe(tmp_path, sample_ids, *, set_key="probe", species_key=None):
    """A `{set}_close_seq.parquet` with a chosen vocabulary, in the real directory shape.

    ⚠ `set_key` is a parameter because the cohort declaration is keyed on `(species, set_key)`:
    `probe` is deliberately UNDECLARED, so a load of it must name its own key and label, while
    `ecoli` resolves through `published_cohorts` and is how the declared path is exercised.
    """
    root = tmp_path / "data"
    (root / "proc" / "mmseqs").mkdir(parents=True)
    pandas.DataFrame(
        {"sample_id": sample_ids, "flat_index": list(range(len(sample_ids)))}
    ).to_parquet(root / "proc" / "mmseqs" / f"{set_key}_close_seq.parquet")
    return CatalogueArtifacts(
        data_root=root,
        set_key=set_key,
        model_label=f"{set_key}_model",
        run_id=f"{set_key}_run",
        species_key=species_key,
    )


# ── the vocabulary ─────────────────────────────────────────────────────────────────────────────
def test_the_vocabulary_is_the_UNIVERSE_sorted_and_it_matches_the_published_counts(
    ecoli_artifacts, published_ecoli_site_catalogue
):
    """⛔ `meta.genomes` is a published contract. This is the list `arr.gid` indexes into."""
    samples, n_genes = read_genome_vocabulary(ecoli_artifacts)
    assert samples == published_ecoli_site_catalogue["meta"]["genomes"]
    assert len(samples) == published_ecoli_site_catalogue["meta"]["n_genomes"] == 100
    assert n_genes == published_ecoli_site_catalogue["meta"]["n_genes"] == 489_146


def test_the_vocabulary_is_sorted_and_distinct_and_that_is_asserted_not_assumed(ecoli_artifacts):
    samples, _ = read_genome_vocabulary(ecoli_artifacts)
    assert samples == sorted(samples)
    assert len(set(samples)) == len(samples)


# ── the two silent failures ────────────────────────────────────────────────────────────────────
def test_a_roster_genome_with_no_genome_row_is_REFUSED_rather_than_skipped(session, tmp_path):
    """A hole would shift every membership after it by one, and nothing downstream could tell."""
    species_ids = ingest_pathogen_species(session)
    session.flush()
    _stage_genome_rows(session, species_ids["ecoli"], ["SAM_A", "SAM_B"])
    artifacts = _fake_universe(tmp_path, ["SAM_A", "SAM_B", "SAM_MISSING"])
    with pytest.raises(RosterError, match="no `genome` row"):
        ingest_genome_collection(
            session, artifacts, pathogen_species_id=species_ids["ecoli"],
            collection_key="probe100", collection_label="E. coli probe",
        )


def test_ordinals_with_a_hole_in_them_are_refused_by_the_READER_too(session, seeded):
    """Written straight to the table, because that is the state a partial delete would leave."""
    collection = session.get(GenomeCollection, seeded["collection"].genome_collection_id)
    session.add(
        GenomeCollectionMembership(
            genome_collection_id=collection.genome_collection_id,
            genome_id=seeded["genome"].genome_id,
            collection_genome_ordinal=7,
        )
    )
    session.flush()
    with pytest.raises(RosterError, match="not 0"):
        genome_id_by_ordinal(session, collection.genome_collection_id)


# ── the ordering, end to end ───────────────────────────────────────────────────────────────────
def test_the_ordinal_follows_the_SORTED_vocabulary_and_not_the_file_order(session, tmp_path):
    """The parquet is written deliberately out of order; the ordinals must not inherit it."""
    species_ids = ingest_pathogen_species(session)
    session.flush()
    _stage_genome_rows(session, species_ids["ecoli"], ["SAM_C", "SAM_A", "SAM_B"])
    artifacts = _fake_universe(tmp_path, ["SAM_C", "SAM_A", "SAM_B"])
    collection_id, report = ingest_genome_collection(
        session, artifacts, pathogen_species_id=species_ids["ecoli"],
        collection_key="probe100", collection_label="E. coli probe",
    )
    session.flush()
    assert report.genome_count == 3
    ordered = session.execute(
        select(GenomeCollectionMembership.requested_sample_id)
        .where(GenomeCollectionMembership.genome_collection_id == collection_id)
        .order_by(GenomeCollectionMembership.collection_genome_ordinal)
    ).scalars().all()
    assert ordered == ["SAM_A", "SAM_B", "SAM_C"]


def test_the_real_ecoli_roster_reproduces_meta_genomes_and_reloads_without_duplicating(
    session, ecoli_artifacts, published_ecoli_site_catalogue
):
    species_ids = ingest_pathogen_species(session)
    session.flush()
    samples, _ = read_genome_vocabulary(ecoli_artifacts)
    _stage_genome_rows(session, species_ids["ecoli"], samples)

    first, report = ingest_genome_collection(
        session, ecoli_artifacts, pathogen_species_id=species_ids["ecoli"]
    )  # ⭐ declared: `(ecoli, ecoli)` resolves through `published_cohorts`, no flags needed
    session.flush()
    assert (report.genome_count, report.members_written, report.unresolved_sample_ids) == (100, 100, [])

    second, _ = ingest_genome_collection(
        session, ecoli_artifacts, pathogen_species_id=species_ids["ecoli"]
    )
    session.flush()
    assert second == first, "a re-ingest must keep the collection id — pangenomes point at it"

    ordered = session.execute(
        select(GenomeCollectionMembership.requested_sample_id)
        .where(GenomeCollectionMembership.genome_collection_id == first)
        .order_by(GenomeCollectionMembership.collection_genome_ordinal)
    ).scalars().all()
    assert len(ordered) == 100, "membership must be REPLACED wholesale, never merged"
    assert ordered == published_ecoli_site_catalogue["meta"]["genomes"]
    assert len(genome_id_by_ordinal(session, first)) == 100


# ── the cohort's LABEL — what the picker shows now the cohort is the axis ─────────────────────
#
# ⛔ These exist because the label is the first field of this table a reader ever sees, and the two
# ways it can be wrong are both silent: it can be absent (and the picker shows a machine token), or
# it can state a genome count that the roster it names contradicts.


def _labelled(session, tmp_path, samples, **kwargs):
    """A collection over a synthetic roster — the label is the subject, so keep the roster small."""
    species_ids = ingest_pathogen_species(session)
    session.flush()
    _stage_genome_rows(session, species_ids["ecoli"], samples)
    collection_id, _ = ingest_genome_collection(
        session, _fake_universe(tmp_path, samples), pathogen_species_id=species_ids["ecoli"],
        **{"collection_key": "probe100", "collection_label": "E. coli probe", **kwargs},
    )
    session.flush()
    return collection_id


def test_an_UNDECLARED_cohort_with_no_key_and_no_label_is_REFUSED_rather_than_invented(
    session, tmp_path
):
    """⛔ The premise of the test this replaced is now false, and the replacement is stricter.

    It used to assert that an unlabelled cohort wears its KEY — visible, unfinished, and refused
    by `offer_catalogue` on the way into the picker. Two things were wrong with relying on that.
    The key itself defaulted to `artifacts.set_key`, so a no-flags load wrote the catalogue key
    `{species}-{set_key}` — `ecoli-ecoli` — which passes every guard and is nobody's cohort. And
    the label is assigned UNCONDITIONALLY on an existing row, so a name a human chose was put back
    to the machine token by the next load that omitted the flag, while `offer_catalogue`'s guard
    fires only on the way IN and an already-offered catalogue never re-checks.

    Both defaults now come from `published_cohorts`, and an undeclared load with neither flag is
    refused by name instead of being invented.
    """
    from bacatlas_backend.ingest.published_cohorts import UndeclaredCohort

    species_ids = ingest_pathogen_species(session)
    session.flush()
    _stage_genome_rows(session, species_ids["ecoli"], ["SAM_A", "SAM_B"])
    with pytest.raises(UndeclaredCohort, match="no cohort declared"):
        ingest_genome_collection(
            session,
            _fake_universe(tmp_path, ["SAM_A", "SAM_B"]),
            pathogen_species_id=species_ids["ecoli"],
        )


def test_a_DECLARED_cohort_takes_the_label_from_the_declaration_and_NOT_its_key(session, tmp_path):
    """⭐ The other side, and the reason the declaration is code rather than a comment.

    `(ecoli, ecoli)` is declared, so a load that passes no flags at all gets the words David chose
    — and gets them again on every later re-ingest, which is what makes a chosen label survive.
    """
    species_ids = ingest_pathogen_species(session)
    session.flush()
    _stage_genome_rows(session, species_ids["ecoli"], ["SAM_A", "SAM_B"])
    collection_id, report = ingest_genome_collection(
        session,
        _fake_universe(tmp_path, ["SAM_A", "SAM_B"], set_key="ecoli", species_key="ecoli"),
        pathogen_species_id=species_ids["ecoli"],
    )
    session.flush()
    collection = session.get(GenomeCollection, collection_id)
    assert (collection.collection_key, collection.collection_label) == ("probe100", "E. coli AMR")
    assert collection.collection_label != collection.collection_key
    assert report.collection_label == "E. coli AMR", "and the operator is told the words, not the key"
    # ⚠ Still no count in the label — the page composes that from `genome_count` beside it.
    assert str(collection.genome_count) not in collection.collection_label


def test_a_supplied_label_is_stored_VERBATIM_and_carries_no_count(session, tmp_path):
    collection = session.get(
        GenomeCollection,
        _labelled(
            session,
            tmp_path,
            ["SAM_A", "SAM_B", "SAM_C"],
            collection_key="probe100",
            collection_label="E. coli probe",
        ),
    )
    assert (collection.collection_key, collection.collection_label) == ("probe100", "E. coli probe")
    assert collection.genome_count == 3
    # ⚠ The count is NOT in the label. The page renders `{label} — {count} genomes` from the two
    # fields, so a label carrying its own number could contradict the roster it names.
    assert str(collection.genome_count) not in collection.collection_label


def test_a_re_ingest_can_CHANGE_the_label_without_renumbering_one_genome(session, tmp_path):
    """Correcting a word on the page must not look like a roster change.

    ⚠ The roster guard refuses a membership reorder while a pangenome depends on the ordinals, so
    if the label rode on the membership write a typo would become unfixable on a live catalogue.
    """
    samples = ["SAM_A", "SAM_B", "SAM_C"]
    species_ids = ingest_pathogen_species(session)
    session.flush()
    _stage_genome_rows(session, species_ids["ecoli"], samples)
    # ⚠ Built ONCE: `_fake_universe` creates its own directory tree and cannot be called twice on
    # one tmp_path. Reusing it is also the honest shape — the roster file did not change.
    artifacts = _fake_universe(tmp_path, samples)
    first, _ = ingest_genome_collection(
        session, artifacts, pathogen_species_id=species_ids["ecoli"],
        collection_key="probe100", collection_label="first",
    )
    session.flush()
    members = lambda cid: session.execute(  # noqa: E731
        select(GenomeCollectionMembership.genome_id, GenomeCollectionMembership.collection_genome_ordinal)
        .where(GenomeCollectionMembership.genome_collection_id == cid)
        .order_by(GenomeCollectionMembership.collection_genome_ordinal)
    ).all()
    before = members(first)

    second, _ = ingest_genome_collection(
        session, artifacts, pathogen_species_id=species_ids["ecoli"],
        collection_key="probe100", collection_label="second",
    )
    session.flush()

    assert second == first, "the same roster must stay the same collection"
    assert session.get(GenomeCollection, second).collection_label == "second"
    assert members(second) == before, "changing a word on the page must not renumber a genome"
    assert len(before) == len(samples) > 0, "and that comparison is over a roster, not two empty lists"
