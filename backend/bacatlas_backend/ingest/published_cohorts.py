"""The cohorts the picker offers — and the one table that tells the FIVE key vocabularies apart.

⛔⛔ **FIVE KEYS, FIVE JOBS, AND NONE IS DERIVABLE FROM ANOTHER.** `published_catalogues.py` warns
about three; the cohort axis adds a fourth and a fifth, and the 71k cohort is the row on which all
five differ:

    pathogen_species.species_key    kp         the organism
    CatalogueArtifacts.set_key      kp71k      nuna's dataset token IN FILENAMES
    genome_collection.collection_key  71k      the roster, scoped to its species
    pangenome.catalogue_key         kp-71k     the service's address, what a URL carries
    PublishedCatalogue.dset         (none)     the static page — only the two nuna4 rows have one

`set_key` is nuna's own (`build_kp_full_rosters.py`: `KP71K, KP13K, HI_SET = "kp71k", "kp13k",
"kp_HI_kp1000"`), and it is the token every artifact path is built from. Collapsing any two of the
five is how a load silently reads another cohort's files.

⛔ **THIS IS LOAD-BEARING, NOT DOCUMENTATION — and that is the whole reason it exists rather than a
comment.** `ingest_genome_collection` used to default `collection_label` to the collection KEY, and
it assigns unconditionally on an existing row. So a label a human chose survived exactly until the
next `--stage pangenome` run that omitted `--collection-label`, which would put the machine token
back into the picker — and `offer_catalogue`'s guard fires only on the way IN, so an
already-published catalogue would never re-check and the token would reach a reader with nothing
refusing it. The default now comes from here.

⛔ **And it is why `--collection-key` can be omitted safely.** That flag used to default to
`artifacts.set_key`, which made the no-flags catalogue key `{species}-{set_key}` — `ecoli-ecoli`,
`kp-kp` — legal under all three of `catalogue_key_for`'s guards and silently wrong. The checked-in
ingest recipe in `backend/README.md` passes neither flag, so following it verbatim produced exactly
that. Resolving `(species_key, set_key)` through this table instead means an undeclared cohort is
REFUSED rather than invented.

⚠ **Separate from `published_catalogues.py` on purpose.** That records two EXPORTS THAT HAPPENED,
every field of which is checkable against files on this laptop and is enforced by two tests. This
records cohorts we INTEND TO SERVE, one of which has no run, no model, no artifacts and no oracle
until the gene layer is compressed. Merging them means weakening those tests or carrying a row that
cannot satisfy its own dataclass. They also move at different times: a published triple must never
change — it is the rollback and the parity oracle — while a cohort's label changes the moment
better words are chosen.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Cohort:
    """One roster the browser offers, and every key that addresses it."""

    #: ⭐ What the picker SHOWS. David, 2026-10-27: *"E.coli AMR 100, Kp AMR 100, Kp 71,500 will
    #: do"* — stored without the counts, because the page renders `"{label} — {count:,} genomes"`
    #: from this and `genome_collection.genome_count`, which is *passed, never inferred*. Stored
    #: verbatim, "Kp 71,500" would render "Kp 71,500 — 71,435 genomes": a label contradicting the
    #: roster it names, with nothing able to say which half was right.
    label: str

    species_key: str
    #: Scoped to the species, so two species may both have a `probe100`.
    collection_key: str
    #: ⛔ nuna's dataset token in FILENAMES, which is not the collection key and not the species
    #: key. `kp71k` for the 71k cohort — a fourth spelling of the same cohort.
    set_key: str
    #: `{species_key}-{collection_key}` — kept explicit rather than composed, because
    #: `catalogue_key_for` admits a deliberate override and a rule that is usually right is the
    #: kind this repo's docstrings keep warning about.
    catalogue_key: str
    #: ⚠ `None` where the roster has not been counted into a database yet. **Not zero** — a cohort
    #: of unknown size and a cohort of no genomes are different claims, and only the first is true
    #: of `kp-71k` here. `reconcile_cohort_genome_count` compares where there is something to
    #: compare and says so where there is not.
    genome_count: int | None


#: ⚠ `kp-71k`'s 71,435 is nuna's roster count and is **not yet verified against any database**: at
#: the measured per-row cost the gene layer alone is ~340 GB at that size, against an 80 GB
#: instance, so it cannot be ingested until the compression is finished. It is carried here so the
#: key vocabulary is complete and the number has somewhere to be reconciled against later.
PUBLISHED_COHORTS: tuple[Cohort, ...] = (
    Cohort(
        label="E. coli AMR",
        species_key="ecoli",
        collection_key="probe100",
        set_key="ecoli",
        catalogue_key="ecoli-probe100",
        genome_count=100,
    ),
    Cohort(
        label="Kp AMR",
        species_key="kp",
        collection_key="probe100",
        set_key="kp",
        catalogue_key="kp-probe100",
        genome_count=100,
    ),
    Cohort(
        label="Kp",
        species_key="kp",
        collection_key="71k",
        set_key="kp71k",
        catalogue_key="kp-71k",
        genome_count=None,
    ),
)


class UndeclaredCohort(LookupError):
    """A load named a (species, dataset) pair that no cohort here describes."""


def cohort_for(species_key: str, set_key: str) -> Cohort:
    """The cohort a load of `(species_key, set_key)` belongs to, or a failure naming what is known.

    ⛔ **Refuses rather than inventing one.** The alternative — falling back to the set key as the
    collection key — is what produced `ecoli-ecoli`: a catalogue key that passes every guard, reads
    like a typo nobody made, and addresses a cohort that does not exist as a cohort.
    """
    for cohort in PUBLISHED_COHORTS:
        if cohort.species_key == species_key and cohort.set_key == set_key:
            return cohort
    known = [(c.species_key, c.set_key) for c in PUBLISHED_COHORTS]
    raise UndeclaredCohort(
        f"no cohort declared for species {species_key!r} + dataset {set_key!r}; known: {known}. "
        f"Add it to `published_cohorts.py` — a cohort reaches a reader by its label, and a load "
        f"that invents one puts a machine token in the picker. Pass --collection-key and "
        f"--collection-label explicitly only to load something that is deliberately not offered."
    )


def cohort_by_collection(species_key: str, collection_key: str) -> Cohort | None:
    """The declared cohort with this roster key, or `None`. ⚠ `None` is a legitimate answer here.

    Used where a cohort may legitimately be undeclared — a staged roster loaded by hand with
    `--collection-key` — so it returns rather than raises, and the caller decides.
    """
    for cohort in PUBLISHED_COHORTS:
        if cohort.species_key == species_key and cohort.collection_key == collection_key:
            return cohort
    return None
