"""Flip a species onto a pangenome — the one-row update that decides what the service serves.

⭐ **This is what makes a re-ingest safe on a running service.** Ingest builds a whole generation
alongside the live one and never touches the pointer; this flips it, in its own transaction, once the
new catalogue is complete. A rollback is the same update backwards, which is why the old generation
is dropped only after the new one has served.

⛔ **It VERIFIES before it flips, and the verification is not a formality.** Publishing an incomplete
catalogue is the one failure this design can produce that a reader cannot see: every page renders,
every number is in range, and a block is simply missing. The checks below are each a thing that has
either gone wrong before or would be invisible if it did — a locus count that disagrees with the
row count, a landing locus that was never written (the render-time mutation that shipped five days
of pages opening on locus 0), a species whose loci belong to a different species, a geometry that
covers one representation and not both.

⚠ **`--force` exists and says what it is.** A catalogue can be published with a known gap — the local
mirror has no gate ledger, for instance — but the gap is named in the return value rather than
silently tolerated.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from sqlalchemy import exists, func, select, update
from sqlalchemy.orm import Session

from bacatlas_backend.models.enumerations import EmbeddingRepresentation
from bacatlas_backend.models.genome_collection import GenomeCollection, GenomeCollectionMembership
from bacatlas_backend.models.locus import Locus
from bacatlas_backend.models.locus_arrangement import LocusArrangement
from bacatlas_backend.models.locus_offset_occupant import LocusOffsetOccupant
from bacatlas_backend.models.locus_similarity import (
    LocusSimilarity,
    PangenomeSimilarityBaseline,
)
from bacatlas_backend.models.pangenome import Pangenome
from bacatlas_backend.models.pangenome_genome_locus_count import PangenomeGenomeLocusCount
from bacatlas_backend.models.pathogen_species import PathogenSpecies


class PublishRefused(RuntimeError):
    """The catalogue is not complete enough to serve, and the failures are named."""


@dataclass
class PublishReport:
    """What was checked, what it found, and what the pointer now is."""

    run_id: str
    pangenome_id: int
    species_key: str
    previous_pangenome_id: int | None
    checks_passed: list = field(default_factory=list)
    failures: list = field(default_factory=list)
    published: bool = False

    def render(self) -> str:
        """Every check by name — a report that says only "published" has told you nothing."""
        lines = [f"publish {self.species_key} → pangenome {self.pangenome_id} ({self.run_id})"]
        for check in self.checks_passed:
            lines.append(f"  ✓ {check}")
        for failure in self.failures:
            lines.append(f"  ✗ {failure}")
        lines.append(
            f"  pointer: {self.previous_pangenome_id} → "
            f"{self.pangenome_id if self.published else self.previous_pangenome_id}"
        )
        return "\n".join(lines)


def verify_pangenome_is_servable(session: Session, pangenome: Pangenome) -> tuple[list[str], list[str]]:
    """`(passed, failed)` — every check by name, so a partial catalogue cannot pass quietly."""
    passed: list[str] = []
    failed: list[str] = []

    def check(name: str, condition: bool, detail: str = "") -> None:
        (passed if condition else failed).append(name if condition else f"{name}: {detail}")

    locus_count = session.execute(
        select(func.count()).select_from(Locus).where(Locus.pangenome_id == pangenome.pangenome_id)
    ).scalar_one()
    check(
        f"{locus_count:,} loci, matching the recorded locus_count",
        locus_count == pangenome.locus_count and locus_count > 0,
        f"the table holds {locus_count:,} and the row claims {pangenome.locus_count:,}",
    )

    # ⛔ Every locus must belong to the SAME species as the pangenome. A mismatch here would serve
    # one species' loci under another's key, and every one of them would render perfectly.
    foreign = session.execute(
        select(func.count())
        .select_from(Locus)
        .where(
            Locus.pangenome_id == pangenome.pangenome_id,
            Locus.pathogen_species_id != pangenome.pathogen_species_id,
        )
    ).scalar_one()
    check(
        "every locus belongs to this pangenome's species",
        foreign == 0,
        f"{foreign:,} loci carry a different pathogen_species_id",
    )

    # ⚠ The render-time mutation that shipped five days of pages opening on locus 0.
    check(
        "a landing locus is written",
        pangenome.landing_locus_id is not None,
        "landing_locus_id is NULL — the page would open on whatever came first",
    )
    check(
        "example loci are written",
        bool(pangenome.example_locus_ids),
        "example_locus_ids is empty — the chips would be blank",
    )

    arrangements = session.execute(
        select(func.count())
        .select_from(LocusArrangement)
        .where(LocusArrangement.pangenome_id == pangenome.pangenome_id)
    ).scalar_one()
    check("the joint view is populated", arrangements > 0, "no arrangements — the track cannot draw")

    occupants = session.execute(
        select(func.count())
        .select_from(LocusOffsetOccupant)
        .where(LocusOffsetOccupant.pangenome_id == pangenome.pangenome_id)
    ).scalar_one()
    check("the marginal view is populated", occupants > 0, "no offset occupants — no bars under the track")

    # ⭐ BOTH representations, or neither half of the card means what it says: ESM is homology and
    # Bacformer is context, and they deliberately disagree about which loci are weak. Measured: of
    # the 1,059 ecoli loci Bacformer flags, ESM rescues 817 (77 %) — so a missing ESM baseline is not
    # half a card, it is the half that says the other half is wrong.
    baselines = {
        row.representation: row
        for row in session.execute(
            select(PangenomeSimilarityBaseline).where(
                PangenomeSimilarityBaseline.pangenome_id == pangenome.pangenome_id
            )
        ).scalars()
    }
    check(
        "both similarity representations are present",
        set(baselines) == set(EmbeddingRepresentation),
        f"only {sorted(r.value for r in baselines)} — the other half of the card would be empty",
    )

    # ⭐ And a FLOOR for each, because without it a cosine has no scale: ESM's random gene pairs sit
    # at ~0.742 and Bacformer's at ~0.059, so the same 0.41 reads oppositely in the two. A card that
    # cannot draw the baseline is a card printing numbers a reader has no way to size.
    for representation, baseline in sorted(baselines.items(), key=lambda item: item[0].value):
        check(
            f"the {representation.value} baseline has a floor",
            baseline.floor_median is not None,
            "no random gene-pair median — the card would print cosines against nothing",
        )

    # ⛔ The denominator is checked against the ROWS, not taken on trust. "p12 of 12,104 loci" is two
    # assertions, and the second one is this number; a baseline that outlived a partial similarity
    # load would print a percentile over a population that is not there.
    for representation, baseline in sorted(baselines.items(), key=lambda item: item[0].value):
        measurable = session.execute(
            select(func.count())
            .select_from(LocusSimilarity)
            .join(Locus, Locus.locus_id == LocusSimilarity.locus_id)
            .where(
                Locus.pangenome_id == pangenome.pangenome_id,
                LocusSimilarity.representation == representation,
                LocusSimilarity.separation_percentile.is_not(None),
            )
        ).scalar_one()
        check(
            f"the {representation.value} denominator matches the rows behind it",
            measurable == baseline.measurable_locus_count,
            f"{measurable:,} ranked loci against a stated {baseline.measurable_locus_count:,}",
        )

    unnamed = session.execute(
        select(func.count())
        .select_from(Locus)
        .where(Locus.pangenome_id == pangenome.pangenome_id, Locus.display_name.is_(None))
    ).scalar_one()
    check("every locus has a display name", unnamed == 0, f"{unnamed:,} loci have none")

    # ⭐ The anchor picker lists genomes FROM this table, so a collection genome with no row is not
    # "a genome in nothing" — it is a genome the picker silently cannot offer, and a stray row is a
    # genome it offers that this catalogue never saw. Both directions, in one statement over the
    # collection, not over the 412 M membership rows the counts were derived from.
    counted = select(PangenomeGenomeLocusCount.genome_id).where(
        PangenomeGenomeLocusCount.pangenome_id == pangenome.pangenome_id
    )
    members = select(GenomeCollectionMembership.genome_id).where(
        GenomeCollectionMembership.genome_collection_id == pangenome.genome_collection_id
    )
    uncounted, stray = session.execute(
        select(
            select(func.count())
            .select_from(GenomeCollectionMembership)
            .where(
                GenomeCollectionMembership.genome_collection_id == pangenome.genome_collection_id,
                ~exists(
                    counted.where(
                        PangenomeGenomeLocusCount.genome_id == GenomeCollectionMembership.genome_id
                    )
                ),
            )
            .scalar_subquery(),
            select(func.count())
            .select_from(PangenomeGenomeLocusCount)
            .where(
                PangenomeGenomeLocusCount.pangenome_id == pangenome.pangenome_id,
                ~exists(
                    members.where(
                        GenomeCollectionMembership.genome_id == PangenomeGenomeLocusCount.genome_id
                    )
                ),
            )
            .scalar_subquery(),
        )
    ).one()
    check(
        "every collection genome has a locus count row",
        uncounted == 0 and stray == 0,
        f"{uncounted:,} collection genomes have no row and {stray:,} rows name a genome outside the "
        "collection — the anchor picker would drop or invent them",
    )
    return passed, failed


def offer_catalogue(session: Session, *, catalogue_key: str, offered: bool = True) -> str:
    """Show a catalogue in the picker, or take it out — WITHOUT touching any default pointer.

    ⭐ **The verb that `publish` is not.** Publishing decides what a bare species request serves;
    this decides what the picker offers. They were one act only while a species had one catalogue.
    A model can now be loaded, addressed by key for review, and shown to readers as three separate
    decisions — which is what lets nuna4 stay the default and the rollback while nuna5 is selectable
    beside it.

    ⚠ It deliberately cannot make a catalogue the default. Promoting is `publish_pangenome`, which
    verifies the whole catalogue first; this writes one boolean and checks **one** thing, because a
    catalogue that is already loaded and addressable has otherwise been verified by whatever put it
    there.

    ⛔ **The one check: a cohort may not reach the picker wearing its key.** The picker's options
    are cohorts now, so `collection_label` is read by a person — and it defaults to
    `collection_key`, which is a machine token with no space in it. Offering such a catalogue would
    put `kp-71k` in a menu where "Kp" belongs, and nothing downstream could tell that from a label
    somebody meant. Loading and addressing it by key stay allowed; only showing it is refused.
    """
    pangenome = session.execute(
        select(Pangenome).where(Pangenome.catalogue_key == catalogue_key)
    ).scalar_one_or_none()
    if pangenome is None:
        raise PublishRefused(f"no catalogue {catalogue_key!r} — nothing to offer.")
    if offered:
        collection = session.get(GenomeCollection, pangenome.genome_collection_id)
        if collection is not None and collection.collection_label == collection.collection_key:
            raise PublishRefused(
                f"{catalogue_key!r} sits on cohort {collection.collection_key!r}, whose label is "
                f"still its own key. The picker shows the LABEL, so this would offer a machine "
                f"token where a cohort name belongs. Re-ingest with "
                f"`--collection-label '<what a reader should see>'`, then offer it. "
                f"⚠ Do not put the genome count in the label — the page renders that from the "
                f"roster's own count ({collection.genome_count:,})."
            )
    was = pangenome.is_published
    session.execute(
        update(Pangenome)
        .where(Pangenome.pangenome_id == pangenome.pangenome_id)
        .values(is_published=offered)
    )
    verb = "offered in the picker" if offered else "hidden from the picker"
    return f"{catalogue_key}: {verb}" + ("" if was != offered else " (already was)")


def publish_pangenome(
    session: Session, *, run_id: str, ingest_generation: int = 1, force: bool = False
) -> PublishReport:
    """Verify the catalogue, then point its species at it. Refuses on any failed check.

    ⛔ **Its own transaction, and only the pointer.** Nothing else is written here, so a publish can
    never be half-applied and a rollback is one row.
    """
    pangenome = session.execute(
        select(Pangenome).where(
            Pangenome.run_id == run_id, Pangenome.ingest_generation == ingest_generation
        )
    ).scalar_one_or_none()
    if pangenome is None:
        raise PublishRefused(f"no pangenome {run_id!r} at generation {ingest_generation}")

    species = session.get(PathogenSpecies, pangenome.pathogen_species_id)
    report = PublishReport(
        run_id=run_id,
        pangenome_id=pangenome.pangenome_id,
        species_key=species.species_key,
        previous_pangenome_id=species.default_pangenome_id,
    )
    report.checks_passed, report.failures = verify_pangenome_is_servable(session, pangenome)

    if report.failures and not force:
        raise PublishRefused(
            f"{species.species_key}: {len(report.failures)} check(s) failed and the pointer was NOT "
            f"moved:\n  " + "\n  ".join(report.failures)
        )

    session.execute(
        update(PathogenSpecies)
        .where(PathogenSpecies.pathogen_species_id == species.pathogen_species_id)
        .values(default_pangenome_id=pangenome.pangenome_id)
    )
    session.execute(
        update(Pangenome)
        .where(Pangenome.pangenome_id == pangenome.pangenome_id)
        .values(is_published=True)
    )
    # ⛔ **THE OUTGOING DEFAULT IS NOT UN-PUBLISHED, AND THAT IS THE WHOLE POINT OF 2026-09-25.**
    # This used to clear `is_published` on whatever the species pointer previously named. That was
    # harmless while the column had no reader and a species had exactly one catalogue: it retired a
    # superseded GENERATION of the same model. It is wrong now, and wrong in the one way that defeats
    # the feature it sits inside.
    #
    # `is_published` means *offer this catalogue in the picker*, and `default_pangenome_id` means
    # *serve this one when asked for the bare species*. They are different facts. The pointer's
    # previous value is now typically a DIFFERENT MODEL — nuna4 — not an older generation of this
    # one. Clearing it made `GET /api/v1/catalogues` drop nuna4 the moment nuna5 was published, so
    # publishing the second model silently deleted the comparator the picker exists to show, with
    # nothing in any output saying so. Reproduced end to end before this was changed.
    #
    # ⚠ **Publishing therefore adds, and never removes.** Hiding a catalogue is a separate, explicit
    # act — it is not a side effect of promoting another one. Rollback is unaffected: it was always
    # the pointer that made it one update, and the pointer still moves.
    report.published = True
    return report
