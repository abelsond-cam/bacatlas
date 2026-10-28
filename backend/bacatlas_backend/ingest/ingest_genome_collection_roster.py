"""The roster a pangenome was built over — `genome_collection` and its ordinals.

⭐ **The ordinal is the whole point of this table.** `meta.genomes` is the genome VOCABULARY that
`arr.gid` indexes into, and it is `sorted(universe["sample_id"].unique())` — `export_payload.main`'s
own line. Every arrangement's membership is a list of positions in that list, so getting the order
wrong renames every genome on the page and nothing downstream can detect it.

⛔ **Stored, never re-derived by sorting at read time.** It is a contract with payloads that are
already published, and re-deriving an ordering that was once written down is exactly the implicit
ordering this project punishes. A locale-sensitive sort, a different pandas version's `unique()`
order, or a roster that later grows by one genome would each silently renumber it.

⚠ **N comes from the roster, not from the rows present.** `genome_count` is the prevalence
denominator *and* the exclusivity denominator, and `profile_from_pairs`' rule is that a genome
contributing no gene to a locus still counts against it. So it is `len(samples)` — which
`build_payload` itself refuses to let disagree with `n_genomes`.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

from sqlalchemy import select
from sqlalchemy.orm import Session

from bacatlas_backend.ingest.artifact_locator import CatalogueArtifacts
from bacatlas_backend.ingest.published_cohorts import cohort_for
from bacatlas_backend.models.enumerations import RosterLineage
from bacatlas_backend.models.genome import Genome
from bacatlas_backend.models.genome_collection import GenomeCollection, GenomeCollectionMembership
from bacatlas_backend.models.pangenome import Pangenome


class RosterError(RuntimeError):
    """A roster that cannot be turned into a collection without inventing something."""


@dataclass
class RosterReport:
    """What the roster load did, in numbers that can be reconciled against the store."""

    collection_key: str
    #: ⭐ What the picker will SHOW for this cohort — reported so the operator sees the words a
    #: reader will see, not just the machine key they were resolved from.
    collection_label: str
    genome_count: int
    members_written: int
    genes_in_universe: int
    #: ⛔ Named, not counted. A roster genome with no `genome` row means the genome layer is short,
    #: and every arrangement referring to it would silently address the wrong ordinal.
    unresolved_sample_ids: list[str]


def read_genome_vocabulary(artifacts: CatalogueArtifacts) -> tuple[list[str], int]:
    """`(samples, n_genes)` from the gene universe — `export_payload.main`'s two lines, verbatim.

    ⚠ The universe, not the assignment: *"a genome with no assigned gene still owns an index"*. The
    two differ exactly where a genome contributed nothing, and that is the case the ordinal must
    survive.
    """
    import pandas  # the `ingest` extra, deliberately not a serving dependency

    universe = pandas.read_parquet(artifacts.gene_universe, columns=["sample_id"])
    samples = sorted(str(value) for value in universe["sample_id"].unique())
    return samples, int(len(universe))


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def ingest_genome_collection(
    session: Session,
    artifacts: CatalogueArtifacts,
    *,
    pathogen_species_id: int,
    collection_key: str | None = None,
    collection_label: str | None = None,
    roster_lineage: RosterLineage = RosterLineage.PROBE_BIOSAMPLE,
) -> tuple[int, RosterReport]:
    """Ensure the collection and its ordered membership exist; return `(id, report)`.

    Idempotent by `(species, collection_key)`, and the membership is **replaced wholesale** rather
    than merged: an ordinal is a position in a list, so a partial update could leave two genomes
    claiming one position while every individual row still looked valid.

    ⛔⛔ **BOTH THE KEY AND THE LABEL DEFAULT FROM `published_cohorts.py`, NOT FROM THE ARTIFACTS,
    and that is what makes a chosen label survive the next load.** Two defaults were wrong in the
    same direction:

    * `collection_key` defaulted to `artifacts.set_key`, so a no-flags load wrote the catalogue key
      `{species}-{set_key}` — `ecoli-ecoli`, `kp-kp` — which passes all three of
      `catalogue_key_for`'s guards and is silently not a cohort anybody named. The checked-in
      recipe in `backend/README.md` passes neither flag, so following it produced exactly that.
    * `collection_label` defaulted to the KEY, and this function **assigns it unconditionally on an
      existing row**. So a label a human chose survived only until the next `--stage pangenome` run
      that omitted `--collection-label`, which put the machine token back; and `offer_catalogue`'s
      guard fires only on the way IN, so an already-offered catalogue never re-checks and the token
      reaches a reader with nothing refusing it. A write-once-per-invocation field defaulting to a
      machine token, on the axis the whole picker is built around.

    ⚠ An explicit `collection_key` / `collection_label` still wins, which is how a roster is staged
    by hand without being declared. An undeclared load with NEITHER is refused rather than invented
    — `UndeclaredCohort` names what is known.
    """
    samples, n_genes = read_genome_vocabulary(artifacts)
    declared = None
    if collection_key is None or collection_label is None:
        # ⚠ Looked up even when only ONE of the two was passed: the other still needs its default,
        # and taking it from the artifacts is what this change exists to stop.
        declared = cohort_for(artifacts.species_key, artifacts.set_key)
    key = collection_key or declared.collection_key
    label = collection_label or declared.label

    # ⛔ **THE COLLECTION KEY MUST BE ABLE TO BECOME A CATALOGUE KEY, and it is now half of one.**
    # `catalogue_key_for(species, cohort)` builds `{species}-{cohort}`, and its guards exist
    # because nuna's payload lookup globs `locus_browser_<key>_*.json` and takes the newest match —
    # so an underscore anywhere in a key lets a shorter key capture a longer one's payloads and
    # serve another catalogue's data in silence. Checked HERE, where it fails before the collection
    # row and its whole membership are written, rather than twenty lines later at the pangenome.
    # ⚠ Called for its guards, not its value: one implementation of the rule, not two.
    from bacatlas_backend.ingest.ingest_pangenome_run import catalogue_key_for

    catalogue_key_for(artifacts.species_key, key)

    collection = session.execute(
        select(GenomeCollection).where(
            GenomeCollection.pathogen_species_id == pathogen_species_id,
            GenomeCollection.collection_key == key,
        )
    ).scalar_one_or_none()
    if collection is None:
        collection = GenomeCollection(pathogen_species_id=pathogen_species_id, collection_key=key)
        session.add(collection)
    collection.roster_lineage = roster_lineage
    collection.collection_label = label
    collection.genome_count = len(samples)
    collection.roster_source_path = str(artifacts.gene_universe)
    collection.roster_source_sha256 = _file_sha256(artifacts.gene_universe)
    session.flush()

    genome_ids = {
        sample_id: genome_id
        for sample_id, genome_id in session.execute(
            select(Genome.sample_id, Genome.genome_id).where(Genome.sample_id.in_(samples))
        ).all()
    }
    unresolved = [sample_id for sample_id in samples if sample_id not in genome_ids]
    if unresolved:
        raise RosterError(
            f"{len(unresolved)} of {len(samples)} roster genomes have no `genome` row "
            f"({unresolved[:5]}). The ordinal is a position in this exact list, so loading a "
            "collection with holes in it would address every arrangement's membership one genome "
            "off. Load the genome layer for this species first."
        )

    # ⛔ **THE ONE ROW TWO MODELS SHARE.** `genome_collection` is unique on
    # `(pathogen_species_id, collection_key)`, so nuna4 and nuna5 over the same 100 genomes share a
    # single collection — and the membership below is deleted and rewritten on EVERY pangenome
    # ingest, with `collection_genome_ordinal = enumerate(samples)`. That ordinal is what
    # `meta.genomes` is ordered by, and `arr.gid` indexes into `meta.genomes`. So if a second
    # model's roster lists the same genomes in a different ORDER, this rewrite silently renames
    # every genome on every arrangement of the model that was already loaded — and every name it
    # then shows is a real genome, so nothing downstream can detect it. The existing parity suite
    # says exactly this: "one element out of place renames every genome on every arrangement."
    #
    # The check above already refuses HOLES. This refuses REORDERING, which only became reachable
    # when a species was allowed to hold more than one catalogue (2026-09-25).
    established = [
        sample_id
        for (sample_id,) in session.execute(
            select(GenomeCollectionMembership.requested_sample_id)
            .where(GenomeCollectionMembership.genome_collection_id == collection.genome_collection_id)
            .order_by(GenomeCollectionMembership.collection_genome_ordinal)
        ).all()
    ]
    if established and established != samples:
        dependants = session.execute(
            select(Pangenome.catalogue_key)
            .where(Pangenome.genome_collection_id == collection.genome_collection_id)
            .order_by(Pangenome.catalogue_key)
        ).scalars().all()
        if dependants:
            first = next(
                (i for i, (a, b) in enumerate(zip(established, samples, strict=False)) if a != b),
                min(len(established), len(samples)),
            )
            raise RosterError(
                f"collection {key!r} already holds a DIFFERENT roster, and "
                f"{len(dependants)} catalogue(s) are built on the order it has: "
                f"{', '.join(dependants)}.\n"
                f"  first difference at ordinal {first}: "
                f"{established[first] if first < len(established) else '(end)'} → "
                f"{samples[first] if first < len(samples) else '(end)'}\n"
                f"  ({len(established)} genomes stored, {len(samples)} offered)\n"
                "  An ordinal is a position in this exact list, and `arr.gid` indexes into it, so "
                "rewriting it would rename every genome on every arrangement of those catalogues "
                "while leaving each name a real genome. Load this model under its own "
                "--collection-key, or reconcile the rosters first."
            )

    session.query(GenomeCollectionMembership).filter(
        GenomeCollectionMembership.genome_collection_id == collection.genome_collection_id
    ).delete(synchronize_session=False)
    for ordinal, sample_id in enumerate(samples):
        session.add(
            GenomeCollectionMembership(
                genome_collection_id=collection.genome_collection_id,
                genome_id=genome_ids[sample_id],
                collection_genome_ordinal=ordinal,
                requested_sample_id=sample_id,
                matched_by="sample_id",
            )
        )
    session.flush()
    return collection.genome_collection_id, RosterReport(
        collection_key=key,
        collection_label=label,
        genome_count=len(samples),
        members_written=len(samples),
        genes_in_universe=n_genes,
        unresolved_sample_ids=unresolved,
    )


def genome_id_by_ordinal(session: Session, genome_collection_id: int) -> list[int]:
    """`genome_id` at each `collection_genome_ordinal`, dense — what `arr.gid` resolves through.

    ⛔ Returned as a list indexed by ordinal, and it **raises on a hole**: `arr.gid` is a positional
    index, so a gap would shift every membership after it by one and produce a page that names the
    wrong genomes with no error anywhere.
    """
    rows = session.execute(
        select(
            GenomeCollectionMembership.collection_genome_ordinal,
            GenomeCollectionMembership.genome_id,
        )
        .where(GenomeCollectionMembership.genome_collection_id == genome_collection_id)
        .order_by(GenomeCollectionMembership.collection_genome_ordinal)
    ).all()
    ordinals = [ordinal for ordinal, _ in rows]
    if ordinals != list(range(len(ordinals))):
        missing = sorted(set(range(len(ordinals))) - set(ordinals))
        raise RosterError(
            f"collection {genome_collection_id} has {len(ordinals)} members but their ordinals are "
            f"not 0..{len(ordinals) - 1} (missing {missing[:5]}). Every arrangement membership is a "
            "position in this list."
        )
    return [genome_id for _, genome_id in rows]
