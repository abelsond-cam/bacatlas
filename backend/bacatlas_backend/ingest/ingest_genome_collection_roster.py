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

    ⚠ **`collection_label` defaults to the KEY, and that default is meant to be visible.** The label
    is what the picker shows now that the cohort is the browser's axis, and a key in a picker reads
    as obviously unfinished — where a blank or a species name would read as finished and wrong.
    `publish_pangenome` refuses to offer a catalogue whose cohort is still wearing its key, so the
    default can never reach a reader; it only survives on a catalogue that is loaded and not offered.
    """
    samples, n_genes = read_genome_vocabulary(artifacts)
    key = collection_key or artifacts.set_key

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
    collection.collection_label = collection_label or key
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
