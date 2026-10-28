"""Every SQLAlchemy model, imported here so Alembic autogenerate sees the whole metadata.

⛔ A model not imported here is **invisible to autogenerate**, which produces a migration that
silently omits its table — and the failure surfaces much later, as a missing relation on one
endpoint. Adding a model means adding it to this file.
"""

from bacatlas_backend.models.gene import (  # noqa: F401
    Gene,
    GeneFunctionalAnnotation,
    GeneLocusMembership,
    GenomeNoncodingFeature,
)
from bacatlas_backend.models.genome import Genome, GenomeAssembly, GenomeContig  # noqa: F401
from bacatlas_backend.models.genome_collection import (  # noqa: F401
    GenomeCollection,
    GenomeCollectionBuildReport,
    GenomeCollectionMembership,
)
from bacatlas_backend.models.intergenic_gap import IntergenicGap, IntergenicGapFeature  # noqa: F401
from bacatlas_backend.models.locus import Locus  # noqa: F401
from bacatlas_backend.models.locus_annotation import (  # noqa: F401
    LocusAnnotationEntry,
    LocusUnirefFamilyCrosstab,
)
from bacatlas_backend.models.locus_arrangement import LocusArrangement  # noqa: F401
from bacatlas_backend.models.locus_offset_occupant import LocusOffsetOccupant  # noqa: F401
from bacatlas_backend.models.locus_similarity import (  # noqa: F401
    LocusNearestLocus,
    LocusSimilarity,
    PangenomeSimilarityBaseline,
)
from bacatlas_backend.models.pangenome import (  # noqa: F401
    Pangenome,
    PangenomeEvaluation,
    PangenomeInputEdge,
    PangenomeStep,
)
from bacatlas_backend.models.pangenome_genome_locus_count import (  # noqa: F401
    PangenomeGenomeLocusCount,
)
from bacatlas_backend.models.pangenome_propagation_rate import (  # noqa: F401
    PangenomePropagationRate,
)
from bacatlas_backend.models.pathogen_species import PathogenSpecies  # noqa: F401
from bacatlas_backend.models.projected_genome import (  # noqa: F401
    ProjectedGenePlacement,
    ProjectedGenome,
)
from bacatlas_backend.models.reference_vocabulary import PfamFamily  # noqa: F401

__all__ = [
    "Gene",
    "GeneFunctionalAnnotation",
    "GeneLocusMembership",
    "Genome",
    "GenomeAssembly",
    "GenomeCollection",
    "GenomeCollectionBuildReport",
    "GenomeCollectionMembership",
    "GenomeContig",
    "GenomeNoncodingFeature",
    "IntergenicGap",
    "IntergenicGapFeature",
    "Locus",
    "LocusAnnotationEntry",
    "LocusArrangement",
    "LocusSimilarity",
    "LocusNearestLocus",
    "PangenomeSimilarityBaseline",
    "LocusOffsetOccupant",
    "LocusUnirefFamilyCrosstab",
    "Pangenome",
    "PangenomeEvaluation",
    "PangenomeGenomeLocusCount",
    "PangenomeInputEdge",
    "PangenomeStep",
    "PathogenSpecies",
    "ProjectedGenePlacement",
    "ProjectedGenome",
]
