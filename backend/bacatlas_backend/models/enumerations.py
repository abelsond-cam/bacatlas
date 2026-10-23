"""Every enumerated vocabulary in the schema, in one place.

⛔ **Three of these are enums specifically because a nullable boolean or a bare string would let a
real distinction be lost silently.** Each carries the failure it prevents.
"""

from __future__ import annotations

import enum


class GateVerdict(enum.Enum):
    """pass / fail / **incomplete** — and `incomplete` is a VALUE, never a NULL.

    ⛔ `nuna.eval.gates` returns ``passed=None`` for a metric a run could not produce, and *"an
    unmeasured gate must not masquerade as a pass or a fail"*. Modelled as a nullable boolean,
    ``bool_and`` over SQL NULLs silently converts *incomplete* into *not counted*, which is exactly
    the failure that module exists to prevent.
    """

    PASS = "pass"
    FAIL = "fail"
    INCOMPLETE = "incomplete"


class ExclusivityForm(enum.Enum):
    """The genome-exclusivity edge weight a step used.

    ⛔ The on-disk tokens are FROZEN and ``-excl`` is a **prefix** of ``-exclLOGP``. Four readers
    got this wrong silently in one day. This column is populated by
    ``nuna.tl.cluster.exclusivity.form_from_run_id`` and by nothing else; after ingest no query ever
    looks at a run_id string again.
    """

    #: ``f = min(1, s·−log₁₀P₀)`` — the standard form. Token ``-exclLOGP``.
    EXCLUSION = "exclusion"
    #: ``h = s·(1−P₀)`` — superseded; it SATURATES at N=100 for a=b≥40. Token ``-excl``.
    DAMPED_EXCLUSION = "damped_exclusion"
    #: The step does not use an exclusivity weight at all.
    NONE = "none"


class ExclusivityFormSource(enum.Enum):
    """How the exclusivity form was determined — so a reader knows how much to trust it."""

    RUN_MANIFEST = "run_manifest"
    RUN_ID_TOKEN = "run_id_token"
    MODEL_REGISTRY = "model_registry"


class PrevalenceBand(enum.Enum):
    """``nuna.eval.prevalence.categorise``'s vocabulary, unchanged.

    core ≥ 0.99N · soft_core ≥ 0.95N · shell ≥ 0.15N · cloud ≥ 2 · rare = singleton.
    """

    CORE = "core"
    SOFT_CORE = "soft_core"
    SHELL = "shell"
    CLOUD = "cloud"
    RARE = "rare"


class RhoRule(enum.Enum):
    """The ρ rail a step ran under. **Per step, never per model.**"""

    OFF = "off"
    CEILING = "ceiling"
    PAIRWISE_MAX = "pairwise_max"


class EmbeddingRepresentation(enum.Enum):
    """The two representations. ESM is sequence-only; Bacformer is a transformer over a genome's."""

    ESM = "esm"
    BACFORMER = "bacformer"


class RosterLineage(enum.Enum):
    """Which of the two incompatible roster constructions a collection came from.

    They do not share a schema and never have: the probe lineage is BioSample-keyed with AMR
    phenotype records behind it, the Klebsiella lineage is `metadata_v2`-derived with sublineage,
    LIN code and assembly-quality flags. Naming the lineage is what stops a loader assuming either.
    """

    PROBE_BIOSAMPLE = "probe_biosample"
    KLEBSIELLA_METADATA_V2 = "klebsiella_metadata_v2"


class SampleIdentifierKind(enum.Enum):
    """What `genome.sample_id` actually is, because it is not one thing.

    ⚠ BioSample-keyed short-read assemblies (`SAMEA103923484`) have ONE form. RefSeq/GenBank
    assembly stems have TWO — `GCF_000512165.1` and `GCF_000512165.1_ASM51216v1_genomic` — because
    NCBI writes the submitter's strain name in the tail (`_INF156_genomic`, `_KSB2_2B_genomic`, two
    underscores). An `_ASM`-only matching rule silently lost 27 genomes.
    """

    BIOSAMPLE = "biosample"
    ASSEMBLY_STEM = "assembly_stem"


class GeneOntologyAgreementVerdict(enum.Enum):
    """The Pfam ladder plus the not-judgeable case, so one verdict vocabulary serves both.

    ⛔ ``no_coverage`` is a **value**, not a NULL: fewer than two annotated members is neither
    agreement nor disagreement, and must never be counted as either.
    """

    NO_COVERAGE = "no_coverage"
    SINGLE = "single"
    SAME_DOMAINS = "same_domains"
    NESTED = "nested"
    OVERLAPPING = "overlapping"
    DISJOINT = "disjoint"


#: ⛔ **The index↔name contract for `locus_annotation_entry.gene_ontology_namespace`**, which is
#: stored as 0/1/2 because one list serves all three namespaces.
#:
#: ⚠ **Verified against the DATA, not against a comment** (2026-09-14, the dev database — then
#: the ingest database): namespace 0
#: is dominated by *nucleic acid binding*, *metal ion binding*, *ion binding*; 1 by *transport*,
#: *metabolic process*, *biosynthetic process*; 2 by *cytoplasm*, *plasma membrane*, *membrane*.
#: Each namespace also carries its own root term as a slim class — literally `molecular_function`,
#: `biological_process`, `cellular_component` — so the ordering identifies itself. A transposed
#: table here would file every molecular function under cellular component and the page would look
#: entirely normal.
GENE_ONTOLOGY_NAMESPACE_NAMES: tuple[str, str, str] = (
    "molecular_function",
    "biological_process",
    "cellular_component",
)


def gene_ontology_namespace_name(index: int | None) -> str | None:
    """0/1/2 → the namespace name the API speaks, or `None` for a non-GO row.

    ⛔ The API must never emit the index. `coverage.go_annotated_gene_count` is already keyed by
    NAME, so an entry keyed by number makes one response speak two dialects — and a client that
    groups entries by a key the coverage block does not use renders three GO cards whose term lists
    are all empty while their coverage lines say otherwise.
    """
    if index is None:
        return None
    try:
        return GENE_ONTOLOGY_NAMESPACE_NAMES[index]
    except IndexError:  # pragma: no cover - a fourth namespace would be a data error
        raise ValueError(f"no GO namespace with index {index}") from None


class EvaluationKind(enum.Enum):
    """Which arm of the funnel an evaluation row belongs to.

    ``merge_evidence`` is reserved: it is another agent's live strand (`docs/merge_evidence.md`),
    and it is named here so it lands as a third kind rather than a fourth table.
    """

    COPY_STRUCTURE_GATE = "copy_structure_gate"
    ACCESSORY_AUDIT = "accessory_audit"
    MERGE_EVIDENCE = "merge_evidence"


class AnnotationKind(enum.Enum):
    """Which vocabulary a `locus_annotation_entry` row draws its term from."""

    GENE_SYMBOL = "gene_symbol"
    PROTEIN_PRODUCT = "protein_product"
    PFAM_ARCHITECTURE = "pfam_architecture"
    COG_ORTHOGROUP = "cog_orthogroup"
    GENE_ONTOLOGY_SLIM = "gene_ontology_slim"
    EC_NUMBER = "ec_number"
    KEGG_ORTHOLOGY = "kegg_orthology"
