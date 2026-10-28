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


# ⛔⛔ **`ExclusivityForm` and `ExclusivityFormSource` ARE DELETED (2026-10-27), with the four
# columns that used them: `pangenome.exclusivity_form` / `_source`, `nuna_model.exclusivity_form`
# and `nuna_model_step.provenance_source`. Both PG types are dropped by `b0c7d4ea9f31`.**
#
# ⚠ The VALUES `exclusion` and `damped_exclusion` are on the never-rename list
# (`nuna/PROJECT_STATE.md`), and this is not a rename: nothing in this service reads a run_id
# token any more, so there is no decoded form to spell one way or another. The tokens themselves
# are untouched — they are nuna's, they are still emitted, and the run_id that carries them is
# stored here verbatim. What is gone is this repo's second, derived copy of that fact.
#
# ⚠ The enum's own docstring was ALSO out of date before it was deleted: it named
# `nuna.tl.cluster.exclusivity.form_from_run_id` as its populator, and nuna renamed that module to
# `hypergeometric` on 2026-10-23.


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


class PropagationVerdict(enum.Enum):
    """Whether a node's own annotated genes agree — the five states, one per locus per vocabulary.

    ⭐ **This is `within_node_propagation`'s decision, named and stored.** It is not a new question:
    `compute_propagation` has always sorted every node into exactly these five cases and then thrown
    four of the labels away, keeping only the per-band tally. Storing the label per locus makes the
    tally a `GROUP BY` over a column instead of a scan over 1.3 M gene rows.

    ⛔ **Five values and not three, because two of them are different kinds of silence.**
    ``no_coverage`` is *no member gene carries this vocabulary at all*; ``not_comparable`` is *two or
    more do and no rung exists that every one of them states*. Both are counted in NEITHER the
    `checkable` nor the `one_gene` column — and collapsing them would lose the only evidence that
    distinguishes a vocabulary this node is silent about from one it speaks incomparably. That is the
    same distinction `no_homology` was retired for having blurred.

    ⚠ ``one_gene`` is not a defect and not a failure: it is the state the whole propagation argument
    is *about* — one annotated gene covering the node's other members. It is unanimous by
    construction, which is exactly why it is kept out of `unanimous`.
    """

    NO_COVERAGE = "no_coverage"
    ONE_GENE = "one_gene"
    NOT_COMPARABLE = "not_comparable"
    UNANIMOUS = "unanimous"
    DISCORDANT = "discordant"
