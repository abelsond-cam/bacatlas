"""David's 2026-11-26 audit redesign, end to end on the backend: the refusal, the read, the lists, the rank.

⭐ The ruling: at 71,435 genomes the audit's per-locus verdicts were decided by stray genes — kp-71k locus
4008 (Wzi, PF14052 on 69,825 of 69,884 genes) read "grouped on context alone" and "conflicting
architecture". The redesigned audit states numbers, and this service must (1) refuse an audit that never
wrote them, (2) read them verbatim, (3) list the footer's residuals from them, and (4) bar no locus from
the showcase.

⛔ Every refusal here is made to FIRE, and every fixture keeps apart the values whose difference is under
test: a measured empty domain set against NULL, `False` against NULL, a singleton against a locus with
members, a locus with no Pfam against one whose architecture is not shared.
"""

from __future__ import annotations

import importlib.util
import json
import math

import pandas
import pytest

from bacatlas_backend.ingest import catalogue_frames
from bacatlas_backend.ingest.artifact_locator import CatalogueArtifacts
from bacatlas_backend.ingest.catalogue_frames import (
    PAGE_STATED_AUDIT_THRESHOLDS,
    REDESIGNED_CONCORDANCE_COLUMNS,
    REDESIGNED_EVIDENCE_COLUMNS,
    PreRedesignAudit,
    refuse_unless_redesigned_audit,
)
from bacatlas_backend.ingest.derive_locus_ranking import (
    BARRED_SENTINEL,
    SHOWCASE_POLICY,
    WEIGHT_DEPTH,
    InterestInputs,
    interest_score,
)
from bacatlas_backend.ingest.ingest_locus_catalogue import _clean
from bacatlas_backend.models.enumerations import PrevalenceBand
from bacatlas_backend.services.audit_residual_service import load_audit_residuals
from bacatlas_backend.services.species_catalogue_service import _load_catalogue_shell

from .conftest import (
    PUBLISHED_SITE_CATALOGUE_DIR,
    artifacts_for,
    make_locus,
    redesigned_audit_artifacts,
    requires_artifacts,
)

#: ⚠ The audit's evidence frame is read through nuna's `load_audit_evidence` — private, and absent in CI —
#: so every test that reaches `_load_audit_frames` says so and skips there rather than failing.
requires_nuna = pytest.mark.skipif(
    importlib.util.find_spec("nuna") is None, reason="nuna is not importable — the audit evidence is read through it"
)

EVIDENCE = ["node", "tier", "resolved_threshold", *REDESIGNED_EVIDENCE_COLUMNS]
CONCORDANCE = ["node", "class_clan", "n_arch", *REDESIGNED_CONCORDANCE_COLUMNS]


def _refuse(evidence=EVIDENCE, source="x_cluster_table.parquet", concordance=CONCORDANCE, summary=None):
    refuse_unless_redesigned_audit(
        evidence_columns=evidence,
        evidence_source=source,
        concordance_columns=concordance,
        concordance_name="x_pfam_concordance.tsv",
        summary=dict(PAGE_STATED_AUDIT_THRESHOLDS) if summary is None else summary,
        summary_name="x_audit_summary.json",
    )


# ── the refusal ────────────────────────────────────────────────────────────────────────────────
def test_a_redesigned_audit_passes():
    _refuse()


@pytest.mark.parametrize("column", list(REDESIGNED_EVIDENCE_COLUMNS))
def test_an_evidence_frame_missing_ANY_new_column_is_refused_by_name(column):
    with pytest.raises(PreRedesignAudit, match=rf"predates the 2026-11-26 redesign(.|\n)*lacks {column}"):
        _refuse(evidence=[c for c in EVIDENCE if c != column])


@pytest.mark.parametrize("column", list(REDESIGNED_CONCORDANCE_COLUMNS))
def test_a_concordance_missing_ANY_new_column_is_refused_by_name(column):
    with pytest.raises(PreRedesignAudit, match=rf"x_pfam_concordance.tsv lacks {column}"):
        _refuse(concordance=[c for c in CONCORDANCE if c != column])


def test_an_ABSENT_concordance_is_refused_rather_than_read_as_no_pfam_facts():
    with pytest.raises(PreRedesignAudit, match="x_pfam_concordance.tsv is absent, so dominant_domains"):
        _refuse(concordance=None)


@pytest.mark.parametrize("key", list(PAGE_STATED_AUDIT_THRESHOLDS))
def test_a_summary_missing_ANY_threshold_is_refused_by_name(key):
    summary = {k: v for k, v in PAGE_STATED_AUDIT_THRESHOLDS.items() if k != key}
    with pytest.raises(PreRedesignAudit, match=rf"x_audit_summary.json lacks {key}"):
        _refuse(summary=summary)


def test_EVERY_problem_is_named_in_one_refusal_not_the_first():
    with pytest.raises(PreRedesignAudit) as caught:
        _refuse(evidence=["node", "tier"], concordance=["node"], summary={})
    message = str(caught.value)
    for name in (*REDESIGNED_EVIDENCE_COLUMNS, *REDESIGNED_CONCORDANCE_COLUMNS, *PAGE_STATED_AUDIT_THRESHOLDS):
        assert name in message, name


def test_a_WATERFALL_source_names_the_nuna_filter_a_present_column_can_be_lost_in():
    with pytest.raises(PreRedesignAudit, match="EVIDENCE_COLS"):
        _refuse(evidence=["node", "tier"], source="x_homology_waterfall.csv")


@pytest.mark.parametrize(
    ("key", "value"), [("mmseqs_group_min_share", 0.9), ("pfam_shared_min_share", 0.8), ("esm_quoted_quantile", 0.1)]
)
def test_a_threshold_the_page_does_NOT_state_is_refused_naming_both_values(key, value):
    with pytest.raises(PreRedesignAudit, match=rf"{key} = {value}, where the page states"):
        _refuse(summary={**PAGE_STATED_AUDIT_THRESHOLDS, key: value})


def test_a_threshold_written_as_its_own_float_rounding_still_passes():
    _refuse(summary={**PAGE_STATED_AUDIT_THRESHOLDS, "pfam_shared_min_share": 0.95000000000001})


@pytest.mark.parametrize("value", [1, 0, "yes", "", "TRUE "])
def test_a_pfam_shared_cell_that_is_neither_true_nor_false_is_REFUSED_not_coerced(value):
    """⛔ `bool("False")` is True: a lenient reading would mark every locus shared."""
    with pytest.raises(catalogue_frames.CatalogueFrameError, match="neither true nor false"):
        catalogue_frames._audit_bool(value)


@pytest.mark.parametrize(("value", "expected"), [(True, True), ("False", False), ("true", True), (math.nan, None)])
def test_a_pfam_shared_cell_reads_as_what_it_says(value, expected):
    assert catalogue_frames._audit_bool(value) is expected


# ── the read, on a tiny audit tree ─────────────────────────────────────────────────────────────
def _audit_tree(tmp_path, *, redesigned: bool) -> CatalogueArtifacts:
    artifacts = CatalogueArtifacts(data_root=tmp_path, set_key="kp", model_label="kp_m", run_id="r")
    artifacts.cluster_table.parent.mkdir(parents=True)
    evidence = pandas.DataFrame(
        {"node": ["1", "2", "3"], "tier": ["mmseq@0.98", "pfam_shared", "esm_only"],
         "bucket": ["mmseq@0.98", "pfam_shared", "esm_only"], "resolved_threshold": [0.98, math.nan, math.nan]}
    )  # fmt: skip
    concordance = pandas.DataFrame({"node": ["1", "2", "3", "4"], "class_clan": ["single", "disjoint", "x", "x"],
                                    "n_arch": [1, 2, 3, 0]})  # fmt: skip
    summary = {"label": "kp_m"}
    if redesigned:
        evidence["resolved_share"] = [0.96, math.nan, math.nan]
        evidence["permissive_share"] = [1.0, 0.4, 0.2]
        evidence["esm_member_sim_p05"] = [0.99, 0.81, 0.7123]
        evidence["esm_member_sim_median"] = [0.995, 0.9, 0.88]
        concordance["dominant_domains"] = ["PF00126,PF03466", "PF14052", "PF00001", ""]
        concordance["n_dominant_domain_genes"] = [97, 69_825, 41, 0]
        concordance["n_other_arch_genes"] = [1, 3, 30, 0]
        concordance["dominant_domain_share"] = [0.97, 0.99916, 0.41, 0.0]
        concordance["pfam_shared"] = [True, True, False, False]
        summary.update(PAGE_STATED_AUDIT_THRESHOLDS)
    evidence.to_parquet(artifacts.cluster_table)
    concordance.to_csv(artifacts.pfam_concordance, sep="\t", index=False)
    artifacts.audit_summary.write_text(json.dumps(summary))
    return artifacts


@requires_nuna
def test_a_PRE_REDESIGN_audit_tree_is_refused_at_load(tmp_path):
    with pytest.raises(PreRedesignAudit, match="kp_m_cluster_table.parquet lacks resolved_share"):
        catalogue_frames._load_audit_frames(_audit_tree(tmp_path, redesigned=False))


@requires_nuna
def test_a_redesigned_audit_tree_reads_and_a_MEASURED_empty_domain_set_stays_a_string(tmp_path):
    """⛔ `""` (no gene annotated — measured) must not become NaN (not measured) on the way in."""
    audit = catalogue_frames._load_audit_frames(_audit_tree(tmp_path, redesigned=True))
    concordance = audit["concordance"]
    assert concordance.loc["4", "dominant_domains"] == ""
    assert concordance.loc["1", "dominant_domains"] == "PF00126,PF03466"
    assert [catalogue_frames._audit_bool(v) for v in concordance["pfam_shared"]] == [True, True, False, False]
    assert audit["evidence"].loc["3", "esm_member_sim_p05"] == pytest.approx(0.7123)


# ── the read, through the WHOLE build, on the real probe ───────────────────────────────────────
@requires_nuna
@requires_artifacts
def test_the_REAL_mirrored_probe_audit_is_refused_BEFORE_the_per_genome_read(monkeypatch):
    """⭐ The refusal on a real pre-redesign audit — and before the expensive half of the build, which at
    71,435 genomes is hours of per-genome reads that would otherwise be spent and discarded."""

    def must_not_be_reached(*args, **kwargs):
        raise AssertionError("the per-genome sources were read before the audit was refused")

    monkeypatch.setattr(catalogue_frames, "_load_source_frames", must_not_be_reached)
    with pytest.raises(PreRedesignAudit, match="predates the 2026-11-26 redesign") as caught:
        catalogue_frames.build_catalogue_frames(artifacts_for("kp"))
    for column in (*REDESIGNED_EVIDENCE_COLUMNS, *REDESIGNED_CONCORDANCE_COLUMNS, *PAGE_STATED_AUDIT_THRESHOLDS):
        assert column in str(caught.value), column


@requires_nuna
@requires_artifacts
@pytest.mark.parametrize("stage", ["pangenome", "all"])
def test_the_CLI_refuses_a_pre_redesign_audit_BEFORE_anything_is_written(monkeypatch, capsys, stage):
    """⛔ Exit 1 (a refusal), the reason on stderr — and before `all` reaches the genome layer, which at 71,435
    genomes is hours, or either stage opens the database at all."""
    from bacatlas_backend.ingest import ingest_command_line

    def reached(*args, **kwargs):
        raise AssertionError("the load went past the audit check")

    for name in ("create_engine", "load_reference_vocabularies", "load_genome_layer", "load_pangenome_layer"):
        monkeypatch.setattr(ingest_command_line, name, reached)
    real = artifacts_for("kp")
    code = ingest_command_line.main([
        "--database-url", "postgresql+psycopg://nobody@localhost:1/never_opened", "--stage", stage,
        "--data-root", str(real.data_root), "--set-key", real.set_key, "--species-key", real.species_key,
        "--model-label", real.model_label, "--run-id", real.run_id,
    ])  # fmt: skip
    assert code == 1
    err = capsys.readouterr().err
    assert "AUDIT REFUSED: this catalogue's audit predates the 2026-11-26 redesign" in err
    assert "lacks resolved_share" in err


@requires_nuna
@requires_artifacts
def test_the_whole_build_READS_each_new_column_into_its_locus_column(tmp_path):
    """⛔ Three real loci given three different sets of values, so no two columns can be swapped unseen.

    The copy's placeholders are overwritten for nodes chosen from the real catalogue — one MMseqs-grouped,
    two not — and everything else keeps its placeholder, which is itself checked to survive as written.
    """
    artifacts = redesigned_audit_artifacts("kp", tmp_path / "root")
    evidence = pandas.read_parquet(artifacts.cluster_table)
    evidence["node"] = evidence["node"].astype(str)
    resolved = evidence.loc[evidence["resolved_threshold"].notna(), "node"].iloc[0]
    unresolved = evidence.loc[evidence["resolved_threshold"].isna(), "node"].iloc[:2].tolist()
    chosen = [resolved, *unresolved]
    values = {
        chosen[0]: dict(resolved_share=0.96, esm_member_sim_p05=0.991, esm_member_sim_median=0.997),
        chosen[1]: dict(resolved_share=math.nan, esm_member_sim_p05=0.812, esm_member_sim_median=0.93),
        chosen[2]: dict(resolved_share=math.nan, esm_member_sim_p05=math.nan, esm_member_sim_median=math.nan),
    }
    for node, row in values.items():
        for column, value in row.items():
            evidence.loc[evidence["node"] == node, column] = value
    evidence.to_parquet(artifacts.cluster_table)
    concordance = pandas.read_csv(artifacts.pfam_concordance, sep="\t", dtype={"node": str},
                                  converters={"dominant_domains": lambda v: v})  # fmt: skip
    pfam = {
        chosen[0]: ("PF00126,PF03466", 5, 1, True),
        chosen[1]: ("PF14052", 7, 2, False),
        chosen[2]: ("", 0, 0, False),
    }
    for node, (domains, dominant, other, shared) in pfam.items():
        mask = concordance["node"] == node
        concordance.loc[mask, ["dominant_domains", "n_dominant_domain_genes", "n_other_arch_genes"]] = [
            domains, dominant, other
        ]
        concordance.loc[mask, "pfam_shared"] = shared
    concordance.to_csv(artifacts.pfam_concordance, sep="\t", index=False)

    loci = catalogue_frames.build_catalogue_frames(artifacts).loci.set_index("node")
    columns = ["mmseqs_group_share", "esm_member_similarity_p05", "esm_member_similarity_median",
               "pfam_dominant_domains", "pfam_dominant_domain_gene_count", "pfam_other_architecture_gene_count",
               "pfam_architecture_shared"]  # fmt: skip
    # ⚠ Through the row builder's own `_clean`, which is what turns a frame's NaN into the NULL it stores.
    got = {node: tuple(_clean(value) for value in loci.loc[node, columns]) for node in chosen}
    assert got == {
        chosen[0]: (0.96, 0.991, 0.997, "PF00126,PF03466", 5, 1, True),
        chosen[1]: (None, 0.812, 0.93, "PF14052", 7, 2, False),
        chosen[2]: (None, None, None, "", 0, 0, False),
    }
    # Coverage, stated: every locus carries the columns, and the untouched ones their placeholder.
    assert len(loci) > 1000 and set(columns) <= set(loci.columns)
    untouched = loci.drop(index=chosen)
    assert (untouched["pfam_dominant_domains"] == "").all()
    assert (untouched["pfam_dominant_domain_gene_count"] == 0).all()


# ── the footer's two lists and the shell's counts of them ──────────────────────────────────────
def _residual_fixture(session, seeded):
    """Eight loci, each separating one condition of the two lists from its neighbour."""
    rows = {
        # label: (tier, pfam_shared, member_gene_count, pfam_annotated_member_count)
        "esm+notshared": ("esm_only", False, 40, 30),      # both lists
        "esm+nopfam": ("esm_only", False, 40, 0),          # list 1 only: no gene carries a domain
        "mmseq+notshared": ("mmseq@0.98", False, 40, 30),  # list 2 only
        "shared": ("pfam_shared", True, 40, 39),           # neither: shared
        "singleton": (None, False, 1, 1),                  # neither: nothing to share
        "unmeasured": ("synteny_only", None, 40, 30),      # neither: audited before 2026-11-26
        "align": ("alignable_not_grouped", True, 40, 40),  # neither
        "esm+unmeasured": ("esm_only", None, 40, 30),      # list 1 only: pfam_shared not measured
    }
    for ordinal, (label, (tier, shared, size, annotated)) in enumerate(rows.items()):
        session.add(
            make_locus(
                seeded, ordinal=400 + ordinal, label=label, collapse_tier=tier, pfam_architecture_shared=shared,
                member_gene_count=size, member_genome_count=min(size, 97), pfam_annotated_member_count=annotated,
                pfam_dominant_domain_gene_count=0 if annotated is not None else None,
                esm_member_similarity_p05=0.8 if tier == "esm_only" else None,
                prevalence_band=PrevalenceBand.RARE if size == 1 else PrevalenceBand.CORE,
            )
        )  # fmt: skip
    session.flush()


def test_the_two_lists_name_exactly_the_loci_their_definitions_name(session, seeded):
    _residual_fixture(session, seeded)
    residuals = load_audit_residuals(session, pangenome_id=seeded["pangenome"].pangenome_id)
    assert [row.label for row in residuals.esm_only] == ["esm+notshared", "esm+nopfam", "esm+unmeasured"]
    assert [row.label for row in residuals.pfam_architecture_not_shared] == ["esm+notshared", "mmseq+notshared"]
    assert residuals.esm_only[0].esm_member_similarity_p05 == pytest.approx(0.8)


def test_the_shell_counts_each_list_from_the_SAME_clauses(session, seeded):
    """⛔ The summary a closed list shows must be the length of the list it opens to."""
    _residual_fixture(session, seeded)
    catalogue = _load_catalogue_shell(session, seeded["species"], seeded["pangenome"])
    residuals = load_audit_residuals(session, pangenome_id=seeded["pangenome"].pangenome_id)
    assert len(residuals.esm_only) == 3 and len(residuals.pfam_architecture_not_shared) == 2
    assert catalogue.residual_locus_counts == {"esm_only": 3, "pfam_architecture_not_shared": 2}


# ── the showcase bars nothing ──────────────────────────────────────────────────────────────────
def _row(tier, pfclass="single", resolved=None):
    return InterestInputs(
        syntenic_a5=0.9, uniref50_major_family_count=2, collapse_tier=tier, resolved_threshold=resolved,
        prevalence_band_index=2, member_gene_count=40, concordance=0.5, pfam_concordance_class=pfclass,
    )  # fmt: skip


@pytest.mark.parametrize("tier", ["esm_only", "synteny_only", "no_homology", "pfam_shared", "mmseq@0.4"])
def test_SHOWCASE_POLICY_bars_no_locus_on_its_pfam_class_or_its_tier(tier):
    """⛔ The two old bars — a `disjoint` clan class, a "failure" tier — are exactly what was ruled a bug."""
    assert interest_score(_row(tier, "disjoint", 0.4 if tier.startswith("mmseq") else None), SHOWCASE_POLICY) > (
        BARRED_SENTINEL
    )


def test_pfam_shared_takes_the_full_rescue_depth_and_esm_only_takes_none():
    shared = interest_score(_row("pfam_shared"), SHOWCASE_POLICY)
    esm_only = interest_score(_row("esm_only"), SHOWCASE_POLICY)
    assert shared - esm_only == pytest.approx(WEIGHT_DEPTH)


def test_a_policy_must_SAY_what_it_bars_and_a_missing_key_is_refused():
    """⛔ The old `.get("contested_pfclass", "disjoint")` barred by omission; now a missing key raises."""
    with pytest.raises(KeyError):
        interest_score(_row("esm_only"), {"failure_tiers": (), "rescue_tiers": ()})


def test_the_rebuild_s_PINNED_policy_is_the_frozen_pages_own():
    """⭐ The payload rebuild reproduces the frozen pages, so its policy is theirs — read from both pages."""
    from bacatlas_backend.instruments.payload_serialiser import INGESTED_PAYLOAD_POLICY

    for species in ("ecoli", "kp"):
        path = PUBLISHED_SITE_CATALOGUE_DIR / f"{species}.json"
        if not path.exists():
            pytest.skip(f"the published site catalogue is not present at {path}")
        assert json.loads(path.read_text())["meta"]["policy"] == INGESTED_PAYLOAD_POLICY, species
