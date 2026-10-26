"""Parity suite T6 — the Sequence tab, the frozen page against the API, gene by gene.

⛔ **What T6 exists to catch: a flank taken from the wrong end.** On a minus-strand gene — about half
of them — the upstream flank sits at HIGHER contig coordinates and is reverse-complemented with the
gene. Slicing `start - flank` unconditionally returns the downstream flank instead and renders as a
perfectly plausible 100 bases. So every flank is compared as a string, with its contig span and its
reverse-complement mark, and a swap is NAMED when it happens.

**The oracle is the page's own code** — `frozen_page_sequence_recorder.js` runs the published
`app.js`'s `seqGene` over the committed `.nseq` + `.loci` files and reads back what the tab showed;
the whole page is also booted and walked over a sample, which proves the lifted `seqGene` draws what
the page draws. **The other side is the endpoint's own path** — `load_gene_sequences` →
`serialise_gene_sequence` — slicing the original gzipped GFF. Both are recorded by
`sequence_parity_comparison`, which reads the page's DISPLAY back into values rather than re-rendering
the API's into sentences.

**Two modes.**

- *Default* (~50 s, both species): every genome of both published catalogues is touched — every gene's placement
  (coordinates, strand, contig, locus) compared from the page's own decoded tables; the rendered tab
  for up to four whole genomes per species chosen to stress it (most ρ > 1 loci, most contig-edge
  genes, most contigs whose name is not their index, the fewest contigs — three distinct on both
  species today), plus a stratified ~50 loci of every other genome with minus-strand, contig-edge
  and exact-truncation-boundary genes over-represented, and two near-core loci it LACKS.
- ``BACATLAS_SEQUENCE_PARITY_EVERY_GENE=1``: the rendered tab for EVERY gene of both catalogues,
  **1,021,997**. ⛔ **Measured 2026-10-26 on 8 workers: 3,335 s for ecoli and 4,367 s for kp —
  7,713 s, 2 h 09 m**, and ~3 h 28 m of wall clock including setup. This docstring claimed "~3.5
  min" until that was timed; budget the real number, because it is the gate the gene-layer redesign
  runs at its stages 3 and 4 and planning two of them as a coffee break is how a stage stalls. The
  cost is ~500 k ``api_genes_at`` calls per species, one per (genome, locus) pair — node is spawned
  only twice per genome, so the subprocesses are not the bottleneck.

⛔ **Coverage is asserted BEFORE any difference is reported**, by every test, through one helper —
genes compared of the database's own count, minus-strand genes, ρ > 1 copies, truncated flanks,
genes on contigs whose name is not their index. `4ab35ca`: a comparison that skips what it cannot
compare reports "0 differ" while skipping exactly what changed.

⚠ Skips, with the reason, when node, the nuna checkout (`app.js`, `dom_shim.js`), the committed
`.nseq` files, `BACATLAS_ROOT_GFF` or the loaded database is absent. Never passes vacuously.
"""

from __future__ import annotations

import json
import multiprocessing
import os
import shutil
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

import pytest
from sqlalchemy import create_engine, select, text
from sqlalchemy.orm import Session

from bacatlas_backend.models.locus import Locus
from bacatlas_backend.models.pangenome import Pangenome
from bacatlas_backend.models.pathogen_species import PathogenSpecies
from bacatlas_backend.services.gene_sequence_service import FLANK_LENGTH
from tests.conftest import PUBLISHED_SITE_CATALOGUE_DIR
from tests.known_parity_exceptions import CONTIG_ROW_SHOWS_NAME_NOT_SEQID
from tests.sequence_parity_comparison import (
    GENE_TABLE_FIELDS,
    RENDERED_FIELDS,
    GenomeReport,
    GenomeTask,
    absent_labels,
    compare_genome,
    compare_rendered,
    genome_gene_rows,
    run_recorder,
    sample_labels,
)

SPECIES_KEYS = ("ecoli", "kp")

#: ⭐ Named for what it does: render and compare the Sequence tab for every gene, not a sample.
EVERY_GENE = os.environ.get("BACATLAS_SEQUENCE_PARITY_EVERY_GENE", "").lower() in {"1", "true", "yes"}
WORKERS = int(os.environ.get("BACATLAS_SEQUENCE_PARITY_WORKERS", min(8, os.cpu_count() or 1)))

FROZEN_APP_JS = Path(
    os.environ.get("BACATLAS_FROZEN_APP_JS", Path.home() / "developer/nuna/src/nuna/tl/locus_browser/app.js")
)
FROZEN_DOM_SHIM = Path(os.environ.get("BACATLAS_FROZEN_DOM_SHIM", Path.home() / "developer/nuna/tests/js/dom_shim.js"))
SEQUENCE_DIR = PUBLISHED_SITE_CATALOGUE_DIR / "seq"
RECORDER = Path(__file__).with_name("frozen_page_sequence_recorder.js")

#: ⛔ Every difference name a report can carry, by the test that owns it. A name outside every group
#: fails `test_T6_every_difference_the_comparison_can_find_is_owned_by_exactly_one_test`, so no kind
#: of difference can be found and then reported by nobody.
FIELD_GROUPS = {
    "placement": tuple(f"table.{name}" for name in GENE_TABLE_FIELDS),
    "flanks": (
        "upstream_flank_sequence",
        "upstream_flank_span",
        "upstream_flank_is_reverse_complemented",
        "downstream_flank_sequence",
        "downstream_flank_span",
        "downstream_flank_is_reverse_complemented",
        "flanks_swapped_end_for_end",
    ),
    "coding_and_protein": (
        "coding_sequence",
        "coding_sequence_span",
        "coding_sequence_is_reverse_complemented",
        "coding_sequence_heading_length",
        "protein_sequence",
        "protein_heading_length",
        "length_stat",
        "gc_percent",
        "gc_percent_shown",
        "is_five_prime_partial",
    ),
    "copies": ("copies_at_locus", "copy_ordinal", "copy_count", "flat_index"),
    "truncation": (
        "upstream_flank_is_truncated_by_contig_end",
        "downstream_flank_is_truncated_by_contig_end",
        "edge_contig",
    ),
    "contig_and_strand": ("seqid", "start_position", "end_position", "strand", "direction"),
    "booted_page": ("booted.error", "booted.differs_from_lifted_seqGene", "booted.head"),
}


# ── the species report — every genome compared once, shared by every test ─────────────────────
@dataclass
class SpeciesReport:
    species_key: str
    every_gene: bool
    genome_count: int
    gene_count: int  # the database's own count for the published pangenome
    whole_genomes: tuple[str, ...]
    genomes_compared: int = 0
    loci_compared: int = 0
    genomes_with_copies_above_one: int = 0
    contig_label_shape_holds: int = 0
    coverage: Counter = field(default_factory=Counter)
    differences: Counter = field(default_factory=Counter)
    differences_by_strand: Counter = field(default_factory=Counter)
    examples: dict = field(default_factory=lambda: defaultdict(list))
    elapsed_seconds: float = 0.0

    def absorb(self, genome: GenomeReport) -> None:
        self.genomes_compared += 1
        self.loci_compared += genome.loci_compared
        self.genomes_with_copies_above_one += genome.coverage["loci_with_copies_above_one"] > 0
        self.contig_label_shape_holds += genome.contig_label_shape_holds
        self.coverage.update(genome.coverage)
        self.differences.update(genome.differences)
        self.differences_by_strand.update(genome.differences_by_strand)
        for name, found in genome.examples.items():
            self.examples[name].extend(found[: 3 - len(self.examples[name])])

    def summary(self) -> str:
        c = self.coverage
        whole = ", ".join(self.whole_genomes)
        mode = "EVERY GENE" if self.every_gene else f"default — whole genomes {whole} + a sample"
        return (
            f"T6 {self.species_key} ({mode}, {self.elapsed_seconds:.0f} s): "
            f"genomes {self.genomes_compared} of {self.genome_count} · "
            f"placement compared on {c['gene_table']:,} of {self.gene_count:,} genes · "
            f"rendered tab compared on {c['genes']:,} of {self.gene_count:,} genes over {self.loci_compared:,} "
            f"(genome, locus) pairs, {c['loci_with_no_gene_on_either_side']:,} of them 'no gene here' on both "
            f"sides — minus-strand {c['minus_strand_genes']:,}, "
            f"ρ>1 {c['loci_with_copies_above_one']:,} pairs / {c['genes_in_copies_above_one']:,} genes in "
            f"{self.genomes_with_copies_above_one} genomes, upstream truncated {c['upstream_flank_truncated']:,} "
            f"(minus-strand {c['minus_strand_upstream_flank_truncated']:,}), downstream truncated "
            f"{c['downstream_flank_truncated']:,}, an empty flank {c['an_empty_flank']:,}, an N in the "
            f"sequence {c['containing_an_ambiguous_base']:,}, on a contig whose name ≠ index+1 "
            f"{c['genes_on_a_contig_whose_name_is_not_index_plus_one']:,} (placement: "
            f"{c['gene_table_contig_name_is_not_index_plus_one']:,}) · booted page: {c['booted_visits']:,} "
            f"visits, {c['booted_genes']:,} genes, {c['booted_no_gene_answers']:,} 'no gene here', "
            f"{c['booted_visits_with_copies_above_one']:,} with ρ>1"
            # ⚠ Named here too, so a coverage failure caused by a defect (a dropped copy shortens the
            # count) still shows the defect rather than only the shortfall.
            + (f" · DIFFERENCES FOUND: {dict(self.differences)}" if self.differences else "")
        )

    def differences_in(self, group: str) -> str | None:
        """A message naming every difference in `group`, or `None` if there are none."""
        found = {name: self.differences[name] for name in FIELD_GROUPS[group] if self.differences[name]}
        if not found:
            return None
        lines = [f"{self.summary()}", f"differences in {group}:"]
        for name, count in found.items():
            by_strand = {
                strand or "n/a": n
                for (field_name, strand), n in self.differences_by_strand.items()
                if field_name == name
            }
            lines.append(f"  {name}: {count:,} ({by_strand})")
            for example in self.examples[name]:
                lines.append(f"      {example}")
        swapped = self.differences["flanks_swapped_end_for_end"]
        if swapped:
            lines.append(
                f"  ⛔ on {swapped:,} genes the API's UPSTREAM flank is the page's DOWNSTREAM and vice versa — "
                "the flanks are being taken from the wrong ends"
            )
        return "\n".join(lines)


def _require(condition, reason: str) -> None:
    if not condition:
        pytest.skip(reason)


@pytest.fixture(scope="module", params=SPECIES_KEYS)
def species_setup(request):
    """Everything both sides need, checked — or a skip naming exactly what is missing."""
    species_key = request.param
    node = shutil.which("node")
    _require(node, "node is not on PATH — the frozen page cannot be run")
    _require(FROZEN_APP_JS.is_file(), f"the frozen page is not at {FROZEN_APP_JS} (BACATLAS_FROZEN_APP_JS)")
    _require(FROZEN_DOM_SHIM.is_file(), f"nuna's DOM shim is not at {FROZEN_DOM_SHIM} (BACATLAS_FROZEN_DOM_SHIM)")
    payload_path = PUBLISHED_SITE_CATALOGUE_DIR / f"{species_key}.json"
    _require(payload_path.is_file(), f"the published catalogue is not at {payload_path}")
    _require(SEQUENCE_DIR.is_dir(), f"the retired .nseq files are not at {SEQUENCE_DIR}")
    url = os.environ.get("BACATLAS_DATABASE_URL")
    _require(url, "BACATLAS_DATABASE_URL is not set — T6 runs against the loaded database")
    gff_root = os.environ.get("BACATLAS_ROOT_GFF")
    _require(gff_root, "BACATLAS_ROOT_GFF is not set — the API reads the original GFFs")
    _require(Path(gff_root).is_dir(), f"BACATLAS_ROOT_GFF={gff_root} is not a directory")

    payload = json.loads(payload_path.read_text())
    engine = create_engine(url, future=True)
    request.addfinalizer(engine.dispose)
    try:
        with engine.connect() as connection:
            connection.execute(text("SELECT 1"))
    except Exception as error:  # noqa: BLE001
        pytest.skip(f"the database at {url} is unreachable: {error}")

    with Session(engine) as session:
        species = session.execute(
            select(PathogenSpecies).where(PathogenSpecies.species_key == species_key)
        ).scalar_one_or_none()
        _require(species is not None and species.default_pangenome_id, f"{species_key} has no published pangenome")
        pangenome = session.get(Pangenome, species.default_pangenome_id)
        # ⛔ The page and the database must describe the SAME model, or every locus compares two answers
        # to different questions.
        assert pangenome.run_id == payload["meta"]["model_id"], (
            f"the database publishes {pangenome.run_id} and the page {payload['meta']['model_id']}"
        )
        labels = list(
            session.execute(
                select(Locus.node_label)
                .where(Locus.pangenome_id == pangenome.pangenome_id)
                .order_by(Locus.catalogue_ordinal)
            ).scalars()
        )
        assert labels == [str(label) for label in payload["nodes"]["label"]], (
            "the locus ORDER differs between the database and the page — every `.loci` index would name a "
            "different locus"
        )
        genomes = session.execute(
            text(
                "SELECT g.sample_id, g.genome_id FROM genome g WHERE EXISTS (SELECT 1 FROM gene_locus_membership m "
                "WHERE m.genome_id = g.genome_id AND m.pangenome_id = :p) ORDER BY g.sample_id"
            ),
            {"p": pangenome.pangenome_id},
        ).all()
        gene_count = session.execute(
            text("SELECT count(*) FROM gene_locus_membership WHERE pangenome_id = :p"), {"p": pangenome.pangenome_id}
        ).scalar_one()
        whole_genomes = _whole_genomes(session, pangenome.pangenome_id)

    samples = sorted(sample for sample, _ in genomes)
    assert sorted(payload["meta"]["genomes"]) == samples, "the page and the database hold different genomes"
    missing = [
        sample
        for sample in samples
        if not all((SEQUENCE_DIR / f"{sample}{suffix}").is_file() for suffix in (".nseq", ".loci"))
    ]
    assert not missing, f"{len(missing)} genomes have no .nseq/.loci in {SEQUENCE_DIR}: {missing[:5]}"
    assert payload["meta"]["seq"]["flank"] == FLANK_LENGTH, (
        f"the page's flank is {payload['meta']['seq']['flank']} and the API's {FLANK_LENGTH}"
    )
    return {
        "species_key": species_key,
        "node": node,
        "url": url,
        "gff_root": gff_root,
        "payload_path": payload_path,
        "pangenome_id": pangenome.pangenome_id,
        "genomes": genomes,
        "gene_count": gene_count,
        "whole_genomes": whole_genomes,
    }


def _whole_genomes(session, pangenome_id: int) -> tuple[str, ...]:
    """The genomes default mode renders WHOLE, one per thing that stresses the tab — the most ρ > 1
    loci, the most genes within a flank of a contig end, the most genes on a contig whose name is not
    its index + 1, and the fewest contigs. Deterministic: ties go to the higher sample id. Distinct,
    so a genome that wins twice is rendered once."""
    rows = session.execute(
        text(
            """
            WITH genes AS (
                SELECT g.genome_id, gm.sample_id, gene.start_position, gene.end_position, c.length_bases,
                       c.contig_index, c.contig_name
                FROM gene_locus_membership g
                JOIN gene ON gene.genome_id = g.genome_id AND gene.flat_index = g.flat_index
                JOIN genome gm ON gm.genome_id = g.genome_id
                JOIN genome_contig c ON c.genome_id = gene.genome_id AND c.contig_index = gene.contig_index
                WHERE g.pangenome_id = :p
            ), copies AS (
                SELECT genome_id, count(*) AS pairs FROM (
                    SELECT genome_id, locus_id FROM gene_locus_membership WHERE pangenome_id = :p
                    GROUP BY genome_id, locus_id HAVING count(*) > 1
                ) AS multi GROUP BY genome_id
            )
            SELECT genes.sample_id,
                   coalesce(max(copies.pairs), 0) AS copy_pairs,
                   count(*) FILTER (WHERE start_position <= :flank OR end_position + :flank > length_bases) AS edge,
                   count(*) FILTER (WHERE contig_name <> 'contig' || lpad((contig_index + 1)::text, 5, '0')) AS renamed,
                   count(DISTINCT contig_index) AS contigs
            FROM genes LEFT JOIN copies ON copies.genome_id = genes.genome_id
            GROUP BY genes.sample_id
            """
        ),
        {"p": pangenome_id, "flank": FLANK_LENGTH},
    ).all()
    chosen = [
        max(rows, key=lambda r: (r.copy_pairs, r.sample_id)).sample_id,
        max(rows, key=lambda r: (r.edge, r.sample_id)).sample_id,
        max(rows, key=lambda r: (r.renamed, r.sample_id)).sample_id,
        min(rows, key=lambda r: (r.contigs, r.sample_id)).sample_id,
    ]
    return tuple(dict.fromkeys(chosen))


def _task(setup, sample_id: str, genome_id: int, every_locus: bool) -> GenomeTask:
    return GenomeTask(
        species_key=setup["species_key"],
        sample_id=sample_id,
        genome_id=genome_id,
        pangenome_id=setup["pangenome_id"],
        database_url=setup["url"],
        gff_root=setup["gff_root"],
        app_js=str(FROZEN_APP_JS),
        payload=str(setup["payload_path"]),
        seq_dir=str(SEQUENCE_DIR),
        recorder=str(RECORDER),
        node=setup["node"],
        dom_shim=str(FROZEN_DOM_SHIM),
        flank_length=FLANK_LENGTH,
        every_locus=every_locus,
    )


@pytest.fixture(scope="module")
def species_report(species_setup) -> SpeciesReport:
    """Every genome of the species, compared in parallel worker processes, then summed."""
    setup = species_setup
    report = SpeciesReport(
        species_key=setup["species_key"],
        every_gene=EVERY_GENE,
        genome_count=len(setup["genomes"]),
        gene_count=setup["gene_count"],
        whole_genomes=setup["whole_genomes"],
    )
    tasks = [
        _task(setup, sample, genome_id, EVERY_GENE or sample in setup["whole_genomes"])
        for sample, genome_id in setup["genomes"]
    ]
    started = time.monotonic()
    # ⚠ `spawn`, not `fork`: the parent holds open database connections, and a forked child inheriting
    # them would share a socket with its parent.
    with ProcessPoolExecutor(WORKERS, mp_context=multiprocessing.get_context("spawn")) as pool:
        for genome in pool.map(compare_genome, tasks):
            report.absorb(genome)
    report.elapsed_seconds = time.monotonic() - started
    print("\n" + report.summary())
    return report


def _assert_coverage(report: SpeciesReport) -> None:
    """⛔ Called FIRST by every test below: what was compared, before anything about what differed."""
    c = report.coverage
    assert report.genomes_compared == report.genome_count, report.summary()
    # Every gene's placement, in both modes — the page's own decoded tables against the database.
    assert c["gene_table"] == report.gene_count, report.summary()
    if report.every_gene:
        assert c["genes"] == report.gene_count, report.summary()
    else:
        assert c["genes"] >= 4_000, f"the default sample rendered only {c['genes']:,} genes — {report.summary()}"
    # About half of all genes are on the minus strand, and that is where a flank goes wrong.
    assert c["minus_strand_genes"] >= 0.4 * c["genes"], report.summary()
    assert c["plus_strand_genes"] > 0, report.summary()
    assert c["loci_with_copies_above_one"] >= 100, report.summary()
    assert c["upstream_flank_truncated"] > 0 and c["downstream_flank_truncated"] > 0, report.summary()
    assert c["minus_strand_upstream_flank_truncated"] > 0, report.summary()
    assert c["an_empty_flank"] > 0, report.summary()
    assert c["gene_table_contig_name_is_not_index_plus_one"] > 0, report.summary()
    assert c["genes_on_a_contig_whose_name_is_not_index_plus_one"] > 0, report.summary()
    assert c["booted_visits"] >= report.genome_count, report.summary()
    assert c["booted_no_gene_answers"] >= report.genome_count, report.summary()
    assert c["booted_visits_with_copies_above_one"] > 0, report.summary()


# ── the tests — each asserts coverage first ───────────────────────────────────────────────────
def test_T6_every_genome_and_every_genes_placement_is_compared_before_anything_is_reported(species_report):
    _assert_coverage(species_report)


def test_T6_every_gene_sits_on_the_same_contig_coordinates_strand_and_locus_on_both_sides(species_report):
    """Every gene, both modes: the page's decoded `.nseq` gene table and `.loci` against the database."""
    _assert_coverage(species_report)
    assert species_report.differences_in("placement") is None, species_report.differences_in("placement")


def test_T6_both_flanks_are_read_in_the_genes_reading_direction_on_both_strands(species_report):
    """⛔⛔ The defect of record: on a minus gene the upstream flank is the HIGHER coordinates, reversed."""
    _assert_coverage(species_report)
    assert species_report.differences_in("flanks") is None, species_report.differences_in("flanks")


def test_T6_the_coding_sequence_protein_and_gc_are_the_strings_the_page_showed(species_report):
    """⭐ GC is compared at full precision: the page's own `gcPct` on the CDS the page DISPLAYED."""
    _assert_coverage(species_report)
    assert species_report.differences_in("coding_and_protein") is None, species_report.differences_in(
        "coding_and_protein"
    )


def test_T6_a_genome_at_rho_above_one_gets_every_copy_with_the_pages_ordinals(species_report):
    _assert_coverage(species_report)
    assert species_report.differences_in("copies") is None, species_report.differences_in("copies")


def test_T6_a_flank_cut_short_by_a_contig_end_is_reported_per_gene_exactly_as_the_page_reported_it(species_report):
    _assert_coverage(species_report)
    assert species_report.differences_in("truncation") is None, species_report.differences_in("truncation")


def test_T6_the_contig_the_page_named_is_the_apis_seqid_on_every_gene(species_report):
    """⛔ The page printed the GFF seqid from the `.nseq` header — never `contig_index + 1`."""
    _assert_coverage(species_report)
    assert species_report.differences_in("contig_and_strand") is None, species_report.differences_in(
        "contig_and_strand"
    )


def test_T6_the_new_cards_contig_row_differs_only_by_the_named_exception_and_on_every_gene(species_report):
    """⚠ A NAMED display difference, asserted exactly: the card prints `contig_name`, the page printed
    the seqid, and `seqid == f"{sample}.{contig_name}"` on every gene compared — not on most."""
    _assert_coverage(species_report)
    assert species_report.species_key in CONTIG_ROW_SHOWS_NAME_NOT_SEQID.species_keys
    assert species_report.contig_label_shape_holds == species_report.coverage["genes"], species_report.summary()


def test_T6_the_booted_page_draws_exactly_what_the_lifted_seqGene_draws_and_says_no_gene_where_the_api_does(
    species_report,
):
    """The lifted functions are only an oracle if they draw what the page draws — shown on every genome."""
    _assert_coverage(species_report)
    assert species_report.differences_in("booted_page") is None, species_report.differences_in("booted_page")


def test_T6_every_difference_the_comparison_can_find_is_owned_by_exactly_one_test(species_report):
    """⛔ A difference found and then reported by NO test is the quiet version of "0 differ"."""
    owned = [name for names in FIELD_GROUPS.values() for name in names]
    assert len(owned) == len(set(owned)), "a field is owned by two groups"
    rendered = set(RENDERED_FIELDS) | {"copies_at_locus", "flanks_swapped_end_for_end"}
    assert rendered | {f"table.{n}" for n in GENE_TABLE_FIELDS} | set(FIELD_GROUPS["booted_page"]) == set(owned)
    unowned = set(species_report.differences) - set(owned)
    assert not unowned, f"differences no test reports: {unowned}"


# ── the endpoint itself, over HTTP ────────────────────────────────────────────────────────────
@pytest.fixture(scope="module")
def http_client():
    from bacatlas_backend.application_factory import create_application
    from bacatlas_backend.configuration import Configuration

    if not os.environ.get("BACATLAS_DATABASE_URL") or not os.environ.get("BACATLAS_ROOT_GFF"):
        pytest.skip("BACATLAS_DATABASE_URL and BACATLAS_ROOT_GFF are both needed for the endpoint")
    return create_application(Configuration.from_environment()).test_client()


def test_T6_the_http_endpoint_serves_what_the_page_showed_including_no_gene_here(species_setup, http_client):
    """The sweep calls the service; this is the same comparison through the URL, so the wiring is covered.

    Two stress genomes, their stratified sample of loci, and two loci each genome LACKS — which must be
    a 200 with `genes: []`, the API's version of the page's "has no gene at".
    """
    setup = species_setup
    report = GenomeReport(sample_id="", every_locus=False)
    no_gene_answers = 0
    genome_ids = dict(setup["genomes"])
    engine = create_engine(setup["url"], future=True)
    stress_genomes = setup["whole_genomes"][:2]
    try:
        for sample_id in stress_genomes:
            report.sample_id = sample_id
            with Session(engine) as session:
                rows = genome_gene_rows(session, genome_id=genome_ids[sample_id], pangenome_id=setup["pangenome_id"])
                lacking = absent_labels(session, genome_id=genome_ids[sample_id], pangenome_id=setup["pangenome_id"])
            labels = sample_labels(rows, FLANK_LENGTH, sample_id)
            page = run_recorder(
                setup["node"],
                str(RECORDER),
                "render",
                {
                    "app_js": str(FROZEN_APP_JS),
                    "payload": str(setup["payload_path"]),
                    "seq_dir": str(SEQUENCE_DIR),
                    "sample": sample_id,
                    "labels": labels + lacking,
                },
            )
            by_label = defaultdict(list)
            for record in page["rendered"]:
                if not record.get("absent"):
                    by_label[record["label"]].append(record)
            for label in labels + lacking:
                response = http_client.get(
                    f"/api/v1/species/{setup['species_key']}/genomes/{sample_id}/loci/{label}/sequence"
                )
                assert response.status_code == 200, (label, response.get_json())
                genes = response.get_json()["genes"]
                if label in lacking:
                    assert genes == [] and not by_label[label], f"{sample_id} ·{label}: the page and the API disagree"
                    no_gene_answers += 1
                compare_rendered(report, label=label, page_genes=by_label[label], api_genes=genes)
    finally:
        engine.dispose()

    c = report.coverage
    assert c["genes"] >= 60 and c["minus_strand_genes"] > 0 and c["loci_with_copies_above_one"] > 0, dict(c)
    assert c["upstream_flank_truncated"] > 0 and c["downstream_flank_truncated"] > 0, dict(c)
    assert no_gene_answers == 2 * len(stress_genomes), no_gene_answers
    assert not report.differences, {name: report.examples[name] for name in report.differences}
