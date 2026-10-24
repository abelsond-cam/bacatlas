"""Parity suite T6's machinery — the frozen page's Sequence tab beside the API's, one genome at a time.

Kept out of the test module because it runs in WORKER PROCESSES: a full sweep is ~1 M genes, each of
which the service answers with one SQL statement, so the genomes are compared in parallel and each
worker must be able to import this without importing pytest's collection.

**The two sides.**

- *Before* — `frozen_page_sequence_recorder.js`, which runs the published `app.js`'s own `seqGene`
  over the committed `.nseq` + `.loci` files and reads back what the tab DISPLAYED.
- *After* — `load_gene_sequences` → `serialise_gene_sequence`, exactly what the endpoint returns,
  slicing the original gzipped Bakta GFF under `BACATLAS_ROOT_GFF`.

**The page's display is PARSED, never re-rendered.** Formatting the API's numbers into the page's
sentences would be a second implementation of the page, and the one thing a parity suite must not
contain. So every displayed string is read back into a value by a strict pattern, a string that does
not parse is itself a difference, and the values are compared with the API's fields.

⛔ **Coverage is counted here and asserted by the caller BEFORE any difference is reported** — genes,
minus-strand genes, ρ > 1 copies, truncated flanks, genes on contigs whose name is not their index.
This repo's own scar (`4ab35ca`): a diff loop that skips what it cannot compare reports "0 differ"
while silently skipping exactly what changed.
"""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
import zlib
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from decimal import ROUND_HALF_UP, Decimal
from functools import lru_cache
from pathlib import Path

#: How many examples of each kind of difference a report keeps. The COUNT is always exact.
EXAMPLES_KEPT = 3

#: Default-mode sample per genome: at most this many loci from each stratum. ⚠ Minus-strand and
#: contig-edge genes are deliberately over-represented — they are where a flank goes wrong.
SAMPLE_PER_STRATUM = {
    "copies_above_one": 8,
    "edge_low_plus": 4,  # a + gene near a contig START: its UPSTREAM flank is short
    "edge_high_plus": 4,  # a + gene near a contig END: its DOWNSTREAM flank is short
    "edge_low_minus": 4,  # a − gene near a contig START: its DOWNSTREAM flank is short
    "edge_high_minus": 4,  # a − gene near a contig END: its UPSTREAM flank is short
    "minus_interior": 10,
    "plus_interior": 6,
    "contig_name_is_not_index_plus_one": 4,
    # ⚠ A flank that reaches EXACTLY to a contig end, or falls one base short of it: the only genes on
    # which `< 1` and `<= 1` disagree. Without this stratum an off-by-one truncation flag was caught
    # by ONE gene of 19,041 in the default sample.
    "flank_ends_at_the_contig_end_or_one_short": 4,
}


@dataclass(frozen=True)
class GenomeTask:
    """Everything a worker needs to compare one genome. Plain values, so it pickles."""

    species_key: str
    sample_id: str
    genome_id: int
    pangenome_id: int
    database_url: str
    gff_root: str
    app_js: str
    payload: str
    seq_dir: str
    recorder: str
    node: str
    #: `nuna/tests/js/dom_shim.js` — the DOM the whole page is booted into.
    dom_shim: str
    flank_length: int
    #: `True` renders EVERY locus the genome sits at; `False` renders the stratified sample.
    every_locus: bool


@dataclass
class GenomeReport:
    """What one genome contributed: exact counts, and a few examples of each difference."""

    sample_id: str
    every_locus: bool
    #: `gene_table` counts genes whose placement was compared (every gene, both modes).
    #: The rest count genes whose RENDERED tab was compared.
    coverage: Counter = field(default_factory=Counter)
    differences: Counter = field(default_factory=Counter)
    #: Differences split by strand, so "every minus gene and no plus gene" is visible at once.
    differences_by_strand: Counter = field(default_factory=Counter)
    examples: dict = field(default_factory=lambda: defaultdict(list))
    #: `(genome, locus)` pairs compared — including those where BOTH sides say "no gene here".
    loci_compared: int = 0
    #: ⭐ The named display exception's exact shape, counted: the page printed the seqid, the new card
    #: prints `contig_name`, and `seqid == f"{sample}.{contig_name}"` on every one.
    contig_label_shape_holds: int = 0


# ── the recorder ──────────────────────────────────────────────────────────────────────────────
def run_recorder(node: str, recorder: str, mode: str, request: dict) -> dict:
    """Run the Node recorder with `request` on stdin; its stdout is one JSON document."""
    completed = subprocess.run(
        [node, recorder, mode, "-"],
        input=json.dumps(request),
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"the frozen-page recorder failed ({mode}):\n{completed.stderr[-4000:]}")
    return json.loads(completed.stdout)


def sequence_digest(sequence: str) -> dict:
    """The recorder's own `{length, digest, head}`, computed on the API's string."""
    return {
        "length": len(sequence),
        "digest": hashlib.sha1(sequence.encode()).hexdigest()[:16],
        "head": sequence[:12],
    }


def js_to_fixed_1(value: float) -> str:
    """`Number.prototype.toFixed(1)` for a non-negative double: exact binary value, ties AWAY from zero.

    ⚠ Not Python's `f"{x:.1f}"`, which rounds a tie to even — 50.25 is exactly representable and the
    two disagree on it. The page printed with `toFixed`, and so does the new card.
    """
    return str(Decimal(value).quantize(Decimal("0.1"), rounding=ROUND_HALF_UP))


# ── reading the page's display back into values ───────────────────────────────────────────────
_COPY_HEADING = re.compile(r"^Copy (\d+) of (\d+) at this locus$")
_SPAN_STAT = re.compile(r"^([\d,]+)–([\d,]+) \(1-based, inclusive\)$")
_LENGTH_STAT = re.compile(r"^([\d,]+) bp · ([\d,]+) aa$")
_GC_STAT = re.compile(r"^(\d+\.\d)%$")
_FLANK_NOTE = re.compile(r"^Contig ([\d,]+)–([\d,]+)(, reverse-complemented)?$")
_CDS_NOTE = re.compile(r"^Contig ([\d,]+)–([\d,]+)(, reverse-complemented)?\. Includes the stop codon\.$")
_CDS_HEADING = re.compile(r"^The gene — ([\d,]+) bases$")
_PROTEIN_HEADING = re.compile(r"^Protein — ([\d,]+) residues$")
_EDGE = re.compile(
    r"^This gene sits within (\d+) bases of the end of (\S+), so the (upstream|downstream|upstream and downstream) "
    r"flanks? (?:is|are) short\."
)
_STRAND_STAT = {"− (reverse)": "-", "+ (forward)": "+"}
_DIRECTION_STAT = {
    "read as the reverse complement of the contig": "-",
    "read as the contig is written": "+",
}
_PARTIAL_NOTE = "5′-partial — this CDS does not begin at a start codon"
_NONE_IN_ASSEMBLY = "None in the assembly."
_NOT_ACGT = re.compile("[^ACGT]")


def _number(text: str) -> int:
    return int(text.replace(",", ""))


class _Unparsed:
    """A displayed string no pattern could read — itself a difference, never silently skipped."""

    def __init__(self, text):
        self.text = text

    def __eq__(self, other):
        return False

    def __hash__(self):
        return hash(self.text)

    def __repr__(self):
        return f"<unparsed {self.text!r}>"


def _match(pattern, text):
    found = pattern.match(text) if text is not None else None
    return found if found else None


def read_page_gene(record: dict) -> dict:
    """One rendered `.seq-gene`, as values in the API's vocabulary."""
    stats = dict(record["stats"])
    out: dict = {}

    heading = record["copy_heading"]
    if heading is None:
        # The page draws no heading for a single copy — so `None` MEANS "1 of 1", and says so.
        out["copy_ordinal"], out["copy_count"] = 1, 1
    elif found := _match(_COPY_HEADING, heading):
        out["copy_ordinal"], out["copy_count"] = int(found[1]), int(found[2])
    else:
        out["copy_ordinal"] = out["copy_count"] = _Unparsed(heading)
    out["flat_index"] = record["flat_index"]

    out["seqid"] = stats.get("Contig", _Unparsed(None))
    if found := _match(_SPAN_STAT, stats.get("Span")):
        out["start_position"], out["end_position"] = _number(found[1]), _number(found[2])
    else:
        out["start_position"] = out["end_position"] = _Unparsed(stats.get("Span"))
    out["strand"] = _STRAND_STAT.get(stats.get("Strand"), _Unparsed(stats.get("Strand")))
    out["direction"] = _DIRECTION_STAT.get(stats.get("Direction"), _Unparsed(stats.get("Direction")))
    if found := _match(_LENGTH_STAT, stats.get("Length")):
        out["length_stat"] = (_number(found[1]), _number(found[2]))
    else:
        out["length_stat"] = _Unparsed(stats.get("Length"))
    out["gc_percent_shown"] = (
        found[1] if (found := _match(_GC_STAT, stats.get("GC content"))) else _Unparsed(stats.get("GC content"))
    )
    note = stats.get("Note")
    out["is_five_prime_partial"] = True if note == _PARTIAL_NOTE else (False if note is None else _Unparsed(note))
    out["gc_percent"] = record["gc_of_shown_cds"]

    blocks = record["blocks"]
    for kind, name in (("up", "upstream_flank"), ("down", "downstream_flank")):
        block = blocks[kind]
        out[f"{name}_sequence"] = _block_sequence(block)
        if block["note"] is None:
            out[f"{name}_span"] = None
            out[f"{name}_is_reverse_complemented"] = None
        elif found := _match(_FLANK_NOTE, block["note"]):
            out[f"{name}_span"] = (_number(found[1]), _number(found[2]))
            out[f"{name}_is_reverse_complemented"] = found[3] is not None
        else:
            out[f"{name}_span"] = out[f"{name}_is_reverse_complemented"] = _Unparsed(block["note"])
    cds = blocks["cds"]
    out["coding_sequence"] = _block_sequence(cds)
    if found := _match(_CDS_NOTE, cds["note"]):
        out["coding_sequence_span"] = (_number(found[1]), _number(found[2]))
        out["coding_sequence_is_reverse_complemented"] = found[3] is not None
    else:
        out["coding_sequence_span"] = out["coding_sequence_is_reverse_complemented"] = _Unparsed(cds["note"])
    out["coding_sequence_heading_length"] = (
        _number(found[1]) if (found := _match(_CDS_HEADING, cds["heading"])) else _Unparsed(cds["heading"])
    )
    protein = blocks["aa"]
    out["protein_sequence"] = _block_sequence(protein)
    out["protein_heading_length"] = (
        _number(found[1]) if (found := _match(_PROTEIN_HEADING, protein["heading"])) else _Unparsed(protein["heading"])
    )

    edge = record["edge"]
    if edge is None:
        out["upstream_flank_is_truncated_by_contig_end"] = False
        out["downstream_flank_is_truncated_by_contig_end"] = False
        out["edge_contig"] = None
    elif found := _match(_EDGE, edge):
        out["upstream_flank_is_truncated_by_contig_end"] = "upstream" in found[3]
        out["downstream_flank_is_truncated_by_contig_end"] = "downstream" in found[3]
        out["edge_contig"] = found[2]
    else:
        out["upstream_flank_is_truncated_by_contig_end"] = _Unparsed(edge)
        out["downstream_flank_is_truncated_by_contig_end"] = _Unparsed(edge)
        out["edge_contig"] = _Unparsed(edge)
    return out


def _block_sequence(block: dict):
    """The block's bases as a digest, or `""`'s digest where the page said "None in the assembly."."""
    if block["sequence"] is not None:
        return block["sequence"]
    if block["none"] == _NONE_IN_ASSEMBLY:
        return sequence_digest("")
    return _Unparsed(block["none"])


def read_api_gene(gene: dict) -> dict:
    """One serialised API gene, in the same vocabulary as `read_page_gene`."""
    minus = gene["strand"] == "-"
    truncated_up = gene["upstream_flank_is_truncated_by_contig_end"]
    truncated_down = gene["downstream_flank_is_truncated_by_contig_end"]
    return {
        "copy_ordinal": gene["copy_ordinal"],
        "copy_count": gene["copy_count"],
        "flat_index": gene["flat_index"],
        "seqid": gene["seqid"],
        "start_position": gene["start_position"],
        "end_position": gene["end_position"],
        "strand": gene["strand"],
        "direction": gene["strand"],
        "length_stat": (len(gene["coding_sequence"]), len(gene["protein_sequence"])),
        "gc_percent_shown": js_to_fixed_1(gene["gc_percent"]),
        "is_five_prime_partial": gene["is_five_prime_partial"],
        "gc_percent": gene["gc_percent"],
        "upstream_flank_sequence": sequence_digest(gene["upstream_flank_sequence"]),
        "upstream_flank_span": tuple(gene["upstream_flank_span"]) if gene["upstream_flank_span"] else None,
        "upstream_flank_is_reverse_complemented": minus if gene["upstream_flank_span"] else None,
        "downstream_flank_sequence": sequence_digest(gene["downstream_flank_sequence"]),
        "downstream_flank_span": tuple(gene["downstream_flank_span"]) if gene["downstream_flank_span"] else None,
        "downstream_flank_is_reverse_complemented": minus if gene["downstream_flank_span"] else None,
        "coding_sequence": sequence_digest(gene["coding_sequence"]),
        "coding_sequence_span": (gene["start_position"], gene["end_position"]),
        "coding_sequence_is_reverse_complemented": minus,
        "coding_sequence_heading_length": len(gene["coding_sequence"]),
        "protein_sequence": sequence_digest(gene["protein_sequence"]),
        "protein_heading_length": len(gene["protein_sequence"]),
        "upstream_flank_is_truncated_by_contig_end": truncated_up,
        "downstream_flank_is_truncated_by_contig_end": truncated_down,
        # ⚠ The page named the contig in the edge sentence by its seqid; the API's seqid is that value.
        "edge_contig": gene["seqid"] if (truncated_up or truncated_down) else None,
    }


#: Every field compared on a rendered gene. ⛔ Asserted to be the key set of BOTH readers before any
#: comparison, so a field one side stops producing cannot drop out of the comparison unnoticed.
RENDERED_FIELDS = (
    "copy_ordinal",
    "copy_count",
    "flat_index",
    "seqid",
    "start_position",
    "end_position",
    "strand",
    "direction",
    "length_stat",
    "gc_percent_shown",
    "is_five_prime_partial",
    "gc_percent",
    "upstream_flank_sequence",
    "upstream_flank_span",
    "upstream_flank_is_reverse_complemented",
    "downstream_flank_sequence",
    "downstream_flank_span",
    "downstream_flank_is_reverse_complemented",
    "coding_sequence",
    "coding_sequence_span",
    "coding_sequence_is_reverse_complemented",
    "coding_sequence_heading_length",
    "protein_sequence",
    "protein_heading_length",
    "upstream_flank_is_truncated_by_contig_end",
    "downstream_flank_is_truncated_by_contig_end",
    "edge_contig",
)

#: The fields that make up "the page's placement of every gene" — compared on EVERY gene in both modes.
GENE_TABLE_FIELDS = (
    "gene_count",
    "start_position",
    "end_position",
    "contig_index",
    "seqid",
    "contig_length",
    "strand",
    "gff_phase",
    "is_five_prime_partial",
    "catalogue_ordinal",
)


def compare_rendered(
    report: GenomeReport,
    *,
    label: str,
    page_genes: list[dict],
    api_genes: list[dict],
    renamed_contig_genes: frozenset[int] = frozenset(),
) -> None:
    """Compare one `(genome, locus)` field by field, counting coverage as it goes.

    `renamed_contig_genes` — the `flat_index` of every gene on a contig whose name is not
    `contig{index + 1}` — is only counted, so the report can say the rendered comparison reached the
    28.6 % of contigs where naming by index would have been wrong.
    """
    report.loci_compared += 1
    report.coverage["loci_with_no_gene_on_either_side"] += not page_genes and not api_genes
    if len(page_genes) != len(api_genes):
        _differ(report, "copies_at_locus", None, label, None, len(page_genes), len(api_genes))
    for page_record, api_gene in zip(page_genes, api_genes, strict=False):
        page, api = read_page_gene(page_record), read_api_gene(api_gene)
        assert set(page) == set(api) == set(RENDERED_FIELDS), (
            f"the two readers disagree on WHICH fields they compare: page-only {set(page) - set(api)}, "
            f"api-only {set(api) - set(page)} — a field missing from one side would drop out unnoticed"
        )
        strand = api_gene["strand"]
        _count_coverage(report, page, api_gene)
        report.coverage["genes_on_a_contig_whose_name_is_not_index_plus_one"] += (
            api_gene["flat_index"] in renamed_contig_genes
        )
        if api_gene["seqid"] == f"{report.sample_id}.{api_gene['contig_name']}":
            report.contig_label_shape_holds += 1
        for name in RENDERED_FIELDS:
            if page[name] != api[name]:
                _differ(report, name, strand, label, api_gene["flat_index"], page[name], api[name])
        # ⛔ The failure T6 exists for, named when it happens rather than left to be inferred: the API's
        # upstream is the page's downstream and vice versa.
        if (
            page["upstream_flank_sequence"] != api["upstream_flank_sequence"]
            and page["upstream_flank_sequence"] == api["downstream_flank_sequence"]
            and page["downstream_flank_sequence"] == api["upstream_flank_sequence"]
        ):
            _differ(
                report, "flanks_swapped_end_for_end", strand, label, api_gene["flat_index"],
                f"up {page['upstream_flank_sequence']['head']}… down {page['downstream_flank_sequence']['head']}…",
                f"up {api['upstream_flank_sequence']['head']}… down {api['downstream_flank_sequence']['head']}…",
            )


def _count_coverage(report: GenomeReport, page: dict, api_gene: dict) -> None:
    """⚠ Counted from the PAGE's side, which is frozen: an API defect must show up as a DIFFERENCE,
    never as a coverage shortfall that hides it (an API that stopped flagging truncation would
    otherwise read as "this sample has no contig-edge genes")."""
    coverage = report.coverage
    minus = page["strand"] == "-"
    copies = page["copy_count"] if isinstance(page["copy_count"], int) else 0
    truncated_up = page["upstream_flank_is_truncated_by_contig_end"] is True
    coverage["genes"] += 1
    coverage["minus_strand_genes"] += minus
    coverage["plus_strand_genes"] += page["strand"] == "+"
    coverage["genes_in_copies_above_one"] += copies > 1
    coverage["loci_with_copies_above_one"] += copies > 1 and page["copy_ordinal"] == 1
    coverage["upstream_flank_truncated"] += truncated_up
    coverage["downstream_flank_truncated"] += page["downstream_flank_is_truncated_by_contig_end"] is True
    coverage["minus_strand_upstream_flank_truncated"] += minus and truncated_up
    coverage["an_empty_flank"] += page["upstream_flank_span"] is None or page["downstream_flank_span"] is None
    # Informational only — the page's digests cannot say which letters they hold, so this one is read
    # from the API's strings, and nothing asserts on it.
    coverage["containing_an_ambiguous_base"] += any(
        _NOT_ACGT.search(api_gene[name])
        for name in ("upstream_flank_sequence", "coding_sequence", "downstream_flank_sequence")
    )


def _differ(report, name, strand, label, flat_index, page_value, api_value) -> None:
    report.differences[name] += 1
    report.differences_by_strand[(name, strand)] += 1
    if len(report.examples[name]) < EXAMPLES_KEPT:
        report.examples[name].append(
            {
                "sample_id": report.sample_id,
                "label": label,
                "flat_index": flat_index,
                "strand": strand,
                "page": _short(page_value),
                "api": _short(api_value),
            }
        )


def _short(value) -> str:
    text = repr(value)
    return text if len(text) <= 160 else text[:157] + "..."


# ── the database side of the gene table ───────────────────────────────────────────────────────
_GENE_TABLE_SQL = """
SELECT g.flat_index, g.start_position, g.end_position, g.contig_index, g.strand, g.gff_phase,
       g.is_five_prime_partial, c.seqid, c.contig_name, c.length_bases, l.node_label, l.catalogue_ordinal
FROM gene g
JOIN genome_contig c ON c.genome_id = g.genome_id AND c.contig_index = g.contig_index
LEFT JOIN gene_locus_membership m
       ON m.genome_id = g.genome_id AND m.flat_index = g.flat_index AND m.pangenome_id = :pangenome_id
LEFT JOIN locus l ON l.locus_id = m.locus_id
WHERE g.genome_id = :genome_id
ORDER BY g.flat_index
"""


def compare_gene_table(report: GenomeReport, page_genome: dict, rows) -> None:
    """EVERY gene of the genome: where the page put it against where the database does."""
    n = page_genome["n_genes"]
    if n != len(rows):
        _differ(report, "table.gene_count", None, None, None, n, len(rows))
    names = page_genome["contig_names"]
    lengths = page_genome["contig_lengths"]
    for row in rows[:n]:
        j = row.flat_index
        if j >= n:
            continue
        report.coverage["gene_table"] += 1
        flags = page_genome["flags"][j]
        ci = page_genome["contig_index"][j]
        page_locus = page_genome["locus_index"][j]
        pairs = (
            ("start_position", page_genome["start"][j], row.start_position),
            ("end_position", page_genome["end"][j], row.end_position),
            ("contig_index", ci, row.contig_index),
            ("seqid", names[ci] if ci < len(names) else None, row.seqid),
            ("contig_length", lengths[ci] if ci < len(lengths) else None, row.length_bases),
            ("strand", "-" if flags & 1 else "+", row.strand),
            ("gff_phase", flags >> 2, row.gff_phase),
            ("is_five_prime_partial", bool(flags & 2), row.is_five_prime_partial),
            ("catalogue_ordinal", None if page_locus < 0 else page_locus, row.catalogue_ordinal),
        )
        for name, page_value, db_value in pairs:
            if page_value != db_value:
                _differ(report, f"table.{name}", row.strand, row.node_label, j, page_value, db_value)
        report.coverage["gene_table_contig_name_is_not_index_plus_one"] += (
            row.contig_name != f"contig{row.contig_index + 1:05d}"
        )


# ── choosing what to render ───────────────────────────────────────────────────────────────────
def sample_labels(rows, flank_length: int, sample_id: str) -> list[str]:
    """A stratified, deterministic sample of this genome's loci. Minus-strand and contig-edge genes
    over-represented; every ρ > 1 locus up to the cap. Ordered by a hash, not by position, so the
    sample spreads across contigs instead of taking the first ones."""
    placed = [row for row in rows if row.node_label is not None]
    copies = Counter(row.node_label for row in placed)
    strata: dict[str, list[str]] = defaultdict(list)
    for row in sorted(placed, key=lambda r: zlib.crc32(f"{sample_id}:{r.flat_index}".encode())):
        low = row.start_position <= flank_length
        high = row.end_position + flank_length > row.length_bases
        minus = row.strand == "-"
        if copies[row.node_label] > 1:
            strata["copies_above_one"].append(row.node_label)
        if low:
            strata["edge_low_minus" if minus else "edge_low_plus"].append(row.node_label)
        if high:
            strata["edge_high_minus" if minus else "edge_high_plus"].append(row.node_label)
        if not (low or high):
            strata["minus_interior" if minus else "plus_interior"].append(row.node_label)
        if row.contig_name != f"contig{row.contig_index + 1:05d}":
            strata["contig_name_is_not_index_plus_one"].append(row.node_label)
        if row.start_position - flank_length in (0, 1) or row.end_position + flank_length - row.length_bases in (0, 1):
            strata["flank_ends_at_the_contig_end_or_one_short"].append(row.node_label)
    chosen: dict[str, None] = {}
    for stratum, cap in SAMPLE_PER_STRATUM.items():
        for label in list(dict.fromkeys(strata[stratum]))[:cap]:
            chosen.setdefault(label, None)
    return list(chosen)


_ABSENT_LOCI_SQL = """
SELECT l.node_label FROM locus l
WHERE l.pangenome_id = :pangenome_id
  AND NOT EXISTS (
      SELECT 1 FROM gene_locus_membership m
      WHERE m.locus_id = l.locus_id AND m.genome_id = :genome_id AND m.pangenome_id = :pangenome_id
  )
ORDER BY l.member_genome_count DESC, l.catalogue_ordinal
LIMIT :limit
"""

#: Loci each genome visits where it has NO gene — the near-core ones it lacks, so "has no gene here"
#: is exercised where a reader would least expect it.
ABSENT_LOCI_PER_GENOME = 2


# ── the booted page ───────────────────────────────────────────────────────────────────────────
_BOOTED_PRESENT = re.compile(r"^(\S+) at .+ ·(\S+?)(?: — (\d+) copies here)?$")
_BOOTED_ABSENT = re.compile(r"^(\S+) has no gene at .+ ·(\S+)\. Walk to a locus it carries, or clear the anchor\.$")
#: Keys the render mode adds for identity and the booted DOM cannot show — dropped before comparing.
_RENDER_ONLY_KEYS = ("label", "locus_index", "flat_index")


def compare_booted(report: GenomeReport, booted: dict, page_by_label: dict, api_count_by_label: dict) -> None:
    """The whole booted page against the lifted `seqGene`, and its "no gene here" answer against the API."""
    for visit in booted["visits"]:
        label = visit["label"]
        report.coverage["booted_visits"] += 1
        if visit["error"] is not None:
            _differ(report, "booted.error", None, label, None, visit["error"], None)
            continue
        lifted = [{k: v for k, v in r.items() if k not in _RENDER_ONLY_KEYS} for r in page_by_label.get(label, [])]
        report.coverage["booted_genes"] += len(visit["genes"])
        if visit["genes"] != lifted:
            _differ(
                report, "booted.differs_from_lifted_seqGene", None, label, None,
                f"{len(visit['genes'])} genes drawn", f"{len(lifted)} genes from the lifted seqGene",
            )
        head = visit["head"] or ""
        api_count = api_count_by_label[label]
        if found := _BOOTED_ABSENT.match(head):
            report.coverage["booted_no_gene_answers"] += 1
            shown = (found[1], found[2], 0)
        elif found := _BOOTED_PRESENT.match(head):
            shown = (found[1], found[2], int(found[3] or 1))
            report.coverage["booted_visits_with_copies_above_one"] += found[3] is not None
        else:
            _differ(report, "booted.head", None, label, None, head, None)
            continue
        # ⛔ "has no gene here" is an ANSWER, and it must be the API's answer too: `genes: []`.
        if shown != (report.sample_id, label, api_count):
            _differ(report, "booted.head", None, label, None, shown, (report.sample_id, label, api_count))


# ── one genome, end to end ────────────────────────────────────────────────────────────────────
@lru_cache(maxsize=1)
def _engine(url: str):
    from sqlalchemy import create_engine

    return create_engine(url, future=True, pool_size=1, max_overflow=0)


def genome_gene_rows(session, *, genome_id: int, pangenome_id: int):
    """Every gene of one genome, with its contig and its locus in this pangenome (NULL if none)."""
    from sqlalchemy import text

    return session.execute(text(_GENE_TABLE_SQL), {"genome_id": genome_id, "pangenome_id": pangenome_id}).all()


def absent_labels(session, *, genome_id: int, pangenome_id: int, limit: int = ABSENT_LOCI_PER_GENOME) -> list[str]:
    """The most widely shared loci this genome has NO gene at."""
    from sqlalchemy import text

    return list(
        session.execute(
            text(_ABSENT_LOCI_SQL), {"genome_id": genome_id, "pangenome_id": pangenome_id, "limit": limit}
        ).scalars()
    )


def api_genes_at(session, *, pangenome_id: int, label: str, sample_id: str, gff_root: Path, flank_length: int):
    """⭐ Exactly the endpoint's path, minus HTTP: the service, then the serialiser.

    ⛔ **A HELD-ONLY provider, deliberately.** This is the parity oracle against the frozen page,
    whose bases came from the GFF tree — so it must read the GFF tree. Letting it fall through to
    ATB would let a parity run pass on bytes the page never saw, which is the one thing a parity
    oracle must not do.
    """
    from bacatlas_backend.sequence_source.contig_sequence_provider import ContigSequenceProvider
    from bacatlas_backend.serialisers.locus_serialiser import serialise_gene_sequence
    from bacatlas_backend.services.gene_sequence_service import load_gene_sequences

    return [
        serialise_gene_sequence(row)
        for row in load_gene_sequences(
            session,
            pangenome_id=pangenome_id,
            node_label=label,
            sample_id=sample_id,
            sequence_provider=ContigSequenceProvider(gff_root=gff_root, atb_store=None),
            flank_length=flank_length,
        )
    ]


def compare_genome(task: GenomeTask) -> GenomeReport:
    """Record the page's tab for this genome, ask the service for the same loci, compare.

    Both modes also boot the WHOLE page over the sampled loci (plus two the genome lacks), which is
    what shows the lifted `seqGene` draws exactly what the page draws.
    """
    from sqlalchemy.orm import Session

    report = GenomeReport(sample_id=task.sample_id, every_locus=task.every_locus)
    with Session(_engine(task.database_url)) as session:
        rows = genome_gene_rows(session, genome_id=task.genome_id, pangenome_id=task.pangenome_id)
        visited = sample_labels(rows, task.flank_length, task.sample_id) + absent_labels(
            session, genome_id=task.genome_id, pangenome_id=task.pangenome_id
        )
        page_request = {"app_js": task.app_js, "payload": task.payload, "seq_dir": task.seq_dir}
        page = run_recorder(
            task.node,
            task.recorder,
            "render",
            {**page_request, "sample": task.sample_id, "labels": None if task.every_locus else visited},
        )
        booted = run_recorder(
            task.node,
            task.recorder,
            "booted",
            {**page_request, "dom_shim": task.dom_shim, "visits": [[task.sample_id, label] for label in visited]},
        )
        if page["flank"] != task.flank_length:
            raise AssertionError(f"the page's flank is {page['flank']} and the API's {task.flank_length}")
        compare_gene_table(report, page["genome"], rows)

        page_by_label: dict[str, list[dict]] = defaultdict(list)
        for record in page["rendered"]:
            if not record.get("absent"):
                page_by_label[record["label"]].append(record)
        if task.every_locus:
            # ⛔ The UNION: a locus only one side places this genome at is a difference, not a skip.
            labels = sorted(
                {row.node_label for row in rows if row.node_label is not None} | set(page_by_label) | set(visited)
            )
        else:
            labels = visited
        renamed = frozenset(
            row.flat_index for row in rows if row.contig_name != f"contig{row.contig_index + 1:05d}"
        )
        api_count_by_label = {}
        for label in labels:
            api_genes = api_genes_at(
                session,
                pangenome_id=task.pangenome_id,
                label=label,
                sample_id=task.sample_id,
                gff_root=Path(task.gff_root),
                flank_length=task.flank_length,
            )
            api_count_by_label[label] = len(api_genes)
            compare_rendered(
                report,
                label=label,
                page_genes=page_by_label.get(label, []),
                api_genes=api_genes,
                renamed_contig_genes=renamed,
            )
        compare_booted(report, booted, page_by_label, api_count_by_label)
    return report
