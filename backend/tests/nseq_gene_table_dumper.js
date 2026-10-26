/* Dump what the PAGE's own decoder makes of a `.nseq` gene table and a `.loci` vector.
 *
 * ⭐ **A reference recorder, not a test.** `gene_spans` claims to be byte-identical to the `.nseq`
 * gene table and `gene_locus_ordinals` to `.loci`. The only thing that can check a claim like that
 * is a second decoder written by someone else, so this lifts `decodeNseq` straight out of the
 * published `app.js` — the precedent is `frozen_page_sequence_recorder.js`, which lifts the whole
 * Sequence tab the same way — and reports the fields it decoded. `test_gene_blocks.py` compares.
 *
 * ⛔ **Nothing is re-implemented here.** The bit positions of the flags word, the stride, the
 * header offsets: all of it comes from `decodeNseq`. The only code written in this file is the
 * canonical rendering of the decoded fields, which must stay identical to the Python side's — it is
 * the comparison, not the measurement.
 *
 * Usage:  node nseq_gene_table_dumper.js <request.json | ->
 *   request: {app_js, genomes: [{sample, nseq, loci}]}
 * Prints one JSON document on stdout: per genome, a digest over every gene's decoded fields, the
 * coverage counters that prove the digest is not over an empty table, and the first rows verbatim
 * so a mismatch is diagnosable rather than merely reported.
 */
"use strict";

const crypto = require("node:crypto");
const fs = require("node:fs");

const requestPath = process.argv[2];
const request = JSON.parse(fs.readFileSync(requestPath === "-" ? 0 : requestPath, "utf8"));
const APP_SOURCE = fs.readFileSync(request.app_js, "utf8");

/* ---- lifting code out of the shipped file, by name (same rule as the sequence recorder) -------- */

function closingBrace(open) {
    let depth = 0,
        quote = null;
    for (let k = open; k < APP_SOURCE.length; k++) {
        const ch = APP_SOURCE[k];
        if (quote) {
            if (ch === "\\") k++;
            else if (ch === quote) quote = null;
        } else if (ch === '"' || ch === "'" || ch === "`") quote = ch;
        else if (ch === "{") depth++;
        else if (ch === "}" && --depth === 0) return k;
    }
    throw new Error("unbalanced braces");
}

/** `    function NAME(` … its matching brace. Four spaces: the body of `window.nunaBoot`. */
function liftFunction(name) {
    const at = APP_SOURCE.indexOf(`\n    function ${name}(`);
    if (at < 0) throw new Error(`function ${name} not found in ${request.app_js}`);
    return APP_SOURCE.slice(at + 1, closingBrace(APP_SOURCE.indexOf("{", at)) + 1);
}

// `utf8` is what `decodeNseq` reads the magic and the JSON header with.
const LIFTED = [liftFunction("utf8"), liftFunction("decodeNseq")].join("\n");
const decodeNseq = new Function(`${LIFTED}\nreturn decodeNseq;`)();

/* ---- the canonical rendering — IDENTICAL to the Python side's ---------------------------------- */
/* One line per gene: start,end,contig,strand,partial5,phase. The flag BITS are read here exactly as
   `app.js` reads them elsewhere in the file, which is the one thing this comparison is about. */
const FLAG_MINUS = 1 << 0;
const FLAG_PARTIAL5 = 1 << 1;

function renderGenes(decoded) {
    const n = decoded.head.n_genes;
    const g = decoded.genes;
    const lines = new Array(n);
    const counters = {
        genes: n,
        minus_strand: 0,
        plus_strand: 0,
        five_prime_partial: 0,
        phase: [0, 0, 0, 0],
        flag_bits_outside_the_defined_four: 0,
    };
    for (let j = 0; j < n; j++) {
        const start = g[4 * j],
            end = g[4 * j + 1],
            contig = g[4 * j + 2],
            flags = g[4 * j + 3];
        const minus = (flags & FLAG_MINUS) !== 0;
        const partial5 = (flags & FLAG_PARTIAL5) !== 0;
        const phase = (flags >> 2) & 3;
        if (minus) counters.minus_strand++;
        else counters.plus_strand++;
        if (partial5) counters.five_prime_partial++;
        counters.phase[phase]++;
        if (flags & ~0xf) counters.flag_bits_outside_the_defined_four++;
        lines[j] = `${start},${end},${contig},${minus ? "-" : "+"},${partial5 ? 1 : 0},${phase}`;
    }
    return { lines, counters };
}

function digest(text) {
    return crypto.createHash("sha256").update(text, "utf8").digest("hex");
}

function fileArrayBuffer(path) {
    // ⚠ `new Uint8Array(buffer)` COPIES into a fresh ArrayBuffer of exactly the file's length.
    // Passing a node Buffer's own `.buffer` would hand `decodeNseq` the whole pooled slab.
    return new Uint8Array(fs.readFileSync(path)).buffer;
}

const out = [];
for (const genome of request.genomes) {
    const decoded = decodeNseq(fileArrayBuffer(genome.nseq));
    const { lines, counters } = renderGenes(decoded);
    const row = {
        sample: genome.sample,
        header_sample: decoded.head.sample,
        n_genes: decoded.head.n_genes,
        n_bases: decoded.head.n_bases,
        contig_count: decoded.head.cname.length,
        contig_names_head: decoded.head.cname.slice(0, 3),
        ambiguous_run_count: (decoded.head.ambiguous || []).length,
        /* ⛔ Digests over the WHOLE tables, not a head and a count. A reversed `clen` moves every
           contig base offset, and a wrong `ambiguous` pair paints the N back in the wrong place —
           both of which a first-3 sample and a length agree on perfectly. `coff` is included even
           though it derives from `clen`, because it is what `contigSpan` actually indexes with. */
        contig_names_digest: digest(decoded.head.cname.join("\n")),
        contig_lengths_digest: digest(decoded.head.clen.join("\n")),
        contig_offsets_digest: digest(Array.prototype.join.call(decoded.coff, "\n")),
        ambiguous_digest: digest(
            (decoded.head.ambiguous || []).map(function (run) { return run[0] + "," + run[1]; }).join("\n")
        ),
        genes_digest: digest(lines.join("\n")),
        genes_head: lines.slice(0, 3),
        counters: counters,
    };
    if (genome.loci) {
        const loci = new Int32Array(fileArrayBuffer(genome.loci));
        const text = Array.prototype.join.call(loci, "\n");
        row.loci_length = loci.length;
        row.loci_digest = digest(text);
        row.loci_head = Array.prototype.slice.call(loci, 0, 3);
        let unassigned = 0;
        for (let k = 0; k < loci.length; k++) if (loci[k] === -1) unassigned++;
        row.loci_unassigned = unassigned;
    }
    out.push(row);
}
process.stdout.write(JSON.stringify({ genomes: out }));
