"""Where the database is KNOWN to differ from the frozen page, each one named and explained.

⛔ **This file exists so that "parity" never becomes a tolerance.** A suite that passes with a
percentage bound has stopped being a parity test: it cannot distinguish two loci that moved for a
recorded reason from two hundred that moved because an ingest is wrong. So every exception here is
an explicit set of locus labels with the reason, and a suite that finds a difference OUTSIDE these
sets fails.

⚠ Each entry also asserts its own size. If an exception ever covers more loci than it was written
for, that is a new difference wearing an old label, and it must fail rather than be absorbed.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ParityException:
    """One recorded difference between the frozen page and what the current artifacts produce."""

    species_key: str
    column: str
    node_labels: frozenset[str]
    frozen_value: str
    current_value: str
    reason: str


# ⛔ **`AUDIT_TIER_RETIREMENT` and `AUDIT_TIER_RETIREMENT_KP` were DELETED on 2026-10-03, and that is
# the result they were written to produce.** They recorded that the published pages (2026-08-25, git
# `41e94f4`) predated the retirement of the tier name `no_homology`, so 2 ecoli and 6 kp loci read
# `no_homology` on the page and `synteny_only` in the database. Both pages were **re-exported at
# schema 16 on 2026-10-03** from the current artifacts, so the page and the database now agree on the
# tier of every locus — 17,531 ecoli and 15,670 kp, with non-null counts matching at 12,104 / 11,272
# so the agreement is not a None-vs-string comparison quietly passing.
#
# ⚠ **They were deleted, not narrowed to `frozenset()`.** An empty exception leaves its test asserting
# `set() == frozenset()`, which passes forever and gates nothing — precisely the "tolerance" this
# module's docstring exists to forbid. The tests that consulted them now assert the POSITIVE
# invariant instead: the tier agrees on every locus, over a stated count, and `RETIRED_TIER` appears
# in neither the string pool nor the database. `RETIRED_TIER` below is kept as that needle.

# ── the allele-variant symbol fold ──────────────────────────────────────────────────
@dataclass(frozen=True)
class SymbolFoldException:
    """The loci an allele-variant fold renames, and the payload blocks that carries.

    ⭐ **Why the ingest folds at all.** A locus is named for its commonest Bakta symbol counted as an
    EXACT STRING, and a substitution-tagged symbol is a different string from its own gene: kp 3878
    read `ompC` 28, `ompK36_T333N` 26, `ompK36` 10 and was named **`ompC` on a plurality of 28 while
    36 genes said OmpK36**. `ingest.allele_variant_symbols` folds the tag away before the count.

    ⚠ **Four loci in 33,201, and they are the only four**: every tagged symbol in both catalogues.
    Verified against the pre-fold database on 2026-09-24 — `display_name`, `display_name_source` and
    `best_product` over all 33,201 loci, four rows differing and no `best_product` moving at all.
    """

    species_key: str
    #: `(node_label, what the frozen page shows, what the fold makes it)` — values, never a count.
    renamed: tuple[tuple[str, str, str], ...]
    #: ⚠ Measured per species and NOT the same set. See `SYMBOL_FOLD_PAYLOAD_BLOCKS` below.
    payload_blocks: frozenset[str]
    reason: str

    @property
    def node_labels(self) -> frozenset[str]:
        return frozenset(label for label, _, _ in self.renamed)


_FOLD_REASON = (
    "An allele-variant symbol (`ompK36_T333N`) is folded onto the gene it is an allele of before "
    "the name vote, so one gene cannot split its own vote. The published page predates the fold."
)

#: ⛔ **`tufA` and `rpsJ` are STRAIGHT suffix strips** — the locus keeps the gene it always named and
#: loses a substitution tag that was never a locus-level fact. **`lon` is the same**, and additionally
#: merges `lon_P403L` 92 with `lon` 8 into one row of 100. **kp 3878 is the only one that changes
#: which GENE is named**, and everything else about that locus already said OmpK36: its modal product
#: (`OmpK36`, 34), its `best_product` (*outer membrane porin OmpK36*) and both major UniRef50
#: families' modal symbols.
SYMBOL_FOLD_ECOLI = SymbolFoldException(
    species_key="ecoli",
    renamed=(("2740", "lon_P403L", "lon"), ("17511", "tufA_E379K", "tufA")),
    # ⚠ **`nodes.name` is ABSENT here and present for kp**, and that is measured, not an oversight:
    # both *E. coli* loci keep the same pool INDEX while the string at that index changes, so the
    # index array is byte-identical and the decoded name is not. A suite asserting `nodes.name`
    # identical for ecoli would pass and mean nothing — which is why the decoded check exists.
    payload_blocks=frozenset({"lists.sym", "lists.u50", "strings.sym"}),
    reason=_FOLD_REASON,
)

SYMBOL_FOLD_KP = SymbolFoldException(
    species_key="kp",
    renamed=(("2512", "rpsJ_V57L", "rpsJ"), ("3878", "ompC", "ompK36")),
    payload_blocks=frozenset({"lists.sym", "lists.u50", "strings.sym", "nodes.name"}),
    reason=_FOLD_REASON,
)

SYMBOL_FOLDS: dict[str, SymbolFoldException] = {
    "ecoli": SYMBOL_FOLD_ECOLI,
    "kp": SYMBOL_FOLD_KP,
}

#: ⚠ **The raw index arrays differ far more widely than the four loci, and that is arithmetic.**
#: The `sym` string pool loses the tagged names (5,261 → 5,260 ecoli; 3,942 → 3,941 kp), so every
#: later index shifts by one; and ecoli 2740's symbol list loses a ROW (`lon_P403L` + `lon` become
#: one `lon`), which shifts every later element of the flat CSR arrays. Measured: 5,823 of 8,799
#: `lists.sym.idx` elements "differ" in *E. coli* while **two** loci actually changed.
#: ⛔ So the suite compares the symbol blocks DECODED, never as indices. A byte comparison here
#: could only ever be a tolerance.
SYMBOL_FOLD_IS_INDEX_SHIFT_NOT_CONTENT = (
    "`strings.sym` shrinks by one and one CSR row is removed, so thousands of indices move while "
    "the decoded content moves on four loci. Compare decoded."
)


def symbol_fold_for(species_key: str) -> SymbolFoldException:
    """The fold exception for one species. ⚠ Both species have one; there is no `None` case."""
    return SYMBOL_FOLDS[species_key]


#: The same fold, as the `bakta_gene_symbol` column — so the interned-string comparison consults one
#: registry rather than growing a second way to be excused. ⚠ `frozen_value`/`current_value` name the
#: RULE here because the four loci carry four different pairs; `SYMBOL_FOLDS` holds those pairs and
#: `test_T5_the_symbol_fold_renames_EXACTLY_the_loci_it_names` asserts every one of them.
SYMBOL_FOLD_COLUMN_ECOLI = ParityException(
    species_key="ecoli",
    column="bakta_gene_symbol",
    node_labels=SYMBOL_FOLD_ECOLI.node_labels,
    frozen_value="the allele-tagged symbol",
    current_value="the gene it is an allele of",
    reason=_FOLD_REASON,
)

SYMBOL_FOLD_COLUMN_KP = ParityException(
    species_key="kp",
    column="bakta_gene_symbol",
    node_labels=SYMBOL_FOLD_KP.node_labels,
    frozen_value="the allele-tagged symbol",
    current_value="the gene it is an allele of",
    reason=_FOLD_REASON,
)

#: Everything, indexed for a suite to consult.
KNOWN_PARITY_EXCEPTIONS: tuple[ParityException, ...] = (
    SYMBOL_FOLD_COLUMN_ECOLI,
    SYMBOL_FOLD_COLUMN_KP,
)


def exceptions_for(species_key: str, column: str) -> ParityException | None:
    """The recorded exception for one (species, column), or `None` if there is none."""
    for exception in KNOWN_PARITY_EXCEPTIONS:
        if exception.species_key == species_key and exception.column == column:
            return exception
    return None


#: ⭐⭐ **The same retirement, seen from the PAYLOAD.** Rebuilding the published catalogue from the
#: database reproduces 114 columns and 3.84 M elements byte-for-byte on both species; what does not
#: reproduce is four blocks, and all four have this one cause. Recorded as an explicit set so the
#: reproduction suite asserts *these four and no others* rather than counting to four.
#:
#: ⛔ **Two of the four went on 2026-10-03 and two did NOT, and the two that stayed are different in
#: kind from the two that left.**
#:
#: Gone, because the re-export at schema 16 made the page and the database agree:
#: * ~~`nodes.tier`~~   — the 2 (ecoli) / 6 (kp) loci now read `synteny_only` on both sides.
#: * ~~`strings.tier`~~ — `no_homology` was interned LAST, so dropping it shifted no other index.
#:   That was luck, not design, which is why the suite checks the pool rather than assuming it.
#:
#: Still here, and **neither is a re-export away** — both are things the DATABASE cannot say:
#: * `meta.audit` — ⛔ **not the six headline keys any more; they stopped moving.** What differs now
#:   is `sources`, and it is a defect on the database side, not the page's: `payload_serialiser`
#:   HARD-CODES the list as `[{label}_pfam_concordance.tsv, {label}_audit_summary.json]`, and no
#:   column anywhere on `pangenome` records which files the audit numbers came from. The 2026-08-25
#:   export named those two because `{label}_cluster_table.parquet` did not exist yet (it is dated
#:   2026-09-04); the 2026-10-03 export read the parquet and names three — correctly, since the
#:   exporter builds the list from files it actually opened. **The ingest read the same three**
#:   (`catalogue_frames.load_audit_evidence` prefers the parquet), so the literal under-reports the
#:   ingest's own inputs. ⚠ Until a column stores them, this is a standing gap in the claim that the
#:   database is a lossless superset of the catalogue — recorded here so it is visible rather than
#:   papered over by widening the literal to three names, which would hide it again until the audit
#:   directory next changes. **A schema decision for David, not one to take here.**
#: * `meta.omitted` — the audit ran with `--skip-seqid-to-medoid`, which the database records and the
#:   payload still has no way to say. (`seqid_coverage` is NULL on all 17,531 ecoli loci, measured.)
AUDIT_RERUN_PAYLOAD_BLOCKS: frozenset[str] = frozenset({"meta.audit", "meta.omitted"})

#: ⛔ **The medoid geometry's retirement, 2026-09-24 — a THIRD recorded cause, and a temporary one.**
#:
#: The published site catalogue is FROZEN at schema 14 and the database is now at 16, so three blocks
#: differ for a reason that is not a defect: the payload carries `map_reps` and `null` that the
#: database no longer has anything to build, and the database carries `sim` that the payload has no
#: key for. `nodes` differs too, on the four `*_d_*` columns.
#:
#: ✅ **EMPTIED 2026-10-03, which is what it was written to do.** Both species were re-exported at
#: schema 16 and republished, so the payload and the database describe the same world: `map_reps` and
#: `null` are gone from both, the four `*_d_*` columns are gone from both, and `sim` is present in
#: both — **byte-identical**, 30 columns over 175,334 ecoli / 156,724 kp elements.
#:
#: ⛔ **Kept as an empty frozenset rather than deleted, and the difference matters.** The reproduction
#: suite subtracts this set from what it found; an empty set subtracts nothing, so any one of these
#: seven lines coming back is a FAILURE naming itself. Deleting the constant would delete that guard
#: along with it. The day one of these blocks reappears is the day a payload and a database have
#: drifted apart again, and it should look exactly like it did on 2026-09-24.
#:
#: ⚠ The entries were difference LINES verbatim, not block names: a block absent on one side is
#: reported as added or removed rather than as changed.
#: ⛔ `schema` was NOT among them, and that sentence used to end "because `INGESTED_PAYLOAD_SCHEMA`
#: pins the rebuild to the schema its rows came from" — which was the excuse that hid a real defect.
#: The pin had said 14 since the 2026-09-25 re-ingest made the rows schema 16, so the rebuild stamped
#: a schema-16-shaped payload `"schema": 14` and the frozen schema-14 oracle agreed with it. Two
#: stale numbers in the same direction. Re-exporting the oracle separated them and `schema` failed
#: immediately; see `payload_serialiser.INGESTED_PAYLOAD_SCHEMA`.
MEDOID_RETIREMENT_PAYLOAD_BLOCKS: frozenset[str] = frozenset()

#: ✅ **EMPTIED 2026-10-03.** These six were arithmetic consequences of the 2/6 loci that moved tier,
#: so they went when the tier difference did: all 21 headline keys now hold equal values on both
#: species. Kept empty for the same reason as `MEDOID_RETIREMENT_PAYLOAD_BLOCKS` — the suite asserts
#: the moved set EQUALS this, so a key that starts moving again fails by name.
#:
#: ⚠ The three `no_homology_*` keys are still EMITTED, at 0, even though the tier name itself is gone
#: from the pool. A reader looking for the retirement in the headline will not find it there.
AUDIT_RERUN_HEADLINE_KEYS: frozenset[str] = frozenset()

#: The tier name that no longer exists, and therefore no longer enters the pool.
RETIRED_TIER = "no_homology"


# ── T6 · the Sequence tab ───────────────────────────────────────────────────────────────────────
@dataclass(frozen=True)
class DisplayedValueException:
    """A value the frozen page PRINTED differently from the new client, where both are right.

    Not a data difference — the API carries the page's exact value — so it is keyed by the rule that
    relates the two, and the suite asserts that rule on every gene it compares rather than counting
    to a number.
    """

    species_keys: tuple[str, ...]
    field: str
    frozen_value: str
    current_value: str
    rule: str
    reason: str


#: ⚠ **The Sequence tab's "Contig" row, and the contig named in its short-flank sentence.** The frozen
#: page printed `cname[contig_index]` from the `.nseq` header — the GFF's full seqid,
#: `SAMEA103923484.contig00324` — and `GeneSequenceCard.vue` prints `contig_name`, `contig00324`.
#:
#: ⛔ **The page never printed `contig_index + 1`.** `.nseq` shipped the names in `contig_idx` order
#: precisely so it would not have to (`genome_sequence.schema.md`, trap 1), and T6 measured the
#: page's contig equal to the API's `seqid` on **every one of 1,021,997 genes** (both species,
#: `BACATLAS_SEQUENCE_PARITY_EVERY_GENE=1`, 2026-09-18). The 28.6 % `index + 1` divergence is real in
#: the database but was never on the page.
#:
#: Measured on the dev database (`syntitude_dev` when measured, `bacatlas_dev` since 2026-09-26 —
#: renamed in place, so the measurement stands): `seqid == f"{sample_id}.{contig_name}"` on all 26,878 contigs, which is
#: the rule below. The API serves both fields, so this is the card's choice of which to print — a
#: difference a reader sees, recorded here rather than absorbed, and a decision for the owner.
CONTIG_ROW_SHOWS_NAME_NOT_SEQID = DisplayedValueException(
    species_keys=("ecoli", "kp"),
    field="Contig",
    frozen_value="SAMEA103923484.contig00324 — the GFF seqid, the `.nseq` header's `cname`",
    current_value="contig00324 — `genome_contig.contig_name`",
    rule='seqid == f"{sample_id}.{contig_name}" on every gene compared',
    reason=(
        "The published page printed the seqid; the new card prints the part after the sample prefix. Both name "
        "the same contig on every gene, and the API returns both."
    ),
)
