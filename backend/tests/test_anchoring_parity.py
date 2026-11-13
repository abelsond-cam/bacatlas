"""Parity suite T4 — ANCHORING: the frozen page's own answers against the new stack's, both species.

**The "before" side is the page itself.** `record_frozen_page_anchoring.js` boots the published
`app.js` under node, exactly as nuna's T2 recorder does, and asks the page's own `anchorRanks`,
`membershipComplete` and `genomeCounts` — under the page's own `setAnchor` — for every genome at every
locus: 1.75 M (ecoli) and 1.57 M (kp) grid cells in under a second. It is run at test time, so it
always grades against the `app.js` and payload that exist, and it SKIPS with the reason where node, the
nuna checkout or the payload is absent. (A committed fixture of the same answers would be ~4 MB per
species, above the ~3.2 MB the T2 fixtures set as the ceiling, and would go stale silently.)

**The "after" side is the new stack's own code path**, never a re-derivation here:

* `anchor_arrangement_ranks` — the function `load_locus_detail` calls — for EVERY (locus, genome) pair
  where the genome has a gene, over the arrangements loaded once per species;
* `load_listed_arrangements` — the statement the locus view issues — for EVERY pair whose anchored rank
  lies past the display cap, which is where the OR that offers it is the whole rule;
* `membership_is_complete` — the function the serialiser calls — at every locus;
* `GET /species/{k}/genomes` for every genome's counts, and the ingest's own count statement;
* and the Flask test client for a sample built to contain every category the rule has to get right.

⚠ **The rank convention, stated once.** The page's `anchorRanks(i)` returns `k`, a 0-based index into
the locus's arrangement list as exported in rank order, and prints `#k+1`; the API's `rank` is
`rank_within_locus`, also 0-based. They are the same integer only because every locus's ranks are
exactly `0..n-1` in both — which the first test asserts rather than assumes.

⛔ **Every test asserts its coverage before it reports a difference** (`4ab35ca`): how many pairs of the
grid it compared, and how many of them were ρ > 1, past the cap, or at an incomplete locus.
"""

from __future__ import annotations

import importlib.util
import json
import os
import random
import shutil
import subprocess
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import pytest
from sqlalchemy import create_engine, select, text
from sqlalchemy.orm import Session

from bacatlas_backend.application_factory import create_application
from bacatlas_backend.configuration import Configuration
from bacatlas_backend.ingest.ingest_genome_locus_counts import GENOME_LOCUS_COUNT_SELECT
from bacatlas_backend.ingest.published_catalogues import catalogue as published_catalogue
from bacatlas_backend.models.genome import Genome
from bacatlas_backend.models.genome_collection import GenomeCollectionMembership
from bacatlas_backend.models.locus import Locus
from bacatlas_backend.models.locus_arrangement import LocusArrangement
from bacatlas_backend.models.pangenome import Pangenome
from bacatlas_backend.services.locus_detail_service import (
    ARRANGEMENT_PAGE_SIZE,
    anchor_arrangement_ranks,
    load_listed_arrangements,
    membership_is_complete,
)
from tests.conftest import PUBLISHED_SITE_CATALOGUE_DIR
from tests.payload_oracle import load_catalogue

SPECIES_KEYS = ("ecoli", "kp")

RECORDER = Path(__file__).with_name("record_frozen_page_anchoring.js")

#: The display cap the locus view applies — read from the service, never restated.
DISPLAY_CAP = ARRANGEMENT_PAGE_SIZE

#: The HTTP/DOM sample, per category. Small enough to run in seconds; every category is required
#: to be non-empty, so a catalogue that stopped exercising one fails rather than passing thinly.
SAMPLE_SIZES = {
    "past_the_cap": 40,
    "rho_straddling_the_cap": 10,
    "rho_above_one": 30,
    "gene_with_no_window": 30,
    "no_gene_at_a_complete_locus": 20,
    "no_gene_at_an_incomplete_locus": 20,
    "ordinary": 30,
}


def _nuna_root() -> Path | None:
    """The nuna checkout the frozen page lives in, found through the importable package."""
    spec = importlib.util.find_spec("nuna")
    if spec is None or not spec.origin:
        return None
    return Path(spec.origin).resolve().parents[2]


@dataclass
class AnchoringParity:
    """One species: the page's recorded answers beside the loaded database, index-aligned."""

    species_key: str
    engine: object
    pangenome: Pangenome
    catalogue: object
    loci: list
    genome_id_by_ordinal: list[int]
    sample_id_by_ordinal: list[str]
    arrangements_by_ordinal: list[list]
    #: ⚠ From `locus.present_genome_ids`: every genome with ≥ 1 gene at the locus. The payload cannot
    #: name these — `arr.gid` omits a genome whose gene reached no window — so the pairs where the
    #: two differ are exactly the "no recorded neighbourhood" cases.
    present_by_ordinal: list[set[int]]
    recording: dict
    #: `(locus ordinal, genome ordinal) → ranks`, non-empty only — the page's sparse answer.
    page_ranks: dict = field(default_factory=dict)
    sample: dict = field(default_factory=dict)

    @property
    def ordinal_by_genome_id(self) -> dict[int, int]:
        return {genome_id: k for k, genome_id in enumerate(self.genome_id_by_ordinal)}

    def page_complete(self, locus_ordinal: int) -> bool:
        return self.recording["membership_complete"][locus_ordinal] == "1"


def _decode_page_ranks(recording: dict) -> dict:
    out = {}
    for k, flat in enumerate(recording["anchor_ranks_by_genome"]):
        position = 0
        while position < len(flat):
            locus, count = flat[position], flat[position + 1]
            out[(locus, k)] = tuple(flat[position + 2 : position + 2 + count])
            position += 2 + count
    return out


def _choose_sample(parity: AnchoringParity, seed: int) -> dict[tuple[int, int], set[str]]:
    """Pairs for the HTTP and DOM checks, drawn per category with a named seed.

    ⭐ Deliberately not uniform: a uniform draw of 200 from ~490k pairs would contain ~6 past the cap
    and ~0 at ρ > 1 — exactly the cases this suite exists for. The deepest past-cap rank and the widest
    ρ > 1 pair are always included, so the extremes are never left to the draw.
    """
    rng = random.Random(seed)
    ranks = parity.page_ranks
    n_loci, n_genomes = len(parity.loci), len(parity.genome_id_by_ordinal)
    by_category: dict[str, list[tuple[int, int]]] = {name: [] for name in SAMPLE_SIZES}

    ordered = sorted(ranks)
    past = [pair for pair in ordered if max(ranks[pair]) >= DISPLAY_CAP]
    rho = [pair for pair in ordered if len(ranks[pair]) > 1]
    straddling = [pair for pair in rho if min(ranks[pair]) < DISPLAY_CAP <= max(ranks[pair])]
    ordinary = [pair for pair in ordered if len(ranks[pair]) == 1 and ranks[pair][0] < DISPLAY_CAP]
    no_window = sorted(
        (i, parity.ordinal_by_genome_id[genome_id])
        for i, present in enumerate(parity.present_by_ordinal)
        for genome_id in present
        if (i, parity.ordinal_by_genome_id[genome_id]) not in ranks
    )

    def draw(pool, size, always=()):
        chosen = list(dict.fromkeys(always))
        rest = [pair for pair in pool if pair not in chosen]
        chosen += rng.sample(rest, min(len(rest), max(0, size - len(chosen))))
        return chosen

    by_category["past_the_cap"] = draw(
        past, SAMPLE_SIZES["past_the_cap"], [max(past, key=lambda pair: max(ranks[pair]))] if past else []
    )
    by_category["rho_straddling_the_cap"] = draw(straddling, SAMPLE_SIZES["rho_straddling_the_cap"])
    by_category["rho_above_one"] = draw(
        rho, SAMPLE_SIZES["rho_above_one"], [max(rho, key=lambda pair: len(ranks[pair]))] if rho else []
    )
    by_category["gene_with_no_window"] = draw(no_window, SAMPLE_SIZES["gene_with_no_window"])
    by_category["ordinary"] = draw(ordinary, SAMPLE_SIZES["ordinary"])

    for name, wanted_complete in (
        ("no_gene_at_a_complete_locus", True),
        ("no_gene_at_an_incomplete_locus", False),
    ):
        loci = [i for i in range(n_loci) if parity.page_complete(i) == wanted_complete]
        tries = 0
        while len(by_category[name]) < SAMPLE_SIZES[name] and tries < 100_000:
            tries += 1
            i, k = rng.choice(loci), rng.randrange(n_genomes)
            if parity.genome_id_by_ordinal[k] not in parity.present_by_ordinal[i]:
                if (i, k) not in by_category[name]:
                    by_category[name].append((i, k))

    tagged: dict[tuple[int, int], set[str]] = {}
    for name, pairs in by_category.items():
        for pair in pairs:
            tagged.setdefault(pair, set()).add(name)
    return tagged


@pytest.fixture(scope="module", params=SPECIES_KEYS)
def parity(request, tmp_path_factory) -> AnchoringParity:
    species_key = request.param
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not on PATH — T4's 'before' side is the frozen page, run under node")
    nuna_root = _nuna_root()
    if nuna_root is None:
        pytest.skip("the nuna package is not importable — the frozen app.js lives in that checkout")
    app_js = nuna_root / "src" / "nuna" / "tl" / "locus_browser" / "app.js"
    dom_shim = nuna_root / "tests" / "js" / "dom_shim.js"
    for path in (app_js, dom_shim):
        if not path.exists():
            pytest.skip(f"the frozen page's {path.name} is not present at {path}")
    payload = PUBLISHED_SITE_CATALOGUE_DIR / f"{species_key}.json"
    if not payload.exists():
        pytest.skip(f"the published site catalogue is not present at {payload}")
    url = os.environ.get("BACATLAS_DATABASE_URL")
    if not url:
        pytest.skip("BACATLAS_DATABASE_URL is not set — T4 compares the page with the loaded database")

    entry = published_catalogue(species_key)
    engine = create_engine(url, future=True)
    with Session(engine) as session:
        pangenome = session.execute(
            select(Pangenome).where(Pangenome.run_id == entry.run_id)
        ).scalar_one_or_none()
        if pangenome is None:
            pytest.skip(f"run {entry.run_id} is not loaded in {url}")
        session.expunge(pangenome)
        loci = list(
            session.execute(
                select(Locus)
                .where(Locus.pangenome_id == pangenome.pangenome_id)
                .order_by(Locus.catalogue_ordinal)
            ).scalars()
        )
        roster = session.execute(
            select(GenomeCollectionMembership.collection_genome_ordinal, Genome.genome_id, Genome.sample_id)
            .join(Genome, Genome.genome_id == GenomeCollectionMembership.genome_id)
            .where(GenomeCollectionMembership.genome_collection_id == pangenome.genome_collection_id)
            .order_by(GenomeCollectionMembership.collection_genome_ordinal)
        ).all()
        ordinal_by_locus_id = {locus.locus_id: locus.catalogue_ordinal for locus in loci}
        arrangements_by_ordinal: list[list] = [[] for _ in loci]
        for arrangement in session.execute(
            select(LocusArrangement)
            .where(LocusArrangement.pangenome_id == pangenome.pangenome_id)
            .order_by(LocusArrangement.locus_id, LocusArrangement.rank_within_locus)
        ).scalars():
            arrangements_by_ordinal[ordinal_by_locus_id[arrangement.locus_id]].append(arrangement)
        present_by_ordinal: list[set[int]] = [set() for _ in loci]
        for locus in loci:
            present_by_ordinal[locus.catalogue_ordinal] = set(locus.present_genome_ids)
        session.expunge_all()

    catalogue = load_catalogue(payload)
    # ⛔ Coverage of the ALIGNMENT, before anything is compared: same loci in the same order, the
    # same genome vocabulary at the same ordinals. Every comparison below indexes by both.
    assert [locus.node_label for locus in loci] == [str(label) for label in catalogue.nodes["label"]]
    assert [ordinal for ordinal, _, _ in roster] == list(range(len(roster)))
    assert [sample for _, _, sample in roster] == list(catalogue.meta["genomes"])

    result = AnchoringParity(
        species_key=species_key,
        engine=engine,
        pangenome=pangenome,
        catalogue=catalogue,
        loci=loci,
        genome_id_by_ordinal=[genome_id for _, genome_id, _ in roster],
        sample_id_by_ordinal=[sample for _, _, sample in roster],
        arrangements_by_ordinal=arrangements_by_ordinal,
        present_by_ordinal=present_by_ordinal,
        recording={},
    )

    # The page's answers without the DOM first, so the sample can be drawn from them…
    grid = _record(node, app_js, dom_shim, payload, None)
    result.recording = grid
    result.page_ranks = _decode_page_ranks(grid)
    result.sample = _choose_sample(result, seed=20260918)

    # …then the sampled pairs and the dropdown queries driven through the rendered page.
    drive = tmp_path_factory.mktemp(f"anchoring_{species_key}") / "drive.json"
    pairs = sorted(result.sample, key=lambda pair: (pair[1], pair[0]))
    drive.write_text(json.dumps({"pairs": pairs, "queries": _queries(result.sample_id_by_ordinal)}))
    rendered = _record(node, app_js, dom_shim, payload, drive)
    # ⛔ The two boots must be the same page answering the same questions.
    assert rendered["anchor_ranks_by_genome"] == grid["anchor_ranks_by_genome"]
    result.recording = rendered
    return result


def _record(node, app_js, dom_shim, payload, drive) -> dict:
    command = [node, str(RECORDER), str(app_js), str(dom_shim), str(payload)]
    if drive is not None:
        command.append(str(drive))
    completed = subprocess.run(command, capture_output=True, text=True, timeout=600, check=False)
    # ⛔ A FAILURE, not a skip: node and the page are both here, so a recorder that cannot run means
    # the page moved under it (the injection point is checked to occur exactly once).
    assert completed.returncode == 0, completed.stderr[-2000:]
    return json.loads(completed.stdout)


def _queries(samples: list[str]) -> list[str]:
    """Dropdown queries: every shape `anchorSearch` distinguishes, drawn from the real accessions."""
    probe, other = samples[37], samples[len(samples) // 2]
    return [
        "",
        "   ",
        probe,
        probe.lower(),
        probe[3:9],
        f"  {probe[3:9].lower()}  ",
        other[-4:],
        probe[:3],
        "NOT-A-GENOME",
        "%",
        "_",
    ]


# ── T4.0 · the two sides describe the same thing ───────────────────────────────────────────────
def test_T4_the_page_and_the_database_share_loci_genomes_and_ONE_rank_space(parity):
    recording, catalogue = parity.recording, parity.catalogue
    n_loci, n_genomes = len(parity.loci), len(parity.genome_id_by_ordinal)
    assert recording["locus_count"] == n_loci == catalogue.n_loci
    assert recording["genome_count"] == n_genomes == 100
    assert recording["grid_pairs_evaluated"] == n_loci * n_genomes
    assert recording["genomes"] == parity.sample_id_by_ordinal
    # ⚠ The rank convention: the page's index into its list IS the API's rank, only because both are
    # exactly 0..n-1 over the same arrangements.
    bad = [
        parity.loci[i].node_label
        for i in range(n_loci)
        if [a.rank_within_locus for a in parity.arrangements_by_ordinal[i]]
        != list(range(len(catalogue.arrangements(i))))
    ]
    assert not bad, bad[:5]


# ── T4.1 · the rank LIST, every genome × every locus it is in ──────────────────────────────────
def test_T4_every_genome_carries_the_SAME_RANK_LIST_at_every_locus_it_is_in(parity):
    """⛔ A list: ρ > 1 puts one genome in two arrangements at one locus, up to 14 at once."""
    ordinal = parity.ordinal_by_genome_id
    genomes_per_locus = parity.catalogue.nodes["genomes"]

    # Coverage first: the database's present genomes are the payload's own per-locus genome count,
    # at every locus, and every pair the page anchors is one the database says is present.
    miscounted = [
        (parity.loci[i].node_label, len(present), genomes_per_locus[i])
        for i, present in enumerate(parity.present_by_ordinal)
        if len(present) != genomes_per_locus[i]
    ]
    assert not miscounted, miscounted[:5]
    outside = [
        pair
        for pair in parity.page_ranks
        if parity.genome_id_by_ordinal[pair[1]] not in parity.present_by_ordinal[pair[0]]
    ]
    assert not outside, f"the page anchors {len(outside)} pairs the database has no gene for"

    examined, multi, past_cap, no_window, differing = 0, 0, 0, 0, []
    for i, present in enumerate(parity.present_by_ordinal):
        arrangements = parity.arrangements_by_ordinal[i]
        for genome_id in present:
            k = ordinal[genome_id]
            ours = anchor_arrangement_ranks(arrangements, genome_id)
            theirs = list(parity.page_ranks.get((i, k), ()))
            examined += 1
            multi += len(theirs) > 1
            past_cap += bool(theirs) and max(theirs) >= DISPLAY_CAP
            no_window += not theirs
            if ours != theirs:
                differing.append((parity.loci[i].node_label, parity.sample_id_by_ordinal[k], ours, theirs))

    assert examined == sum(genomes_per_locus), f"compared {examined:,} (locus, genome) pairs"
    assert examined - no_window == len(parity.page_ranks) == parity.recording["recorded_pairs"]
    assert multi == parity.recording["multi_rank_pairs"] > 0, "no ρ > 1 pair — the LIST is untested"
    assert past_cap > 0 and no_window > 0
    assert not differing, differing[:5]


# ── T4.2 · membership completeness, every locus ────────────────────────────────────────────────
def test_T4_membership_completeness_selects_the_same_sentence_at_every_locus(parity):
    """⛔ *"has no gene at this locus"* is FALSE wherever a present genome reached no window."""
    assert len(parity.recording["membership_complete"]) == len(parity.loci)
    differing = [
        locus.node_label
        for i, locus in enumerate(parity.loci)
        if membership_is_complete(locus) != parity.page_complete(i)
    ]
    incomplete = [i for i in range(len(parity.loci)) if not parity.page_complete(i)]
    assert len(incomplete) > 0, "no incomplete locus — the second sentence is untested"
    assert not differing, differing[:5]

    # ⭐ And at each incomplete locus, the genomes the database holds present-but-unarranged are
    # exactly as many as the page knows are missing from its union — so the page's "incomplete" is
    # the same genomes, not merely the same verdict.
    arranged_on_page: dict[int, set[int]] = {}
    for i, k in parity.page_ranks:
        arranged_on_page.setdefault(i, set()).add(k)
    for i in incomplete:
        missing_on_page = parity.catalogue.nodes["genomes"][i] - len(arranged_on_page.get(i, ()))
        unarranged_here = [
            g for g in parity.present_by_ordinal[i]
            if parity.ordinal_by_genome_id[g] not in arranged_on_page.get(i, ())
        ]
        assert missing_on_page == len(unarranged_here) > 0, parity.loci[i].node_label


# ── T4.3 · the anchored arrangement is OFFERED past the display cap ────────────────────────────
def test_T4_every_anchored_arrangement_PAST_THE_CAP_is_offered_by_the_views_own_statement(parity):
    """⛔ `arrShown`: *"otherwise the reader is told their genome sits in #37 and has no button"*.

    Every pair, not a sample: `load_listed_arrangements` is the one statement the locus view uses to
    pick its rows, driven here once per anchored pair past the cap.
    """
    past = sorted(pair for pair, ranks in parity.page_ranks.items() if max(ranks) >= DISPLAY_CAP)
    independent = {
        (i, parity.ordinal_by_genome_id[genome_id])
        for i, arrangements in enumerate(parity.arrangements_by_ordinal)
        for arrangement in arrangements
        if arrangement.rank_within_locus >= DISPLAY_CAP
        for genome_id in arrangement.member_genome_ids
    }
    assert set(past) == independent and len(past) > 0, f"{len(past):,} past-cap pairs"

    differing = []
    with Session(parity.engine) as session:
        for i, k in past:
            genome_id = parity.genome_id_by_ordinal[k]
            listed = load_listed_arrangements(
                session,
                locus_id=parity.loci[i].locus_id,
                arrangement_limit=DISPLAY_CAP,
                anchor_genome_id=genome_id,
            )
            listed_ranks = [arrangement.rank_within_locus for arrangement in listed]
            expected = list(parity.page_ranks[(i, k)])
            commonest = list(range(min(DISPLAY_CAP, len(parity.arrangements_by_ordinal[i]))))
            if (
                anchor_arrangement_ranks(listed, genome_id) != expected
                or not set(expected) <= set(listed_ranks)
                or listed_ranks[: len(commonest)] != commonest
            ):
                differing.append((parity.loci[i].node_label, parity.sample_id_by_ordinal[k], listed_ranks, expected))
    assert not differing, differing[:5]


# ── T4.4 · the per-genome locus count ──────────────────────────────────────────────────────────
@pytest.fixture(scope="module")
def client():
    url = os.environ.get("BACATLAS_DATABASE_URL")
    if not url:
        pytest.skip("BACATLAS_DATABASE_URL is not set")
    return create_application(Configuration(database_url=url)).test_client()


def test_T4_each_genomes_locus_count_is_the_dropdowns_and_is_COUNT_DISTINCT(parity, client):
    """⛔ Per LOCUS, not per arrangement: a genome at ρ > 1 is counted once (`app.js:3925`)."""
    recording = parity.recording
    served = client.get(
        f"/api/v1/species/{parity.species_key}/genomes", query_string={"limit": 1000}
    ).get_json()["genomes"]
    assert [g["sample_id"] for g in served] == recording["genomes"], "the vocabularies differ"
    assert len(served) == len(recording["genome_counts_from_function"]) == 100

    # The hook's counts are the ones the dropdown PRINTS, so the page's own two readings agree first.
    assert recording["genome_counts_from_dropdown"] == recording["genome_counts_from_function"]
    page = recording["genome_counts_from_function"]

    # ⭐ The served number against the page's, genome for genome.
    assert [g["arrangement_locus_count"] for g in served] == page

    # ⭐ And the ingest's own statement, recomputed now, against the page — so the definition that
    # writes the next catalogue is graded, not only the rows the migration backfilled.
    with parity.engine.connect() as connection:
        recomputed = {
            genome_id: (present, arranged)
            for genome_id, present, arranged in connection.execute(
                text(GENOME_LOCUS_COUNT_SELECT), {"pangenome_id": parity.pangenome.pangenome_id}
            )
        }
    assert [recomputed[genome_id][1] for genome_id in parity.genome_id_by_ordinal] == page

    # ⛔ `locus_count` has no counterpart on the page — `arr.gid` cannot see a gene with no window —
    # so it is held to the payload's own identities instead. Summed over genomes it is the per-locus
    # genome count summed over loci; per genome it is the page's count plus the loci where the
    # genome is present and the page anchors it nowhere.
    genomes_per_locus = parity.catalogue.nodes["genomes"]
    assert sum(g["locus_count"] for g in served) == sum(genomes_per_locus)
    unanchored = Counter(
        k
        for i, present in enumerate(parity.present_by_ordinal)
        for genome_id in present
        for k in [parity.ordinal_by_genome_id[genome_id]]
        if (i, k) not in parity.page_ranks
    )
    expected_present = [page[k] + unanchored[k] for k in range(len(page))]
    assert [g["locus_count"] for g in served] == expected_present
    # …and the same of the ingest statement's own `locus_count`, for the same reason as above.
    assert [recomputed[genome_id][0] for genome_id in parity.genome_id_by_ordinal] == expected_present
    assert all(g["locus_count"] <= len(parity.loci) for g in served)
    assert sum(unanchored.values()) > 0, "no present-but-unanchored pair — the two counts would be one"


# ── T4.5 · the rendered page, and the wire, on a sample of every category ──────────────────────
def _dom_ranks(rendered: dict) -> list[int]:
    """What the rendered switcher SAYS the anchored genome carries.

    ⚠ A locus with one arrangement draws no buttons (`renderArrangements`, `tot <= 1`), only the line;
    there the line's "Anchored to …" is the page saying rank 0, the only one there is.
    """
    if rendered["option_ranks"]:
        return sorted(rendered["anchored_option_ranks"])
    return [0] if rendered["anchor_line"] and not rendered["anchor_line_is_muted"] else []


def test_T4_the_hook_reports_exactly_what_the_page_RENDERS(parity):
    """⭐ The recorder's hook is checked, not trusted: each sampled pair driven through the real control."""
    rendered = {(row["locus_index"], row["genome_ordinal"]): row for row in parity.recording["dom_pairs"]}
    assert set(rendered) == set(parity.sample), "the DOM drive did not cover the sample"
    assert not [row for row in rendered.values() if "error" in row]
    differing, silent = [], 0
    for (i, k), row in rendered.items():
        hooked = list(parity.page_ranks.get((i, k), ()))
        sentence_complete = None
        if not parity.catalogue.arrangements(i):
            # ⚠ A locus with NO arrangement at all — every member gene alone on its contig (622 ecoli,
            # 225 kp) — draws no switcher and therefore NO anchor line (`renderArrangements` returns
            # on an empty list). The page says nothing there, rather than either sentence; the API
            # still answers `[]` and `membership_is_complete: false`, and a client must not turn
            # that into a sentence the page never printed.
            silent += 1
            if hooked or row["option_ranks"] or row["anchor_line"] is not None:
                differing.append((parity.loci[i].node_label, parity.sample_id_by_ordinal[k], row, hooked))
            continue
        if not hooked:
            # The anchor line's two sentences, as the page wrote them.
            sentence_complete = "has no gene at this locus" in (row["anchor_line"] or "")
            assert row["anchor_line_is_muted"] is True
        if _dom_ranks(row) != hooked or (
            sentence_complete is not None and sentence_complete != parity.page_complete(i)
        ):
            differing.append((parity.loci[i].node_label, parity.sample_id_by_ordinal[k], row, hooked))
    assert silent < len(rendered), "every sampled locus was one the page draws nothing at"
    assert not differing, differing[:3]


def test_T4_the_WIRE_agrees_with_the_page_in_every_category(parity, client):
    """The locus response for each sampled (locus, anchored genome), against the page's answer."""
    counts = Counter(name for names in parity.sample.values() for name in names)
    assert set(counts) == set(SAMPLE_SIZES), f"categories never sampled: {set(SAMPLE_SIZES) - set(counts)}"
    assert all(counts[name] > 0 for name in SAMPLE_SIZES), counts

    differing = []
    for (i, k), names in sorted(parity.sample.items()):
        label, sample_id = parity.loci[i].node_label, parity.sample_id_by_ordinal[k]
        response = client.get(
            f"/api/v1/species/{parity.species_key}/loci/{label}", query_string={"anchor": sample_id}
        )
        assert response.status_code == 200, (label, sample_id, response.status_code)
        body = response.get_json()
        expected = list(parity.page_ranks.get((i, k), ()))
        listed = [row["rank"] for row in body["arrangements"]["listed"]]
        problems = []
        if body["anchor"] != {"is_anchored": True, "arrangement_ranks": expected}:
            problems.append(("anchor", body["anchor"], expected))
        if not set(expected) <= set(listed):
            problems.append(("not offered", listed, expected))
        if body["arrangements"]["membership_is_complete"] != parity.page_complete(i):
            problems.append(("completeness", body["arrangements"]["membership_is_complete"]))
        if body["arrangements"]["total"] != parity.catalogue.total_arrangement_count(i):
            problems.append(("total", body["arrangements"]["total"]))
        if "past_the_cap" in names and not any(rank >= DISPLAY_CAP for rank in listed):
            problems.append(("past-cap row absent", listed))
        if problems:
            differing.append((label, sample_id, sorted(names), problems))
    assert not differing, differing[:3]


# ── T4.6 · the picker's query means what `anchorSearch` meant ──────────────────────────────────
def test_T4_the_genome_query_lists_exactly_what_the_dropdown_listed(parity, client):
    queries = parity.recording["dom_queries"]
    assert len(queries) == len(_queries(parity.sample_id_by_ordinal))
    sizes = {len(entry["listed"]) for entry in queries}
    assert 0 in sizes and 100 in sizes and any(0 < size < 100 for size in sizes), sizes
    differing = []
    for entry in queries:
        served = client.get(
            f"/api/v1/species/{parity.species_key}/genomes",
            query_string={"q": entry["query"], "limit": 1000},
        ).get_json()
        listed = [genome["sample_id"] for genome in served["genomes"]]
        if listed != entry["listed"] or served["matched_genome_count"] != len(entry["listed"]):
            differing.append((entry["query"], len(listed), len(entry["listed"])))
    assert not differing, differing
