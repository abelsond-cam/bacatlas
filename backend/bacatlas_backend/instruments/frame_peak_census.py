"""frame_peak_census — what the ingest's two big in-memory frames actually cost, measured in a clean process.

⛔ **Why this is not `DataFrame.memory_usage(deep=True)`.** That counts a shared Python string once **per
cell**, so an `object` column of 372 M rows drawn from 71,435 distinct genome ids is reported at ~71,000×
its true cost. Every B/row figure in the 71k plan that came from it is an upper bound of the wrong order.
The kernel's high-water mark is not foolable that way.

**Two numbers, because they answer different questions** (both from nuna's `resource_log`, which the model
steps already use, so the figures join):

* ``VmHWM`` — peak resident set, what SLURM checks ``--mem`` against. ⚠ Monotone: it never falls, so it
  reports the peak *so far* and a later cheaper phase shows a delta of 0. This is the number that sizes a job.
* ``RssAnon`` — anonymous heap only, page cache excluded, and it **does** fall when a frame is freed. This is
  the number for *"what does this frame cost"* and for comparing two builds of it.

⭐ **The two frames are measured in ONE process and in the ingest's own order**, because that is how the
ingest holds them: `_annotation_frame` is built first and is still resident when `oriented_windows` runs.
Measuring either alone understates the peak that matters. ⛔ And one process measures ONE configuration —
`VmHWM` cannot be reset, so `categorical=True` and `categorical=False` are separate invocations.

What it answers, in the plan's terms: the measured B/row for *this* window frame (replacing a 215 B/row
inference from a differently-shaped frame), whether `categorical=True` is worth its ~1.5× build time at a
given cohort size, whether the per-genome fold is needed at all, and — since T2 was withdrawn —
`_annotation_frame`'s own peak, which is now the open term.

    python -m bacatlas_backend.instruments.frame_peak_census --data-root … --set-key kp --model-label … \
        --run-id … [--categorical] [--roster roster.txt] [--json out.json]
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Mark:
    """One named point in the run: the monotone peak, and the heap as it stands."""

    label: str
    hwm_gb: float
    anon_gb: float | None
    rows: int | None = None
    seconds: float | None = None
    note: str = ""


@dataclass
class Census:
    """Every mark, plus what the caller needs to reproduce it."""

    categorical: bool
    n_genomes: int
    n_genes: int
    n_loci: int
    marks: list[Mark] = field(default_factory=list)

    def add(self, label: str, *, rows: int | None = None, seconds: float | None = None, note: str = "") -> Mark:
        """Read both kernel counters NOW and file them under `label`. The read itself is two file reads."""
        from nuna.tl.cluster.resource_log import peak_rss_gb, rss_anon_gb

        mark = Mark(label, round(peak_rss_gb(), 3), None, rows, seconds, note)
        anon = rss_anon_gb()
        mark.anon_gb = None if anon is None else round(anon, 3)
        self.marks.append(mark)
        return mark

    def delta_anon(self, before: str, after: str) -> float | None:
        """Heap added between two marks — the frame's own cost, where RssAnon is available."""
        got = {m.label: m.anon_gb for m in self.marks}
        if got.get(before) is None or got.get(after) is None:
            return None
        return round(got[after] - got[before], 3)

    def to_dict(self) -> dict:
        """The census as JSON — so a number quoted from it names the file it was read from."""
        return {
            "categorical": self.categorical,
            "n_genomes": self.n_genomes,
            "n_genes": self.n_genes,
            "n_loci": self.n_loci,
            "marks": [vars(m) for m in self.marks],
        }


def _roster_filter(source: dict, roster: Path | None) -> dict:
    """Restrict every per-gene frame to a roster's genomes — the slice knob, applied ONCE and to all of them.

    ⛔ Filtering only some frames would make the measurement meaningless rather than wrong-looking: a
    `coords` of 1,000 genomes merged against an `assign` of 71,435 still yields a 1,000-genome window frame,
    and the peak would silently carry the big frame too.
    """
    if roster is None:
        return source
    wanted = {line.strip() for line in roster.read_text().splitlines() if line.strip()}
    missing = wanted - set(source["samples"])
    if missing:
        raise SystemExit(
            f"{len(missing):,} of {len(wanted):,} roster genomes are not in this catalogue's universe "
            f"(e.g. {sorted(missing)[:3]}) — the slice would measure a different cohort"
        )
    source["samples"] = [s for s in source["samples"] if s in wanted]
    for key in ("assign", "products", "coords", "noncoding", "dbxrefs", "pfam"):
        frame = source.get(key)
        if frame is not None and len(frame) and "sample_id" in frame.columns:
            source[key] = frame[frame["sample_id"].astype(str).isin(wanted)].reset_index(drop=True)
    source["n_genes"] = int(len(source["assign"]))
    return source


def run(artifacts, *, categorical: bool, roster: Path | None = None, allow_missing_strand: bool = False) -> Census:
    """Load the sources, build both frames in the ingest's order, and mark the kernel's counters throughout."""
    import time

    from nuna.eval.metrics.syntology import NEAR_SHELLS
    from nuna.tl.locus_browser.export_payload import oriented_windows

    from bacatlas_backend.ingest.catalogue_frames import _annotation_frame, _load_source_frames

    census = Census(categorical=categorical, n_genomes=0, n_genes=0, n_loci=0)
    census.add("start")

    t0 = time.perf_counter()
    source = _roster_filter(_load_source_frames(artifacts, allow_missing_strand=allow_missing_strand), roster)
    assign = source["assign"]
    census.n_genomes = len(source["samples"])
    census.n_genes = int(len(assign))
    census.n_loci = int(assign["node"].nunique())
    census.add("sources_loaded", rows=census.n_genes, seconds=round(time.perf_counter() - t0, 1))

    # ── the wide annotation frame, first and still resident, exactly as the ingest holds it ──────
    t0 = time.perf_counter()
    annotation = _annotation_frame(assign, source["products"], source["dbxrefs"], source["pfam"], source["numpy"])
    census.add(
        "annotation_frame", rows=int(len(annotation)), seconds=round(time.perf_counter() - t0, 1),
        note="T2 withdrawn 2026-10-30, so this frame stays per-gene and its peak is the open term",
    )  # fmt: skip

    # ── the ±5 neighbourhood: 2*NEAR_SHELLS rows per gene ────────────────────────────────────────
    nodes = sorted(assign["node"].astype(str).unique())
    node_position = {node: index for index, node in enumerate(nodes)}
    assign_coords = source["coords"].merge(assign, on=["sample_id", "flat_index"], how="inner")
    census.add("assign_coords", rows=int(len(assign_coords)))

    t0 = time.perf_counter()
    windows = oriented_windows(assign_coords, node_position, radius=NEAR_SHELLS, categorical=categorical)
    census.add(
        "oriented_windows", rows=int(len(windows)), seconds=round(time.perf_counter() - t0, 1),
        note=f"radius={NEAR_SHELLS}, categorical={categorical}",
    )  # fmt: skip
    del windows
    census.add("windows_freed", note="RssAnon here vs oriented_windows is the frame's own resident cost")
    return census


def report(census: Census) -> str:
    """One block a reader can act on: the peak that sizes `--mem`, and the B/row that sizes the frame."""
    lines = [
        f"=== frame peak census | categorical={census.categorical} | {census.n_genomes:,} genomes, "
        f"{census.n_genes:,} genes, {census.n_loci:,} loci ===",
        f"{'mark':<20} {'VmHWM GB':>10} {'RssAnon GB':>11} {'rows':>14} {'s':>7}  note",
    ]
    for m in census.marks:
        lines.append(
            f"{m.label:<20} {m.hwm_gb:>10.3f} {('-' if m.anon_gb is None else f'{m.anon_gb:.3f}'):>11} "
            f"{('-' if m.rows is None else f'{m.rows:,}'):>14} {('-' if m.seconds is None else f'{m.seconds:.1f}'):>7}"
            f"  {m.note}"
        )
    rows = {m.label: m.rows for m in census.marks}
    for frame, before, after in (
        ("annotation_frame", "sources_loaded", "annotation_frame"),
        ("oriented_windows", "assign_coords", "oriented_windows"),
    ):
        delta, n = census.delta_anon(before, after), rows.get(frame)
        if delta is None:
            lines.append(f"  {frame}: RssAnon unavailable (not Linux) — B/row not measurable here")
        elif not n:
            lines.append(f"  {frame}: no rows")
        else:
            lines.append(f"  {frame}: +{delta:.3f} GB heap over {n:,} rows = {delta * (1 << 30) / n:.1f} B/row")
    freed = census.delta_anon("oriented_windows", "windows_freed")
    if freed is not None:
        lines.append(f"  freeing the window frame returned {-freed:.3f} GB")
    peak = max((m.hwm_gb for m in census.marks), default=0.0)
    lines.append(f"  ⇒ VmHWM over the whole run: {peak:.3f} GB — this is what --mem is checked against")
    return "\n".join(lines)


def main() -> None:
    """CLI entry: ONE configuration, in a clean process, because VmHWM cannot be reset."""
    from bacatlas_backend.ingest.artifact_locator import CatalogueArtifacts

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root", type=Path, required=True)
    ap.add_argument("--set-key", required=True)
    ap.add_argument("--model-label", required=True)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--species-key", default=None)
    ap.add_argument("--gene-universe", type=Path, default=None)
    ap.add_argument("--assignment", type=Path, default=None)
    ap.add_argument("--categorical", action="store_true", help="build the window frame with categorical ids")
    ap.add_argument("--roster", type=Path, default=None, help="one genome per line; the slice to measure")
    ap.add_argument("--allow-missing-strand", action="store_true")
    ap.add_argument("--json", type=Path, default=None, help="write the census here as well as printing it")
    args = ap.parse_args()

    artifacts = CatalogueArtifacts(
        data_root=args.data_root,
        set_key=args.set_key,
        model_label=args.model_label,
        run_id=args.run_id,
        species_key=args.species_key,
        gene_universe_path=args.gene_universe,
        assignment_path=args.assignment,
    )
    census = run(
        artifacts,
        categorical=args.categorical,
        roster=args.roster,
        allow_missing_strand=args.allow_missing_strand,
    )
    print(report(census), flush=True)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(census.to_dict(), indent=2) + "\n")
        print(f"  → {args.json}")


if __name__ == "__main__":
    main()
