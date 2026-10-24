"""Where a genome's bases come from — the one seam the Sequence tab reads through.

⭐ **The seam already existed; this only widens it.** `models/gene.py` states the design: *"Everything
the Sequence tab shows is a column here EXCEPT the base letters themselves."* Coordinates, strand,
phase, partiality, GC, protein length, symbol and product are all `gene` columns, deliberately split
from `gene_locus_membership` so a model re-run never rewrites 412 M coordinate rows. Only the base
letters came off the disk, through exactly one expression. So the change is to replace *one
resolver* — "build a path from `bakrep_dataset_id`, open it" — with "get this genome's contigs, from
wherever", keeping the `dict[seqid, str]` contract. `read_gene_sequence` is untouched, and so are
the schema, the API shape and the translation rule.

⛔ **The order is HELD-FIRST, and that is what makes this safe to turn on.** A genome whose Bakta
GFF is on this server is read from it exactly as before, so enabling ATB streaming is purely
*additive*: it can only change the answer for genomes that previously returned
``sequence_unavailable``. There is no cohort for which this swaps a local answer for a remote one.

```
1. the local Bakta GFF        held genomes — unchanged behaviour, ~30 ms
2. the ATB disk cache         previously fetched — ~30 ms
3. AllTheBacteria             one GET, measured mean 0.241 s / p95 0.330 s
4. SequenceUnavailable        named, and distinguished from "no gene here"
```

⛔ **Three answers that must never look alike**, which is the published page's own rule
(`app.js:4604`) and the reason this module raises rather than returning an empty map:

| the situation | what the caller does |
|---|---|
| the genome has no gene at this locus | returns `[]` — an **answer**, rendered as one |
| the gene is there, the bases are not readable | raises, and the endpoint says why |
| the locus or genome does not exist | 404, upstream of here |

⚠ **The parsed cache is shared across sources deliberately.** Whichever route supplied a genome, it
is cached the same way and the next locus on the same genome is free. That is also why the parse
lives here rather than in the store: there is one `lru_cache`, not one per source.
"""

from __future__ import annotations

import logging
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path

from bacatlas_backend.gff.gff_cds_parser import parse_genome_annotation
from bacatlas_backend.sequence_source.atb_assembly_store import AssemblyUnavailable, AtbAssemblyStore

logger = logging.getLogger(__name__)

#: How many genomes' contig sequences to keep parsed, per process. ⚠ A bacterial genome is ~5.3 MB
#: as Python strings (measured on SAMN04393397, 189 contigs, 5,328,536 bases), so this is a bound in
#: the tens of MB per worker — deliberately small, because a reader walking a track hits ONE genome.
PARSED_GENOME_CACHE_SIZE = 4


class SequenceSourceUnavailable(LookupError):
    """No source could supply this genome's bases. ``is_retryable`` carries whether to try again.

    ⚠ Mapped onto `SequenceUnavailable` by the service, which is what the endpoint catches — this
    type exists so the resolver can be tested without importing the service layer.
    """

    def __init__(self, message: str, *, is_retryable: bool = False) -> None:
        super().__init__(message)
        self.is_retryable = is_retryable


@dataclass(frozen=True)
class GenomeSequenceRequest:
    """Which genome, and everything needed to locate its held copy if there is one."""

    sample_id: str
    #: The BakRep dataset directory. ⚠ `None` is ordinary now, not a failure: a streamed genome has
    #: no local directory at all, where before a missing value meant the genome was unreadable.
    bakrep_dataset_id: str | None


class ContigSequenceProvider:
    """Resolve one genome's `{seqid: bases}` through the held → cached → remote chain."""

    def __init__(
        self,
        *,
        gff_root: Path | None,
        atb_store: AtbAssemblyStore | None,
        cache_size: int = PARSED_GENOME_CACHE_SIZE,
    ) -> None:
        self._gff_root = gff_root
        self._atb_store = atb_store
        # ⚠ An `OrderedDict` rather than `functools.lru_cache` on a method: the decorator keys on
        # `self` too, so every provider instance would keep its own entries alive for the lifetime
        # of the process — a leak that shows only under the test suite, which builds many. And an
        # OrderedDict rather than a list of keys because `move_to_end` is the LRU; "append on every
        # touch and dedupe" is the version of this that silently evicts the hot entry.
        self._parsed: OrderedDict[str, dict[str, str]] = OrderedDict()
        self._cache_size = max(1, cache_size)

    @property
    def can_stream(self) -> bool:
        """Whether a genome we do not hold can be fetched at all."""
        return self._atb_store is not None

    def contig_sequences(self, request: GenomeSequenceRequest) -> dict[str, str]:
        """The genome's contigs, by seqid. Raises `SequenceSourceUnavailable` with a reason."""
        hit = self._parsed.get(request.sample_id)
        if hit is not None:
            self._parsed.move_to_end(request.sample_id)
            return hit

        local_path = self._held_path(request)
        if local_path is not None and local_path.exists():
            sequences = self._parse_held(local_path)
        elif self._atb_store is not None:
            sequences = self._fetch_streamed(request.sample_id)
        else:
            raise SequenceSourceUnavailable(
                self._nothing_to_read_message(request, local_path), is_retryable=False
            )

        self._remember(request.sample_id, sequences)
        return sequences

    # ------------------------------------------------------------------- held

    def _held_path(self, request: GenomeSequenceRequest) -> Path | None:
        if self._gff_root is None or not request.bakrep_dataset_id:
            return None
        return self._gff_root / request.bakrep_dataset_id / request.sample_id / f"{request.sample_id}.bakta.gff3.gz"

    def _parse_held(self, path: Path) -> dict[str, str]:
        parsed = parse_genome_annotation(path, want_sequence=True)
        if not parsed.carries_sequence:
            # ⚠ Checked rather than assumed. All 280 probe GFFs carry a `##FASTA` block (measured
            # 2026-09-04), and one that does not must fail with a sentence rather than a KeyError.
            raise SequenceSourceUnavailable(
                f"{path.name} carries no ##FASTA block, so it holds no bases", is_retryable=False
            )
        return parsed.contig_sequences

    # --------------------------------------------------------------- streamed

    def _fetch_streamed(self, sample_id: str) -> dict[str, str]:
        assert self._atb_store is not None  # guarded by the caller
        try:
            return self._atb_store.contig_sequences(sample_id)
        except AssemblyUnavailable as error:
            raise SequenceSourceUnavailable(str(error), is_retryable=error.is_retryable) from error

    def _nothing_to_read_message(self, request: GenomeSequenceRequest, local_path: Path | None) -> str:
        """⚠ The message names WHICH of the three ways this can happen, because the fix differs."""
        if self._gff_root is None:
            return (
                "this server has no annotation store configured and no assembly source enabled, so "
                "no sequence can be read (set BACATLAS_ROOT_GFF or BACATLAS_ATB_STREAMING)"
            )
        if local_path is None:
            return (
                f"{request.sample_id} has no recorded annotation directory and assembly streaming is "
                "disabled, so its bases cannot be located"
            )
        return (
            f"{request.sample_id}'s annotation file is not on this server and assembly streaming is "
            "disabled, so its bases cannot be read"
        )

    # ------------------------------------------------------------ parsed cache

    def _remember(self, sample_id: str, sequences: dict[str, str]) -> None:
        self._parsed[sample_id] = sequences
        self._parsed.move_to_end(sample_id)
        while len(self._parsed) > self._cache_size:
            self._parsed.popitem(last=False)

    def clear(self) -> None:
        """Drop the parsed genomes. For tests, and for a process that wants its memory back."""
        self._parsed.clear()
