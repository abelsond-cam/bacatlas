"""Where a genome's bases come from — the one seam the Sequence tab reads through.

⭐ **The seam already existed; this only widens it.** Everything the Sequence tab shows is in the
database EXCEPT the base letters themselves — coordinates, strand, phase and partiality in
`genome_genes`, kept apart from the model's `genome_locus_assignment` so a model re-run never
rewrites them. Only the base letters came off the disk, through exactly one expression. So the change is to replace *one
resolver* — "build a path from `bakrep_dataset_id`, open it" — with "get this genome's contigs, from
wherever", keeping the `dict[seqid, str]` contract. `read_gene_sequence` is untouched, and so are
the schema, the API shape and the translation rule.

⛔ **The order is HELD-FIRST, and that is what makes this safe to turn on.** A genome whose Bakta
GFF is on this server is read from it exactly as before, so enabling ATB streaming is purely
*additive*: it can only change the answer for genomes that previously returned
``sequence_unavailable``. There is no cohort for which this swaps a local answer for a remote one.

```
1. the local Bakta GFF        held genomes — unchanged behaviour, ~30 ms
2. the disk cache             previously fetched — ~30 ms (ATB's and NCBI's, side by side)
3a. AllTheBacteria            BioSample-keyed genomes: one GET, measured mean 0.241 s / p95 0.330 s
3b. NCBI Datasets             RefSeq/GenBank-keyed genomes: one GET, measured 1.08 s, rate-limited
4. SequenceUnavailable        named, and distinguished from "no gene here"
```

⛔ **3a and 3b are ONE station chosen by key kind, never a chain of two.** `genome.sample_id_kind` says
what the genome's key IS: a BioSample id (`SAMN…`) names an AllTheBacteria object and is nothing to NCBI;
an assembly stem (`GCF_…_genomic`) names an NCBI assembly, whose `bare_accession` is the thing to ask for,
and the ATB bucket holds it under no name (measured 2026-10-29: the cohort's 5,389 `GCF_`/`GCA_` stems were
found **0 times** in the bucket). Asking the wrong source costs a request, a rate token, and a 404 that the
negative cache would then remember for ten minutes — so the kind decides, and only that source is asked.

⚠ **An NCBI assembly answers in its own namespace.** Its deflines carry `NZ_CP003999.1` where our
coordinates (built from a GenBank-style GFF) say `CP003999.1`. The request therefore carries the contig
names the caller wants, and for an NCBI-keyed genome this provider returns the map RE-KEYED to them via
`contig_namespace.resolve_contig_names` — exact first, else one `NZ_` prefix rule proved on the whole
genome, else a named refusal. The parsed cache holds the assembly as NCBI keyed it; the re-keying is per
call and costs a set comparison over a few hundred names. BioSample genomes are untouched: the held GFF
and the ATB object write the same first token (280/280 measured), so the service's own lookup stays the
one that speaks for them.

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
from bacatlas_backend.models.enumerations import SampleIdentifierKind
from bacatlas_backend.sequence_source.atb_assembly_store import AssemblyUnavailable, AtbAssemblyStore
from bacatlas_backend.sequence_source.contig_namespace import ContigNamespaceMismatch, resolve_contig_names
from bacatlas_backend.sequence_source.ncbi_assembly_store import NcbiAssemblyStore

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
    #: ⭐ Which remote station may be asked when the genome is not held — see the module docstring.
    #: `None` means the caller did not say, and the chain behaves as it did before NCBI existed: ATB.
    sample_id_kind: SampleIdentifierKind | None = None
    #: `GCF_000512165.1` for an assembly stem — what NCBI is asked for. `None` for a BioSample genome.
    bare_accession: str | None = None
    #: The contig names the caller's coordinates are stored under. For an NCBI-keyed genome the answer
    #: is re-keyed to these (the `NZ_` rule); `None` returns the assembly exactly as its source keyed it.
    wanted_seqids: frozenset[str] | None = None


class ContigSequenceProvider:
    """Resolve one genome's `{seqid: bases}` through the held → cached → remote chain."""

    def __init__(
        self,
        *,
        gff_root: Path | None,
        atb_store: AtbAssemblyStore | None,
        ncbi_store: NcbiAssemblyStore | None = None,
        cache_size: int = PARSED_GENOME_CACHE_SIZE,
    ) -> None:
        self._gff_root = gff_root
        self._atb_store = atb_store
        self._ncbi_store = ncbi_store
        # ⚠ An `OrderedDict` rather than `functools.lru_cache` on a method: the decorator keys on
        # `self` too, so every provider instance would keep its own entries alive for the lifetime
        # of the process — a leak that shows only under the test suite, which builds many. And an
        # OrderedDict rather than a list of keys because `move_to_end` is the LRU; "append on every
        # touch and dedupe" is the version of this that silently evicts the hot entry.
        self._parsed: OrderedDict[str, dict[str, str]] = OrderedDict()
        self._cache_size = max(1, cache_size)

    @property
    def can_stream(self) -> bool:
        """Whether a genome we do not hold can be fetched at all, from either remote station."""
        return self._atb_store is not None or self._ncbi_store is not None

    @property
    def streams_atb(self) -> bool:
        """Whether BioSample-keyed genomes we do not hold are fetched from AllTheBacteria."""
        return self._atb_store is not None

    @property
    def streams_ncbi(self) -> bool:
        """Whether RefSeq/GenBank-keyed genomes are fetched from NCBI Datasets."""
        return self._ncbi_store is not None

    @property
    def stations(self) -> tuple[str, ...]:
        """The configured sources, in chain order — what `/health` prints, so nobody has to infer it."""
        names: list[str] = []
        if self._gff_root is not None:
            names.append("held_annotation_tree")
        if self._atb_store is not None:
            names.append("allthebacteria")
        if self._ncbi_store is not None:
            names.append("ncbi_datasets")
        return tuple(names)

    def contig_sequences(self, request: GenomeSequenceRequest) -> dict[str, str]:
        """The genome's contigs, by seqid. Raises `SequenceSourceUnavailable` with a reason."""
        sequences = self._parsed.get(request.sample_id)
        if sequences is not None:
            self._parsed.move_to_end(request.sample_id)
        else:
            sequences = self._read_through_the_chain(request)
            self._remember(request.sample_id, sequences)
        return self._in_the_callers_namespace(request, sequences)

    def _read_through_the_chain(self, request: GenomeSequenceRequest) -> dict[str, str]:
        local_path = self._held_path(request)
        if local_path is not None and local_path.exists():
            return self._parse_held(local_path)
        if request.sample_id_kind is SampleIdentifierKind.ASSEMBLY_STEM:
            # ⛔ NCBI or nothing — never ATB, which holds no stem under any name (see the docstring).
            if self._ncbi_store is None:
                raise SequenceSourceUnavailable(
                    f"{request.sample_id} is a RefSeq/GenBank assembly that is not on this server, and "
                    "NCBI streaming is disabled, so its bases cannot be read (set BACATLAS_NCBI_STREAMING)",
                    is_retryable=False,
                )
            return self._fetch_from_ncbi(request)
        if self._atb_store is not None:
            return self._fetch_streamed(request.sample_id)
        raise SequenceSourceUnavailable(self._nothing_to_read_message(request, local_path), is_retryable=False)

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

    def _fetch_from_ncbi(self, request: GenomeSequenceRequest) -> dict[str, str]:
        assert self._ncbi_store is not None  # guarded by the caller
        if not request.bare_accession:
            # ⚠ A data gap named as one. `build_genome_set.bare_accession` fills this for every stem, so
            # an empty value is an ingest fault — not something to repair by parsing the stem here.
            raise SequenceSourceUnavailable(
                f"{request.sample_id} has no bare accession recorded, so NCBI cannot be asked for it",
                is_retryable=False,
            )
        try:
            return self._ncbi_store.contig_sequences(request.bare_accession)
        except AssemblyUnavailable as error:
            raise SequenceSourceUnavailable(str(error), is_retryable=error.is_retryable) from error

    # ------------------------------------------------------------- namespace

    @staticmethod
    def _in_the_callers_namespace(request: GenomeSequenceRequest, sequences: dict[str, str]) -> dict[str, str]:
        """For an NCBI-keyed genome, `{wanted seqid: bases}`; for everything else, the map as keyed."""
        if request.wanted_seqids is None or request.sample_id_kind is not SampleIdentifierKind.ASSEMBLY_STEM:
            return sequences
        try:
            resolved = resolve_contig_names(request.wanted_seqids, sequences)
        except ContigNamespaceMismatch as error:
            # ⛔ Permanent, and about THIS genome's coordinates against THIS assembly — never a reason to
            # ask again, and never to be read as a missing file.
            raise SequenceSourceUnavailable(f"{request.sample_id}: {error}", is_retryable=False) from error
        if any(wanted != header for wanted, header in resolved.items()):
            logger.info("%s: contig names resolved by the NZ_ prefix rule", request.sample_id)
        return {wanted: sequences[header] for wanted, header in resolved.items()}

    def _nothing_to_read_message(self, request: GenomeSequenceRequest, local_path: Path | None) -> str:
        """⚠ The message names WHICH of the three ways this can happen, because the fix differs."""
        if self._gff_root is None:
            return (
                "this server has no annotation store configured and no assembly source enabled, so "
                "no sequence can be read (set BACATLAS_ROOT_GFF, BACATLAS_ATB_STREAMING or "
                "BACATLAS_NCBI_STREAMING)"
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
