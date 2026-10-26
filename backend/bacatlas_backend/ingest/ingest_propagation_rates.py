"""Write one catalogue's unanimity tally — twenty rows that replace a 1.3 M-row scan per request.

⛔ **Counted from the SAME per-locus verdicts that are written to `locus`, deliberately, and the
docstring says so rather than claiming an independence it does not have.** The tally and the
verdict columns are two writes of one in-memory list, so a mistake in the DECISION moves both. What
the pair does catch is everything between the decision and the database: a column projected into the
wrong position, an enum coerced to its value instead of its name, a partial `COPY`, a band string
that arrived as something `PrevalenceBand` does not accept. `GROUP BY` the stored verdicts and
compare — `tests/test_locus_annotation_rollups.py` does — and the write path is proved.

⭐ The decision itself is gated differently and elsewhere: against
`within_node_propagation.compute_propagation` reading the gene rows, which is a different reader of
a different store. Neither gate covers the other, which is why both exist.
"""

from __future__ import annotations

from collections import defaultdict

from sqlalchemy.orm import Session

from bacatlas_backend.ingest.derive_locus_propagation import PROPAGATION_VOCABULARIES
from bacatlas_backend.instruments.within_node_propagation import VERDICT_TALLY
from bacatlas_backend.models.enumerations import PrevalenceBand, PropagationVerdict
from bacatlas_backend.models.pangenome_propagation_rate import PangenomePropagationRate


def write_propagation_rates(
    session: Session,
    *,
    pangenome_id: int,
    bands: list[str],
    verdicts: dict[str, list[PropagationVerdict | None]],
) -> int:
    """Tally the per-locus verdicts into (vocabulary × band) rows and insert them. Returns the count.

    `bands[i]` and `verdicts[suffix][i]` describe the same locus. A vocabulary whose verdicts are all
    `None` was never measured and gets **no rows at all** — the same statement its NULL column makes,
    and the same thing `compute_propagation` does today when a catalogue carries no dbxrefs.

    ⚠ **A band present in the catalogue gets a row even when every number in it is zero.** RARE is
    that band in every catalogue measured so far: zero checkable nodes beside hundreds of one-gene
    ones. `services/function_inference_service.propagation_for` also substitutes zeros for a band it
    cannot find, which is belt and braces on purpose — a *missing key* reaching the page where it
    expects *not measurable* is the three-states-into-two collapse `no_homology` was retired for.
    """
    rows = []
    for suffix, annotation_kind in PROPAGATION_VOCABULARIES.items():
        column = verdicts[suffix]
        if all(verdict is None for verdict in column):
            continue
        tally: dict[str, list[int]] = defaultdict(lambda: [0, 0, 0])
        for band, verdict in zip(bands, column, strict=True):
            counts = tally[band]
            # ⛔ `VERDICT_TALLY` and not a hand-written `>= 2`: `checkable` means *agreed OR
            # disagreed*, so a `not_comparable` node belongs in neither column. Deriving it as
            # "two or more annotated genes" quietly moves 81 ecoli GO loci into the denominator and
            # drops the published rate with no test turning red.
            for index, increment in enumerate(VERDICT_TALLY[verdict]):
                counts[index] += increment
        for band, (checkable, unanimous, one_gene) in tally.items():
            rows.append(
                PangenomePropagationRate(
                    pangenome_id=pangenome_id,
                    annotation_kind=annotation_kind,
                    prevalence_band=PrevalenceBand(band) if not isinstance(band, PrevalenceBand) else band,
                    checkable=checkable,
                    unanimous=unanimous,
                    one_gene=one_gene,
                )
            )
    session.add_all(rows)
    session.flush()
    return len(rows)
