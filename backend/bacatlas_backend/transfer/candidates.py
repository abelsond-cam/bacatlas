"""Candidate donors for a node — ⛔ pooled across both representations, with NO donor chosen.

⭐ **This is the module that fixes the defect the old design was built on.** The old scheme took the
ESM-nearest labelled node and the Bacformer-nearest labelled node — *different nodes, 24 % of the
time* — and fed one covariate from each into one regression. David, 2026-10-12: *"the nearest in these
two may not be the same node … That makes little sense in doing the regression as we are not doing the
slope on same samples."*

So nothing here chooses. Every annotated neighbour either representation offers becomes a candidate,
and each carries its similarity under **both** representations measured on **its own pair**. The
disagreement between the two representations is preserved as data rather than resolved by a rule,
because it is evidence:

> *"if N1 and N2 and N3 have relative prediction probability of 0.95, 0.9 and 0.85 and all agree,
> this is surely a much stronger evidence than N1 and N2 having different annotations."* — David

⚠ **A candidate's similarity under the OTHER representation is often unknown today**, because the
database stores five neighbours per representation and the lists differ. That is the gap S2 closes by
rebuilding at depth 100; until then `missing_pair_similarity` counts it, and ⛔ it is counted rather
than filled, because substituting a different pair's number is precisely the defect above.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from sqlalchemy import select
from sqlalchemy.orm import Session

from bacatlas_backend.models.enumerations import EmbeddingRepresentation
from bacatlas_backend.models.locus import Locus
from bacatlas_backend.models.locus_similarity import LocusNearestLocus
from bacatlas_backend.transfer.claims import Catalogue

REPRESENTATIONS = ("ESM", "BACFORMER")


@dataclass(frozen=True)
class Candidate:
    """One possible donor for one recipient — a node, and the geometry of that one pair.

    ⚠ `similarity[rep]` is `None` where that representation's stored shortlist does not contain this
    pair. `None` means **not stored**, which is a different statement from *not similar*; treating it
    as a low similarity would manufacture a correlation with shortlist depth.
    """

    donor: int
    similarity: dict[str, float | None]
    found_by: frozenset[str]

    def ranked_under(self, representation: str) -> bool:
        """Did this representation's own shortlist offer this candidate?"""
        return representation in self.found_by


def within_species(session: Session, catalogue: Catalogue) -> dict[int, list[tuple[str, int, int, float]]]:
    """Every stored neighbour of every node, as `(representation, rank, neighbour, similarity)`.

    ⚠ The stored similarity is the **median** cosine over gene pairs (`locus_nearest_locus`'s own
    column comment). The max is computed upstream and dropped at export — the asymmetry S2 fixes.
    """
    walks: dict[int, list[tuple[str, int, int, float]]] = defaultdict(list)
    for representation in (EmbeddingRepresentation.ESM, EmbeddingRepresentation.BACFORMER):
        for node, rank, neighbour, similarity in session.execute(
            select(
                LocusNearestLocus.locus_id,
                LocusNearestLocus.rank,
                LocusNearestLocus.neighbour_locus_id,
                LocusNearestLocus.cross_similarity,
            )
            .join(Locus, Locus.locus_id == LocusNearestLocus.locus_id)
            .where(
                Locus.pangenome_id == catalogue.pangenome_id,
                LocusNearestLocus.representation == representation,
                LocusNearestLocus.cross_similarity.is_not(None),
            )
            .order_by(LocusNearestLocus.locus_id, LocusNearestLocus.rank)
        ).all():
            walks[node].append((representation.name, rank, neighbour, similarity))
    return dict(walks)


def read_cross_map(path: str, *, donor: str, recipient: str) -> tuple[dict[str, list[dict]], dict]:
    """A cross-species map keyed by recipient node LABEL, plus its provenance.

    ⛔⛔ **THE DIRECTION IS NOT IN THE HEADER.** A map in either direction carries both `*_node`
    columns, so a header check reads an `ecoli→kp` map as a `kp→ecoli` one without complaint and every
    rank then means the opposite of what it says. The direction comes from the provenance JSON, which
    records it explicitly; the header is a second, independent check.

    ⚠ The sidecar is REQUIRED. It carries the shortlist's recall audit, and a map whose recall was
    never measured is not something to report transfer rates from.
    """
    sidecar = Path(path).with_suffix(".json")
    if not sidecar.exists():
        raise SystemExit(f"FATAL: no provenance beside the map — expected {sidecar}")
    provenance = json.loads(sidecar.read_text())
    stated = (provenance.get("donor"), provenance.get("recipient"))
    if stated != (donor, recipient):
        raise SystemExit(
            f"FATAL: {sidecar.name} says donor={stated[0]} recipient={stated[1]}; "
            f"this run was told {donor} → {recipient}"
        )
    rows: dict[str, list[dict]] = defaultdict(list)
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = (f"{recipient}_node", f"{donor}_node", "cross", "cross_max", "rank_cross", "rank_cross_max")
        missing = [column for column in required if column not in (reader.fieldnames or ())]
        if missing:
            raise SystemExit(f"FATAL: {path} lacks {missing}; columns are {reader.fieldnames}")
        for row in reader:
            rows[row[f"{recipient}_node"]].append(
                {
                    "donor_label": row[f"{donor}_node"],
                    "cross": float(row["cross"]),
                    "cross_max": float(row["cross_max"]),
                    "rank_cross": int(row["rank_cross"]),
                    "rank_cross_max": int(row["rank_cross_max"]),
                    "n_pairs": int(row["n_pairs"]),
                    "n_donor_genes": int(row[f"n_{donor}_genes"]),
                }
            )
    return dict(rows), provenance
