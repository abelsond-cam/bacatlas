"""Export the node universe and every claim it states — ⛔ the ONE database touch in this method.

David, 2026-10-23: *"Keep all of this running on CSD3. It doesn't need to involve bacatlas browser at
all. That was a mistake in the plan."* So this script is the whole of bacatlas's involvement: it runs
once on the laptop, against the loaded build database, and writes two small files that go to CSD3
beside the neighbour artifacts. ⛔ Nothing downstream of it opens a connection.

⭐ **It goes through `transfer/claims.py`**, which wraps `gene_ontology_fold` / `go_rung` / `ec_levels`
/ `cog_levels` rather than reimplementing them. Re-deriving a claim in a second place is how the two
halves of one answer come to disagree, and the GO closure in particular cost a page that read
*"classes differ"* for a term and its own parent.

    export BACATLAS_DATABASE_URL="postgresql+psycopg://$USER@localhost:5432/bacatlas_ingest"
    .venv/bin/python -m bacatlas_backend.scripts.export_claim_table --out-dir /tmp/claim_table

Reader: `nuna.tl.annotation_transfer.claim_table.load_claim_table`.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from bacatlas_backend.transfer import claims as C

SEPARATOR = "|"  # ⛔ must match nuna.tl.annotation_transfer.claim_table.SEPARATOR
CATALOGUES = ("ecoli-nuna4", "kp-nuna4")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--catalogues", nargs="+", default=list(CATALOGUES))
    ap.add_argument("--database-url", default=os.environ.get("BACATLAS_DATABASE_URL"))
    args = ap.parse_args()
    if not args.database_url:
        raise SystemExit("FATAL: set BACATLAS_DATABASE_URL (the build database is `bacatlas_ingest`)")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    node_rows: list[tuple[str, str, int]] = []
    claim_rows: list[tuple[str, str, str, str]] = []
    report: dict = {"catalogues": {}, "separator": SEPARATOR, "database_url_name": args.database_url.rsplit("/", 1)[-1]}

    with Session(create_engine(args.database_url)) as session:
        for catalogue in args.catalogues:
            loaded = C.load(session, catalogue)
            # ⛔ EVERY node, annotated or not — the unannotated ones are the population this work
            # exists to annotate, so they are the subject of the file and not an omission from it.
            for node, size in loaded.gene_count.items():
                node_rows.append((catalogue, loaded.node_label[node], int(size)))
            per_claim = {}
            for kind in C.CLAIM_KINDS:
                carried = loaded.claims[kind.key]
                for node, classes in carried.items():
                    # sorted, so the file is byte-stable across runs and a diff means a real change
                    claim_rows.append(
                        (catalogue, kind.key, loaded.node_label[node], SEPARATOR.join(sorted(classes)))
                    )
                per_claim[kind.key] = {
                    "label": kind.label,
                    "nodes_carrying": len(carried),
                    "node_share": round(len(carried) / loaded.node_total, 4),
                }
                print(f"  {catalogue:<12} {kind.key:<4} {len(carried):>6,} nodes "
                      f"({100 * len(carried) / loaded.node_total:.1f}%)")  # fmt: skip
            report["catalogues"][catalogue] = {
                "node_total": loaded.node_total,
                "gene_total": loaded.gene_total,
                "claims": per_claim,
            }

    nodes = args.out_dir / "nodes.tsv"
    claims = args.out_dir / "claims.tsv"
    with nodes.open("w") as fh:
        fh.write("catalogue\tnode_label\tn_genes\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in node_rows)
    with claims.open("w") as fh:
        fh.write("catalogue\tclaim\tnode_label\tclasses\n")
        fh.writelines("\t".join(r) + "\n" for r in claim_rows)
    (args.out_dir / "claim_table.json").write_text(json.dumps(report, indent=2))
    print(f"\nwrote {nodes.name} ({len(node_rows):,} nodes), {claims.name} ({len(claim_rows):,} rows), claim_table.json")


if __name__ == "__main__":
    main()
