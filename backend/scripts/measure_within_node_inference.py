"""How many GENES gain a function because their node carries one — across every axis, not just COG.

⛔ **This is the question `measure_cog_function_inference.py` answers for two vocabularies only.**
That script was written to assess the NEIGHBOUR ladder, and its within-node section inherited the
same two axes. But a node is the unit of annotation: if a locus carries a GO term on 3 of its 100
genes, the clustering has given a function to all 100. Counting only COG and EC understates what the
model delivers, and GO is the axis it understates most — 351,559 of 489,146 *E. coli* genes carry a
GO term against 313,617 for COG.

Definitions, stated because every number below depends on them:
  * a node **carries** an axis when at least one member gene has it;
  * a gene **gains** the axis when it is in such a node and has none of its own;
  * the node is **checkable** when at least two member genes carry it — one is unanimous by
    construction, which is the absence of evidence (518 of 4,993 *E. coli* COG nodes).
⚠ The gain is what the CLUSTERING delivers. It is not an independent prediction: the evidence is
that these genes were grouped together, which is the model's claim, not a second opinion on it.
"""

from __future__ import annotations

import os
import sys

from sqlalchemy import create_engine, text
from sqlalchemy.orm import Session

CATALOGUES = ("ecoli-nuna4", "kp-nuna4")

#: (label, the gene-level predicate that means "this gene carries it")
AXES = (
    ("COG orthogroup", "f.cog_id is not null"),
    ("GO term", "f.gene_ontology_terms is not null"),
    ("EC number", "f.ec_numbers is not null"),
    ("KEGG orthology", "f.kegg_orthology_id is not null"),
    ("any of the four", "(f.cog_id is not null or f.gene_ontology_terms is not null"
                        " or f.ec_numbers is not null or f.kegg_orthology_id is not null)"),
    # ⚠ The NAME, which is a different kind of claim and reported apart from the four. Bakta emits a
    # genome-private locus tag where it has no symbol, so the predicate asks for a real symbol.
    ("gene symbol", "g.bakta_gene_symbol is not null"),
)

QUERY = """
with member as (
    select m.locus_id,
           ({predicate}) as carries
      from gene_locus_membership m
      join locus l on l.locus_id = m.locus_id
      join pangenome p on p.pangenome_id = l.pangenome_id
      join gene g on g.genome_id = m.genome_id and g.flat_index = m.flat_index
      left join gene_functional_annotation f
             on f.genome_id = m.genome_id and f.flat_index = m.flat_index
     where p.catalogue_key = :catalogue
),
node as (
    select locus_id,
           count(*) as genes,
           count(*) filter (where carries) as carrying
      from member group by locus_id
)
select
    (select sum(genes) from node)                                              as total_genes,
    (select count(*) from node)                                                as total_nodes,
    coalesce(sum(carrying), 0)                                                 as genes_already,
    count(*) filter (where carrying > 0)                                       as nodes_carrying,
    coalesce(sum(genes - carrying) filter (where carrying >= 2), 0)            as gained_checkable,
    coalesce(sum(genes - carrying) filter (where carrying = 1), 0)             as gained_one_gene,
    count(*) filter (where carrying = 1)                                       as nodes_one_gene
  from node
"""


def main() -> None:
    """Print, per catalogue and axis, what the clustering hands to genes that had nothing."""
    url = os.environ.get("BACATLAS_DATABASE_URL")
    if not url:
        sys.exit("BACATLAS_DATABASE_URL is not set — see backend/README.md")
    with Session(create_engine(url)) as session:
        for catalogue in CATALOGUES:
            print(f"\n{'=' * 108}\n{catalogue}")
            print(f"  {'axis':<18}{'nodes':>8}{'genes have':>12}{'→ %':>7}"
                  f"{'GAINED':>10}{'→ % after':>11}{'(unverifiable)':>16}{'nodes w/ 1':>12}")
            print("  " + "-" * 94)
            for label, predicate in AXES:
                row = session.execute(
                    text(QUERY.format(predicate=predicate)), {"catalogue": catalogue}
                ).one()
                before = 100 * row.genes_already / row.total_genes
                after = 100 * (row.genes_already + row.gained_checkable) / row.total_genes
                print(
                    f"  {label:<18}{row.nodes_carrying:>8,}{row.genes_already:>12,}{before:>6.1f}%"
                    f"{row.gained_checkable:>10,}{after:>10.1f}%"
                    f"{row.gained_one_gene:>16,}{row.nodes_one_gene:>12,}"
                )
            print(f"  total: {row.total_genes:,} genes in {row.total_nodes:,} nodes")

            # ⛔ The three framings answer three different questions and must never be conflated.
            # Per-axis gains CANNOT be added: a gene lacking a COG usually already has a GO, so the
            # sum counts the same gene repeatedly. The union is the honest "is this gene annotated
            # at all" figure, and it is SMALLER than the COG column alone.
            totals = session.execute(
                text("""
                    with member as (
                        select m.locus_id,
                               (f.cog_id is not null) as cog,
                               (f.gene_ontology_terms is not null) as go,
                               (f.ec_numbers is not null) as ec,
                               (f.kegg_orthology_id is not null) as kegg
                          from gene_locus_membership m
                          join locus l on l.locus_id = m.locus_id
                          join pangenome p on p.pangenome_id = l.pangenome_id
                          left join gene_functional_annotation f
                                 on f.genome_id = m.genome_id and f.flat_index = m.flat_index
                         where p.catalogue_key = :catalogue
                    ),
                    node as (
                        select locus_id, count(*) genes,
                               count(*) filter (where cog) c, count(*) filter (where go) g,
                               count(*) filter (where ec) e, count(*) filter (where kegg) k,
                               count(*) filter (where cog or go or ec or kegg) any_of
                          from member group by locus_id
                    )
                    select
                      sum(case when c >= 2 then genes - c else 0 end)
                      + sum(case when g >= 2 then genes - g else 0 end)
                      + sum(case when e >= 2 then genes - e else 0 end)
                      + sum(case when k >= 2 then genes - k else 0 end)   as assignments_created,
                      sum(genes) filter (where any_of = 0)                as genes_in_dark_nodes,
                      count(*) filter (where any_of = 0)                  as dark_nodes
                      from node
                """),
                {"catalogue": catalogue},
            ).one()
            print(f"  (gene, axis) ASSIGNMENTS created across the four: "
                  f"{totals.assignments_created:,}  — ⚠ not genes; the axes overlap")
            print(f"  genes in nodes where NO member carries any of the four: "
                  f"{totals.genes_in_dark_nodes:,} "
                  f"({100 * totals.genes_in_dark_nodes / row.total_genes:.1f} %) "
                  f"in {totals.dark_nodes:,} nodes")


if __name__ == "__main__":
    main()
