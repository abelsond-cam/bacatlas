"""Regenerate the front end's real-bytes fixtures from the live API — BOTH species.

Writes, per species, `frontend/tests/fixtures/api_responses_{species}.json`. (It also wrote the two
catalogue sprites until the page stopped drawing them, 2026-09-22.) Run from `backend/` with
`BACATLAS_DATABASE_URL` set, both catalogues loaded and published, and `BACATLAS_ROOT_GFF` set (the
sequence endpoint reads the GFFs):

    python scripts/dump_api_locus_fixture.py

⭐ **Every case is chosen for what it EXERCISES, by a query that says so** — never "the first locus".
The front-end suite asserts each case still exercises its property, so a regenerate that picked an
ordinary locus for a discriminating case fails loudly rather than quietly testing nothing:

- `ordinary` — every arrangement listed, both member remainders zero, and a MEASURED ZERO ESM
  within-medoid distance (the case a truthiness test drops).
- `over_cap` — more arrangements than the API lists, so members sit past the cap.
- `no_window` — the most members with no recorded neighbourhood at all.
- `function_rich` — EC and KEGG and all three GO namespaces, the parts a typical locus lacks.
- `inferred_function` — **no COG of its own**, and a COG-bearing ESM neighbour close enough to quote,
  so the inferred-function card renders on real bytes rather than only on a hand-built object. ⚠ Its
  first annotated donor must sit past rank 1, because the walk reporting the empty ranks before the
  donor is the part a fixture with a rank-1 donor would never exercise.
- sequences: a MINUS-strand gene well inside its contig (where the flank orientation can be wrong and
  still look plausible), and a genome with NO gene at that locus (an answer, not a failure).
"""

from __future__ import annotations

import json
import pathlib

from sqlalchemy import select, text
from sqlalchemy.orm import Session

from bacatlas_backend.application_factory import create_application
from bacatlas_backend.configuration import Configuration
from bacatlas_backend.models.gene import Gene, GeneLocusMembership
from bacatlas_backend.models.genome import Genome
from bacatlas_backend.models.genome_collection import GenomeCollectionMembership
from bacatlas_backend.models.locus import Locus
from bacatlas_backend.models.pangenome import Pangenome
from bacatlas_backend.models.pathogen_species import PathogenSpecies

SPECIES_KEYS = ("ecoli", "kp")
FIXTURES = pathlib.Path(__file__).resolve().parents[2] / "frontend" / "tests" / "fixtures"
#: The API's arrangement display cap (`locus_detail_service`). `ordinary` must fit inside it.
ARRANGEMENT_CAP = 8


def _one(session, statement, what):
    value = session.execute(statement).scalar_one_or_none()
    if value is None:
        raise SystemExit(f"no locus in this catalogue exercises `{what}` — the fixture cannot be built")
    return value


def choose_cases(session: Session, pangenome_id: int) -> dict[str, str]:
    """One locus label per case, each selected by the property it exists to exercise."""
    in_catalogue = Locus.pangenome_id == pangenome_id
    return {
        "ordinary": _one(
            session,
            select(Locus.node_label)
            .where(
                in_catalogue,
                Locus.total_arrangement_count.between(2, 5),
                Locus.member_gene_count == Locus.arrangement_member_gene_count,
            )
            .order_by(Locus.member_gene_count.desc(), Locus.catalogue_ordinal)
            .limit(1),
            "ordinary",
        ),
        "over_cap": _one(
            session,
            select(Locus.node_label)
            .where(in_catalogue, Locus.total_arrangement_count > ARRANGEMENT_CAP + 4)
            .order_by(Locus.total_arrangement_count.desc(), Locus.catalogue_ordinal)
            .limit(1),
            "over_cap",
        ),
        "no_window": _one(
            session,
            select(Locus.node_label)
            .where(in_catalogue, Locus.member_gene_count > Locus.arrangement_member_gene_count)
            .order_by(
                (Locus.member_gene_count - Locus.arrangement_member_gene_count).desc(), Locus.catalogue_ordinal
            )
            .limit(1),
            "no_window",
        ),
        "function_rich": _one(
            session,
            select(Locus.node_label)
            .where(
                in_catalogue,
                Locus.ec_annotated_member_count > 0,
                Locus.kegg_annotated_member_count > 0,
                Locus.go_annotated_member_count_molecular_function > 0,
                Locus.go_annotated_member_count_biological_process > 0,
                Locus.go_annotated_member_count_cellular_component > 0,
            )
            .order_by(Locus.member_gene_count.desc(), Locus.catalogue_ordinal)
            .limit(1),
            "function_rich",
        ),
        # ⭐ A node with NO COG whose first COG-bearing ESM neighbour sits past rank 1 — so the
        # fixture exercises both the suggestion and the reported walk over the ranks that carried
        # nothing. For *E. coli* COG only 1,679 of 4,143 transfers come from rank 1, so a rank-1
        # donor would be the unrepresentative case.
        # ⚠ Plain SQL: the "no NEARER neighbour carries one" clause correlates a second scan of
        # `locus_nearest_locus` to the first, which reads clearly here and not in the ORM.
        "inferred_function": _one(
            session,
            text("""
                select focal.node_label
                  from locus focal
                  join locus_nearest_locus n
                       on n.locus_id = focal.locus_id and n.representation = 'ESM'
                  join locus donor on donor.locus_id = n.neighbour_locus_id
                 where focal.pangenome_id = :pangenome_id
                   and focal.cog_distinct_id_count = 0
                   and donor.cog_distinct_id_count > 0
                   and n.cross_similarity >= 0.98
                   and n.rank >= 2
                   and not exists (
                         select 1
                           from locus_nearest_locus earlier
                           join locus e on e.locus_id = earlier.neighbour_locus_id
                          where earlier.locus_id = focal.locus_id
                            and earlier.representation = 'ESM'
                            and earlier.rank < n.rank
                            and e.cog_distinct_id_count > 0)
                 order by focal.member_gene_count desc, focal.catalogue_ordinal
                 limit 1
            """).bindparams(pangenome_id=pangenome_id),
            "inferred_function",
        ),
    }


def choose_sequences(session: Session, pangenome_id: int) -> dict[str, tuple[str, str]]:
    """A minus-strand gene well inside its contig, and a genome with no gene at that locus."""
    minus = session.execute(
        select(Genome.sample_id, Locus.node_label)
        .select_from(Gene)
        .join(
            GeneLocusMembership,
            (GeneLocusMembership.genome_id == Gene.genome_id) & (GeneLocusMembership.flat_index == Gene.flat_index),
        )
        .join(Genome, Genome.genome_id == Gene.genome_id)
        .join(Locus, Locus.locus_id == GeneLocusMembership.locus_id)
        .where(
            GeneLocusMembership.pangenome_id == pangenome_id,
            Gene.strand == "-",
            Gene.start_position > 500,
            Gene.length_nt > 900,
            # ⚠ A locus SOME collection genome lacks, so the "no gene here" case below exists at the
            # same locus — on kp the first minus-strand gene sat at a locus all 100 genomes carry.
            Locus.member_genome_count
            < select(Pangenome.genome_count).where(Pangenome.pangenome_id == pangenome_id).scalar_subquery(),
        )
        .order_by(Gene.genome_id, Gene.flat_index)
        .limit(1)
    ).one()
    members = (
        select(GeneLocusMembership.genome_id)
        .join(Locus, Locus.locus_id == GeneLocusMembership.locus_id)
        .where(Locus.node_label == minus.node_label, Locus.pangenome_id == pangenome_id)
    )
    # ⛔ A genome of THIS COLLECTION that lacks a gene here — not merely one of the species. The
    # database holds more genomes per species than the pangenome modelled, and one outside the
    # collection is a different claim (404, "not in this catalogue"), which the recorded "no gene here"
    # answer would then have silently been.
    collection = select(Pangenome.genome_collection_id).where(Pangenome.pangenome_id == pangenome_id)
    absent = session.execute(
        select(Genome.sample_id)
        .join(GenomeCollectionMembership, GenomeCollectionMembership.genome_id == Genome.genome_id)
        .where(
            GenomeCollectionMembership.genome_collection_id == collection.scalar_subquery(),
            Genome.genome_id.not_in(members),
        )
        .order_by(Genome.sample_id)
        .limit(1)
    ).scalar_one()
    return {"minus_strand": (minus.sample_id, minus.node_label), "no_gene_here": (absent, minus.node_label)}


def record_species(client, session: Session, species_key: str) -> dict:
    """Every response one species' fixture carries, each asserted to have been answered."""
    pangenome_id = session.execute(
        select(PathogenSpecies.default_pangenome_id).where(PathogenSpecies.species_key == species_key)
    ).scalar_one()

    def get(path, **query):
        reply = client.get(path, query_string=query)
        assert reply.status_code == 200, (path, reply.status_code)
        return reply

    out: dict = {"recorded_from": "the API test client", "species_key": species_key, "loci": {}, "sequences": {}}
    # ⭐ The species response too: the locus card reads its null baselines and measurable-locus counts
    # from it and everything else from the locus endpoint, so only the recorded pair can show the two
    # disagreeing.
    out["species"] = get(f"/api/v1/species/{species_key}").get_json()
    for name, label in choose_cases(session, pangenome_id).items():
        # ⛔ The function block for EVERY case: a second endpoint describing the same locus, and only
        # a recorded pair can show the two disagreeing.
        out["loci"][name] = {
            "label": label,
            "response": get(f"/api/v1/species/{species_key}/loci/{label}").get_json(),
            "function": get(f"/api/v1/species/{species_key}/loci/{label}/function").get_json(),
        }
    for name, (sample, label) in choose_sequences(session, pangenome_id).items():
        response = get(f"/api/v1/species/{species_key}/genomes/{sample}/loci/{label}/sequence").get_json()
        out["sequences"][name] = {"sample_id": sample, "locus_label": label, "response": response}
    out["residuals"] = get(f"/api/v1/species/{species_key}/audit/residual-loci").get_json()
    # The anchor picker's list, cut short on purpose so the `truncated` path is recorded too.
    out["genomes"] = get(f"/api/v1/species/{species_key}/genomes", limit=5).get_json()
    out["search"] = {"ligase": get(f"/api/v1/species/{species_key}/search", q="ligase", limit=40).get_json()}
    return out


def main() -> None:
    """Write both species' fixtures and print what each case exercises."""
    application = create_application(Configuration.from_environment())
    client = application.test_client()
    engine = application.extensions["bacatlas_database"].engine
    with Session(engine) as session:
        for species_key in SPECIES_KEYS:
            recorded = record_species(client, session, species_key)
            target = FIXTURES / f"api_responses_{species_key}.json"
            target.write_text(json.dumps(recorded, indent=1, sort_keys=True) + "\n")
            print(f"wrote {target.name}")
            for name, entry in recorded["loci"].items():
                arrangements = entry["response"]["arrangements"]
                print(
                    f"  {name:14s} locus {entry['label']:>6s} listed={len(arrangements['listed'])} "
                    f"total={arrangements['total']} past_cap={arrangements['members_in_arrangements_not_listed']} "
                    f"no_window={arrangements['members_without_a_neighbourhood']}"
                )
            for name, entry in recorded["sequences"].items():
                print(f"  sequence {name:14s} {entry['sample_id']} @ {entry['locus_label']} -> "
                      f"{len(entry['response']['genes'])} gene(s)")


if __name__ == "__main__":
    main()
