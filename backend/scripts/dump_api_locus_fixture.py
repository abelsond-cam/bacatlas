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
- `predicted_function` — several annotation-transfer v0.1 predictions (a COG category and a GO class
  with its lineage), at least one from a donor in ANOTHER catalogue, so the predicted-function card
  renders on real bytes rather than only on a hand-built object. ⚠ Needs the served predictions
  loaded (`ingest --stage served-annotations`).
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
from bacatlas_backend.gene_blocks import GeneSpanTable, decode_int32_vector
from bacatlas_backend.models.genome import Genome
from bacatlas_backend.models.genome_collection import GenomeCollectionMembership
from bacatlas_backend.models.genome_genes import GenomeGenes
from bacatlas_backend.models.genome_locus_assignment import GenomeLocusAssignment
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
        # ⭐ A locus annotation-transfer v0.1 makes SEVERAL predictions for — a COG category, and a GO
        # class with the lineage it implies — at least one copied from a donor in ANOTHER catalogue, so
        # the predicted-function card renders every part of itself on real bytes: the stated class and
        # its "within", a letter named on the page, and a cross-catalogue link. ⚠ Needs the catalogue's
        # predictions loaded (`ingest --stage served-annotations`); a database without them cannot
        # build this case, and `_one` says so.
        "predicted_function": _one(
            session,
            text("""
                select focal.node_label
                  from locus focal
                  join served_locus_annotation s on s.locus_id = focal.locus_id
                  join locus donor on donor.locus_id = s.donor_locus_id
                 where focal.pangenome_id = :pangenome_id
                 group by focal.locus_id, focal.node_label, focal.member_gene_count, focal.catalogue_ordinal
                having bool_or(s.annotation_kind = 'COG_FUNCTIONAL_CATEGORY')
                   and bool_or(s.annotation_kind = 'GENE_ONTOLOGY_SLIM' and cardinality(s.term_values) >= 2)
                   and bool_or(donor.pangenome_id <> focal.pangenome_id)
                 order by count(*) desc, focal.member_gene_count desc, focal.catalogue_ordinal
                 limit 1
            """).bindparams(pangenome_id=pangenome_id),
            "predicted_function",
        ),
    }


def choose_sequences(session: Session, pangenome_id: int) -> dict[str, tuple[str, str]]:
    """A minus-strand gene well inside its contig, and a genome with no gene at that locus.

    Read from the genome-grain blocks: each genome's `gene_locus_ordinals` against its `gene_spans`,
    in `(genome_id, flat_index)` order, the first gene that qualifies.
    """
    genome_count = session.execute(
        select(Pangenome.genome_count).where(Pangenome.pangenome_id == pangenome_id)
    ).scalar_one()
    # ⚠ A locus SOME collection genome lacks, so the "no gene here" case below exists at the same
    # locus — on kp the first minus-strand gene sat at a locus all 100 genomes carry.
    label_of_ordinal = dict(
        session.execute(
            select(Locus.catalogue_ordinal, Locus.node_label).where(
                Locus.pangenome_id == pangenome_id, Locus.member_genome_count < genome_count
            )
        ).all()
    )
    blocks = session.execute(
        select(
            Genome.sample_id,
            GenomeGenes.gene_count,
            GenomeGenes.gene_spans,
            GenomeLocusAssignment.gene_locus_ordinals,
        )
        .join(GenomeGenes, GenomeGenes.genome_id == GenomeLocusAssignment.genome_id)
        .join(Genome, Genome.genome_id == GenomeLocusAssignment.genome_id)
        .where(GenomeLocusAssignment.pangenome_id == pangenome_id)
        .order_by(GenomeLocusAssignment.genome_id)
    )
    minus = None
    for sample_id, gene_count, spans, ordinal_block in blocks:
        table = GeneSpanTable(spans, gene_count=gene_count)
        ordinals = decode_int32_vector(ordinal_block, expected_count=gene_count)
        flat_index = next(
            (
                index
                for index, ordinal in enumerate(ordinals)
                if ordinal in label_of_ordinal
                and table.is_minus_strand(index)
                and table.start_position(index) > 500
                and table.length_nt(index) > 900
            ),
            None,
        )
        if flat_index is not None:
            minus = (sample_id, label_of_ordinal[ordinals[flat_index]])
            break
    if minus is None:
        raise SystemExit("no minus-strand gene qualifies for the `minus_strand` sequence case")
    present = session.execute(
        select(Locus.present_genome_ids).where(Locus.node_label == minus[1], Locus.pangenome_id == pangenome_id)
    ).scalar_one()
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
            Genome.genome_id.not_in(present),
        )
        .order_by(Genome.sample_id)
        .limit(1)
    ).scalar_one()
    return {"minus_strand": minus, "no_gene_here": (absent, minus[1])}


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
