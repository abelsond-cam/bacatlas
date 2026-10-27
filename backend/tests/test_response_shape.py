"""Every cacheable response's SHAPE, written down — so a field cannot move without a bump.

⛔⛔ **Nothing in this repo asserted the shape of any JSON response, on either side.** Swept before
this file was written: `backend/tests/` contains no `set(payload)`, no `payload.keys()`, no
`sorted(payload)` — not one key-set assertion on any served body. The frontend has exactly one
(`apiContract.test.ts:163`, over the locus detail) and nothing equivalent for the shell. So the six
field moves of 2026-10-27 — the picker's `model` block replaced by `cohort` + `model_label`, the
shell's `exclusivity_form` gone and `model_label` + `cohort` arrived, and its top-level `model`,
`steps` and `provenance_rows` deleted — would each have been invisible to all 1,400-odd tests.

⛔ **And the one test that LOOKS like it guards `RESPONSE_SHAPE_VERSION` cannot.**
`test_every_catalogue_response_is_cacheable_and_tagged_by_its_pangenome_AND_its_shape` imports the
constant and asserts the header agrees with it — true by construction, for any value. Its own
docstring describes the 2026-10-25 deploy the bump exists to prevent, and it would not have caught
a forgotten one. **The literal below is what closes that**: the version is written out as a number
here, so changing the code's constant without editing this file turns a test red, and changing a
response's shape without editing this file turns a *different* test red.

⭐ **The walk is RECURSIVE, because the regression that created this mechanism was a nested field.**
`genome_source_groups` was added one level down, inside `pangenome`. A top-level key set would not
have seen it, and would not see any of this stage's changes either.

⚠ **Data-keyed maps are allowlisted, not walked.** In `prevalence_census`, `audit_headline`,
`genome_source_groups` and their kin the KEYS are values — band names, metric names, source groups
— so walking into them would pin data as though it were shape, and the test would go red whenever a
catalogue gained a metric. They are asserted to be *present and to be objects*, and their contents
are the business of the tests that measure them.
"""

from __future__ import annotations

import os

import pytest
from sqlalchemy import create_engine, func, select
from sqlalchemy.orm import Session

from bacatlas_backend.application_factory import Configuration, create_application
from bacatlas_backend.models.pathogen_species import PathogenSpecies

#: ⛔ A LITERAL, deliberately not imported. Importing the constant is what makes the existing
#: assertion vacuous. This number has to be edited by hand when a shape moves, and that edit is
#: the whole mechanism.
EXPECTED_RESPONSE_SHAPE_VERSION = 4

#: Key paths whose CHILDREN are data rather than shape — walked to, never walked into.
#: ⚠ Each is a map whose keys are values: band names, metric names, anchor source groups. Walking
#: in would make this test fail the day a catalogue gains a metric, which is not a shape change.
DATA_KEYED = frozenset({
    "pangenome.genome_source_groups",
    "pangenome.omitted_sections",
    "prevalence_census",
    "prevalence_gene_census",
    "audit_headline",
})


@pytest.fixture(scope="module")
def client():
    url = os.environ.get("BACATLAS_DATABASE_URL")
    if not url:
        pytest.skip("BACATLAS_DATABASE_URL is not set — the shape tests need a loaded database")
    engine = create_engine(url, future=True)
    with Session(engine) as session:
        published = session.execute(
            select(func.count())
            .select_from(PathogenSpecies)
            .where(PathogenSpecies.default_pangenome_id.is_not(None))
        ).scalar_one()
    if published < 1:
        pytest.skip(f"nothing is published in {url}; run the loader with --publish")
    return create_application(Configuration(database_url=url)).test_client()


def _key_paths(value, prefix: str = "") -> set[str]:
    """Every dotted key path in a JSON body, with list elements collapsed to their first.

    ⚠ Lists collapse because their LENGTH is data — two catalogues today, three tomorrow — while
    the shape of an element is not. The first element stands for all of them, which is sound only
    because every list here is homogeneous; `test_every_element_of_a_list_has_the_SAME_shape`
    below is what stops that assumption rotting.
    """
    paths: set[str] = set()
    if isinstance(value, dict):
        for key, child in value.items():
            path = f"{prefix}.{key}" if prefix else key
            paths.add(path)
            if path not in DATA_KEYED:
                paths |= _key_paths(child, path)
    elif isinstance(value, list) and value:
        paths |= _key_paths(value[0], prefix)
    return paths


#: ⛔ The shape of record for each cacheable list route, written out. An addition, a removal or a
#: rename anywhere in these bodies fails here and nowhere else.
SPECIES_LIST = {
    "species",
    "species.key", "species.scientific_name", "species.ncbi_taxonomy_id",
    "species.published", "species.genome_count", "species.gene_count",
    "species.locus_count", "species.catalogue_key",
}  # fmt: skip

CATALOGUE_LIST = {
    "catalogues",
    "catalogues.key",
    "catalogues.species", "catalogues.species.key", "catalogues.species.scientific_name",
    # ⭐ The axis. `model: {key,label,step_count,exclusivity_form}` stood here until 2026-10-27.
    "catalogues.cohort", "catalogues.cohort.key", "catalogues.cohort.label",
    "catalogues.cohort.genome_count",
    "catalogues.model_label",
    "catalogues.is_default", "catalogues.genome_count", "catalogues.gene_count",
    "catalogues.locus_count", "catalogues.run_id",
}  # fmt: skip


def test_the_species_list_has_exactly_this_shape(client):
    assert _key_paths(client.get("/api/v1/species").get_json()) == SPECIES_LIST


def test_the_catalogue_list_has_exactly_this_shape(client):
    """⛔ The picker's own route, which had NO test of any kind before this one.

    No HTTP test in the backend suite, no entry in the deploy smoke test's sixteen checks, no
    recorded fixture, and the component that renders it is mounted by nothing. It could have 500'd
    in production with every gate green.
    """
    assert _key_paths(client.get("/api/v1/catalogues").get_json()) == CATALOGUE_LIST


def test_the_shell_carries_the_cohort_and_the_model_label_and_NOT_the_model_block(client):
    """The shell is too large to pin whole; these are the six fields this stage moved.

    ⚠ Asserted in BOTH directions. A test that only checks the new fields arrived would pass with
    the old ones still being served beside them — which is the state a half-applied change leaves,
    and the one a `RESPONSE_SHAPE_VERSION` bump is least able to describe.
    """
    paths = _key_paths(client.get("/api/v1/species/ecoli").get_json())

    assert "pangenome.model_label" in paths
    assert {"pangenome.cohort.key", "pangenome.cohort.label",
            "pangenome.cohort.genome_count"} <= paths  # fmt: skip

    # ⚠ `model`, `steps` and `provenance_rows` were TOP-LEVEL keys of the shell, siblings of
    # `pangenome` and not children of it.
    gone = {"pangenome.exclusivity_form", "model", "steps", "provenance_rows"}
    assert not (gone & paths), f"still served: {sorted(gone & paths)}"


def test_every_element_of_a_list_has_the_SAME_shape(client):
    """⚠ `_key_paths` collapses a list to its first element. This is what makes that honest.

    Without it, a catalogue whose second row carried an extra field — or lacked one — would be
    invisible to every assertion above, and the shape tests would read as exhaustive while
    examining one row in N.
    """
    for route, block in (("/api/v1/species", "species"), ("/api/v1/catalogues", "catalogues")):
        rows = client.get(route).get_json()[block]
        assert len(rows) >= 2, f"{route} returned {len(rows)} row(s) — too few to compare shapes"
        shapes = {frozenset(_key_paths(row)) for row in rows}
        assert len(shapes) == 1, (
            f"{route} returns rows of {len(shapes)} different shapes; the first-element collapse "
            f"in `_key_paths` is then measuring one row and reporting on all of them"
        )


def test_the_shape_version_is_the_number_this_file_says_it_is():
    """⛔ The literal is the mechanism. Importing it, as the older test does, asserts nothing.

    `test_every_catalogue_response_is_cacheable_and_tagged_by_its_pangenome_AND_its_shape` compares
    the header against the imported constant, so it is true for any value and cannot detect a
    forgotten bump. This one fails until somebody edits this file, and the shape assertions above
    fail until somebody edits this file, so the two edits happen together or neither does.
    """
    from bacatlas_backend.api.blueprint_species import RESPONSE_SHAPE_VERSION

    assert RESPONSE_SHAPE_VERSION == EXPECTED_RESPONSE_SHAPE_VERSION, (
        "the served shape version moved. If a response's shape changed, update the key-path sets "
        "in this file in the same commit; if it did not, the bump is spurious."
    )


def test_the_two_LIST_routes_now_carry_a_validator_and_it_tracks_the_BODY(client):
    """⛔ They carried no ETag at all — worse than the incident the mechanism was built for.

    `_immutable(payload, None)` skipped the header entirely, so `/species` and `/catalogues` went
    out with `max-age=86400` and nothing to revalidate against. ⚠ And a shape-only tag would have
    been actively wrong here: these bodies change when a catalogue is published or hidden, or when
    a cohort is relabelled — database edits with no deploy and no shape change. The tag is a digest
    of the body for that reason.
    """
    tags = {}
    for route in ("/api/v1/species", "/api/v1/catalogues"):
        response = client.get(route)
        assert response.headers["Cache-Control"] == "public, max-age=86400"
        tag = response.headers["ETag"]
        digest, _, shape = tag.strip('"').partition("-")
        assert shape == f"v{EXPECTED_RESPONSE_SHAPE_VERSION}", tag
        assert len(digest) == 16 and all(c in "0123456789abcdef" for c in digest), tag
        tags[route] = tag

    # ⛔ Two routes, two bodies, two tags. One shared string would make a hand-written
    # `If-None-Match` against the wrong route look valid, and would say nothing about WHICH route's
    # content had moved.
    assert tags["/api/v1/species"] != tags["/api/v1/catalogues"]

    # ⚠ And it is stable: a digest that moved per request would defeat the point entirely.
    assert client.get("/api/v1/catalogues").headers["ETag"] == tags["/api/v1/catalogues"]
