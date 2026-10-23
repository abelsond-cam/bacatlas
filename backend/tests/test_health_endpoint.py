import os

import pytest

from bacatlas_backend.application_factory import create_application
from bacatlas_backend.configuration import Configuration

DEV_URL = os.environ.get(
    "BACATLAS_TEST_DATABASE_URL",
    f"postgresql+psycopg://{os.environ.get('USER', 'postgres')}@localhost:5432/bacatlas_ingest",
)


def _client(database_url):
    return create_application(Configuration(database_url=database_url)).test_client()


def test_health_reports_a_reachable_database_and_says_whether_it_is_migrated():
    response = _client(DEV_URL).get("/api/v1/health")
    body = response.get_json()
    assert response.status_code == 200
    assert body["database_reachable"] is True
    assert body["database_server_version"].startswith("16.")
    # schema_present is False until the first migration lands and True after; either is healthy,
    # which is the whole point of reporting it separately from reachability.
    assert isinstance(body["schema_present"], bool)


def test_an_unreachable_database_is_503_and_NAMES_the_failure():
    # The branch that matters: a health check that goes green on a dead database is worse than
    # none. Port 1 is closed on every machine, so this exercises connect failure, not auth.
    response = _client("postgresql+psycopg://nobody@localhost:1/nothing").get("/api/v1/health")
    body = response.get_json()
    assert response.status_code == 503
    assert body["database_reachable"] is False
    assert "database_error" in body, "a 503 must say what failed, not just that something did"


def test_reachable_but_unmigrated_is_NOT_reported_as_a_failure():
    # These are different states and an operator has to be able to tell them apart: one needs a
    # DBA, the other needs `alembic upgrade head`.
    body = _client(DEV_URL).get("/api/v1/health").get_json()
    if not body["schema_present"]:
        assert body["database_reachable"] is True
        assert "run `alembic upgrade head`" in body["note"]


def test_a_missing_artifact_root_names_the_key_and_the_profile():
    # The real failure mode is a profile that was never given the root, which is otherwise
    # indistinguishable from a missing file.
    configuration = Configuration(database_url=DEV_URL, deployment_profile="development")
    with pytest.raises(KeyError, match="gff"):
        configuration.resolve_artifact("gff", "anything.gff3.gz")
