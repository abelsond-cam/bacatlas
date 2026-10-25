"""``GET /health`` — liveness, and whether the database behind it is actually usable.

⚠ A health check that only proves the process is running is worse than none: it goes green
while every request 500s on a database that is unreachable, or reachable and empty. This one
reports the connection, the server version and whether the schema has been migrated, and it
says which of those failed rather than returning a bare 503.
"""

from __future__ import annotations

from flask import Blueprint, current_app, jsonify
from sqlalchemy import text
from sqlalchemy.exc import SQLAlchemyError

from bacatlas_backend.sequence_source.provider_factory import SEQUENCE_PROVIDER_KEY

health_blueprint = Blueprint("health", __name__)

#: Present once the first migration has run. Until then the API is up and has nothing to serve,
#: which is a distinct state from "the database is down" and must read as one.
_SCHEMA_PROBE = text("SELECT to_regclass('public.pangenome') IS NOT NULL AS schema_present")


@health_blueprint.get("/health")
def report_health():
    """Report process liveness and database reachability. 200 healthy, 503 otherwise."""
    database = current_app.extensions["bacatlas_database"]
    report: dict[str, object] = {
        "service": "bacatlas-backend",
        "deployment_profile": current_app.config["BACATLAS"].deployment_profile,
        "database_reachable": False,
        "schema_present": False,
    }
    try:
        with database.session() as session:
            report["database_server_version"] = session.execute(text("SHOW server_version")).scalar_one()
            report["database_reachable"] = True
            report["schema_present"] = bool(session.execute(_SCHEMA_PROBE).scalar_one())
    except SQLAlchemyError as error:
        # The class name and message, not a traceback: this endpoint is public.
        report["database_error"] = f"{type(error).__name__}: {error}"
        return jsonify(report), 503

    report["sequence_source"] = _sequence_source_report()

    if not report["schema_present"]:
        report["note"] = "database reachable but not migrated — run `alembic upgrade head`"
    return jsonify(report), 200


def _sequence_source_report() -> dict[str, object]:
    """⭐ **Where the Sequence tab's bases come from, answerable from a browser.**

    ⛔ Added because the only other evidence was a `logger.info`, and gunicorn configures its own
    loggers rather than the root one — so a module logger sits at WARNING and that line is
    swallowed. Someone asked to confirm they were testing the ATB pull and, as shipped, could not:
    a page that reads correctly looks identical whether the bytes came off the local tree or the
    wire. This states which, and `cached_assemblies` moves when a fetch really happens.

    ⚠ Deliberately NOT gated differently from the rest of `/health`: it names no path a reader
    could not already infer from the deployment, carries no credential, and the whole origin sits
    behind basic auth anyway.
    """
    configuration = current_app.config["BACATLAS"]
    provider = current_app.extensions.get(SEQUENCE_PROVIDER_KEY)
    gff_root = configuration.artifact_roots.get("gff")
    cached: object = "no cache directory configured"
    if configuration.atb_cache_directory is not None:
        try:
            cached = sum(1 for _ in configuration.atb_cache_directory.glob("*.fa.gz"))
        except OSError as error:
            cached = f"unreadable: {type(error).__name__}"
    return {
        # ⚠ "held tree present" is about CONFIGURATION, not about any one genome: a genome absent
        # from a configured tree still streams. It answers "could a local file have served this?"
        "held_annotation_tree": str(gff_root) if gff_root else None,
        "assembly_streaming": bool(provider is not None and provider.can_stream),
        "cached_assemblies": cached,
    }
