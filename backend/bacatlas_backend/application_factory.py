"""``create_application(configuration) -> Flask`` — registers blueprints and nothing else.

The factory decides nothing. Which endpoints exist is the blueprint list below; what they
return is the services layer; where the data is, is `Configuration`. Keeping it that thin is
what lets a test build an app against a throwaway database in two lines.
"""

from __future__ import annotations

import logging

from flask import Flask, jsonify, request
from werkzeug.exceptions import HTTPException, InternalServerError

from bacatlas_backend.api.blueprint_health import health_blueprint
from bacatlas_backend.api.blueprint_species import species_blueprint
from bacatlas_backend.configuration import Configuration
from bacatlas_backend.database import Database
from bacatlas_backend.sequence_source.provider_factory import (
    SEQUENCE_PROVIDER_KEY,
    build_sequence_provider,
)

logger = logging.getLogger(__name__)

#: Every blueprint the application serves, in registration order. One line per resource.
BLUEPRINTS = (health_blueprint, species_blueprint)

#: All API routes live under this prefix. The version is in the path rather than a header so a
#: breaking change can run alongside its predecessor instead of replacing it.
API_PREFIX = "/api/v1"


def create_application(configuration: Configuration | None = None) -> Flask:
    """Build the app. Pass a `Configuration` in tests; omit it to read the environment."""
    configuration = configuration or Configuration.from_environment()
    application = Flask(__name__)
    application.config["BACATLAS"] = configuration
    application.extensions["bacatlas_database"] = Database(configuration)
    application.extensions[SEQUENCE_PROVIDER_KEY] = build_sequence_provider(configuration)

    for blueprint in BLUEPRINTS:
        application.register_blueprint(blueprint, url_prefix=API_PREFIX)

    # Read-only service, immutable responses: JSON key order carries no meaning and sorting it
    # makes two responses byte-comparable, which the parity suites rely on.
    application.json.sort_keys = True
    _register_json_errors(application)
    return application


def _register_json_errors(application: Flask) -> None:
    """Every error under `/api/` is JSON, in the same `{error, detail}` shape as a named 404.

    ⚠ Without this, a route that does not exist answered with Flask's HTML 404 and an unhandled
    exception — e.g. a database that was never migrated — with its HTML 500, and the client reported
    "Unexpected token '<'" instead of the status. The routes' own named 404s are untouched; this is
    only what happens when nothing more specific did.

    ⛔ A 500 NEVER carries the exception's text: it can name tables, paths and SQL, and the reader of
    this response is a browser, not whoever runs the server. The exception goes to the log.
    """

    def is_api() -> bool:
        return request.path.startswith(API_PREFIX + "/") or request.path == API_PREFIX

    @application.errorhandler(HTTPException)
    def http_error(error: HTTPException):
        if not is_api():
            return error
        return jsonify({"error": (error.name or "error").lower().replace(" ", "_"), "detail": error.description}), (
            error.code or 500
        )

    @application.errorhandler(InternalServerError)
    def server_error(error: InternalServerError):
        if not is_api():
            return error
        original = getattr(error, "original_exception", None)
        if original is not None:
            logging.getLogger(__name__).error("unhandled error on %s", request.path, exc_info=original)
        return jsonify({"error": "server_error", "detail": "the server could not answer this request"}), 500
