# syntax=docker/dockerfile:1
#
# The BacAtlas API — gunicorn over `create_application()`, with the SERVING dependencies only.
#
# ⛔ No `ingest` extra, and never `nuna`. `nuna` is a private repo the server has no access to and
# does not need: only `ingest/` imports it, ingest runs on a machine that already has the science
# package, and what it produces reaches the server as a `pg_dump`. Installing `ingest` here would
# add pandas, pyarrow and scikit-learn to an image that never calls them.
#
# Built by `compose.yaml` at the repo root; the context is `backend/`.

FROM python:3.13-slim

ENV UV_PROJECT_ENVIRONMENT=/opt/venv \
    UV_PYTHON_DOWNLOADS=never \
    UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy \
    PATH=/opt/venv/bin:$PATH \
    PYTHONUNBUFFERED=1

WORKDIR /srv/bacatlas/backend

# Dependencies first, from the lock alone, so a source edit does not reinstall them.
# ⛔ `--locked`, not `--frozen`: a `pyproject.toml` edited without re-running `uv lock` FAILS the
# build here. `--frozen` would install the stale lock and ship a server that resolved something
# other than what the tests ran against — the lock exists so the med school box gets what we tested.
#
# uv (pinned) is bind-mounted for the RUN only and its cache lives in a BuildKit cache mount, so
# neither ships: together they were ~100 MB of layers in an image that never runs either.
COPY pyproject.toml uv.lock ./
RUN --mount=from=ghcr.io/astral-sh/uv:0.11.8,source=/uv,target=/bin/uv \
    --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-dev --no-install-project --extra serve

# `alembic.ini` + the package (which carries `migrations/`), so that
# `docker compose run --rm api alembic upgrade head` works against a fresh database.
COPY alembic.ini ./
COPY bacatlas_backend ./bacatlas_backend
RUN --mount=from=ghcr.io/astral-sh/uv:0.11.8,source=/uv,target=/bin/uv \
    --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-dev --extra serve

# Not root: the process only ever reads — the database over a socket, the GFF tree through a
# read-only bind mount. ⚠ On a Linux host that tree must therefore be readable by uid 10001 (or
# world-readable); a permission error surfaces as the sequence endpoint's named 503, not at start-up.
RUN useradd --system --uid 10001 --no-create-home bacatlas
# ⛔ The assembly cache directory is created HERE, owned by the app user, and that is not tidiness.
# Docker seeds a fresh named volume from whatever the image has at the mount path, ownership
# included — so if the path does not exist in the image the volume is created root-owned and uid
# 10001 cannot write to it. The store degrades to "no cache" on an unwritable directory rather than
# failing, so the symptom would be a silent re-fetch on every single request, not an error.
RUN install -d -o bacatlas -g bacatlas -m 0755 /srv/bacatlas/seqcache
USER bacatlas

EXPOSE 8000

# `create_application()` with no argument reads `BACATLAS_*` from the environment and raises at
# start-up on anything missing, so a misconfigured container dies here, visibly, not on a request.
# Workers come from `WEB_CONCURRENCY`, which gunicorn reads itself; `compose.yaml` sets it.
# `--no-control-socket`: gunicorn 26's runtime-control socket lives in $HOME, which this user does not
# have, so it logs an ERROR at every start — and the container is managed by Docker, not by `gunicornc`.
CMD ["gunicorn", \
     "--bind", "0.0.0.0:8000", \
     "--access-logfile", "-", \
     "--no-control-socket", \
     "bacatlas_backend.application_factory:create_application()"]
