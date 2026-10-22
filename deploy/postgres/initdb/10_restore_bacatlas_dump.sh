#!/usr/bin/env bash
# Restore the BacAtlas database from `/restore/bacatlas.dump`, or start empty if there is none.
#
# Run by the official postgres image's entrypoint, and ⚠ ONLY on the first start against an EMPTY
# data volume — that is the image's rule, not ours. A new dump is therefore never picked up by a
# restart: a deploy is `docker compose down -v` (drops the volume) and `docker compose up -d`.
#
# ⛔ If this script fails, the image has ALREADY run `initdb`, so the volume is no longer empty: the
# container restarts (`restart: unless-stopped`), SKIPS this script, and comes up healthy on whatever
# half-restored state it left — measured with a truncated dump: healthy within a second, 0 tables.
# So a failure leaves a marker in the data directory, and the compose healthcheck refuses to report
# healthy while it exists — the api never starts, and `up --wait` fails instead of serving an empty
# catalogue. The only way out is the right one: `docker compose down -v`, fix the dump, `up` again.
#
# ⚠ Keep this file EXECUTABLE. The entrypoint runs an executable `.sh` in its own process but SOURCES
# a non-executable one into the entrypoint's shell, where `set -u` and any early exit would leak.
set -Eeuo pipefail

dump=/restore/bacatlas.dump
failed_marker="$PGDATA/BACATLAS_RESTORE_FAILED"

trap 'echo "bacatlas: RESTORE FAILED — the database is incomplete. Run: docker compose down -v" >&2;
      echo "restore of $dump failed at $(date -u +%FT%TZ); docker compose down -v" > "$failed_marker"' ERR

# `pg_trgm` whether or not a dump follows. The dump creates it too (`IF NOT EXISTS`), but a database
# that starts EMPTY needs it before `alembic upgrade head` can build `ix_locus__search_text_trigram`:
# the migrations do not create it — `backend/README.md` has it as a provisioning step — and without
# it the upgrade dies on `operator class "gin_trgm_ops" does not exist`.
psql -v ON_ERROR_STOP=1 --username "$POSTGRES_USER" --dbname "$POSTGRES_DB" \
  -c 'CREATE EXTENSION IF NOT EXISTS pg_trgm'

if [ ! -f "$dump" ]; then
  echo "bacatlas: no dump at $dump — the database starts EMPTY (no schema). See README.md."
else
  echo "bacatlas: restoring $(ls -lh "$dump" | awk '{print $5}') from $dump into $POSTGRES_DB"
  started=$(date +%s)
  # `--no-owner --no-privileges`: the dump was taken on the Mac as that user, who does not exist
  # here; every object is owned by the role that restores it instead.
  pg_restore --username "$POSTGRES_USER" --dbname "$POSTGRES_DB" \
    --no-owner --no-privileges --exit-on-error \
    --jobs "${BACATLAS_RESTORE_JOBS:-4}" \
    "$dump"
  # ⚠ pg_restore restores rows, not planner statistics. Until autovacuum gets round to it, every
  # query is planned against a table the planner believes is empty.
  psql -v ON_ERROR_STOP=1 --username "$POSTGRES_USER" --dbname "$POSTGRES_DB" -c 'ANALYZE'
  echo "bacatlas: restore + ANALYZE finished in $(( $(date +%s) - started )) s"
fi
