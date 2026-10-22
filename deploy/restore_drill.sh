#!/usr/bin/env bash
# The restore drill — the S3 exit criterion, as one command the owner runs unaided.
#
#   deploy/restore_drill.sh --yes-delete-the-database
#
# It DELETES the running database volume, restores it from the dump, times that, then checks the
# stack end to end THROUGH THE PUBLISHED PORT — the path a reader's browser takes — and prints one
# PASS/FAIL line per check and a summary. Exit status 0 means every check passed.
#
# ⛔ Destructive by design: `docker compose down -v` is how a new dump is loaded (the postgres image
# restores only into an EMPTY volume), so the drill IS the upgrade procedure. It refuses to run
# without the flag above, and it refuses before deleting anything if the dump is missing — a drill
# that deletes the live database and then finds no dump has made an outage out of a rehearsal.
#
# Needs: Docker Engine 25+ with Compose v2, curl, and `.env` beside compose.yaml (see .env.example).
# ⚠ Nothing else on the host — no Python: JSON is read with the api container's own interpreter.
set -Eeuo pipefail

cd "$(dirname "$0")/.."

if [ "${1:-}" != "--yes-delete-the-database" ]; then
  echo "usage: deploy/restore_drill.sh --yes-delete-the-database"
  echo "  deletes the database volume, restores it from the dump, and checks the stack end to end."
  exit 2
fi

failures=0
pass() { printf '  PASS  %s\n' "$1"; }
fail() { printf '  FAIL  %s\n' "$1"; failures=$((failures + 1)); }

# ── preflight: nothing is deleted until all of this holds ─────────────────────────────────────────
[ -f .env ] || { echo "no .env beside compose.yaml — copy .env.example and fill it in"; exit 2; }
set -a; . ./.env; set +a
dump_dir="${BACATLAS_HOST_DUMP_DIR:-./deploy/dump}"
dump="$dump_dir/bacatlas.dump"
[ -f "$dump" ] || { echo "no dump at $dump — refusing to delete a database it cannot replace"; exit 2; }
command -v curl >/dev/null || { echo "curl is required"; exit 2; }
docker compose version >/dev/null 2>&1 || { echo "Docker Compose v2 is required"; exit 2; }
engine=$(docker version --format '{{.Server.Version}}' 2>/dev/null || echo 0)
[ "${engine%%.*}" -ge 25 ] 2>/dev/null || { echo "Docker Engine 25+ is required (found $engine)"; exit 2; }

# ⛔ THE DUMP MUST MATCH THE CODE, and this is checked BEFORE anything is deleted.
# A dump four months behind the code restored cleanly, reported `schema_present: true`, passed the
# health check — and then 500'd EVERY locus request on `relation "enzyme_class" does not exist`
# (2026-10-22). Nothing in the stack noticed, because a restore cannot tell that the schema it just
# laid down is older than the code reading it. So read the vintage off the dump itself.
#
# Both halves need nothing new on the host:
#   - the dump's revision: `pg_restore` to stdout, in the postgres image the stack already uses. The
#     archive is never restored, only read.
#   - the code's revision: the migration that no other migration names as its parent. No alembic, no
#     database, just the files on disk.
versions=backend/bacatlas_backend/migrations/versions
dump_rev=$(docker run --rm -v "$(cd "$dump_dir" && pwd)":/r:ro postgres:16 \
             pg_restore -f - --data-only --table=alembic_version "/r/$(basename "$dump")" 2>/dev/null \
           | grep -Eoi '\b[0-9a-f]{12}\b' | tail -1)
code_rev=$(comm -23 \
  <(grep -h '^revision' "$versions"/*.py | grep -Eo "'[0-9a-f]{12}'" | tr -d "'" | sort -u) \
  <(grep -h '^down_revision' "$versions"/*.py | grep -Eo "'[0-9a-f]{12}'" | tr -d "'" | sort -u))
[ -n "$dump_rev" ] || { echo "could not read alembic_version from $dump — is it a \`pg_dump -Fc\` archive?"; exit 2; }
[ "$(printf '%s' "$code_rev" | grep -c .)" = "1" ] || {
  echo "the migrations directory has no single head, so there is nothing to match against:"
  printf '  %s\n' $code_rev; exit 2; }
if [ "$dump_rev" != "$code_rev" ]; then
  echo "⛔ THE DUMP IS THE WRONG VINTAGE — refusing, and nothing has been deleted."
  echo "     the dump carries  $dump_rev"
  echo "     this code needs   $code_rev"
  echo "   Restoring it would leave a database that answers /health happily and 500s on real queries."
  echo "   Take a fresh one on the machine that ran the ingest, then re-run this drill:"
  echo "     pg_dump -Fc -d bacatlas_dev -f $dump"
  exit 2
fi

base="${BACATLAS_PUBLIC_BASE:-/}"
host="${BACATLAS_HTTP_BIND:-127.0.0.1}"
[ "$host" = "0.0.0.0" ] && host=127.0.0.1
url="http://$host:${BACATLAS_HTTP_PORT:-8080}${base}"
api="${url}api/v1"

echo "BacAtlas restore drill — $(date -u +%FT%TZ)"
checksum=$( (sha256sum "$dump" 2>/dev/null || shasum -a 256 "$dump") | cut -c1-16)
echo "  dump    $dump ($(du -h "$dump" | cut -f1), sha256 ${checksum}…)"
echo "  schema  $dump_rev (matches the code)"
echo "  serves  $url"
echo "  engine  Docker $engine"
echo

# ── the drill ─────────────────────────────────────────────────────────────────────────────────────
echo "deleting the database volume and restoring from the dump…"
docker compose down -v --remove-orphans >/dev/null 2>&1
started=$(date +%s)
if docker compose up -d --wait >/dev/null 2>&1; then
  elapsed=$(( $(date +%s) - started ))
  pass "stack healthy ${elapsed} s after the volume was deleted"
else
  elapsed=$(( $(date +%s) - started ))
  fail "stack did not come up healthy (${elapsed} s) — docker compose ps; docker compose logs db"
  docker compose ps
  exit 1
fi
restore_line=$(docker compose logs db 2>/dev/null | grep -o 'restore + ANALYZE finished in [0-9]* s' | tail -1 || true)
[ -n "$restore_line" ] && pass "$restore_line (inside the db container)" || fail "no restore line in the db log — did the dump restore?"

# ── end to end, through the published port ────────────────────────────────────────────────────────
# ⛔ The origin may sit behind an HTTP-Basic gate (one shared lab credential in front of BOTH the app
# and the API). Without a credential every check below 401s and the drill reports a broken deployment
# when nothing is broken — which is what happened the first time this ran against a gated origin.
# Set `BACATLAS_DRILL_AUTH=user:password` in `.env` (gitignored) and the drill authenticates. Leave it
# unset and nothing changes, so an ungated deployment behaves exactly as before.
auth=()
[ -n "${BACATLAS_DRILL_AUTH:-}" ] && auth=(-u "$BACATLAS_DRILL_AUTH")
get() { curl -fsS --max-time 20 "${auth[@]}" "$1"; }
status() { curl -s -o /dev/null -w '%{http_code}' --max-time 20 "${auth[@]}" "$@"; }
# Deliberately WITHOUT the credential: used to prove the gate is actually closed.
status_anon() { curl -s -o /dev/null -w '%{http_code}' --max-time 20 "$@"; }
json() { docker compose exec -T api python -c "import json,sys; d=json.load(sys.stdin); print($1)"; }

# The gate, if there is one. A rebuilt `web` that silently lost its auth config would otherwise pass
# every check below and sit open on the internet — a failure no other check here can see.
if [ -n "${BACATLAS_DRILL_AUTH:-}" ]; then
  [ "$(status_anon "$url")" = "401" ] && pass "gate: the page refuses an anonymous request (401)" \
    || fail "gate: the page answered $(status_anon "$url") to NO credential — the gate is OPEN"
  [ "$(status_anon "$api/health")" = "401" ] && pass "gate: the API refuses an anonymous request (401)" \
    || fail "gate: the API answered $(status_anon "$api/health") to NO credential — the API is OPEN"
  [ "$(status_anon "${url}healthz")" = "200" ] && pass "gate: /healthz stays open, for the container healthcheck" \
    || fail "gate: /healthz answered $(status_anon "${url}healthz") — the web healthcheck will fail"
fi

health=$(get "$api/health" || true)
if echo "$health" | json 'd.get("schema_present") is True' 2>/dev/null | grep -q True; then
  pass "health: database reachable, schema present"
else
  fail "health: $health"
fi

species=$(get "$api/species" | json '" ".join(s["key"] for s in d["species"] if s["published"])' 2>/dev/null || true)
[ -n "$species" ] && pass "published species: $species" || fail "no published species in the restored database"

[ "$(status "$url")" = "200" ] && pass "the page itself: 200" || fail "the page itself: $(status "$url")"

for key in $species; do
  catalogue=$(get "$api/species/$key" || true)
  landing=$(echo "$catalogue" | json 'd["landing_locus"]' 2>/dev/null || true)
  loci=$(echo "$catalogue" | json 'd["pangenome"]["locus_count"]' 2>/dev/null || true)
  [ -n "$landing" ] && pass "$key: catalogue of $loci loci, landing locus $landing" || { fail "$key: no catalogue"; continue; }

  [ "$(status "$api/species/$key/loci/$landing")" = "200" ] && pass "$key: landing locus 200" \
    || fail "$key: landing locus $(status "$api/species/$key/loci/$landing")"

  # ⛔ The catalogue-sprite checks were REMOVED 2026-10-22. They read
  # `map_projections[].scatter_sprite` -> `scatter.png` -> its 304, but `locus_map_scatter_sprite` is
  # created by migration 858439e8fd7c and DROPPED by d7c2b6e4a915. No sprite table exists at head, so
  # the checks could not pass on any host — which made this drill, the S3 exit criterion, incapable of
  # ever returning 0. A test that cannot pass reports nothing. Restore them only with the table.

  genome=$(get "$api/species/$key/genomes?limit=1" | json 'd["genomes"][0]["sample_id"]' 2>/dev/null || true)
  [ -n "$genome" ] && pass "$key: genome list answers ($genome first)" || fail "$key: genome list"

  if [ -n "${BACATLAS_HOST_GFF_DIR:-}" ] && [ -n "$genome" ]; then
    code=$(status "$api/species/$key/genomes/$genome/loci/$landing/sequence")
    [ "$code" = "200" ] && pass "$key: sequence read from the GFF tree" \
      || fail "$key: sequence endpoint $code — is the GFF tree readable by uid 10001?"
  fi
done

echo
if [ "$failures" -eq 0 ]; then
  echo "DRILL PASSED — restored and serving in ${elapsed} s."
else
  echo "DRILL FAILED — $failures check(s) failed. The stack is left running for inspection."
  exit 1
fi
