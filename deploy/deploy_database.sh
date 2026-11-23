#!/usr/bin/env bash
# Deploy a database to this host: RESTORE, then SMOKE TEST. One command, two phases.
#
#   deploy/deploy_database.sh --yes-delete-the-database --expect-published FILE
#
#   FILE        the pangenome rows the dump must hold, and nothing else: catalogue key, run id, model
#               label, published. Written AT DEPLOY TIME from what is to be served — nuna's
#               `python -m nuna.config --write-published-json FILE` — and carried beside the dump. This
#               repo never reads nuna; it reads the file. (launching_bacatlas_org.md §7b, in nuna.)
#   RESTORE     delete the database volume and reload it from the snapshot, timed.
#   SMOKE TEST  check the stack end to end THROUGH THE PUBLISHED PORT — the path a reader's browser
#               takes. "Smoke test" in the usual sense: turn it on and see whether smoke comes out.
#               It proves the deployment basically works. It does NOT prove the data is correct.
#
# One PASS/FAIL line per check, then a summary. Exit status 0 means every check passed.
#
# ⚠ Renamed from `restore_drill.sh` (2026-10-22, David). A "restore drill" is a real term — rehearsing
# a restore you hope never to need — but it named one of this script's three jobs and the least used:
# day to day this is not a rehearsal, it is how new data reaches production.
#
# ⛔ Destructive by design: `docker compose down -v` is how a new snapshot is loaded (the postgres
# image restores only into an EMPTY volume), so this script IS the upgrade procedure, and running it
# against the CURRENT snapshot is also how you rehearse a recovery. It refuses without the flag
# above, and refuses before deleting anything if the snapshot is missing or the wrong vintage — a
# deploy that deletes the live database and then finds no snapshot has made an outage out of a test.
#
# ⛔ TWO SITES SHARE THE HOST — main (/srv/bacatlas/app, bacatlas.org) and beta (/srv/bacatlas/beta,
# bacatlas.org/beta_testing/), each a checkout with its own `.env` and its own database. This deploys the
# site its own `.env` names (BACATLAS_SITE), into the compose project it names (COMPOSE_PROJECT_NAME), and
# refuses — before anything is deleted — a `.env` that could reach the other site's volume or images, and a
# rows FILE written for the other site.
#
# Needs: Docker Engine 25+ with Compose v2, curl, and `.env` beside compose.yaml (see .env.example).
# ⚠ Nothing else on the host — no Python: JSON is read with the api container's own interpreter.
set -Eeuo pipefail

usage() {
  echo "usage: deploy/deploy_database.sh --yes-delete-the-database --expect-published FILE"
  echo "  deletes the database volume, restores it from the dump, and checks the stack end to end."
  echo "  FILE: the pangenome rows the dump must hold, written at deploy time by nuna's"
  echo "        \`python -m nuna.config --write-published-json FILE\` and carried beside the dump."
  exit 2
}
confirmed=""
expect=""
while [ $# -gt 0 ]; do
  case "$1" in
    --yes-delete-the-database) confirmed=1; shift ;;
    --expect-published) [ $# -ge 2 ] || usage; expect="$2"; shift 2 ;;
    --expect-published=*) expect="${1#--expect-published=}"; shift ;;
    *) usage ;;
  esac
done
[ -n "$confirmed" ] && [ -n "$expect" ] || usage
# Resolved BEFORE the `cd` below, so a relative path means what the operator typed.
[ -f "$expect" ] || { echo "no expected-rows file at $expect — refusing, and nothing has been deleted"; exit 2; }
expect="$(cd "$(dirname "$expect")" && pwd)/$(basename "$expect")"

cd "$(dirname "$0")/.."

failures=0
checks=0
# ⚠ `set -e` exits on an unguarded failure and prints NOTHING; after `down -v` that is an outage with
# no explanation. Say where it stopped.
# ⚠ Main shell only: `set -E` carries the trap into command substitutions, where a failing curl is
# often CONSUMED (`$(status …)` inside a `[ … ] && pass || fail`) and the deploy carries on.
trap 'rc=$?; [ "$BASH_SUBSHELL" -eq 0 ] && echo "  FAIL  unexpected error at line $LINENO (exit $rc) — the deploy stopped here; docker compose ps; docker compose logs db" >&2' ERR
pass() { printf '  PASS  %s\n' "$1"; checks=$((checks + 1)); }
fail() { printf '  FAIL  %s\n' "$1"; failures=$((failures + 1)); checks=$((checks + 1)); }

# ── preflight: nothing is deleted until all of this holds ─────────────────────────────────────────
[ -f .env ] || { echo "no .env beside compose.yaml — copy .env.example and fill it in"; exit 2; }
# ⛔ WHICH SITE, AND WHICH VOLUME — from `.env` ALONE, before docker is called at all. `down -v` below deletes
# the volume of whatever compose project COMPOSE_PROJECT_NAME names, and Compose lets the invoking SHELL
# override `.env`: a name exported for the other site would aim this deploy at that site's database. So the
# three are cleared first, `.env` must set them, and the pairs that belong to the live site are refused here.
unset COMPOSE_PROJECT_NAME BACATLAS_SITE BACATLAS_IMAGE_TAG
set -a; . ./.env; set +a
site="${BACATLAS_SITE:-}"
project="${COMPOSE_PROJECT_NAME:-}"
image_tag="${BACATLAS_IMAGE_TAG:-local}"
refuse_env() { echo "$1 — refusing, and nothing has been deleted"; exit 2; }
if [ -z "$project" ]; then
  refuse_env "COMPOSE_PROJECT_NAME is not set in .env (main: bacatlas, beta: bacatlas-beta), and it names the volume this deletes"
fi
case "$site" in
  main|beta) ;;
  *) refuse_env "BACATLAS_SITE is '$site' in .env, not main or beta" ;;
esac
if [ "$site" = beta ] && [ "$project" = bacatlas ]; then
  refuse_env "BACATLAS_SITE=beta with COMPOSE_PROJECT_NAME=bacatlas, the LIVE site's project, whose database down -v would delete"
fi
if [ "$site" = beta ] && [ "$image_tag" = local ]; then
  refuse_env "BACATLAS_SITE=beta with BACATLAS_IMAGE_TAG=local (or unset), the LIVE site's images (set BACATLAS_IMAGE_TAG=beta)"
fi
echo "  site    $site — compose project $project, images bacatlas-{api,web}:$image_tag; down -v deletes ${project}_pgdata"
dump_dir="${BACATLAS_HOST_DUMP_DIR:-./deploy/dump}"
dump="$dump_dir/bacatlas.dump"
[ -f "$dump" ] || { echo "no dump at $dump — refusing to delete a database it cannot replace"; exit 2; }
command -v curl >/dev/null || { echo "curl is required"; exit 2; }
docker compose version >/dev/null 2>&1 || { echo "Docker Compose v2 is required"; exit 2; }
# ⚠ `up --wait` with no bound waits FOREVER on a service whose health stays `starting` — a crash-looping
# api resets it on every restart — so bound it where this Compose can. Decided here, before anything is
# deleted: an unknown flag passed AFTER `down -v` would turn a deploy into an outage.
# ⚠ Read into a variable first: `docker compose up --help | grep -q` under pipefail reads as "absent"
# whenever grep exits early and the writer dies of SIGPIPE.
wait_bound=()
up_help=$(docker compose up --help 2>/dev/null || true)
grep -q -- '--wait-timeout' <<<"$up_help" && wait_bound=(--wait-timeout "${BACATLAS_UP_WAIT_SECONDS:-900}")
engine=$(docker version --format '{{.Server.Version}}' 2>/dev/null || echo 0)
[ "${engine%%.*}" -ge 25 ] 2>/dev/null || { echo "Docker Engine 25+ is required (found $engine)"; exit 2; }

# ⛔ THE DB CONTAINER MUST BE ABLE TO READ THE DUMP, and this too is checked BEFORE anything is deleted.
# The restore (deploy/postgres/initdb) runs pg_restore as the image's `postgres` user, but every check below
# reads the dump in a `docker run` that runs as ROOT — so a dump the postgres user cannot open passes them all,
# `down -v` deletes the live database, and the restore then fails on "Permission denied": an outage.
# Measured 2026-11-23: a dump scp'd from CSD3 arrived `-rw-r----- ubuntu` (RDS writes 660 + an ACL). Docker
# Desktop on the Mac maps file modes away, so the laptop drill of the SAME file passed; only the server sees it.
if ! docker run --rm --user postgres -v "$(cd "$dump_dir" && pwd)":/r:ro postgres:16 \
       test -r "/r/$(basename "$dump")"; then
  echo "⛔ THE DB CONTAINER CANNOT READ THE DUMP — refusing, and nothing has been deleted."
  echo "   $(ls -l "$dump")"
  echo "   The restore runs as the postgres image's \`postgres\` user, not as root. Make it readable:"
  echo "     chmod 644 $dump   (and its directory searchable: chmod o+x on each parent)"
  exit 2
fi

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
#
# ⛔ **THE PATTERN MATCHES THE ID, NOT ITS QUOTES — and the quotes are what made this gate wrong.**
# It used to grep `'[0-9a-f]{12}'`, single-quoted, because alembic's template writes single quotes
# and all 24 migrations then on disk did. `b1f39c8b40e9` was hand-written with DOUBLE quotes, so the
# pattern missed it, `comm -23` returned its PARENT `c99e925f7341` as the head, and the gate refused
# every dump taken after it — telling the operator to re-take a dump that would mismatch identically
# forever. One migration later it fails differently and no less fatally: with two unmatched ids the
# head is not unique and the gate exits on "no single head".
# ⚠ The fix belongs HERE, not in the migration: `backend/pyproject.toml` exempts
# `migrations/versions/*` from every ruff rule precisely so a migration's bytes never change after it
# has run, and rewriting one to suit a grep is the tail wagging the dog.
# ⭐ `tests/test_deploy_vintage_gate.py` now runs this exact pipeline and compares its answer with
# alembic's own `ScriptDirectory`, so a quoting style a future migration happens to use cannot move
# the answer again without a test going red.
versions=backend/bacatlas_backend/migrations/versions
dump_rev=$(docker run --rm -v "$(cd "$dump_dir" && pwd)":/r:ro postgres:16 \
             pg_restore -f - --data-only --table=alembic_version "/r/$(basename "$dump")" 2>/dev/null \
           | grep -Eoi '\b[0-9a-f]{12}\b' | tail -1 || true)
code_rev=$(comm -23 \
  <(grep -h '^revision' "$versions"/*.py | grep -Eo '\b[0-9a-f]{12}\b' | sort -u) \
  <(grep -h '^down_revision' "$versions"/*.py | grep -Eo '\b[0-9a-f]{12}\b' | sort -u))
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
  echo "     pg_dump -Fc -d bacatlas_ingest -f $dump"
  exit 2
fi

# ⛔ THE DUMP MUST SERVE WHAT IS MEANT TO BE SERVED, and this too is checked BEFORE anything is deleted.
# Nothing above can tell one well-formed dump from another: a dump of the wrong build — the laptop's
# drill database, a catalogue loaded for testing, last month's run — restores, passes the smoke test
# and serves the wrong model under the right name. What is meant to be served is `[published]` in nuna's
# config/config.toml (private), so nuna writes it to FILE at deploy time and this compares the DUMP with
# it: every pangenome row, published or not, by catalogue key, run id, model label and published flag —
# never a path. Equal in both directions; a row missing, a row extra, or a row that differs refuses.
# Nothing new on the host: the dump is read with the postgres image (as above), and the rows are compared
# by the api IMAGE's own Python — not `docker compose exec`, which needs a running stack that a first
# deploy or a broken one does not have. ⭐ `tests/test_deploy_published_gate.py` runs this program, lifted
# verbatim, on a real `pg_dump` of a real pangenome table.
api_image="bacatlas-api:$image_tag"   # compose.yaml tags it so (BACATLAS_IMAGE_TAG); `docker compose build` (§7e) makes it
docker image inspect "$api_image" >/dev/null 2>&1 \
  || { echo "no $api_image image on this host — run \`docker compose build\` first; nothing has been deleted"; exit 2; }
pangenome_copy=$(docker run --rm -v "$(cd "$dump_dir" && pwd)":/r:ro postgres:16 \
                   pg_restore -f - --data-only --table=pangenome "/r/$(basename "$dump")" 2>/dev/null) \
  || { echo "could not read the pangenome table from $dump — nothing has been deleted"; exit 2; }
IFS= read -r -d '' published_check_py <<'PY' || true
import json
import sys

FIELDS = ("catalogue_key", "run_id", "model_label", "is_published")
#: Each site's table in nuna's config/config.toml — named in the refusals, as nuna names it.
TABLES = {"main": "[published]", "beta": "[beta]"}
ESCAPES = {"b": "\b", "f": "\f", "n": "\n", "r": "\r", "t": "\t", "v": "\v"}


def refuse(message):
    print(f"  ✗ {message}")
    sys.exit(1)


def copy_value(text):
    """One field of COPY's text format: \\N is NULL, and a backslash escapes the character after it."""
    if text == "\\N":
        return None
    out, chars = [], iter(text)
    for char in chars:
        if char == "\\":
            char = next(chars, "")
            char = ESCAPES.get(char, char)
        out.append(char)
    return "".join(out)


def held_rows(lines):
    """The pangenome table's COPY block, as dicts keyed by the header's own column names."""
    columns = None
    rows = []
    for line in lines:
        if columns is None:
            if line.startswith("COPY ") and line.endswith(") FROM stdin;"):
                table, _, names = line[len("COPY "):-len(") FROM stdin;")].partition(" (")
                if table.rsplit(".", 1)[-1].strip('"') == "pangenome":
                    columns = [name.strip().strip('"') for name in names.split(",")]
            continue
        if line == "\\.":
            return columns, rows
        values = [copy_value(value) for value in line.split("\t")]
        if len(values) != len(columns):
            refuse(f"a pangenome row has {len(values)} fields; the COPY header names {len(columns)}")
        rows.append(dict(zip(columns, values)))
    if columns is None:
        refuse("the dump holds no data for the pangenome table")
    refuse("the pangenome table's data has no end marker — is the dump truncated?")


try:
    expected = json.loads(sys.argv[1], parse_constant=lambda name: refuse(f"the expected rows hold {name}"))
except ValueError as error:
    refuse(f"the expected-rows file is not JSON: {error}")
if not isinstance(expected, dict) or expected.get("schema") != 1 or not isinstance(expected.get("rows"), list):
    refuse("the expected-rows file is not schema 1 with a list of rows — write it again with nuna's --write-published-json")
# ⛔ The file names the site it was written for, and it must be THIS checkout's (BACATLAS_SITE, argv[2]): each
# site serves its own table, so the other site's rows would refuse this dump for the wrong reason — or, if the
# two tables happen to agree, pass it. Exit 3, not 1: the file is wrong, not the dump.
site = sys.argv[2]
if expected.get("site") != site:
    named = f"is for site {expected['site']!r}" if "site" in expected else "names no site"
    print(f"  ✗ the expected-rows file {named}; this checkout's .env says BACATLAS_SITE={site}")
    sys.exit(3)
for row in expected["rows"]:
    if not isinstance(row, dict) or sorted(row) != sorted(FIELDS) or not isinstance(row["is_published"], bool):
        refuse(f"an expected row is not {list(FIELDS)} with a true/false is_published: {row!r}")
want = {row["catalogue_key"]: row for row in expected["rows"]}
if len(want) != len(expected["rows"]) or not want:
    refuse(f"the expected rows are empty or name a catalogue twice: {[row['catalogue_key'] for row in expected['rows']]}")

# ⚠ Split on "\n" only: COPY escapes a newline inside a value, but not the other characters that
# `str.splitlines` also breaks on (\x1c, \x85, \u2028 and others).
lines = sys.stdin.read().split("\n")
if lines and lines[-1] == "":
    lines.pop()  # the newline that ENDS the last line opens no row
columns, rows = held_rows(lines)
absent = [field for field in FIELDS if field not in columns]
if absent:
    refuse(f"the dump's pangenome table has no column {absent}")
held = []
for row in rows:
    flag = {"t": True, "f": False}.get(row["is_published"])
    if flag is None:
        refuse(f"{row['catalogue_key']}: is_published is {row['is_published']!r}, not t or f")
    held.append({"catalogue_key": row["catalogue_key"], "run_id": row["run_id"],
                 "model_label": row["model_label"], "is_published": flag})

# ⛔ The same comparison, and the same words, as nuna's `config.pangenome_row_differences` — the build job's
# check of its database before the dump ships. An empty dump is a difference: no rows read is not a match.
problems = []
if not held:
    problems.append("no pangenome rows were read — nothing to compare")
else:
    have = {}
    for row in held:
        if row["catalogue_key"] in have:
            problems.append(f"{row['catalogue_key']}: held twice")
        have[row["catalogue_key"]] = row
    for key in sorted(set(want) | set(have)):
        if key not in have:
            problems.append(f"{key}: expected, and not held")
        elif key not in want:
            problems.append(f"{key}: held (is_published={have[key]['is_published']}), and not in {TABLES[site]}")
        else:
            for field in FIELDS[1:]:
                if have[key][field] != want[key][field]:
                    problems.append(f"{key}: {field} is {have[key][field]!r}, expected {want[key][field]!r}")
for problem in problems:
    print(f"  ✗ {problem}")
if problems:
    sys.exit(1)
print(f"  rows    the dump's {len(held)} pangenome rows are exactly the expected ({', '.join(sorted(want))}; "
      f"{site} site; config sha256 {str(expected.get('config_sha256', '?'))[:16]}…)")
PY
if published_verdict=$(printf '%s\n' "$pangenome_copy" \
      | docker run --rm -i --pull never --network none --entrypoint python "$api_image" \
          -c "$published_check_py" "$(cat "$expect")" "$site" 2>&1); then
  printf '%s\n' "$published_verdict"
  published_rows="rows ${expect##*/} ✓"   # named in the closing line: a script without this gate never prints it
else
  published_rc=$?
  printf '%s\n' "$published_verdict"
  # The program exits 1 on a dump it refuses, 3 on a rows file for another site, and only then; anything else
  # is the container failing to run it.
  if [ "$published_rc" -eq 3 ]; then
    echo "⛔ THE EXPECTED ROWS ARE NOT FOR THIS SITE — refusing, and nothing has been deleted."
    echo "     this checkout   BACATLAS_SITE=$site (compose project $project)"
    echo "     expected rows   $expect"
    echo "   Write them for this site, from the commit that serves it:"
    echo "     python -m nuna.config --write-published-json FILE --site $site"
    exit 2
  fi
  if [ "$published_rc" -ne 1 ]; then
    echo "⛔ COULD NOT RUN THE ROW CHECK in $api_image (exit $published_rc) — refusing, and nothing has been deleted."
    echo "   This says nothing about the dump: rebuild the image (docker compose build) and run this again."
    exit 2
  fi
  echo "⛔ THE DUMP DOES NOT HOLD WHAT IS TO BE SERVED — refusing, and nothing has been deleted."
  echo "     dump            $dump"
  echo "     expected rows   $expect"
  echo "   Either this is not the dump meant for this deploy, or the file was not written from the commit"
  echo "   that serves it: a deploy moves the catalogue into this site's table in nuna's config (main: [published],"
  echo "   beta: [beta]) in the same commit, and the file is written from that commit"
  echo "   (python -m nuna.config --write-published-json FILE --site $site)."
  exit 2
fi

# ⛔ THE DISK MUST HOLD THE RESTORE, and this too is checked BEFORE anything is deleted.
# `docker compose down -v` deletes the serving database first and restores second, so a restore that
# runs out of disk destroys the live database and leaves nothing — the vintage gate above cannot see
# that, and neither can the health check until it is too late. Measured on the 100-genome probe: the
# restored database is ~4x the compressed dump (`pg_dump -Fc` emits bytea as hex and omits every
# index), and a restore peaks at ~1.6x the finished database while WAL and the index builds overlap.
# So the free space on the volume's disk must cover the dump itself plus FACTOR x its size, where the
# factor defaults to 8 (= 1 + 4 x 1.6, rounded up) and is overridable for a cohort measured otherwise.
factor="${BACATLAS_RESTORE_SIZE_FACTOR:-8}"
docker_root=$(docker info --format '{{.DockerRootDir}}' 2>/dev/null || echo /var/lib/docker)
[ -d "$docker_root" ] || docker_root=/
dump_kb=$(du -k "$dump" | cut -f1)
free_kb=$(df -Pk "$docker_root" | awk 'NR==2 {print $4}')
need_kb=$(( dump_kb * (factor + 1) ))
if [ "${free_kb:-0}" -lt "$need_kb" ]; then
  echo "⛔ NOT ENOUGH DISK FOR THE RESTORE — refusing, and nothing has been deleted."
  echo "     dump            $(( dump_kb / 1024 )) MB at $dump"
  echo "     needed          $(( need_kb / 1024 )) MB  (dump x (1 + $factor); BACATLAS_RESTORE_SIZE_FACTOR=$factor)"
  echo "     free            $(( ${free_kb:-0} / 1024 )) MB on $docker_root"
  echo "   A restore that fills the disk destroys the serving database and restores nothing. Free space"
  echo "   first (old dumps under $dump_dir, docker system prune), or lower the factor only with a measured"
  echo "   database:dump ratio for THIS cohort."
  exit 2
fi
echo "  disk    $(( ${free_kb:-0} / 1024 )) MB free on $docker_root for a $(( dump_kb / 1024 )) MB dump (x$((factor + 1)) needed)"

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

# ── RESTORE ──────────────────────────────────────────────────────────────────────────────────────
echo "RESTORE — deleting the $site site's database volume (${project}_pgdata) and reloading it from the snapshot…"
docker compose down -v --remove-orphans >/dev/null 2>&1
started=$(date +%s)
# ⛔ WAIT FOR THE RESTORE ITSELF, NOT FOR THE HEALTHCHECK. While the first-start restore runs, the
# image's temporary server listens on the socket only, so the TCP healthcheck fails throughout, and
# compose counts those failures once `start_period` (10m) is over: five more and the db is marked
# unhealthy, `up --wait` gives up, and the api is never started — though the restore carries on and
# finishes. A 100-genome dump restores in under a minute; the 71,435-genome kp-71k dump (5.6 GiB, a
# 16 GB database) does not, and how long it takes depends on the host. So start the db alone, watch
# its log for the init script's own verdict however long that takes, and only then wait for health.
# A failed restore is refused the moment its line appears, never after a timeout.
docker compose up -d db >/dev/null 2>&1 || { fail "the db container did not start — docker compose logs db"; exit 1; }
db_id=$(docker compose ps -q db 2>/dev/null || true)
# ⛔ The verdict line can also NEVER come, and a loop that waits only for it would print "still
# restoring" forever after the volume was deleted. The init script runs only on the first start of an
# EMPTY volume, so a db container that restarted mid-restore (killed, OOM, a dockerd restart) skips it,
# says "Skipping initialization", and comes up on a PARTIAL database. That, a stopped or restarting
# container, and a ceiling (BACATLAS_RESTORE_TIMEOUT_MINUTES, default 360) each end the wait with a FAIL.
deadline=$(( started + ${BACATLAS_RESTORE_TIMEOUT_MINUTES:-360} * 60 ))
noted=$started
while :; do
  db_log=$(docker compose logs db 2>/dev/null || true)
  if grep -q 'RESTORE FAILED' <<<"$db_log"; then
    fail "the restore FAILED after $(( $(date +%s) - started )) s — docker compose logs db"
    exit 1
  fi
  grep -qE 'restore \+ ANALYZE finished in [0-9]+ s|the database starts EMPTY' <<<"$db_log" && break
  db_state=$(docker inspect -f '{{.State.Status}} {{.RestartCount}}' "$db_id" 2>/dev/null || echo "missing ?")
  if grep -q 'Skipping initialization' <<<"$db_log" || [ "$db_state" != "running 0" ]; then
    fail "the db container stopped or restarted mid-restore (state: $db_state) — its database is PARTIAL; docker compose down -v, then re-run"
    exit 1
  fi
  if [ "$(date +%s)" -ge "$deadline" ]; then
    fail "no restore verdict after $(( ($(date +%s) - started) / 60 )) min — it may still be running; docker compose logs db"
    exit 1
  fi
  if [ $(( $(date +%s) - noted )) -ge 300 ]; then
    noted=$(date +%s)
    echo "  … still restoring, $(( (noted - started) / 60 )) min since the volume was deleted"
  fi
  sleep 10
done
# The image restarts the server with TCP after the init script; a db already marked unhealthy flips
# back to healthy on its next passing check (interval 10s), and `up --wait` refuses an unhealthy one.
for _ in $(seq 1 60); do
  [ "$(docker inspect -f '{{if .State.Health}}{{.State.Health.Status}}{{end}}' "$db_id" 2>/dev/null || true)" = healthy ] && break
  sleep 5
done
up_log=$(mktemp)
if docker compose up -d --wait ${wait_bound[@]+"${wait_bound[@]}"} >"$up_log" 2>&1; then
  elapsed=$(( $(date +%s) - started ))
  pass "stack healthy ${elapsed} s after the volume was deleted"
  rm -f "$up_log"
else
  elapsed=$(( $(date +%s) - started ))
  fail "stack did not come up healthy (${elapsed} s) — docker compose ps; docker compose logs db"
  tail -20 "$up_log"; rm -f "$up_log"
  docker compose ps
  exit 1
fi
restore_line=$(docker compose logs db 2>/dev/null | grep -o 'restore + ANALYZE finished in [0-9]* s' | tail -1 || true)
[ -n "$restore_line" ] && pass "$restore_line (inside the db container)" || fail "no restore line in the db log — did the dump restore?"

restore_checks=$checks

# ── SMOKE TEST — end to end, through the published port ──────────────────────────────────────────
echo
echo "SMOKE TEST — through the published port, the path a reader's browser takes…"
# ⛔ The origin may sit behind BASIC AUTH (HTTP Basic authentication: one shared credential, in front
# of BOTH the app and the API — no accounts, no sessions, no logout). Without a credential every check
# below 401s and the smoke test reports a broken deployment when nothing is broken, which is what
# happened the first time it ran against a protected origin. Set
# `BACATLAS_SMOKE_TEST_AUTH=user:password` in `.env` (gitignored) and the smoke test authenticates.
# Leave it unset and nothing changes, so an unprotected deployment behaves exactly as before.
auth=()
[ -n "${BACATLAS_SMOKE_TEST_AUTH:-}" ] && auth=(-u "$BACATLAS_SMOKE_TEST_AUTH")
# ⚠ `${auth[@]+"${auth[@]}"}`, not `"${auth[@]}"`: bash before 4.4 — macOS's /bin/bash is 3.2 —
# treats an EMPTY array as unbound under `set -u`, so with no credential set the first request died
# with "auth[@]: unbound variable" (the laptop drill of the kp-71k dump, 2026-11-14).
get() { curl -fsS --max-time 20 ${auth[@]+"${auth[@]}"} "$1"; }
status() { curl -s -o /dev/null -w '%{http_code}' --max-time 20 ${auth[@]+"${auth[@]}"} "$@"; }
# Deliberately WITHOUT the credential: used to prove basic auth is actually in force.
status_anon() { curl -s -o /dev/null -w '%{http_code}' --max-time 20 "$@"; }
json() { docker compose exec -T api python -c "import json,sys; d=json.load(sys.stdin); print($1)"; }

# Basic auth, if this origin uses it. A rebuilt `web` that silently lost its auth configuration would
# otherwise pass every check below while sitting open on the internet — a failure no other check here
# can see, because every other one arrives authenticated.
if [ -n "${BACATLAS_SMOKE_TEST_AUTH:-}" ]; then
  [ "$(status_anon "$url")" = "401" ] && pass "basic auth: the page refuses an anonymous request (401)" \
    || fail "basic auth: the page answered $(status_anon "$url") to NO credential — it is OPEN"
  [ "$(status_anon "$api/health")" = "401" ] && pass "basic auth: the API refuses an anonymous request (401)" \
    || fail "basic auth: the API answered $(status_anon "$api/health") to NO credential — the API is OPEN"
  [ "$(status_anon "${url}healthz")" = "200" ] && pass "basic auth: /healthz stays exempt, for the container healthcheck" \
    || fail "basic auth: /healthz answered $(status_anon "${url}healthz") — the web healthcheck will fail"
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

  # ⭐ The Function tab, once per WORKER, all at once. Its agreement ladder is computed on first use
  # and kept for the life of the process (`function_inference_service.calibration_for`): 5.75 s per
  # worker on the 71,435-genome catalogue on the laptop, and slower on a 2-vCPU host. A sync gunicorn
  # worker serves one request at a time, so as many concurrent requests as there are workers land on
  # every worker — which both proves the route answers and warms each worker before a reader arrives.
  # ⛔ Its own, longer time limit: the 20 s the other checks use would fail a cold ladder that is
  # merely slow. The real bounds are gunicorn's `--timeout` (backend/Dockerfile) and nginx's
  # `proxy_read_timeout` — and on the gated instance nginx serves the BIND-MOUNTED copy of its
  # template, which a `git pull` does not update (launching_bacatlas_org.md §7a): until that copy is
  # re-rendered, nginx's default 60 s is what bounds this.
  workers="${BACATLAS_API_WORKERS:-4}"
  fn_started=$(date +%s)
  fn_codes=$(for _ in $(seq 1 "$workers"); do
               curl -s -o /dev/null -w '%{http_code} ' --max-time 170 ${auth[@]+"${auth[@]}"} \
                 "$api/species/$key/loci/$landing/function" &
             done; wait)
  fn_seconds=$(( $(date +%s) - fn_started ))
  fn_bad=0
  fn_count=0
  for code in $fn_codes; do
    fn_count=$((fn_count + 1))
    [ "$code" = "200" ] || fn_bad=$((fn_bad + 1))
  done
  if [ "$fn_bad" -eq 0 ] && [ "$fn_count" -eq "$workers" ]; then
    pass "$key: function tab 200 on $workers concurrent requests, ${fn_seconds} s — every worker's ladder is warm"
  else
    fail "$key: function tab answered [${fn_codes% }] to $workers concurrent requests (${fn_seconds} s)"
  fi

  genome=$(get "$api/species/$key/genomes?limit=1" | json 'd["genomes"][0]["sample_id"]' 2>/dev/null || true)
  [ -n "$genome" ] && pass "$key: genome list answers ($genome first)" || fail "$key: genome list"

  # ⛔ **The sequence checks pick a genome from the SOURCE GROUP each station serves.** The 71,435-genome
  # catalogue holds two populations — BioSample-keyed AllTheBacteria assemblies and RefSeq/GenBank stems
  # that only NCBI can serve — and the first genome in the list may be either. A check labelled
  # "AllTheBacteria" that took it blind could pass by reading NCBI, or 503 on a stem with NCBI off: a
  # gate running a configuration its label does not name. `source_group` is the API's own derivation.
  atb_genome=$(get "$api/species/$key/genomes?source_group=atb&limit=1" | json 'd["genomes"][0]["sample_id"]' 2>/dev/null || true)
  ncbi_genome=""
  for group in ncbi_complete ncbi_draft; do
    [ -n "$ncbi_genome" ] && break
    ncbi_genome=$(get "$api/species/$key/genomes?source_group=$group&limit=1" | json 'd["genomes"][0]["sample_id"]' 2>/dev/null || true)
  done

  # One locus, one genome, through the published port. ⚠ A 200 with NO genes is the answer "this genome
  # has no gene here" — and the provider is never asked for such a genome, so the FETCH was not
  # exercised. That is said, not passed off as a read: silence is not success.
  # ⛔ **The sequence checks ask at the catalogue's MOST PREVALENT locus, not its landing locus.** On
  # kp-71k the landing locus (12517) is in 249 of 71,435 genomes and the four example loci in 9-10,
  # so the first genome of each source group had no gene at any of them: both checks answered "200,
  # but the fetch was NOT exercised" — a deploy whose two remote stations were never touched, passing.
  # The most prevalent locus is in (nearly) every genome. It is read from the database because no
  # route lists a genome's loci; the request itself still goes through the published port.
  seq_locus=$(printf '%s\n' "select l.node_label from locus l join pathogen_species s on s.default_pangenome_id = l.pangenome_id where s.species_key = :'key' order by l.member_genome_count desc, l.locus_id limit 1;" \
                | docker compose exec -T db sh -c 'psql -U "$POSTGRES_USER" -d "$POSTGRES_DB" -At -v key="$1"' sh "$key" 2>/dev/null || true)
  if [ -z "$seq_locus" ]; then
    seq_locus="$landing"
    echo "  note  $key: could not read the most prevalent locus from the db — probing the landing locus, which may not exercise a fetch"
  fi

  check_sequence() {  # <source label> <sample_id>
    local source_name="$1" probe="$2" code genes
    # ⚠ The FIRST request is a real remote fetch, always cold (`down -v` drops the assembly cache), and
    # the NCBI station's own timeout is 30 s — so 60 s here, and `|| true` so a curl that gives up
    # reads as code 000 through the named FAIL below rather than ending the smoke test under `set -e`.
    code=$(status --max-time 60 "$api/species/$key/genomes/$probe/loci/$seq_locus/sequence" || true)
    if [ "$code" != "200" ]; then
      fail "$key: sequence endpoint $code for $probe from $source_name — GFF readable by uid 10001, or $source_name reachable?"
      return
    fi
    genes=$(get "$api/species/$key/genomes/$probe/loci/$seq_locus/sequence" | json 'len(d["genes"])' 2>/dev/null || true)
    if [ "${genes:-0}" -gt 0 ] 2>/dev/null; then
      pass "$key: sequence read from $source_name ($probe, $genes gene(s) at $seq_locus)"
    else
      pass "$key: sequence endpoint 200 for $probe, but it has no gene at $seq_locus — the $source_name fetch was NOT exercised"
    fi
  }

  # ⛔ **Gated on EITHER source, and it says which.** This was gated on the GFF tree alone, so the
  # moment the tree was cleared to serve sequence from AllTheBacteria the check vanished — the drill
  # went on returning 0 while the one test that proves sequence works stopped running, precisely
  # when sequence had moved to a remote source and needed testing most. A check that silently skips
  # reports nothing, which is the same defect the removed sprite checks had in the other direction.
  if [ -n "${BACATLAS_HOST_GFF_DIR:-}" ] || [ -n "${BACATLAS_ATB_STREAMING:-}" ]; then
    if [ -n "${BACATLAS_HOST_GFF_DIR:-}" ]; then source_name="the GFF tree"; else source_name="AllTheBacteria"; fi
    if [ -n "$atb_genome" ]; then
      check_sequence "$source_name" "$atb_genome"
    else
      pass "$key: no BioSample-keyed genome in the catalogue — $source_name not exercised"
    fi
  else
    # ⚠ Not silence: a stack serving no sequence at all is a legitimate configuration and must say so.
    pass "$key: sequence not configured (no GFF tree, no ATB streaming) — AllTheBacteria not exercised"
  fi

  # ⭐ The fourth station: a RefSeq/GenBank-keyed genome, whose bases only NCBI can supply.
  if [ -n "${BACATLAS_NCBI_STREAMING:-}" ]; then
    if [ -n "$ncbi_genome" ]; then
      check_sequence "NCBI" "$ncbi_genome"
    else
      pass "$key: no NCBI-keyed genome in the catalogue — NCBI streaming not exercised"
    fi
  elif [ -n "$ncbi_genome" ]; then
    # ⚠ Named: this catalogue HAS genomes only NCBI can serve, and the station is off. Their Sequence
    # tab answers a named 503 — a legitimate configuration, but one the operator should see stated.
    pass "$key: NCBI streaming off, and $ncbi_genome is NCBI-keyed — its sequence is not served"
  else
    pass "$key: NCBI streaming off and no NCBI-keyed genome in the catalogue — nothing to exercise"
  fi
done

echo
smoke_checks=$(( checks - restore_checks ))
if [ "$failures" -eq 0 ]; then
  echo "DATABASE DEPLOYED — ${published_rows}, restore ${restore_checks}/${restore_checks}, smoke test ${smoke_checks}/${smoke_checks}, serving in ${elapsed} s."
else
  echo "DEPLOY FAILED — $failures of $checks check(s) failed. The stack is left running for inspection."
  exit 1
fi
