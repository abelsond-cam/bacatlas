"""`deploy_database.sh --expect-published FILE` refuses a dump that does not hold exactly the expected rows.

⛔ **Why the gate exists.** Every other check in the deploy is blind to WHICH build it is restoring: a dump of
the laptop's drill database, a catalogue loaded for testing, or last month's run restores cleanly, passes the
vintage gate and the smoke test, and serves the wrong model under the right name. What is meant to be served
is written down in one place — `[published]` in nuna's config/config.toml (private) — so nuna writes the
expected rows to a file at deploy time, and the deploy compares the DUMP with it before `down -v`.

Three layers, each run for real:

* the comparison program, **lifted verbatim** out of the script (a copy is how two things drift — the reason
  `test_deploy_vintage_gate.py` lifts its pipeline too), fed by `pg_restore` of a real `pg_dump -Fc`;
* the same program on hand-written COPY text, for the malformed inputs a real dump cannot be made to hold;
* the SCRIPT itself, run in a scratch tree with `docker` stubbed: a refused dump stops before `down -v`, an
  accepted one reaches it — and no real docker, database or `.env` is ever touched.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from sqlalchemy import delete
from sqlalchemy.orm import Session

from bacatlas_backend.models.enumerations import RosterLineage
from bacatlas_backend.models.genome_collection import GenomeCollection
from bacatlas_backend.models.pangenome import Pangenome
from bacatlas_backend.models.pathogen_species import PathogenSpecies

BACKEND_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = BACKEND_ROOT.parent
GATE = REPO_ROOT / "deploy" / "deploy_database.sh"
RESTORE_FLAGS = ["-f", "-", "--data-only", "--table=pangenome"]

pytestmark = pytest.mark.skipif(not GATE.exists(), reason=f"{GATE} is not in this checkout")

needs_pg_tools = pytest.mark.skipif(
    not (shutil.which("pg_dump") and shutil.which("pg_restore")), reason="pg_dump / pg_restore not on PATH"
)


def _gate_program() -> str:
    """The script's comparison program, verbatim. Raises rather than testing nothing if it was rewritten."""
    # ⚠ split("\n"), not splitlines(): the latter also breaks on characters bash does not ( , \x1c, …),
    # and a comment holding one turned the lifted program into a SyntaxError the script never had.
    lines = GATE.read_text(encoding="utf-8").split("\n")
    opener = "IFS= read -r -d '' published_check_py <<'PY' || true"
    assert lines.count(opener) == 1, f"re-anchor this test: no single `{opener}` in {GATE}"
    start = lines.index(opener) + 1
    end = lines.index("PY", start)
    return "\n".join(lines[start:end]) + "\n"


def test_the_script_reads_the_dump_with_the_flags_this_test_uses():
    text = GATE.read_text(encoding="utf-8")
    assert "pg_restore " + " ".join(RESTORE_FLAGS) + ' "/r/$(basename "$dump")"' in text


def _expected(*rows: tuple[str, str, str, bool]) -> dict:
    return {
        "schema": 1,
        "config_sha256": "0123456789abcdef" * 4,
        "rows": [{"catalogue_key": k, "run_id": r, "model_label": m, "is_published": p} for k, r, m, p in rows],
    }


def _check(copy_text: str, expected: dict | str) -> subprocess.CompletedProcess:
    argument = expected if isinstance(expected, str) else json.dumps(expected)
    return subprocess.run(
        [sys.executable, "-c", _gate_program(), argument],
        input=copy_text, capture_output=True, text=True, check=False,
    )  # fmt: skip


# ── a real dump of a real pangenome table ────────────────────────────────────────────────────────────────
#: A label Postgres must ESCAPE in COPY text — a tab, a newline, a backslash — so the program's unescaping
#: is checked against what pg_restore actually writes, not against this file's idea of it.
AWKWARD = "lab\tel\nwith \\ in it"
ROWS = [
    ("ecoli-probe100", "ecoli_run_res0.1_seed0", "ecoli_label", True),
    ("kp-71k", "kp71k_run_res0.6_seed0", AWKWARD, True),
]


@pytest.fixture(scope="module")
def owners(engine):
    with Session(engine) as session:
        species = PathogenSpecies(species_key="gate", scientific_name="Klebsiella pneumoniae")
        session.add(species)
        session.flush()
        collection = GenomeCollection(
            pathogen_species_id=species.pathogen_species_id,
            collection_key="gate100",
            collection_label="Gate probe",
            roster_lineage=RosterLineage.PROBE_BIOSAMPLE,
            genome_count=100,
        )
        session.add(collection)
        session.flush()
        ids = species.pathogen_species_id, collection.genome_collection_id
        session.commit()
    return ids


@pytest.fixture()
def dump_of(engine, owners, tmp_path):
    """`dump_of(rows)` → the path of a `pg_dump -Fc` whose pangenome table holds exactly `rows`."""
    url = engine.url.set(drivername="postgresql").render_as_string(hide_password=False)

    def make(rows) -> Path:
        with Session(engine) as session:
            session.execute(delete(Pangenome))
            for key, run_id, label, published in rows:
                session.add(Pangenome(
                    run_id=run_id, catalogue_key=key, model_label=label, is_published=published,
                    ingest_generation=1, pathogen_species_id=owners[0], genome_collection_id=owners[1],
                    genome_count=100, gene_count=1, locus_count=1,
                    # JSONB with escapes of its own, so a sibling column cannot shift the fields
                    run_manifest_json={"note": "a\ttab, a\nnewline, a \\ backslash"},
                ))  # fmt: skip
            session.commit()
        dump = tmp_path / f"gate-{len(list(tmp_path.glob('*.dump')))}.dump"
        subprocess.run(["pg_dump", "-Fc", "-t", "pangenome", "-f", str(dump), url], check=True)
        return dump

    return make


def _restore(dump: Path) -> str:
    return subprocess.run(["pg_restore", *RESTORE_FLAGS, str(dump)], capture_output=True, text=True,
                          check=True).stdout  # fmt: skip


@needs_pg_tools
def test_a_dump_holding_exactly_the_expected_rows_passes(dump_of):
    got = _check(_restore(dump_of(ROWS)), _expected(*ROWS))
    assert got.returncode == 0, got.stdout + got.stderr
    assert "the dump's 2 pangenome rows are exactly the expected (ecoli-probe100, kp-71k" in got.stdout


MISMATCHES = [
    pytest.param(
        [ROWS[0], ("kp-71k", "kp71k_OTHER_run", AWKWARD, True)],
        "kp-71k: run_id is 'kp71k_OTHER_run', expected 'kp71k_run_res0.6_seed0'",
        id="another-run",
    ),  # fmt: skip
    pytest.param(
        [ROWS[0], ("kp-71k", ROWS[1][1], "kp71k_hypergeo", True)],
        "kp-71k: model_label is 'kp71k_hypergeo', expected",
        id="another-label",
    ),  # fmt: skip
    pytest.param([ROWS[0], (*ROWS[1][:3], False)], "kp-71k: is_published is False, expected True", id="hidden"),
    pytest.param([ROWS[0]], "kp-71k: expected, and not held", id="missing"),
    pytest.param(
        [*ROWS, ("kp-probe100", "kp_run", "kp_label", True)],
        "kp-probe100: held (is_published=True), and not in [published]",
        id="extra-published",
    ),  # fmt: skip
    pytest.param(
        [*ROWS, ("kp-nuna5", "kp5_run", "kp5_label", False)],
        "kp-nuna5: held (is_published=False), and not in [published]",
        id="extra-hidden",
    ),  # fmt: skip
    pytest.param([], "no pangenome rows were read — nothing to compare", id="empty-table"),
]


@needs_pg_tools
@pytest.mark.parametrize("held, reason", MISMATCHES)
def test_a_dump_that_differs_in_any_way_is_refused(dump_of, held, reason):
    got = _check(_restore(dump_of(held)), _expected(*ROWS))
    assert got.returncode == 1, got.stdout + got.stderr
    assert f"  ✗ {reason}" in got.stdout


# ── inputs a real dump cannot be made to hold ────────────────────────────────────────────────────────────
HEADER = "COPY public.pangenome (pangenome_id, run_id, catalogue_key, model_label, is_published) FROM stdin;"
ONE = _expected(("a-1", "r", "m", True))

MALFORMED = [
    pytest.param("", ONE, "the dump holds no data for the pangenome table", id="no-copy-block"),
    pytest.param(
        "COPY public.genome (genome_id) FROM stdin;\n1\n\\.\n",
        ONE,
        "the dump holds no data for the pangenome table",
        id="another-table-only",
    ),  # fmt: skip
    pytest.param(f"{HEADER}\n1\tr\ta-1\tm\tt\n", ONE, "has no end marker", id="truncated"),
    pytest.param(
        f"{HEADER}\n1\tr\ta-1\tm\n\\.\n", ONE, "a pangenome row has 4 fields; the COPY header names 5", id="short-row"
    ),  # fmt: skip
    pytest.param(
        f"{HEADER}\n1\tr\ta-1\tm\t\\N\n\\.\n", ONE, "a-1: is_published is None, not t or f", id="null-flag"
    ),  # fmt: skip
    pytest.param(
        "COPY public.pangenome (pangenome_id, run_id, catalogue_key) FROM stdin;\n\\.\n",
        ONE,
        "the dump's pangenome table has no column ['model_label', 'is_published']",
        id="no-column",
    ),  # fmt: skip
    pytest.param(
        f"{HEADER}\n1\tr\ta-1\tm\tt\n\\.\n", "{not json", "the expected-rows file is not JSON", id="not-json"
    ),  # fmt: skip
    pytest.param(f"{HEADER}\n1\tr\ta-1\tm\tt\n\\.\n", {**ONE, "schema": 2}, "is not schema 1", id="schema-2"),
    pytest.param(
        f"{HEADER}\n1\tr\ta-1\tm\tt\n\\.\n",
        '{"schema": 1, "rows": [], "x": NaN}',
        "the expected rows hold NaN",
        id="nan",
    ),  # fmt: skip
    pytest.param(f"{HEADER}\n1\tr\ta-1\tm\tt\n\\.\n", _expected(), "the expected rows are empty", id="no-rows"),
    pytest.param(
        f"{HEADER}\n1\tr\ta-1\tm\tt\n\\.\n",
        _expected(("a-1", "r", "m", True), ("a-1", "r", "m", True)),
        "name a catalogue twice",
        id="a-key-twice",
    ),  # fmt: skip
    pytest.param(
        f"{HEADER}\n1\tr\ta-1\tm\tt\n\\.\n",
        {"schema": 1, "rows": [{"catalogue_key": "a-1", "run_id": "r", "model_label": "m", "is_published": "t"}]},
        "with a true/false is_published",
        id="flag-as-text",
    ),  # fmt: skip
    pytest.param(f"{HEADER}\n1\tr\ta-1\tm\tt\n2\tr\ta-1\tm\tt\n\\.\n", ONE, "a-1: held twice", id="held-twice"),
]


@pytest.mark.parametrize("copy_text, expected, reason", MALFORMED)
def test_malformed_input_is_refused_by_name(copy_text, expected, reason):
    got = _check(copy_text, expected)
    assert got.returncode == 1, got.stdout + got.stderr
    assert reason in got.stdout, got.stdout


def test_the_header_is_read_by_name_not_by_position():
    reordered = (
        "COPY public.pangenome (is_published, model_label, catalogue_key, run_id) FROM stdin;\nt\tm\ta-1\tr\n\\.\n"
    )
    got = _check(reordered, ONE)
    assert got.returncode == 0, got.stdout + got.stderr


# ── the script itself, with docker stubbed ───────────────────────────────────────────────────────────────
STUB = r"""#!/usr/bin/env bash
printf '%s\n' "$*" >> "$STUB_LOG"
case "$1 ${2:-}" in
  "compose version") echo "Docker Compose version v2.29.0"; exit 0 ;;
  "compose up") [ "${3:-}" = "--help" ] && { echo "  --wait-timeout int"; exit 0; } ;;
  "compose down") exit 97 ;;   # ⛔ the destructive step: the gate let the deploy through, and the stub stops it here
  "version --format") echo "27.3.1"; exit 0 ;;
  "info --format") echo "$STUB_DOCKER_ROOT"; exit 0 ;;
  "image inspect") [ -z "${STUB_NO_API_IMAGE:-}" ]; exit ;;
esac
if [ "$1" = run ]; then
  case "$*" in
    *--table=alembic_version*)
      printf 'COPY public.alembic_version (version_num) FROM stdin;\n%s\n\\.\n' "$STUB_CODE_HEAD"; exit 0 ;;
    *--table=pangenome*) exec pg_restore -f - --data-only --table=pangenome "$STUB_DUMP" ;;
    *bacatlas-api:local*)
      if [ -n "${STUB_API_EXIT:-}" ]; then
        echo "docker: exec: \"python\": executable file not found" >&2; exit "$STUB_API_EXIT"
      fi
      # as docker does: without -i the container gets NO stdin, whatever the pipe feeding `docker run` holds
      interactive=""
      while [ "$1" != "bacatlas-api:local" ]; do [ "$1" = -i ] && interactive=1; shift; done; shift
      if [ -n "$interactive" ]; then exec "$STUB_PYTHON" "$@"; else exec "$STUB_PYTHON" "$@" </dev/null; fi ;;
  esac
fi
echo "stub docker: unexpected call: $*" >&2
exit 98
"""


def _code_head() -> str:
    got = subprocess.run([sys.executable, "-m", "alembic", "heads"], cwd=BACKEND_ROOT, capture_output=True,
                         text=True, check=False)  # fmt: skip
    if got.returncode != 0:
        pytest.skip(f"alembic heads failed in this interpreter: {got.stderr.strip()}")
    (head,) = re.findall(r"\b[0-9a-f]{12}\b", got.stdout)
    return head


@pytest.fixture()
def scratch_tree(tmp_path):
    """A copy of the script in a tree of its own: no real `.env`, dump, compose stack or docker."""
    root = tmp_path / "app"
    (root / "deploy").mkdir(parents=True)
    shutil.copy2(GATE, root / "deploy" / "deploy_database.sh")
    (root / "backend").symlink_to(BACKEND_ROOT)  # the migrations the vintage gate reads
    (root / "dump").mkdir()
    (root / ".env").write_text(f"BACATLAS_HOST_DUMP_DIR={root / 'dump'}\n")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "docker").write_text(STUB)
    (bin_dir / "docker").chmod(0o755)
    log = tmp_path / "docker.log"
    env = {
        "PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}",
        "HOME": str(tmp_path),
        "STUB_LOG": str(log),
        "STUB_DOCKER_ROOT": str(tmp_path),
        "STUB_PYTHON": sys.executable,
    }
    return root, env, log


def _deploy(root: Path, env: dict, *args: str, cwd: Path | None = None) -> subprocess.CompletedProcess:
    return subprocess.run(["bash", str(root / "deploy" / "deploy_database.sh"), *args], cwd=cwd or root,
                          env=env, capture_output=True, text=True, check=False)  # fmt: skip


@pytest.mark.parametrize(
    "args",
    [
        pytest.param((), id="no-flags"),
        pytest.param(("--yes-delete-the-database",), id="no-expected-rows"),
        pytest.param(("--expect-published", "x.json"), id="not-confirmed"),
        pytest.param(("--yes-delete-the-database", "--expect-published"), id="flag-without-file"),
        pytest.param(("--yes-delete-the-database", "--expect-published", "x.json", "--extra"), id="unknown-flag"),
    ],
)
def test_the_script_requires_both_flags_and_touches_nothing_without_them(scratch_tree, args):
    root, env, log = scratch_tree
    got = _deploy(root, env, *args)
    assert got.returncode == 2
    assert "usage: deploy/deploy_database.sh --yes-delete-the-database --expect-published FILE" in got.stdout
    assert not log.exists(), "docker was called before the arguments were checked"


def test_a_missing_expected_rows_file_is_refused_before_anything(scratch_tree):
    root, env, log = scratch_tree
    got = _deploy(root, env, "--yes-delete-the-database", "--expect-published", "nowhere.json")
    assert got.returncode == 2 and "no expected-rows file at nowhere.json" in got.stdout
    assert not log.exists()


@needs_pg_tools
@pytest.mark.parametrize(
    "held, reaches_down",
    [pytest.param(ROWS, True, id="matching-dump-reaches-down"),
     pytest.param([ROWS[0]], False, id="mismatched-dump-refused-before-down")],
)  # fmt: skip
def test_the_script_gates_the_destructive_step_on_the_rows(scratch_tree, dump_of, tmp_path, held, reaches_down):
    root, env, log = scratch_tree
    dump = dump_of(held)
    shutil.copy2(dump, root / "dump" / "bacatlas.dump")
    expect = tmp_path / "operator" / "expected.json"
    expect.parent.mkdir()
    expect.write_text(json.dumps(_expected(*ROWS)))
    env = {**env, "STUB_DUMP": str(dump), "STUB_CODE_HEAD": _code_head()}
    # relative to where the operator stands, not to the repo the script cd's into
    got = _deploy(root, env, "--yes-delete-the-database", "--expect-published", "expected.json", cwd=expect.parent)
    calls = log.read_text().splitlines()
    if reaches_down:
        assert got.returncode == 97, got.stdout + got.stderr  # the stub's `compose down`: the gate let it through
        assert "the dump's 2 pangenome rows are exactly the expected" in got.stdout
        assert calls[-1].startswith("compose down -v")
    else:
        assert got.returncode == 2, got.stdout + got.stderr
        assert "  ✗ kp-71k: expected, and not held" in got.stdout
        assert "THE DUMP DOES NOT HOLD WHAT IS TO BE SERVED — refusing, and nothing has been deleted." in got.stdout
        assert not any(c.startswith("compose down") for c in calls)
    # the comparison ran in the api IMAGE, offline, never through the running stack
    (api_call,) = [c for c in calls if "bacatlas-api:local" in c and c.startswith("run ")]
    assert "run --rm -i --pull never --network none --entrypoint python" in api_call


@needs_pg_tools
def test_a_container_that_cannot_run_the_check_is_not_reported_as_a_wrong_dump(scratch_tree, dump_of, tmp_path):
    root, env, log = scratch_tree
    dump = dump_of(ROWS)
    shutil.copy2(dump, root / "dump" / "bacatlas.dump")
    expect = tmp_path / "expected.json"
    expect.write_text(json.dumps(_expected(*ROWS)))
    env = {**env, "STUB_DUMP": str(dump), "STUB_CODE_HEAD": _code_head(), "STUB_API_EXIT": "127"}
    got = _deploy(root, env, "--yes-delete-the-database", "--expect-published", str(expect))
    assert got.returncode == 2, got.stdout + got.stderr
    assert "COULD NOT RUN THE ROW CHECK in bacatlas-api:local (exit 127)" in got.stdout
    assert "THE DUMP DOES NOT HOLD" not in got.stdout  # the dump is right; the operator must not be told otherwise
    assert not any(c.startswith("compose down") for c in log.read_text().splitlines())


def test_no_api_image_is_refused_before_anything_is_deleted(scratch_tree, tmp_path):
    root, env, log = scratch_tree
    (root / "dump" / "bacatlas.dump").write_bytes(b"not read: the vintage gate is stubbed")
    expect = tmp_path / "expected.json"
    expect.write_text(json.dumps(ONE))
    env = {**env, "STUB_NO_API_IMAGE": "1", "STUB_CODE_HEAD": _code_head()}
    got = _deploy(root, env, "--yes-delete-the-database", "--expect-published", str(expect))
    assert got.returncode == 2 and "no bacatlas-api:local image on this host" in got.stdout
    assert not any(c.startswith("compose down") for c in log.read_text().splitlines())
