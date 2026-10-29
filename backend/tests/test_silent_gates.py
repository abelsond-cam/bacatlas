"""⛔ A check that SKIPS has not passed, and a run that prints `36 skipped` has not said which 36.

This matters more the further the real work moves off this machine. The 71k programme puts the
ingest on CSD3 and reduces the laptop to a conduit for a dump — and that laptop still runs this
suite, which can report green having checked the parity oracles, the route costs and the
column-meaning audit **not at all**.

Two mechanisms, and this file is what stops either becoming decoration:

* **`require_oracle` FAILS** where an independent oracle is absent. An oracle is the only
  independent witness a piece of code has; `test_gene_blocks`' differential against the published
  `app.js` is the gene-block codec's *only* one, and it skipped silently when `node` was missing.
* **the skip census** prints every skip, grouped by reason, on **every** run. `-rs` would do it on
  request, and nobody requests it.

⚠ `BACATLAS_ALLOW_UNORACLED=1` is the typed opt-out, and it is here so a machine that genuinely
cannot run an oracle says so once rather than every suite everywhere quietly measuring less.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from tests.conftest import require_oracle

BACKEND = Path(__file__).resolve().parents[1]

#: ⭐ The ONE case that is run end to end, as a subprocess, because the behaviour under test is
#: what pytest reports — and a test asserting on `pytest.fail` from inside pytest asserts on the
#: exception, not on the report.
GENE_BLOCK_ORACLE = "tests/test_gene_blocks.py"


def test_a_present_oracle_is_a_no_op():
    require_oracle(True, "nothing is missing", supplies="nothing")


def test_an_absent_oracle_FAILS_and_names_what_it_would_have_checked():
    # ⚠ `pytest.fail.Exception`, not `Exception`: both `Failed` and `Skipped` derive from
    # BaseException, so `pytest.raises(Exception)` catches NEITHER — and a test written that way
    # does not fail, it is reported as skipped by the outcome it failed to catch. Exactly the
    # silent gate this file exists to close, reproduced inside the file closing it.
    with pytest.raises(pytest.fail.Exception) as refusal:
        require_oracle(False, "node is not on PATH", supplies="the codec's differential oracle")
    message = str(refusal.value)
    assert "ORACLE UNAVAILABLE" in message
    assert "node is not on PATH" in message
    assert "the codec's differential oracle" in message, "a refusal that does not say what was lost"
    assert "BACATLAS_ALLOW_UNORACLED" in message, "a refusal with no way past it is a wall, not a gate"


def test_the_opt_out_is_TYPED_and_defaults_off(monkeypatch):
    """⛔ Default off. A default-on escape hatch is the silent skip with extra steps."""
    import tests.conftest as conftest

    monkeypatch.setattr(conftest, "ALLOW_UNORACLED", True)
    with pytest.raises(pytest.skip.Exception):
        conftest.require_oracle(False, "node is not on PATH", supplies="x")

    monkeypatch.setattr(conftest, "ALLOW_UNORACLED", False)
    with pytest.raises(pytest.fail.Exception):
        conftest.require_oracle(False, "node is not on PATH", supplies="x")


def _run(*args: str, **env: str) -> subprocess.CompletedProcess:
    import os

    return subprocess.run(
        [sys.executable, "-m", "pytest", GENE_BLOCK_ORACLE, "-q", "-p", "no:cacheprovider"],
        cwd=BACKEND, capture_output=True, text=True, check=False,
        # ⛔ A PATH with no `node`. This is the real condition, not a monkeypatch of `shutil.which`:
        # the thing being tested is whether the SUITE goes red, and only a real run reports that.
        env={**os.environ, "PATH": "/usr/bin:/bin", **env},
    )


def test_the_suite_goes_RED_when_node_is_removed_from_PATH():
    """⭐ The gate, stated as the plan states it. Before this, it went green."""
    result = _run()
    assert result.returncode != 0, "the suite passed without its only independent oracle"
    assert "ORACLE UNAVAILABLE" in result.stdout
    assert " 0 skipped" not in result.stdout


def test_the_opt_out_turns_it_back_into_a_skip_and_the_census_NAMES_it():
    """⚠ Skipping is allowed — silently skipping is not. The census is what makes the difference."""
    result = _run(BACATLAS_ALLOW_UNORACLED="1")
    assert result.returncode == 0, result.stdout[-2000:]
    assert "what this run did NOT check" in result.stdout, "the census did not print"
    assert "ORACLE UNAVAILABLE" in result.stdout
