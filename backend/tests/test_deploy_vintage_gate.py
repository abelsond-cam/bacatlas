"""`deploy_database.sh` reads the migration head out of the FILES. This is the test that it can.

⛔ **The gate was wrong for five days and nothing could see it.** `deploy_database.sh` refuses to
restore a dump whose `alembic_version` differs from the code's head — the fix for the 2026-10-22
incident where a four-month-old dump restored cleanly, answered `/health` happily, and 500'd every
locus request. It computes the code's head with `grep` over the migration files rather than with
alembic, deliberately: the gate runs on the instance before anything is installed, where there is
no alembic and no database to point it at.

That grep matched **single-quoted** revision ids, because alembic's template writes single quotes
and all 24 migrations then on disk did. `b1f39c8b40e9` was hand-written with double quotes. The
pattern missed it, `comm -23` returned its PARENT as the head, and the gate refused every dump taken
afterwards — instructing the operator to take a fresh dump that would mismatch identically forever.
One more migration and it fails differently: two unmatched ids make the head non-unique and the gate
exits on "no single head". Either way the only path that can deploy data is closed, and
`test_migration_round_trip` cannot tell, because it drives alembic and never reads this shell.

⭐ **So this module runs the gate's OWN pipeline, lifted verbatim out of the script**, rather than a
second copy of it. A copy is how the two drift apart again — the same reason `cohort_label_refusal`
is one function with two callers rather than two conditions.

Two tests, and they fail for different reasons on purpose:

* against the real tree, the gate must agree with `alembic heads`;
* against a synthetic tree written in three quoting styles, it must find the head anyway. The second
  is the one that would have caught this before `b1f39c8b40e9` was written, because it does not
  depend on any style happening to be present today.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

BACKEND_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = BACKEND_ROOT.parent
GATE = REPO_ROOT / "deploy" / "deploy_database.sh"

pytestmark = pytest.mark.skipif(
    not GATE.exists(), reason=f"{GATE} is not in this checkout"
)


def _gate_head_pipeline() -> tuple[str, str]:
    """The `versions=` and `code_rev=` statements of the gate, verbatim.

    ⚠ Lifted from the script rather than retyped. If the gate is rewritten in a form these two
    anchors cannot find, this raises instead of silently testing nothing — *not looked at* and
    *no difference* must not look the same (`a-zero-needs-its-denominator`).
    """
    lines = GATE.read_text(encoding="utf-8").splitlines()

    versions = next((line for line in lines if line.startswith("versions=")), None)
    assert versions is not None, (
        f"no `versions=` assignment in {GATE}. The gate has been rewritten and this test is now "
        "measuring nothing — re-anchor it on whatever names the migrations directory."
    )

    start = next((i for i, line in enumerate(lines) if line.startswith("code_rev=")), None)
    assert start is not None, (
        f"no `code_rev=` assignment in {GATE}. The gate has been rewritten and this test is now "
        "measuring nothing — re-anchor it on whatever computes the code's head."
    )
    statement = [lines[start]]
    while statement[-1].rstrip().endswith("\\"):
        statement.append(lines[start + len(statement)])
    return versions, "\n".join(statement)


def _run_gate_head(versions_statement: str) -> str:
    """Run the gate's head computation with `versions=` pointed wherever the caller says."""
    _, code_rev = _gate_head_pipeline()
    script = f'set -u\n{versions_statement}\n{code_rev}\nprintf "%s" "$code_rev"\n'
    # ⚠ bash, not sh: `comm -23 <(…) <(…)` is process substitution, which dash does not have. The
    # instance runs this file through `bash`, so testing it under anything else tests another
    # program.
    finished = subprocess.run(
        ["bash", "-c", script], cwd=REPO_ROOT, capture_output=True, text=True, check=False
    )
    assert finished.returncode == 0, f"the gate's own pipeline failed:\n{finished.stderr}"
    return finished.stdout.strip()


def _alembic_head() -> str:
    finished = subprocess.run(
        [sys.executable, "-m", "alembic", "heads"],
        cwd=BACKEND_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    if finished.returncode != 0:
        pytest.skip(f"alembic heads failed in this interpreter: {finished.stderr.strip()}")
    heads = re.findall(r"\b[0-9a-f]{12}\b", finished.stdout)
    assert len(heads) == 1, (
        f"alembic reports {len(heads)} heads ({heads}); the gate assumes one and so does this test."
    )
    return heads[0]


def test_the_deploy_gate_computes_the_same_head_alembic_does():
    """The one assertion the 2026-10-27 bug violated: grep's answer == alembic's answer."""
    versions, _ = _gate_head_pipeline()
    assert _run_gate_head(versions) == _alembic_head()


#: `(name, how the revision line is written)` — the three forms a migration file has used or could.
#: ⚠ The bare form is not hypothetical politeness: it is what a hand-edit that drops the quotes
#: leaves behind, and the gate must survive it for the same reason it must survive double quotes.
QUOTING_STYLES = (
    ("single, alembic's own template", "'{revision}'"),
    ("double, as b1f39c8b40e9 was hand-written", '"{revision}"'),
    ("bare, as a hand-edit can leave it", "{revision}"),
)


@pytest.mark.parametrize(("style", "template"), QUOTING_STYLES, ids=[s for s, _ in QUOTING_STYLES])
def test_the_gate_finds_the_head_whatever_quotes_the_revision_wears(tmp_path, style, template):
    """⭐ The test that would have caught this BEFORE the double-quoted migration was written.

    A two-migration chain whose head is written in `style` and whose parent is written in alembic's
    own. Pin the pattern's generality, not the styles that happen to be on disk today — a test over
    the current tree goes green again the moment the offending file is reformatted, which is exactly
    the fix this comment argues against.
    """
    parent, head = "aaaaaaaaaaaa", "bbbbbbbbbbbb"
    versions = tmp_path / "versions"
    versions.mkdir()
    (versions / f"{parent}_the_parent.py").write_text(
        f"revision: str = '{parent}'\ndown_revision: str | None = None\n", encoding="utf-8"
    )
    (versions / f"{head}_the_head.py").write_text(
        f"revision: str = {template.format(revision=head)}\n"
        f"down_revision: str | None = {template.format(revision=parent)}\n",
        encoding="utf-8",
    )
    assert _run_gate_head(f"versions={versions}") == head, (
        f"the gate cannot see a revision written with {style}. It would then report this chain's "
        f"head as {parent!r} — the PARENT — and refuse every dump taken from a database at "
        f"{head!r}, telling the operator to take a fresh dump that mismatches identically forever."
    )
