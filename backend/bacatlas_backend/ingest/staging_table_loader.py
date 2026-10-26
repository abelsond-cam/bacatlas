"""Bulk loading — ``COPY``, and the value coercions that must happen at the boundary.

⛔ **``COPY``, never ``executemany``, and never the ORM.** 1.02 M gene rows at 100 genomes becomes
412 M at 80,000, and the difference between the three is two orders of magnitude, not a constant
factor. The ORM path also builds a Python object per row, which is where an ingest quietly acquires
a memory ceiling.

⚠ **Text-format COPY, deliberately, for now.** Binary is roughly 2–3× faster but requires the
column type OIDs to be resolved and pinned per table — including the eleven user-defined enum types
— which is a plumbing surface that fails silently by writing a plausible wrong value. Text format
lets Postgres parse each field into the column's own type, so an enum, an array and a JSONB column
all arrive correctly with no type table to keep in step. Measured cost is what should change this,
not preference: if a rung of the scale ladder shows the load is COPY-bound, add a binary path
*behind this same function* and cross-check the two on one table.

⛔ **NaN is coerced here, at the boundary.** Every ``measurement()`` column has a CHECK forbidding
NaN, and the CHECK is the backstop — this is the intent. Postgres defines ``NaN = NaN`` as TRUE, so
a NaN that got in would compare equal to itself and be indistinguishable from a real value in every
later query.
"""

from __future__ import annotations

import enum
import math
from collections.abc import Iterable, Iterator, Sequence
from typing import Any

from sqlalchemy import LargeBinary, Table, delete
from sqlalchemy.orm import Session

from bacatlas_backend.instruments.sql_cost_oracle import record_driver_statement

#: How many rows to hand the driver between flushes. Large enough that the per-batch overhead is
#: noise, small enough that a failing row is found in seconds rather than after a whole genome.
COPY_BATCH_ROWS = 50_000


def coerce_value(value: Any) -> Any:
    """One row value → something the driver can write, or ``None``.

    Four coercions, each of which is a trap somewhere else in this codebase:

    - **NaN → NULL.** *Not measured* and *measured zero* are different findings.
    - **A Python enum → its NAME**, because SQLAlchemy's ``Enum`` persists ``.name`` by default and
      the CHECK constraints were written against that (``exclusivity_form = 'EXCLUSION'``). Passing
      ``.value`` writes ``'exclusion'``, which the enum type rejects — loudly, fortunately.
    - **A numpy scalar → its Python equivalent.** ``numpy.int64`` has no text dumper registered, and
      the failure is an adaptation error a long way from the column that caused it.
    - **A list is passed through**, so psycopg dumps it to an array literal. ⚠ An empty list and
      ``None`` are *different* — ``{}`` is "annotated with nothing", ``NULL`` is "not annotated" —
      so nothing here collapses one into the other.
    """
    if value is None:
        return None
    if isinstance(value, float):
        return None if math.isnan(value) else value
    if isinstance(value, enum.Enum):
        return value.name
    if isinstance(value, (list, tuple)):
        return [coerce_value(item) for item in value]
    # numpy scalars and pandas NA, without importing numpy into the serving path.
    item = getattr(value, "item", None)
    if item is not None and type(value).__module__.startswith("numpy"):
        return coerce_value(item())
    if value is not value:  # pandas NaT / NA compare unequal to themselves
        return None
    return value


def _refuse_text_in_a_binary_column(
    table: Table, columns: Sequence[str], binary_positions: Sequence[int], row: Sequence[Any]
) -> None:
    r"""⛔ A `str` handed to a `bytea` column is reinterpreted, not stored. Measured 2026-10-26.

    `coerce_value` passes a `str` straight through, psycopg dumps it as text, and bytea's input
    function then reads it as an escaped literal: 34 characters of ``\x4141…`` become **16 bytes**,
    a legal block length, so the column's own CHECK is satisfied and nothing raises. That is the
    exact failure mode `copy_rows`' docstring rejects the binary path for — *"fails silently by
    writing a plausible wrong value"* — reached through the text path instead.

    ⚠ **Checked on the FIRST row only, deliberately.** This is a programming error, uniform across
    rows by construction, and the gene layer will push 380 M rows through here at 72,801 genomes;
    a per-row `isinstance` on that is a cost with no new information after the first row.
    """
    for position in binary_positions:
        value = row[position]
        if isinstance(value, str):
            raise TypeError(
                f"{table.name}.{columns[position]} is a bytea column and was given a str "
                f"({value[:24]!r}…). Pass `bytes`: a str is reinterpreted by bytea's input "
                "function, so it is stored as DIFFERENT bytes at a plausible length and no "
                "constraint can see it."
            )


def _batched(rows: Iterable[Sequence[Any]], size: int) -> Iterator[list[Sequence[Any]]]:
    batch: list[Sequence[Any]] = []
    for row in rows:
        batch.append(row)
        if len(batch) >= size:
            yield batch
            batch = []
    if batch:
        yield batch


def copy_rows(
    session: Session,
    table: Table,
    columns: Sequence[str],
    rows: Iterable[Sequence[Any]],
    *,
    batch_rows: int = COPY_BATCH_ROWS,
) -> int:
    """``COPY`` rows into ``table``, returning how many were written.

    ``rows`` is an iterable of positional tuples matching ``columns``. It is consumed lazily, so a
    generator over a parquet file never materialises the whole frame.

    ⚠ The count returned is what this function WROTE, not what the table now holds. Ingest verifies
    against the source's own row count; a loader that reported the table's total would confirm
    itself.
    """
    unknown = [name for name in columns if name not in table.c]
    if unknown:
        raise KeyError(
            f"{table.name} has no column(s) {unknown!r}. Known columns: {sorted(table.c.keys())}. "
            "A COPY into a mis-named column is a silent no-op in some drivers, so this is checked "
            "before a single row is sent."
        )

    column_list = ", ".join(f'"{name}"' for name in columns)
    statement = f'COPY "{table.name}" ({column_list}) FROM STDIN'

    # ⛔ Straight to the DBAPI cursor, because `COPY … FROM STDIN` has no SQLAlchemy expression —
    # which means SQLAlchemy's `after_cursor_execute` never fires for it. The cost oracle would
    # therefore score a 5,000-row load as zero statements, so this path REPORTS ITSELF. See
    # `instruments/sql_cost_oracle.py` — the blind spot is named there, not left to be found.
    # Which of `columns` are bytea, resolved once so the first-row check costs one tuple walk.
    binary_positions = tuple(
        index for index, name in enumerate(columns) if isinstance(table.c[name].type, LargeBinary)
    )

    driver_connection = session.connection().connection.driver_connection
    written = 0
    first_row_checked = False
    with driver_connection.cursor() as cursor, cursor.copy(statement) as copy:
        for batch in _batched(rows, batch_rows):
            for row in batch:
                if binary_positions and not first_row_checked:
                    _refuse_text_in_a_binary_column(table, columns, binary_positions, row)
                    first_row_checked = True
                copy.write_row(tuple(coerce_value(value) for value in row))
            written += len(batch)
    record_driver_statement("copy", statement, written)
    return written


def replace_rows_for(
    session: Session,
    table: Table,
    columns: Sequence[str],
    rows: Iterable[Sequence[Any]],
    *,
    where,
    batch_rows: int = COPY_BATCH_ROWS,
) -> tuple[int, int]:
    """Delete the rows matching ``where``, then COPY the replacements. Returns ``(deleted, written)``.

    ⭐ **This is what makes an ingest re-runnable at a useful granularity.** A load that fails on
    genome 200 of 280 must be resumable without dropping the database, and a load that is re-run
    after a bug fix must not double every row. Scoping delete-then-load to one genome (or one
    pangenome) keeps the blast radius equal to the unit of work, and the caller's transaction makes
    the pair atomic: a crash between the two leaves the previous rows intact.

    ⚠ ``where`` must actually be scoped. An unbounded predicate here silently becomes a table wipe,
    so a ``where`` of ``None`` is refused rather than treated as "everything".
    """
    if where is None:
        raise ValueError(
            "replace_rows_for needs a scoped predicate. An unscoped delete here would empty the "
            "table, and the caller would see only a successful load."
        )
    deleted = session.execute(delete(table).where(where)).rowcount
    written = copy_rows(session, table, columns, rows, batch_rows=batch_rows)
    return deleted, written
