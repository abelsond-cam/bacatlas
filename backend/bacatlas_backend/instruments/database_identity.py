"""Are two BacAtlas databases the same database? Table by table, every row, as a multiset.

⛔ **This is the gate for a change that claims to move no number.** A memory fix, a re-plumbed
loader, a chunked fold: each is argued bit-exact on the frames it touches, and the argument stops at
the function boundary. What reaches a reader is rows, so the claim is checked where it lands — build
the database twice by the SAME sequence of ingest commands, once with the untouched code and once
with the change, and compare everything.

⭐ **Everything means everything.** Every table in the schema, every column that is not a timestamp,
surrogate ids included: two builds that insert in the same order assign the same ids, so an id that
differs is a difference worth seeing rather than noise to be joined away. Rows are compared as a
multiset (an md5 per row, aggregated in sorted order), never by physical position — the heap order
of the big locus tables varies between two runs of the SAME code, which was measured, not assumed.

⚠ **It asserts its own coverage before it reports agreement.** Same tables, same columns and types,
and the tables a loaded catalogue cannot be without are not empty. *Not looked at* and *no
difference* must not print the same line: two empty databases agree perfectly.

⚠ **What is excluded is excluded BY NAME and printed.** Timestamp columns by type, and
`BUILD_STAMPS` — columns that record which commit ran the build and so differ between any two
builds by construction. Nothing else is skipped, and a column added to that set is a decision.

    python -m bacatlas_backend.instruments.database_identity <reference url> <candidate url>

exits 0 when no table differs and no required table is empty, 1 otherwise, 2 when the two
databases do not even hold the same tables and columns.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field

from sqlalchemy import create_engine, text
from sqlalchemy.engine import Engine

#: Never compared: alembic's own bookkeeping is one row naming the migration, and two databases at
#: different migrations fail the column comparison long before it could matter.
SKIPPED_TABLES = frozenset({"alembic_version"})

#: Columns that record WHICH BUILD wrote the row. `pangenome.git_sha` is nuna's HEAD at ingest
#: time, so it differs between a build before a commit and a build after it — which is every use
#: this instrument has.
BUILD_STAMPS = frozenset({("pangenome", "git_sha")})

#: Tables a database holding a published catalogue with a projection cannot have empty. A
#: comparison of two databases that both failed to load would otherwise report perfect agreement.
MUST_HOLD_ROWS = frozenset({
    "genome", "genome_genes", "genome_contig", "genome_locus_assignment", "locus", "locus_arrangement",
    "locus_offset_occupant", "locus_annotation_entry", "locus_similarity", "locus_nearest_locus",
    "intergenic_gap", "pangenome", "pangenome_propagation_rate", "pangenome_genome_locus_count",
    "projected_genome", "projected_gene_placement",
})  # fmt: skip

#: How many differing rows of each side are quoted per table — enough to see what kind of
#: difference it is, never the whole table.
EXAMPLES_PER_SIDE = 2


class CoverageMismatch(RuntimeError):
    """The two databases do not hold the same tables and columns, so no row comparison means anything."""


@dataclass
class TableComparison:
    """One table, both sides."""

    table: str
    reference_rows: int
    candidate_rows: int
    compared_columns: tuple[str, ...]
    excluded_columns: tuple[str, ...]
    same: bool
    only_in_reference: int = 0
    only_in_candidate: int = 0
    reference_examples: tuple[str, ...] = ()
    candidate_examples: tuple[str, ...] = ()


@dataclass
class IdentityReport:
    """Every table's comparison, and the two ways the whole can fail."""

    tables: list[TableComparison] = field(default_factory=list)
    column_count: int = 0
    empty_required: list[str] = field(default_factory=list)

    @property
    def differing(self) -> list[TableComparison]:
        """The tables whose row multisets are not equal."""
        return [table for table in self.tables if not table.same]

    @property
    def identical(self) -> bool:
        """No table differs AND no required table is empty — agreement alone is not enough."""
        return not self.differing and not self.empty_required

    def render(self) -> str:
        """One line per table, the differing ones with examples, then the count that decides."""
        lines = [
            f"coverage: {len(self.tables)} tables, {self.column_count} columns, "
            "identical table and column lists in both"
        ]
        for table in self.tables:
            minus = f", minus {','.join(table.excluded_columns)}" if table.excluded_columns else ""
            if table.same:
                lines.append(
                    f"  same    {table.table:34s} {table.reference_rows:>11,} rows  "
                    f"({len(table.compared_columns)} columns{minus})"
                )
                continue
            lines.append(
                f"  DIFFERS {table.table:34s} {table.reference_rows:,} vs {table.candidate_rows:,} rows; "
                f"{table.only_in_reference:,} only in the reference, {table.only_in_candidate:,} only in the candidate"
            )
            lines += [f"      reference: {row[:300]}" for row in table.reference_examples]
            lines += [f"      candidate: {row[:300]}" for row in table.candidate_examples]
        lines += [f"  EMPTY   {name}: a table that must hold rows" for name in self.empty_required]
        total = sum(table.reference_rows for table in self.tables)
        lines.append(
            f"\n{len(self.tables)} tables, {total:,} reference rows: {len(self.differing)} table(s) differ, "
            f"{len(self.empty_required)} required table(s) empty"
        )
        return "\n".join(lines)


def _columns(engine: Engine, schema: str) -> dict[str, list[tuple[str, str]]]:
    """`table → [(column, type)]` in declaration order, for every table of one schema."""
    with engine.connect() as connection:
        rows = connection.execute(
            text(
                "select table_name, column_name, data_type from information_schema.columns "
                "where table_schema = :schema order by table_name, ordinal_position"
            ),
            {"schema": schema},
        ).all()
    found: dict[str, list[tuple[str, str]]] = {}
    for table, column, data_type in rows:
        if table not in SKIPPED_TABLES:
            found.setdefault(table, []).append((column, data_type))
    return found


def _row_expression(columns: tuple[str, ...]) -> str:
    return "md5(row(" + ", ".join(f'"{column}"' for column in columns) + ")::text)"


def _digest(engine: Engine, schema: str, table: str, columns: tuple[str, ...]) -> tuple[int, str]:
    """Row count, and one md5 over the sorted per-row md5s — equal iff the two multisets are."""
    query = (
        "select count(*), md5(coalesce(string_agg(h, '' order by h), '')) "
        f'from (select {_row_expression(columns)} as h from "{schema}"."{table}") rows'
    )
    with engine.connect() as connection:
        count, digest = connection.execute(text(query)).one()
    return int(count), digest


def _rows_by_hash(engine: Engine, schema: str, table: str, columns: tuple[str, ...]) -> dict[str, tuple[int, str]]:
    """`row md5 → (how many rows carry it, the row as text)` — only read for a table that differs."""
    listed = ", ".join(f'"{column}"' for column in columns)
    query = (
        f"select {_row_expression(columns)} as h, count(*), min(row({listed})::text) "
        f'from "{schema}"."{table}" group by h'
    )
    with engine.connect() as connection:
        return {digest: (int(count), row) for digest, count, row in connection.execute(text(query)).all()}


def compare_databases(
    reference: Engine,
    candidate: Engine,
    *,
    reference_schema: str = "public",
    candidate_schema: str = "public",
    must_hold_rows: frozenset[str] = MUST_HOLD_ROWS,
) -> IdentityReport:
    """Compare every table of two schemas; raise `CoverageMismatch` if they are not the same shape.

    ⚠ A duplicated row is a difference: the multiset counts it. The per-row listing of a differing
    table therefore counts multiplicity too, so `only_in_*` is rows, not distinct rows.
    """
    ours, theirs = _columns(reference, reference_schema), _columns(candidate, candidate_schema)
    if set(ours) != set(theirs):
        raise CoverageMismatch(f"the two databases hold different tables: {sorted(set(ours) ^ set(theirs))}")
    if not ours:
        raise CoverageMismatch(f"schema {reference_schema!r} holds no table at all — nothing was compared")
    for table in ours:
        if ours[table] != theirs[table]:
            raise CoverageMismatch(
                f"`{table}` has different columns or types: {sorted(set(ours[table]) ^ set(theirs[table]))}"
            )

    report = IdentityReport(column_count=sum(len(columns) for columns in ours.values()))
    for table in sorted(ours):
        compared = tuple(
            column for column, data_type in ours[table]
            if "timestamp" not in data_type and (table, column) not in BUILD_STAMPS
        )  # fmt: skip
        excluded = tuple(column for column, _ in ours[table] if column not in compared)
        if not compared:
            raise CoverageMismatch(f"`{table}` has no comparable column — every one is a timestamp or a build stamp")
        reference_rows, reference_digest = _digest(reference, reference_schema, table, compared)
        candidate_rows, candidate_digest = _digest(candidate, candidate_schema, table, compared)
        if table in must_hold_rows and (reference_rows == 0 or candidate_rows == 0):
            report.empty_required.append(table)
        same = reference_rows == candidate_rows and reference_digest == candidate_digest
        comparison = TableComparison(
            table=table, reference_rows=reference_rows, candidate_rows=candidate_rows,
            compared_columns=compared, excluded_columns=excluded, same=same,
        )  # fmt: skip
        if not same:
            mine = _rows_by_hash(reference, reference_schema, table, compared)
            yours = _rows_by_hash(candidate, candidate_schema, table, compared)
            surplus_mine = {h: n - yours.get(h, (0, ""))[0] for h, (n, _) in mine.items()}
            surplus_yours = {h: n - mine.get(h, (0, ""))[0] for h, (n, _) in yours.items()}
            comparison.only_in_reference = sum(n for n in surplus_mine.values() if n > 0)
            comparison.only_in_candidate = sum(n for n in surplus_yours.values() if n > 0)
            comparison.reference_examples = tuple(
                mine[h][1] for h in sorted(h for h, n in surplus_mine.items() if n > 0)[:EXAMPLES_PER_SIDE]
            )
            comparison.candidate_examples = tuple(
                yours[h][1] for h in sorted(h for h, n in surplus_yours.items() if n > 0)[:EXAMPLES_PER_SIDE]
            )
        report.tables.append(comparison)
    return report


def main(argv: list[str] | None = None) -> int:
    """Compare two databases by URL; 0 identical, 1 differing or empty, 2 not comparable."""
    arguments = sys.argv[1:] if argv is None else argv
    if len(arguments) != 2:
        usage = "python -m bacatlas_backend.instruments.database_identity <reference url> <candidate url>"
        print(f"usage: {usage}", file=sys.stderr)
        return 2
    try:
        report = compare_databases(create_engine(arguments[0], future=True), create_engine(arguments[1], future=True))
    except CoverageMismatch as error:
        print(f"NOT COMPARED: {error}", file=sys.stderr)
        return 2
    print(f"excluded by name: {sorted(BUILD_STAMPS)} (build stamps) and every timestamp column")
    print(report.render())
    return 0 if report.identical else 1


if __name__ == "__main__":
    raise SystemExit(main())
