"""What every column actually HOLDS, measured — so a claim can be checked against its contents.

⛔ **This exists because the defects this project keeps producing are not bugs.** They are mismatches
between what a column *means* and what it *holds*, and no assertion knows the difference: `cog_cat`
held `CP` in a column that fits it, `bakta_gene_symbol` lost 125 of 1,180,647 symbols, a Pfam
architecture was truncated at 256 characters into a valid-looking different architecture. Every one
passed every test, because every test compared a value to a value.

⭐ **So measure the contents and put them beside the declaration.** A `String(16)` whose longest
value is 4 characters is fine; one whose longest is 16 is one genome away from silent truncation. A
column documented as a scalar whose values contain a separator is a set. A `measurement()` column
with no NULLs has either never been absent or is being filled with a zero. None of those is decidable
from the code alone, and all three are obvious from a census.

⚠ **A census FLAGS; it does not judge.** Every row it emits is a question, and the answer is in the
column's docstring or in the artifact it was read from. A finding is only a finding once it has been
verified against Postgres or the real files — which is the rule the pass this serves runs under.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from sqlalchemy import Enum as SqlEnum
from sqlalchemy import func, inspect, select, text
from sqlalchemy.orm import Session

from bacatlas_backend.database import Base

#: A `String(n)` whose longest observed value is at least this share of `n` is flagged: the next
#: cohort could exceed it, and the failure mode is a truncation that produces a plausible value.
WIDTH_HEADROOM_WARNING = 0.75

#: Characters that make a scalar string look like a concatenated SET. ⚠ `,` and `;` are the usual
#: joiners; a bare run of single uppercase letters is how Bakta writes COG categories, which is the
#: case that had no separator at all and was therefore invisible.
SET_LIKE_SEPARATORS = (",", ";", "|")


@dataclass
class ColumnObservation:
    """One column, as declared and as filled."""

    table: str
    column: str
    declared_type: str
    nullable: bool
    row_count: int = 0
    null_count: int = 0
    distinct_count: int | None = None
    #: ⚠ How many values the TYPE declares, for an enum only — the denominator `distinct_count` had
    #: no meaning without. A count of distinct values says nothing on its own; against the size of
    #: the declared vocabulary it says whether the column is using it.
    declared_value_count: int | None = None
    minimum: object = None
    maximum: object = None
    maximum_text_length: int | None = None
    minimum_text_length: int | None = None
    declared_text_length: int | None = None
    separator_bearing_count: int = 0
    zero_count: int | None = None
    negative_count: int | None = None
    questions: list = field(default_factory=list)

    @property
    def is_empty(self) -> bool:
        """Whether nothing was examined.

        ⚠ A column in an empty table says nothing, and must not read as a clean bill.
        """
        return self.row_count == 0

    def render(self) -> str:
        """One line: the declaration, what fills it, and every question that raises."""
        parts = [f"{self.table}.{self.column} [{self.declared_type}]"]
        if self.is_empty:
            return parts[0] + "  ⚠ NOT EXAMINED — the table is empty"
        parts.append(f"rows={self.row_count:,}")
        parts.append(f"null={self.null_count:,}")
        if self.maximum_text_length is not None:
            width = f"/{self.declared_text_length}" if self.declared_text_length else ""
            parts.append(f"maxlen={self.maximum_text_length}{width}")
        if self.minimum is not None or self.maximum is not None:
            parts.append(f"range=[{self.minimum}, {self.maximum}]")
        line = "  ".join(parts)
        for question in self.questions:
            line += f"\n      ⚠ {question}"
        return line


def _string_length(column) -> int | None:
    return getattr(column.type, "length", None)


def census_column(session: Session, table, column) -> ColumnObservation:
    """Measure one column. One statement, whatever the type."""
    observation = ColumnObservation(
        table=table.name,
        column=column.name,
        declared_type=str(column.type),
        nullable=column.nullable,
        declared_text_length=_string_length(column),
    )
    observation.row_count = session.execute(select(func.count()).select_from(table)).scalar_one()
    if not observation.row_count:
        return observation

    observation.null_count = session.execute(
        select(func.count()).select_from(table).where(column.is_(None))
    ).scalar_one()

    kind = str(column.type).upper()
    # ⚠ An enum is a user-defined Postgres type, so `length()` does not apply to it — and its
    # interesting fact is its VALUE SET, not its width. Checked first, because SQLAlchemy renders
    # several enum spellings that would otherwise match the text branch below.
    if isinstance(column.type, SqlEnum):
        observation.distinct_count = session.execute(
            select(func.count(func.distinct(column))).select_from(table)
        ).scalar_one()
        observation.declared_value_count = len(column.type.enums or ())
        observation.declared_text_length = None
        return observation
    is_array = "[]" in kind or kind.startswith("ARRAY")
    if is_array:
        # An array's own length distribution is the interesting fact, not its values.
        lengths = session.execute(
            select(
                func.min(func.cardinality(column)),
                func.max(func.cardinality(column)),
            ).select_from(table)
        ).one()
        observation.minimum, observation.maximum = lengths
        return observation

    # ⚠ A blob has no length, no distinct set worth counting and no min/max — but "we did not look
    # at it" and "there was nothing to see" must not read the same, which is this instrument's whole
    # premise. Its size distribution is the one fact that can be wrong: a sprite of 300 bytes is a
    # blank square, and nothing else in the schema would say so.
    if "BYTEA" in kind or "BLOB" in kind:
        observation.minimum, observation.maximum = session.execute(
            select(func.min(func.octet_length(column)), func.max(func.octet_length(column)))
            .select_from(table)
        ).one()
        return observation

    if any(token in kind for token in ("CHAR", "TEXT")):
        observation.minimum_text_length, observation.maximum_text_length = session.execute(
            select(func.min(func.length(column)), func.max(func.length(column))).select_from(table)
        ).one()
        observation.distinct_count = session.execute(
            select(func.count(func.distinct(column))).select_from(table)
        ).scalar_one()
        separators = " OR ".join(
            f"{column.name} LIKE '%{separator}%'" for separator in SET_LIKE_SEPARATORS
        )
        observation.separator_bearing_count = session.execute(
            text(f"SELECT count(*) FROM {table.name} WHERE {separators}")  # noqa: S608
        ).scalar_one()
    elif any(token in kind for token in ("INT", "NUMERIC", "FLOAT", "DOUBLE", "REAL", "SMALL")):
        observation.minimum, observation.maximum = session.execute(
            select(func.min(column), func.max(column)).select_from(table)
        ).one()
        observation.zero_count = session.execute(
            select(func.count()).select_from(table).where(column == 0)
        ).scalar_one()
        observation.negative_count = session.execute(
            select(func.count()).select_from(table).where(column < 0)
        ).scalar_one()
    return observation


def ask_the_questions(observation: ColumnObservation) -> list[str]:
    """Every way this column's contents could disagree with its declaration, as questions."""
    questions: list[str] = []
    if observation.is_empty:
        return questions

    filled = observation.row_count - observation.null_count
    if observation.declared_text_length and observation.maximum_text_length:
        headroom = observation.maximum_text_length / observation.declared_text_length
        # ⚠ A value of FIXED width that exactly fills its column is not near an overflow — a sha256
        # is 64 characters by construction and a strand is 1. Only a column whose values VARY in
        # length and reach the cap is one genome away from a silent truncation, so the flag needs
        # both facts. Without this the rule fires on every hash and every flag character and its
        # signal is lost in its own noise.
        varies = observation.minimum_text_length != observation.maximum_text_length
        if headroom >= WIDTH_HEADROOM_WARNING and varies:
            questions.append(
                f"WIDTH: values range {observation.minimum_text_length}–"
                f"{observation.maximum_text_length} against a declared "
                f"{observation.declared_text_length} ({headroom:.0%}). A longer one truncates "
                "silently into a valid-looking different value."
            )
    if observation.separator_bearing_count and filled:
        share = observation.separator_bearing_count / filled
        questions.append(
            f"SET: {observation.separator_bearing_count:,} of {filled:,} values ({share:.1%}) "
            "contain a separator. Is this column a scalar that holds a set?"
        )
    # ⛔ **An enum column was unauditable by construction until this question existed.** The census
    # measured `distinct_count` for an enum and nothing read it, so a column declaring five states
    # and holding one — a whole vocabulary written as a constant, or never written at all — raised
    # nothing from the instrument built to catch exactly that. It is the `uniref50_impurity` shape
    # (100 % of one value, indistinguishable from the documented meaning) in a different type.
    if observation.declared_value_count and filled:
        if observation.distinct_count == 1:
            questions.append(
                f"ONE VALUE: declares {observation.declared_value_count} states and {filled:,} rows "
                "hold exactly one of them. Is this a constant, or a column nothing writes?"
            )
        elif observation.distinct_count and observation.distinct_count < observation.declared_value_count:
            questions.append(
                f"UNUSED STATES: {observation.distinct_count} of "
                f"{observation.declared_value_count} declared states occur. Which cases never arise, "
                "and is that a fact about the data or about the writer?"
            )
    if not observation.nullable and observation.null_count:
        questions.append("NULL: declared NOT NULL and holds NULLs — impossible; re-read the schema.")
    # ⚠ "Nullable and never NULL" alone fires on almost every optional column and says nothing. The
    # case that matters is a column that is nullable, never NULL, **and holds zeros** — because that
    # is what an absence written as a zero looks like from here, and it is the exact shape of the
    # NULL-vs-0.0 trap this schema documents everywhere.
    if (
        observation.nullable
        and observation.null_count == 0
        and filled
        and observation.zero_count
    ):
        share = observation.zero_count / filled
        questions.append(
            f"ABSENCE: nullable, never NULL, and {observation.zero_count:,} of {filled:,} values "
            f"({share:.1%}) are zero. Is a zero here *measured zero* or an absence written as one?"
        )
    if observation.negative_count:
        questions.append(
            f"SIGN: {observation.negative_count:,} negative values. Is a negative a measurement "
            "(a signed gap length), a sentinel (-1, -32768, -1e12), or a fault?"
        )
    return questions


def run_column_census(session: Session, *, tables: list[str] | None = None) -> list[ColumnObservation]:
    """Every column of every table, measured — and the coverage is the return length.

    ⛔ Returns observations for EMPTY tables too, flagged as not examined. A census that silently
    omitted them would report full coverage over the half of the schema that happens to be loaded.
    """
    inspector = inspect(session.get_bind())
    present = set(inspector.get_table_names())
    out: list[ColumnObservation] = []
    for table in Base.metadata.sorted_tables:
        if tables and table.name not in tables:
            continue
        if table.name not in present:
            continue
        for column in table.columns:
            observation = census_column(session, table, column)
            observation.questions = ask_the_questions(observation)
            out.append(observation)
    return out


def render_census(observations: list[ColumnObservation], *, only_questions: bool = True) -> str:
    """The census as lines, with its own coverage stated first."""
    examined = [row for row in observations if not row.is_empty]
    empty_tables = sorted({row.table for row in observations if row.is_empty})
    lines = [
        f"column census: {len(observations):,} columns over "
        f"{len({row.table for row in observations})} tables; "
        f"{len(examined):,} examined against real rows",
    ]
    if empty_tables:
        lines.append(f"  ⚠ NOT EXAMINED (empty tables): {', '.join(empty_tables)}")
    flagged = [row for row in observations if row.questions]
    lines.append(f"  {len(flagged):,} columns raise a question:")
    for row in observations if not only_questions else flagged:
        lines.append(f"    {row.render()}")
    return "\n".join(lines)
