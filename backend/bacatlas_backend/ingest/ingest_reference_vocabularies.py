"""Load the public reference vocabularies the catalogue joins against.

⛔ **Global, not per pangenome.** A Pfam family is the same family in every species and every model,
so this runs once and is idempotent — nothing here is keyed by `pangenome_id`, and re-running it
after a new catalogue lands is a no-op rather than a duplicate.

⚠ **Read through `nuna.tl.locus_browser.vendor_reference`, never by parsing the gzip here.** That
module is deliberately the single reader of `pfam_names.tsv.gz`: three call sites wanted the table
and three private parsers is exactly the drift that lets two of them disagree about which column the
clan is in. Ingest becomes a fourth reader only through the same door.
"""

from __future__ import annotations

from dataclasses import dataclass

from sqlalchemy import func, select
from sqlalchemy.dialects.postgresql import insert
from sqlalchemy.orm import Session

from bacatlas_backend.models.reference_vocabulary import EnzymeClass, KeggOrthology, PfamFamily

#: ⛔ The vendored table's field name → our column. Read POSITIONALLY would be a silent
#: transposition waiting to happen — `clan` and `clan_id` are adjacent, one is an accession and the
#: other a readable id, and swapping them produces a page that looks entirely plausible. The
#: positions are taken from `vendor_reference.FIELDS` at load time, so a reordered table follows.
PFAM_COLUMN_FOR_FIELD = {
    "short_name": "short_name",
    "description": "description",
    "interpro": "interpro_accession",
    "interpro_name": "interpro_name",
    "clan": "clan_accession",
    "clan_id": "clan_name",
}

#: Every column the upsert writes — the mapping's values, in one place, so the SET clause below and
#: the row builder cannot drift apart.
PFAM_COLUMNS = tuple(PFAM_COLUMN_FOR_FIELD.values())


@dataclass(frozen=True)
class VocabularyReport:
    """What was loaded, for the reconciliation the CLI prints."""

    pfam_families_read: int
    pfam_families_in_table: int
    #: ⚠ Fields named above that the vendored table no longer has — reported rather than swallowed,
    #: because a column that quietly becomes empty everywhere reads as "Pfam has no clans".
    unmapped_fields: tuple[str, ...] = ()
    enzyme_classes_read: int = 0
    enzyme_classes_in_table: int = 0
    kegg_orthologies_read: int = 0
    kegg_orthologies_in_table: int = 0

    def render(self) -> str:
        """One line per vocabulary, naming any Pfam field the vendored table did not map.

        ⛔ **A zero is printed, never skipped.** Each of these references is optional — the page
        falls back to a bare accession — so an absent vendored table loads nothing and raises
        nothing. Printing the count is the only thing standing between that and a reader inferring
        "EC has no names" from a blank column three screens later.
        """
        lines = [
            f"pfam reference: {self.pfam_families_read:,} families read, "
            f"{self.pfam_families_in_table:,} in the table"
        ]
        if self.unmapped_fields:
            lines[0] += f" — ⚠ absent from the vendored table: {', '.join(self.unmapped_fields)}"
        lines.append(
            f"ec reference: {self.enzyme_classes_read:,} codes read, "
            f"{self.enzyme_classes_in_table:,} in the table"
            + ("" if self.enzyme_classes_read else " — ⚠ no vendored table; EC codes stay unnamed")
        )
        lines.append(
            f"kegg reference: {self.kegg_orthologies_read:,} KOs read, "
            f"{self.kegg_orthologies_in_table:,} in the table"
            + (
                " — ⛔ licence unresolved; internal deployments only"
                if self.kegg_orthologies_read
                else " — ⚠ no vendored table; KO ids stay linked and unnamed"
            )
        )
        return "\n".join(lines)


def load_pfam_reference(session: Session) -> VocabularyReport:
    """Upsert every Pfam-A family from the vendored reference.

    Returns counts rather than raising on an absent table: the reference is **optional everywhere**
    in nuna — the page falls back to bare accessions and the audit's clan map falls back to identity
    — so an absent table must degrade the chips, not fail the load. ⚠ It still reports zero, so a
    silently missing reference is visible in the run output rather than inferred from blank chips
    three screens later.
    """
    from nuna.tl.locus_browser.vendor_reference import FIELDS, pfam_reference

    # ⛔ Resolve each field's POSITION from nuna's own FIELDS tuple, by name. A field nuna has
    # dropped simply yields an empty column; one it has added is ignored until named above.
    position_of = {field: FIELDS.index(field) for field in PFAM_COLUMN_FOR_FIELD if field in FIELDS}
    missing = sorted(set(PFAM_COLUMN_FOR_FIELD) - set(position_of))

    reference = pfam_reference()
    rows = [
        {
            "pfam_accession": accession,
            # Padded by the reader, so a short row degrades to empty strings rather than raising.
            **{
                column: _text(fields, position_of.get(field, -1))
                for field, column in PFAM_COLUMN_FOR_FIELD.items()
            },
        }
        for accession, fields in reference.items()
    ]

    if rows:
        for start in range(0, len(rows), 5_000):
            batch = rows[start : start + 5_000]
            statement = insert(PfamFamily).values(batch)
            session.execute(
                statement.on_conflict_do_update(
                    index_elements=[PfamFamily.pfam_accession],
                    set_={column: statement.excluded[column] for column in PFAM_COLUMNS},
                )
            )

    in_table = session.execute(select(func.count()).select_from(PfamFamily)).scalar_one()
    return VocabularyReport(
        pfam_families_read=len(rows),
        pfam_families_in_table=int(in_table),
        unmapped_fields=tuple(missing),
    )


def _text(fields: tuple[str, ...], index: int) -> str:
    """One field, as text, tolerating an absent field (`-1`) or a row padded shorter than expected."""
    return "" if index < 0 or index >= len(fields) else (fields[index] or "")


def load_enzyme_reference(session: Session) -> tuple[int, int]:
    """Upsert every EC code ExPASy names, at every level. Returns (read, in table).

    ⭐ **Levels 1-3 matter as much as level 4.** 126 of the 1,410 distinct EC codes in the two
    published catalogues state fewer than four fields, so a level-4-only load would leave exactly
    the vaguest codes unnamed. `2.7.10.-` is a row here, named *Protein-tyrosine kinases*.

    ⚠ Attribution is the whole licence condition (CC BY 4.0) and it is the **page's** job — this
    function only moves rows. See `EnzymeClass` and `nuna.tl.locus_browser.vendor_ec`.
    """
    from nuna.tl.locus_browser.vendor_ec import enzyme_reference

    rows = [{"ec_code": code, "name": name} for code, name in enzyme_reference().items()]
    _upsert(session, EnzymeClass, rows, key="ec_code", columns=("name",))
    in_table = session.execute(select(func.count()).select_from(EnzymeClass)).scalar_one()
    return len(rows), int(in_table)


def load_kegg_reference(session: Session) -> tuple[int, int]:
    """Upsert every KEGG KO name. Returns (read, in table).

    ⛔⛔ **Loading this is a licence decision, not a default.** KEGG requires an academic service
    provider licence from anyone offering services. `kegg_reference()` returns `{}` where the
    vendored table is absent, which is the correct state for any deployment that has not confirmed
    the subscription — the card then falls back to a bare linked accession, exactly as the published
    static pages do permanently.
    """
    from nuna.tl.locus_browser.vendor_kegg import kegg_reference

    rows = [
        {"ko_id": ko, "symbol": symbol, "definition": definition}
        for ko, (symbol, definition) in kegg_reference().items()
    ]
    _upsert(session, KeggOrthology, rows, key="ko_id", columns=("symbol", "definition"))
    in_table = session.execute(select(func.count()).select_from(KeggOrthology)).scalar_one()
    return len(rows), int(in_table)


def load_reference_vocabularies(session: Session) -> VocabularyReport:
    """Every global reference the catalogue joins against, in one idempotent pass.

    ⛔ **One entry point, so a new vocabulary is added in one place.** The two `--stage` call sites
    in `ingest_command_line` each used to name `load_pfam_reference` directly; a third vocabulary
    would then have had to be added twice, and the second site is the one that gets forgotten.
    """
    pfam = load_pfam_reference(session)
    ec_read, ec_in = load_enzyme_reference(session)
    kegg_read, kegg_in = load_kegg_reference(session)
    return VocabularyReport(
        pfam_families_read=pfam.pfam_families_read,
        pfam_families_in_table=pfam.pfam_families_in_table,
        unmapped_fields=pfam.unmapped_fields,
        enzyme_classes_read=ec_read,
        enzyme_classes_in_table=ec_in,
        kegg_orthologies_read=kegg_read,
        kegg_orthologies_in_table=kegg_in,
    )


def _upsert(session: Session, model: type, rows: list[dict], *, key: str, columns: tuple[str, ...]) -> None:
    """Batched `ON CONFLICT DO UPDATE` — the shape `load_pfam_reference` established.

    ⚠ Batched at 5,000 because a single `VALUES` of 28,516 KEGG rows exceeds what psycopg will
    bind comfortably, and the failure is a driver-level error rather than anything readable.
    """
    for start in range(0, len(rows), 5_000):
        statement = insert(model).values(rows[start : start + 5_000])
        session.execute(
            statement.on_conflict_do_update(
                index_elements=[key],
                set_={column: statement.excluded[column] for column in columns},
            )
        )
