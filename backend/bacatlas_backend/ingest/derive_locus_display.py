r"""What a locus is CALLED and what it is DESCRIBED as — `app.js`'s rules, moved to ingest.

⭐ **This is the fan-out fix.** `displayOf` falls through the modal Bakta symbol → the locus's *whole*
Pfam list plus the Pfam name table → its *whole* product list. Rendering one locus touches a median
15–19 other loci and up to 303, because a neighbour's NAME is that transitive read. Computed once
into a NOT NULL column, the fan-out becomes `WHERE locus_id = ANY($1)` over ~20 ids.

⛔ **Never let an inferred name reach a graded statistic.** `app.js` says it in as many words:
*"NEVER let this reach n_named or any graded statistic — it is a display fallback, not an
annotation."* `named_member_count` and `bakta_gene_symbol` come from the assignment's own symbols and
are untouched by anything here.

⚠ **Three regex fidelity notes, because JS and Python `re` are not the same engine.**

1. **`re.ASCII` on every `\\b` pattern.** JavaScript's `\\b` is ASCII-only (`\\w` is `[A-Za-z0-9_]`);
   Python's is Unicode-aware by default. A Bakta product carrying an accented character would put a
   word boundary in a different place, and the swap guard would admit a candidate the page rejects.
2. **`fullmatch` where the JS pattern is anchored `^…$`.** `re.match` alone would accept a prefix.
3. **`re.escape` for the symbol**, which is a superset of the JS escape class and never a subset.

⚠ **The share denominator is the LISTED products, not the locus.** `bestProduct` sums `p[j].n` over
the entries the payload ships — capped at `TOP_PRODUCTS = 6` — so a locus with a long tail has a
denominator smaller than its member count and its 10 % floor is correspondingly easier to clear.
That is the published behaviour; computing the share against `member_gene_count` would silently
raise the bar and change what the card says.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

#: `app.js::GENE_SHAPE` — the shape a Pfam short name's last token must have to be read as a symbol.
GENE_SHAPE = re.compile(r"[A-Z][a-z]{2}[A-Z0-9]?", re.ASCII)

#: `app.js::VAGUE`. A candidate matching this *"says less than almost anything, so it may never WIN a
#: swap"*. ⛔ Deliberately NOT applied to the modal product — replacing a vague modal with something
#: specific is the entire point of the swap.
VAGUE = re.compile(
    r"^(hypothetical|uncharacteri[sz]ed|putative|conserved|predicted)?\s*"
    r"(protein|membrane protein|lipoprotein|exported protein|secreted protein)"
    r"(\s+[A-Za-z]{3}[A-Za-z0-9]*)?$",
    re.IGNORECASE | re.ASCII,
)

#: `app.js::inferredName`'s product pattern — a symbol named right before `family` or `-like`.
PRODUCT_SYMBOL = re.compile(r"\b([A-Z][a-z]{2}[A-Z0-9]?)\b(?=\s+family\b|-like\b)", re.ASCII)

#: The minimum share of the LISTED product counts a candidate must hold to be allowed to win.
PRODUCT_SWAP_MINIMUM_SHARE = 0.1

#: `locus.display_name_source`'s vocabulary — so the page can say WHY it calls a locus what it does.
SOURCE_BAKTA_SYMBOL = "bakta_symbol"
SOURCE_PFAM_ARCHITECTURE = "pfam_architecture"
SOURCE_PRODUCT = "product"
SOURCE_LABEL = "label"
#: A symbol Bakta gave some member genes, chosen over the commonest one because other annotation
#: names it — see `corroborated_symbol`. `source_accession` holds the Pfam or COG id that did.
SOURCE_CORROBORATED_SYMBOL = "corroborated_symbol"

#: Below this share of a locus's genes, its commonest Bakta symbol is too thin a vote to name the
#: locus on its own. ⭐ Set by census over all four catalogues in `bacatlas_ingest` (2026-11-18):
#: anywhere from 2 % to 20 % the rule renames exactly kp 722 (`gumC` → `wzc`) and ecoli 5168
#: (`orfA` → `tnp`); at 50 % it renames 5 kp and 11 ecoli loci.
WEAK_SYMBOL_SHARE = 0.2

#: A gene name as it appears INSIDE annotation text — `wzc`/`gumC` (symbol case) or `Wzc`/`GumC`/`Wzc1`
#: (protein-name case). All-caps tokens (`ABC`, `HTH`, `AAA`) are deliberately not names.
NAME_IN_TEXT = re.compile(r"[a-z]{3}[A-Z0-9]?|[A-Z][a-z]{2}[A-Z0-9]?", re.ASCII)
WORD = re.compile(r"[A-Za-z0-9]+", re.ASCII)


@dataclass(frozen=True)
class DisplayName:
    """What to call a locus, and on what evidence."""

    name: str
    source: str
    source_accession: str | None = None


def sole_architecture(pfam_entries: list[tuple[str, int]], pfam_names: dict) -> str | None:
    """`app.js::soleArch` — the commonest architecture's accession, but only if it has exactly one.

    ⚠ *Sole* refers to the architecture having a single accession, not to the locus having a single
    architecture. A two-domain architecture names no gene, so it is not a source for a symbol.
    """
    if not pfam_entries:
        return None
    # The list is already rank-ordered by count descending, and `app.js` keeps the first on a tie
    # (its loop advances only on a strict `>`), so rank 0 IS the answer.
    label = str(pfam_entries[0][0])
    accessions = label.split(",")
    if len(accessions) != 1:
        return None
    return accessions[0] if accessions[0] in pfam_names else None


def infer_name(
    pfam_entries: list[tuple[str, int]],
    product_entries: list[tuple[str, int]],
    pfam_names: dict,
) -> DisplayName | None:
    """`app.js::inferredName` — a gene name for a locus Bakta never named, or `None`.

    Two levels, in order. **Level one needs `pfam_names`**, which is attached at RENDER time and is
    NOT in the exported payload: skipping it silently demotes 370 *E. coli* loci to `locus <label>`,
    which looks like a data gap rather than a missing lookup table.
    """
    accession = sole_architecture(pfam_entries, pfam_names)
    if accession:
        short_name = str(pfam_names[accession][0])
        token = short_name.split("_")[-1]
        if GENE_SHAPE.fullmatch(token):
            return DisplayName(token, SOURCE_PFAM_ARCHITECTURE, accession)
    if product_entries:
        match = PRODUCT_SYMBOL.search(str(product_entries[0][0]))
        if match:
            return DisplayName(match.group(1), SOURCE_PRODUCT, None)
    return None


def names_in_text(text: str | None) -> frozenset[str]:
    """The gene names a piece of annotation text mentions, casefolded — `Wzc1` also yields `wzc`."""
    keys: set[str] = set()
    for token in WORD.findall(text or ""):
        if NAME_IN_TEXT.fullmatch(token):
            key = token.casefold()
            keys.add(key)
            if key[-1].isdigit():
                keys.add(key[:-1])
    return frozenset(keys)


def corroborated_symbol(
    symbol_entries: list[tuple[str, int]],
    member_gene_count: int,
    pfam_entries: list[tuple[str, int]],
    product_entries: list[tuple[str, int]],
    cog_entries: list[tuple[str, int]],
    pfam_names: dict,
    cog_names: dict,
) -> DisplayName | None:
    """Where the symbol vote is too thin to decide, the Bakta symbol other annotation also names.

    ⭐ **The `gumC` case** (kp 722): 100 member genes, two named — `gumC` ×1, `wzc` ×1 — and the
    tie went to the alphabet. Every one of the 100 carries Pfam PF23607 *"Tyrosine-protein kinase
    Wzc, N-terminal domain"*, so `wzc` is the name the locus's own annotation supports.

    Fires only when the commonest symbol is on fewer than `WEAK_SYMBOL_SHARE` of the member genes.
    Each of the locus's listed symbols is then scored by the most genes any ONE source names it on
    — its own Bakta count, the architectures with a domain whose Pfam short name or description
    names it, its COG names, its products — and the highest wins; ties go to the Bakta count, then
    the alphabet. Returns `None` when the winner is the commonest symbol anyway, so every locus
    this does not rename keeps `SOURCE_BAKTA_SYMBOL` exactly as before.

    ⛔ **Candidates are only symbols Bakta gave this locus's genes — never a name read out of the
    other sources.** Pfam, COG and product text name a FAMILY by its founding member (LysR, DeoR,
    Wzz, Wzc). Letting them propose a name, or override a decisive vote, was measured on the same
    catalogues: `yhaJ` 95 → `lysR` 5, `glpR` → `deoR`, `etk` → `wzc`, and 600+ loci named `the`
    and `and`. They may only choose between names the genes themselves carry.
    """
    if not symbol_entries or member_gene_count <= 0:
        return None
    modal, modal_count = str(symbol_entries[0][0]), int(symbol_entries[0][1])
    if modal_count >= WEAK_SYMBOL_SHARE * member_gene_count:
        return None

    # ⚠ Per source, so one gene is never counted twice for one name; the score is the max.
    bakta: dict[str, int] = {}
    for symbol, count in symbol_entries:
        bakta[str(symbol).casefold()] = bakta.get(str(symbol).casefold(), 0) + int(count)
    pfam: dict[str, int] = {}
    pfam_by_accession: dict[str, dict[str, int]] = {}
    for label, count in pfam_entries:
        named_here: set[str] = set()
        for accession in set(str(label).split(",")):
            fields = pfam_names.get(accession)
            if not fields:
                continue
            keys = names_in_text(str(fields[0]).replace("_", " ")) | names_in_text(str(fields[1]))
            named_here |= keys
            for key in keys:
                by_accession = pfam_by_accession.setdefault(key, {})
                by_accession[accession] = by_accession.get(accession, 0) + int(count)
        for key in named_here:
            pfam[key] = pfam.get(key, 0) + int(count)
    cog: dict[str, int] = {}
    cog_by_id: dict[str, dict[str, int]] = {}
    for cog_id, count in cog_entries:
        for key in names_in_text(cog_names.get(str(cog_id))):
            cog[key] = cog.get(key, 0) + int(count)
            by_id = cog_by_id.setdefault(key, {})
            by_id[str(cog_id)] = by_id.get(str(cog_id), 0) + int(count)
    product: dict[str, int] = {}
    for text, count in product_entries:
        for key in names_in_text(str(text)):
            product[key] = product.get(key, 0) + int(count)

    def support(key: str) -> int:
        return max(source.get(key, 0) for source in (bakta, pfam, cog, product))

    own_count = {str(symbol): int(count) for symbol, count in reversed(symbol_entries)}
    winner = min(own_count, key=lambda s: (-support(s.casefold()), -own_count[s], s))
    if winner == modal:
        return None
    key = winner.casefold()
    # Name the accession that carried the most genes for the winner, so the page can say WHICH.
    evidence = None
    if pfam.get(key, 0) >= max(cog.get(key, 0), product.get(key, 0)) and key in pfam_by_accession:
        evidence = min(pfam_by_accession[key].items(), key=lambda item: (-item[1], item[0]))[0]
    elif cog.get(key, 0) >= product.get(key, 0) and key in cog_by_id:
        evidence = min(cog_by_id[key].items(), key=lambda item: (-item[1], item[0]))[0]
    return DisplayName(winner, SOURCE_CORROBORATED_SYMBOL, evidence)


def display_name(
    node_label: str,
    modal_symbol: str | None,
    pfam_entries: list[tuple[str, int]],
    product_entries: list[tuple[str, int]],
    pfam_names: dict,
    *,
    symbol_entries: list[tuple[str, int]],
    member_gene_count: int,
    cog_entries: list[tuple[str, int]],
    cog_names: dict,
) -> DisplayName:
    """`app.js::displayOf` — symbol, else inferred, else `locus <label>`. Never NULL.

    ⭐ One departure from `app.js`, by decision (David, 2026-11-18): a thin symbol vote is settled by
    `corroborated_symbol` before the commonest symbol is taken. ⚠ The keyword arguments are required
    so no caller can run the old rule by leaving them out.
    """
    if modal_symbol and symbol_entries and str(symbol_entries[0][0]) != str(modal_symbol):
        raise ValueError(
            f"locus {node_label}: modal symbol {modal_symbol!r} is not the first listed symbol "
            f"{symbol_entries[0][0]!r} — the two must come from one vote"
        )
    corroborated = corroborated_symbol(
        symbol_entries, member_gene_count, pfam_entries, product_entries, cog_entries,
        pfam_names, cog_names,
    )
    if corroborated is not None:
        return corroborated
    if modal_symbol:
        return DisplayName(str(modal_symbol), SOURCE_BAKTA_SYMBOL, None)
    inferred = infer_name(pfam_entries, product_entries, pfam_names)
    if inferred is not None:
        return inferred
    return DisplayName(f"locus {node_label}", SOURCE_LABEL, None)


def best_product(
    modal_symbol: str | None, product_entries: list[tuple[str, int]]
) -> str | None:
    """`app.js::bestProduct` — the most SPECIFIC product, not the commonest.

    ⛔ **Not `products[0]`.** At `rfaL` the modal is *"Ligase"* while the same locus also holds
    *"O-antigen ligase RfaL"*. A candidate may win only if it **extends** the modal string or
    **names this locus's own symbol** on a word boundary the modal does not — either way it cannot
    become a different protein, which plain "longest wins" did (`araB` matched inside *"L-arabinose
    isomerase"*).
    """
    if not product_entries:
        return None
    modal = str(product_entries[0][0])
    total = sum(int(count) for _, count in product_entries)
    symbol_pattern = (
        re.compile(rf"\b{re.escape(str(modal_symbol))}\b", re.IGNORECASE | re.ASCII)
        if modal_symbol
        else None
    )
    modal_lower = modal.lower()
    best, best_extends = None, False
    for value, count in product_entries:
        candidate = str(value)
        if candidate.lower() == modal_lower:
            continue  # case-only churn
        if total and count / total < PRODUCT_SWAP_MINIMUM_SHARE:
            continue
        if VAGUE.fullmatch(candidate.strip()):
            continue  # never swap INTO a vague product
        extends = modal_lower in candidate.lower()
        names = bool(
            symbol_pattern
            and symbol_pattern.search(candidate)
            and not symbol_pattern.search(modal)
        )
        if not extends and not names:
            continue
        if (
            best is None
            or (extends and not best_extends)
            or (extends == best_extends and len(candidate) > len(best))
        ):
            best, best_extends = candidate, extends
    return best or modal


def search_text(
    node_label: str,
    modal_symbol: str | None,
    product_entries: list[tuple[str, int]],
    uniref50_accessions: list[str],
) -> str:
    """`app.js::HAY` — symbol, label, every listed product and every listed UniRef50, lowercased.

    ⛔ **Pfam, COG, GO, EC and KEGG are deliberately absent**, because that is what the page searches
    today. Adding them is a product change and not an implementation detail: a query for `PF00593`
    would start matching, which is either an improvement or a surprise, and it is not this port's
    call to make.

    ⚠ The API must escape `%` and `_` before interpolating a query into `ILIKE '%'||q||'%'`. They
    are LIKE wildcards and plain literals to `indexOf`, so the two are **not** equivalent for a
    query containing either — and product strings contain both.
    """
    parts = [
        str(modal_symbol or ""),
        str(node_label),
        " ".join(str(value) for value, _ in product_entries),
        " ".join(str(accession) for accession in uniref50_accessions),
    ]
    return " ".join(parts).lower()
