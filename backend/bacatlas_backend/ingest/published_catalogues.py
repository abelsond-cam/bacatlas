"""The catalogues the live pages serve — the `(species, model_label, run_id)` triples, checked in.

⛔ **The run_id is NOT the model label with the species swapped, and guessing it fails silently in
the worst way: it resolves.** The two published runs are

    ecoli   …_preclusterstrict98pm3b-3b0.5-excl_…
    kp      …_preclusterkp98pm3b-3b0.5-excl_…

`strict98` against `kp98`. Substituting `kp` for `ecoli` in the *E. coli* run_id yields a path that
does not exist — which is the lucky case. A run_id that differs only in a **default-emitted token**
would resolve to a different run's assignment and load a whole catalogue that is wrong about which
model produced it, because *"run ids emit tokens non-default-only"*: two models differing only in a
ρ rule have produced byte-identical provenance.

So the triples live here, in the repo, rather than being reconstructed by a rule. `meta.model_id` in
a published payload is the only other place they are written down, and that is an output.

⚠ Both are **`nuna4` on the superseded DAMPED weight** — the `-excl` token, decoded by
`exclusivity.form_from_run_id`, is `damped_exclusion`. `-excl` is a *prefix* of `-exclLOGP`, so
never test it with `in` / `endswith` / `.replace`.

⛔ **`nuna5` is deliberately ABSENT, and that absence is the decision — not an omission to be fixed.**
This file names what the STATIC site serves, and the static site is `nuna4` only (David, 2026-09-25):
it is the parity oracle and the rollback, so it must not move when the service does. `nuna5` is
ingested into the database and reached through the service, which holds every catalogue and lets the
display pick one; it has no static page and no shipped payload for these triples to reconcile against.
Adding a row here would assert a publication that does not exist.

⚠ So this file and `pangenome.is_published` now mean **different things**, on purpose: `is_published`
is picker visibility inside the service (four catalogues are visible), while a row here is a page on
GitHub Pages (two). The day a nuna5 catalogue is published statically it gains a row and its own
`dset` — which is why that key is carried explicitly rather than derived from the species.

⛔ **It used to say "three keys, three jobs". It is FIVE** — the cohort axis added two, and the
table now lives in `published_cohorts.py`, which owns the vocabulary. The three this file is about:

    pathogen_species.species_key    ecoli          the organism
    PublishedCatalogue.dset         ecoli          the static page (`ecoli.html`), = `--dset`
    pangenome.catalogue_key         ecoli-nuna4    the service's address, what the picker sends

`dset` and `catalogue_key` differ for nuna4 and always will: the service reads the hyphen to tell
"pin this catalogue" from "follow this species' default", so it cannot spell nuna4 as `ecoli`; and
`ecoli.html` is a published url, so the static site cannot spell it `ecoli-nuna4`.

⚠ **`dset` is free to equal `catalogue_key` for any FUTURE catalogue and frozen for these two.**
GitHub Pages was deleted 2026-10-22, so nothing new will ever get a static page — but six test
modules read `data/{dset}.json` as the parity oracle, and `ecoli.html` / `kp.html` / `index.html`'s
redirect pin these two spellings.

⚠ **`PublishedCatalogue.model_label` is the ARTIFACT PATH COMPONENT** (`ecoli_nuna4_g2_0.98_…`),
consumed by `CatalogueArtifacts(model_label=…)`. `pangenome.model_label`, added 2026-10-27, is the
free-text string the footer prints. They hold the same characters today and are different facts —
do not make a third thing of that name.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PublishedCatalogue:
    """One species' live catalogue, as `published.tsv` and the shipped payload agree it is."""

    #: ⛔ The STATIC PUBLICATION key — the `published.tsv` row and the page it renders to, which is
    #: what nuna's exporter takes as `--dset`. It equals `species_key` for both rows here because
    #: nuna4 keeps the legacy urls `ecoli.html` / `kp.html`, and it must NOT be assumed to in general.
    #: ⚠ This is NOT `pangenome.catalogue_key`, which addresses the SERVICE and is `ecoli-nuna4`. The
    #: service reads the hyphen to tell "pin this catalogue" from "follow this species' default", so
    #: it cannot spell nuna4 as `ecoli`; the static site cannot spell it `ecoli-nuna4` without
    #: breaking a published url. Two vocabularies, both right, and they meet only here.
    dset: str
    species_key: str
    set_key: str
    model_label: str
    run_id: str
    #: What the shipped payload's `meta` says, so an ingest can reconcile against it rather than
    #: against a number someone remembered.
    genome_count: int
    gene_count: int
    locus_count: int


PUBLISHED_CATALOGUES: tuple[PublishedCatalogue, ...] = (
    PublishedCatalogue(
        dset="ecoli",
        species_key="ecoli",
        set_key="ecoli",
        model_label="ecoli_nuna4_g2_0.98_3b0.5rhoPAIRMAX_step4g0.1rhoCEIL",
        run_id="ecoli_bacformer_clever_exploded_preclusterstrict98pm3b-3b0.5-excl_k100_g100_res0.1_seed0",
        genome_count=100,
        gene_count=489_146,
        locus_count=17_531,
    ),
    PublishedCatalogue(
        dset="kp",
        species_key="kp",
        set_key="kp",
        model_label="kp_nuna4_g2_0.98_3b0.5rhoPAIRMAX_step4g0.1rhoCEIL",
        # ⛔ `preclusterkp98pm3b`, NOT `preclusterstrict98pm3b`. See the module docstring.
        run_id="kp_bacformer_clever_exploded_preclusterkp98pm3b-3b0.5-excl_k100_g100_res0.1_seed0",
        genome_count=100,
        gene_count=532_851,
        locus_count=15_670,
    ),
)


def catalogue(dset: str) -> PublishedCatalogue:
    """The statically published catalogue with this `published.tsv` key, or a failure naming what is known.

    ⛔ Keyed on `dset`, not `species_key`. Scanning on the species cannot express "this species has two
    published catalogues" — it returns whichever row comes first, silently — and that is precisely the
    shape the rebuild removed everywhere else. The two keys coincide for both rows TODAY; the lookup is
    written so that stops being load-bearing.
    """
    for entry in PUBLISHED_CATALOGUES:
        if entry.dset == dset:
            return entry
    known = [entry.dset for entry in PUBLISHED_CATALOGUES]
    raise KeyError(f"no published catalogue for {dset!r}; known: {known}")
