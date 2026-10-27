# `backend/` — the BacAtlas API

Postgres + Flask + SQLAlchemy. **Read-only for users**: no login, no accounts, no user writes.
Everything that writes is an offline ingest step.

Design of record: **`docs/design/serving_from_a_database.md`** in the `nuna` repo. Build order and
status: **`PROJECT_STATE.md` §6 and §Layer 6** there. Nothing in this directory carries a status
block.

## Bring it up on a Mac

```bash
brew install postgresql@16 && brew services start postgresql@16
export PATH="/opt/homebrew/opt/postgresql@16/bin:$PATH"
createdb bacatlas_ingest
psql -d bacatlas_ingest -c "CREATE EXTENSION IF NOT EXISTS pg_trgm;"

uv venv --python 3.13 .venv                    # from the repo root
uv pip install --python .venv/bin/python -e backend

export BACATLAS_DATABASE_URL="postgresql+psycopg://$USER@localhost:5432/bacatlas_ingest"
.venv/bin/python -m bacatlas_backend.serve --debug
curl -s localhost:5001/api/v1/health | python3 -m json.tool
```

`pg_trgm` is not optional — it is what preserves the page's **exact substring** search semantics
(`ligase` finds *O-antigen ligase RfaL* mid-string) rather than approximating them with prefixes.

## Configuration

Everything comes from the environment, resolved once at start-up, and a missing value is an
exception then rather than a 404 on one endpoint three weeks later.

| variable | meaning |
|---|---|
| `BACATLAS_DATABASE_URL` | required |
| `BACATLAS_PROFILE` | `development` (default) or `production` |
| `BACATLAS_SQL_ECHO` | `1` to log SQL |
| `BACATLAS_ROOT_GFF` | where the gzipped Bakta GFFs live — the **sequence source** |
| `BACATLAS_ROOT_ASSEMBLIES` | assembly FASTAs, the fallback if a GFF carries no `##FASTA` |
| `BACATLAS_ROOT_EMBEDDINGS` | ESM / Bacformer `.npy`; offline jobs only, never a request |
| `BACATLAS_ROOT_ANALYSIS` | audit tables, catalogue maps, assignment TSVs |

An `artifact_pointer` row stores a root **key** plus a relative path, so moving the store is a
config change and never a database migration.

## Layout

| path | what |
|---|---|
| `application_factory.py` | `create_application(configuration) -> Flask`; registers blueprints and decides nothing |
| `configuration.py` | the only module that reads `os.environ` |
| `database.py` | engine, session, declarative base — **and the constraint naming convention** |
| `api/` | one blueprint per resource |
| `services/` | the work; no Flask imports, so it is testable without a request context |
| `serialisers/` | dicts only — no ORM object escapes a serialiser |
| `models/` | SQLAlchemy models, one module per entity group |
| `ingest/` | the **only** writer. Offline. The one place that imports `nuna` |
| `gff/` | the one place a GFF is opened or parsed |
| `queries/` | hand-written SQL for hot paths, each with its `EXPLAIN` checked in |
| `migrations/` | alembic |

## Two rules that are load-bearing

**The serving install must not need `nuna`.** It is a private repo; the med school server has no
access to it and does not need any of it. Only `ingest/` imports it, and only on a machine that
already has the science package. That is why `nuna` is not in `dependencies`.

**The constraint naming convention in `database.py` had to exist before the first migration.**
Without it Postgres invents names, Alembic autogenerate produces a different one on each machine,
and a migration that drops a constraint by name works only for whoever wrote it.

## Loading data

Writes are ours alone and offline; nothing on a request path inserts. The loader needs the cluster
artifacts mirrored locally and `nuna` importable — the serving install needs neither.

```bash
uv pip install --python .venv/bin/python -e backend           # from the repo root
uv pip install --python .venv/bin/python -e ~/developer/nuna  # ingest only; never on the server
uv pip install --python .venv/bin/python scikit-learn         # nuna.eval.metrics imports it at module scope
```

⭐ **The pangenome layer imports `nuna` and the genome layer does not**, and the split is deliberate.
The genome layer needed one 20-line GFF reader, so it vendors it and cross-checks against the
original. The pangenome layer needs `node_order`, `oriented_windows`, `arrangements`,
`neighbour_counts`, `gap_table` and the vendored Pfam/GO/COG tables — the science itself — and a
second implementation of any of those is exactly what this design exists to avoid.

```bash
export BACATLAS_DATABASE_URL="postgresql+psycopg://$USER@localhost:5432/bacatlas_ingest"

# ONE genome layer for both species: contigs, genes, functional labels. Model-INDEPENDENT — a new
# pangenome must never rewrite it. ~180 s for all 280 genomes.
.venv/bin/python -m bacatlas_backend.ingest --stage genomes \
  --model-label ecoli_nuna4_g2_0.98_3b0.5rhoPAIRMAX_step4g0.1rhoCEIL \
  --run-id ecoli_bacformer_clever_exploded_preclusterstrict98pm3b-3b0.5-excl_k100_g100_res0.1_seed0

# then ONE pangenome layer per species: the model registry, the roster, the run and the catalogue.
# ~60 s each, re-runnable — the blast radius is exactly one pangenome.
.venv/bin/python -m bacatlas_backend.ingest --stage pangenome --set-key ecoli --species-key ecoli \
  --model-label ecoli_nuna4_g2_0.98_3b0.5rhoPAIRMAX_step4g0.1rhoCEIL \
  --run-id ecoli_bacformer_clever_exploded_preclusterstrict98pm3b-3b0.5-excl_k100_g100_res0.1_seed0

.venv/bin/python -m bacatlas_backend.ingest --stage pangenome --set-key kp --species-key kp \
  --model-label kp_nuna4_g2_0.98_3b0.5rhoPAIRMAX_step4g0.1rhoCEIL \
  --run-id kp_bacformer_clever_exploded_preclusterkp98pm3b-3b0.5-excl_k100_g100_res0.1_seed0
```

Then **publish** — a separate, explicit step, because the pointer is what the service reads and a
load that is not published changes nothing a reader can see. `--publish` runs every check in
`verify_pangenome_is_servable` and refuses to move the pointer if any fails: locus count against the
row count, the landing locus, both map representations, both **catalogue scatter sprites** and that
each accounts for every locus, every locus named, **every collection genome has its per-genome locus
count row**, and **that the cohort carries a label a reader can read** rather than its own key.
⚠ The checks are **named in the report**, and are deliberately not counted here — this sentence said
"twelve" while the function already ran thirteen, which is what a prose count does.

```bash
.venv/bin/python -m bacatlas_backend.ingest --stage pangenome --publish  ...  # as above
```

## The API

| endpoint | replaces |
|---|---|
| `GET /api/v1/species` | `published.tsv` |
| `GET /api/v1/species/{key}` | `meta`, the census, the footer provenance, `meta.landing`/`examples` |
| `GET /api/v1/species/{key}/loci/{label}` | ⭐ the hot path — `show(i)` minus the function tab |
| `GET …/loci/{label}/arrangements` | `focalCard`'s full scroller, paged |
| `GET …/loci/{label}/function` | `renderFunction`, on tab open |
| `GET /api/v1/species/{key}/search` | `search()` and the 3.0 MB `HAY` |
| `GET /api/v1/species/{key}/genomes` | `meta.genomes` + `anchorSearch` + `GENOME_N` — the anchor picker |
| `GET …/genomes/{bs}/loci/{label}/sequence` | ⭐ `withSeq` + `decodeNseq` + `genesAtLocus` + `translate` |
| `GET /api/v1/species/{key}/map/{rep}/scatter.png` | ⭐ `drawGlobal`'s catalogue dust — the one non-JSON response |

Measured on the loaded *E. coli* catalogue: the landing locus is **11.8 kB in 36 ms**, and a locus
view is **8 statements whether it resolves 2 neighbours or 49** — the property
`tests/test_api_endpoints.py` asserts, rather than a recorded budget that could be re-baselined.

⭐ **The anchor picker's counts are read, never aggregated per request.** Each genome carries TWO
numbers, and they are different facts: `locus_count` (loci where it has a gene) and
`arrangement_locus_count` (loci where it sits in an arrangement — exactly what the published dropdown
printed). A gene alone on its contig reaches no window, so the two differ for every probe genome
(ecoli: up to 286 loci, median 17.5; kp: up to 45, median 14.5). Both are `COUNT(DISTINCT locus)`,
written by the pangenome layer into `pangenome_genome_locus_count`; the aggregate behind them is ~412 M
membership rows at the design target, so a keystroke costs two statements instead.

⭐ **The scatter sprite is a picture because it is the one part of the map that is O(catalogue).**
Everything else the map needs is a handful of numbers; the dust behind the six dots is 889,160
positions at the design target, sent on every page load to draw a texture no reader reads a value
out of. Rendered once at ingest into `locus_map_scatter_sprite` and served with an **ETag that is
the content digest** — not the pangenome id like every JSON response here — because an image is
cached hard by caches we do not control, so a re-render has to be able to invalidate it.

⚠ **It does not win at today's size, and that is fine.** Measured: 17,531 loci is ~118 kB as a
sprite against ~70 kB of raw positions; the crossover is near 30k loci, and at 889,160 the positions
are 3.56 MB while the sprite is ~0.3 MB. A sprite is bounded by 1200² pixels; a coordinate array is
bounded by nothing.

⛔ **The sprite's four viewport columns are a contract, not metadata.** The published page could
project the six foreground loci itself because it held every coordinate; a client holding only the
picture cannot. Re-derive the transform anywhere and the focal dot lands *beside* its own speck
rather than on it — and the result still looks exactly like a scatter plot with a highlighted point.
`test_a_REAL_locus_projected_with_the_SERVED_viewport_lands_on_LIT_dust` checks 2,000 real loci
against the actual PNG bytes for precisely that.

⭐ **The parity suites need a loaded database.** `tests/test_catalogue_parity.py` (T1, T3a, T5, T7)
reads `BACATLAS_DATABASE_URL` and compares every locus against the published catalogue in
`data/{species}.json`. Without that variable, or with the run not loaded, it **skips with the
reason** rather than passing.

⛔ **Do not reconstruct a `run_id`.** The two published ones differ by `strict98` against `kp98`,
not by their species token, and a run id that differed only in a *default-emitted* token would
resolve to a different run and load a catalogue that is wrong about which model produced it. The
triples are checked in at `bacatlas_backend/ingest/published_catalogues.py`, with a test that both
resolve against the store.

⭐ **The Sequence tab reads the original GFF, and `.nseq` is retired rather than ported.** That
custom 2-bit-packed format existed for one reason: GitHub Pages applies `Range` to the *compressed*
stream, so byte offsets return plausible wrong bytes with a 206 and no error — the browser therefore
had to fetch whole files and decode DNA itself. A server has no such constraint. It slices the file
Bakta wrote (~5 kB out, where the page pulled 1.3 MB in), and the endpoint needs `BACATLAS_ROOT_GFF`.

⛔ **Flanks are in the GENE's reading direction, not the contig's.** On a minus-strand gene the
upstream flank sits at HIGHER contig coordinates and is reverse-complemented with it. Slicing
`start - flank` unconditionally returns the *downstream* flank for about half of all genes and
renders as a perfectly plausible 100 bases of DNA, so `test_gene_sequence_endpoint.py` goes back to
the contig and reads the bases itself.

⚠ **Three answers, and none of them may look like another:** `genes: []` with a 200 is *this genome
has no gene at this locus* (an answer); a 503 with a named reason is *the bases could not be read*;
a 404 is *no such genome or locus*.

| variable | meaning |
|---|---|
| `BACATLAS_ROOT_GFF` | the annotation store — **required by the sequence endpoint**, which is the only one that opens a file |
| `BACATLAS_NUNA_DATA_ROOT` | the local artifact mirror (default `~/developer/nuna/data`) |
| `BACATLAS_FULL_COHORT=1` | run the alignment gate over all 280 genomes (~85 s) instead of 30 |

The loader reports counts and then **asks the database what it holds**, because a loader that
reports its own writes back confirms itself. A refusal names the genome and the check that caught
it, rolls that genome back, and the run continues.

## Tests

```bash
.venv/bin/python -m pytest backend/tests -q
```

They need the `bacatlas_ingest` database running; override with `BACATLAS_TEST_DATABASE_URL`.

### ⚠ The `nuna` cross-check runs from nuna's venv, not this one

`gff/` carries **pinned copies** of two rules that live in `nuna` — `open_text`'s gzip-by-magic-number
and `genome_sequence.translate_cds` — because the serving install must not depend on a private repo.
A copy is only safe while something fails when it diverges, and that something is
`tests/test_gff_reader_against_nuna.py`. It is **not** run by the command above: `nuna/__init__.py`
pulls anndata, which pulls sklearn, and mirroring that chain into a serving package's test extra to
prove one 20-line function still matches would be the wrong trade. Run it where nuna already lives:

```bash
cd ~/developer/nuna
PYTHONPATH=~/developer/bacatlas/backend .venv/bin/python -m pytest \
  ~/developer/bacatlas/backend/tests/test_gff_reader_against_nuna.py -q
```

It skips — saying which is missing — if `nuna` is not importable or the probe GFFs have not been
pulled to `~/developer/nuna/data/raw/gff`. **A skip here is not a pass**; if the cross-check cannot
run, the vendored copies are unverified.

## Linting — `ruff check` is enforced, `ruff format` is NOT

```bash
uvx ruff check backend/bacatlas_backend backend/tests    # must pass; run before every commit
```

⚠ **Do not run `ruff format` over the tree.** Both are configured in `backend/pyproject.toml`, but only
`check` has ever been applied: **54 of 88 files would be reformatted**, so a formatting pass produces a
diff that is almost entirely unrelated to whatever you were doing and buries it. Format the lines you
touch by hand, to the configured `line-length = 120`.

⚠ **Fix lint findings by hand rather than with `--fix`.** This is the sibling repo's rule and it applies
here for the same reason: `--fix` mutates the tree during a commit, and it once deleted three re-exports
that were pinned by an identity test, turning ten tests red in a commit whose hooks all passed.
