# Claude Code instructions for this repo

This repo is a **generalized dataset-construction pipeline**, not a
one-off analysis. It exists so that turning a new binary PPI positive set
into a V3-style training/eval dataset never again has to be re-derived by
hand -- it should always mean "run these 8 scripts (or `run_pipeline.py`)
against a new input file with a new `--outdir`."

Read `README.md` first -- it documents the pipeline stages, directory
layout, seed table, and the explicit "Design decisions" this generalized
version made where the source projects (MED4, Synechocystis) didn't have a
single settled answer.

## Ground rules specific to this repo

- **Outputs never go in this repo.** `--outdir` should always point outside
  the repository (e.g. under `~/Dropbox/<project>/`). This repo holds only
  the pipeline code and docs. Never commit a generated dataset here, even a
  small one for "example" purposes -- if an example is needed, generate it
  on demand from the synthetic-data pattern below rather than checking in
  CSVs.
- **LF line endings are load-bearing, not a style preference.** Every CSV
  write must go through `common.write_rows_lf` / `common.open_lf`. If you
  add a new script that writes CSV, do not use a bare `csv.writer` on a
  normally-opened file -- `csv.writer` defaults to CRLF, and a stray `\r`
  is an out-of-vocabulary character for the char-level tokenizers these
  datasets feed. This exact bug has broken a build in this dataset lineage
  before; don't reintroduce it.
- **Pair identity is always orientation-independent.** Use
  `common.canonical_key(a, b)` (a `frozenset`), never a raw tuple, for any
  de-duplication, membership check, or depletion logic. A pair that leaks
  through its swapped orientation is a silent, hard-to-detect bug.
- **Every source of randomness is a named, defaulted, overridable CLI
  seed.** No bare `random.random()` / unseeded `random.Random()` anywhere.
  If you add a new randomized step, give it a seed flag, document the
  default in both the script docstring and the README seed table, and pick
  a default that doesn't collide with the existing seed ranges (42-46 for
  the parent-set stage, 1000+k/2000+k/.../5000+k for MCCV replicate k).
- **Global uniqueness of generated (fake) sequences is enforced, not
  assumed.** Route new decoy-generation code through `common.UniqueSeqPool`
  (or extend it) rather than calling `random.choice` in a loop directly --
  it's what tracks collisions against real sequences and against every
  other generated sequence in the same build.
- **Verify by re-reading the file you just wrote, not by trusting the
  in-memory object that wrote it.** Every existing verification step (e.g.
  the MCCV leak check in `05_make_mccv_replicates.py`) opens the written
  file fresh and checks it independently of the generator's internal
  state. Follow this pattern for any new check you add -- it's what has
  caught real bugs in this dataset lineage (e.g. a nondeterministic
  configuration-model build caused by iterating an unsorted Python `set`,
  found only because a rerun was diffed against the first run).
- **Positives (and anything built row-for-row from them) carry both
  orientations.** `01_make_positives.py` writes every non-homodimer pair as
  both (A,B) and (B,A); `03_make_random_substitution_sets.py` inherits this
  for free by iterating the positives file as-is. If you add a new step
  that builds a negative set from scratch (not by substitution on an
  existing orientation-symmetric file), it must also emit both orientations
  -- an orientation-asymmetric negative set next to an orientation-symmetric
  positive set turns "which column is this sequence in" into a usable
  shortcut for the label.
- **Step 2 and step 5 must be run with matching `--mode`/`--proteome-*`
  flags.** Step 5 re-derives step 2's paralog-forbidden set from scratch
  (see "Verify by re-reading" above) rather than trusting inheritance --
  that only works if it's told the same pool step 2 used. `run_pipeline.py`
  handles this via one shared `--rrs-mode` flag; if you add a third mode or
  another mode-dependent step, keep that single-source-of-truth pattern
  rather than letting each script default independently.
- **The paralog filter is a cheap proxy, not a homology search.** See the
  "Paralog filter" section of README.md before changing
  `build_close_match_map`/`paralog_forbidden_keys` in `common.py` -- it's
  k-mer Jaccard similarity, deliberately alignment-free so the pipeline
  stays dependency-free. If a change needs real homology calls, wire them
  in via `configuration_model_negative_pairs`'s `extra_forbidden` parameter
  rather than trying to make the k-mer proxy more sophisticated in place.

## When adding a new pipeline step

1. New scripts go in `scripts/`, numbered to reflect where they sit in the
   pipeline (or as an unnumbered variant of an existing numbered step, e.g.
   `05b_...py`, if it's an alternative rather than a replacement).
2. Import shared logic from `common.py`; do not duplicate CSV I/O,
   sequence-cleaning, or pair-identity logic in a new script.
3. Every script takes `--outdir` and reads/writes the fixed relative paths
   documented in the README's directory-layout section. If your script
   introduces a new canonical file that should be format-converted (steps
   06/07) or manifested (step 08), add it to `common.iter_canonical_files`
   and to the `gather_files` glob list in `08_make_manifest_and_provenance.py`
   -- both are the single source of truth for "what counts as pipeline
   output," and adding a file elsewhere silently excludes it from
   conversion and the manifest.
4. Update `README.md`'s pipeline diagram, directory layout, and (if
   applicable) seed table and design-decisions section in the same change.
   Docs drifting from the actual scripts is exactly the failure mode this
   repo exists to avoid.
5. Test on a small synthetic dataset before running on anything real (see
   below) -- it's cheap and catches most bugs (degree-matching failures,
   off-by-one row counts, format mismatches) in seconds instead of minutes.

## Testing changes

There is no fixed example dataset in the repo (see "outputs never go in
this repo" above). To smoke-test a change, generate a tiny synthetic
positive set on the fly and run the full pipeline against it:

```python
import random, csv
rng = random.Random(7)
AA = 'ACDEFGHIKLMNPQRSTVWY'
proteins = {f'P{i}': ''.join(rng.choice(AA) for _ in range(rng.randint(20, 60))) for i in range(30)}
names = list(proteins.keys())
pairs, rows = set(), []
while len(rows) < 80:
    a, b = rng.sample(names, 2)
    key = frozenset((a, b))
    if key in pairs:
        continue
    pairs.add(key)
    rows.append((a, proteins[a], b, proteins[b]))
with open('/tmp/toy_positives.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['idA', 'seqA', 'idB', 'seqB'])
    w.writerows(rows)
```

Then:

```bash
cd scripts
python run_pipeline.py --input /tmp/toy_positives.csv \
    --seq1-col seqA --seq2-col seqB \
    --outdir /tmp/toy_output --n-replicates 3 --frac 0.10
```

Check: the run completes with `total leak-check failures = 0`, `total
orientation-leak failures = 0`, and `total RRS paralog collisions = 0`;
`find /tmp/toy_output -name '*.csv' | xargs grep -lU $'\r'` returns nothing
(no CR anywhere); and re-running into a second `--outdir` produces
byte-identical `.csv` files (only the `*_stats.json` files should differ,
and only in the embedded output-path strings). Also spot-check that
`Parent_sequences/positives.csv` is byte-identical to `01_positives.csv`,
and that `01_positives.csv` has (very close to) 2x the unique-pair count in
rows, since every non-homodimer pair should appear twice.

## What "provenance" means in this project

Every real (non-toy) run should end with step 08 run against it, producing
`manifest.json` (machine-readable: per-file row counts, sha256, LF check)
and `provenance.md` (human-readable: pipeline diagram with actual counts,
full seed table, verification summary) in the output directory --
alongside the data, not in this repo. This is the standing convention
across the whole MED4/Synechocystis dataset lineage: a dataset without a
provenance record next to it is treated as suspect until one is generated.
