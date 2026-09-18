# ppi-v3-dataset-pipeline

Generalized, from-scratch pipeline for turning **any binary protein-protein
interaction (PPI) positive set** into a full V3-style training/evaluation
dataset: degree-matched real-protein negatives, three partner-randomized
decoy controls, ten Monte-Carlo cross-validation (MCCV) replicates with
matched holdouts, and two output encodings (ppiDCE/ppiBTEP and ppiGPLM).

This generalizes the dataset-construction conventions developed across the
MED4 and Synechocystis interactome projects (V2/V3 training sets, PRS/RRS
evaluation, partner-randomized decoy augmentation) into a reusable,
organism-agnostic tool, so a new positive interactome (any species, any
assay) can be turned into the same family of training artifacts without
re-deriving the recipe by hand each time.

## Input assumption

The only thing this pipeline needs to start is **a binary PPI positive
set**: a CSV with a header row containing (at minimum) two columns holding
the amino-acid sequences of each pair's two partners. Everything else
(negatives, decoys, replicates, formats, manifest) is derived from that.

## Pipeline

```
<your positive set>.csv  (supplied, any schema with 2 sequence columns)
        |  01_make_positives.py           --seq1-col / --seq2-col
        v
01_positives.csv                          canonical SEQ1,SEQ2,label=1 -- BOTH orientations
Parent_sequences/positives.csv            (byte-identical copy, lives with the other parent sets)
        |  02_make_random_pairs.py     (--mode in-distribution|out-of-distribution; seed 42, paralog-filtered)
        v
Parent_sequences/random_pairs.csv       real-protein negatives, both orientations (see "RRS mode" below)
        |  03_make_random_substitution_sets.py         (seeds 43 / 44 / 45)
        v
Parent_sequences/ps1_random.csv                        first partner replaced by random seq
Parent_sequences/ps2_random.csv                         second partner replaced by random seq
Parent_sequences/both_random.csv                        both partners replaced (the "floor")
        |  04_make_master_mix.py                       (seed 46)
        v
master_training_mix.csv                   all 5 sets concatenated + shuffled, 1 pos : 4 neg
        |  05_make_mccv_replicates.py       (seeds 1000s-6000s + k; final shuffle;
        v                                    orientation-depletion + RRS paralog checks)
MCCV/   10 replicates x {PRS, RRS, 3 random controls, depleted training set}
        |  06_convert_format_ppiDCE_ppiBTEP.py
        |  07_convert_format_ppiGPLM.py
        v
formats/ppiDCE_ppiBTEP/   formats/ppiGPLM/     same tree, two encodings
        |  08_make_manifest_and_provenance.py
        v
manifest.json + provenance.md
```

Run it all with one command:

```bash
cd scripts
python run_pipeline.py \
    --input /path/to/positive_set.csv \
    --seq1-col protein1_sequence --seq2-col protein2_sequence \
    --outdir /path/to/output_dir
```

Or, for an **out-of-distribution** random-pair (RRS) set instead of the
default in-distribution one:

```bash
python run_pipeline.py \
    --input /path/to/positive_set.csv \
    --seq1-col protein1_sequence --seq2-col protein2_sequence \
    --outdir /path/to/output_dir \
    --rrs-mode out-of-distribution \
    --proteome-csv /path/to/full_proteome.csv \
    --proteome-id-col orf_id --proteome-seq-col sequence \
    --proteome-filter-col in_scope_flag --proteome-filter-value True
```

Or run any step standalone -- every script takes `--outdir` and reads/writes
the fixed relative paths shown in the diagram above (see each script's
docstring/`--help` for its full flag list, including every seed).

## What each step does

**01 -- normalize.** Cleans sequences (strip whitespace, uppercase, drop one
trailing stop-codon `*`), drops rows with empty or non-standard-amino-acid
sequences, deduplicates by unordered pair, and reports measured diagnostics
(unique pairs, homodimers, distinct proteins) rather than assuming them.
Every non-homodimer pair is then written out in **both orientations** --
(A,B) and (B,A) -- so position (which column a protein lands in) never
correlates with label on its own; every downstream negative set inherits
this same structure automatically. Saves the canonical file at
`01_positives.csv` and an identical copy at `Parent_sequences/positives.csv`.

**02 -- random pairs (RRS), two modes.** The "ordinary" real-protein negative
set. Both modes write the same file (`Parent_sequences/random_pairs.csv`,
both orientations, paralog-filtered -- see below) but draw candidate
proteins from different pools, and answer different questions:

- **`--mode in-distribution`** (default). Real proteins from the **positive
  pool only**, paired at random via a **configuration model**
  (stub-shuffling) so the negative degree distribution matches the positive
  graph's. This is the point of the exercise -- without a degree match, a
  model can learn to read the label off how often it has seen a protein
  rather than off the pair. Every candidate protein has interacted with
  *something*; the question this set poses is "does it interact with this
  specific other one?" Reports the achieved Pearson r between positive and
  negative degree.
- **`--mode out-of-distribution`**. Real proteins drawn **uniformly at
  random from an external full-proteome pool** (`--proteome-csv`), not
  restricted to the positive pool -- most candidates here were never tested
  in the underlying screen at all. There's no positive-graph degree
  structure to match against an essentially-arbitrary external pool, so this
  is plain rejection-sampled uniform pairing rather than a configuration
  model (Pearson r is reported as N/A). This is the easier negative,
  answering "are two arbitrary proteome proteins a pair?" rather than "does
  this known interactor pair with this other one?" -- useful as a second,
  distinct difficulty tier alongside the in-distribution set. `--proteome-csv`
  needs `--proteome-id-col`/`--proteome-seq-col` (defaults: `orf_id`,
  `sequence`, matching `Yeast_proteome_reference_SGD.csv`'s schema) and
  optionally `--proteome-filter-col`/`--proteome-filter-value` to restrict
  the pool (e.g. to `in_paper_scope_2026_manuscript == True`).

Both modes: self-pairs, duplicate negative pairs, and accidental positive
pairs are rejected and re-drawn. A **paralog filter** additionally rejects
any candidate pair (X,Y) where X is a close sequence match to one partner of
some positive pair and Y is a close match to the other partner -- such a
pair could plausibly be a true, conserved interaction rather than a real
negative. Closeness is k-mer Jaccard similarity (`--paralog-k`,
`--paralog-threshold`; alignment-free, no external tools -- see "Paralog
filter" below for what this proxy does and doesn't guarantee). In
out-of-distribution mode the close-match map spans the positive pool *and*
the proteome pool, so a proteome protein that happens to be a near-duplicate
of a positive-pair partner is still caught. Output is written in both
orientations, matching step 1.

**03 -- random-substitution decoys.** Three "hard decoy" parent sets, one
substituted row per original positive row: `ps1_random` (first partner
replaced), `ps2_random` (second partner replaced), `both_random` (both
replaced -- the floor: no real protein content survives, only length and
composition). Each set uses its own seed for maximum diversity of the fake
proteins used. All random sequences are drawn uniformly from the 20
standard amino acids and are enforced **globally unique** -- never equal to
a real sequence or to any other generated sequence anywhere in the build.

**04 -- master training mix.** Concatenates the positive set and all four
parent negative sets, shuffles once. `both_random` is included in the bulk
mix (not held out for eval only), giving a 1:4 positive:negative ratio --
see "Design decisions" below.

**05 -- MCCV replicates.** Ten (default) independent Monte-Carlo
cross-validation replicates. Each replicate independently samples 10%
(default `--frac`) from **each of the five parent pools** (positives,
random pairs, and the three decoy sets) with its own seed,
writes those as that replicate's PRS / RRS / random-control evaluation
files, then removes every row matching one of those selected pairs
(orientation-independent, by canonical key -- so both the (A,B) and (B,A)
row are removed in one pass) from the master mix, and **re-shuffles what's
left** (its own seed, `--seed-base-shuffle` + k) to produce that replicate's
depleted training set. Three independent checks run automatically and the
script fails loudly if any of them find a problem:

1. **leak check** -- re-read the written depleted training set, confirm no
   held-out pair's canonical key survives.
2. **orientation-depletion check** -- separately confirm, for every held-out
   pair, that neither its (A,B) nor its (B,A) literal row survived (a second
   way of checking what (1) already guarantees by construction, rather than
   trusting the mechanism that did the removing).
3. **RRS paralog-collision check** -- re-derive the same paralog-forbidden
   set step 2 used and confirm this replicate's RRS holdout contains none of
   it (should be automatic, since RRS is subsampled from an
   already-filtered parent set -- checked anyway). **Pass the same
   `--mode`/`--proteome-csv` (and proteome column flags) to this script that
   step 2 was run with** -- otherwise the re-derived forbidden set is built
   over the wrong pool and this check isn't meaningful (`run_pipeline.py`
   handles this automatically via `--rrs-mode`).

MCCV, not k-fold: replicates draw independently, so their evaluation sets
may overlap between replicates. What's guaranteed is that each replicate's
training set is clean of *that replicate's own* held-out pairs.

**06 / 07 -- format conversion.** Same content, two encodings:

```
ppiDCE / ppiBTEP:  training  SEQ1,SEQ2,label        (3 col, label 0/1)
                    eval      SEQ1,SEQ2               (2 col, no label)

ppiGPLM (alt1):    training  $,SEQ1,!,SEQ2,<label>   (5 col, label <0>/<1>)
                    eval      $,SEQ1,!,SEQ2,<         (5 col, bare trailing marker --
                                                        the model completes the prompt)
```

Both are headerless, LF-line-ending only (see "Conventions").

**08 -- manifest and provenance.** Walks every canonical and converted
file, records row counts, sha256, and an LF-only check into `manifest.json`,
and writes a human-readable `provenance.md` with the pipeline diagram, the
full seed table, and a verification summary (leak checks, collision counts).
Run this last, after every other step has completed for a given `--outdir`.

## Directory layout produced under `--outdir`

```
01_positives.csv
01_positives_stats.json
Parent_sequences/
    positives.csv                       identical copy of 01_positives.csv
    random_pairs.csv                    real-protein negatives (in-distribution or out-of-distribution)
    ps1_random.csv
    ps2_random.csv
    both_random.csv
    *_stats.json
master_training_mix.csv
master_training_mix_stats.json
MCCV/
    training_sets/depleted_training_set-V3-{1..10}.csv
    PRS-RRS/PRS-V3-{1..10}.csv, RRS-V3-{1..10}.csv
    random_controls/ps1random-V3-{1..10}.csv, ps2random-V3-{1..10}.csv, bothrandom-V3-{1..10}.csv
    mccv_stats.json
formats/
    ppiDCE_ppiBTEP/   (same tree as above, converted)
    ppiGPLM/          (same tree as above, converted)
manifest.json
provenance.md
```

## Conventions (enforced by `scripts/common.py`, don't bypass)

- **LF line endings only, everywhere.** `csv.writer` defaults to CRLF; a
  stray `\r` is an out-of-vocabulary character for a char-level tokenizer
  and has broken a previous build in this dataset lineage. Always write
  through `common.write_rows_lf`.
- **Orientation-independent pair identity.** De-duplication and MCCV
  depletion always key on `frozenset({seq1, seq2})`, never on `(seq1, seq2)`
  order, so a pair can never leak into training via its swapped orientation.
- **Positives (and the in-distribution real-protein negatives) are stored in
  both orientations.** Every non-homodimer pair appears as both (A,B) and
  (B,A), so which column a protein lands in never correlates with label.
- **Global random-sequence uniqueness.** Every decoy sequence generated
  anywhere in a build is checked against every real sequence and every
  previously-generated decoy sequence. Collision counts are reported, not
  assumed to be zero.
- **Every seed is an explicit, documented CLI argument.** Nothing is left
  to an unseeded RNG. See the seed table in `provenance.md` after a run, or
  the table below for the pipeline defaults.

## Paralog filter -- what it is and isn't

`build_close_match_map` / `paralog_forbidden_keys` in `common.py` use k-mer
(default k=4) Jaccard similarity between sequences as a fast, dependency-free
proxy for "these two proteins are close enough that pairing one with a true
interactor's partner would risk fabricating a plausible real interaction as
a negative." This is **not** a real homology search -- no alignment, no
E-values, no gap handling. It is a cheap pre-filter, tuned via
`--paralog-k`/`--paralog-threshold`, that catches near-identical and
highly-repetitive sequences reasonably well and nothing more subtle. If you
have actual BLAST/HMMER results (or any real alignment-based homology calls)
for your protein set, prefer building `extra_forbidden` from those instead
of (or in addition to) this filter -- `configuration_model_negative_pairs`
in `common.py` takes `extra_forbidden` as a plain set of `canonical_key()`s,
so any external homology source can be wired in without touching the
sampling logic. Disable entirely with `--no-paralog-filter` on steps 02 and
05 (both must agree, or step 05's RRS paralog-verification will fail
against a forbidden set step 02 never applied).

## Seed table (defaults)

| stage | seed(s) |
|---|---|
| random pairs (either RRS mode) | 42 |
| ps1_random parent set | 43 |
| ps2_random parent set | 44 |
| both_random parent set | 45 |
| master mix shuffle | 46 |
| MCCV PRS selection, replicate k | 1000+k |
| MCCV RRS selection, replicate k | 2000+k |
| MCCV ps1-control selection, replicate k | 3000+k |
| MCCV ps2-control selection, replicate k | 4000+k |
| MCCV both-control selection, replicate k | 5000+k |
| MCCV final training-set shuffle, replicate k | 6000+k |

## Design decisions

These are defaults chosen so this generalized pipeline has *a* sane answer
for every question the source projects (MED4 V2/V3, Synechocystis V3) also
had to answer -- flag any of these if you want them changed:

- **RRS mode defaults to in-distribution.** Pass `--rrs-mode out-of-distribution
  --proteome-csv <path>` to draw the random-pair negative set from an
  external full-proteome pool instead of the positive pool. There is no
  "more correct" default -- in-distribution asks whether a known interactor
  pairs with *this* other protein; out-of-distribution asks whether two
  largely-untested proteome proteins pair at all. They're different
  difficulty tiers, not a strict/loose pair; pick (or build both and
  compare) based on what the downstream model needs to be evaluated against.
- **`both_random` is in the training mix**, following the Synechocystis
  build (1 pos : 4 neg) rather than the original MED4 V3 build, which used
  its equivalent for evaluation only (1:3). Drop `both_random.csv` from the
  sources list in `04_make_master_mix.py` and re-shuffle to restore 1:3.
- **MCCV replicates sample independently from all five parent pools**
  (positives + 4 negative types), rather than deriving the decoy controls
  by re-randomizing the selected PRS/RRS rows on the fly. This means a
  replicate's `ps1random`/`ps2random`/`bothrandom` control files are not
  guaranteed to correspond to the *same* underlying positive pairs as that
  replicate's PRS -- they are independent 10% draws from the same-sized
  parent pools. If you need row-for-row correspondence between a
  replicate's PRS and its decoy controls, sample by row index instead of
  by independent draw (a straightforward change to `05_make_mccv_replicates.py`).
- **Decoys derived from held-out positives are not depleted from training
  by default** (`--deplete-decoys-from-heldout` turns this on). A
  `ps1_random`/`ps2_random`/`both_random` row built from a positive pair
  that later gets held out as part of some replicate's PRS still appears in
  that replicate's training set, because its own canonical key (containing
  a fabricated sequence) doesn't match the held-out pair's key. This
  mirrors a known, explicitly-flagged gap in the Synechocystis V3 build;
  the opt-in flag closes it via row-alignment tracking through step 3.
- **MCCV fraction is 10% of unique positive pairs per replicate**
  (`--frac`), rather than a fixed absolute count (the source projects used
  a fixed 100 or 200). Pick whichever suits your dataset's scale --
  fixed-count evaluation sets are more comparable across replicates of
  different sizes; fraction-based sets scale with the input.

## Reproducing a run

Every script is deterministic given its seeds and inputs; re-running
reproduces byte-identical CSVs (verified in testing: two independent runs
of the full pipeline on the same input produced identical file content for
every `.csv` output). Check `manifest.json`'s `sha256` field per file to
confirm.

## Repository layout

```
scripts/
    common.py                              shared utilities (LF I/O, sequence
                                            cleaning, configuration-model
                                            negative generation, MCCV sampling,
                                            manifest helpers)
    01_make_positives.py
    02_make_random_pairs.py                (in-distribution or out-of-distribution RRS)
    03_make_random_substitution_sets.py
    04_make_master_mix.py
    05_make_mccv_replicates.py
    06_convert_format_ppiDCE_ppiBTEP.py
    07_convert_format_ppiGPLM.py
    08_make_manifest_and_provenance.py
    run_pipeline.py                        orchestrates 01-08
README.md
Claude-Code-instruction.md
LICENSE                                    PolyForm Noncommercial License 1.0.0
```

No third-party dependencies -- pure Python 3 standard library throughout.

## License

[PolyForm Noncommercial License 1.0.0](LICENSE) -- free for any noncommercial
purpose (research, education, personal use); commercial use requires a
separate license from the copyright holder.
