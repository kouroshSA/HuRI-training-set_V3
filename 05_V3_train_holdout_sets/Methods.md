# Methods — HuRI V3 Training/Holdout Set Construction

## Purpose

Turn the standardized HuRI positive interactome
(`huri_interactome_uniprotID_seqs_cleaned_deduped_standard.csv`, 50,176
positive pairs / 8,008 proteins, from `../01_cleaning_dedup_standardization/`)
into a full V3-style training/evaluation dataset: a real-protein random
negative set (RRS), three partner-randomized decoy controls, and 10
Monte-Carlo cross-validation (MCCV) replicates each with its own held-out
Positive Reference Set (PRS), RRS, and decoy controls — using the existing,
previously-built
[`ppi-v3-dataset-pipeline`](https://github.com/kouroshSA/ppi-v3-dataset-pipeline)
rather than re-deriving this from scratch.

## Reviewing the pipeline before use

The repository was cloned fresh and read in full (`README.md`,
`Claude-Code-instruction.md`, `scripts/common.py`, and every numbered
script) to confirm current behavior before running it. As of this build,
the repository was at:

```
commit d41e5083ad233be011beef95d32990f75b231bd1
date:   2026-08-21 13:15:17 +0400
subject: Add out-of-distribution RRS mode alongside in-distribution
```

Its history shows two prior commits (`f1101ca` initial generalized
pipeline, `3aef188` orientation-duplicate positives / paralog filter / final
shuffle / license) — the "out-of-distribution RRS mode" in the latest commit
is the modification most likely to be new relative to any earlier use of
this repo.

**A pinned snapshot of the pipeline code (this exact commit) is vendored
into `pipeline_src/`** in this folder, so this dataset build stays
reproducible even if the upstream repository changes later. It is a plain
file copy, not a git checkout — see `pipeline_src/VENDORED_COMMIT.txt`.

## What the pipeline does (as reviewed)

Full detail is in `pipeline_src/README.md`; summary of the terminology used
throughout this build:

- **PRS (Positive Reference Set)**: a held-out sample of true interacting
  pairs.
- **RRS (Random Reference Set)**: a held-out sample of real-protein pairs
  that are *not* known to interact, built by degree-matched random pairing
  within the positive protein pool (**in-distribution** mode — see below),
  paralog-filtered so a pair can't be a plausible-but-untested true
  interaction via conserved paralogous binding.
- **Random-substitution decoys** (`ps1_random`, `ps2_random`,
  `both_random`): one partner (or both) of every positive pair replaced with
  a same-length, uniformly random amino-acid sequence, so no real protein
  content survives on the substituted side.
- **MCCV (Monte-Carlo cross-validation)**: 10 independent replicates, each
  sampling its own 10% holdout from every one of the five pools (positives,
  RRS, and the three decoy sets), then removing every held-out pair from the
  master training mix and re-shuffling what's left.

## Command run

```bash
bash scripts/01_run_v3_pipeline.sh \
    <positives.csv> <outdir> pipeline_src
```

which invokes:

```bash
python3 pipeline_src/scripts/run_pipeline.py \
    --input <positives.csv> \
    --seq1-col Sequence_1 --seq2-col Sequence_2 \
    --outdir <outdir> \
    --n-replicates 10 --frac 0.10 --version-tag V3 \
    --rrs-mode in-distribution
```

(All other flags — seeds, paralog-filter k/threshold, `both_random`
inclusion in the training mix — left at the pipeline's documented defaults;
see the seed table in `pipeline_src/README.md`.)

**RRS mode: in-distribution** (the pipeline's default) was used — real
proteins drawn from the same 8,008-protein HuRI pool, paired via a
degree-matched configuration model. This asks "does this known interactor
pair with this *other* specific protein?", the harder and more informative
negative for this dataset, matching how the separately-built
`huri_negative_pairs_indistribution.csv` (stage `../03_negative_pairs/`) was
also framed as in-distribution.

### Relationship to the earlier custom negative-pair set

Stage `../03_negative_pairs/` already built a rigorous in-distribution
negative set using a real mmseqs2 alignment-based homology table
(`protein_homology_table.csv`) rather than the pipeline's lighter-weight
k-mer Jaccard paralog proxy. That earlier set was **not** substituted in as
this build's RRS. The pipeline's own `Parent_sequences/random_pairs.csv` was
generated fresh via its own `02_make_random_pairs.py`, using its own
paralog filter, because:

1. The intent here was to use the pipeline as reviewed, producing its
   standard PRS/RRS/decoy artifacts end to end.
2. `common.py`'s `configuration_model_negative_pairs` exposes an
   `extra_forbidden` parameter specifically for wiring in an external,
   alignment-based homology source instead of (or alongside) the built-in
   k-mer filter, but no CLI flag currently loads one from a file — using it
   would mean extending the pipeline's CLI surface, which was out of scope
   for this build.

Wiring the alignment-based `protein_homology_table.csv` into
`02_make_random_pairs.py`/`05_make_mccv_replicates.py` via `extra_forbidden`
remains a documented, available future enhancement if a stricter paralog
guarantee than the k-mer proxy is ever needed for this dataset.

## An expected data-quality finding: 538 additional duplicate pairs

Step 1 of the pipeline (`01_make_positives.py`) reported:

```
dropped_duplicate_unordered_pair: 538
n_unique_unordered_pairs: 49638   (input had 50,176 rows / 49,706 non-self ID pairs)
n_homodimers: 466
```

This is expected, not a bug: stage `../01_cleaning_dedup_standardization/`
dedupes by **Protein_ID pair**, while this pipeline dedupes by
**sequence-content pair** (`canonical_key(seq1, seq2)`, i.e. the two
amino-acid strings themselves, since the pipeline only ever sees sequences,
not IDs — it was run with `--seq1-col Sequence_1 --seq2-col Sequence_2` and
never sees the ID columns). 538 pairs of distinct Protein_ID pairs turned
out to share identical sequence content on at least one side (e.g. an
isoform ID and its base accession mapping to the same string, or two
accessions for what is functionally the same sequence) and collapsed under
the pipeline's sequence-level identity — a stricter, independent uniqueness
check that stage 1's ID-based dedup could not have caught. `n_homodimers`
(466) is correspondingly slightly lower than stage 1's earlier count of 470
self-pairs by ID, for the same reason.

## Label convention

Inherited unchanged from the pipeline, and unambiguous throughout every
output file (`SEQ1,SEQ2,label`, no header):

- **`1`** = a true interacting pair (positives / PRS)
- **`0`** = every negative type: RRS (real-protein random pairs) **and**
  all three random-substitution decoys (`ps1_random`, `ps2_random`,
  `both_random`)

## Results

| Stage | Count |
|---|---|
| Input positive rows | 50,176 |
| Unique sequence-level positive pairs | 49,638 (466 homodimers) |
| `01_positives.csv` rows (both orientations) | 98,810 |
| RRS (`random_pairs.csv`) unique pairs / rows | 49,638 / 99,276 |
| RRS degree-match Pearson r (vs. positive graph) | **0.9965** |
| Paralog-shadow pairs excluded from RRS candidates | 54,379 (out of 7,991 distinct proteins checked; 244 proteins had ≥1 close match) |
| `ps1_random` / `ps2_random` / `both_random` rows | 98,810 each |
| Random decoy sequences generated | 395,240 (0 collisions with real or other generated sequences) |
| `master_training_mix.csv` rows | 494,516 (98,810 pos : 395,706 neg = **1 : 4.00**) |
| MCCV replicates | 10, each a 10% holdout of every one of the 5 parent pools |
| Per-replicate PRS / RRS size | 4,964 pairs each |
| Per-replicate decoy-control size (each of the 3) | 9,881 pairs each |
| Per-replicate depleted training set size | ~445,055–445,072 rows (varies slightly per replicate's holdout draw) |

### Verification (from `<outdir>/provenance.md` / `manifest.json` in your own run)

- **199 files** produced and checked; **199/199 pass the LF-only
  line-ending check** (no `\r` anywhere).
- **MCCV leak check** (a held-out pair re-appearing in that replicate's own
  depleted training set), summed over all 10 replicates: **0**.
- **Orientation-depletion check** (both the `(A,B)` and `(B,A)` row of a
  held-out pair independently confirmed absent from the depleted set),
  summed over all 10 replicates: **0**.
- **RRS paralog-collision check** (a replicate's held-out RRS pair matching
  the paralog-forbidden set), summed over all 10 replicates: **0**.
- Every `.csv` output re-read from disk and its SHA-256 recorded in
  `manifest.json`.

## Per-replicate folder reorganization

The pipeline's native output groups files by *type* across all replicates
(`MCCV/training_sets/depleted_training_set-V3-{1..10}.csv`,
`MCCV/PRS-RRS/PRS-V3-{1..10}.csv`, etc.). An additional reorganization step
(`scripts/09_organize_per_replicate_folders.py`) copies each replicate's
six files into one self-contained folder per replicate — `HuRI-V3-1/` ..
`HuRI-V3-10/` — each holding:

```
training_set.csv           depleted, shuffled training mix for this replicate (label 1/0, 1:4 ratio)
PRS.csv                    held-out true positives (label 1)
RRS.csv                    held-out real-protein random pairs (label 0)
ps1_random_control.csv     held-out: partner 1 randomized, partner 2 native (label 0)
ps2_random_control.csv     held-out: partner 1 native, partner 2 randomized (label 0)
both_random_control.csv    held-out: both partners randomized (label 0)
replicate_summary.json     this replicate's row counts and seeds (from mccv_stats.json)
```

Every copy is verified byte-identical (SHA-256) against its source at copy
time. The canonical `<outdir>/MCCV/` tree is left intact alongside these —
this is an additional view, not a replacement.

## Bonus artifacts from running the full pipeline

Since `run_pipeline.py` runs steps 1–8 end to end, this build also produces
`formats/ppiDCE_ppiBTEP/` and `formats/ppiGPLM/` — the same full tree
converted to two alternative headerless encodings (3-column
`SEQ1,SEQ2,label` and a 5-column prompt-style format respectively; see
`pipeline_src/README.md` "06 / 07 — format conversion"), plus
`manifest.json` and `provenance.md` with full per-file row counts, SHA-256
hashes, and the verification summary quoted above.

**None of these data outputs (positives, RRS, decoys, master mix, MCCV
replicates, per-replicate folders, format conversions) are committed to
this repo** — see the repo root README for why. Only the code that produces
them (this folder) is tracked.

See `recipe.md` for exact reproduction commands.
