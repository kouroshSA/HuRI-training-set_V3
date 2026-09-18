# HuRI-traing-set_V3

Code, methods, and recipes for turning the raw HuRI human binary
protein-protein interaction table into a full V3-style training/holdout
dataset: cleaning and standardization, a single-protein homology table,
homology-safe in-distribution negative pairs, and 10 Monte-Carlo
cross-validation (MCCV) replicates each with its own held-out PRS, RRS, and
three partner-randomized decoy controls.

**This repo holds the starting data plus every script/Methods/recipe needed
to reproduce the full workflow — it does not hold the final train/holdout
datasets themselves.** Those are large, fully regeneratable from the code
here, and are excluded on purpose (see "What's not here" below).

## Workflow

```
data/raw/huri_interactome_uniprotID_seqs.csv.xz     the starting HuRI set (included, compressed)
        |
        v
01_cleaning_dedup_standardization/    drop missing/duplicate pairs, standardize
        |                             non-standard amino acids -> the positive
        |                             interaction table every later stage uses
        v
02_homology_table/                    mmseqs2 all-vs-all search over every
        |                             unique protein -> a %identity / aligned-
        |                             length table for any pair with detectable
        |                             similarity
        v
03_negative_pairs/                    degree-matched, homology-safe in-
        |                             distribution negative (non-interacting)
        |                             pairs, same count as the positives
        v
04_master_assembly/                   stage positives + negatives together,
        |                             kept as two separate files (not merged
        |                             or labeled yet)
        v
05_V3_train_holdout_sets/             the vendored ppi-v3-dataset-pipeline:
                                       RRS, 3 random-substitution decoys,
                                       label 1/0, shuffle, 10 MCCV replicates
                                       each with PRS/RRS/decoy holdouts ->
                                       HuRI-V3-1 .. HuRI-V3-10
```

Each stage folder has its own `Methods.md` (what was done and why, with
exact counts from the reference run), `recipe.md` (exact commands to
reproduce it), and `scripts/` (the code). Read them in order — each one
says what the next stage needs from it.

## What's here

- `data/raw/` — the starting HuRI interactome table (compressed).
- `01_cleaning_dedup_standardization/` … `05_V3_train_holdout_sets/` —
  scripts, Methods, and recipes for every stage.
- `05_V3_train_holdout_sets/pipeline_src/` — a pinned snapshot of the
  separately-maintained
  [`ppi-v3-dataset-pipeline`](https://github.com/kouroshSA/ppi-v3-dataset-pipeline)
  repo (see `REFERENCES.md`), used unmodified for stage 5.
- `REFERENCES.md` — citation for the underlying HuRI dataset and the
  vendored pipeline.

## What's not here

The actual generated datasets — cleaned/standardized positives, the
homology table, negative pairs, the master training mix, and the final
`HuRI-V3-1` .. `HuRI-V3-10` train/holdout folders — are **not** committed.
They are large (the full V3 build is on the order of tens of GB across all
10 replicates and format conversions) and fully reproducible from the code
and starting data in this repo. Run each stage's `recipe.md` in order to
regenerate them locally.

## Requirements

- Python 3.9+, `pandas`
- `mmseqs2` for stage 2 (`conda create -n homology -c bioconda -c conda-forge mmseqs2 -y`)
- bash, `xz`, `sha256sum` (standard on Linux/macOS)
- No third-party dependencies for the vendored pipeline in stage 5 (pure
  Python standard library)

## Portability note

Every script here takes explicit `--input`/`--outdir`-style arguments (or
positional equivalents for the bash scripts) — there are no
machine-specific absolute paths anywhere in this repo. Each stage's
`recipe.md` gives the exact commands assuming you run them from the repo
root, writing outputs into an `output/` folder inside that stage (ignored by
git — see `.gitignore`).

## License / reuse

Private repository. The vendored pipeline in
`05_V3_train_holdout_sets/pipeline_src/` carries its own license (PolyForm
Noncommercial 1.0.0 — see that folder's `LICENSE`); everything else in this
repo is unlicensed by default (all rights reserved) unless/until a license
is added.
