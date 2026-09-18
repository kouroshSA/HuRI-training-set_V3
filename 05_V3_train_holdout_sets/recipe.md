# Recipe — building the HuRI V3 train/holdout sets

## Requirements

- Python 3.9+ (pure standard library — the vendored pipeline has no
  third-party dependencies)
- bash

## Input

```
../04_master_assembly/output/huri_interactome_uniprotID_seqs_cleaned_deduped_standard.csv
```

(Only the positive table is needed here — the pipeline builds its own RRS
and decoy negatives from it; see `Methods.md` for why the earlier
`huri_negative_pairs_indistribution.csv` isn't used as an input to this
stage.)

## Step-by-step

Run from the repo root:

```bash
POS=04_master_assembly/output/huri_interactome_uniprotID_seqs_cleaned_deduped_standard.csv
OUT=05_V3_train_holdout_sets/output/build

# 1) Run the vendored ppi-v3-dataset-pipeline end to end (steps 01-08):
#    positives -> RRS -> 3 random-substitution decoys -> master mix ->
#    10 MCCV replicates (PRS/RRS/decoy holdouts + depleted training sets) ->
#    2 format conversions -> manifest + provenance.
#    Takes several minutes: steps 02 and 05 each run an O(n^2) k-mer
#    paralog-similarity scan over ~8,000 proteins in pure Python.
bash 05_V3_train_holdout_sets/scripts/01_run_v3_pipeline.sh \
    "$POS" "$OUT" 05_V3_train_holdout_sets/pipeline_src

# 2) Reorganize each MCCV replicate's 6 files into one folder per replicate
#    (HuRI-V3-1 .. HuRI-V3-10).
python3 05_V3_train_holdout_sets/scripts/09_organize_per_replicate_folders.py \
    --build-dir "$OUT" \
    --outdir 05_V3_train_holdout_sets/output
```

## Output

```
05_V3_train_holdout_sets/output/
├── build/                                     canonical pipeline output (see pipeline_src/README.md)
│   ├── 01_positives.csv, 01_positives_stats.json
│   ├── Parent_sequences/{positives,random_pairs,ps1_random,ps2_random,both_random}.csv (+ *_stats.json)
│   ├── master_training_mix.csv, master_training_mix_stats.json
│   ├── MCCV/training_sets/depleted_training_set-V3-{1..10}.csv
│   ├── MCCV/PRS-RRS/{PRS,RRS}-V3-{1..10}.csv
│   ├── MCCV/random_controls/{ps1random,ps2random,bothrandom}-V3-{1..10}.csv
│   ├── MCCV/mccv_stats.json
│   ├── formats/ppiDCE_ppiBTEP/   formats/ppiGPLM/     (same tree, 2 alternate encodings)
│   └── manifest.json, provenance.md
└── HuRI-V3-1/ .. HuRI-V3-10/                   one self-contained folder per MCCV replicate
    ├── training_set.csv
    ├── PRS.csv
    ├── RRS.csv
    ├── ps1_random_control.csv
    ├── ps2_random_control.csv
    ├── both_random_control.csv
    └── replicate_summary.json
```

None of this is committed to this repo (it is all derived, regeneratable
data/training material) — see the repo root README.

## Verifying a run

```bash
OUT=05_V3_train_holdout_sets/output/build

# no CR anywhere
find "$OUT" -name '*.csv' | xargs grep -lU $'\r'   # expect no output

# leak / orientation-leak / RRS-paralog-collision totals (should all be 0)
python3 -c "
import json
d = json.load(open('$OUT/MCCV/mccv_stats.json'))
r = d['replicates']
print('leak:', sum(x['leak_check_pairs_present_in_depleted_set'] for x in r))
print('orient leak:', sum(x['orientation_leak_check']['forward_A-B'] + x['orientation_leak_check']['reverse_B-A'] for x in r))
print('rrs paralog collisions:', sum(x['rrs_paralog_collisions'] for x in r))
"

# label sanity per replicate folder (PRS all 1, everything else all 0)
for f in 05_V3_train_holdout_sets/output/HuRI-V3-1/*.csv; do
  echo -n "$f: "; awk -F',' '{c[$NF]++} END {for (k in c) printf "%s=%d ", k, c[k]; print ""}' "$f"
done
```

## Reusing the vendored pipeline directly

Every step is also runnable standalone; see `pipeline_src/README.md` and
`pipeline_src/Claude-Code-instruction.md` for the full flag reference
(every seed is an explicit, overridable CLI argument). To pick a different
RRS mode (out-of-distribution, against an external proteome pool), a
different MCCV fraction, or a different replicate count, edit
`scripts/01_run_v3_pipeline.sh` or call `pipeline_src/scripts/run_pipeline.py`
directly with different flags.

## Notes for reuse

- **Re-running is deterministic**: every seed is fixed (see the seed table
  in `Methods.md` / `pipeline_src/README.md`), so re-running
  `01_run_v3_pipeline.sh` against the same input reproduces byte-identical
  `.csv` output (verify via the `sha256` field in `manifest.json`).
- **`pipeline_src/` is a pinned snapshot, not a live clone** — it will not
  pick up upstream changes. To rebuild against a newer version of
  `ppi-v3-dataset-pipeline`, re-clone it, re-copy `scripts/` (+ docs) over
  `pipeline_src/`, update `pipeline_src/VENDORED_COMMIT.txt`, and re-run.
- **Script 09 depends only on `<outdir>/MCCV/`** and can be re-run any time
  after step 1 without redoing the (slow) pipeline run, e.g. if the
  per-replicate folder layout itself needs to change.
