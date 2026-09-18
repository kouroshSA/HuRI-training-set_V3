# Recipe — building the single-protein homology table

## Requirements

- Python 3.9+, `pandas`
- `mmseqs2` (install once into a dedicated conda environment):

```bash
conda create -n homology -c bioconda -c conda-forge mmseqs2 -y
```

## Input

```
../01_cleaning_dedup_standardization/output/huri_interactome_uniprotID_seqs_cleaned_deduped_standard.csv
```

## Step-by-step

Run from the repo root:

```bash
POS=01_cleaning_dedup_standardization/output/huri_interactome_uniprotID_seqs_cleaned_deduped_standard.csv
OUT=02_homology_table/output

# 3) Pool + dedupe every (Protein_ID, Sequence) referenced anywhere in the
#    interaction table -> proteins/unique_proteins.fasta (+ .csv lookup)
python3 02_homology_table/scripts/03_extract_unique_proteins.py \
    --input "$POS" --outdir "$OUT"

# 4) All-vs-all similarity search with mmseqs2 -> mmseqs/all_vs_all.m8
#    (activates the "homology" conda env itself)
bash 02_homology_table/scripts/04_run_all_vs_all_search.sh \
    "$OUT/proteins/unique_proteins.fasta" "$OUT"

# 5) Collapse directed hits into one row per unordered protein pair,
#    with %identity and aligned-stretch length -> protein_homology_table.csv
python3 02_homology_table/scripts/05_build_homology_table.py --outdir "$OUT"
```

Each script prints a summary to stdout and also writes it to
`$OUT/logs/`. Step 4's mmseqs2 log is verbose (parameter dump); the
important lines are the final hit count and `Wrote raw hits to ...`.

## Output

| File | Description |
|---|---|
| `proteins/unique_proteins.fasta` | The 8,008 unique proteins referenced in the interaction table, as FASTA |
| `proteins/unique_proteins.csv` | Protein_ID, Length lookup for the same set |
| `mmseqs/all_vs_all.m8` | Raw mmseqs2 hits, BLAST-tab format, one row per detected directed alignment |
| `protein_homology_table.csv` | **The deliverable** — one row per unordered pair of proteins with detectable similarity, with percent identity, aligned-stretch length, protein lengths, coverage, gap/mismatch counts, E-value, bit score |

None of the above are committed to this repo — regenerate them with the
commands above.

## Notes for reuse / re-running

- No identity or coverage threshold has been applied to
  `protein_homology_table.csv` — it intentionally includes everything mmseqs2
  detected down to its `-e 10` E-value cutoff (percent identities from 14.2%
  up to 100% in this run), so a threshold can be chosen by inspecting the
  table rather than by re-running the search each time a different cutoff is
  tried.
- To actually remove paralog pairs from the interaction table once a
  threshold is chosen: pick a `Percent_Identity` (and optionally
  `Coverage_Shorter_Protein`) cutoff, take the set of `Protein_A`/`Protein_B`
  IDs in rows above that cutoff, and decide a policy (e.g. drop one member of
  each paralogous cluster from the training set; or ensure no two proteins
  above the cutoff end up on opposite sides of a train/test split). That
  policy step is deliberately not implemented here — this deliverable is the
  table needed to make that decision.
- Steps 3–5 are cheap to re-run (a few seconds to ~20 seconds total); if the
  underlying interaction table changes, just re-run all three in order.
- `Gap_Openings`/`Mismatches` require `-a 1` (backtrace) in step 4; removing
  that flag will silently make both columns 0 for every row (mmseqs2 still
  computes correct `Percent_Identity`/`Alignment_Length` without it, but not
  gap/mismatch counts).

## Next stage

`protein_homology_table.csv` (this deliverable) is used as-is — with no
threshold applied — by `../03_negative_pairs/scripts/06_build_negative_pairs.py`
to build a homology-safe, degree-matched set of in-distribution negative
pairs.
