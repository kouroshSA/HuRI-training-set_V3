# Recipe — building the in-distribution negative pair set

## Requirements

- Python 3.9+, `pandas` (no mmseqs2/conda env needed for this step — it only
  reads `protein_homology_table.csv`, already computed by stage 2)

## Inputs

```
../01_cleaning_dedup_standardization/output/huri_interactome_uniprotID_seqs_cleaned_deduped_standard.csv
../02_homology_table/output/protein_homology_table.csv
../02_homology_table/output/proteins/unique_proteins.fasta
```

## Step-by-step

Run from the repo root:

```bash
python3 03_negative_pairs/scripts/06_build_negative_pairs.py \
    --positives 01_cleaning_dedup_standardization/output/huri_interactome_uniprotID_seqs_cleaned_deduped_standard.csv \
    --homology  02_homology_table/output/protein_homology_table.csv \
    --fasta     02_homology_table/output/proteins/unique_proteins.fasta \
    --outdir    03_negative_pairs/output
```

Prints a round-by-round convergence log (phase 1) followed by the phase-2
top-up summary and a final verification pass, and writes the same summary to
`03_negative_pairs/output/logs/06_build_negative_pairs.log`.

## Output

| File | Description |
|---|---|
| `huri_negative_pairs_indistribution.csv` | **The deliverable** — 50,176 negative (non-interacting) pairs, same schema as the positive table (`Protein_ID_1, Sequence_1, Protein_ID_2, Sequence_2`), drawn from the same 8,008 proteins, degree-matched, and free of any positive-pair, duplicate, or homology-rule violation. |

Not committed to this repo — regenerate with the command above.

## Design choices worth knowing before reusing this

- **"Significant identity" = any row in `protein_homology_table.csv`.** No
  extra numeric %identity/coverage cutoff is layered on top. If you want a
  looser rule (e.g. only reject above 30% identity), filter
  `protein_homology_table.csv` down to that subset before running this
  script — the script does not take a threshold argument, it treats the
  entire table as "significant" by design (the strictest option the data
  supports).
- **Self-interacting proteins are excluded from the negative pool** — no
  `X-X` negative pairs are generated. Their positive self-loops still count
  toward that protein's target usage (a self-loop contributes 2 to a
  protein's stub count), so they still influence how often that protein is
  paired with *other* proteins in the negative set.
- **Exact degree-matching is a soft target, not a hard guarantee.** The
  script does its best via phase 1 (strict) and phase 2 (top-up against the
  full pool), and logs the final deviation per protein. In this run, 75.4%
  of proteins matched exactly; the rest deviate by a handful of pairs at
  most (see `Methods.md` for the exact numbers). If a stricter guarantee is
  ever needed, the alternative discussed but not used here is to remove the
  highest-identity/highest-degree positive pairs first, shrinking the
  hardest cases before sampling — deliberately not the default, since it
  touches the positive table.
- Uses a fixed random seed (`RANDOM_SEED = 20260918` in the script) for
  reproducibility; change it to get a different draw.

## Next stage

`huri_negative_pairs_indistribution.csv` (this deliverable) is staged
alongside the positive pair table by `../04_master_assembly/`.
