# Recipe — assembling the master/train staging folder

## Requirements

- bash, `sha256sum` (standard on Linux/macOS)

## Inputs

```
../01_cleaning_dedup_standardization/output/huri_interactome_uniprotID_seqs_cleaned_deduped_standard.csv
../03_negative_pairs/output/huri_negative_pairs_indistribution.csv
```

## Step-by-step

Run from the repo root:

```bash
bash 04_master_assembly/scripts/07_assemble_master_and_train.sh \
    01_cleaning_dedup_standardization/output/huri_interactome_uniprotID_seqs_cleaned_deduped_standard.csv \
    03_negative_pairs/output/huri_negative_pairs_indistribution.csv \
    04_master_assembly/output
```

Copies both tables into `04_master_assembly/output/`, verifies each copy is
byte-identical to its source (SHA-256 + row count), and writes the
verification output to `04_master_assembly/output/logs/07_assemble_master_and_train.log`.
Fails loudly (non-zero exit) if any checksum or row count doesn't match.

## Output

```
04_master_assembly/output/
├── huri_interactome_uniprotID_seqs_cleaned_deduped_standard.csv   (positive pairs, 50,176 rows)
├── huri_negative_pairs_indistribution.csv                         (negative pairs, 50,176 rows)
└── logs/07_assemble_master_and_train.log
```

The two pair tables are **not merged** — no label column, no shuffling, no
split. Not committed to this repo — regenerate with the command above. See
`Methods.md` for why, and for what the next stage does with them.

## Next stage

`../05_V3_train_holdout_sets/` turns these two tables into the full V3
training/holdout dataset (RRS, random-substitution decoys, labeled and
shuffled MCCV replicates).
