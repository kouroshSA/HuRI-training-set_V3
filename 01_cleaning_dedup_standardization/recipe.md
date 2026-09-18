# Recipe — cleaning, deduplicating, standardizing the raw HuRI table

## Requirements

- Python 3.9+
- `pandas`

## Input

The raw HuRI CSV at `data/raw/huri_interactome_uniprotID_seqs.csv.xz` in this
repo. Decompress it first:

```bash
xz -dk data/raw/huri_interactome_uniprotID_seqs.csv.xz
# -> data/raw/huri_interactome_uniprotID_seqs.csv
```

## Step-by-step

Run from the repo root:

```bash
OUT=01_cleaning_dedup_standardization/output

# 1) Remove missing-sequence rows, strip stop codons, dedupe exact and
#    swapped-orientation pairs -> huri_interactome_uniprotID_seqs_cleaned_deduped.csv
python3 01_cleaning_dedup_standardization/scripts/01_clean_and_dedup.py \
    --input data/raw/huri_interactome_uniprotID_seqs.csv \
    --outdir "$OUT"

# 2) Scan for non-standard amino acid symbols, replace them, and write a
#    manifest of every change -> huri_interactome_uniprotID_seqs_cleaned_deduped_standard.csv
#    + huri_interactome_uniprotID_seqs_cleaned_deduped_standard_manifest.csv
python3 01_cleaning_dedup_standardization/scripts/02_standardize_amino_acids.py \
    --input "$OUT/huri_interactome_uniprotID_seqs_cleaned_deduped.csv" \
    --outdir "$OUT"
```

Each script prints a summary to stdout and also writes it to
`$OUT/logs/01_clean_and_dedup.log` / `02_standardize_amino_acids.log`.

## Output

| File | Description |
|---|---|
| `huri_interactome_uniprotID_seqs_cleaned_deduped.csv` | Output of script 1: missing-sequence rows removed, stop codons stripped, exact and swapped-orientation duplicate pairs removed. 50,176 rows. |
| `huri_interactome_uniprotID_seqs_cleaned_deduped_standard.csv` | Output of script 2: same rows as above, with every non-standard amino acid symbol replaced by a standard residue. 50,176 rows. **This is the positive pair table every later stage builds on.** |
| `huri_interactome_uniprotID_seqs_cleaned_deduped_standard_manifest.csv` | One row per residue substitution made in script 2 (row index, Protein_ID, sequence column, position, original/replacement residue, reason). 384 rows. |

None of the above `.csv` outputs are committed to this repo (see the repo
root README for why) — regenerate them with the commands above.

## Design choices worth knowing before reusing this on another table

- **"Missing sequence"** = blank field or the literal string `NAN`
  (case-insensitive), in either `Sequence_1` or `Sequence_2`. A row is
  dropped if *either* side is missing.
- **Deduplication** is ID-based, not sequence-based: two rows are considered
  the same interaction if they share the same pair of Protein_IDs, regardless
  of order. Step 3 catches same-order repeats; Step 4 additionally catches
  `(A,B)`/`(B,A)` swaps. Self-interactions (`ID_1 == ID_2`) are left alone —
  they are not duplicates.
- **Non-standard amino acid replacement** only ever touches `Sequence_1` /
  `Sequence_2`, never the ID columns. The replacement table (see
  `Methods.md`) uses a chemically closest standard residue when one clearly
  exists (`U → C`, `O → K`), and Alanine (`A`) as the ambiguous default for
  anything else (`B`, `Z`, `J`, `X`, or any unforeseen symbol).
- Both scripts are idempotent and safe to re-run — given the same
  `--input`/`--outdir`, they overwrite their own output files.

## Next stage

Scripts in `../02_homology_table/` build a single-protein homology table
from `huri_interactome_uniprotID_seqs_cleaned_deduped_standard.csv`, for
removing paralog/high-identity proteins before making train/validation/test
or in-/out-of-distribution splits.
