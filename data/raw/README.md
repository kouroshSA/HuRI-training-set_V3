# Raw starting data

`huri_interactome_uniprotID_seqs.csv.xz` — the starting HuRI interactome
table this whole workflow is built from. See `../../REFERENCES.md` for the
citation and where the underlying interaction data comes from.

- Columns: `Protein_ID_1, Sequence_1, Protein_ID_2, Sequence_2`
- Uncompressed: 171,543 rows, 141 MB, SHA-256
  `ff858265c5c54af4f1ec64f97e554364055b5672272b2062e30cdbb1447573a5`
- Compressed with `xz -9` (3.5 MB) rather than gzip (79 MB) — the file has
  massive long-range sequence redundancy (most of the ~8,000 distinct
  proteins recur across many of the 171,543 rows) that only a large-window
  compressor like LZMA captures.

## Decompress before use

```bash
xz -dk huri_interactome_uniprotID_seqs.csv.xz
# -> huri_interactome_uniprotID_seqs.csv

sha256sum huri_interactome_uniprotID_seqs.csv
# should print ff858265c5c54af4f1ec64f97e554364055b5672272b2062e30cdbb1447573a5
```

This is the `--input` for `../../01_cleaning_dedup_standardization/scripts/01_clean_and_dedup.py`
— see that stage's `recipe.md` to continue the workflow.
