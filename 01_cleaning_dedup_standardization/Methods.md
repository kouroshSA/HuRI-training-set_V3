# Methods — HuRI Interactome Cleaning, Deduplication, Standardization

## Source data

`data/raw/huri_interactome_uniprotID_seqs.csv` at the repo root (see
`../REFERENCES.md` for where the underlying HuRI interaction data comes
from), a table of human binary protein-protein interactions with
UniProt/Ensembl identifiers and amino acid sequences for both interaction
partners.

- Columns: `Protein_ID_1, Sequence_1, Protein_ID_2, Sequence_2`
- Rows (interactions): 171,543
- Sequences contain no characters outside the 20 standard amino acids plus
  `U` (Selenocysteine, 3 occurrences) and `X` (unknown residue, 1,241
  occurrences) across the whole file; no `*` stop-codon symbols were present
  anywhere in the raw file.
- A missing sequence is encoded as the literal string `NAN`; 5,008 rows had
  `NAN` in `Sequence_1` (or an outright blank field) and 3,943 in `Sequence_2`.
- The raw table turned out to contain heavy redundancy: only 52,929 of the
  171,543 rows represent distinct undirected interaction pairs once orientation
  (A,B) vs (B,A) is ignored — most pairs are listed 3 times, consistent with
  HuRI being a union of three separate yeast two-hybrid screening stages that
  was not deduplicated before export.

## Step 1 — Remove rows with a missing sequence

A sequence field was treated as missing if, after trimming whitespace, it was
an empty string or equal to `NAN` (case-insensitive). Any row missing either
`Sequence_1` or `Sequence_2` was dropped, since a usable training pair needs
both partner sequences.

- **8,903 rows removed** (171,543 → 162,640).

## Step 2 — Remove stop codons

Any `*` stop-codon symbol occurring anywhere in a sequence was stripped.

- **0 sequences affected** — no `*` characters were present in this dataset.
  The step was still run (and is left in the pipeline for reuse on other
  interactome tables where it may not be a no-op).

## Step 3 — Deduplicate repeated pairs (same orientation)

Rows with an identical `(Protein_ID_1, Protein_ID_2)` pair in the same order
were collapsed, keeping the first occurrence. (All duplicate rows for a given
ordered ID pair carried identical sequences, so no information was lost by
keeping only the first.)

- **86,730 rows removed** (162,640 → 75,910).

## Step 4 — Deduplicate swapped-orientation pairs

Because HuRI interactions are undirected, a pair listed as `(A, B)` in one row
and `(B, A)` in another row represents the same interaction. Pairs were
canonicalized by sorting the two IDs alphabetically and dropping rows whose
canonical pair had already been seen, keeping the first occurrence
encountered in the file.

- **25,734 rows removed** (75,910 → 50,176).
- 470 self-interaction rows (`Protein_ID_1 == Protein_ID_2`, i.e. homodimers)
  remained and were **not** removed — they are not duplicates, only their
  reported self-pairing.

## Step 5 — Save cleaned/deduplicated table

Result: `huri_interactome_uniprotID_seqs_cleaned_deduped.csv` —
**50,176 rows**, same 4 columns as the source file, no missing sequences, no
stop-codon symbols, no exact or swapped-orientation duplicate pairs.

## Step 6 — Standardize non-standard amino acid symbols

The two sequence columns (never the ID columns) of the cleaned/deduplicated
table were scanned character-by-character against the 20 standard amino acid
letters (`ACDEFGHIKLMNPQRSTVWY`). Non-standard symbols found were replaced
according to the rule requested: an unambiguous, chemically-closest standard
residue where one exists, otherwise Alanine (`A`) as the ambiguous default.

| Symbol | Meaning | Replacement | Rationale |
|---|---|---|---|
| `U` | Selenocysteine | `C` (Cysteine) | Direct structural analogue of cysteine (Se → S) |
| `O` | Pyrrolysine | `K` (Lysine) | Direct structural analogue of lysine |
| `B` | Asx (Asp *or* Asn) | `A` (Alanine) | Ambiguous between two residues — no single closest one |
| `Z` | Glx (Glu *or* Gln) | `A` (Alanine) | Ambiguous between two residues — no single closest one |
| `J` | Xle (Leu *or* Ile) | `A` (Alanine) | Ambiguous between two residues — no single closest one |
| `X` | Unknown residue | `A` (Alanine) | Unknown identity |
| any other symbol | — | `A` (Alanine) | Fallback for any unforeseen non-standard character |

Only `U` and `X` were actually encountered in the cleaned/deduplicated table:

- **384 total residue substitutions** across 384 individual sequence
  positions (each substitution is in a distinct sequence/row, i.e. no sequence
  needed more than one fix):
  - `U → C`: 1 occurrence
  - `X → A`: 383 occurrences
- 28 unique Protein_IDs were affected.

Every substitution is logged in
`huri_interactome_uniprotID_seqs_cleaned_deduped_standard_manifest.csv`, one
row per residue change, with the row index (into the standardized/cleaned
file), the affected Protein_ID, which sequence column, the 1-based position
in the sequence, the original and replacement residue, and the reason.

Result: `huri_interactome_uniprotID_seqs_cleaned_deduped_standard.csv` — same
50,176 rows, sequences now restricted to the 20 standard amino acid letters
only. **This is the positive interaction table used by every later stage.**

## Row-count summary

| Stage | Rows |
|---|---|
| Raw HuRI table | 171,543 |
| After removing missing-sequence rows (Step 1) | 162,640 |
| After stop-codon stripping (Step 2) | 162,640 (no change) |
| After same-orientation dedup (Step 3) | 75,910 |
| After swapped-orientation dedup (Step 4) | 50,176 |
| Final cleaned/deduped file (Step 5) | 50,176 |
| Final standardized file (Step 6) | 50,176 (384 residues changed) |

## Verification performed

- Confirmed 0 residues outside the 20 standard amino acids remain in the
  final standardized file.
- Confirmed 0 rows with a missing/blank sequence, 0 rows containing `*`, and
  0 duplicate or swapped-orientation pairs remain in the final files.
- Spot-checked manifest entries against both the pre- and
  post-standardization sequences to confirm the logged position and
  substitution match the actual file contents.

## Files produced (not committed to this repo — see repo root README)

- `huri_interactome_uniprotID_seqs_cleaned_deduped.csv`
- `huri_interactome_uniprotID_seqs_cleaned_deduped_standard.csv`
- `huri_interactome_uniprotID_seqs_cleaned_deduped_standard_manifest.csv`

## Next stage

`huri_interactome_uniprotID_seqs_cleaned_deduped_standard.csv` (this stage's
output) is used by `../02_homology_table/` to build a single-protein
homology table, for removing paralog/high-identity proteins before making
train/validation/test or in-/out-of-distribution splits.

See `recipe.md` for the exact commands and `scripts/` for the code.
