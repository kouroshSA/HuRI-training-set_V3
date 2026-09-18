# Methods — Master/Train Staging Assembly

## Purpose

Stage the two labeled-by-construction pair tables built so far — the
standardized positive (interacting) pairs (`../01_cleaning_dedup_standardization/`)
and the in-distribution negative (random, non-interacting) pairs
(`../03_negative_pairs/`) — together in one location, so downstream
training/evaluation work has a single place to look. The two tables are
copied there **as separate files, not merged**: no label column, shuffling,
or train/test split has been applied yet. That is a deliberate, distinct
next step (see stage `../05_V3_train_holdout_sets/`), not an oversight.

## Inputs

| File | Role | Rows | Provenance |
|---|---|---|---|
| `huri_interactome_uniprotID_seqs_cleaned_deduped_standard.csv` | Positive (true interacting) pairs | 50,176 | Output of `../01_cleaning_dedup_standardization/` |
| `huri_negative_pairs_indistribution.csv` | Negative (in-distribution random, non-interacting) pairs | 50,176 | Output of `../03_negative_pairs/` |

Both tables share the same schema (`Protein_ID_1, Sequence_1, Protein_ID_2,
Sequence_2`) and are drawn from the same 8,008-protein universe, so they are
directly comparable/combinable later, but are kept apart for now.

## What was done

Each source file was copied (not moved — the original stage outputs stay in
place for provenance/re-run purposes) into the master staging folder, and
the copy was verified byte-for-byte against the source via SHA-256 checksum
and row count.

Both checksums matched exactly between source and destination in this
build; see the script's log output for the exact hashes and row counts of
your own run.

## Why kept separate

- Positive and negative pairs may need different treatment before they can
  be safely combined (e.g. a stratified train/validation/test split that
  keeps homologous protein clusters — per `protein_homology_table.csv` —
  on the same side of the split, which is easier to reason about with the
  two classes still distinguishable rather than a single loose CSV of
  pairs).
- Keeping them as two files with no label column avoids silently baking in
  an assumption about class balance, split ratios, or shuffling seed before
  those decisions are actually made.

## Not yet done (handled by the next stage)

Combining these into an actual training file needs, at minimum: adding a
label column (e.g. `Label = 1` / `0`), deciding the positive:negative ratio
to use, shuffling, and splitting into train/validation/(test) sets in a way
that respects the homology table (no near-duplicate protein pair straddling
a split boundary). That is exactly what `../05_V3_train_holdout_sets/` does.

## Files produced (not committed to this repo)

- `huri_interactome_uniprotID_seqs_cleaned_deduped_standard.csv` (copy)
- `huri_negative_pairs_indistribution.csv` (copy)

See `recipe.md` for the exact command.
