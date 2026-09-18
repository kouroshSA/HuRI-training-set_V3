# Methods — Single-Protein Homology Table

## Purpose

`huri_interactome_uniprotID_seqs_cleaned_deduped_standard.csv` (from
`../01_cleaning_dedup_standardization/`) contains 50,176 interaction pairs
built from 8,008 unique human proteins. Before this table can be split into
train/validation/test sets (or used to build in-distribution and
out-of-distribution controls), we need to know which *individual proteins*
are paralogs or otherwise highly similar to each other. Without removing
near-duplicate proteins across splits, a held-out "test" protein can simply
be a near-identical copy of a protein the model saw in training, which
inflates apparent performance.

This step produces a **single-protein** homology table (not a table of
interaction pairs): for every pair of the 8,008 proteins that show any
detectable sequence similarity, it reports the percent identity and the
length of the aligned (similar) stretch, plus enough surrounding context
(protein lengths, alignment coverage, E-value) to choose a removal threshold
afterward. No identity/coverage threshold has been applied yet — that
decision is deliberately left for downstream analysis.

## Why an alignment tool was needed

There are 8,008 unique proteins, i.e. ~32.05 million possible pairs. A naive
all-pairs Smith-Waterman in Python is not tractable at that scale.
**mmseqs2** is the tool used here — a fast, actively-maintained local
alignment search tool (k-mer-indexed prefiltering followed by a gapped local
alignment), avoiding the full quadratic blow-up while still reporting
BLAST-style percent identity / alignment length / E-value per hit. Install it
into a dedicated conda environment:

```bash
conda create -n homology -c bioconda -c conda-forge mmseqs2 -y
```

## Step 1 — Extract the unique protein set

Every `(Protein_ID, Sequence)` appearing in either column-pair
(`Protein_ID_1`/`Sequence_1`, `Protein_ID_2`/`Sequence_2`) of the
cleaned/standardized interaction table was pooled and deduplicated by ID.

- **8,008 unique proteins**, no ID mapped to more than one distinct sequence.
- Sequence lengths range 25–6,885 residues (mean 450.9).
- Written to `proteins/unique_proteins.fasta` (for mmseqs2) and
  `proteins/unique_proteins.csv` (ID/length lookup used later for sanity
  checks).

## Step 2 — All-vs-all similarity search (mmseqs2)

```bash
mmseqs easy-search unique_proteins.fasta unique_proteins.fasta all_vs_all.m8 tmp \
    -s 7.5 -e 10 --max-seqs 2000 -a 1 \
    --format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov
```

Key parameter choices:

- **`-s 7.5`** (maximum sensitivity) so weak/remote paralogs aren't missed by
  the k-mer prefilter.
- **`-e 10`** — a deliberately permissive E-value *detection* threshold. This
  is not the paralog-removal threshold; it just controls which candidate
  hits are computed and reported at all. Being permissive here means the
  table can be filtered by %identity/coverage downstream without having to
  re-run the search.
- **`--max-seqs 2000`** — generous per-query hit cap so promiscuous protein
  families aren't truncated (only 8,008 sequences total, so this is cheap).
- **`-a 1`** — compute the full alignment backtrace. Without this flag,
  mmseqs2 still reports correct (gap-aware) percent identity and alignment
  length, but leaves `mismatch`/`gapopen` at 0 for every hit because those
  two columns are only populated from an actual traced alignment path. `-a 1`
  was added specifically so the mismatch/gap-open counts in the final table
  are real, not placeholders.

Runtime: ~16 seconds total on 28 threads for all 8,008×8,008 comparisons.
Self-hits (a protein vs. itself) were excluded downstream, not by the search
itself (`--add-self-matches 0` is the mmseqs default, but 8,003/8,008
proteins still produced a query==target row during the search and were
filtered out explicitly in step 3 below).

Raw output: `mmseqs/all_vs_all.m8` — 285,603 directed hit rows (a pair can
appear as both `A→B` and `B→A`, each potentially with slightly different
alignment statistics since the search is heuristic, not exhaustive).

## Step 3 — Build the single-protein homology table

For each unordered pair of proteins `(A, B)` (`A < B` alphabetically):

1. Self-hits (`query == target`) were dropped (8,003 rows).
2. If both directions `A→B` and `B→A` were found (111,236 of the 166,364
   pairs), the direction with the higher percent identity was kept (ties
   broken by bit score) — i.e. each pair is reported once, using its
   best-scoring alignment.
3. Protein lengths, coverage, and gap/mismatch counts were re-expressed
   relative to `Protein_A`/`Protein_B` rather than the winning direction's
   query/target labels.

Output: `protein_homology_table.csv` — **166,364 unordered protein pairs**
with detectable similarity (out of ~32.05 million possible pairs; mmseqs2's
k-mer prefilter is what makes searching only the plausible subset feasible).

### Columns

| Column | Meaning |
|---|---|
| `Protein_A`, `Protein_B` | The two protein IDs (alphabetically ordered; this is a single canonical row per unordered pair) |
| `Percent_Identity` | % identical residues over the aligned region (mmseqs `pident`) |
| `Alignment_Length` | Length, in residues, of the aligned (similar) stretch — this is the "length of the stretch with similarity" |
| `Length_A`, `Length_B` | Full length of each protein |
| `Coverage_Shorter_Protein` | `Alignment_Length / min(Length_A, Length_B)` — how much of the smaller protein the similar stretch covers. Can slightly exceed 1.0 for alignments with insertions (the alignment "length" counts gap columns too), e.g. up to 1.61 in this dataset for a handful of heavily-gapped hits |
| `Coverage_Longer_Protein` | `Alignment_Length / max(Length_A, Length_B)` |
| `Mismatches`, `Gap_Openings` | From the traced alignment (`-a 1`) |
| `E_value`, `Bit_Score` | mmseqs2 statistical significance / score for the reported alignment |

## Sanity checks performed

- 200 randomly sampled rows had their `Length_A`/`Length_B` cross-checked
  against the independently built `proteins/unique_proteins.csv` lookup —
  0 mismatches.
- Confirmed `Gap_Openings` is no longer trivially zero everywhere after
  adding `-a 1` (141,696 of 166,364 pairs have at least one gap opening;
  before adding `-a 1`, this column was uniformly 0, which is a known mmseqs2
  behavior when the backtrace isn't requested — see Step 2 notes).

## Identity distribution (for choosing a threshold)

| Percent identity band | Pairs |
|---|---|
| 90–100% | 223 |
| 70–<90% | 778 |
| 50–<70% | 11,989 |
| 30–<50% | 91,807 |
| <30% (down to the 14.2% floor mmseqs2 still reported at `-e 10`) | 61,567 |

How many *proteins* would be touched by a paralog-removal step at a few
candidate thresholds (i.e. how many of the 8,008 proteins appear in at least
one pair at or above that identity):

| Threshold (≥) | Pairs at/above | Unique proteins involved |
|---|---|---|
| 30% | 104,797 | 7,561 |
| 40% | 43,417 | 6,222 |
| 50% | 12,990 | 4,202 |
| 70% | 1,001 | 1,175 |
| 90% | 223 | 258 |

These numbers are only descriptive — no rows have been removed from any
interaction table at this stage. A common convention for "safe" train/test
protein splits in interaction-prediction work is to remove one side of any
pair above ~25–40% sequence identity with ≥~50% coverage of the shorter
protein (to avoid counting two unrelated proteins that merely share a short
motif as paralogs); the `Coverage_Shorter_Protein` column is included
specifically so that combined identity+coverage rule can be applied, rather
than identity alone.

## Files produced (not committed to this repo)

- `proteins/unique_proteins.fasta`, `proteins/unique_proteins.csv`
- `mmseqs/all_vs_all.m8` (raw directed hits)
- `protein_homology_table.csv` (final single-protein homology table)

## Next stage

`protein_homology_table.csv` (this deliverable) is used as-is — with no
threshold applied — by `../03_negative_pairs/` to build a homology-safe,
degree-matched set of in-distribution negative pairs.

See `recipe.md` for exact reproduction commands and `scripts/` for the code.
