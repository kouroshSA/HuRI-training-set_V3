#!/usr/bin/env python3
"""
Clean and deduplicate the raw HuRI interactome table.

Input: a CSV with columns Protein_ID_1, Sequence_1, Protein_ID_2, Sequence_2
(the raw HuRI export used by this project lives at data/raw/ in this repo).

Steps (in order):
  1. Drop any row missing a sequence for either partner (literal "NAN" string,
     case-insensitive, or a blank/whitespace-only field).
  2. Strip stop-codon symbols ("*") from the end (or anywhere) of a sequence,
     if present.
  3. Drop exact repeated pairs (same Protein_ID_1, Protein_ID_2 in the same order).
  4. Drop pairs that appear again with the two IDs swapped (A,B) vs (B,A).
  5. Write the cleaned, deduplicated table.

Usage:
  python3 01_clean_and_dedup.py --input <raw.csv> --outdir <dir>
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

ID_COLS = ["Protein_ID_1", "Protein_ID_2"]
SEQ_COLS = ["Sequence_1", "Sequence_2"]


def is_missing(series: pd.Series) -> pd.Series:
    """A sequence field counts as missing if it's blank or the literal string
    'NAN' (any case), which is how this file marks an unmapped sequence."""
    return series.fillna("").str.strip().str.upper().isin(["NAN", ""])


def strip_stop_codons(series: pd.Series):
    """Remove '*' stop-codon symbols from sequences. Returns the cleaned
    series and the number of sequences that were changed."""
    changed = series.str.contains(r"\*", regex=True)
    cleaned = series.str.replace(r"\*", "", regex=True)
    return cleaned, int(changed.sum())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="path to the raw HuRI CSV")
    ap.add_argument("--outdir", required=True, help="directory to write outputs into")
    args = ap.parse_args()

    input_csv = Path(args.input)
    outdir = Path(args.outdir)
    output_csv = outdir / "huri_interactome_uniprotID_seqs_cleaned_deduped.csv"
    log_file = outdir / "logs" / "01_clean_and_dedup.log"

    lines = []

    def log(msg: str) -> None:
        print(msg)
        lines.append(msg)

    log(f"Reading {input_csv}")
    df = pd.read_csv(input_csv, dtype=str, keep_default_na=False)
    n_start = len(df)
    log(f"Rows read: {n_start}")

    # ---- Step 1: drop rows with a missing sequence on either side ----
    missing_mask = is_missing(df["Sequence_1"]) | is_missing(df["Sequence_2"])
    n_missing = int(missing_mask.sum())
    df = df.loc[~missing_mask].copy()
    log(f"Step 1 - removed {n_missing} rows missing Sequence_1 and/or Sequence_2 "
        f"(literal 'NAN' or blank). Rows remaining: {len(df)}")

    # ---- Step 2: strip stop-codon symbols ----
    n_star_total = 0
    for col in SEQ_COLS:
        df[col], n_changed = strip_stop_codons(df[col])
        n_star_total += n_changed
    log(f"Step 2 - sequences containing a '*' stop-codon symbol that were stripped: "
        f"{n_star_total} (0 expected for this dataset; symbols are removed if ever present)")

    # ---- Step 3: drop exact repeated pairs (same order) ----
    n_before = len(df)
    df = df.drop_duplicates(subset=ID_COLS, keep="first").copy()
    log(f"Step 3 - removed {n_before - len(df)} exact repeated (Protein_ID_1, Protein_ID_2) "
        f"rows, same order. Rows remaining: {len(df)}")

    # ---- Step 4: drop swapped-orientation duplicates (A,B) == (B,A) ----
    n_before = len(df)
    pair_key = df.apply(lambda r: tuple(sorted((r["Protein_ID_1"], r["Protein_ID_2"]))), axis=1)
    df["_pair_key"] = pair_key
    df = df.drop_duplicates(subset=["_pair_key"], keep="first").copy()
    df = df.drop(columns=["_pair_key"])
    log(f"Step 4 - removed {n_before - len(df)} rows that repeated an already-seen pair "
        f"in swapped orientation (A,B) vs (B,A). Rows remaining: {len(df)}")

    n_self = int((df["Protein_ID_1"] == df["Protein_ID_2"]).sum())
    log(f"Note - {n_self} self-interaction (homodimer) rows remain (Protein_ID_1 == "
        f"Protein_ID_2); these were not removed, only true duplicate/swapped rows were.")

    # ---- Step 5: save ----
    outdir.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_csv, index=False)
    log(f"Step 5 - wrote {len(df)} rows to {output_csv}")
    log(f"Total rows removed: {n_start - len(df)} ({n_start} -> {len(df)})")

    log_file.parent.mkdir(parents=True, exist_ok=True)
    log_file.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
