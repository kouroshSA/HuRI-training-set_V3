#!/usr/bin/env python3
"""
Extract the set of unique single proteins (Protein_ID + Sequence) referenced
anywhere in the cleaned/standardized HuRI table, and write them as a FASTA
file for all-vs-all homology searching.

Input:
  huri_interactome_uniprotID_seqs_cleaned_deduped_standard.csv (output of
  stage 01_cleaning_dedup_standardization)
Output:
  proteins/unique_proteins.fasta
  proteins/unique_proteins.csv   (ID, Length) lookup table

Usage:
  python3 03_extract_unique_proteins.py --input <standardized.csv> --outdir <dir>
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="path to huri_interactome_uniprotID_seqs_cleaned_deduped_standard.csv")
    ap.add_argument("--outdir", required=True, help="directory to write outputs into")
    args = ap.parse_args()

    input_csv = Path(args.input)
    outdir = Path(args.outdir)
    proteins_dir = outdir / "proteins"
    fasta_out = proteins_dir / "unique_proteins.fasta"
    table_out = proteins_dir / "unique_proteins.csv"
    log_file = outdir / "logs" / "03_extract_unique_proteins.log"

    lines = []

    def log(msg: str) -> None:
        print(msg)
        lines.append(msg)

    log(f"Reading {input_csv}")
    df = pd.read_csv(input_csv, dtype=str, keep_default_na=False)
    log(f"Rows read: {len(df)}")

    left = df[["Protein_ID_1", "Sequence_1"]].rename(
        columns={"Protein_ID_1": "Protein_ID", "Sequence_1": "Sequence"}
    )
    right = df[["Protein_ID_2", "Sequence_2"]].rename(
        columns={"Protein_ID_2": "Protein_ID", "Sequence_2": "Sequence"}
    )
    proteins = pd.concat([left, right], ignore_index=True).drop_duplicates()

    n_conflict = proteins.groupby("Protein_ID")["Sequence"].nunique()
    conflicts = n_conflict[n_conflict > 1]
    if len(conflicts):
        log(f"WARNING: {len(conflicts)} Protein_IDs map to more than one distinct "
            f"sequence; keeping the first sequence seen for each. IDs: "
            f"{list(conflicts.index)[:10]}{'...' if len(conflicts) > 10 else ''}")
        proteins = proteins.drop_duplicates(subset=["Protein_ID"], keep="first")

    proteins = proteins.sort_values("Protein_ID").reset_index(drop=True)
    proteins["Length"] = proteins["Sequence"].str.len()
    log(f"Unique proteins: {len(proteins)}")
    log(f"Sequence length range: {proteins['Length'].min()}-{proteins['Length'].max()} "
        f"(mean {proteins['Length'].mean():.1f})")

    proteins_dir.mkdir(parents=True, exist_ok=True)
    with open(fasta_out, "w") as fh:
        for pid, seq in zip(proteins["Protein_ID"], proteins["Sequence"]):
            fh.write(f">{pid}\n{seq}\n")
    log(f"Wrote FASTA ({len(proteins)} sequences) to {fasta_out}")

    proteins[["Protein_ID", "Length"]].to_csv(table_out, index=False)
    log(f"Wrote ID/length lookup table to {table_out}")

    log_file.parent.mkdir(parents=True, exist_ok=True)
    log_file.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
