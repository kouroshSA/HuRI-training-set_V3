#!/usr/bin/env python3
"""
Scan the cleaned/deduplicated HuRI table for non-standard amino acid symbols
in the sequence columns (never the ID columns) and replace them with a
standard residue, logging every change to a manifest.

Replacement rule:
  - If the symbol has one unambiguous closest standard amino acid (e.g.
    Selenocysteine, U, which is a Cys analogue), replace with that residue.
  - If the symbol is ambiguous / unknown (could be more than one residue, or
    is a catch-all placeholder), replace with Alanine (A).

Replacement table used:
  U (Selenocysteine)      -> C (Cysteine)   - direct structural analogue of Cys
  O (Pyrrolysine)         -> K (Lysine)     - direct structural analogue of Lys
  B (Asx: Asp or Asn)     -> A (Alanine)    - ambiguous, no single closest residue
  Z (Glx: Glu or Gln)     -> A (Alanine)    - ambiguous, no single closest residue
  J (Xle: Leu or Ile)     -> A (Alanine)    - ambiguous, no single closest residue
  X (unknown residue)     -> A (Alanine)    - unknown, no closest residue
  any other non-standard  -> A (Alanine)    - fallback for anything unforeseen

Input:
  huri_interactome_uniprotID_seqs_cleaned_deduped.csv (output of 01_clean_and_dedup.py)
Outputs:
  huri_interactome_uniprotID_seqs_cleaned_deduped_standard.csv
  huri_interactome_uniprotID_seqs_cleaned_deduped_standard_manifest.csv

Usage:
  python3 02_standardize_amino_acids.py --input <cleaned_deduped.csv> --outdir <dir>
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

REPLACEMENT = {
    "U": ("C", "Selenocysteine -> closest standard residue Cysteine"),
    "O": ("K", "Pyrrolysine -> closest standard residue Lysine"),
    "B": ("A", "Asx (Asp/Asn ambiguity) -> ambiguous, defaulted to Alanine"),
    "Z": ("A", "Glx (Glu/Gln ambiguity) -> ambiguous, defaulted to Alanine"),
    "J": ("A", "Xle (Leu/Ile ambiguity) -> ambiguous, defaulted to Alanine"),
    "X": ("A", "Unknown residue -> ambiguous, defaulted to Alanine"),
}
FALLBACK = ("A", "Unrecognized non-standard symbol -> defaulted to Alanine")

ID_COLS = {"Sequence_1": "Protein_ID_1", "Sequence_2": "Protein_ID_2"}


def standardize_sequence(seq: str, row_idx: int, id_col: str, id_val: str, seq_col: str, manifest_rows: list) -> str:
    if not seq:
        return seq
    chars = list(seq)
    for pos, ch in enumerate(chars):
        if ch in STANDARD_AA:
            continue
        replacement, reason = REPLACEMENT.get(ch, FALLBACK)
        manifest_rows.append(
            {
                "Row_Index": row_idx,
                "Protein_ID": id_val,
                "Sequence_Column": seq_col,
                "Position_1based": pos + 1,
                "Original_Residue": ch,
                "Replacement_Residue": replacement,
                "Reason": reason,
            }
        )
        chars[pos] = replacement
    return "".join(chars)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="path to huri_interactome_uniprotID_seqs_cleaned_deduped.csv")
    ap.add_argument("--outdir", required=True, help="directory to write outputs into")
    args = ap.parse_args()

    input_csv = Path(args.input)
    outdir = Path(args.outdir)
    output_csv = outdir / "huri_interactome_uniprotID_seqs_cleaned_deduped_standard.csv"
    manifest_csv = outdir / "huri_interactome_uniprotID_seqs_cleaned_deduped_standard_manifest.csv"
    log_file = outdir / "logs" / "02_standardize_amino_acids.log"

    lines = []

    def log(msg: str) -> None:
        print(msg)
        lines.append(msg)

    log(f"Reading {input_csv}")
    df = pd.read_csv(input_csv, dtype=str, keep_default_na=False)
    log(f"Rows read: {len(df)}")

    manifest_rows: list = []
    for seq_col, id_col in ID_COLS.items():
        df[seq_col] = [
            standardize_sequence(seq, idx, id_col, pid, seq_col, manifest_rows)
            for idx, (seq, pid) in enumerate(zip(df[seq_col], df[id_col]))
        ]

    manifest = pd.DataFrame(
        manifest_rows,
        columns=[
            "Row_Index",
            "Protein_ID",
            "Sequence_Column",
            "Position_1based",
            "Original_Residue",
            "Replacement_Residue",
            "Reason",
        ],
    )

    log(f"Total residue substitutions made: {len(manifest)}")
    if len(manifest):
        counts = manifest.groupby(["Original_Residue", "Replacement_Residue"]).size()
        for (orig, repl), n in counts.items():
            log(f"  {orig} -> {repl}: {n} occurrence(s)")
        log(f"Sequences affected: {manifest[['Row_Index', 'Sequence_Column']].drop_duplicates().shape[0]}")
        log(f"Unique Protein_IDs affected: {manifest['Protein_ID'].nunique()}")

    outdir.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_csv, index=False)
    log(f"Wrote standardized table ({len(df)} rows) to {output_csv}")

    manifest.to_csv(manifest_csv, index=False)
    log(f"Wrote manifest ({len(manifest)} rows) to {manifest_csv}")

    log_file.parent.mkdir(parents=True, exist_ok=True)
    log_file.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
