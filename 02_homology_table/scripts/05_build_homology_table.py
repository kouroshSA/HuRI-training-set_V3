#!/usr/bin/env python3
"""
Turn the raw mmseqs2 all-vs-all hits into a single-protein homology table:
one row per unordered pair of proteins (A, B) that showed any detectable
local similarity, with the percent identity and the length of the aligned
(similar) stretch, plus enough context (protein lengths, coverage, E-value)
to later choose an identity/coverage threshold for removing paralog pairs
or for building in-/out-of-distribution controls.

Input:
  <outdir>/mmseqs/all_vs_all.m8   (from 04_run_all_vs_all_search.sh; columns:
    query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,
    evalue,bits,qlen,tlen,qcov,tcov)
  <outdir>/proteins/unique_proteins.csv  (Protein_ID, Length lookup, for sanity check)

Output:
  <outdir>/protein_homology_table.csv

Each unordered pair (A, B) can appear in the raw hits as up to two directed
rows (A as query vs B as target, and vice versa). Local alignment is not
always perfectly symmetric under a heuristic search, so for each pair we
keep the single best-scoring direction (highest percent identity, ties
broken by bit score) and report that alignment's identity and length.

Usage:
  python3 05_build_homology_table.py --outdir <dir>
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

COLUMNS = [
    "query", "target", "pident", "alnlen", "mismatch", "gapopen",
    "qstart", "qend", "tstart", "tend", "evalue", "bits", "qlen", "tlen",
    "qcov", "tcov",
]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True,
                     help="directory containing mmseqs/all_vs_all.m8 and proteins/unique_proteins.csv; "
                          "also where protein_homology_table.csv is written")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    m8_file = outdir / "mmseqs" / "all_vs_all.m8"
    lengths_file = outdir / "proteins" / "unique_proteins.csv"
    output_table = outdir / "protein_homology_table.csv"
    log_file = outdir / "logs" / "05_build_homology_table.log"

    lines = []

    def log(msg: str) -> None:
        print(msg)
        lines.append(msg)

    log(f"Reading {m8_file}")
    hits = pd.read_csv(m8_file, sep="\t", header=None, names=COLUMNS)
    log(f"Raw directed hits: {len(hits)}")

    # Drop self-hits (a protein "matching itself"); mmseqs was run with
    # --add-self-matches 0 so these shouldn't be present, but filter
    # defensively in case of upstream ID collisions.
    n_before = len(hits)
    hits = hits[hits["query"] != hits["target"]].copy()
    log(f"Removed {n_before - len(hits)} self-hit rows (query == target)")

    # Canonicalize each pair regardless of query/target order.
    pair_a = hits[["query", "target"]].min(axis=1)
    pair_b = hits[["query", "target"]].max(axis=1)
    hits["Protein_A"] = pair_a
    hits["Protein_B"] = pair_b

    n_directed = len(hits)
    n_pairs = hits[["Protein_A", "Protein_B"]].drop_duplicates().shape[0]
    n_both_dir = n_directed - n_pairs
    log(f"Unordered protein pairs with detectable similarity: {n_pairs} "
        f"({n_both_dir} of them were seen in both directions)")

    # For each unordered pair, keep the best-scoring direction: highest
    # percent identity, ties broken by bit score.
    hits = hits.sort_values(["pident", "bits"], ascending=False)
    best = hits.drop_duplicates(subset=["Protein_A", "Protein_B"], keep="first").copy()

    # Re-express length/coverage columns relative to Protein_A / Protein_B
    # rather than query/target (whichever happened to win as best hit).
    is_a_query = best["query"] == best["Protein_A"]
    best["Length_A"] = best["qlen"].where(is_a_query, best["tlen"])
    best["Length_B"] = best["tlen"].where(is_a_query, best["qlen"])

    best["Alignment_Length"] = best["alnlen"]
    best["Percent_Identity"] = best["pident"]
    shorter = best[["Length_A", "Length_B"]].min(axis=1)
    longer = best[["Length_A", "Length_B"]].max(axis=1)
    best["Coverage_Shorter_Protein"] = (best["Alignment_Length"] / shorter).round(4)
    best["Coverage_Longer_Protein"] = (best["Alignment_Length"] / longer).round(4)

    out = best[[
        "Protein_A", "Protein_B",
        "Percent_Identity", "Alignment_Length",
        "Length_A", "Length_B",
        "Coverage_Shorter_Protein", "Coverage_Longer_Protein",
        "mismatch", "gapopen", "evalue", "bits",
    ]].rename(columns={
        "mismatch": "Mismatches",
        "gapopen": "Gap_Openings",
        "evalue": "E_value",
        "bits": "Bit_Score",
    })
    out = out.sort_values(["Percent_Identity", "Alignment_Length"], ascending=False)
    out = out.reset_index(drop=True)

    out.to_csv(output_table, index=False)
    log(f"Wrote homology table ({len(out)} pairs) to {output_table}")

    # Sanity check against the independently computed length lookup.
    lengths = pd.read_csv(lengths_file, dtype={"Protein_ID": str}).set_index("Protein_ID")["Length"]
    mismatched_len = 0
    for _, row in out.sample(min(200, len(out)), random_state=0).iterrows():
        if lengths.get(row["Protein_A"]) != row["Length_A"] or lengths.get(row["Protein_B"]) != row["Length_B"]:
            mismatched_len += 1
    log(f"Spot-check: {mismatched_len}/200 sampled rows had a length mismatch against "
        f"the independent unique_proteins.csv lookup (expect 0)")

    # A quick distribution summary to help pick a threshold later.
    for lo, hi in [(90, 100), (70, 90), (50, 70), (30, 50), (0, 30)]:
        n = ((out["Percent_Identity"] >= lo) & (out["Percent_Identity"] < hi if hi < 100 else out["Percent_Identity"] <= hi)).sum()
        log(f"Pairs with {lo}-{hi}% identity: {n}")

    log_file.parent.mkdir(parents=True, exist_ok=True)
    log_file.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
