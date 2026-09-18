#!/usr/bin/env python3
"""
Step 1 -- normalize an arbitrary binary PPI positive set into the pipeline's
canonical working format: headerless, LF-terminated CSV with three fields

    SEQ1,SEQ2,label      label is always 1 here

Input is any CSV **with a header row** that has two columns holding the two
partners' amino-acid sequences; point --seq1-col / --seq2-col at them. No
column-name guessing is done on purpose -- the input schema varies by source
(YeRI ships protein1_sequence/protein2_sequence; other interactomes may not).

Cleaning applied per sequence: strip whitespace, uppercase, drop one trailing
stop-codon '*' if present. Rows are dropped (and counted) if, after cleaning,
either sequence is empty or contains a non-standard-amino-acid character --
every downstream random-sequence generator assumes the clean 20-letter
alphabet, so this is where that guarantee is established for the whole run.

**Orientation duplication.** Each distinct unordered pair is deduplicated to
one representative row, then written out in BOTH orientations -- (A,B) and
(B,A) -- unless it is a homodimer (A==B), which is written once (a literal
swap of a homodimer is byte-identical and would just be a duplicate row).
This matters beyond completeness: every downstream negative set (built by
substitution or degree-matched re-pairing) inherits this same row structure,
so position (SEQ1 vs SEQ2) never correlates with label on its own -- without
it, "which column has the more label-1-typical residue composition" could
become a shortcut a model learns instead of judging the pair.

Writes:
  <outdir>/01_positives.csv                canonical 3-field positive set (both orientations)
  <outdir>/Parent_sequences/positives.csv  identical copy, so the positive set lives alongside
                                            the other Parent_sequences/ sets it was built from
  <outdir>/01_positives_stats.json         measured diagnostics (not assumed)
"""
import argparse
import csv
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from common import clean_sequence, is_standard_aa, canonical_key, write_rows_lf, dump_json


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="input CSV with a header row")
    ap.add_argument("--seq1-col", required=True)
    ap.add_argument("--seq2-col", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    with open(args.input, newline="") as f:
        reader = csv.DictReader(f)
        rows_in = list(reader)

    n_input = len(rows_in)
    dropped_empty = 0
    dropped_nonstandard = 0
    dropped_duplicate_unordered = 0

    seen_keys = set()
    base_pairs = []  # one representative row per unique unordered pair
    for r in rows_in:
        s1 = clean_sequence(r[args.seq1_col])
        s2 = clean_sequence(r[args.seq2_col])
        if not s1 or not s2:
            dropped_empty += 1
            continue
        if not (is_standard_aa(s1) and is_standard_aa(s2)):
            dropped_nonstandard += 1
            continue
        key = canonical_key(s1, s2)
        if key in seen_keys:
            dropped_duplicate_unordered += 1
            continue
        seen_keys.add(key)
        base_pairs.append((s1, s2))

    positives = []
    n_homodimers = 0
    for s1, s2 in base_pairs:
        positives.append((s1, s2, 1))
        if s1 == s2:
            n_homodimers += 1
        else:
            positives.append((s2, s1, 1))  # explicit swap-orientation duplicate

    out_path = f"{args.outdir}/01_positives.csv"
    write_rows_lf(out_path, positives)
    parent_copy_path = f"{args.outdir}/Parent_sequences/positives.csv"
    write_rows_lf(parent_copy_path, positives)

    proteins = set()
    degree = {}
    for a, b, _ in positives:
        proteins.add(a)
        proteins.add(b)
        degree[a] = degree.get(a, 0) + 1
        degree[b] = degree.get(b, 0) + 1
    max_degree = max(degree.values()) if degree else 0

    stats = {
        "input_file": args.input,
        "n_input_rows": n_input,
        "dropped_empty_sequence": dropped_empty,
        "dropped_nonstandard_aa": dropped_nonstandard,
        "dropped_duplicate_unordered_pair": dropped_duplicate_unordered,
        "n_unique_unordered_pairs": len(base_pairs),
        "n_homodimers": n_homodimers,
        "n_positive_rows_written": len(positives),
        "orientation_duplication": "every non-homodimer pair written as both (A,B) and (B,A)",
        "n_distinct_proteins": len(proteins),
        "max_degree": max_degree,
        "output_file": out_path,
        "parent_sequences_copy": parent_copy_path,
    }
    dump_json(f"{args.outdir}/01_positives_stats.json", stats)

    print(f"[01] wrote {len(positives)} positive rows ({len(base_pairs)} unique unordered pairs, "
          f"{n_homodimers} homodimers, both orientations included) -> {out_path}")
    print(f"[01] also saved a copy -> {parent_copy_path}")
    if dropped_empty or dropped_nonstandard or dropped_duplicate_unordered:
        print(f"[01] dropped: {dropped_empty} empty, {dropped_nonstandard} non-standard-AA, "
              f"{dropped_duplicate_unordered} duplicate-unordered-pair rows")


if __name__ == "__main__":
    main()
