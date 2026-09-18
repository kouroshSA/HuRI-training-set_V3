#!/usr/bin/env python3
"""
Step 4 -- master training mix: concatenate the positive set and all four
parent negative sets (random_pairs -- in-distribution or out-of-distribution,
whichever mode step 2 was run in -- ps1_random, ps2_random, both_random),
then shuffle once with a fixed seed (default 46).

both_random is included in the bulk training mix (not held out for eval
only) -- this is a deliberate convention choice carried over from the most
recent generalization of this pipeline (the Synechocystis build), and it is
why the resulting class ratio is 1 positive : 4 negative rather than 1:3.
If you want both_random eval-only instead, drop it before calling this
script and the ratio reverts to 1:3.

Writes:
  <outdir>/master_training_mix.csv
  <outdir>/master_training_mix_stats.json
"""
import argparse
import random
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from common import read_rows, write_rows_lf, dump_json


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=46)
    args = ap.parse_args()

    sources = {
        "positives": f"{args.outdir}/01_positives.csv",
        "random_pairs": f"{args.outdir}/Parent_sequences/random_pairs.csv",
        "ps1_random": f"{args.outdir}/Parent_sequences/ps1_random.csv",
        "ps2_random": f"{args.outdir}/Parent_sequences/ps2_random.csv",
        "both_random": f"{args.outdir}/Parent_sequences/both_random.csv",
    }

    all_rows = []
    counts = {}
    for name, path in sources.items():
        rows = read_rows(path)
        counts[name] = len(rows)
        all_rows.extend(rows)

    rng = random.Random(args.seed)
    rng.shuffle(all_rows)

    out_path = f"{args.outdir}/master_training_mix.csv"
    write_rows_lf(out_path, all_rows)

    n_pos = counts["positives"]
    n_neg = sum(v for k, v in counts.items() if k != "positives")
    stats = {
        "seed": args.seed,
        "sources": sources,
        "counts_per_source": counts,
        "n_total_rows": len(all_rows),
        "n_positive": n_pos,
        "n_negative": n_neg,
        "class_ratio_pos_to_neg": f"1:{n_neg / n_pos:.2f}" if n_pos else None,
        "output_file": out_path,
    }
    dump_json(f"{args.outdir}/master_training_mix_stats.json", stats)

    print(f"[04] wrote master_training_mix.csv: {len(all_rows)} rows "
          f"({n_pos} pos : {n_neg} neg = 1:{n_neg / n_pos:.2f}), seed={args.seed}")


if __name__ == "__main__":
    main()
