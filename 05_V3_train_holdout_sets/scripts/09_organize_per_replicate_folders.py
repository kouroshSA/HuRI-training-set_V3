#!/usr/bin/env python3
"""
Reorganize the pipeline's canonical MCCV/ output (grouped by file TYPE across
all 10 replicates) into one self-contained folder PER REPLICATE --
HuRI-V3-1 .. HuRI-V3-10 -- each holding that replicate's depleted training
set plus its own PRS, RRS, and three random-substitution decoy holdouts.

The canonical <build>/MCCV/ tree is left untouched; this only adds a second,
replicate-centric view of the same files (each file is copied, then
verified byte-identical via SHA-256 against its source).

Label convention inherited unchanged from the pipeline (see Methods.md):
  1 = true interacting pair (PRS / positives)
  0 = every negative type -- RRS (real-protein random pairs) AND all three
      random-substitution decoys (ps1_random, ps2_random, both_random)

Per-replicate folder contents:
  training_set.csv           depleted_training_set-V3-{k}.csv (shuffled, 1:4 pos:neg)
  PRS.csv                    held-out true positives (label 1)
  RRS.csv                    held-out real-protein random pairs (label 0)
  ps1_random_control.csv     held-out: partner 1 replaced by random seq, partner 2 native (label 0)
  ps2_random_control.csv     held-out: partner 1 native, partner 2 replaced by random seq (label 0)
  both_random_control.csv    held-out: both partners replaced by random seqs (label 0)
  replicate_summary.json     this replicate's row counts / seeds, pulled from mccv_stats.json

Usage:
  python3 09_organize_per_replicate_folders.py --build-dir <build> --outdir <dir>
"""
import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

N_REPLICATES = 10
TAG = "V3"

FILES = [
    ("training_sets", f"depleted_training_set-{TAG}-{{k}}.csv", "training_set.csv"),
    ("PRS-RRS", f"PRS-{TAG}-{{k}}.csv", "PRS.csv"),
    ("PRS-RRS", f"RRS-{TAG}-{{k}}.csv", "RRS.csv"),
    ("random_controls", f"ps1random-{TAG}-{{k}}.csv", "ps1_random_control.csv"),
    ("random_controls", f"ps2random-{TAG}-{{k}}.csv", "ps2_random_control.csv"),
    ("random_controls", f"bothrandom-{TAG}-{{k}}.csv", "both_random_control.csv"),
]


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build-dir", required=True, help="pipeline output root (contains MCCV/)")
    ap.add_argument("--outdir", required=True, help="directory under which HuRI-V3-1..10 folders are written")
    args = ap.parse_args()

    build = Path(args.build_dir)
    outdir = Path(args.outdir)
    mccv = build / "MCCV"

    mccv_stats = json.loads((mccv / "mccv_stats.json").read_text())
    stats_by_replicate = {r["replicate"]: r for r in mccv_stats["replicates"]}

    for k in range(1, N_REPLICATES + 1):
        dest_dir = outdir / f"HuRI-{TAG}-{k}"
        dest_dir.mkdir(parents=True, exist_ok=True)

        for subdir, src_pattern, dest_name in FILES:
            src = mccv / subdir / src_pattern.format(k=k)
            dst = dest_dir / dest_name
            shutil.copyfile(src, dst)
            src_sum, dst_sum = sha256_of(src), sha256_of(dst)
            if src_sum != dst_sum:
                sys.exit(f"CHECKSUM MISMATCH copying {src} -> {dst}")

        summary = dict(stats_by_replicate[k])
        summary["source_build_dir"] = str(build)
        summary["label_convention"] = "1 = true interacting pair (PRS); 0 = RRS or any random-substitution decoy"
        (dest_dir / "replicate_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

        print(f"[09] HuRI-{TAG}-{k}: training_set={summary['n_depleted_training_rows']} rows, "
              f"PRS={summary['n_prs']} RRS={summary['n_rrs']} "
              f"ps1ctrl={summary['n_ps1_control']} ps2ctrl={summary['n_ps2_control']} "
              f"bothctrl={summary['n_both_control']} -- all copies verified byte-identical")

    print(f"[09] done: {N_REPLICATES} per-replicate folders written under {outdir}")


if __name__ == "__main__":
    main()
