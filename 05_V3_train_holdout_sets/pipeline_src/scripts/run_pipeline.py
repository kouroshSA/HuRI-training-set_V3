#!/usr/bin/env python3
"""
Orchestrator -- runs steps 01-08 end to end with one command.

Example:
  python run_pipeline.py \\
      --input /path/to/positive_set.csv \\
      --seq1-col protein1_sequence --seq2-col protein2_sequence \\
      --outdir /path/to/output_dir

Every step is also runnable standalone (see README.md) -- this just chains
them with the documented default seeds. Pass --help to see every knob;
anything not overridden here uses the per-script default.
"""
import argparse
import subprocess
import sys

HERE = __file__.rsplit("/", 1)[0]


def run(cmd):
    print(f"\n$ {' '.join(cmd)}")
    subprocess.run([sys.executable] + cmd, check=True, cwd=HERE)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True)
    ap.add_argument("--seq1-col", required=True)
    ap.add_argument("--seq2-col", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--n-replicates", type=int, default=10)
    ap.add_argument("--frac", type=float, default=0.10)
    ap.add_argument("--version-tag", default="V3")
    ap.add_argument("--deplete-decoys-from-heldout", action="store_true")
    ap.add_argument("--seed-rand-pairs", type=int, default=42)
    ap.add_argument("--seed-ps1", type=int, default=43)
    ap.add_argument("--seed-ps2", type=int, default=44)
    ap.add_argument("--seed-both", type=int, default=45)
    ap.add_argument("--seed-master-shuffle", type=int, default=46)
    ap.add_argument("--seed-base-shuffle", type=int, default=6000,
                     help="base seed for the final per-replicate depleted-training-set shuffle (uses base+k)")
    ap.add_argument("--paralog-k", type=int, default=4)
    ap.add_argument("--paralog-threshold", type=float, default=0.5)
    ap.add_argument("--no-paralog-filter", action="store_true")
    # RRS mode (step 2's "ordinary" real-protein negative set)
    ap.add_argument("--rrs-mode", choices=["in-distribution", "out-of-distribution"], default="in-distribution",
                     help="in-distribution: degree-matched pairing within the positive pool (default). "
                          "out-of-distribution: uniform random pairing from --proteome-csv, an external full-proteome pool.")
    ap.add_argument("--proteome-csv", default=None, help="required if --rrs-mode out-of-distribution")
    ap.add_argument("--proteome-id-col", default="orf_id")
    ap.add_argument("--proteome-seq-col", default="sequence")
    ap.add_argument("--proteome-filter-col", default=None)
    ap.add_argument("--proteome-filter-value", default="True")
    args = ap.parse_args()

    if args.rrs_mode == "out-of-distribution" and not args.proteome_csv:
        raise SystemExit("--proteome-csv is required for --rrs-mode out-of-distribution")

    run(["01_make_positives.py", "--input", args.input,
         "--seq1-col", args.seq1_col, "--seq2-col", args.seq2_col, "--outdir", args.outdir])

    rand_cmd = ["02_make_random_pairs.py",
                "--positives", f"{args.outdir}/01_positives.csv", "--outdir", args.outdir,
                "--mode", args.rrs_mode,
                "--seed", str(args.seed_rand_pairs),
                "--paralog-k", str(args.paralog_k), "--paralog-threshold", str(args.paralog_threshold)]
    if args.rrs_mode == "out-of-distribution":
        rand_cmd += ["--proteome-csv", args.proteome_csv,
                     "--proteome-id-col", args.proteome_id_col, "--proteome-seq-col", args.proteome_seq_col]
        if args.proteome_filter_col:
            rand_cmd += ["--proteome-filter-col", args.proteome_filter_col,
                         "--proteome-filter-value", args.proteome_filter_value]
    if args.no_paralog_filter:
        rand_cmd.append("--no-paralog-filter")
    run(rand_cmd)

    run(["03_make_random_substitution_sets.py",
         "--positives", f"{args.outdir}/01_positives.csv", "--outdir", args.outdir,
         "--seed-ps1", str(args.seed_ps1), "--seed-ps2", str(args.seed_ps2), "--seed-both", str(args.seed_both)])

    run(["04_make_master_mix.py", "--outdir", args.outdir, "--seed", str(args.seed_master_shuffle)])

    mccv_cmd = ["05_make_mccv_replicates.py", "--outdir", args.outdir,
                "--n-replicates", str(args.n_replicates), "--frac", str(args.frac),
                "--version-tag", args.version_tag,
                "--seed-base-shuffle", str(args.seed_base_shuffle),
                "--mode", args.rrs_mode,
                "--paralog-k", str(args.paralog_k), "--paralog-threshold", str(args.paralog_threshold)]
    if args.rrs_mode == "out-of-distribution":
        mccv_cmd += ["--proteome-csv", args.proteome_csv,
                     "--proteome-id-col", args.proteome_id_col, "--proteome-seq-col", args.proteome_seq_col]
        if args.proteome_filter_col:
            mccv_cmd += ["--proteome-filter-col", args.proteome_filter_col,
                         "--proteome-filter-value", args.proteome_filter_value]
    if args.deplete_decoys_from_heldout:
        mccv_cmd.append("--deplete-decoys-from-heldout")
    if args.no_paralog_filter:
        mccv_cmd.append("--no-paralog-filter")
    run(mccv_cmd)

    run(["06_convert_format_ppiDCE_ppiBTEP.py", "--outdir", args.outdir])
    run(["07_convert_format_ppiGPLM.py", "--outdir", args.outdir])
    run(["08_make_manifest_and_provenance.py", "--outdir", args.outdir, "--source-positive-file", args.input])

    print(f"\n[run_pipeline] complete -> {args.outdir}")


if __name__ == "__main__":
    main()
