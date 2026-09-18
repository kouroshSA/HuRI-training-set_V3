#!/usr/bin/env python3
"""
Step 2 -- the "ordinary" real-protein negative set (RRS). Two modes:

  --mode in-distribution (default)
      Real proteins from the POSITIVE POOL, paired at random via a
      configuration model (stub-shuffling) so the negative degree
      distribution matches the positive graph's -- without a degree match, a
      model can learn to read the label off how often it has seen a protein
      rather than off the pair. Every candidate pool member has already
      interacted with something in the positive set; this set asks "does
      this specific OTHER pairing also interact?"

  --mode out-of-distribution
      Real proteins drawn uniformly at random from an EXTERNAL full-proteome
      pool (--proteome-csv), not restricted to proteins that appear in any
      positive pair. Most candidate proteins here have never been tested in
      the underlying screen at all -- this set asks a different question,
      "are two essentially-arbitrary proteins from the proteome a pair?",
      which is a much easier negative (no shared-degree structure to match,
      so no configuration model -- plain rejection-sampled uniform pairing).

**Paralog filter (both modes).** A random pair (X, Y) is rejected if it is a
"paralog shadow" of some positive pair (A, B) -- i.e. X is a close sequence
match to A (or B) and Y is a close match to B (or A). Such a pair could
plausibly be a true (but unassayed) interaction via conserved paralogous
binding, which would make it a bad negative example. Closeness is k-mer
Jaccard similarity (--paralog-k, --paralog-threshold; alignment-free, no
external tools -- see common.build_close_match_map for caveats). In
out-of-distribution mode the close-match map is built over the UNION of the
positive-pool proteins and the proteome pool, so a proteome protein that
happens to be a near-duplicate of a positive-pair partner is still caught.
Disable with --no-paralog-filter.

**Orientation duplication.** Matches step 1: the n_target unique unordered
negative pairs are each written out in both orientations, so this set has
the same row-count scale and position-symmetry as 01_positives.csv (default
n_target = n unique positive pairs, so output is 2 x n_target rows).

Default seed is 42, as specified for this step.

Writes:
  <outdir>/Parent_sequences/random_pairs.csv   (SEQ1,SEQ2,label=0, both orientations)
  <outdir>/Parent_sequences/random_pairs_stats.json
"""
import argparse
import csv
import random
import sys
import time

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from common import (read_rows, write_rows_lf, dump_json, canonical_key,
                     configuration_model_negative_pairs, build_close_match_map,
                     paralog_forbidden_keys)


def load_proteome_pool(path, id_col, seq_col, filter_col, filter_value):
    pool = {}  # orf_id -> sequence
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            if filter_col and row.get(filter_col) != filter_value:
                continue
            pool[row[id_col].strip()] = row[seq_col].strip()
    return pool


def uniform_random_negative_pairs(pool_sequences, n_target, seed, pos_keys, extra_forbidden,
                                   max_attempts_multiplier=50):
    """
    Rejection-sampled uniform random pairing over an arbitrary protein pool
    (no degree structure to match -- unlike configuration_model_negative_pairs,
    which is for pairing WITHIN the positive pool). Rejects self-pairs,
    duplicate negative pairs, positive-pair matches, and paralog-forbidden
    pairs. Raises if n_target isn't reached within a generous attempt budget
    (never silently returns short).
    """
    rng = random.Random(seed)
    pool = sorted(pool_sequences)  # sorted: deterministic despite str-hash salting
    neg_pairs = set()
    max_attempts = max(n_target * max_attempts_multiplier, 10000)
    attempts = 0
    while len(neg_pairs) < n_target and attempts < max_attempts:
        attempts += 1
        a, b = rng.sample(pool, 2)
        key = canonical_key(a, b)
        if key in pos_keys or key in neg_pairs or key in extra_forbidden:
            continue
        neg_pairs.add(key)
    if len(neg_pairs) < n_target:
        raise RuntimeError(
            f"uniform_random_negative_pairs: only reached {len(neg_pairs)}/{n_target} "
            f"unique negative pairs after {max_attempts} attempts (pool too small relative "
            f"to n_target, or --paralog-threshold too strict -- try a larger pool, lower "
            f"threshold, or --no-paralog-filter)"
        )
    result = sorted(tuple(sorted(k)) for k in neg_pairs)
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--positives", required=True, help="path to 01_positives.csv")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--mode", choices=["in-distribution", "out-of-distribution"], default="in-distribution")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--n-target", type=int, default=None,
                     help="unique negative pairs to generate; default = n unique positive pairs (1:1)")
    ap.add_argument("--paralog-k", type=int, default=4, help="k-mer size for paralog closeness")
    ap.add_argument("--paralog-threshold", type=float, default=0.5,
                     help="k-mer Jaccard similarity at/above which two proteins count as close matches")
    ap.add_argument("--no-paralog-filter", action="store_true")
    # out-of-distribution mode only
    ap.add_argument("--proteome-csv", default=None,
                     help="[out-of-distribution mode] CSV with the full candidate protein pool")
    ap.add_argument("--proteome-id-col", default="orf_id")
    ap.add_argument("--proteome-seq-col", default="sequence")
    ap.add_argument("--proteome-filter-col", default=None,
                     help="optional column to filter the proteome pool on, e.g. in_paper_scope_2026_manuscript")
    ap.add_argument("--proteome-filter-value", default="True",
                     help="value --proteome-filter-col must equal to keep a row (string-compared, as read from CSV)")
    args = ap.parse_args()

    if args.mode == "out-of-distribution" and not args.proteome_csv:
        raise SystemExit("--proteome-csv is required for --mode out-of-distribution")

    pos_rows = read_rows(args.positives)
    pos_pairs_all_rows = [(r[0], r[1]) for r in pos_rows]  # both orientations, as written by step 1
    unique_pos_pairs = sorted(set(tuple(sorted(p)) for p in pos_pairs_all_rows))
    pos_keys = set(canonical_key(a, b) for a, b in unique_pos_pairs)
    positive_proteins = set()
    for a, b in unique_pos_pairs:
        positive_proteins.add(a)
        positive_proteins.add(b)

    n_target = args.n_target if args.n_target is not None else len(unique_pos_pairs)

    proteome_pool = {}
    if args.mode == "out-of-distribution":
        proteome_pool = load_proteome_pool(args.proteome_csv, args.proteome_id_col, args.proteome_seq_col,
                                            args.proteome_filter_col, args.proteome_filter_value)

    # --- paralog filter setup ---
    forbidden = set()
    close_map_stats = {}
    if not args.no_paralog_filter:
        t0 = time.time()
        if args.mode == "in-distribution":
            paralog_check_pool = set(positive_proteins)
        else:
            paralog_check_pool = set(positive_proteins) | set(proteome_pool.values())
        close_map = build_close_match_map(paralog_check_pool, k=args.paralog_k, threshold=args.paralog_threshold)
        forbidden = paralog_forbidden_keys(unique_pos_pairs, close_map)
        elapsed = time.time() - t0
        group_sizes = [len(v) for v in close_map.values()]
        close_map_stats = {
            "paralog_k": args.paralog_k,
            "paralog_threshold": args.paralog_threshold,
            "n_distinct_proteins_checked": len(paralog_check_pool),
            "n_proteins_with_a_close_match_besides_self": sum(1 for s in group_sizes if s > 1),
            "max_close_group_size": max(group_sizes) if group_sizes else 0,
            "n_forbidden_paralog_shadow_pairs": len(forbidden),
            "build_time_seconds": round(elapsed, 2),
        }

    # --- negative-pair generation ---
    if args.mode == "in-distribution":
        neg_pairs, pearson_r = configuration_model_negative_pairs(
            unique_pos_pairs, n_target, seed=args.seed, extra_forbidden=forbidden)
        pool_size_stats = {"n_pool_proteins": len(positive_proteins)}
    else:
        pool_sequences = set(proteome_pool.values())
        neg_pairs = uniform_random_negative_pairs(pool_sequences, n_target, seed=args.seed,
                                                    pos_keys=pos_keys, extra_forbidden=forbidden)
        pearson_r = None  # not meaningful: most OOD-pool proteins have positive-degree 0 by design
        pool_size_stats = {
            "n_pool_orf_ids_after_filter": len(proteome_pool),
            "n_pool_distinct_sequences": len(pool_sequences),
            "n_pool_proteins_also_in_positive_pool": len(set(proteome_pool.values()) & positive_proteins),
        }

    rows = []
    for a, b in neg_pairs:
        rows.append((a, b, 0))
        rows.append((b, a, 0))  # orientation duplication, matching step 1

    out_path = f"{args.outdir}/Parent_sequences/random_pairs.csv"
    write_rows_lf(out_path, rows)

    # independent re-verification: re-read the written file, confirm 0 rows
    # collide with the forbidden paralog-shadow set (not just trust the
    # generator's internal bookkeeping)
    written_keys = set(canonical_key(r[0], r[1]) for r in read_rows(out_path))
    paralog_collisions_found = len(written_keys & forbidden) if forbidden else 0
    positive_collisions_found = len(written_keys & pos_keys)

    stats = {
        "mode": args.mode,
        "source_positives": args.positives,
        "seed": args.seed,
        "n_unique_positive_pairs": len(unique_pos_pairs),
        "n_target_negative_pairs": n_target,
        "n_negative_pairs_written": len(neg_pairs),
        "n_negative_rows_written": len(rows),
        "degree_match_pearson_r": pearson_r,
        "pool": pool_size_stats,
        "proteome_csv": args.proteome_csv,
        "paralog_filter_enabled": not args.no_paralog_filter,
        "paralog_filter": close_map_stats,
        "paralog_collisions_found_on_reread": paralog_collisions_found,
        "positive_pair_collisions_found_on_reread": positive_collisions_found,
        "output_file": out_path,
    }
    dump_json(f"{args.outdir}/Parent_sequences/random_pairs_stats.json", stats)

    pr_str = f"{pearson_r:.4f}" if pearson_r is not None else "N/A (out-of-distribution mode)"
    print(f"[02] mode={args.mode}: wrote {len(rows)} random pair rows ({len(neg_pairs)} unique pairs x2 "
          f"orientations, seed={args.seed}, degree-match Pearson r={pr_str})")
    if not args.no_paralog_filter:
        print(f"[02] paralog filter: {len(forbidden)} shadow pairs forbidden "
              f"(k={args.paralog_k}, threshold={args.paralog_threshold}); "
              f"{paralog_collisions_found} collisions found on re-read")
    print(f"[02] -> {out_path}")
    if paralog_collisions_found:
        raise SystemExit(f"[02] PARALOG-SHADOW LEAK: {paralog_collisions_found} pair(s) in the written "
                          f"file match the forbidden set -- this should not be possible, investigate")
    if positive_collisions_found:
        raise SystemExit(f"[02] POSITIVE-PAIR LEAK: {positive_collisions_found} pair(s) in the written "
                          f"file match a real positive pair -- this should not be possible, investigate")


if __name__ == "__main__":
    main()
