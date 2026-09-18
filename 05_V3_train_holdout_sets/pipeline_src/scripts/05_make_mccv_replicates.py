#!/usr/bin/env python3
"""
Step 5 -- ten (default) Monte-Carlo cross-validation (MCCV) replicates.

MCCV, not k-fold: each replicate draws its own independent random sample, so
evaluation sets from different replicates may overlap -- what matters is that
replicate k's training set is clean of *that replicate's own* held-out pairs.

For each replicate k = 1..N (--n-replicates, default 10), a fraction
(--frac, default 0.10 = 10%) of each of the five parent pools is sampled
independently (own seed per pool per replicate, see seed table in README):

  PRS-V{tag}-{k}.csv          10% of unique positive pairs
  RRS-V{tag}-{k}.csv          10% of the random-pair set (step 2's random_pairs.csv --
                               in-distribution or out-of-distribution, whichever mode was used)
  random_controls/ps1random-V{tag}-{k}.csv    10% of the ps1_random parent set
  random_controls/ps2random-V{tag}-{k}.csv    10% of the ps2_random parent set
  random_controls/bothrandom-V{tag}-{k}.csv   10% of the both_random parent set

**Depletion, both orientations.** depleted_training_set-V{tag}-{k}.csv is the
master training mix with every row whose canonical key (orientation-
independent: {A,B} and {B,A} are the same key) matches one of this
replicate's selected pairs removed. Because canonical_key() is used, this
removes BOTH the (A,B) row and the (B,A) row from the master mix in one
pass -- an explicit check re-verifies this for every replicate (see
"orientation_leak_check" in mccv_stats.json) rather than trusting that by
construction.

**Final shuffle.** After depletion, each depleted training set is re-shuffled
with its own seed (--seed-base-shuffle + k, default 6000+k) before being
written -- an explicit mixing step, not just "whatever order survived
depletion of an already-shuffled master file."

**Paralog verification for RRS.** Because the random-pair parent set (step 2)
is already paralog-filtered against the full positive set, any subsample of
it -- including each replicate's RRS-V{tag}-{k} -- inherits that property.
This script re-verifies it explicitly per replicate (re-derives the same
forbidden-pair set step 2 used, from the same positives file and the same
--paralog-k/--paralog-threshold, and checks RRS-V{tag}-{k} against it)
rather than assuming inheritance is bug-free. If step 2 was run in
--mode out-of-distribution, pass the SAME --mode out-of-distribution and
--proteome-csv (and proteome column args) here too, so the re-derived
forbidden set is built over the same pool step 2 used -- otherwise this
verification checks against the wrong (in-distribution-only) close-match
map and isn't meaningful. Disable with --no-paralog-filter to skip (matches
step 2's flag; skip only if step 2 was also run with it).

Optional --deplete-decoys-from-heldout: also removes, from the depleted
training set, ps1_random/ps2_random/both_random rows in the master mix whose
*originating* positive pair (the row of 01_positives.csv they were built
from, tracked by row alignment through step 3) is among this replicate's
held-out PRS pairs. Off by default, matching this pipeline's documented
lineage (the Synechocystis build left this undone and flagged it explicitly
as a known, accepted gap) -- turn it on for a stricter leakage guarantee.

Writes, under <outdir>/MCCV/:
  training_sets/depleted_training_set-V{tag}-{k}.csv
  PRS-RRS/PRS-V{tag}-{k}.csv, RRS-V{tag}-{k}.csv
  random_controls/ps1random-V{tag}-{k}.csv, ps2random-V{tag}-{k}.csv, bothrandom-V{tag}-{k}.csv
  mccv_stats.json   (per-replicate counts, seeds, leak-check and paralog-check results)
"""
import argparse
import csv
import random
import sys
import time

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from common import (read_rows, write_rows_lf, dump_json, canonical_key, mccv_sample,
                     build_close_match_map, paralog_forbidden_keys)


def load_proteome_pool(path, id_col, seq_col, filter_col, filter_value):
    pool = {}
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            if filter_col and row.get(filter_col) != filter_value:
                continue
            pool[row[id_col].strip()] = row[seq_col].strip()
    return pool


def load_pairs(path):
    return [(r[0], r[1], int(r[2])) for r in read_rows(path)]


def unique_by_key(rows):
    """First-occurrence de-dup by canonical key; returns list of full rows."""
    seen = set()
    out = []
    for row in rows:
        k = canonical_key(row[0], row[1])
        if k not in seen:
            seen.add(k)
            out.append(row)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--n-replicates", type=int, default=10)
    ap.add_argument("--frac", type=float, default=0.10)
    ap.add_argument("--version-tag", default="V3")
    ap.add_argument("--deplete-decoys-from-heldout", action="store_true")
    # seed bases -- one addend per pool type, replicate k uses base+k
    ap.add_argument("--seed-base-prs", type=int, default=1000)
    ap.add_argument("--seed-base-rrs", type=int, default=2000)
    ap.add_argument("--seed-base-ps1ctrl", type=int, default=3000)
    ap.add_argument("--seed-base-ps2ctrl", type=int, default=4000)
    ap.add_argument("--seed-base-bothctrl", type=int, default=5000)
    ap.add_argument("--seed-base-shuffle", type=int, default=6000)
    # paralog re-verification (must match step 2's settings to be meaningful)
    ap.add_argument("--paralog-k", type=int, default=4)
    ap.add_argument("--paralog-threshold", type=float, default=0.5)
    ap.add_argument("--no-paralog-filter", action="store_true")
    # must match whatever step 2 was run with -- see docstring
    ap.add_argument("--mode", choices=["in-distribution", "out-of-distribution"], default="in-distribution")
    ap.add_argument("--proteome-csv", default=None)
    ap.add_argument("--proteome-id-col", default="orf_id")
    ap.add_argument("--proteome-seq-col", default="sequence")
    ap.add_argument("--proteome-filter-col", default=None)
    ap.add_argument("--proteome-filter-value", default="True")
    args = ap.parse_args()

    if args.mode == "out-of-distribution" and not args.no_paralog_filter and not args.proteome_csv:
        raise SystemExit("--proteome-csv is required when --mode out-of-distribution and the paralog "
                          "filter is enabled (pass --no-paralog-filter to skip RRS paralog verification instead)")

    tag = args.version_tag
    positives = load_pairs(f"{args.outdir}/01_positives.csv")
    rand_pairs = load_pairs(f"{args.outdir}/Parent_sequences/random_pairs.csv")
    ps1r = load_pairs(f"{args.outdir}/Parent_sequences/ps1_random.csv")
    ps2r = load_pairs(f"{args.outdir}/Parent_sequences/ps2_random.csv")
    bothr = load_pairs(f"{args.outdir}/Parent_sequences/both_random.csv")
    master = load_pairs(f"{args.outdir}/master_training_mix.csv")

    pos_unique = unique_by_key(positives)
    rand_unique = unique_by_key(rand_pairs)
    # ps1/ps2/both decoy rows are already unique by construction (fake side is
    # globally unique), but de-dup defensively for consistency.
    ps1_unique = unique_by_key(ps1r)
    ps2_unique = unique_by_key(ps2r)
    both_unique = unique_by_key(bothr)

    # row-alignment map for the optional strict-depletion flag: decoy row
    # content -> originating positive canonical key (step 3 preserves order).
    orig_key_by_content = {}
    if args.deplete_decoys_from_heldout:
        for decoy_rows in (ps1r, ps2r, bothr):
            for i, row in enumerate(decoy_rows):
                orig_a, orig_b = positives[i][0], positives[i][1]
                orig_key_by_content[(row[0], row[1])] = canonical_key(orig_a, orig_b)

    # re-derive the paralog-forbidden set for RRS verification (mirrors step 2,
    # including its mode -- see docstring on why --mode/--proteome-csv must match)
    paralog_forbidden = set()
    if not args.no_paralog_filter:
        t0 = time.time()
        pos_pairs_for_paralog = [(r[0], r[1]) for r in pos_unique]
        distinct_proteins = set()
        for a, b in pos_pairs_for_paralog:
            distinct_proteins.add(a)
            distinct_proteins.add(b)
        if args.mode == "out-of-distribution":
            proteome_pool = load_proteome_pool(args.proteome_csv, args.proteome_id_col, args.proteome_seq_col,
                                                args.proteome_filter_col, args.proteome_filter_value)
            distinct_proteins |= set(proteome_pool.values())
        close_map = build_close_match_map(distinct_proteins, k=args.paralog_k, threshold=args.paralog_threshold)
        paralog_forbidden = paralog_forbidden_keys(pos_pairs_for_paralog, close_map)
        print(f"[05] re-derived paralog-forbidden set for RRS verification (mode={args.mode}): "
              f"{len(paralog_forbidden)} pairs ({time.time()-t0:.1f}s)")

    per_replicate_stats = []

    for k in range(1, args.n_replicates + 1):
        prs = mccv_sample(pos_unique, args.frac, seed=args.seed_base_prs + k)
        rrs = mccv_sample(rand_unique, args.frac, seed=args.seed_base_rrs + k)
        ps1_ctrl = mccv_sample(ps1_unique, args.frac, seed=args.seed_base_ps1ctrl + k)
        ps2_ctrl = mccv_sample(ps2_unique, args.frac, seed=args.seed_base_ps2ctrl + k)
        both_ctrl = mccv_sample(both_unique, args.frac, seed=args.seed_base_bothctrl + k)

        write_rows_lf(f"{args.outdir}/MCCV/PRS-RRS/PRS-{tag}-{k}.csv", prs)
        write_rows_lf(f"{args.outdir}/MCCV/PRS-RRS/RRS-{tag}-{k}.csv", rrs)
        write_rows_lf(f"{args.outdir}/MCCV/random_controls/ps1random-{tag}-{k}.csv", ps1_ctrl)
        write_rows_lf(f"{args.outdir}/MCCV/random_controls/ps2random-{tag}-{k}.csv", ps2_ctrl)
        write_rows_lf(f"{args.outdir}/MCCV/random_controls/bothrandom-{tag}-{k}.csv", both_ctrl)

        # paralog re-verification on this replicate's RRS, independent re-read
        rrs_keys = set(canonical_key(r[0], r[1]) for r in read_rows(f"{args.outdir}/MCCV/PRS-RRS/RRS-{tag}-{k}.csv"))
        rrs_paralog_collisions = len(rrs_keys & paralog_forbidden) if paralog_forbidden else 0

        selected_keys = set()
        for rows in (prs, rrs, ps1_ctrl, ps2_ctrl, both_ctrl):
            for row in rows:
                selected_keys.add(canonical_key(row[0], row[1]))

        depleted = []
        removed = 0
        removed_decoy_by_origin = 0
        for row in master:
            key = canonical_key(row[0], row[1])
            if key in selected_keys:
                removed += 1
                continue
            if args.deplete_decoys_from_heldout:
                orig_key = orig_key_by_content.get((row[0], row[1]))
                if orig_key is not None and orig_key in {canonical_key(r[0], r[1]) for r in prs}:
                    removed_decoy_by_origin += 1
                    continue
            depleted.append(row)

        # final shuffle -- explicit mixing step, own seed, not inherited order
        random.Random(args.seed_base_shuffle + k).shuffle(depleted)

        depleted_path = f"{args.outdir}/MCCV/training_sets/depleted_training_set-{tag}-{k}.csv"
        write_rows_lf(depleted_path, depleted)

        # leak check -- independent re-read from what was just written
        depleted_rows_reread = read_rows(depleted_path)
        depleted_keys = set(canonical_key(r[0], r[1]) for r in depleted_rows_reread)
        leaked = len(depleted_keys & selected_keys)

        # explicit both-orientation depletion check: for every selected pair
        # (a,b), verify NEITHER the (a,b) row NOR the (b,a) row survives,
        # checked as literal ordered tuples (not just via canonical key,
        # which is what actually did the removal above -- this re-derives
        # the same guarantee a different way as an independent check).
        depleted_ordered = set((r[0], r[1]) for r in depleted_rows_reread)
        leak_forward = 0
        leak_reverse = 0
        for row in prs + rrs + ps1_ctrl + ps2_ctrl + both_ctrl:
            a, b = row[0], row[1]
            if (a, b) in depleted_ordered:
                leak_forward += 1
            if a != b and (b, a) in depleted_ordered:
                leak_reverse += 1

        per_replicate_stats.append({
            "replicate": k,
            "seeds": {
                "prs": args.seed_base_prs + k,
                "rrs": args.seed_base_rrs + k,
                "ps1_control": args.seed_base_ps1ctrl + k,
                "ps2_control": args.seed_base_ps2ctrl + k,
                "both_control": args.seed_base_bothctrl + k,
                "final_shuffle": args.seed_base_shuffle + k,
            },
            "n_prs": len(prs),
            "n_rrs": len(rrs),
            "n_ps1_control": len(ps1_ctrl),
            "n_ps2_control": len(ps2_ctrl),
            "n_both_control": len(both_ctrl),
            "n_master_rows": len(master),
            "n_removed_by_selected_keys": removed,
            "n_removed_decoys_by_origin": removed_decoy_by_origin,
            "n_depleted_training_rows": len(depleted),
            "leak_check_pairs_present_in_depleted_set": leaked,
            "orientation_leak_check": {"forward_A-B": leak_forward, "reverse_B-A": leak_reverse},
            "rrs_paralog_collisions": rrs_paralog_collisions,
        })

        print(f"[05] replicate {k}: PRS={len(prs)} RRS={len(rrs)} "
              f"ps1ctrl={len(ps1_ctrl)} ps2ctrl={len(ps2_ctrl)} bothctrl={len(both_ctrl)} "
              f"depleted_train={len(depleted)} leak={leaked} "
              f"orient_leak(fwd/rev)={leak_forward}/{leak_reverse} rrs_paralog_collisions={rrs_paralog_collisions}")

    stats = {
        "n_replicates": args.n_replicates,
        "frac": args.frac,
        "version_tag": tag,
        "deplete_decoys_from_heldout": args.deplete_decoys_from_heldout,
        "paralog_filter_enabled": not args.no_paralog_filter,
        "paralog_k": args.paralog_k,
        "paralog_threshold": args.paralog_threshold,
        "n_paralog_forbidden_pairs": len(paralog_forbidden),
        "n_unique_positive_pairs": len(pos_unique),
        "rrs_mode": args.mode,
        "n_unique_random_pairs": len(rand_unique),
        "n_ps1_random_rows": len(ps1_unique),
        "n_ps2_random_rows": len(ps2_unique),
        "n_both_random_rows": len(both_unique),
        "replicates": per_replicate_stats,
    }
    dump_json(f"{args.outdir}/MCCV/mccv_stats.json", stats)

    total_leak = sum(r["leak_check_pairs_present_in_depleted_set"] for r in per_replicate_stats)
    total_orient_leak = sum(r["orientation_leak_check"]["forward_A-B"] + r["orientation_leak_check"]["reverse_B-A"]
                             for r in per_replicate_stats)
    total_paralog_collisions = sum(r["rrs_paralog_collisions"] for r in per_replicate_stats)
    print(f"[05] done: {args.n_replicates} replicates written; "
          f"total leak-check failures = {total_leak}; "
          f"total orientation-leak failures = {total_orient_leak}; "
          f"total RRS paralog collisions = {total_paralog_collisions}")
    if total_leak or total_orient_leak:
        raise SystemExit(f"[05] LEAK DETECTED -- see MCCV/mccv_stats.json")
    if total_paralog_collisions:
        raise SystemExit(f"[05] PARALOG-SHADOW LEAK in RRS -- see MCCV/mccv_stats.json")


if __name__ == "__main__":
    main()
