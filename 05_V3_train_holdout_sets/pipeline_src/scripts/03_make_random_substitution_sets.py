#!/usr/bin/env python3
"""
Step 3 -- the three "hard decoy" parent sets: partner-matched negatives built
by substituting one or both partners of every positive pair with a random,
same-length sequence. These stop a model from learning "protein X is always
an interactor" instead of judging the pair.

  ps1_random.csv    first partner (SEQ1) replaced by a random sequence
  ps2_random.csv    second partner (SEQ2) replaced by a random sequence
  both_random.csv   both partners replaced (the "floor" control -- no real
                     protein content survives, only length/composition)

Each of the three sets is generated with its own seed (--seed-ps1,
--seed-ps2, --seed-both -- defaults 43/44/45) for maximum diversity of the
random proteins used across sets. Random sequences are drawn uniformly from
the 20 standard amino acids and are enforced globally unique -- against every
real sequence in the positive set AND against every random sequence
generated anywhere else in this step (across all three files), never just
within one file.

Saved separately, one file per variant, in <outdir>/Parent_sequences/ --
these are the full-scale ("parent") sets, one substituted row per original
positive row (including both (A,B) and (B,A) orientation rows, since
01_positives.csv now carries both -- this script iterates positives as-is
and needs no special-casing to inherit that symmetry), before any
Monte-Carlo subsampling.

Writes:
  <outdir>/Parent_sequences/ps1_random.csv
  <outdir>/Parent_sequences/ps2_random.csv
  <outdir>/Parent_sequences/both_random.csv
  <outdir>/Parent_sequences/random_substitution_stats.json
"""
import argparse
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from common import read_rows, write_rows_lf, dump_json, UniqueSeqPool


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--positives", required=True, help="path to 01_positives.csv")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed-ps1", type=int, default=43)
    ap.add_argument("--seed-ps2", type=int, default=44)
    ap.add_argument("--seed-both", type=int, default=45)
    args = ap.parse_args()

    pos_rows = read_rows(args.positives)
    pairs = [(r[0], r[1]) for r in pos_rows]

    real_sequences = set()
    for a, b in pairs:
        real_sequences.add(a)
        real_sequences.add(b)

    generated_so_far = set()  # carried across all three variants -> global uniqueness

    def new_pool(seed):
        pool = UniqueSeqPool(seed)
        for s in real_sequences:
            pool.register_real(s)
        for s in generated_so_far:
            pool.register_real(s)
        return pool

    # --- ps1_random ---
    pool = new_pool(args.seed_ps1)
    ps1_rows = []
    for a, b in pairs:
        fake = pool.draw(len(a))
        generated_so_far.add(fake)
        ps1_rows.append((fake, b, 0))
    write_rows_lf(f"{args.outdir}/Parent_sequences/ps1_random.csv", ps1_rows)
    ps1_collisions = pool.collisions

    # --- ps2_random ---
    pool = new_pool(args.seed_ps2)
    ps2_rows = []
    for a, b in pairs:
        fake = pool.draw(len(b))
        generated_so_far.add(fake)
        ps2_rows.append((a, fake, 0))
    write_rows_lf(f"{args.outdir}/Parent_sequences/ps2_random.csv", ps2_rows)
    ps2_collisions = pool.collisions

    # --- both_random ---
    pool = new_pool(args.seed_both)
    both_rows = []
    for a, b in pairs:
        fake_a = pool.draw(len(a))
        generated_so_far.add(fake_a)
        fake_b = pool.draw(len(b))
        generated_so_far.add(fake_b)
        both_rows.append((fake_a, fake_b, 0))
    write_rows_lf(f"{args.outdir}/Parent_sequences/both_random.csv", both_rows)
    both_collisions = pool.collisions

    stats = {
        "source_positives": args.positives,
        "n_rows_per_set": len(pairs),
        "seed_ps1_random": args.seed_ps1,
        "seed_ps2_random": args.seed_ps2,
        "seed_both_random": args.seed_both,
        "n_random_sequences_generated_total": len(generated_so_far),
        "collision_retries": {
            "ps1_random": ps1_collisions,
            "ps2_random": ps2_collisions,
            "both_random": both_collisions,
        },
        "global_uniqueness": "enforced: no generated sequence equals any real sequence "
                              "or any other generated sequence across all three files",
    }
    dump_json(f"{args.outdir}/Parent_sequences/random_substitution_stats.json", stats)

    print(f"[03] wrote ps1_random.csv, ps2_random.csv, both_random.csv "
          f"({len(pairs)} rows each; {len(generated_so_far)} unique random sequences total, "
          f"0 collisions with real or prior-generated sequences)")


if __name__ == "__main__":
    main()
