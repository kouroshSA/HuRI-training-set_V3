#!/usr/bin/env python3
"""
Build an in-distribution, homology-safe set of negative (non-interacting)
protein pairs from the same protein universe as the cleaned HuRI
interactome, for use as random/negative controls alongside the positive
interaction table.

Design goals (from the request this implements):
  - Same number of negative pairs as positive pairs.
  - A negative pair (X, Y) must not itself be a true interacting pair, in
    either orientation.
  - Each protein should be used as a negative-pair endpoint about as often as
    it is used as a positive-pair endpoint (degree-preserving), so the
    negative set draws on the same proteins with the same usage frequency as
    the positive set, rather than under-using low-degree proteins or
    over-using a handful of unconstrained ones.
  - "As stringent as possible" homology safety: a candidate pair (X, Y) is
    rejected if there is ANY detectable sequence similarity (i.e. it shows up
    at all in protein_homology_table.csv, which was itself built at a
    permissive detection threshold - see the Methods.md in stage
    02_homology_table) between:
      (a) X and Y themselves,
      (b) Y and any true interacting partner of X, or
      (c) X and any true interacting partner of Y.
    Rule (b)/(c) is the direct implementation of the requested check ("C of
    A-C has no significant identity with B of A-B", generalized to *all* of
    A's true partners, not just one).

Method: degree-preserving random pairing via a stub/configuration-model
approach. A "stub list" is built where each protein appears once per positive
interaction it participates in (twice for a self-interacting/homodimer
protein, matching standard graph degree bookkeeping for self-loops). The
stub list is repeatedly shuffled and paired up; pairs that violate any rule
above are put back into a pool and re-shuffled in the next round, until the
full stub list is consumed (or no further progress can be made).

Input:
  huri_interactome_uniprotID_seqs_cleaned_deduped_standard.csv (stage 01 output)
  protein_homology_table.csv (stage 02 output)
  proteins/unique_proteins.fasta (stage 02 output; sequence lookup)
Output:
  huri_negative_pairs_indistribution.csv
  logs/06_build_negative_pairs.log

Usage:
  python3 06_build_negative_pairs.py --positives <standardized.csv> \
      --homology <protein_homology_table.csv> --fasta <unique_proteins.fasta> \
      --outdir <dir>
"""

import argparse
import sys
import random
from collections import defaultdict
from pathlib import Path

import pandas as pd

RANDOM_SEED = 20260918
MAX_ROUNDS = 2000


def load_fasta(path: Path) -> dict:
    seqs = {}
    pid, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if pid is not None:
                    seqs[pid] = "".join(chunks)
                pid, chunks = line[1:].strip(), []
            else:
                chunks.append(line)
        if pid is not None:
            seqs[pid] = "".join(chunks)
    return seqs


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--positives", required=True, help="path to huri_interactome_uniprotID_seqs_cleaned_deduped_standard.csv")
    ap.add_argument("--homology", required=True, help="path to protein_homology_table.csv")
    ap.add_argument("--fasta", required=True, help="path to unique_proteins.fasta")
    ap.add_argument("--outdir", required=True, help="directory to write outputs into")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    output_csv = outdir / "huri_negative_pairs_indistribution.csv"
    log_file = outdir / "logs" / "06_build_negative_pairs.log"

    rng = random.Random(RANDOM_SEED)
    lines = []

    def log(msg: str) -> None:
        print(msg)
        lines.append(msg)

    log(f"Reading {args.positives}")
    pos = pd.read_csv(args.positives, dtype=str, keep_default_na=False)
    log(f"Positive pairs: {len(pos)}")

    log(f"Reading {args.homology}")
    hom = pd.read_csv(args.homology, dtype={"Protein_A": str, "Protein_B": str})
    log(f"Homology pairs (any detectable similarity): {len(hom)}")

    log(f"Reading {args.fasta}")
    seqs = load_fasta(Path(args.fasta))
    log(f"Sequences loaded: {len(seqs)}")

    # ---- positive-graph structures ----
    pos_set = set()          # frozenset({id1,id2}) for non-self positive pairs
    partners = defaultdict(set)  # true interacting partners of each protein
    stub_count = defaultdict(int)
    for a, b in zip(pos["Protein_ID_1"], pos["Protein_ID_2"]):
        if a == b:
            stub_count[a] += 2
            partners[a].add(a)
        else:
            stub_count[a] += 1
            stub_count[b] += 1
            partners[a].add(b)
            partners[b].add(a)
            pos_set.add(frozenset((a, b)))

    total_stubs = sum(stub_count.values())
    log(f"Unique proteins (degree > 0): {len(stub_count)}; total stubs: {total_stubs} "
        f"(target negative pairs = total_stubs / 2 = {total_stubs // 2})")
    target_n = total_stubs // 2

    # ---- homology adjacency (any row = "significant" by construction; see
    #      stage 02_homology_table Methods.md - the table only contains
    #      detected hits at a permissive but real detection threshold) ----
    homolog_of = defaultdict(set)
    for a, b in zip(hom["Protein_A"], hom["Protein_B"]):
        homolog_of[a].add(b)
        homolog_of[b].add(a)

    def violates_homology_rule(x: str, y: str) -> bool:
        # (a) X and Y themselves detectably similar
        if y in homolog_of.get(x, ()):
            return True
        # (b) Y resembles a true partner of X
        if not partners[x].isdisjoint(homolog_of.get(y, ())):
            return True
        # (c) X resembles a true partner of Y
        if not partners[y].isdisjoint(homolog_of.get(x, ())):
            return True
        return False

    # ---- build stub list ----
    stubs = []
    for pid, n in stub_count.items():
        stubs.extend([pid] * n)
    log(f"Stub list length: {len(stubs)}")

    neg_set = set()
    round_no = 0
    stuck_rounds = 0
    prev_len = None
    while stubs and len(neg_set) < target_n and round_no < MAX_ROUNDS:
        round_no += 1
        rng.shuffle(stubs)
        leftover = []
        n_pairs_this_round = len(stubs) // 2
        for i in range(n_pairs_this_round):
            x, y = stubs[2 * i], stubs[2 * i + 1]
            if x == y:
                leftover.extend([x, y])
                continue
            key = frozenset((x, y))
            if key in pos_set or key in neg_set:
                leftover.extend([x, y])
                continue
            if violates_homology_rule(x, y):
                leftover.extend([x, y])
                continue
            neg_set.add(key)
        if len(stubs) % 2 == 1:
            leftover.append(stubs[-1])
        stubs = leftover

        if prev_len is not None and len(stubs) >= prev_len:
            stuck_rounds += 1
        else:
            stuck_rounds = 0
        prev_len = len(stubs)

        if round_no <= 10 or round_no % 25 == 0:
            log(f"Round {round_no}: accepted so far = {len(neg_set)}, "
                f"remaining unmatched stubs = {len(stubs)}")

        if stuck_rounds >= 15:
            log(f"No progress for {stuck_rounds} consecutive rounds "
                f"({len(stubs)} stubs remain unmatched) - stopping.")
            break

    log(f"Phase 1 (degree-preserving stub matching) finished after {round_no} rounds: "
        f"{len(neg_set)} negative pairs accepted, {len(stubs)} stubs left unmatched, "
        f"target was {target_n}")

    if stubs:
        from collections import Counter
        leftover_counts = Counter(stubs)
        log(f"Unmatched stub composition (protein: leftover stub count), "
            f"top 10: {leftover_counts.most_common(10)}")

    # ---- Phase 2: top-up ----
    # Phase 1 only re-pairs a protein against OTHER proteins that are also
    # still short of their target degree, which can deadlock a handful of
    # hub/heavily-paralogous proteins against each other even though they
    # each individually have plenty of valid (non-positive, homology-safe)
    # candidates in the FULL protein pool. Phase 2 resolves each remaining
    # unmatched stub against the full pool: it first tries other still-needy
    # proteins (keeping degree-preservation as tight as possible), and only
    # falls back to an already-satisfied protein (accepting a small, logged
    # "overshoot" of +1 to that protein's negative-degree) if that fails.
    all_proteins = list(stub_count.keys())
    remaining_capacity = defaultdict(int)
    for pid in stubs:
        remaining_capacity[pid] += 1
    overshoot_count = 0
    overshoot_uses = defaultdict(int)

    if stubs:
        log(f"Phase 2 (full-pool top-up) starting: {len(stubs)} stubs "
            f"({len(set(stubs))} unique proteins) still need a partner")

    pending = list(stubs)
    rng.shuffle(pending)
    i = 0
    max_attempts_per_stub = 500
    while i < len(pending) and len(neg_set) < target_n:
        x = pending[i]
        if remaining_capacity[x] <= 0:
            i += 1
            continue
        placed = False
        # First choice: another still-needy protein.
        needy_candidates = [p for p in remaining_capacity if remaining_capacity[p] > 0 and p != x]
        rng.shuffle(needy_candidates)
        for y in needy_candidates[:max_attempts_per_stub]:
            key = frozenset((x, y))
            if key in pos_set or key in neg_set:
                continue
            if violates_homology_rule(x, y):
                continue
            neg_set.add(key)
            remaining_capacity[x] -= 1
            remaining_capacity[y] -= 1
            placed = True
            break
        if not placed:
            # Fall back to the full pool; accept a degree overshoot on the partner.
            for _ in range(max_attempts_per_stub):
                y = rng.choice(all_proteins)
                if y == x:
                    continue
                key = frozenset((x, y))
                if key in pos_set or key in neg_set:
                    continue
                if violates_homology_rule(x, y):
                    continue
                neg_set.add(key)
                remaining_capacity[x] -= 1
                overshoot_count += 1
                overshoot_uses[y] += 1
                placed = True
                break
        if not placed:
            log(f"WARNING: could not find any valid partner for {x} after "
                f"{max_attempts_per_stub} attempts in both needy and full pools")
        i += 1

    log(f"Phase 2 finished: {len(neg_set)} negative pairs total "
        f"(target {target_n}); {overshoot_count} pairs used a full-pool "
        f"overshoot partner (protein used beyond its exact positive degree)")
    if overshoot_uses:
        top_overshoot = sorted(overshoot_uses.items(), key=lambda kv: -kv[1])[:10]
        log(f"Top overshoot partners (protein: extra times used beyond its "
            f"positive degree): {top_overshoot}")

    # ---- verification pass: re-check every accepted pair against all rules ----
    violations = 0
    for pair in neg_set:
        x, y = tuple(pair)
        if x == y or pair in pos_set or violates_homology_rule(x, y):
            violations += 1
    log(f"Verification: {violations} of {len(neg_set)} accepted pairs violate a rule "
        f"on re-check (expect 0)")

    # ---- write output table, same schema as the positive interactome file ----
    rows = []
    for pair in neg_set:
        x, y = tuple(pair)
        rows.append((x, seqs[x], y, seqs[y]))
    out = pd.DataFrame(rows, columns=["Protein_ID_1", "Sequence_1", "Protein_ID_2", "Sequence_2"])
    out = out.sort_values(["Protein_ID_1", "Protein_ID_2"]).reset_index(drop=True)
    outdir.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_csv, index=False)
    log(f"Wrote {len(out)} negative pairs to {output_csv}")

    # ---- degree-preservation diagnostic ----
    neg_degree = defaultdict(int)
    for pair in neg_set:
        x, y = tuple(pair)
        neg_degree[x] += 1
        neg_degree[y] += 1
    diffs = [neg_degree.get(pid, 0) - stub_count[pid] for pid in stub_count]
    import statistics
    log(f"Degree-match check (negative usage minus positive degree per protein): "
        f"mean diff = {statistics.mean(diffs):.4f}, "
        f"proteins with 0 difference = {sum(1 for d in diffs if d == 0)}/{len(diffs)}, "
        f"max |diff| = {max(abs(d) for d in diffs)}")

    log_file.parent.mkdir(parents=True, exist_ok=True)
    log_file.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
