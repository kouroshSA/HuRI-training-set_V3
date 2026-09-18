"""
Shared utilities for the PPI V3 training-set pipeline.

Conventions enforced project-wide (see README.md "Conventions"):
  - All CSV output is headerless, comma-separated, LF line endings ONLY.
    (csv.writer defaults to CRLF -- that has bitten this project before:
    a stray \\r is an out-of-vocabulary character for the char-level
    ppiGPLM tokenizer. Always go through write_rows_lf() / open_lf().)
  - Random sequences are drawn uniformly from the 20 standard amino acids
    and must be globally unique within a build (tracked via UniqueSeqPool).
  - Pair identity for de-duplication / depletion is orientation-independent:
    canonical_key(a, b) = frozenset({a, b}).
  - Every seed used anywhere in the pipeline is an explicit CLI argument
    with a documented default -- nothing is left to an unseeded RNG.
"""

import argparse
import csv
import hashlib
import json
import os
import random

STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"
STANDARD_AA_SET = set(STANDARD_AA)


# --------------------------------------------------------------------------
# I/O
# --------------------------------------------------------------------------

def open_lf(path, mode="w"):
    """Open a file guaranteeing LF-only line endings on write."""
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)
    return open(path, mode, newline="\n")


def write_rows_lf(path, rows):
    """Write an iterable of row-tuples/lists as headerless CSV, LF endings."""
    with open_lf(path) as f:
        w = csv.writer(f, lineterminator="\n")
        for row in rows:
            w.writerow(row)


def read_rows(path):
    with open(path, newline="") as f:
        return list(csv.reader(f))


def check_lf(path):
    """Return True if the file contains no bare \\r anywhere."""
    with open(path, "rb") as f:
        return b"\r" not in f.read()


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------
# Sequence handling
# --------------------------------------------------------------------------

def clean_sequence(seq):
    """Strip whitespace and a single trailing stop-codon '*', uppercase."""
    s = seq.strip().upper()
    if s.endswith("*"):
        s = s[:-1]
    return s


def is_standard_aa(seq):
    return len(seq) > 0 and set(seq) <= STANDARD_AA_SET


def canonical_key(seq1, seq2):
    """Orientation-independent identity for a pair."""
    return frozenset((seq1, seq2))


class UniqueSeqPool:
    """
    Generates random same-length amino-acid sequences that are (a) globally
    unique against every real sequence ever passed to register_real(), and
    (b) globally unique against every random sequence this pool has already
    emitted -- across the *whole* pipeline run, not just one file. Matches
    the collision-free guarantee documented throughout the MED4/Synecho
    lineage of this pipeline (0 collisions, enforced not assumed).
    """

    def __init__(self, seed):
        self._rng = random.Random(seed)
        self._used = set()
        self._collisions = 0

    def register_real(self, seq):
        self._used.add(seq)

    def draw(self, length):
        while True:
            s = "".join(self._rng.choice(STANDARD_AA) for _ in range(length))
            if s not in self._used:
                self._used.add(s)
                return s
            self._collisions += 1

    @property
    def collisions(self):
        return self._collisions

    @property
    def n_generated(self):
        return len(self._used)


# --------------------------------------------------------------------------
# Degree-matched negative pairing (configuration model)
# --------------------------------------------------------------------------

def configuration_model_negative_pairs(positive_pairs, n_target, seed, max_rounds=400,
                                        extra_forbidden=None):
    """
    Build n_target unique unordered negative pairs over the same node set as
    `positive_pairs`, with a degree sequence matching the positive graph
    (configuration model / stub-shuffling), rejecting self-pairs, duplicate
    negative pairs, any pair that is itself a positive pair, and (if given)
    any pair whose canonical key is in `extra_forbidden` -- used to also
    reject paralog-shadow pairs, see build_close_match_map /
    paralog_forbidden_keys below.

    positive_pairs: iterable of (seq1, seq2) tuples (one row per pair; if a
        pair should count once, de-duplicate by canonical_key before calling).
    n_target: number of unique unordered negative pairs to produce. If it
        exceeds what the stub multiset can support without exhausting valid
        partners, rounds are repeated with re-shuffled stubs until reached
        or max_rounds is hit (raises RuntimeError on failure, never silently
        returns short -- callers should not have to guess whether they got
        what they asked for).
    seed: RNG seed (deterministic).
    extra_forbidden: optional set of canonical_key()s to also reject, on top
        of self-pairs / duplicates / positive pairs.

    Returns: (list of (seq1, seq2) negative pairs, achieved Pearson r between
        positive and negative degree over the shared node set).
    """
    rng = random.Random(seed)
    extra_forbidden = extra_forbidden or set()

    pos_keys = set(canonical_key(a, b) for a, b in positive_pairs)
    degree = {}
    for a, b in positive_pairs:
        degree[a] = degree.get(a, 0) + 1
        degree[b] = degree.get(b, 0) + 1

    # Build stub list: node repeated `degree` times. Sort first -- Python
    # salts str-hash per process, so iterating an unsorted set/dict would
    # make this non-deterministic across runs despite the fixed seed
    # (this exact bug was found and fixed in the Synecho build).
    nodes_sorted = sorted(degree.keys())

    neg_pairs = set()
    rounds = 0
    while len(neg_pairs) < n_target and rounds < max_rounds:
        rounds += 1
        stubs = []
        for node in nodes_sorted:
            stubs.extend([node] * degree[node])
        rng.shuffle(stubs)
        for i in range(0, len(stubs) - 1, 2):
            a, b = stubs[i], stubs[i + 1]
            if a == b:
                continue
            key = canonical_key(a, b)
            if key in pos_keys or key in neg_pairs or key in extra_forbidden:
                continue
            neg_pairs.add(key)
            if len(neg_pairs) >= n_target:
                break

    if len(neg_pairs) < n_target:
        raise RuntimeError(
            f"configuration_model_negative_pairs: only reached {len(neg_pairs)}/"
            f"{n_target} unique negative pairs after {max_rounds} rounds "
            f"(paralog filtering may be too strict for this node set -- "
            f"try a lower --paralog-threshold, a higher --paralog-k, or "
            f"--no-paralog-filter)"
        )

    result = [tuple(sorted(k)) for k in neg_pairs]
    # deterministic emission order
    result.sort()
    r = _degree_pearson_r(positive_pairs, result)
    return result, r


def _degree_pearson_r(positive_pairs, negative_pairs):
    pos_deg, neg_deg = {}, {}
    for a, b in positive_pairs:
        pos_deg[a] = pos_deg.get(a, 0) + 1
        pos_deg[b] = pos_deg.get(b, 0) + 1
    for a, b in negative_pairs:
        neg_deg[a] = neg_deg.get(a, 0) + 1
        neg_deg[b] = neg_deg.get(b, 0) + 1
    nodes = sorted(pos_deg.keys())
    xs = [pos_deg[n] for n in nodes]
    ys = [neg_deg.get(n, 0) for n in nodes]
    n = len(nodes)
    if n < 2:
        return float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0:
        return float("nan")
    return cov / (vx ** 0.5 * vy ** 0.5)


# --------------------------------------------------------------------------
# Paralog / close-match detection (alignment-free, k-mer Jaccard)
# --------------------------------------------------------------------------
#
# Purpose: a "random" real-protein negative pair (X, Y) is a bad negative if
# X is a near-duplicate (paralog) of some protein A, Y is a near-duplicate of
# some protein B, and (A, B) is a genuine positive pair -- (X, Y) could very
# plausibly be a true (but unassayed) interaction via conserved paralogous
# binding, not a real negative. This is a fast, dependency-free proxy for
# that check (k-mer Jaccard similarity, not a real alignment/BLAST E-value --
# if you have actual homology search results, prefer those; this exists so
# the pipeline needs no external tools).

def kmer_set(seq, k):
    if len(seq) < k:
        return {seq}
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def jaccard(set_a, set_b):
    if not set_a and not set_b:
        return 0.0
    inter = len(set_a & set_b)
    union = len(set_a | set_b)
    return inter / union if union else 0.0


def build_close_match_map(sequences, k=4, threshold=0.5, length_ratio_min=0.5):
    """
    For a collection of distinct protein sequences, return {seq: set(seqs)}
    mapping every sequence to the set of sequences (including itself) whose
    k-mer Jaccard similarity to it is >= threshold. O(n^2) pairwise, with a
    length-ratio pre-filter (pairs whose lengths differ by more than
    length_ratio_min are skipped without computing Jaccard -- they cannot
    plausibly clear a sane threshold, and this is what keeps a few thousand
    proteins tractable in pure Python).
    """
    seqs = sorted(set(sequences))  # sorted: deterministic despite str-hash salting
    kmers = {s: kmer_set(s, k) for s in seqs}
    lengths = {s: len(s) for s in seqs}
    close = {s: {s} for s in seqs}
    n = len(seqs)
    for i in range(n):
        si = seqs[i]
        li = lengths[si]
        ki = kmers[si]
        for j in range(i + 1, n):
            sj = seqs[j]
            lj = lengths[sj]
            ratio = lj / li if li >= lj else li / lj
            if ratio < length_ratio_min:
                continue
            if jaccard(ki, kmers[sj]) >= threshold:
                close[si].add(sj)
                close[sj].add(si)
    return close


def paralog_forbidden_keys(positive_pairs, close_map):
    """
    For every positive pair (A, B), and every (X in close_map[A], Y in
    close_map[B]) or (X in close_map[B], Y in close_map[A]), add
    canonical_key(X, Y) to the forbidden set. Includes the positive pairs
    themselves (X=A, Y=B is always in its own close set).
    """
    forbidden = set()
    for a, b in positive_pairs:
        close_a = close_map.get(a, {a})
        close_b = close_map.get(b, {b})
        for x in close_a:
            for y in close_b:
                forbidden.add(canonical_key(x, y))
        for x in close_b:
            for y in close_a:
                forbidden.add(canonical_key(x, y))
    return forbidden


# --------------------------------------------------------------------------
# Monte-Carlo fractional sampling
# --------------------------------------------------------------------------

def mccv_sample(items, frac, seed):
    """
    Sample round(frac * len(items)) items without replacement, seeded.
    Independent per call -- draws across different `seed`s (e.g. different
    replicates) may overlap; that is the defining property of Monte-Carlo
    cross-validation as opposed to k-fold, and is intentional here.
    """
    rng = random.Random(seed)
    items = list(items)
    n = round(frac * len(items))
    return rng.sample(items, n)


# --------------------------------------------------------------------------
# manifest helpers
# --------------------------------------------------------------------------

def count_rows(path):
    with open(path, newline="") as f:
        return sum(1 for _ in f)


def load_json(path):
    with open(path) as f:
        return json.load(f)


def dump_json(path, obj):
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)
    with open(path, "w", newline="\n") as f:
        json.dump(obj, f, indent=2, sort_keys=True)
        f.write("\n")


def common_arg_parser(description):
    p = argparse.ArgumentParser(description=description)
    return p


# --------------------------------------------------------------------------
# canonical file inventory (used by the two format converters, 06 and 07,
# so both walk exactly the same set of files with the same train/eval
# classification -- adding a new canonical file only needs to change this
# one place)
# --------------------------------------------------------------------------

import glob


def iter_canonical_files(outdir):
    """
    Yield (relative_path, role) for every canonical 3-field pipeline output
    under `outdir`, role in {"train", "eval"}. "train" files keep their
    label column when converted; "eval" files (PRS/RRS/controls -- bare
    prompts the model completes) have the label column stripped.
    """
    train_globs = [
        "01_positives.csv",
        "Parent_sequences/random_pairs.csv",
        "Parent_sequences/ps1_random.csv",
        "Parent_sequences/ps2_random.csv",
        "Parent_sequences/both_random.csv",
        "master_training_mix.csv",
        "MCCV/training_sets/depleted_training_set-*.csv",
    ]
    eval_globs = [
        "MCCV/PRS-RRS/PRS-*.csv",
        "MCCV/PRS-RRS/RRS-*.csv",
        "MCCV/random_controls/*.csv",
    ]
    for pattern in train_globs:
        for path in sorted(glob.glob(os.path.join(outdir, pattern))):
            yield os.path.relpath(path, outdir), "train"
    for pattern in eval_globs:
        for path in sorted(glob.glob(os.path.join(outdir, pattern))):
            yield os.path.relpath(path, outdir), "eval"
