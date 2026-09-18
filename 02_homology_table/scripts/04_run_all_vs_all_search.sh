#!/usr/bin/env bash
# All-vs-all protein similarity search over the unique HuRI proteins, using
# mmseqs2 (fast k-mer-indexed local alignment; a full O(n^2) Smith-Waterman
# in pure Python is not tractable for ~8,000 proteins / ~32M pairs).
#
# Requires the "homology" conda environment (created with:
#   conda create -n homology -c bioconda -c conda-forge mmseqs2 -y
# ), because mmseqs2 is not part of the base environment.
#
# Usage:
#   04_run_all_vs_all_search.sh <unique_proteins.fasta> <outdir>
#
# Input:
#   <unique_proteins.fasta>  (from 03_extract_unique_proteins.py)
# Output:
#   <outdir>/mmseqs/all_vs_all.m8   (BLAST-tab-like hits, one row per detected
#                                     local alignment between two proteins/directions)

set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: $0 <unique_proteins.fasta> <outdir>" >&2
    exit 1
fi

FASTA="$1"
OUTDIR="$2/mmseqs"
TMPDIR="$OUTDIR/tmp"
OUT_M8="$OUTDIR/all_vs_all.m8"

mkdir -p "$OUTDIR"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate homology

# -s 7.5      : maximum sensitivity, so weak/remote paralogs are not missed
# -e 10       : permissive E-value cutoff (detection threshold only) - we
#               deliberately do NOT pre-filter on identity/length here; the
#               resulting table is meant to let a %identity / alignment
#               length threshold be chosen downstream.
# --max-seqs  : generous cap on hits reported per query so promiscuous
#               protein families are not truncated
# -a 1        : compute the full backtrace (real gapped alignment), not just
#               a score/identity estimate - without this, mmseqs reports
#               mismatch/gapopen as 0 for every hit even though pident/alnlen
#               are already gap-aware.
mmseqs easy-search "$FASTA" "$FASTA" "$OUT_M8" "$TMPDIR" \
    -s 7.5 \
    -e 10 \
    --max-seqs 2000 \
    -a 1 \
    --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov" \
    --threads "$(nproc)"

echo "Wrote raw hits to $OUT_M8"
wc -l "$OUT_M8"
