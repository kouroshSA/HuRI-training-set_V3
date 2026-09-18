#!/usr/bin/env bash
# Run the vendored ppi-v3-dataset-pipeline (pinned at commit
# d41e5083ad233be011beef95d32990f75b231bd1, see ../pipeline_src/VENDORED_COMMIT.txt)
# end to end against the standardized HuRI positive interactome, producing
# the canonical V3 pipeline output tree (positives, RRS, the three
# random-substitution decoy sets, master training mix, and 10 MCCV
# replicates with PRS/RRS/decoy holdouts).
#
# Usage:
#   01_run_v3_pipeline.sh <positives.csv> <outdir> <pipeline_src_dir>
#
# RRS mode: in-distribution (default) -- real-protein negatives drawn from
# the same positive pool via a degree-matched configuration model,
# paralog-filtered by the pipeline's own k-mer Jaccard proxy (see
# pipeline_src/README.md "Paralog filter"). This is a different, lighter-
# weight filter than the mmseqs2 alignment-based protein_homology_table.csv
# built in stage 02_homology_table and used for huri_negative_pairs_indistribution.csv
# in stage 03_negative_pairs; see this stage's Methods.md "Relationship to
# the earlier custom negative-pair set" for why the pipeline's own RRS was
# used here instead of substituting that one in.

set -euo pipefail

if [ $# -ne 3 ]; then
    echo "Usage: $0 <positives.csv> <outdir> <pipeline_src_dir>" >&2
    exit 1
fi

INPUT="$1"
OUTDIR="$2"
PIPELINE_SRC="$3"
PIPELINE="$PIPELINE_SRC/scripts/run_pipeline.py"

echo "Input positive set: $INPUT"
echo "Pipeline (vendored): $PIPELINE"
echo "Output dir:          $OUTDIR"
echo

python3 "$PIPELINE" \
    --input "$INPUT" \
    --seq1-col Sequence_1 --seq2-col Sequence_2 \
    --outdir "$OUTDIR" \
    --n-replicates 10 \
    --frac 0.10 \
    --version-tag V3 \
    --rrs-mode in-distribution

echo
echo "Pipeline run complete. See $OUTDIR/provenance.md and $OUTDIR/manifest.json."
