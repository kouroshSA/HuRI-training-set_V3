#!/usr/bin/env bash
# Assemble the standardized positive (interacting) and negative (in-distribution
# random) pair tables into a master/training staging folder, as two SEPARATE
# files - they are deliberately NOT merged/labeled/shuffled into one combined
# table at this stage. That combination step (adding a label column, shuffling,
# and splitting into train/validation/test) is a later, separate step (see
# stage 05_V3_train_holdout_sets).
#
# Usage:
#   07_assemble_master_and_train.sh <positives.csv> <negatives.csv> <destdir>
#
# Source files (already built by stages 01 and 03):
#   huri_interactome_uniprotID_seqs_cleaned_deduped_standard.csv  (positive pairs)
#   huri_negative_pairs_indistribution.csv                        (negative pairs)
#
# Each file is copied (not moved, so the original pipeline outputs stay in
# place for provenance/re-run purposes) and verified byte-for-byte via SHA-256.

set -euo pipefail

if [ $# -ne 3 ]; then
    echo "Usage: $0 <positives.csv> <negatives.csv> <destdir>" >&2
    exit 1
fi

POS_SRC="$1"
NEG_SRC="$2"
DEST="$3"
LOG="$DEST/logs/07_assemble_master_and_train.log"

mkdir -p "$DEST" "$DEST/logs"

{
  echo "Assembling master/train staging folder: $DEST"
  echo "Run date (UTC): $(date -u '+%Y-%m-%d %H:%M:%S')"
  echo

  for SRC in "$POS_SRC" "$NEG_SRC"; do
    BASENAME="$(basename "$SRC")"
    DST="$DEST/$BASENAME"
    cp "$SRC" "$DST"

    SRC_SUM="$(sha256sum "$SRC" | cut -d' ' -f1)"
    DST_SUM="$(sha256sum "$DST" | cut -d' ' -f1)"
    SRC_ROWS="$(wc -l < "$SRC")"
    DST_ROWS="$(wc -l < "$DST")"

    echo "File: $BASENAME"
    echo "  source: $SRC"
    echo "  dest:   $DST"
    echo "  source sha256: $SRC_SUM"
    echo "  dest   sha256: $DST_SUM"
    echo "  rows (incl. header) source/dest: $SRC_ROWS / $DST_ROWS"
    if [ "$SRC_SUM" != "$DST_SUM" ]; then
      echo "  ERROR: checksum mismatch after copy!" >&2
      exit 1
    fi
    if [ "$SRC_ROWS" != "$DST_ROWS" ]; then
      echo "  ERROR: row count mismatch after copy!" >&2
      exit 1
    fi
    echo "  OK: copy verified byte-identical"
    echo
  done

  echo "Done. Both tables are present in $DEST as separate, unmerged files."
} | tee "$LOG"
