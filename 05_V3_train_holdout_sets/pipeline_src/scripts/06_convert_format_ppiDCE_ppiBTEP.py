#!/usr/bin/env python3
"""
Step 6 -- convert every canonical pipeline output into the ppiDCE / ppiBTEP
format (the two models share one format; no marker tokens, sequences and
pairing copied verbatim):

    training : SEQ1,SEQ2,label      label in {0,1}   (3 columns)
    PRS / RRS / controls: SEQ1,SEQ2                  (2 columns, no label)

Headerless, LF line endings. Mirrors the existing project convention seen in
ppiDCE_ppiBTEP_train_PRS-RRS_/.

Writes the whole tree under <outdir>/formats/ppiDCE_ppiBTEP/, same relative
layout as the canonical files (see common.iter_canonical_files).
"""
import argparse
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from common import read_rows, write_rows_lf, dump_json, iter_canonical_files


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True, help="pipeline output root (contains 01_positives.csv etc.)")
    args = ap.parse_args()

    dest_root = f"{args.outdir}/formats/ppiDCE_ppiBTEP"
    converted = []
    for rel_path, role in iter_canonical_files(args.outdir):
        rows = read_rows(f"{args.outdir}/{rel_path}")
        if role == "train":
            out_rows = [(r[0], r[1], r[2]) for r in rows]
        else:
            out_rows = [(r[0], r[1]) for r in rows]
        out_path = f"{dest_root}/{rel_path}"
        write_rows_lf(out_path, out_rows)
        converted.append({"file": rel_path, "role": role, "rows": len(out_rows)})

    dump_json(f"{dest_root}/conversion_manifest.json", {"format": "ppiDCE_ppiBTEP", "files": converted})
    print(f"[06] converted {len(converted)} files -> {dest_root}")


if __name__ == "__main__":
    main()
