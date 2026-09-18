#!/usr/bin/env python3
"""
Step 7 -- convert every canonical pipeline output into ppiGPLM native format,
the "_alt1" token variant (short single-character partner separators):

    training : $,SEQ1,!,SEQ2,<label>     label in {<0>,<1>}   (5 fields)
    PRS / RRS / controls: $,SEQ1,!,SEQ2,<   (trailing open marker, no label
                                              -- the model completes the
                                              prompt; leaving a label in would
                                              hand it the answer)

Headerless, LF line endings. Mirrors ppiGPLM_train_PRS-RRS_alt1/ in the
existing project (originally `<ps1>`/`<ps2>` tokens, rewritten to `$`/`!`;
this script writes the alt1 tokens directly).

Writes the whole tree under <outdir>/formats/ppiGPLM/, same relative layout
as the canonical files (see common.iter_canonical_files).
"""
import argparse
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from common import read_rows, write_rows_lf, dump_json, iter_canonical_files

PS1_TOKEN = "$"
PS2_TOKEN = "!"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True, help="pipeline output root (contains 01_positives.csv etc.)")
    args = ap.parse_args()

    dest_root = f"{args.outdir}/formats/ppiGPLM"
    converted = []
    for rel_path, role in iter_canonical_files(args.outdir):
        rows = read_rows(f"{args.outdir}/{rel_path}")
        out_rows = []
        for r in rows:
            seq1, seq2 = r[0], r[1]
            if role == "train":
                label_tok = f"<{r[2]}>"
                out_rows.append((PS1_TOKEN, seq1, PS2_TOKEN, seq2, label_tok))
            else:
                out_rows.append((PS1_TOKEN, seq1, PS2_TOKEN, seq2, "<"))
        out_path = f"{dest_root}/{rel_path}"
        write_rows_lf(out_path, out_rows)
        converted.append({"file": rel_path, "role": role, "rows": len(out_rows)})

    dump_json(f"{dest_root}/conversion_manifest.json",
              {"format": "ppiGPLM_alt1", "ps1_token": PS1_TOKEN, "ps2_token": PS2_TOKEN, "files": converted})
    print(f"[07] converted {len(converted)} files -> {dest_root}")


if __name__ == "__main__":
    main()
