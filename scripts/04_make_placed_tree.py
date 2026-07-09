#!/usr/bin/env python3
from pathlib import Path
import subprocess, shutil, sys

# Paths relative to repo root; tweak if yours differ
OUT_DIR = Path("../functional_prediction/16s")
INTERMEDIATE = OUT_DIR / "intermediate"
PLACED_TRE = OUT_DIR / "placed_seqs.tre"

# Candidates we might already have from step 03 (GAPPA)
CANDIDATES = [
    INTERMEDIATE / "placed_seqslabelled_tree.newick",
    INTERMEDIATE / "placed_seqs_tree.newick",
    INTERMEDIATE / "placed_seqs.newick",
]

# EPA jplace + reference tree (used if we need to rebuild)
JPLACE = INTERMEDIATE / "epa_out" / "epa_result.jplace"
REF_TRE = Path("../picrust2/picrust2/default_files/bacteria/bac_ref/bac_ref.tre")

def shell(cmd):
    print("➤", " ".join(cmd))
    subprocess.run(cmd, check=True)

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1) If a placed tree already exists, keep it
    if PLACED_TRE.exists() and PLACED_TRE.stat().st_size > 0:
        print(f"✅ Already present: {PLACED_TRE}")
        return

    # 2) Try to copy an existing newick from step 03
    for cand in CANDIDATES:
        if cand.exists() and cand.stat().st_size > 0:
            print(f"📎 Using existing: {cand} → {PLACED_TRE}")
            shutil.copy(cand, PLACED_TRE)
            print("✅ Done.")
            return

    # 3) Rebuild from JPLACE with GAPPA
    if not JPLACE.exists():
        sys.exit(f"⛔ Cannot find {JPLACE}. Run 02_place_with_epa.py first.")
    if not REF_TRE.exists():
        sys.exit(f"⛔ Cannot find reference tree: {REF_TRE}")

    tmp_out = OUT_DIR  # write into the final folder
    # Use gappa to make a placed tree (best-hit placements)
    shell([
        "gappa", "examine", "to-tree",
        "--jplace-path", str(JPLACE),
        "--ref-tree",    str(REF_TRE),
        "--out-dir",     str(tmp_out),
        "--file-prefix", "placed_seqs",
        "--best-hit",
        "--allow-file-overwriting",
    ])

    # GAPPA writes either .tree or .newick; normalize to .tre
    produced = None
    for name in ("placed_seqs.tree", "placed_seqs.newick"):
        cand = OUT_DIR / name
        if cand.exists() and cand.stat().st_size > 0:
            produced = cand
            break
    if not produced:
        sys.exit("⛔ GAPPA did not produce a placed tree.")

    shutil.move(str(produced), str(PLACED_TRE))
    print(f"✅ Wrote {PLACED_TRE}")

if __name__ == "__main__":
    main()
