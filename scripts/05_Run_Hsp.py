#!/usr/bin/env python3
"""
05_run_hsp.py — robust HSP runner for PICRUSt2

• Uses your existing placed tree at functional_prediction/16s/placed_seqs.tre
• If that file is missing or looks like a *labelled* tree, tries to rebuild an
  unlabelled placed tree from the EPA-NG JPLACE using GAPPA.
• Verifies Rscript + R packages are available (HSP calls castor in R).
• Executes HSP via your local picrust2 clone (module form) with the right flags.

Run from anywhere (paths are resolved relative to this script's location):

    conda activate picrust2
    python scripts/05_run_hsp.py               # default: trait=16S, processes=1

Optional flags:
    --processes 4                             # number of CPU processes for HSP
    --trait 16S|KO|EC|...                     # trait table to predict (default 16S)
    --placed-tree /path/to/tree.newick        # override tree path
    --regen-tree                              # force regenerating tree with GAPPA

Outputs (in functional_prediction/16s/):
    marker_predicted.tsv.gz   (and confidence intervals)
"""
from __future__ import annotations
import argparse
import shutil
import subprocess
import sys
import os
from pathlib import Path

# --------------------- CLI ---------------------
parser = argparse.ArgumentParser(description="Run PICRUSt2 HSP with safety checks.")
parser.add_argument("--trait", default="16S",
                    choices=["16S","BIGG","CAZY","EC","GENE_NAMES","GO","KO","PFAM","COG","TIGRFAM","PHENO"],
                    help="Trait table to predict (default: 16S).")
parser.add_argument("--processes", "-p", type=int, default=1,
                    help="Number of CPU processes (default: 1).")
parser.add_argument("--placed-tree", default=None,
                    help="Path to an EXISTING placed tree to use (NEWICK). Overrides autodetect.")
parser.add_argument("--regen-tree", action="store_true",
                    help="Force regenerate an unlabelled placed tree from JPLACE using GAPPA.")
args = parser.parse_args()

# --------------------- Paths (relative to repo root) ---------------------
BASE = Path(__file__).resolve().parents[1]      # repo root (one up from scripts/)
PICRUST2_DIR = BASE / "picrust2"
OUT_DIR = BASE / "functional_prediction/16s"
INTERMEDIATE = OUT_DIR / "intermediate"
PLACED_TRE = OUT_DIR / "placed_seqs.tre"
JPLACE = INTERMEDIATE / "epa_out" / "epa_result.jplace"
REF_TRE = BASE / "picrust2/picrust2/default_files/bacteria/bac_ref/bac_ref.tre"

# --------------------- Helpers ---------------------
def sh(cmd: list[str], check=True, env=None):
    print("➤", " ".join(map(str, cmd)))
    return subprocess.run(cmd, check=check, env=env)

def have(cmd: str) -> bool:
    return shutil.which(cmd) is not None

# --------------------- Checks ---------------------
OUT_DIR.mkdir(parents=True, exist_ok=True)

# 1) Ensure Rscript exists (HSP needs R + castor)
if not have("Rscript"):
    sys.exit("⛔ Rscript not found in PATH. Activate the picrust2 env and install R: \n"
             "   conda install -c conda-forge r-base=4.3 r-castor r-ape r-data.table r-optparse")

# 2) Determine placed tree
if args.placed_tree:
    src = Path(args.placed_tree)
    if not src.exists():
        sys.exit(f"⛔ --placed-tree not found: {src}")
    if src.resolve() != PLACED_TRE.resolve():
        shutil.copy(src, PLACED_TRE)
        print(f"📎 Using provided placed tree → {PLACED_TRE}")

# Rebuild if requested or missing
def regenerate_with_gappa():
    if not have("gappa"):
        sys.exit("⛔ Need gappa to regenerate the placed tree. Install it or run in your qiime2 env.")
    if not JPLACE.exists():
        sys.exit(f"⛔ JPLACE not found: {JPLACE} — run 02_place_with_epa.py first.")
    if not REF_TRE.exists():
        sys.exit(f"⛔ Reference tree not found: {REF_TRE}")
    # Try the common subcommand; some builds use 'to-tree' under 'examine'
    try:
        sh(["gappa","examine","to-tree",
            "--jplace-path", str(JPLACE),
            "--ref-tree", str(REF_TRE),
            "--out-dir", str(OUT_DIR),
            "--file-prefix","placed_seqs",
            "--best-hit","--allow-file-overwriting"])
    except subprocess.CalledProcessError:
        # Fallback: some versions expect 'to-tree' as separate; try --help to confirm
        sys.exit("⛔ gappa examine to-tree failed. Run 'gappa examine to-tree --help' to confirm availability.")
    # Normalize filename to .tre
    produced = None
    for name in ("placed_seqs.tree","placed_seqs.newick"):
        cand = OUT_DIR / name
        if cand.exists() and cand.stat().st_size > 0:
            produced = cand
            break
    if not produced:
        sys.exit("⛔ GAPPA did not produce a placed tree file.")
    shutil.move(str(produced), str(PLACED_TRE))
    print(f"✅ Wrote {PLACED_TRE}")

if args.regen_tree or not PLACED_TRE.exists():
    print("🌳 Building unlabelled placed tree …")
    regenerate_with_gappa()
else:
    print(f"✅ Using existing placed tree: {PLACED_TRE}")

# 3) Run HSP via local picrust2 clone (module form)
py_env = os.environ.copy()
py_env["PYTHONPATH"] = str(PICRUST2_DIR)  # allow local clone import
py_env["PYTHONPATH"] = str(PICRUST2_DIR)

cmd = [
    sys.executable,
    str(PICRUST2_DIR / "scripts" / "hsp.py"),
    "-t", str(PLACED_TRE),
    "-o", str(OUT_DIR),
    "-i", args.trait,
    "-p", str(args.processes),
]

sh(cmd, env=py_env)

# 4) Confirm outputs
marker = OUT_DIR / "marker_predicted.tsv.gz"
if marker.exists() and marker.stat().st_size > 0:
    print(f"✅ HSP complete: {marker}")
else:
    sys.exit("⛔ HSP finished but marker_predicted.tsv.gz not found — check logs above.")
