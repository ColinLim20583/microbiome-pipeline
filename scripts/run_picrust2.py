"""
PICRUSt2 functional prediction — clean, local, single-command version.

This replaces the previous approaches, which did not work reliably:
  * the Galaxy-EU version (bioblend) embedded a hard-coded API key and depended
    on a remote server + exact tool names;
  * the "stratified" version called `sudo` to build a swapfile (hangs with no
    terminal) and relied on the vendored picrust2/ source tree.

Here we call the OFFICIAL `picrust2_pipeline.py` command, which runs sequence
placement, hidden-state prediction, metagenome inference and pathway abundance
in one step. Install it in its own conda env (it conflicts with QIIME2):

    conda create -n picrust2 -c bioconda -c conda-forge picrust2
    conda activate picrust2
    python scripts/run_picrust2.py            # or run this whole script in that env

Inputs (produced by the export step) per marker under data/exported/<marker>/:
    rep_seqs.fna            (representative sequences)
    feature-table.biom      (feature/ASV table)

By default only 16S is processed (PICRUSt2 is validated for prokaryotic 16S).
Set PEANUT_PICRUST2_MARKERS="16s,18s" to override.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE_DIR))
try:
    import config as cfg
    EXPORTED_DIR = cfg.EXPORTED_DIR
except Exception:
    EXPORTED_DIR = BASE_DIR / "data" / "exported"

# PICRUSt2 is only biologically valid for 16S by default.
MARKERS = os.environ.get("PEANUT_PICRUST2_MARKERS", "16s").split(",")
THREADS = os.environ.get("PEANUT_PICRUST2_THREADS", "2")

PICRUST2_CMD = "picrust2_pipeline.py"


def have_picrust2() -> bool:
    return shutil.which(PICRUST2_CMD) is not None


def run_marker(marker: str) -> None:
    folder = EXPORTED_DIR / marker
    rep_fna = folder / "rep_seqs.fna"
    # the export step names it rep_seqs.fna; fall back to dna-sequences.fasta
    if not rep_fna.exists():
        alt = folder / "dna-sequences.fasta"
        rep_fna = alt if alt.exists() else rep_fna
    biom = folder / "feature-table.biom"
    out_dir = folder / "picrust2_output"

    if not rep_fna.exists() or not biom.exists():
        print(f"⏩ Skipping {marker.upper()}: need {rep_fna.name} and {biom.name} "
              f"in {folder} (run the export step first).")
        return

    if out_dir.exists():
        print(f"🧹 Removing previous output: {out_dir}")
        shutil.rmtree(out_dir)

    cmd = [
        PICRUST2_CMD,
        "-s", str(rep_fna),
        "-i", str(biom),
        "-o", str(out_dir),
        "-p", THREADS,
        "--stratified",
        "--verbose",
    ]
    print(f"🔬 Running PICRUSt2 for {marker.upper()}:\n   $ {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
        print(f"✅ PICRUSt2 finished for {marker.upper()} -> {out_dir}")
        print("   Key outputs: KO_metagenome_out/, EC_metagenome_out/, "
              "pathways_out/ (all *.tsv.gz).")
    except subprocess.CalledProcessError as e:
        print(f"❌ PICRUSt2 failed for {marker.upper()} (exit {e.returncode}).")


def main() -> None:
    if not have_picrust2():
        print("❌ 'picrust2_pipeline.py' is not on your PATH.")
        print("   PICRUSt2 must be installed in an active conda env:")
        print("     conda create -n picrust2 -c bioconda -c conda-forge picrust2")
        print("     conda activate picrust2")
        print("   Then re-run this step from that environment.")
        # Exit 0 so the pipeline runner treats this as 'skipped', not 'crashed'.
        return

    for marker in [m.strip() for m in MARKERS if m.strip()]:
        run_marker(marker)
    print("\n🏁 PICRUSt2 step complete (for available markers).")


if __name__ == "__main__":
    main()
