"""
PICRUSt2 functional prediction — one clean, working command.

This single script replaces the entire hand-built PICRUSt2 chain that used to
live here (01_align_seqs, 02_place_with_epa, 03_gappa_to_newick,
04_make_placed_tree, 05_Run_Hsp, picrust2_balance, run_picrust2_stratified,
and the Galaxy version). All of those manually reconstructed — and fragile-ly —
what the official `picrust2_pipeline.py` does in a single validated call:
sequence placement → hidden-state prediction → metagenome inference →
pathway abundance.

Install PICRUSt2 in its OWN conda env (it conflicts with QIIME2):

    conda create -n picrust2 -c bioconda -c conda-forge picrust2
    conda activate picrust2
    python scripts/run_picrust2.py

Inputs (produced by the export step) per marker in data/exported/<marker>/:
    rep_seqs.fna         representative sequences
    feature-table.biom   feature/ASV table

Outputs in data/exported/<marker>/picrust2_output/:
    KO_metagenome_out/  EC_metagenome_out/  pathways_out/   (all *.tsv.gz)
    plus a human-readable top_pathways_summary.csv + .png written by this script.

PICRUSt2 is only biologically valid for prokaryotic 16S, so only 16S runs by
default. Override with  PEANUT_PICRUST2_MARKERS="16s,18s".
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

MARKERS = os.environ.get("PEANUT_PICRUST2_MARKERS", "16s").split(",")
THREADS = os.environ.get("PEANUT_PICRUST2_THREADS", "2")
TOP_N = int(os.environ.get("PEANUT_PICRUST2_TOP_N", "25"))
PICRUST2_CMD = "picrust2_pipeline.py"


def have_picrust2() -> bool:
    return shutil.which(PICRUST2_CMD) is not None


def summarize_pathways(out_dir: Path, marker: str) -> None:
    """Write a human-readable top-pathways summary (CSV + bar plot)."""
    path_abun = out_dir / "pathways_out" / "path_abun_unstrat.tsv.gz"
    if not path_abun.exists():
        print(f"   (no pathways_out/path_abun_unstrat.tsv.gz to summarize for {marker.upper()})")
        return
    try:
        import pandas as pd
        df = pd.read_csv(path_abun, sep="\t", index_col=0)
        totals = df.sum(axis=1).sort_values(ascending=False)
        top = totals.head(TOP_N).rename("total_predicted_abundance").reset_index()
        top.columns = ["pathway", "total_predicted_abundance"]
        csv_fp = out_dir / "top_pathways_summary.csv"
        top.to_csv(csv_fp, index=False)
        print(f"   📄 wrote {csv_fp.name} (top {TOP_N} predicted MetaCyc pathways)")

        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(9, max(4, 0.3 * len(top))))
            ax.barh(top["pathway"][::-1], top["total_predicted_abundance"][::-1], color="#2c7bb6")
            ax.set_xlabel("Total predicted abundance")
            ax.set_title(f"Top {TOP_N} predicted pathways ({marker.upper()})")
            fig.tight_layout()
            png_fp = out_dir / "top_pathways_summary.png"
            fig.savefig(png_fp, dpi=200)
            plt.close(fig)
            print(f"   📈 wrote {png_fp.name}")
        except Exception as e:
            print(f"   (plot skipped: {e})")
    except Exception as e:
        print(f"   ⚠️ could not summarize pathways for {marker.upper()}: {e}")


def run_marker(marker: str) -> None:
    folder = EXPORTED_DIR / marker
    rep_fna = folder / "rep_seqs.fna"
    if not rep_fna.exists():
        alt = folder / "dna-sequences.fasta"
        rep_fna = alt if alt.exists() else rep_fna
    biom = folder / "feature-table.biom"
    out_dir = folder / "picrust2_output"

    if not rep_fna.exists() or not biom.exists():
        print(f"⏩ Skipping {marker.upper()}: need {rep_fna.name} and feature-table.biom "
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
        summarize_pathways(out_dir, marker)
    except subprocess.CalledProcessError as e:
        print(f"❌ PICRUSt2 failed for {marker.upper()} (exit {e.returncode}).")


def main() -> None:
    if not have_picrust2():
        print("❌ 'picrust2_pipeline.py' is not on your PATH.")
        print("   PICRUSt2 must be installed in an active conda env:")
        print("     conda create -n picrust2 -c bioconda -c conda-forge picrust2")
        print("     conda activate picrust2")
        print("   Then re-run:  python scripts/run_picrust2.py")
        return  # exit 0 -> pipeline treats this as 'skipped', not 'crashed'

    for marker in [m.strip() for m in MARKERS if m.strip()]:
        run_marker(marker)
    print("\n🏁 PICRUSt2 step complete (for available markers).")


if __name__ == "__main__":
    main()
