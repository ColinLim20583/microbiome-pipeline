import os
import subprocess
from pathlib import Path

print("🚀 Starting per-region DADA2 denoising with QZV support...")

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR.parent / "data"

# Region-specific DADA2 trimming/truncating settings
DADA2_PARAMS = {
    "16s": {"trunc_len_f": 220, "trunc_len_r": 200, "trim_left_f": 0, "trim_left_r": 0},
    "its": {"trunc_len_f": 200, "trunc_len_r": 220, "trim_left_f": 0, "trim_left_r": 0},
    "18s": {"trunc_len_f": 220, "trunc_len_r": 220, "trim_left_f": 20, "trim_left_r": 20},
}

for region, params in DADA2_PARAMS.items():
    demux_file = DATA_DIR / f"{region}_demux.qza"
    if not demux_file.exists():
        print(f"⚠️ Skipping {region.upper()}: Missing demux file at {demux_file}")
        continue

    print(f"\n🔬 Processing {region.upper()}...")

    # Output paths
    table_qza = DATA_DIR / f"{region}_table.qza"
    rep_seqs_qza = DATA_DIR / f"{region}_rep_seqs.qza"
    stats_qza = DATA_DIR / f"{region}_denoising_stats.qza"

    # Run DADA2
    try:
        subprocess.run([
            "qiime", "dada2", "denoise-paired",
            "--i-demultiplexed-seqs", str(demux_file),
            "--p-trim-left-f", str(params["trim_left_f"]),
            "--p-trim-left-r", str(params["trim_left_r"]),
            "--p-trunc-len-f", str(params["trunc_len_f"]),
            "--p-trunc-len-r", str(params["trunc_len_r"]),
            "--p-max-ee-f", "2",
            "--p-max-ee-r", "2",
            "--p-n-threads", "2",  # WSL-safe
            "--p-n-reads-learn", "100000",  # ✅ correct flag
            "--o-table", str(table_qza),
            "--o-representative-sequences", str(rep_seqs_qza),
            "--o-denoising-stats", str(stats_qza),
            "--verbose"
        ], check=True)

        print(f"✅ DADA2 completed for {region.upper()}")
    except subprocess.CalledProcessError:
        print(f"❌ DADA2 failed for {region.upper()}")
        continue

    # Generate visualizations
    try:
        subprocess.run([
            "qiime", "feature-table", "summarize",
            "--i-table", str(table_qza),
            "--o-visualization", str(DATA_DIR / f"{region}_table.qzv")
        ], check=True)

        subprocess.run([
            "qiime", "feature-table", "tabulate-seqs",
            "--i-data", str(rep_seqs_qza),
            "--o-visualization", str(DATA_DIR / f"{region}_rep_seqs.qzv")
        ], check=True)

        subprocess.run([
            "qiime", "metadata", "tabulate",
            "--m-input-file", str(stats_qza),
            "--o-visualization", str(DATA_DIR / f"{region}_denoising_stats.qzv")
        ], check=True)

        print(f"📊 QZV visualizations created for {region.upper()}")
    except subprocess.CalledProcessError:
        print(f"❌ Failed to generate QZV files for {region.upper()}")

print("\n🏁 All regions processed (if data available).")
