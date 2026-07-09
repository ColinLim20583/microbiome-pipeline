import subprocess
import sys
from pathlib import Path

print("🚀 Starting per-region DADA2 denoising with QZV support...")

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR.parent / "data"

# DADA2 truncation lengths are dataset-dependent and should be tuned per marker
# from the demux.qzv quality profile. They live in config.py (and are editable
# from the Streamlit UI) rather than being hard-coded here.
sys.path.insert(0, str(BASE_DIR.parent))
try:
    import config as cfg
    _D = cfg.DADA2
    MARKERS = cfg.MARKERS
    DADA2_PARAMS = {m: _D[m] for m in MARKERS if m in _D}
    MAX_EE_F = str(_D.get("max_ee_f", 2))
    MAX_EE_R = str(_D.get("max_ee_r", 2))
    MIN_OVERLAP = str(_D.get("min_overlap", 20))
    N_THREADS = str(_D.get("n_threads", 1))
    N_READS_LEARN = str(_D.get("n_reads_learn", 100000))
except Exception:
    DADA2_PARAMS = {
        "16s": {"trunc_len_f": 220, "trunc_len_r": 220, "trim_left_f": 0, "trim_left_r": 0},
        "its": {"trunc_len_f": 200, "trunc_len_r": 220, "trim_left_f": 0, "trim_left_r": 0},
        "18s": {"trunc_len_f": 220, "trunc_len_r": 220, "trim_left_f": 0, "trim_left_r": 0},
    }
    MAX_EE_F = MAX_EE_R = "2"
    MIN_OVERLAP = "20"
    N_THREADS = "1"
    N_READS_LEARN = "100000"

for region, params in DADA2_PARAMS.items():
    demux_file = DATA_DIR / f"{region}_demux.qza"
    if not demux_file.exists():
        print(f"⚠️ Skipping {region.upper()}: Missing demux file at {demux_file}")
        continue

    print(f"\n🔬 Processing {region.upper()}...")

    table_qza = DATA_DIR / f"{region}_table.qza"
    rep_seqs_qza = DATA_DIR / f"{region}_rep_seqs.qza"
    stats_qza = DATA_DIR / f"{region}_denoising_stats.qza"

    try:
        subprocess.run([
            "qiime", "dada2", "denoise-paired",
            "--i-demultiplexed-seqs", str(demux_file),
            "--p-trim-left-f", str(params["trim_left_f"]),
            "--p-trim-left-r", str(params["trim_left_r"]),
            "--p-trunc-len-f", str(params["trunc_len_f"]),
            "--p-trunc-len-r", str(params["trunc_len_r"]),
            "--p-max-ee-f", MAX_EE_F,
            "--p-max-ee-r", MAX_EE_R,
            "--p-min-overlap", MIN_OVERLAP,
            "--p-n-threads", N_THREADS,
            "--p-n-reads-learn", N_READS_LEARN,
            "--o-table", str(table_qza),
            "--o-representative-sequences", str(rep_seqs_qza),
            "--o-denoising-stats", str(stats_qza),
            "--verbose"
        ], check=True)

        print(f"✅ DADA2 completed for {region.upper()}")
    except subprocess.CalledProcessError:
        print(f"❌ DADA2 failed for {region.upper()}")
        continue

    # QZVs
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
