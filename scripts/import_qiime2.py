import os
import subprocess

BASE_DIR = os.path.abspath(os.path.dirname(__file__))
DATA_DIR = os.path.join(BASE_DIR, "../data")
MODES = ["16s", "its", "18s"]

def run_import(mode):
    manifest_path = os.path.join(DATA_DIR, f"{mode}_manifest.tsv")
    qza_path = os.path.join(DATA_DIR, f"{mode}_demux.qza")
    qzv_path = os.path.join(DATA_DIR, f"{mode}_demux.qzv")

    if not os.path.exists(manifest_path):
        print(f"⚠️  Skipping {mode.upper()}: Manifest not found at {manifest_path}")
        return

    print(f"\n📂 Importing {mode.upper()} using manifest: {manifest_path}")
    import_cmd = [
        "qiime", "tools", "import",
        "--type", "SampleData[PairedEndSequencesWithQuality]",
        "--input-path", manifest_path,
        "--output-path", qza_path,
        "--input-format", "PairedEndFastqManifestPhred33V2"
    ]

    result = subprocess.run(import_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"❌ Failed to import {mode.upper()}")
        print(result.stderr)
        return
    else:
        print(f"✅ Imported {mode.upper()} demux to: {qza_path}")

    # Generate QZV summary
    print(f"📊 Summarizing {mode.upper()} demux to: {qzv_path}")
    summary_cmd = [
        "qiime", "demux", "summarize",
        "--i-data", qza_path,
        "--o-visualization", qzv_path
    ]
    subprocess.run(summary_cmd, check=True)
    print(f"✅ Summary written: {qzv_path}")

def main():
    print("📦 Importing FASTQ data using QIIME 2...")
    for mode in MODES:
        run_import(mode)
    print("\n🏁 FASTQ import completed for available regions.")

if __name__ == "__main__":
    main()
