"""
FastQC on raw (uploads/) and trimmed (data/) reads.

Markers come from config.py (lower-case, e.g. 16s/its/18s) so folder names match
what the uploader and priming step actually create - fixes the previous
case-sensitivity bug where a hard-coded "ITS" folder was never found on Linux.
"""
import os
import subprocess
import sys

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

sys.path.insert(0, os.path.abspath(os.path.join(BASE_DIR, "..")))
try:
    import config as cfg
    UPLOAD_FOLDER = str(cfg.UPLOAD_DIR)
    TRIMMED_FOLDER = str(cfg.DATA_DIR)
    MARKERS = list(cfg.MARKERS)
    FASTQ_EXTS = tuple(cfg.FASTQ_EXTS)
except Exception:
    UPLOAD_FOLDER = os.path.join(BASE_DIR, "../uploads")
    TRIMMED_FOLDER = os.path.join(BASE_DIR, "../data")
    MARKERS = ["16s", "its", "18s"]
    FASTQ_EXTS = (".fastq", ".fq", ".fastq.gz", ".fq.gz")

FASTQC_PRE = os.path.join(BASE_DIR, "../results/fastqc_pretrim")
FASTQC_POST = os.path.join(BASE_DIR, "../results/fastqc_posttrim")
os.makedirs(FASTQC_PRE, exist_ok=True)
os.makedirs(FASTQC_POST, exist_ok=True)


def run_fastqc(input_file, outdir):
    if not os.path.exists(input_file):
        print(f"❌ ERROR: File not found - {input_file}")
        return
    print(f"🚀 Running FastQC on {input_file}...")
    try:
        subprocess.run(["fastqc", input_file, "--outdir", outdir], check=True)
        print(f"✅ FastQC complete: {os.path.basename(input_file)}")
    except FileNotFoundError:
        print("❌ 'fastqc' not found on PATH. Install it (conda install -c bioconda fastqc).")
    except subprocess.CalledProcessError as e:
        print(f"❌ FastQC failed for {input_file}: {e}")


def get_fastq_files(folder):
    found = []
    for root, _, files in os.walk(folder):
        for file in files:
            if file.lower().endswith(FASTQ_EXTS):
                found.append(os.path.join(root, file))
    return found


if __name__ == "__main__":
    print("🔹 Running FastQC on pre-trimmed files (uploads)...")
    for marker in MARKERS:
        # try the exact marker name and a couple of case variants, for safety
        candidates = {marker, marker.lower(), marker.upper()}
        path = next((os.path.join(UPLOAD_FOLDER, m)
                     for m in candidates if os.path.isdir(os.path.join(UPLOAD_FOLDER, m))), None)
        if path is None:
            print(f"⏩ Skipping missing folder for {marker.upper()}")
            continue
        raw_files = get_fastq_files(path)
        if not raw_files:
            print(f"⚠️ No FASTQ files in {marker.upper()}")
            continue
        for f in raw_files:
            run_fastqc(f, FASTQC_PRE)

    print("🔹 Running FastQC on post-trimmed files (data)...")
    trimmed_files = get_fastq_files(TRIMMED_FOLDER)
    if not trimmed_files:
        print("⚠️ No trimmed files found in data/")
    else:
        for f in trimmed_files:
            run_fastqc(f, FASTQC_POST)
