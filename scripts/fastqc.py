import os
import subprocess

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

UPLOAD_FOLDER = os.path.join(BASE_DIR, "../uploads")
TRIMMED_FOLDER = os.path.join(BASE_DIR, "../data")

FASTQC_PRE = os.path.join(BASE_DIR, "../results/fastqc_pretrim")
FASTQC_POST = os.path.join(BASE_DIR, "../results/fastqc_posttrim")

# Create result directories
os.makedirs(FASTQC_PRE, exist_ok=True)
os.makedirs(FASTQC_POST, exist_ok=True)

def run_fastqc(input_file, outdir):
    if not os.path.exists(input_file):
        print(f"❌ ERROR: File not found - {input_file}")
        return
    print(f"🚀 Running FastQC on {input_file}...")
    cmd = f"fastqc \"{input_file}\" --outdir=\"{outdir}\""
    try:
        subprocess.run(cmd, shell=True, check=True)
        print(f"✅ FastQC complete: {os.path.basename(input_file)}")
    except subprocess.CalledProcessError as e:
        print(f"❌ FastQC failed for {input_file}: {e}")

def get_fastq_files(folder):
    fastq_files = []
    for root, _, files in os.walk(folder):
        for file in files:
            if file.endswith((".fastq", ".fq", ".fastq.gz")):
                fastq_files.append(os.path.join(root, file))
    return fastq_files

if __name__ == "__main__":
    print("🔹 Running FastQC on pre-trimmed files (uploads)...")

    for region in ["16s", "ITS", "18s"]:
        path = os.path.join(UPLOAD_FOLDER, region)
        if not os.path.exists(path):
            print(f"⏩ Skipping missing folder: {path}")
            continue
        raw_files = get_fastq_files(path)
        if not raw_files:
            print(f"⚠️ No FASTQ files in {region}")
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
