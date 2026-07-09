import subprocess
from pathlib import Path

# ---- CONFIG ----
DATA_DIR = Path("../data/exported/16s")  # Adjust if run from 'scripts/' directory
PICRUST2_DIR = Path("../picrust2")
REP_SEQS = DATA_DIR / "rep_seqs.fna"
OUT_DIR = Path("../functional_prediction/16s")
INTERMEDIATE = OUT_DIR / "intermediate"

# Input files (produced by 01_align_seqs.py)
QUERY_ALN = INTERMEDIATE / "query_align.stockholm"
STUDY_ALN = INTERMEDIATE / "study_seqs_hmmalign.fasta"
REF_MSA = INTERMEDIATE / "ref_seqs_hmmalign.fasta"

# Reference tree and model
REF_DIR = PICRUST2_DIR / "picrust2" / "default_files" / "bacteria" / "bac_ref"
TREE = REF_DIR / "bac_ref.tre"
MODEL = REF_DIR / "bac_ref.model"

# Output directory for EPA-NG
EPA_OUT = INTERMEDIATE / "epa_out"
EPA_OUT.mkdir(parents=True, exist_ok=True)

# ---- Optional: Ensure swap (for systems with low RAM) ----
def ensure_swapfile(swap_path="/swapfile2", size_gb=16):
    try:
        print("🧠 Ensuring swap is active (best-effort)...")
        result = subprocess.run(["swapon", "--show"], capture_output=True, text=True)
        if swap_path not in result.stdout:
            print(f"ℹ️  No swap at {swap_path}. Attempting to create {size_gb}GB swapfile...")
            subprocess.run(["sudo", "swapoff", swap_path], stderr=subprocess.DEVNULL)
            subprocess.run(["sudo", "rm", "-f", swap_path])
            # The next calls may fail on hosts without sudo — that's OK.
            subprocess.run(["sudo", "fallocate", "-l", f"{size_gb}G", swap_path], check=False)
            subprocess.run(["sudo", "chmod", "600", swap_path], check=False)
            subprocess.run(["sudo", "mkswap", swap_path], check=False)
            subprocess.run(["sudo", "swapon", swap_path], check=False)
        print("✅ Continuing regardless of swap status.")
    except Exception as e:
        print(f"⚠️  Swap setup failed (non-fatal): {e}. Continuing.")

# ---- Run EPA-NG placement ----
def run_epa_ng():
    print("📌 Running EPA-NG for phylogenetic placement...")

    cmd = [
        "epa-ng",
        "--tree", str(TREE),
        "--ref-msa", str(REF_MSA),
        "--query", str(STUDY_ALN),
        "--chunk-size", "3000",
        "-T", "2",
        "-m", str(MODEL),
        "-w", str(EPA_OUT),
        "--filter-acc-lwr", "0.99",
        "--filter-max", "100",
        "--redo"
    ]

    try:
        subprocess.run(cmd, check=True)
        print("✅ EPA-NG placement complete.")
    except subprocess.CalledProcessError as e:
        print("❌ EPA-NG failed.")
        print("Command:", " ".join(e.cmd))
        print("Exit Code:", e.returncode)
        raise

# ---- MAIN ----
if __name__ == "__main__":
    ensure_swapfile()
    run_epa_ng()
