import os
import subprocess
from pathlib import Path

# ---- CONFIG ----
MARKERS = {
    "16s": "bacteria/bac_ref",
    "its": "fungi_ITS",
    "its_amf": "fungi_ITS",
    "18s_amf": "fungi_18S",
}

PICRUST2_DIR = Path("../picrust2")
DEFAULTS_DIR = PICRUST2_DIR / "picrust2" / "default_files"
PLACE_SEQS_SCRIPT = PICRUST2_DIR / "scripts/place_seqs.py"

env = os.environ.copy()
env["PYTHONPATH"] = str(PICRUST2_DIR)


def run_place_seqs(marker, ref_subdir):
    print(f"\n🔬 Processing marker: {marker}")
    DATA_DIR = Path(f"../data/exported/{marker}")
    REP_SEQS = DATA_DIR / "rep_seqs.fna"
    OUT_DIR = Path(f"../functional_prediction/{marker}")
    INTERMEDIATE = OUT_DIR / "intermediate"
    OUT_TREE = OUT_DIR / "placed_seqs.tre"

    REF_DIR = DEFAULTS_DIR / ref_subdir
    REF_FASTA = list(REF_DIR.glob("*.fna")) + list(REF_DIR.glob("*.fna.gz"))
    HMM_FILE = list(REF_DIR.glob("*.hmm"))

    # Input checks
    if not REP_SEQS.exists():
        print(f"⚠️ Skipping {marker}: {REP_SEQS} not found.")
        return
    if not REF_FASTA or not HMM_FILE:
        print(f"❌ Skipping {marker}: missing reference files in {REF_DIR}")
        return

    INTERMEDIATE.mkdir(parents=True, exist_ok=True)
    print(f"📁 Created intermediate folder: {INTERMEDIATE}")

    cmd = [
        "python", str(PLACE_SEQS_SCRIPT),
        "-s", str(REP_SEQS),
        "-r", str(REF_DIR),  # Picrust2 script expects directory containing all ref files
        "-o", str(OUT_TREE),
        "--intermediate", str(INTERMEDIATE),
        "--verbose"
    ]

    try:
        subprocess.run(cmd, check=True, env=env)
        print(f"✅ place_seqs.py completed for {marker}.")
    except subprocess.CalledProcessError as e:
        print(f"❌ place_seqs.py failed for {marker}. Continuing to next...")
        print(f"↪️  Command: {' '.join(e.cmd)}")
        print(f"↪️  Exit Code: {e.returncode}")


# ---- MAIN ----
if __name__ == "__main__":
    for marker, ref_path in MARKERS.items():
        run_place_seqs(marker, ref_path)
