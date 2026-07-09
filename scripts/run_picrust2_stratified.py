import os
import subprocess
from pathlib import Path
import shutil

# ---- CONFIG ----
DATA_DIR = Path("../data/exported/16s")
PICRUST2_DIR = Path("../picrust2")
PICRUST2_SCRIPTS = PICRUST2_DIR / "scripts"
REP_SEQS = DATA_DIR / "rep_seqs.fna"
FEATURE_TABLE = DATA_DIR / "feature-table.biom"
OUT_DIR = Path("functional_prediction/16s")

# ---- CLEANUP ----
if OUT_DIR.exists():
    print("🧹 Cleaning up previous output directory...")
    shutil.rmtree(OUT_DIR)
OUT_DIR.mkdir(parents=True)

# Clean EPA output if it exists (to avoid EPA-NG error)
epa_out = OUT_DIR / "intermediate" / "epa_out"
if epa_out.exists():
    print("🧹 Removing old EPA-NG output...")
    shutil.rmtree(epa_out)

# ---- ENVIRONMENT ----
env = os.environ.copy()
env["PYTHONPATH"] = str(PICRUST2_DIR)
env["EPA_NG_EXTRA_ARGS"] = "--chunk-size 3000 -T 2"  # SAFE SETTINGS for EPA-NG

# ---- HELPER FUNCTION ----
def run_step(name, cmd, env=None):
    print(name)
    try:
        result = subprocess.run(
            cmd,
            check=True,
            env=env,
            text=True,
            capture_output=True
        )
    except subprocess.CalledProcessError as e:
        print("❌ Error during subprocess execution.")
        print("Command:", ' '.join(e.cmd))
        print("Exit code:", e.returncode)
        print("\n📤 STDOUT:\n", e.stdout)
        print("\n📥 STDERR:\n", e.stderr)
        raise

def ensure_swapfile(swap_path="/swapfile2", size_gb=16):
    print("🧠 Ensuring swap is active...")
    result = subprocess.run(["swapon", "--show"], capture_output=True, text=True)
    if swap_path not in result.stdout:
        print(f"⚠️  No swap at {swap_path}. Creating {size_gb}GB swapfile...")
        subprocess.run(["sudo", "swapoff", swap_path], stderr=subprocess.DEVNULL)
        subprocess.run(["sudo", "rm", "-f", swap_path])
        subprocess.run(["sudo", "fallocate", "-l", f"{size_gb}G", swap_path], check=True)
        subprocess.run(["sudo", "chmod", "600", swap_path], check=True)
        subprocess.run(["sudo", "mkswap", swap_path], check=True)
        subprocess.run(["sudo", "swapon", swap_path], check=True)
        print("✅ Swap enabled.")
    else:
        print("✅ Swap already active.")


# Call this early in your script
ensure_swapfile()

# ---- STEP 1: Place sequences ----
run_step("🌲 Step 1: Placing sequences...", [
    "python",
    str(PICRUST2_SCRIPTS / "place_seqs.py"),
    "-s", str(REP_SEQS),
    "-o", str(OUT_DIR / "placed_seqs.tre"),
    "-p", "1",
    "--intermediate", str(OUT_DIR / "intermediate"),
    "--verbose"
], env=env)

# ---- STEP 2: HSP ----
run_step("🧠 Step 2: Hidden-state prediction...", [
    "python", str(PICRUST2_SCRIPTS / "hsp.py"),
    "-i", str(OUT_DIR / "placed_seqs.tre"),
    "-o", str(OUT_DIR),
    "--threads", "1",
    "--redo"
], env=env)

# ---- STEP 3: Predict metagenome ----
run_step("📊 Step 3: Predicting metagenome...", [
    "python", str(PICRUST2_SCRIPTS / "metagenome_pipeline.py"),
    "-i", str(FEATURE_TABLE),
    "-m", str(OUT_DIR / "marker_predicted.tsv.gz"),
    "-f", str(OUT_DIR / "functional_predicted.tsv.gz"),
    "-o", str(OUT_DIR / "pred_metagenome_unstrat.tsv.gz"),
    "--threads", "1",
    "--redo"
], env=env)

# ---- STEP 4: Predict pathways ----
run_step("🔬 Step 4: Predicting pathways (stratified)...", [
    "python", str(PICRUST2_SCRIPTS / "pathway_pipeline.py"),
    "-i", str(OUT_DIR / "pred_metagenome_unstrat.tsv.gz"),
    "-o", str(OUT_DIR),
    "--stratified",
    "--no-regroup",
    "--threads", "1",
    "--redo"
], env=env)

print("✅ Stratified metagenome prediction complete!")