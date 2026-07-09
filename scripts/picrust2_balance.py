import os
import subprocess
from pathlib import Path

# ---- CONFIG ----
DATA_DIR = Path("../data/exported/16s")
PICRUST2_DIR = Path("../picrust2")
PICRUST2_SCRIPTS = PICRUST2_DIR / "scripts"
FEATURE_TABLE = DATA_DIR / "feature-table.biom"
OUT_DIR = Path("functional_prediction/16s")
PLACED_TREE = OUT_DIR / "placed_seqs.tre"

# ---- ENVIRONMENT ----
env = os.environ.copy()
env["PYTHONPATH"] = str(PICRUST2_DIR)

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

# ---- STEP 2: HSP ----
run_step("🧠 Step 2: Hidden-state prediction...", [
    "python", str(PICRUST2_SCRIPTS / "hsp.py"),
    "-i", str(PLACED_TREE),
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
