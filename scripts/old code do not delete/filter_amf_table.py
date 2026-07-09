import subprocess
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
data_dir = BASE_DIR / "data"
classification_dir = data_dir / "classification"

# Input files (ITS)
table_in = data_dir / "its_table.qza"
rep_seqs_in = data_dir / "its_rep_seqs.qza"
taxonomy_in = classification_dir / "its_taxonomy.qza"   # ✅ fixed

# Output files
table_out = data_dir / "its_amf_table.qza"
rep_seqs_out = data_dir / "its_amf_rep_seqs.qza"
rep_seqs_vis = data_dir / "its_amf_rep_seqs.qzv"

print("🔍 Filtering AMF features and representative sequences (Glomeromycetes)...")

# sanity checks
for p in [table_in, rep_seqs_in, taxonomy_in]:
    if not p.exists():
        raise FileNotFoundError(f"Missing required input: {p}")

try:
    subprocess.run([
        "qiime", "taxa", "filter-table",
        "--i-table", str(table_in),
        "--i-taxonomy", str(taxonomy_in),
        "--p-include", "Glomeromycetes",
        "--o-filtered-table", str(table_out)
    ], check=True)

    subprocess.run([
        "qiime", "taxa", "filter-seqs",
        "--i-sequences", str(rep_seqs_in),
        "--i-taxonomy", str(taxonomy_in),
        "--p-include", "Glomeromycetes",
        "--o-filtered-sequences", str(rep_seqs_out)
    ], check=True)

    subprocess.run([
        "qiime", "feature-table", "tabulate-seqs",
        "--i-data", str(rep_seqs_out),
        "--o-visualization", str(rep_seqs_vis)
    ], check=True)

    print("✅ AMF filtering complete.")

except subprocess.CalledProcessError:
    print("❌ An error occurred during AMF filtering or visualization.")
