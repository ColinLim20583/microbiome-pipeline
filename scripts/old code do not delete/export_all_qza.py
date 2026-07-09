import subprocess
from pathlib import Path

# Define directories
DATA_DIR = Path("../data")
EXPORT_DIR = DATA_DIR / "exported"
EXPORT_DIR.mkdir(exist_ok=True)

# Define datasets and corresponding filenames
datasets = {
    "16s": {
        "table": DATA_DIR / "16s_table.qza",
        "rep_seqs": DATA_DIR / "16s_rep_seqs.qza"
    },
    "its": {
        "table": DATA_DIR / "its_table.qza",
        "rep_seqs": DATA_DIR / "its_rep_seqs.qza"
    },
    "its_amf": {
        "table": DATA_DIR / "its_amf_table.qza",
        "rep_seqs": DATA_DIR / "its_rep_seqs.qza"  # same rep_seqs as ITS
    },
    "18s_amf": {
        "table": DATA_DIR / "18s_amf_table.qza",
        "rep_seqs": DATA_DIR / "18s_rep_seqs.qza"
    }
}

# Loop through and export
for label, paths in datasets.items():
    export_path = EXPORT_DIR / label
    export_path.mkdir(exist_ok=True)
    print(f"\n📦 Exporting {label.upper()}...")

    # Export feature table
    if paths["table"].exists():
        print(f"🔹 Exporting {label} table...")
        subprocess.run([
            "qiime", "tools", "export",
            "--input-path", str(paths["table"]),
            "--output-path", str(export_path)
        ], check=True)
    else:
        print(f"⚠️ {label} table not found: {paths['table']}")

    # Export rep seqs
    if paths["rep_seqs"].exists():
        print(f"🔹 Exporting {label} representative sequences...")
        subprocess.run([
            "qiime", "tools", "export",
            "--input-path", str(paths["rep_seqs"]),
            "--output-path", str(export_path)
        ], check=True)
    else:
        print(f"⚠️ {label} rep_seqs not found: {paths['rep_seqs']}")

print("\n✅ Export complete. Check the 'exported' folder.")
