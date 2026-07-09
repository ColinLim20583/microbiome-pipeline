import subprocess
from pathlib import Path

print("🔍 Exporting QIIME2 artifacts...\n")

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
EXPORT_DIR = DATA_DIR / "exported"
EXPORT_DIR.mkdir(exist_ok=True)

datasets = {
    "16s": {
        "table": DATA_DIR / "16s_table.qza",
        "rep_seqs": DATA_DIR / "16s_rep_seqs.qza",
    },
    "its": {
        "table": DATA_DIR / "its_table.qza",
        "rep_seqs": DATA_DIR / "its_rep_seqs.qza",
    },
    "18s": {
        "table": DATA_DIR / "18s_table.qza",
        "rep_seqs": DATA_DIR / "18s_rep_seqs.qza",
    },
}

for label, paths in datasets.items():
    out_dir = EXPORT_DIR / label
    out_dir.mkdir(exist_ok=True)

    print(f"\n📦 {label.upper()}")
    print(f"   table:    {paths['table']}")
    print(f"   rep_seqs: {paths['rep_seqs']}")
    print(f"   out:      {out_dir}")

    if paths["table"].exists():
        print("   📤 exporting table...")
        subprocess.run([
            "qiime", "tools", "export",
            "--input-path", str(paths["table"]),
            "--output-path", str(out_dir),
        ], check=True)
    else:
        print("   ⚠️ missing table")

    if paths["rep_seqs"].exists():
        print("   📤 exporting rep_seqs...")
        subprocess.run([
            "qiime", "tools", "export",
            "--input-path", str(paths["rep_seqs"]),
            "--output-path", str(out_dir),
        ], check=True)
    else:
        print("   ⚠️ missing rep_seqs")

print("\n✅ Export complete.")
