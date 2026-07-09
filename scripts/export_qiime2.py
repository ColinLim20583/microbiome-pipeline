# Export only
from pathlib import Path
import subprocess
import shutil

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
EXPORT_DIR = DATA_DIR / "exported"
DATASETS = {
    "16s": {"rep_seqs_qza": DATA_DIR / "16s_rep_seqs.qza", "table_qza": DATA_DIR / "16s_table.qza"},
    "its": {"rep_seqs_qza": DATA_DIR / "its_rep_seqs.qza", "table_qza": DATA_DIR / "its_table.qza"},
    "its_amf": {"rep_seqs_qza": DATA_DIR / "its_amf_rep_seqs.qza", "table_qza": DATA_DIR / "its_amf_table.qza"},
    "18s_amf": {"rep_seqs_qza": DATA_DIR / "18s_rep_seqs.qza", "table_qza": DATA_DIR / "18s_table.qza"},
}

print("🔍 Exporting QIIME2 artifacts...\n")

for label, files in DATASETS.items():
    export_path = EXPORT_DIR / label
    export_path.mkdir(parents=True, exist_ok=True)

    rep_qza = files["rep_seqs_qza"]
    table_qza = files["table_qza"]
    rep_fna = export_path / "rep_seqs.fna"
    table_biom = export_path / "feature-table.biom"

    if rep_qza.exists() and not rep_fna.exists():
        print(f"📤 Exporting rep_seqs for {label.upper()}")
        subprocess.run(["qiime", "tools", "export", "--input-path", str(rep_qza), "--output-path", str(export_path)], check=True)
        (export_path / "dna-sequences.fasta").rename(rep_fna)

    if table_qza.exists() and not table_biom.exists():
        print(f"📤 Exporting table for {label.upper()}")
        subprocess.run(["qiime", "tools", "export", "--input-path", str(table_qza), "--output-path", str(export_path)], check=True)
        (export_path / "feature-table.biom").rename(table_biom)
