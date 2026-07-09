import subprocess
from pathlib import Path

def prepare_taxon_file():
    print("🔧 Preparing GAPPA-compatible taxon file...")
    csv_path = Path("../picrust2/picrust2/default_files/bacteria/bacteria_metadata.csv.gz")
    tsv_path = Path("../picrust2/picrust2/default_files/bacteria/bacteria_metadata.tsv")

    # Only generate if not already present
    if not tsv_path.exists():
        result = subprocess.run(
            f"zcat {csv_path} | tail -n +2 | awk -F',' '{{print $1 \"\\t\" $8}}' > {tsv_path}",
            shell=True,
            executable="/bin/bash"
        )
        if result.returncode != 0:
            print("❌ Failed to generate taxon file.")
            exit(1)
    else:
        print("✅ Taxon file already exists.")

def extract_tree_with_gappa():
    print("🌳 Extracting placed tree using GAPPA...")
    prepare_taxon_file()  # Ensure taxon file is ready

    result = subprocess.run([
        "gappa", "examine", "assign",
        "--jplace-path", "../functional_prediction/16s/intermediate/epa_out/epa_result.jplace",
        "--taxon-file", "../picrust2/picrust2/default_files/bacteria/bacteria_metadata.tsv",
        "--out-dir", "../functional_prediction/16s/intermediate",
        "--file-prefix", "placed_seqs",
        "--allow-file-overwriting"
    ], check=True)

if __name__ == "__main__":
    extract_tree_with_gappa()
