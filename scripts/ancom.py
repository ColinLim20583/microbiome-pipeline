from pathlib import Path
import qiime2
from qiime2 import Artifact, Metadata
from qiime2.plugins import composition
from qiime2.plugins.feature_table.methods import filter_features

# === Define base directories ===
DATA_DIR = Path("../data").resolve()
ANCOM_DIR = Path("../ANCOM").resolve()
CLASS_DIR = Path("../classification").resolve()

# Ensure output folder exists
ANCOM_DIR.mkdir(parents=True, exist_ok=True)

# === Metadata ===
metadata_file = DATA_DIR / 'metadata.tsv'
metadata = Metadata.load(str(metadata_file))
group_column = 'treatment'  # Change if using a different column

# === Dataset keys you want to process ===
datasets = ['16s', 'its', 'its_amf', '18s']

# === Process each dataset ===
for dataset in datasets:
    print(f"\n🔄 Running ANCOM for: {dataset}")

    # Define file paths
    table_path = DATA_DIR / f"{dataset}_table.qza"
    taxonomy_path = CLASS_DIR / f"{dataset}_taxonomy.qza"
    comp_output_path = ANCOM_DIR / f"{dataset}_comp_table.qza"
    ancom_output_path = ANCOM_DIR / f"{dataset}_ancom.qzv"

    # Skip dataset if table doesn't exist
    if not table_path.exists():
        print(f"⚠️ Skipping {dataset}: {table_path.name} not found.")
        continue

    # Load and filter feature table
    original_table = Artifact.load(str(table_path))
    filtered_result = filter_features(table=original_table, min_frequency=10)
    filtered_table = filtered_result.filtered_table

    # Add pseudocount
    pseudo_result = composition.methods.add_pseudocount(table=filtered_table)
    pseudo_result.composition_table.save(str(comp_output_path))

    # Run ANCOM
    try:
        ancom_result = composition.visualizers.ancom(
            table=pseudo_result.composition_table,
            metadata=metadata.get_column(group_column)
        )
        ancom_result.visualization.save(str(ancom_output_path))
        print(f"✅ ANCOM completed for {dataset} → {ancom_output_path.name}")
    except Exception as e:
        print(f"❌ Failed to run ANCOM for {dataset}: {e}")
