import pandas as pd
from pathlib import Path
import subprocess

# === CONFIG ===
markers = {
    "16s": "../data/exported/16s",
    "its": "../data/exported/its",
    "18s": "../data/exported/18s",
}
output_dir = Path("../data/krona_charts")
output_dir.mkdir(parents=True, exist_ok=True)

print("🌈 Generating Krona charts...\n")

krona_inputs = []

for marker, folder in markers.items():
    folder = Path(folder)
    feature_table_file = folder / "cleaned-feature-table.tsv"
    taxonomy_file = folder / "taxonomy.tsv"
    mapping_file = folder / "asv_mapping.tsv"
    krona_input_file = output_dir / f"{marker}_krona_input.txt"
    krona_output_file = output_dir / f"{marker}_krona_chart.html"

    if not feature_table_file.exists() or not taxonomy_file.exists():
        print(f"⚠️ Skipping {marker.upper()} - missing files.\n")
        continue

    # Load feature table and compute total counts
    df_feat = pd.read_csv(feature_table_file, sep='\t', skiprows=1, index_col=0)
    df_feat["Total"] = df_feat.sum(axis=1)

    # Load taxonomy
    df_tax = pd.read_csv(taxonomy_file, sep='\t')

    # Handle optional mapping
    if mapping_file.exists():
        df_map = pd.read_csv(mapping_file, sep='\t', names=["ASV ID", "Feature ID"])
        df_map.set_index("Feature ID", inplace=True)
        df_tax["ASV ID"] = df_tax["Feature ID"].map(df_map["ASV ID"])
        df_tax.dropna(subset=["ASV ID"], inplace=True)
        df_tax.set_index("ASV ID", inplace=True)
    else:
        df_tax.set_index("Feature ID", inplace=True)

    # Match
    matched_ids = df_feat.index.intersection(df_tax.index)
    if matched_ids.empty:
        print(f"⚠️ No overlapping Feature IDs found for {marker.upper()}\n")
        continue

    df_feat = df_feat.loc[matched_ids]
    df_tax = df_tax.loc[matched_ids]

    # Build Krona input
    krona_df = pd.DataFrame({
        "Total": df_feat["Total"],
        "Taxon": df_tax["Taxon"]
    }).dropna()

    # Clean taxonomy strings and prepend dataset marker as root
    krona_df["Taxon"] = krona_df["Taxon"].str.replace(r"k__|d__|p__|c__|o__|f__|g__|s__", "", regex=True)
    krona_df["Taxon"] = marker.upper() + "\t" + krona_df["Taxon"].str.replace(";", "\t")

    # Save input
    krona_df.to_csv(krona_input_file, sep="\t", header=False, index=False)

    # Generate individual Krona chart
    try:
        subprocess.run(["ktImportText", "-n", marker.upper(), str(krona_input_file), "-o", str(krona_output_file)], check=True)
        print(f"✅ Krona chart created for {marker.upper()}: {krona_output_file.name}")
    except subprocess.CalledProcessError:
        print(f"❌ Failed to generate Krona chart for {marker.upper()}")

    # Store only the file path (no -n) for combined chart
    krona_inputs.append((str(krona_input_file),))

# === COMBINED CHART ===
if krona_inputs:
    combined_output = output_dir / "krona_combined_chart.html"
    cmd = ["ktImportText", "-o", str(combined_output)]
    for input_file in krona_inputs:
        cmd += list(input_file)

    try:
        subprocess.run(cmd, check=True)
        print(f"\n🎯 Combined Krona chart created: {combined_output.name}")
    except subprocess.CalledProcessError:
        print("❌ Failed to generate combined Krona chart.")

print(f"\n📁 All Krona charts saved to: {output_dir.resolve()}")
