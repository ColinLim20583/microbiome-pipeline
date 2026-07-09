import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# ------------ CONFIG ------------
DATA_DIR = Path("data/exported/16s")  # change to your dataset folder
TAXONOMY_FILE = DATA_DIR / "taxonomy.tsv"
STRATIFIED_FILE = DATA_DIR / "stratified_metagenome.tsv"
OUTPUT_SUMMARY = DATA_DIR / "top_taxa_function_summary.tsv"
OUTPUT_PLOT = DATA_DIR / "top_taxa_function_plot.png"
TOP_N = 20
PLOT = True

# ------------ LOAD DATA ------------
print("📄 Loading stratified metagenome...")
strat_df = pd.read_csv(STRATIFIED_FILE, sep="\t", comment="#")
strat_df.columns = ["function", "asv_id"] + list(strat_df.columns[2:])

print("📄 Loading taxonomy file...")
tax_df = pd.read_csv(TAXONOMY_FILE, sep="\t")
tax_map = dict(zip(tax_df["Feature ID"], tax_df["Taxon"]))

# ------------ MAP ASV TO TAXONOMY ------------
print("🔁 Mapping ASV IDs to taxonomy...")
strat_df["taxonomy"] = strat_df["asv_id"].map(tax_map).fillna("[unmapped]")

# ------------ SUMMARIZE CONTRIBUTIONS ------------
print("📊 Summarizing total abundance...")
melted = strat_df.melt(id_vars=["function", "taxonomy"], var_name="sample", value_name="abundance")
summary = (
    melted.groupby(["function", "taxonomy"])["abundance"]
    .sum()
    .reset_index()
    .sort_values("abundance", ascending=False)
    .head(TOP_N)
)

# ------------ EXPORT RESULTS ------------
summary.to_csv(OUTPUT_SUMMARY, sep="\t", index=False)
print(f"✅ Summary saved to: {OUTPUT_SUMMARY}")

# ------------ OPTIONAL PLOT ------------
if PLOT:
    print("📈 Creating plot...")
    summary["label"] = summary["taxonomy"] + " → " + summary["function"]
    plt.figure(figsize=(10, 6))
    sns.barplot(data=summary, y="label", x="abundance", palette="viridis")
    plt.xlabel("Total Predicted Abundance")
    plt.ylabel("Taxon → Function")
    plt.title(f"Top {TOP_N} Taxon-Function Contributions")
    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT, dpi=300)
    plt.close()
    print(f"✅ Plot saved to: {OUTPUT_PLOT}")
