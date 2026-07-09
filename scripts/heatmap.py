import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import re

# === CONFIG (paths come from config.py; no hard-coded layout) ===
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
try:
    import config as cfg
    input_excel = cfg.EXPORTED_DIR / "asv_tables_combined.xlsx"
    output_dir = cfg.DATA_DIR / "visualizations"
    TAX_RANKS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
except Exception:
    input_excel = Path("../data/exported/asv_tables_combined.xlsx")
    output_dir = Path("../data/visualizations")
    TAX_RANKS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
output_dir.mkdir(parents=True, exist_ok=True)

# Columns that are NOT sample-abundance columns.
NON_SAMPLE = set(TAX_RANKS) | {"ASV ID", "ASV_ID", "Feature ID", "OTU ID"}


def detect_sample_columns(df):
    """Dynamically find numeric sample columns - do NOT assume a 'sample' prefix."""
    cols = []
    for c in df.columns:
        if str(c) in NON_SAMPLE:
            continue
        if pd.to_numeric(df[c], errors="coerce").notna().any():
            cols.append(c)
    return cols


# === Helper to sort sample names numerically ===
def sort_samples_numerically(sample_names):
    def extract_number(s):
        match = re.search(r'\d+', s)
        return int(match.group()) if match else float('inf')

    return sorted(sample_names, key=extract_number)


# === Plotting Functions ===
def plot_barplot(df, marker):
    order_sums = df.groupby("Order")[sample_cols].sum().sum(axis=1)
    top_orders = order_sums.sort_values(ascending=False).head(10).index
    df_top = df[df["Order"].isin(top_orders)].groupby("Order")[sample_cols].sum()

    # Reorder sample columns numerically
    sorted_samples = sort_samples_numerically(df_top.columns)
    df_top = df_top[sorted_samples]

    df_top.T.plot(kind="bar", stacked=True, figsize=(12, 6))
    plt.title(f"Top 10 Orders - Barplot ({marker.upper()})")
    plt.xlabel("Samples")
    plt.ylabel("Abundance")
    plt.xticks(rotation=45, ha='right')
    plt.legend(title="Order", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(output_dir / f"{marker}_top10_orders_barplot.png")
    plt.close()


def plot_heatmap(df, marker):
    order_sums = df.groupby("Order")[sample_cols].sum().sum(axis=1)
    top_orders = order_sums.sort_values(ascending=False).head(10).index
    heatmap_data = df[df["Order"].isin(top_orders)].groupby("Order")[sample_cols].sum()

    # Reorder sample columns numerically
    sorted_samples = sort_samples_numerically(heatmap_data.columns)
    heatmap_data = heatmap_data[sorted_samples]

    plt.figure(figsize=(10, 6))
    sns.heatmap(heatmap_data, cmap="viridis", annot=False)
    plt.title(f"Top 10 Orders - Heatmap ({marker.upper()})")
    plt.ylabel("Order")
    plt.xlabel("Sample")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_dir / f"{marker}_top10_orders_heatmap.png")
    plt.close()


# === Process Each Sheet ===
excel_file = pd.ExcelFile(input_excel)
for sheet in excel_file.sheet_names:
    df = excel_file.parse(sheet)

    if "Order" not in df.columns:
        print(f"⚠️ Skipping {sheet} - no Order column.")
        continue

    sample_cols = detect_sample_columns(df)
    if not sample_cols:
        print(f"⚠️ Skipping {sheet} - no numeric sample columns detected.")
        continue

    print(f"📊 Generating plots for: {sheet}")
    plot_barplot(df, sheet)
    plot_heatmap(df, sheet)

print(f"\n✅ All plots saved to: {output_dir.resolve()}")
