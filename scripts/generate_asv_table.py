from pathlib import Path
import pandas as pd
import subprocess
import re

# =========================
# CONFIG
# =========================
markers = {
    "16s": "../data/exported/16s",
    "its": "../data/exported/its",
    "18s": "../data/exported/18s",
}

taxonomy_qza_dir = Path("../data/classification")

output_excel = Path("../data/exported/asv_tables_combined.xlsx")
output_excel.parent.mkdir(parents=True, exist_ok=True)

TAXONOMY_NAME = {
    "16s": "16s",
    "its": "its",
    "18s": "18s",
}

rank_map = {
    "d__": "Kingdom",
    "k__": "Kingdom",
    "p__": "Phylum",
    "c__": "Class",
    "o__": "Order",
    "f__": "Family",
    "g__": "Genus",
    "s__": "Species",
}
tax_ranks = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


# =========================
# HELPERS
# =========================
def extract_sample_number(name: str):
    m = re.search(r"\d+", str(name))
    return int(m.group()) if m else float("inf")


def export_taxonomy_qza_to_tsv(taxonomy_qza: Path, out_dir: Path) -> Path:
    """
    Export taxonomy.qza into out_dir, returning the path to taxonomy.tsv.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    taxonomy_tsv = out_dir / "taxonomy.tsv"
    if taxonomy_tsv.exists():
        return taxonomy_tsv

    if not taxonomy_qza.exists():
        return taxonomy_tsv

    subprocess.run(
        ["qiime", "tools", "export", "--input-path", str(taxonomy_qza), "--output-path", str(out_dir)],
        check=True,
    )
    return taxonomy_tsv


def load_asv_mapping(mapping_file: Path) -> pd.DataFrame:
    """
    Load asv_mapping.tsv and return a DataFrame indexed by Feature ID with column 'ASV ID'.
    Supports both orientations:
      ASV_ID <tab> Feature_ID
      Feature_ID <tab> ASV_ID
    """
    df_map = pd.read_csv(mapping_file, sep="\t", header=None, dtype=str)

    first = str(df_map.iloc[0, 0]) if not df_map.empty else ""
    if first.startswith(("16S_", "ITS_", "18S_", "ASV")):
        df_map.columns = ["ASV ID", "Feature ID"]
    else:
        df_map.columns = ["Feature ID", "ASV ID"]

    df_map = df_map.dropna()
    df_map["Feature ID"] = df_map["Feature ID"].astype(str)
    df_map["ASV ID"] = df_map["ASV ID"].astype(str)
    return df_map.set_index("Feature ID")


def load_taxonomy_as_asv_index(taxonomy_tsv: Path, df_map: pd.DataFrame) -> pd.DataFrame:
    """
    Load taxonomy.tsv and return a DataFrame indexed by ASV ID with rank columns.
    """
    df_tax = pd.read_csv(taxonomy_tsv, sep="\t", dtype=str)

    feature_col = None
    for cand in ["Feature ID", "FeatureID", "feature-id", "feature_id"]:
        if cand in df_tax.columns:
            feature_col = cand
            break
    if feature_col is None:
        raise ValueError(f"Can't find Feature ID column in {taxonomy_tsv}. Columns: {list(df_tax.columns)}")

    if "Taxon" not in df_tax.columns:
        raise ValueError(f"Can't find Taxon column in {taxonomy_tsv}. Columns: {list(df_tax.columns)}")

    df_tax = df_tax.set_index(feature_col)
    tax_series = df_tax["Taxon"].fillna("").astype(str)

    tax_expanded = pd.DataFrame(index=tax_series.index, columns=tax_ranks, dtype=str)
    tax_expanded.loc[:, :] = "Unassigned"

    for fid, tax_str in tax_series.items():
        parts = [p.strip() for p in tax_str.split(";") if p.strip()]
        rank_values = {}
        for p in parts:
            for prefix, rank in rank_map.items():
                if p.startswith(prefix):
                    val = p[len(prefix):].strip()
                    rank_values[rank] = val if val else "Unassigned"
                    break
        for r in tax_ranks:
            tax_expanded.loc[fid, r] = rank_values.get(r, "Unassigned")

    def to_asv(fid: str) -> str:
        return df_map.loc[fid, "ASV ID"] if fid in df_map.index else fid

    tax_expanded.index = [to_asv(str(x)) for x in tax_expanded.index]
    return tax_expanded


# =========================
# MAIN
# =========================
writer = pd.ExcelWriter(output_excel, engine="openpyxl")
sheets_written = 0

for marker, folder in markers.items():
    folder_path = Path(folder)
    tsv_file = folder_path / "cleaned-feature-table.tsv"
    mapping_file = folder_path / "asv_mapping.tsv"

    tax_prefix = TAXONOMY_NAME.get(marker, marker)
    taxonomy_qza = taxonomy_qza_dir / f"{tax_prefix}_taxonomy.qza"
    taxonomy_tsv = folder_path / "taxonomy.tsv"

    print(f"\n=== {marker.upper()} ===")
    print(f"Feature table: {tsv_file}")
    print(f"ASV mapping:   {mapping_file}")
    print(f"Taxonomy QZA:  {taxonomy_qza}")

    if not tsv_file.exists():
        print(f"❌ Skipping {marker.upper()} - TSV file not found.")
        continue
    if not mapping_file.exists():
        print(f"❌ Skipping {marker.upper()} - mapping file not found.")
        continue

    # Load feature table (BIOM TSV usually has 1 comment row)
    df = pd.read_csv(tsv_file, sep="\t", skiprows=1, index_col=0)
    df.index = df.index.astype(str)

    # Load mapping
    df_map = load_asv_mapping(mapping_file)

    # Export taxonomy if needed
    if taxonomy_qza.exists():
        try:
            taxonomy_tsv = export_taxonomy_qza_to_tsv(taxonomy_qza, folder_path)
            if taxonomy_tsv.exists():
                print(f"📤 Taxonomy TSV ready: {taxonomy_tsv}")
            else:
                print(f"⚠️ Taxonomy export ran, but taxonomy.tsv not found in {folder_path}")
        except subprocess.CalledProcessError as e:
            print(f"❌ Failed to export taxonomy for {marker.upper()}: {e}")
    else:
        print(f"⚠️ No taxonomy.qza found for {marker.upper()} at: {taxonomy_qza}")

    # Add ASV ID column by mapping feature-table index
    df["ASV ID"] = df.index.map(lambda fid: df_map.loc[fid, "ASV ID"] if fid in df_map.index else fid)
    df = df.set_index("ASV ID")

    # Identify numeric sample columns ONLY
    # Convert everything else to numeric where possible
    df_numeric = df.apply(pd.to_numeric, errors="coerce")
    sample_cols = [c for c in df_numeric.columns if df_numeric[c].notna().any()]

    # Keep only those numeric columns as sample abundance columns
    df = df[sample_cols].copy()

    # Add taxonomy (default Unassigned)
    for rank in tax_ranks:
        df[rank] = "Unassigned"

    if taxonomy_tsv.exists():
        try:
            tax_expanded = load_taxonomy_as_asv_index(taxonomy_tsv, df_map)
            df = df.drop(columns=tax_ranks, errors="ignore").join(tax_expanded, how="left")
            df = df.fillna("Unassigned")
        except Exception as e:
            print(f"⚠️ Could not load/merge taxonomy for {marker.upper()}: {e}")
    else:
        print(f"⚠️ No taxonomy.tsv found for {marker.upper()} (taxonomy will be Unassigned)")

    df = df.reset_index()

    # Sort sample columns nicely
    sample_cols = sorted(sample_cols, key=extract_sample_number)

    final_cols = ["ASV ID"] + sample_cols + tax_ranks
    df = df[final_cols]

    # Write ASV sheet
    sheet_name = marker.upper()
    df.to_excel(writer, sheet_name=sheet_name, index=False)
    print(f"📄 Wrote sheet: {sheet_name}")
    sheets_written += 1

    # Top 10 Orders
    if "Order" in df.columns:
        df_order = df.copy()
        df_order["Order"] = df_order["Order"].fillna("Unassigned")
        order_summary = df_order.groupby("Order")[sample_cols].sum(numeric_only=True)
        order_summary["Total"] = order_summary.sum(axis=1)

        assigned = order_summary.drop(index="Unassigned", errors="ignore")
        top_orders = assigned.sort_values("Total", ascending=False).head(10).drop(columns="Total")
        top_orders = top_orders.reset_index()

        cols_after = [c for c in top_orders.columns if c != "Order"]
        top_orders = top_orders[["Order"] + sorted(cols_after, key=extract_sample_number)]

        top_sheet = f"{marker.upper()}_TOP_ORDERS"
        top_orders.to_excel(writer, sheet_name=top_sheet, index=False)
        print(f"⭐ Wrote sheet: {top_sheet}")

if sheets_written:
    writer.close()
    print(f"\n✅ Combined ASV table saved to: {output_excel.resolve()}")
else:
    print("\n⚠️ No ASV tables written — Excel file was not saved.")
