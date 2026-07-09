import pandas as pd
import os


def load_order_table(xls: pd.ExcelFile, sheet_name: str) -> pd.DataFrame:
    """Load and normalize order-level ASV data from a specific sheet."""
    df = xls.parse(sheet_name).set_index("Order")
    df_pct = df.div(df.sum()) * 100  # Convert to relative abundance (%)
    return df_pct.T  # Transpose to SampleID × Order


def prefix_columns(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    """Prefix each column with the data source (e.g., Bacteria_, Fungi_)."""
    df.columns = [f"{prefix}{col}" for col in df.columns]
    return df


def align_sample_ids(df_its: pd.DataFrame) -> pd.DataFrame:
    """Normalize ITS sample names to match 16S (remove _ITS suffix)."""
    df_its.index = df_its.index.str.replace("_ITS", "", regex=False)
    return df_its


def main():
    # === Config ===
    input_path = os.path.join("..", "data", "exported", "asv_tables_combined.xlsx")
    output_path = os.path.join("..", "data", "exported", "all_orders_relative_abundance.csv")

    # === Load Excel File ===
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")

    xls = pd.ExcelFile(input_path)

    # === Load and Process Tables ===
    df_16s = prefix_columns(load_order_table(xls, "16S_TOP_ORDERS"), "Bacteria_")
    df_its = align_sample_ids(load_order_table(xls, "ITS_TOP_ORDERS"))
    df_its = prefix_columns(df_its, "Fungi_")

    # === Merge and Save ===
    df_combined = pd.concat([df_16s, df_its], axis=1).reset_index()
    df_combined.rename(columns={"index": "SampleID"}, inplace=True)

    df_combined.to_csv(output_path, index=False)
    print(f"✅ File saved successfully:\n{output_path}")


if __name__ == "__main__":
    main()
