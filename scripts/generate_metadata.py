import os
import pandas as pd

BASE_DIR = os.path.abspath(os.path.dirname(__file__))
DATA_DIR = os.path.join(BASE_DIR, "../data")
OUTPUT_FILE = os.path.join(DATA_DIR, "metadata.tsv")

def load_metadata(mode):
    path = os.path.join(DATA_DIR, f"{mode}_metadata.tsv")
    if not os.path.exists(path):
        print(f"⚠️ Metadata file for {mode} not found at {path}")
        return pd.DataFrame()

    df = pd.read_csv(path, sep="\t", dtype=str)

    # Normalize column names
    df.columns = [c.strip().lower().replace("-", "_") for c in df.columns]

    # Sanitize sample id column
    if "sample_id" in df.columns:
        df["sample_id"] = df["sample_id"].astype(str).str.replace("\r", "", regex=False).str.strip()
    elif "sample-id" in df.columns:
        df["sample_id"] = df["sample-id"].astype(str).str.replace("\r", "", regex=False).str.strip()
    else:
        print(f"❌ No sample-id column in {path}")
        return pd.DataFrame()

    df["mode"] = mode
    return df

def strip_marker_suffix(sample_id: str) -> str:
    sample_id = str(sample_id)
    for tag in ["_16s", "_18s", "_its", "_its_amf", "_ITS", "_16S", "_18S", "_ITS_AMF"]:
        if sample_id.endswith(tag):
            return sample_id[: -len(tag)]
    return sample_id

def main():
    print("📄 Generating QIIME metadata (metadata.tsv)...")
    combined = []

    for mode in ["16s", "its", "18s"]:
        df = load_metadata(mode)
        if df.empty:
            continue

        for _, row in df.iterrows():
            sid = row.get("sample_id", "")
            if not sid:
                print(f"⚠️ Missing sample ID in {mode}, skipping row.")
                continue

            base_id = strip_marker_suffix(sid)

            combined.append({
                "#SampleID": sid,
                "treatment": row.get("treatment", ""),
                "timepoint": row.get("timepoint", ""),
                "location": row.get("location", ""),
                "soil": row.get("soil", ""),
                "crop": row.get("crop", ""),
                "ph": row.get("ph", ""),
                "note": row.get("note", "")
            })

    if not combined:
        print("⚠️ No metadata to write. Exiting.")
        return

    final_df = pd.DataFrame(combined)

    # Collapse duplicates (same sample seen in multiple modes)
    final_df = (
        final_df.groupby("#SampleID", as_index=False)
        .agg(lambda x: next((v for v in x if pd.notna(v) and str(v).strip() != ""), ""))
    )

    final_df.to_csv(OUTPUT_FILE, sep="\t", index=False)
    print(f"✅ QIIME metadata written to {OUTPUT_FILE}")
    print(final_df)

if __name__ == "__main__":
    main()
