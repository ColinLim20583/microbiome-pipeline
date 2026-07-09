"""
Filter feature tables down to Arbuscular Mycorrhizal Fungi (AMF, phylum
Glomeromycota) for the ITS / 18S markers.

Previously this file was an accidental copy of generate_asv_table.py and did NOT
filter anything. It now performs a real, config-driven AMF filter:

  * reads each exported marker folder (cleaned-feature-table.tsv + taxonomy.tsv),
  * keeps only features whose taxonomy string contains a known AMF taxon
    (the target list lives in config.py -> amf_target_taxa, so it is editable
    without touching code),
  * writes an AMF-only feature table (TSV) and a combined Excel workbook.

It is intentionally pandas-based (no QIIME dependency) so it can run in the same
environment as the rest of the post-processing and is easy to unit test.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

# ---- config (paths + AMF target taxa) ------------------------------------- #
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
try:
    import config as cfg
    EXPORTED_DIR = cfg.EXPORTED_DIR
    AMF_TARGETS = [t.lower() for t in cfg.get("amf_target_taxa", [])]
    # Which markers can contain AMF (fungi / eukaryotes).
    AMF_SOURCE_MARKERS = ["its", "18s"]
except Exception:
    EXPORTED_DIR = Path("../data/exported")
    AMF_TARGETS = [t.lower() for t in [
        "Glomeromycota", "Glomeromycetes", "Glomerales", "Glomeraceae",
        "Rhizophagus", "Funneliformis", "Claroideoglomus", "Glomus",
        "Gigaspora", "Acaulospora", "Diversispora", "Scutellospora",
        "Paraglomus", "Archaeospora", "Ambispora", "Septoglomus",
    ]]
    AMF_SOURCE_MARKERS = ["its", "18s"]

OUTPUT_XLSX = EXPORTED_DIR / "amf_filtered_tables.xlsx"


def is_amf(taxon_string: str) -> bool:
    """True if the taxonomy string mentions any configured AMF taxon."""
    s = str(taxon_string).lower()
    return any(t in s for t in AMF_TARGETS)


def load_taxonomy(folder: Path) -> pd.DataFrame | None:
    """Load taxonomy.tsv (Feature ID -> Taxon) for a marker if present."""
    tax = folder / "taxonomy.tsv"
    if not tax.exists():
        return None
    df = pd.read_csv(tax, sep="\t", dtype=str)
    # normalise the feature-id column name
    for cand in ["Feature ID", "FeatureID", "feature-id", "feature_id", "#OTU ID"]:
        if cand in df.columns:
            df = df.rename(columns={cand: "Feature ID"})
            break
    return df if {"Feature ID", "Taxon"}.issubset(df.columns) else None


def filter_marker(marker: str) -> pd.DataFrame | None:
    folder = EXPORTED_DIR / marker
    table_tsv = folder / "cleaned-feature-table.tsv"
    if not table_tsv.exists():
        print(f"⏩ {marker.upper()}: no cleaned-feature-table.tsv, skipping.")
        return None

    tax_df = load_taxonomy(folder)
    if tax_df is None:
        print(f"⚠️ {marker.upper()}: taxonomy.tsv missing/invalid, cannot filter AMF.")
        return None

    # BIOM-exported TSVs carry one leading comment row
    table = pd.read_csv(table_tsv, sep="\t", skiprows=1, index_col=0)
    table.index = table.index.astype(str)

    amf_ids = set(tax_df.loc[tax_df["Taxon"].apply(is_amf), "Feature ID"].astype(str))
    kept = table.loc[table.index.intersection(amf_ids)].copy()

    print(f"🍄 {marker.upper()}: {len(kept)} AMF features kept "
          f"out of {len(table)} total.")
    if kept.empty:
        return None

    # attach the taxon string for readability
    tax_lookup = tax_df.set_index("Feature ID")["Taxon"].to_dict()
    kept.insert(0, "Taxon", [tax_lookup.get(i, "") for i in kept.index])
    kept = kept.reset_index().rename(columns={"index": "Feature ID"})
    return kept


def main() -> None:
    print("🍄 Filtering feature tables for AMF (Glomeromycota)...")
    print(f"   AMF target taxa: {', '.join(sorted(set(AMF_TARGETS)))}")

    results = {}
    for marker in AMF_SOURCE_MARKERS:
        out = filter_marker(marker)
        if out is not None and not out.empty:
            results[marker] = out
            # per-marker TSV alongside the source
            per_marker = EXPORTED_DIR / marker / "amf-filtered-feature-table.tsv"
            out.to_csv(per_marker, sep="\t", index=False)
            print(f"   ✅ wrote {per_marker}")

    if not results:
        print("⚠️ No AMF features found in any marker - nothing written.")
        return

    OUTPUT_XLSX.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(OUTPUT_XLSX, engine="openpyxl") as writer:
        for marker, df in results.items():
            df.to_excel(writer, sheet_name=f"{marker.upper()}_AMF"[:31], index=False)
    print(f"\n✅ Combined AMF workbook: {OUTPUT_XLSX.resolve()}")


if __name__ == "__main__":
    main()
