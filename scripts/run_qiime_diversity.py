import subprocess
from pathlib import Path
import shutil
import pandas as pd

# Base directories
DATA_DIR = Path("../data")
EXPORT_DIR = DATA_DIR / "exported"
RESULTS_DIR = EXPORT_DIR / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

metadata_file = DATA_DIR / "metadata.tsv"

# Keep exactly these markers
markers = ["16s", "its", "its_amf", "18s"]


def run_cmd(cmd, label):
    """Run a shell command and fail loudly if it errors."""
    print(f"\n▶️ {label}")
    print("   $ " + " ".join(map(str, cmd)))
    subprocess.run(cmd, check=True)


def get_sampling_depth_and_ids(marker, quantile=0.10):
    """
    Export the QIIME2 feature table to BIOM, convert to TSV, and compute:
      - sampling depth (10th percentile by default; more robust than min)
      - list of sample IDs (table columns)
    """
    table_qza = DATA_DIR / f"{marker}_table.qza"
    if not table_qza.exists():
        print(f"❌ Skipping {marker.upper()}: Missing table {table_qza}")
        return None, None

    temp_dir = RESULTS_DIR / f"tmp_{marker}"
    if temp_dir.exists():
        shutil.rmtree(temp_dir)
    temp_dir.mkdir(parents=True, exist_ok=True)

    run_cmd(
        [
            "qiime", "tools", "export",
            "--input-path", str(table_qza),
            "--output-path", str(temp_dir),
        ],
        f"{marker.upper()} - Exporting table"
    )

    biom_file = temp_dir / "feature-table.biom"
    if not biom_file.exists():
        print(f"❌ {marker.upper()}: Exported BIOM not found at {biom_file}")
        shutil.rmtree(temp_dir)
        return None, None

    tsv_file = temp_dir / "feature-table.tsv"
    run_cmd(
        [
            "biom", "convert",
            "-i", str(biom_file),
            "-o", str(tsv_file),
            "--to-tsv",
        ],
        f"{marker.upper()} - Convert BIOM to TSV"
    )

    # TSV has a comment row, then header row; skiprows=1 is correct here
    df = pd.read_csv(tsv_file, sep="\t", skiprows=1, index_col=0)

    if df.empty or df.shape[1] == 0:
        print(f"⚠️ {marker.upper()}: Table empty or no samples.")
        shutil.rmtree(temp_dir)
        return None, None

    sample_freqs = df.sum(axis=0)
    sample_ids = list(df.columns)

    # Robust depth choice (vendor-like pipelines rarely use the absolute minimum)
    depth = int(sample_freqs.quantile(quantile))
    if depth < 1:
        depth = int(sample_freqs.min())

    print(f"📉 {marker.upper()} sampling depth = {depth} (quantile={quantile})")
    shutil.rmtree(temp_dir)
    return depth, sample_ids


def make_matched_metadata(metadata_path: Path, table_sample_ids, out_path: Path) -> Path:
    """
    Read metadata TSV and write a cleaned + filtered metadata TSV that matches the
    sample IDs in the feature table.

    Fixes common issues:
      - First column header not '#SampleID'
      - Windows CRLF (\r)
      - Leading/trailing whitespace in sample IDs
      - Extra samples not present in the table (filters them out)
    """
    if not metadata_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")

    md = pd.read_csv(metadata_path, sep="\t", dtype=str)

    if md.shape[1] < 1:
        raise ValueError("Metadata file has no columns.")

    # Ensure first column is '#SampleID'
    first_col = md.columns[0]
    if first_col != "#SampleID":
        md = md.rename(columns={first_col: "#SampleID"})

    # Clean sample IDs
    md["#SampleID"] = (
        md["#SampleID"]
        .astype(str)
        .str.replace("\r", "", regex=False)
        .str.strip()
    )

    # Filter metadata to only IDs present in table
    table_set = set(table_sample_ids)
    md_filtered = md[md["#SampleID"].isin(table_set)].copy()

    if md_filtered.empty:
        md_ids_preview = md["#SampleID"].dropna().unique().tolist()[:10]
        table_ids_preview = list(table_sample_ids)[:10]
        raise ValueError(
            "After filtering, metadata has 0 matching sample IDs.\n"
            f"Example metadata IDs: {md_ids_preview}\n"
            f"Example table IDs: {table_ids_preview}\n"
            "IDs likely differ by naming (e.g., '-' vs '_' or suffix differences)."
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    md_filtered.to_csv(out_path, sep="\t", index=False)
    return out_path


for marker in markers:
    print(f"\n🔍 Checking {marker.upper()} files...")

    table = DATA_DIR / f"{marker}_table.qza"
    rep_seqs = DATA_DIR / f"{marker}_rep_seqs.qza"

    if not table.exists() or not rep_seqs.exists():
        print(f"❌ Skipping {marker.upper()}: missing table/rep_seqs.")
        continue

    depth, table_ids = get_sampling_depth_and_ids(marker, quantile=0.10)
    if depth is None or depth < 1 or not table_ids:
        print(f"⚠️ Skipping {marker.upper()}: invalid sampling depth or no table IDs.")
        continue

    marker_out = RESULTS_DIR / f"{marker}_core_metrics_results"
    if marker_out.exists():
        shutil.rmtree(marker_out)
    marker_out.mkdir(parents=True, exist_ok=True)

    # Create metadata file matched to THIS marker's sample IDs
    matched_md = marker_out / "metadata_matched.tsv"
    try:
        make_matched_metadata(metadata_file, table_ids, matched_md)
        print(f"🧾 Using matched metadata: {matched_md}")
    except Exception as e:
        print(f"❌ Skipping {marker.upper()}: metadata mismatch.\n{e}")
        continue

    # Filter low-depth samples BEFORE metrics (keeps results meaningful)
    filtered_table = marker_out / f"{marker}_table_filtered.qza"
    run_cmd(
        [
            "qiime", "feature-table", "filter-samples",
            "--i-table", str(table),
            "--p-min-frequency", str(depth),
            "--o-filtered-table", str(filtered_table),
        ],
        f"{marker.upper()} - Filter low-depth samples (min-frequency={depth})"
    )

    table_for_metrics = filtered_table

    if marker in ["16s", "18s"]:
        # Build phylogeny
        run_cmd(
            [
                "qiime", "phylogeny", "align-to-tree-mafft-fasttree",
                "--i-sequences", str(rep_seqs),
                "--o-alignment", str(marker_out / "aligned-rep-seqs.qza"),
                "--o-masked-alignment", str(marker_out / "masked-aligned-rep-seqs.qza"),
                "--o-tree", str(marker_out / "unrooted-tree.qza"),
                "--o-rooted-tree", str(marker_out / "rooted-tree.qza"),
            ],
            f"{marker.upper()} - Build Phylogeny"
        )

        run_cmd(
            [
                "qiime", "diversity", "core-metrics-phylogenetic",
                "--i-phylogeny", str(marker_out / "rooted-tree.qza"),
                "--i-table", str(table_for_metrics),
                "--p-sampling-depth", str(depth),
                "--m-metadata-file", str(matched_md),
                "--o-rarefied-table", str(marker_out / "rarefied_table.qza"),
                "--o-faith-pd-vector", str(marker_out / "faith_pd_vector.qza"),
                "--o-observed-features-vector", str(marker_out / "observed_features_vector.qza"),
                "--o-shannon-vector", str(marker_out / "shannon_vector.qza"),
                "--o-evenness-vector", str(marker_out / "evenness_vector.qza"),
                "--o-unweighted-unifrac-distance-matrix", str(marker_out / "unweighted_unifrac_distance_matrix.qza"),
                "--o-weighted-unifrac-distance-matrix", str(marker_out / "weighted_unifrac_distance_matrix.qza"),
                "--o-jaccard-distance-matrix", str(marker_out / "jaccard_distance_matrix.qza"),
                "--o-bray-curtis-distance-matrix", str(marker_out / "bray_curtis_distance_matrix.qza"),
                "--o-unweighted-unifrac-pcoa-results", str(marker_out / "unweighted_unifrac_pcoa.qza"),
                "--o-weighted-unifrac-pcoa-results", str(marker_out / "weighted_unifrac_pcoa.qza"),
                "--o-jaccard-pcoa-results", str(marker_out / "jaccard_pcoa.qza"),
                "--o-bray-curtis-pcoa-results", str(marker_out / "bray_curtis_pcoa.qza"),
                "--o-unweighted-unifrac-emperor", str(marker_out / "unweighted_unifrac_emperor.qzv"),
                "--o-weighted-unifrac-emperor", str(marker_out / "weighted_unifrac_emperor.qzv"),
                "--o-jaccard-emperor", str(marker_out / "jaccard_emperor.qzv"),
                "--o-bray-curtis-emperor", str(marker_out / "bray_curtis_emperor.qzv"),
            ],
            f"{marker.upper()} - Core Metrics (Phylo)"
        )

    else:
        run_cmd(
            [
                "qiime", "diversity", "core-metrics",
                "--i-table", str(table_for_metrics),
                "--p-sampling-depth", str(depth),
                "--m-metadata-file", str(matched_md),
                "--o-rarefied-table", str(marker_out / "rarefied_table.qza"),
                "--o-observed-features-vector", str(marker_out / "observed_features_vector.qza"),
                "--o-shannon-vector", str(marker_out / "shannon_vector.qza"),
                "--o-evenness-vector", str(marker_out / "evenness_vector.qza"),
                "--o-jaccard-distance-matrix", str(marker_out / "jaccard_distance_matrix.qza"),
                "--o-bray-curtis-distance-matrix", str(marker_out / "bray_curtis_distance_matrix.qza"),
                "--o-jaccard-pcoa-results", str(marker_out / "jaccard_pcoa.qza"),
                "--o-bray-curtis-pcoa-results", str(marker_out / "bray_curtis_pcoa.qza"),
                "--o-jaccard-emperor", str(marker_out / "jaccard_emperor.qzv"),
                "--o-bray-curtis-emperor", str(marker_out / "bray_curtis_emperor.qzv"),
            ],
            f"{marker.upper()} - Core Metrics"
        )

print("\n✅ All completed (where data available).")
