import subprocess
import os
from pathlib import Path

print("🔍 Starting taxonomic classification...\n")

# ✅ Limit threads to reduce memory usage
os.environ["OMP_NUM_THREADS"] = "1"

# === Directories ===
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
CLASSIFIER_DIR = BASE_DIR / "classifiers"
CLASSIFICATION_DIR = DATA_DIR / "classification"
CLASSIFICATION_DIR.mkdir(exist_ok=True)

# === Parameters for each marker ===
CLASSIFY_PARAMS = {
    "16s": {
        "classifier": CLASSIFIER_DIR / "silva-138-99-nb-classifier.qza",
        "rep_seqs": DATA_DIR / "16s_rep_seqs.qza",
        "label": "16s",
        "method": "sklearn"
    },
    "its": {
        "rep_seqs": DATA_DIR / "its_rep_seqs.qza",
        "label": "its",
        "method": "blast",
        "reference_reads": CLASSIFIER_DIR / "unite-ref-seqs.qza",
        "reference_taxonomy": CLASSIFIER_DIR / "unite-ref-taxonomy.qza"
    },
    "its_amf": {
        "classifier": CLASSIFIER_DIR / "maarjam_its_classifier.qza",
        "rep_seqs": DATA_DIR / "its_amf_rep_seqs.qza",
        "label": "its_amf",
        "method": "sklearn"
    },
    "18s_amf": {
        "classifier": CLASSIFIER_DIR / "maarjam_18s_classifier.qza",
        "rep_seqs": DATA_DIR / "18s_rep_seqs.qza",
        "label": "18s_amf",
        "method": "sklearn"
    },
}

# === Run classification ===
for key, params in CLASSIFY_PARAMS.items():
    label = params["label"]
    taxonomy = CLASSIFICATION_DIR / f"{label}_taxonomy.qza"
    taxonomy_vis = CLASSIFICATION_DIR / f"{label}_taxonomy.qzv"

    rep_seqs = params["rep_seqs"]
    if not rep_seqs.exists():
        print(f"⚠️ Skipping {label.upper()}: Missing representative sequences.")
        continue

    print(f"🧬 Classifying {label.upper()} sequences using {params['method']}...")

    try:
        if params["method"] == "sklearn":
            classifier = params["classifier"]
            if not classifier.exists():
                print(f"⚠️ Skipping {label.upper()}: Missing classifier file.")
                continue
            subprocess.run([
                "qiime", "feature-classifier", "classify-sklearn",
                "--i-classifier", str(classifier),
                "--i-reads", str(rep_seqs),
                "--p-confidence", "0.8",
                "--p-n-jobs", "1",
                "--o-classification", str(taxonomy)
            ], check=True)

        elif params["method"] == "blast":
            reference_reads = params["reference_reads"]
            reference_tax = params["reference_taxonomy"]
            if not reference_reads.exists() or not reference_tax.exists():
                print(f"⚠️ Skipping {label.upper()}: Missing reference reads/taxonomy.")
                continue

            subprocess.run([
                "qiime", "feature-classifier", "classify-consensus-blast",
                "--i-query", str(rep_seqs),
                "--i-reference-reads", str(reference_reads),
                "--i-reference-taxonomy", str(reference_tax),
                "--p-perc-identity", "0.8",
                "--p-num-threads", "1",
                "--o-classification", str(taxonomy),
                "--o-search-results", str(CLASSIFICATION_DIR / f"{label}_blast_hits.qza")
            ], check=True)

        # ✅ Tabulate the result
        subprocess.run([
            "qiime", "metadata", "tabulate",
            "--m-input-file", str(taxonomy),
            "--o-visualization", str(taxonomy_vis)
        ], check=True, env=os.environ.copy())

        print(f"✅ Taxonomy classified for {label.upper()} and saved to: {taxonomy_vis}\n")

    except subprocess.CalledProcessError as e:
        print(f"❌ Classification failed for {label.upper()}")
        print(f"   Error: {e}\n")

print("🏁 Taxonomic classification complete.")
