import os
import subprocess
import sys
from pathlib import Path

print("🔍 Starting taxonomic classification...\n")

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
CLASSIFICATION_DIR = DATA_DIR / "classification"
CLASSIFICATION_DIR.mkdir(parents=True, exist_ok=True)

# Config drives thread limits, confidence, read-orientation and the set of
# markers. Classifier .qza files are DISCOVERED on disk by keyword match
# (config.discover_classifier) instead of being referenced by a hard-coded
# filename - so renaming/upgrading a classifier never silently breaks the run.
sys.path.insert(0, str(BASE_DIR))
try:
    import config as cfg
    CLASSIFY_CFG = cfg.CLASSIFY
    MARKER_LIST = list(cfg.MARKERS) + list(cfg.AMF_MARKERS)
    N_THREADS = str(CLASSIFY_CFG.get("n_jobs", 1))
    CONFIDENCE = str(CLASSIFY_CFG.get("confidence", 0.8))
    READ_ORIENT = str(CLASSIFY_CFG.get("read_orientation", "auto"))
    _discover = cfg.discover_classifier
except Exception:
    cfg = None
    CLASSIFIER_DIR = BASE_DIR / "classifiers"
    MARKER_LIST = ["16s", "its", "18s", "its_amf", "18s_amf"]
    N_THREADS = "1"
    CONFIDENCE = "0.8"
    READ_ORIENT = "auto"

    def _discover(marker):
        # crude fallback: first .qza whose name contains the marker keyword
        key = {"16s": "16s", "its": "its", "18s": "18s",
               "its_amf": "maarjam", "18s_amf": "maarjam"}.get(marker, marker)
        for q in sorted(CLASSIFIER_DIR.glob("*.qza")):
            if key in q.name.lower():
                return q
        return None

os.environ["OMP_NUM_THREADS"] = N_THREADS
ENV = os.environ.copy()

def run(cmd: list[str]) -> None:
    """Run a command with consistent env and nice printing."""
    print("   $ " + " ".join(cmd))
    subprocess.run(cmd, check=True, env=ENV)

# AMF markers reuse their parent region's rep-seqs (ITS / 18S).
REP_SEQS_FOR = {
    "16s": "16s_rep_seqs.qza",
    "its": "its_rep_seqs.qza",
    "18s": "18s_rep_seqs.qza",
    "its_amf": "its_rep_seqs.qza",
    "18s_amf": "18s_rep_seqs.qza",
}

CLASSIFY_PARAMS = {}
for _m in MARKER_LIST:
    CLASSIFY_PARAMS[_m] = {
        "label": _m,
        "rep_seqs": DATA_DIR / REP_SEQS_FOR.get(_m, f"{_m}_rep_seqs.qza"),
        "method": "sklearn",
        "classifier": _discover(_m),   # resolved dynamically from disk
        "confidence": CONFIDENCE,
        "read_orientation": READ_ORIENT,
    }

for key, p in CLASSIFY_PARAMS.items():
    label = p["label"]
    rep_seqs = p["rep_seqs"]

    taxonomy_qza = CLASSIFICATION_DIR / f"{label}_taxonomy.qza"
    taxonomy_qzv = CLASSIFICATION_DIR / f"{label}_taxonomy.qzv"

    if not rep_seqs.exists():
        print(f"⚠️ Skipping {label.upper()}: missing rep_seqs: {rep_seqs}")
        continue

    print(f"🧬 Classifying {label.upper()} using {p['method']}...")

    try:
        if p["method"] == "sklearn":
            classifier = p["classifier"]
            if classifier is None or not Path(classifier).exists():
                print(f"⚠️ Skipping {label.upper()}: no classifier found on disk "
                      f"(looked for keywords for '{label}' in classifiers/).")
                continue

            print(f"   using classifier: {Path(classifier).name}")
            run([
                "qiime", "feature-classifier", "classify-sklearn",
                "--i-classifier", str(classifier),
                "--i-reads", str(rep_seqs),
                "--p-confidence", str(p.get("confidence", "0.8")),
                "--p-read-orientation", str(p.get("read_orientation", "auto")),
                "--p-n-jobs", N_THREADS,
                "--o-classification", str(taxonomy_qza),
            ])

        elif p["method"] == "blast":
            # (kept for completeness; you likely won't need it now)
            ref_reads = p["reference_reads"]
            ref_tax = p["reference_taxonomy"]
            if not ref_reads.exists() or not ref_tax.exists():
                print(f"⚠️ Skipping {label.upper()}: missing BLAST refs.")
                continue

            blast_hits = CLASSIFICATION_DIR / f"{label}_blast_hits.qza"
            run([
                "qiime", "feature-classifier", "classify-consensus-blast",
                "--i-query", str(rep_seqs),
                "--i-reference-reads", str(ref_reads),
                "--i-reference-taxonomy", str(ref_tax),
                "--p-perc-identity", p.get("perc_identity", "0.8"),
                "--p-num-threads", "1",
                "--o-classification", str(taxonomy_qza),
                "--o-search-results", str(blast_hits),
            ])
        else:
            print(f"⚠️ Skipping {label.upper()}: unknown method {p['method']}")
            continue

        # Tabulate results
        run([
            "qiime", "metadata", "tabulate",
            "--m-input-file", str(taxonomy_qza),
            "--o-visualization", str(taxonomy_qzv),
        ])

        print(f"✅ Done: {label.upper()} taxonomy saved:\n   {taxonomy_qza}\n   {taxonomy_qzv}\n")

    except subprocess.CalledProcessError as e:
        print(f"❌ Classification failed for {label.upper()}\n   {e}\n")

print("🏁 Taxonomic classification complete.")
