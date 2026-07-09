"""
Central configuration for the peanut microbiome pipeline.

Everything that used to be hard-coded inside individual scripts (primers, DADA2
truncation lengths, classifier file names, marker lists, folder layout) lives
here and can additionally be overridden at run time from a YAML file or from
environment variables. Nothing downstream should contain literal paths, primer
strings, or magic numbers - import from this module instead.

Resolution order for any setting (highest priority first):
    1. Environment variable  PEANUT_<UPPER_SNAKE_KEY>
    2. config.local.yaml sitting next to this file (git-ignored, user supplied)
    3. The DEFAULTS defined below

Classifiers are *discovered* on disk rather than assumed, so renaming a
classifier file never silently breaks classification - the pipeline picks up
whatever is present in the classifiers/ folder and matches it to a marker.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any

try:
    import yaml  # optional; only needed if config.local.yaml is used
except Exception:  # pragma: no cover
    yaml = None

# --------------------------------------------------------------------------- #
# Folder layout (all relative to the project root = this file's directory)
# --------------------------------------------------------------------------- #
BASE_DIR = Path(__file__).resolve().parent
UPLOAD_DIR = BASE_DIR / "uploads"
DATA_DIR = BASE_DIR / "data"
EXPORTED_DIR = DATA_DIR / "exported"
CLASSIFICATION_DIR = DATA_DIR / "classification"
CLASSIFIER_DIR = BASE_DIR / "classifiers"
SCRIPTS_DIR = BASE_DIR / "scripts"
# Curated reference data (evidence base for microbial interactions) lives with
# the scripts so it is always version-controlled alongside the code.
REFERENCE_DIR = SCRIPTS_DIR
RESULTS_DIR = EXPORTED_DIR / "results"

# --------------------------------------------------------------------------- #
# Defaults
# --------------------------------------------------------------------------- #
DEFAULTS: dict[str, Any] = {
    # Amplicon markers the pipeline understands. Add/remove freely; every
    # script iterates over these rather than a private hard-coded list.
    "markers": ["16s", "its", "18s"],

    # AMF is handled by filtering the standard ITS/18S tables for Glomeromycota
    # (scripts/filter_amf_table.py), so there are no separate AMF sub-markers.
    # The old its_amf / 18s_amf marker path was never completed (DADA2 never
    # produced those tables) and has been removed.
    "amf_markers": [],

    # Vendor primer pairs (forward, reverse). "I" (inosine) is normalised to N.
    "primers": {
        "16s": ["ACTCCTACGGGAGGCAGCAG", "GGACTACHVGGGTWTCTAAT"],   # 338F / 806R
        "its": ["CTTGGTCATTTAGAGGAAGTAA", "GCTGCGTTCTTCATCGATGC"],  # ITS1F / ITS2
        "18s": ["CGWTAACGAACGAGACCT", "AICCATTCAATCGGTAIT"],        # FF390 / FR1
        "amf": ["ATCAACTTTCGATGGTAGGATAGA", "GAACCCAAACACTTTGGTTTCC"],  # AML1 / AML2
    },

    # Cutadapt behaviour.
    "cutadapt": {"error_rate": 0.12, "times": 1, "discard_untrimmed": True},

    # DADA2 truncation / trimming per marker. These are dataset-dependent and
    # should be tuned from each demux.qzv quality profile; exposed here (and in
    # the Streamlit UI) so they are never buried inside the denoise script.
    "dada2": {
        "16s": {"trunc_len_f": 220, "trunc_len_r": 220, "trim_left_f": 0, "trim_left_r": 0},
        "its": {"trunc_len_f": 200, "trunc_len_r": 220, "trim_left_f": 0, "trim_left_r": 0},
        "18s": {"trunc_len_f": 220, "trunc_len_r": 220, "trim_left_f": 0, "trim_left_r": 0},
        "max_ee_f": 2, "max_ee_r": 2, "min_overlap": 20,
        "n_threads": 1, "n_reads_learn": 100000,
    },

    # Taxonomic classification.
    "classify": {"confidence": 0.8, "read_orientation": "auto", "n_jobs": 1},

    # How to match a classifier file on disk to a marker. Each marker lists
    # case-insensitive keyword sets; a classifier file matches if it contains
    # ANY keyword group (all tokens in the group present in the filename).
    "classifier_keywords": {
        "16s": [["silva", "16s"], ["silva", "138"], ["greengenes"]],
        "its": [["unite", "its"], ["unite"]],
        "18s": [["silva", "18s"], ["pr2"]],
        "its_amf": [["maarjam", "its"]],
        "18s_amf": [["maarjam", "18s"], ["maarjam"]],
    },

    # AMF taxa of interest for the AMF filter step (Glomeromycota lineage).
    "amf_target_taxa": [
        "Glomeromycota", "Glomeromycetes", "Glomerales", "Glomeraceae",
        "Rhizophagus", "Funneliformis", "Claroideoglomus", "Glomus",
        "Gigaspora", "Acaulospora", "Diversispora", "Scutellospora",
        "Paraglomus", "Archaeospora", "Ambispora", "Septoglomus",
    ],

    # Diversity settings.
    "diversity": {"sampling_depth_quantile": 0.10, "markers": ["16s", "its", "its_amf", "18s"]},

    # Metadata columns collected in the UI (first must be the sample id column).
    "metadata_columns": [
        "sample-id", "treatment", "timepoint", "location",
        "soil", "crop", "ph", "note",
    ],

    # Accepted FASTQ extensions.
    "fastq_extensions": [".fastq", ".fastq.gz", ".fq", ".fq.gz"],

    # Ordered pipeline steps: (script relative to scripts/, human label).
    "pipeline_steps": [
        ["fastqc.py", "FastQC on raw reads"],
        ["priming.py", "Primer trimming (cutadapt)"],
        ["fastqc.py", "FastQC on trimmed reads"],
        ["generate_metadata.py", "Combine metadata"],
        ["generate_manifest.py", "Generate QIIME2 manifest"],
        ["import_qiime2.py", "Import demux (paired-end)"],
        ["dada2_denoise.py", "DADA2 denoising"],
        ["taxonomic_classification.py", "Taxonomic classification"],
        ["export_all_qza.py", "Export BIOM / FASTA"],
        ["clean_asv_fasta.py", "Standardise ASV IDs"],
        ["filter_amf_table.py", "Filter AMF (Glomeromycota)"],
        ["run_qiime_diversity.py", "Alpha / beta diversity"],
        ["generate_asv_table.py", "Build ASV table"],
        ["heatmap.py", "Heatmap & barplots"],
        ["krona_chart.py", "Krona chart"],
        ["export_qiime2.py", "Export for PICRUSt2"],
        ["run_picrust2.py", "PICRUSt2 functional prediction (needs picrust2 env)"],
        ["generate_functional_order.py", "Functional order (TND)"],
        ["network_analysis.py", "Microbial interaction network"],
    ],

    # Microbial interaction / co-occurrence network settings.
    "network": {
        "min_prevalence": 0.30,        # taxon must be present in >= 30% of samples
        "min_abundance": 0.0,          # minimum total count to keep a taxon
        "clr_pseudocount": 1.0,        # added before centred-log-ratio transform
        "correlation_method": "spearman",  # spearman | pearson
        "fdr_alpha": 0.05,             # BH-adjusted significance threshold
        "min_abs_correlation": 0.60,   # |rho| edge threshold
        "taxonomic_rank": "Genus",     # rank at which interactions are reported
        "reference_db": "microbial_interactions.json",
    },
}

# --------------------------------------------------------------------------- #
# Loader
# --------------------------------------------------------------------------- #
def _load_local_overrides() -> dict[str, Any]:
    local = BASE_DIR / "config.local.yaml"
    if local.exists() and yaml is not None:
        with local.open("r", encoding="utf-8") as fh:
            return yaml.safe_load(fh) or {}
    return {}


def _deep_merge(base: dict, override: dict) -> dict:
    out = dict(base)
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out


_CONFIG = _deep_merge(DEFAULTS, _load_local_overrides())


def get(key: str, default: Any = None) -> Any:
    """Fetch a top-level config value, honouring PEANUT_<KEY> env overrides."""
    env_key = "PEANUT_" + key.upper()
    if env_key in os.environ:
        raw = os.environ[env_key]
        # allow comma lists for simple sequence overrides
        return [p.strip() for p in raw.split(",")] if "," in raw else raw
    return _CONFIG.get(key, default)


# Convenience accessors used across scripts ---------------------------------- #
MARKERS: list[str] = get("markers")
AMF_MARKERS: list[str] = get("amf_markers")
PRIMERS: dict[str, list[str]] = get("primers")
DADA2: dict[str, Any] = get("dada2")
CLASSIFY: dict[str, Any] = get("classify")
NETWORK: dict[str, Any] = get("network")
FASTQ_EXTS: tuple[str, ...] = tuple(get("fastq_extensions"))


def discover_classifier(marker: str) -> Path | None:
    """
    Find a classifier .qza on disk for the given marker by keyword matching,
    so renaming/upgrading a classifier file does not break the pipeline.
    Returns None if nothing suitable is present.
    """
    if not CLASSIFIER_DIR.exists():
        return None
    keyword_groups = get("classifier_keywords", {}).get(marker, [])
    candidates = sorted(CLASSIFIER_DIR.glob("*.qza"))
    for qza in candidates:
        name = qza.name.lower()
        for group in keyword_groups:
            if all(tok.lower() in name for tok in group):
                return qza
    return None


if __name__ == "__main__":
    print("Project root :", BASE_DIR)
    print("Markers      :", MARKERS)
    print("Classifiers found on disk:")
    for m in MARKERS + AMF_MARKERS:
        print(f"   {m:8s} -> {discover_classifier(m)}")