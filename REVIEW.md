# Code Review & Change Log — Microbiome Pipeline

This document records the review of the pipeline and the changes made to turn it
into a clean, evidence-based, general-purpose microbiome platform. Goals:
**(1)** run on GitHub + Streamlit, **(2)** provide an evidence base for microbial
interactions and abundance, and **(3)** remove hard-coding / make it dynamic.

Legend: ✅ done · 🟡 improved / optional follow-up · 🔴 was a real problem.

---

## 1. Evidence layers added (the core value)

The original pipeline produced taxonomy, diversity, heatmaps, Krona and PICRUSt2,
but **nothing interpreted the results**. Two evidence-based layers were built:

- **Microbial interaction network** — `scripts/network_analysis.py` +
  `scripts/microbial_interactions.json`. A compositionally-aware co-occurrence
  network (CLR → Spearman → p-values → Benjamini–Hochberg FDR → effect-size
  threshold) where every significant edge is cross-referenced to published
  literature and labelled `CONSISTENT` / `DISCORDANT` / `CONTEXT` / `NOVEL`.
  Validated numerically: ρ and p-values match `scipy.stats.spearmanr` exactly,
  and BH q-values match `statsmodels` to 1e-16.

- **Taxon insights** — `scripts/taxon_insights.py` +
  `scripts/taxon_insights.json`. For the taxa in a dataset it explains **why they
  are high or low** (drivers such as pH, salinity, nitrogen, organic carbon),
  what that **means**, candidate **interventions**, and **study citations**. This
  runs at any sample size (relative abundance, not correlation).

Both knowledge bases are **general** (soil/rhizosphere microbiomes broadly) and
extensible — add entries to the JSON with no code change. They include some
example organisms relevant to legume/peanut systems, but are not limited to them.

> ⚠️ Literature matches are **candidate annotations to verify** against the cited
> source — not automated causal proof for a specific dataset.

---

## 2. Hard-coding removed → centralised in `config.py` ✅

Values that were scattered/duplicated across scripts now live in one `config.py`
(overridable via `config.local.yaml` or `PEANUT_*` env vars):

| Was hard-coded | Now |
|---|---|
| Primer sequences | `config.primers` |
| DADA2 truncation lengths / max-EE / threads | `config.dada2` (editable in UI) |
| Classifier filenames | **discovered on disk** via `config.discover_classifier` |
| Marker lists | `config.markers` |
| Cutadapt options | `config.cutadapt` |
| AMF target taxa | `config.amf_target_taxa` |
| Pipeline step order | `config.pipeline_steps` |
| Network / analysis thresholds | `config.network` |

🟡 A few remaining scripts still use literal `../data/...` paths; they work but
could import from `config.py` for full consistency.

---

## 3. Real bugs fixed 🔴 → ✅

- **`filter_amf_table.py` was a byte-for-byte copy of `generate_asv_table.py`** —
  the AMF-filtering step did nothing. Rewritten to actually filter features whose
  taxonomy matches `config.amf_target_taxa`.
- **Classifier filenames didn't match the files on disk** — every sklearn
  classification would have skipped. Classifiers are now **discovered by
  keyword**, so present files are picked up automatically.
- **Heatmap sample detection was brittle** (`col.startswith("sample")`) — samples
  not named `sample*` were silently dropped. Replaced with numeric-column
  detection that makes no assumption about names.
- **FastQC case bug** — hard-coded `"ITS"` folder was never found on Linux; now
  reads markers from config and matches case-insensitively.

---

## 4. PICRUSt2 redesigned to actually work ✅

The original had three fragile, non-working approaches: a hand-built placement
chain (`01–05_*`, `picrust2_balance`), a `run_picrust2_stratified.py` that called
`sudo` to build a swapfile, and a Galaxy version with a hard-coded API key. All
replaced by a single `scripts/run_picrust2.py` that calls the official
`picrust2_pipeline.py` and writes a top-pathways summary. The 8 dead scripts +
`galaxy.py` were deleted.

---

## 5. Cleanup ✅

- Removed the orphan `ancom.py` (unwired, wrong path).
- Removed dead `its_amf` / `18s_amf` code branches; AMF is handled only by
  `filter_amf_table.py`.
- Removed the legacy Flask app (`main.py`, `templates/`) — replaced by Streamlit.
- Removed personal notes and an old duplicate `readme.txt`.
- Added a root `.gitignore` excluding FASTQ, `data/`, `uploads/`, `*.qza/qzv/biom`,
  `classifiers/`, vendored `picrust2/`, and `__pycache__/`.
- Split dependencies: `requirements.txt` (pip / app) + `environment.yml` (conda /
  QIIME2). The repo now tracks only ~28 files (~300 KB); multi-GB data stays local.
- No secrets in the repo — scripts read any credentials from the environment.

---

## 6. Architecture: Flask → Streamlit ✅

The old Flask app drove the pipeline via chained HTTP redirects (fragile; long
steps blocked requests). Replaced by `streamlit_app.py` with five sections:
Upload & metadata → Run pipeline (live logs) → Results → Interaction network →
Taxon insights. Verified to import without QIIME2, so the analysis pages deploy
on Streamlit Community Cloud.

---

## 7. Key limitation to remember (data, not code)

A co-occurrence interaction network correlates taxa **across samples**, so it
needs **≥ 3** samples and realistically **15–20+** for meaningful results. The
example dataset here has only 2 samples per marker, so the network correctly
reports "not enough samples" rather than fabricating correlations. Taxonomy,
composition and **taxon insights** work at any sample size; interaction
statistics need a larger sample set.

---

## 8. Optional follow-ups

1. Finish migrating remaining scripts to read paths from `config.py`.
2. Rename the conda env / env-var prefix from peanut-specific names to generic
   ones (`microbiome`, `MB_`) if desired — touches `config.py` + `environment.yml`.
3. Expand `microbial_interactions.json` and `taxon_insights.json` with taxa /
   drivers specific to your study system (no code change needed).
4. Move the git working copy off the network drive to local disk (network mounts
   corrupt the git index).
