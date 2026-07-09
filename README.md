# 🥜 Peanut Microbiome Pipeline

A 16S / ITS / 18S amplicon pipeline (QIIME2 → DADA2 → taxonomy → diversity →
PICRUSt2) with a **Streamlit** front-end and an **evidence-based microbial
interaction network** module. Everything is configuration-driven — no dataset
values are hard-coded.

> **Data note:** the 16S data study peanut rhizosphere response to intercropping
> with quinoa in saline-alkali soil; the ITS1 data study peanut soil fungal
> response to nitrogen fertiliser. The 16S and ITS datasets are **not** from the
> same samples.

---

## Highlights

- **Streamlit UI** (`streamlit_app.py`) — upload FASTQ, edit metadata, run the
  pipeline step-by-step with live logs, browse results, build the interaction
  network, and generate taxon insights. Replaces the old Flask app.
- **Evidence-based interaction network** (`scripts/network_analysis.py`) —
  compositionally-aware (CLR) co-occurrence network with p-values,
  Benjamini–Hochberg FDR control, effect-size thresholds, and every edge
  cross-referenced against a curated literature database
  (`scripts/microbial_interactions.json`).
- **Taxon insights** (`scripts/taxon_insights.py`) — for the taxa in your data,
  explains **why they may be high or low** (drivers), what that **means**,
  candidate **interventions**, and **study citations**, from a general
  literature-backed knowledge base (`scripts/taxon_insights.json`). Works with
  any microbiome dataset, not just peanut.
- **Nothing hard-coded** — markers, primers, DADA2 parameters, classifier files,
  pipeline steps and analysis thresholds all live in `config.py` and can be
  overridden with `config.local.yaml` or `PEANUT_*` environment variables.
- **Classifiers are discovered on disk** by keyword, so renaming a `.qza` file
  never silently breaks classification.
- **Cloud-ready** — the analysis pages (interaction network, taxon insights,
  results) deploy to Streamlit Community Cloud; the QIIME2 processing runs
  locally.

---

## Repository layout

```
config.py                      # central, dynamic configuration (edit here)
streamlit_app.py               # Streamlit UI (run this)
requirements.txt               # pip layer (app + analysis)
environment.yml                # conda layer (QIIME2 / bioinformatics CLIs)
REVIEW.md                      # code review notes & findings
.streamlit/config.toml         # app theme (used locally and on Streamlit Cloud)
scripts/
  network_analysis.py          # microbial interaction network + evidence base
  microbial_interactions.json  # curated literature interaction database
  taxon_insights.py            # why taxa are high/low + causes + solutions + evidence
  taxon_insights.json          # general literature-backed taxon knowledge base
  run_picrust2.py              # PICRUSt2 functional prediction (official CLI)
  priming.py                   # cutadapt primer trimming (config-driven)
  dada2_denoise.py             # DADA2 denoising (config-driven)
  taxonomic_classification.py  # classify-sklearn (classifiers auto-discovered)
  filter_amf_table.py          # AMF (Glomeromycota) filter
  ... (fastqc, manifest, import, export, diversity, heatmap, krona)
classifiers/                   # *.qza classifiers (git-ignored; large)
data/                          # pipeline outputs (git-ignored)
uploads/                       # uploaded FASTQ (git-ignored)
```

---

## Setup

### 1. Clone
```bash
git clone <your-repo-url>.git
cd peanut
```

### 2. Bioinformatics stack (conda)
The QIIME steps need QIIME2. Install the current **QIIME2 amplicon
distribution** from the official docs for your OS, then the extras:
```bash
conda env create -f environment.yml
conda activate peanut-microbiome
```
PICRUSt2 conflicts with QIIME2 — install it in its own env:
```bash
conda create -n picrust2 -c bioconda -c conda-forge picrust2
```

### 3. App / analysis layer (pip)
Already covered by the conda env, or install standalone:
```bash
pip install -r requirements.txt
```

### 4. Classifiers
Drop your trained classifiers into `classifiers/` (e.g. SILVA for 16S/18S,
UNITE for ITS, MaarjAM for AMF). They are matched to markers by keyword in
`config.py → classifier_keywords`, so exact filenames don't matter.

---

## Running

### Streamlit app
```bash
streamlit run streamlit_app.py
```
Then work through the sidebar: **Upload & metadata → Run pipeline → Results →
Interaction network → Taxon insights**.

### Interaction network from the command line
```bash
python scripts/network_analysis.py \
    --input data/exported/asv_tables_combined.xlsx \
    --sheet 16S --rank Genus --fdr-alpha 0.05
```
Outputs (edges with statistics + evidence, node hubs, summary JSON) are written
to `data/exported/results/network/`.

### Taxon insights from the command line
```bash
python scripts/taxon_insights.py \
    --input data/exported/asv_tables_combined.xlsx --sheet 16S
```
For each taxon it reports role, why-high / why-low drivers, implication,
interventions and citations. Unlike the network, this works with any number of
samples (it uses relative abundance, not correlation).

---

## The interaction network method (what makes it evidence-based)

Amplicon data are **compositional**, so naive correlation invents spurious
associations (Friedman & Alm 2012). The module therefore:

1. Collapses ASVs to a taxonomic rank (default Genus).
2. Filters by prevalence / abundance.
3. Applies a **centred-log-ratio (CLR)** transform.
4. Computes **Spearman** correlations with p-values for every taxon pair.
5. Controls the false discovery rate with **Benjamini–Hochberg**.
6. Keeps edges passing both an FDR (`q`) and an effect-size (`|ρ|`) threshold.
7. **Cross-references** each surviving edge against
   `scripts/microbial_interactions.json` and labels it
   `CONSISTENT` / `DISCORDANT` / `CONTEXT` / `NOVEL` relative to the literature,
   attaching the mechanism and citations.

Each reported interaction thus carries **both the maths** (ρ, p, q, co-prevalence)
**and the evidence** (interaction type, mechanism, references). The literature
matches are *candidate* annotations to verify against the primary source — see
the disclaimer in the JSON.

To extend the evidence base, add entries to
`scripts/microbial_interactions.json` (no code change needed).

---

## Taxon insights (why high/low → cause → solution → evidence)

`scripts/taxon_insights.py` reads any ASV table, computes each taxon's mean
relative abundance and prevalence, flags taxa as **High / Medium / Low**, and
annotates them from a general knowledge base (`scripts/taxon_insights.json`):

- **role** — what the taxon does ecologically,
- **why high / why low** — the drivers (pH, salinity, nitrogen, organic carbon,
  intercropping, moisture, …),
- **implication** — what it means agronomically/ecologically,
- **interventions** — candidate management actions,
- **evidence** — study citations.

The knowledge base is general (soil/rhizosphere microbiomes broadly) and
extensible — add taxa or factors to the JSON with no code change.

> ⚠️ These are general, literature-based **candidate explanations to interpret
> and verify** — not automated causal proof for a specific soil. Statistically
> proving which taxa differ between conditions still requires enough samples per
> group.

---

## Deploying to Streamlit Community Cloud

The analysis pages run online for free:

1. Push the repo to GitHub.
2. Go to https://share.streamlit.io → sign in with GitHub → **Create app**.
3. Repository = your repo, branch = `main`, main file = `streamlit_app.py`.
4. Deploy — installs `requirements.txt` and gives a public URL.

The **Interaction network**, **Taxon insights** and **Results** pages work on the
cloud. The **Run pipeline** page (QIIME2/DADA2/PICRUSt2) needs a local
bioinformatics environment and is not expected to run on the free cloud.

---

## Configuration

Edit `config.py`, or override without touching it:

- **YAML:** create `config.local.yaml` (git-ignored) mirroring keys in `DEFAULTS`.
- **Env vars:** e.g. `PEANUT_MARKERS="16s,its"`.

Key knobs: `markers`, `primers`, `dada2`, `classifier_keywords`,
`amf_target_taxa`, `pipeline_steps`, and `network` thresholds.

---

## GitHub notes / gotchas

- `data/`, `uploads/`, `raw Data/`, `*.qza`, and FASTQ files are **git-ignored**
  (large/binary). Share those out-of-band.
- **`picrust2/` ships its own `.git`.** Committed as-is it becomes an embedded
  repo and breaks `git add`. It is git-ignored here; install PICRUSt2 via conda.
  If you must keep it in-tree: `rm -rf picrust2/.git` first, then un-ignore it.
- Secrets: `scripts/galaxy.py` now reads `GALAXY_API_KEY` from the environment —
  never commit API keys.

---

## Credits & method references

- Friedman J, Alm EJ (2012) *Inferring Correlation Networks from Genomic Survey
  Data.* PLoS Comput Biol 8(9):e1002687. doi:10.1371/journal.pcbi.1002687
- Benjamini Y, Hochberg Y (1995) *Controlling the False Discovery Rate.*
  J R Stat Soc B 57(1):289–300.

See `scripts/microbial_interactions.json` for the biological interaction
citations.