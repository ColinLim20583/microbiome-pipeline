# 🧬 Microbiome Pipeline

A configurable **16S / ITS / 18S amplicon microbiome pipeline** (QIIME2 → DADA2 →
taxonomy → diversity → PICRUSt2) with a **Streamlit** interface and two
**evidence-based interpretation layers**:

- a **microbial interaction network** — which taxa co-occur or exclude each other, and
- **taxon insights** — why a taxon is high or low, what it means, and what to do,

with every claim backed by published literature.

Works with **any soil/rhizosphere microbiome dataset** — nothing is hard-coded to
a particular crop or study. (A peanut rhizosphere dataset is used as the built-in
example, but markers, primers, taxa and thresholds are all configurable.)

---

## Highlights

- **Streamlit UI** (`streamlit_app.py`) — upload FASTQ, edit metadata, run the
  pipeline with live logs, browse results, build the interaction network, and
  generate taxon insights.
- **Evidence-based interaction network** (`scripts/network_analysis.py`) —
  compositionally-aware (CLR) co-occurrence network with p-values,
  Benjamini–Hochberg FDR control, effect-size thresholds, and every edge
  cross-referenced against a literature database (`scripts/microbial_interactions.json`).
- **Taxon insights** (`scripts/taxon_insights.py`) — for the taxa in your data,
  explains **why they may be high or low** (drivers), what that **means**,
  candidate **interventions**, and **study citations**, from a general
  literature-backed knowledge base (`scripts/taxon_insights.json`).
- **Nothing hard-coded** — markers, primers, DADA2 parameters, classifier files,
  pipeline steps and analysis thresholds all live in `config.py` and can be
  overridden with `config.local.yaml` or `PEANUT_*` environment variables.
- **Classifiers discovered on disk** by keyword, so renaming a `.qza` file never
  silently breaks classification.
- **Cloud-ready** — the analysis pages deploy to Streamlit Community Cloud; the
  QIIME2 processing runs locally.

---

## Repository layout

```
config.py                      # central, dynamic configuration (edit here)
streamlit_app.py               # Streamlit UI (run this)
requirements.txt               # pip layer (app + analysis)
environment.yml                # conda layer (QIIME2 / bioinformatics CLIs)
.streamlit/config.toml         # app theme (used locally and on Streamlit Cloud)
README.md                      # this file
REVIEW.md                      # code review notes & change log
scripts/
  network_analysis.py          # microbial interaction network + evidence base
  microbial_interactions.json  # literature-backed interaction database
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

## What the pipeline produces

| Stage | Question it answers | Output |
|---|---|---|
| Taxonomy / ASV table | What microbes are present, and how much? | ASV tables, Krona charts |
| Diversity | How diverse / different are the samples? | alpha & beta diversity |
| PICRUSt2 | What functions might the community perform? | KO / EC / pathway predictions |
| **Interaction network** | Which taxa interact (co-occur / exclude)? | edges with ρ, p, q + literature |
| **Taxon insights** | Why is a taxon high/low, and what can be done? | drivers, interventions, citations |

---

## Setup

### 1. Clone
```bash
git clone https://github.com/<you>/microbiome-pipeline.git
cd microbiome-pipeline
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
Sidebar workflow: **Upload & metadata → Run pipeline → Results →
Interaction network → Taxon insights**.

### Interaction network (command line)
```bash
python scripts/network_analysis.py \
    --input data/exported/asv_tables_combined.xlsx \
    --sheet 16S --rank Genus --fdr-alpha 0.05
```
Outputs (edges with statistics + evidence, node hubs, summary JSON) go to
`data/exported/results/network/`.

### Taxon insights (command line)
```bash
python scripts/taxon_insights.py \
    --input data/exported/asv_tables_combined.xlsx --sheet 16S
```
Reports role, why-high / why-low drivers, implication, interventions and
citations per taxon. Works with any number of samples (relative abundance).

### PICRUSt2 (functional prediction)
```bash
conda activate picrust2
python scripts/run_picrust2.py
```
Runs the official `picrust2_pipeline.py` on the exported rep-seqs + feature table
and writes a top-pathways summary.

---

## The interaction network method (why it's evidence-based)

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
   attaching mechanism and citations.

Each reported interaction carries **both the maths** (ρ, p, q, co-prevalence)
**and the evidence** (interaction type, mechanism, references). Literature
matches are *candidate* annotations to verify against the primary source.

> **Sample-size note:** a co-occurrence network correlates taxa across samples,
> so it needs **≥ 3** samples and realistically **15–20+** for meaningful
> results. With fewer, the app reports "not enough samples" rather than
> fabricating correlations.

To extend the interaction database, add entries to
`scripts/microbial_interactions.json` (no code change needed).

---

## Taxon insights (why high/low → cause → solution → evidence)

`scripts/taxon_insights.py` reads any ASV table, computes each taxon's mean
relative abundance and prevalence, flags taxa as **High / Medium / Low**, and
annotates them from a general knowledge base (`scripts/taxon_insights.json`):

- **role** — what the taxon does ecologically,
- **why high / why low** — the drivers (pH, salinity, nitrogen, organic carbon,
  intercropping, moisture, …),
- **implication** — what it means agronomically / ecologically,
- **interventions** — candidate management actions,
- **evidence** — study citations.

The knowledge base is general (soil/rhizosphere microbiomes broadly) and
extensible — add taxa or factors to the JSON with no code change.

> ⚠️ These are general, literature-based **candidate explanations to interpret
> and verify** — not automated causal proof for a specific soil. Statistically
> proving which taxa differ between conditions still requires enough samples per
> group.

---

## Configuration

Edit `config.py`, or override without touching it:

- **YAML:** create `config.local.yaml` (git-ignored) mirroring keys in `DEFAULTS`.
- **Env vars:** e.g. `PEANUT_MARKERS="16s,its"`.

Key knobs: `markers`, `primers`, `dada2`, `classifier_keywords`,
`amf_target_taxa`, `pipeline_steps`, and `network` thresholds.

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

## GitHub notes / gotchas

- `data/`, `uploads/`, `raw Data/`, `*.qza`, FASTQ and `picrust2/` are
  **git-ignored** (large/binary). Only code, config and docs are committed
  (~300 KB); raw data stays local.
- Install PICRUSt2 via conda in its own environment (it conflicts with QIIME2).
- Never commit API keys — pipeline scripts read any credentials from the
  environment.

---

## Credits & method references

- Friedman J, Alm EJ (2012) *Inferring Correlation Networks from Genomic Survey
  Data.* PLoS Comput Biol 8(9):e1002687. doi:10.1371/journal.pcbi.1002687
- Benjamini Y, Hochberg Y (1995) *Controlling the False Discovery Rate.*
  J R Stat Soc B 57(1):289–300.

See `scripts/microbial_interactions.json` and `scripts/taxon_insights.json` for
the biological interaction and taxon-driver citations.
