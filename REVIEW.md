# Code Review — Peanut Microbiome Pipeline

Review of the pipeline against your three goals: **(1)** build/run on GitHub +
Streamlit, **(2)** provide an *evidence base for microbial interaction*, and
**(3)** remove hard-coding / make it dynamic. FASTQ handling was treated as
out-of-scope per your note.

Legend: ✅ fixed in this pass · 🟡 improved, needs your input · 🔴 flagged, action
recommended.

---

## 1. Biggest gap: there was no microbial interaction analysis at all ✅

The 17-step pipeline produced taxonomy, diversity, heatmaps, Krona and PICRUSt2,
but **nothing computed or evidenced microbial interactions**. That is the core
of your request, so it was built from scratch:

- **`scripts/network_analysis.py`** — a compositionally-aware co-occurrence
  network with real statistics and a literature evidence base.
- **`scripts/microbial_interactions.json`** — a curated, citable reference DB of
  peanut/rhizosphere interactions (Bradyrhizobium–peanut N-fixation, AMF
  symbiosis, mycorrhiza helper bacteria, Trichoderma/Bacillus vs *Aspergillus
  flavus*, Pseudomonas vs *Fusarium*, …).

**How it makes interactions evidence-based (all three senses you asked for):**

| Your ask | How it's delivered |
|---|---|
| Statistical rigor | CLR transform → Spearman → p-values → **Benjamini–Hochberg FDR** → effect-size threshold. Avoids the spurious-correlation trap of compositional data (Friedman & Alm 2012). |
| Literature citations | Every significant edge is cross-referenced to the JSON DB and labelled `CONSISTENT` / `DISCORDANT` / `CONTEXT` / `NOVEL`, with mechanism + citations. |
| Show the math/data | Output table surfaces ρ, p, q, co-prevalence, sample count per edge; a "Method & references" panel exposes the parameters used. |

Validated numerically: ρ and p-values match `scipy.stats.spearmanr` exactly, and
BH q-values match `statsmodels` to 1e-16. A planted positive symbiosis
(ρ≈+0.93) and negative antagonism (ρ≈−0.78) were both correctly detected and
labelled `SUPPORTED/CONSISTENT`.

> ⚠️ The literature matches are *candidate* annotations at genus level — a name
> match flags an edge for you to verify against the cited paper; it is not proof
> that the specific ASVs interact. This caveat is stated in the JSON and UI.

---

## 2. Hard-coding removed → centralised in `config.py` ✅ / 🟡

Values were scattered and duplicated across scripts. They now live in one
`config.py` (overridable via `config.local.yaml` or `PEANUT_*` env vars).

| Was hard-coded | Where | Now |
|---|---|---|
| Primer sequences | `priming.py` | `config.primers` ✅ |
| DADA2 truncation lengths, max-EE, threads | `dada2_denoise.py` | `config.dada2` ✅ (editable in UI) |
| Classifier filenames | `taxonomic_classification.py` | **discovered on disk** via `config.discover_classifier` ✅ |
| Marker list `["16s","its","18s"]` | many scripts | `config.markers` ✅ (refactored the key scripts) |
| Cutadapt error rate / options | `priming.py` | `config.cutadapt` ✅ |
| AMF target taxa | (missing) | `config.amf_target_taxa` ✅ |
| Pipeline step order | `main.py` | `config.pipeline_steps` ✅ |
| Network thresholds | (n/a) | `config.network` ✅ |

🟡 Remaining scripts (`import_qiime2.py`, `export_all_qza.py`, `export_qiime2.py`,
`run_qiime_diversity.py`, `generate_asv_table.py`, `krona_chart.py`, the
`generate_functional_*`/`picrust2` helpers, and the `01–05_*` placement scripts)
still contain literal `../data/...` paths and marker lists. They work, but for
full consistency they should import from `config.py` the same way — a
mechanical follow-up. Say the word and I'll convert them too.

---

## 3. Real bugs found

### 3.1 `filter_amf_table.py` was a copy-paste of `generate_asv_table.py` 🔴→✅
The two files were **byte-for-byte identical**. The AMF-filtering step (main.py
step 11) therefore did nothing — it just rebuilt the ASV table. Rewritten to
actually filter features whose taxonomy matches `config.amf_target_taxa`
(Glomeromycota lineage) and emit per-marker + combined AMF tables.

### 3.2 Classifier filenames didn't match the files on disk 🔴→✅
`taxonomic_classification.py` referenced e.g.
`silva-138.1-16s-v3v4-nb-classifier.qza`, `unite-its1-nb-classifier.qza`,
`silva-138.1-18s-V7V8-nb-classifier.qza` — **none of which exist** in
`classifiers/` (which actually holds `silva-138-99-nb-classifier.qza`,
`maarjam_its_classifier.qza`, `maarjam_18s_classifier.qza`). Every sklearn
classification would have skipped with "missing classifier". Now classifiers are
**discovered by keyword**, so the present files are picked up automatically.

### 3.3 Heatmap sample detection was brittle 🔴→✅
`heatmap.py` selected sample columns with
`col.lower().startswith("sample")`. Any sample not literally named `sample*`
(e.g. `P1`, `T3_rep2`) would be silently dropped and the heatmap would be empty
or wrong — directly against your "must be dynamic" requirement. Replaced with
numeric-column detection that ignores taxonomy/ID columns and makes no
assumption about names.

### 3.4 `18s_amf` reused the wrong rep-seqs 🟡
In the original config, `18s_amf` pointed at `18s_rep_seqs.qza` with a comment
admitting it may be wrong. Kept the mapping explicit in `REP_SEQS_FOR` so it's
obvious and easy to correct if you do have separate AMF rep-seqs.

---

## 4. Security & secrets 🔴→✅

- `scripts/galaxy.py` committed a **live-looking Galaxy API key**
  (`26e8755e…`) and `http://localhost:8080`. Now reads `GALAXY_URL` /
  `GALAXY_API_KEY` from the environment and errors if unset. **Rotate that key**
  if it was ever real.
- `main.py` (Flask) used `app.secret_key = "supersecret"`. The Streamlit app
  doesn't need it; if you keep Flask, move it to an env var.

---

## 5. GitHub readiness 🔴→✅

- **No root `.gitignore`** existed (only a PyCharm one under `.idea/`). Added a
  proper one that excludes FASTQ, `data/`, `uploads/`, `*.qza/qzv/biom`, and
  `__pycache__/`.
- **`picrust2/` contains its own `.git`** (a full nested repo, ~hundreds of
  files). Committed as-is it becomes an embedded repo/gitlink and breaks
  `git add`. It's now git-ignored; install PICRUSt2 via conda instead. If you
  must vendor it, `rm -rf picrust2/.git` first.
- A committed `scripts/__pycache__/*.pyc` was present — now ignored.
- **`requirements.txt` listed non-pip packages** (`qiime2`, `krona`, `fastqc`).
  These are conda/bioinformatics tools and `pip install` would fail. Split into
  a realistic `requirements.txt` (pip layer) + `environment.yml` (conda layer).

---

## 6. Architecture: Flask → Streamlit ✅

The Flask app (`main.py`) drove the pipeline by chaining HTTP redirects, one per
step, running scripts via `subprocess`. That's fragile (a long DADA2 step blocks
a request; errors dump raw HTML). Replaced with `streamlit_app.py`:

- dynamic marker tabs and metadata editor (grid),
- step selection + live logs + progress bar,
- results browser, and the interaction-network page.

`main.py` is left in place for reference; you can delete it once you're happy
with the Streamlit version.

---

## 7. Suggested next steps

1. Point me at the remaining scripts (§2 🟡) to finish the config migration.
2. Confirm the `18s_amf` rep-seqs mapping (§3.4).
3. Rotate the Galaxy API key (§4).
4. Decide whether to vendor or conda-install PICRUSt2 (§5), then push.
5. Expand `microbial_interactions.json` with any interactions specific to your
   study system — it's designed to grow without code changes.

---

## 8. Update log (work completed after the initial review)

- **PICRUSt2 redesigned to actually work ✅** — replaced the fragile hand-built
  chain (`01–05_*`, `picrust2_balance`, `run_picrust2_stratified` with its
  `sudo` swapfile, and the Galaxy version with a second exposed API key) with a
  single `scripts/run_picrust2.py` that calls the official `picrust2_pipeline.py`
  and writes a top-pathways summary. Deleted the 8 dead scripts + `galaxy.py`.
- **Removed differential-abundance orphan `ancom.py`** (unwired, wrong path).
- **Cleaned dead `its_amf` / `18s_amf` code branches** across export/diversity/
  ASV-table/krona scripts; AMF is now handled only by `filter_amf_table.py`.
- **FastQC case bug fixed** — markers read from config; ITS no longer skipped on
  Linux.
- **New: Taxon Insights ✅** — `scripts/taxon_insights.py` +
  `scripts/taxon_insights.json` explain **why each taxon is high/low, what it
  means, what to do, and the study evidence**, from a general (any-microbiome)
  literature-backed knowledge base. Exposed as Streamlit page 5. Runs on any
  sample count (relative abundance, not correlation).
- **Repo trimmed & cloud-ready** — non-essential notes/old-code removed from
  tracking; `.streamlit/config.toml` theme added; app verified to import without
  QIIME2 so it deploys on Streamlit Community Cloud.

### Key limitation to remember (data, not code)
The current dataset has **only 2 samples per marker** (`FenceRow`, `WorkArea`).
The interaction network needs ≥3 (ideally 15–20+) samples, so it correctly
reports "not enough samples" until more are added. Taxon insights and taxonomy/
composition views work now; interaction statistics need a larger sample set.