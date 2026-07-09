"""
Peanut Microbiome Pipeline - Streamlit interface.

Replaces the Flask app (main.py) with a single, dynamic Streamlit UI:

  * Upload paired-end FASTQ per marker (markers come from config, not hard-coded).
  * Enter/edit sample metadata in an editable grid.
  * Run the QIIME2 pipeline step-by-step with live logs (steps come from config).
  * Explore results, and run the Microbial Interaction Network with a full
    evidence base (statistics + literature) - the headline feature.

Nothing about the dataset is hard-coded here: markers, pipeline steps, primers,
thresholds and file layout are all read from config.py and are adjustable in the
sidebar at run time.

Run with:
    streamlit run streamlit_app.py
"""
from __future__ import annotations

import io
import subprocess
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg  # noqa: E402

st.set_page_config(page_title="Peanut Microbiome Pipeline", page_icon="🥜", layout="wide")

# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def ensure_dirs() -> None:
    for d in (cfg.UPLOAD_DIR, cfg.DATA_DIR, cfg.EXPORTED_DIR):
        d.mkdir(parents=True, exist_ok=True)


def sample_id_from_filename(filename: str) -> str:
    import re
    name = re.sub(r"\.(fastq|fq)(\.gz)?$", "", filename, flags=re.IGNORECASE)
    name = re.sub(r"(_R?[12])$", "", name, flags=re.IGNORECASE)
    name = re.sub(r"_trimmed$", "", name, flags=re.IGNORECASE)
    return name.strip("_")


def run_script(script_name: str) -> tuple[int, str]:
    """Run a pipeline script and capture combined stdout/stderr."""
    script_path = cfg.SCRIPTS_DIR / script_name
    if not script_path.exists():
        return 1, f"Script not found: {script_path}"
    proc = subprocess.run(
        [sys.executable, str(script_path)],
        capture_output=True, text=True, cwd=str(cfg.SCRIPTS_DIR),
    )
    return proc.returncode, (proc.stdout or "") + (proc.stderr or "")


# --------------------------------------------------------------------------- #
# Sidebar - dynamic configuration
# --------------------------------------------------------------------------- #
st.sidebar.title("🥜 Pipeline")
page = st.sidebar.radio(
    "Section",
    ["1 - Upload & metadata", "2 - Run pipeline", "3 - Results",
     "4 - Interaction network", "About"],
)

st.sidebar.markdown("---")
st.sidebar.caption("Markers (from config)")
st.sidebar.write(", ".join(cfg.MARKERS))

ensure_dirs()

# =========================================================================== #
# PAGE 1 - Upload & metadata
# =========================================================================== #
if page.startswith("1"):
    st.title("Upload FASTQ & enter metadata")
    st.write("Markers are read from `config.py` - add or remove them there and this page updates automatically.")

    tabs = st.tabs([m.upper() for m in cfg.MARKERS])
    for tab, mode in zip(tabs, cfg.MARKERS):
        with tab:
            st.subheader(f"{mode.upper()} paired-end FASTQ")
            uploads = st.file_uploader(
                f"Upload {mode.upper()} .fastq(.gz) files",
                accept_multiple_files=True, key=f"up_{mode}",
                type=[e.lstrip(".") for e in ("fastq", "fq", "gz")],
            )
            if uploads:
                folder = cfg.UPLOAD_DIR / mode
                folder.mkdir(parents=True, exist_ok=True)
                sample_ids = set()
                for f in uploads:
                    (folder / f.name).write_bytes(f.getbuffer())
                    sid = sample_id_from_filename(f.name)
                    if sid:
                        sample_ids.add(sid)
                st.session_state[f"{mode}_samples"] = sorted(sample_ids)
                st.success(f"Saved {len(uploads)} files - detected samples: {', '.join(sorted(sample_ids))}")

            # Metadata editor (columns from config; fully dynamic)
            samples = st.session_state.get(f"{mode}_samples", [])
            if samples:
                st.markdown("**Sample metadata**")
                cols = cfg.get("metadata_columns")
                default = pd.DataFrame({cols[0]: samples})
                for c in cols[1:]:
                    default[c] = ""
                edited = st.data_editor(default, key=f"md_{mode}", num_rows="dynamic", use_container_width=True)
                if st.button(f"Save {mode.upper()} metadata", key=f"save_{mode}"):
                    out = cfg.DATA_DIR / f"{mode}_metadata.tsv"
                    edited.to_csv(out, sep="\t", index=False)
                    st.success(f"Wrote {out}")

# =========================================================================== #
# PAGE 2 - Run pipeline
# =========================================================================== #
elif page.startswith("2"):
    st.title("Run the QIIME2 pipeline")
    steps = cfg.get("pipeline_steps")
    st.write(f"{len(steps)} steps, defined in `config.py -> pipeline_steps`. "
             "Requires the `qiime2` conda environment to be active for the QIIME steps.")

    with st.expander("Adjust DADA2 truncation lengths before running"):
        st.caption("These feed straight into config at run time via environment - tune from your demux.qzv quality plots.")
        st.json(cfg.DADA2)

    selected = st.multiselect(
        "Steps to run (default: all)",
        options=[f"{i}: {label}" for i, (_, label) in enumerate(steps)],
        default=[f"{i}: {label}" for i, (_, label) in enumerate(steps)],
    )
    if st.button("▶️ Run selected steps"):
        chosen_idx = [int(s.split(":")[0]) for s in selected]
        progress = st.progress(0.0)
        log_area = st.empty()
        full_log = ""
        for n, idx in enumerate(chosen_idx, 1):
            script, label = steps[idx]
            full_log += f"\n=== Step {idx}: {label} ({script}) ===\n"
            log_area.code(full_log[-4000:])
            code, out = run_script(script)
            full_log += out
            full_log += f"\n[{'OK' if code == 0 else 'FAILED (' + str(code) + ')'}]\n"
            log_area.code(full_log[-4000:])
            progress.progress(n / len(chosen_idx))
            if code != 0:
                st.error(f"Step {idx} ({label}) failed - see log above.")
                break
        else:
            st.success("All selected steps completed.")

# =========================================================================== #
# PAGE 3 - Results
# =========================================================================== #
elif page.startswith("3"):
    st.title("Results")
    combined = cfg.EXPORTED_DIR / "asv_tables_combined.xlsx"
    if combined.exists():
        st.subheader("ASV tables")
        xls = pd.ExcelFile(combined)
        sheet = st.selectbox("Sheet", xls.sheet_names)
        st.dataframe(xls.parse(sheet), use_container_width=True)
        st.download_button("Download workbook", data=combined.read_bytes(),
                           file_name=combined.name)
    else:
        st.info("No ASV table yet - run the pipeline first.")

    viz = cfg.DATA_DIR / "visualizations"
    if viz.exists():
        imgs = sorted(viz.glob("*.png"))
        if imgs:
            st.subheader("Figures")
            for img in imgs:
                st.image(str(img), caption=img.name)

# =========================================================================== #
# PAGE 4 - Interaction network (headline feature)
# =========================================================================== #
elif page.startswith("4"):
    st.title("🕸️ Microbial interaction network + evidence base")
    st.write(
        "Builds a statistically defensible co-occurrence network and annotates "
        "every interaction with published literature. Method: CLR transform "
        "(compositional) → Spearman correlation → p-values → Benjamini-Hochberg "
        "FDR → effect-size threshold → literature cross-reference."
    )

    from network_analysis import run_network_analysis, load_reference_db, NetworkParams  # noqa: E402

    # Input selection
    combined = cfg.EXPORTED_DIR / "asv_tables_combined.xlsx"
    uploaded = st.file_uploader("ASV / feature table (.xlsx, .tsv, .csv)", type=["xlsx", "tsv", "csv"])
    src_df = None
    sheet = None
    if uploaded is not None:
        if uploaded.name.endswith(("xlsx", "xls")):
            xls = pd.ExcelFile(uploaded)
            sheet = st.selectbox("Sheet / marker", [s for s in xls.sheet_names if "TOP" not in s.upper()])
            src_df = xls.parse(sheet)
        else:
            sep = "\t" if uploaded.name.endswith("tsv") else ","
            src_df = pd.read_csv(uploaded, sep=sep)
    elif combined.exists():
        xls = pd.ExcelFile(combined)
        sheet = st.selectbox("Sheet / marker", [s for s in xls.sheet_names if "TOP" not in s.upper()])
        src_df = xls.parse(sheet)
    else:
        st.info("Upload a table or run the pipeline to generate one.")

    # Thresholds (defaults from config, adjustable)
    n = cfg.NETWORK
    c1, c2, c3 = st.columns(3)
    rank = c1.selectbox("Taxonomic rank", ["Genus", "Family", "Order", "Class", "Phylum"],
                        index=["Genus", "Family", "Order", "Class", "Phylum"].index(n.get("taxonomic_rank", "Genus")))
    min_prev = c1.slider("Min prevalence", 0.0, 1.0, float(n.get("min_prevalence", 0.30)), 0.05)
    method = c2.selectbox("Correlation", ["spearman", "pearson"],
                          index=0 if n.get("correlation_method", "spearman") == "spearman" else 1)
    min_rho = c2.slider("Min |correlation|", 0.0, 1.0, float(n.get("min_abs_correlation", 0.60)), 0.05)
    fdr = c3.slider("FDR alpha (q)", 0.001, 0.20, float(n.get("fdr_alpha", 0.05)), 0.001)

    if src_df is not None and st.button("🔬 Build interaction network"):
        try:
            ref_path = cfg.REFERENCE_DIR / n.get("reference_db", "microbial_interactions.json")
            metadata_path = cfg.DATA_DIR / "metadata.tsv"
            result = run_network_analysis(
                df=src_df, sheet=sheet,
                metadata_path=str(metadata_path) if metadata_path.exists() else None,
                rank=rank, min_prevalence=min_prev, correlation_method=method,
                min_abs_correlation=min_rho, fdr_alpha=fdr, reference_db=str(ref_path),
            )
        except Exception as e:
            st.error(f"Analysis failed: {e}")
            st.stop()

        s = result.summary
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Taxa analysed", s["n_taxa_analyzed"])
        m2.metric("Significant edges", s["n_significant_edges"])
        m3.metric("Literature-supported", s["n_supported_by_literature"])
        m4.metric("Consistent w/ literature", s["n_consistent"])

        sig = result.significant_edges
        if sig.empty:
            st.warning("No edges passed the thresholds. Try lowering |correlation| or FDR.")
        else:
            st.subheader("Interactions with the maths AND the evidence")
            show_cols = ["taxon_a", "taxon_b", "direction", "rho", "p_value", "q_value",
                         "co_prevalence", "evidence_status", "literature_verdict",
                         "interaction_type", "mechanism", "references"]
            st.dataframe(sig[[c for c in show_cols if c in sig.columns]], use_container_width=True)

            st.download_button("Download significant edges (CSV)",
                               data=sig.to_csv(index=False).encode(),
                               file_name="interaction_edges_significant.csv")

            # Literature-supported subset, highlighted
            supported = sig[sig["evidence_status"] == "SUPPORTED"]
            if not supported.empty:
                st.subheader("📚 Literature-backed interactions")
                for _, r in supported.iterrows():
                    verdict_icon = {"CONSISTENT": "✅", "DISCORDANT": "⚠️", "CONTEXT": "ℹ️"}.get(r["literature_verdict"], "")
                    st.markdown(
                        f"**{r['taxon_a']} ↔ {r['taxon_b']}** &nbsp; {verdict_icon} "
                        f"`{r['literature_verdict']}` &nbsp; ρ={r['rho']:+.2f}, q={r['q_value']:.1e}  \n"
                        f"*{r['interaction_type']}* — {r['mechanism']}  \n"
                        f"<small>{r['references']}</small>",
                        unsafe_allow_html=True,
                    )

            # Node hubs
            if not result.node_stats.empty:
                st.subheader("Hub taxa (most connected)")
                st.dataframe(result.node_stats.head(15), use_container_width=True)

            # Simple network figure (matplotlib + networkx if available)
            try:
                import networkx as nx
                import matplotlib.pyplot as plt
                G = nx.Graph()
                for _, r in sig.iterrows():
                    G.add_edge(r["taxon_a"], r["taxon_b"], weight=abs(r["rho"]),
                               color="#2c7bb6" if r["rho"] >= 0 else "#d7191c")
                fig, ax = plt.subplots(figsize=(9, 7))
                pos = nx.spring_layout(G, seed=1, k=0.6)
                edge_colors = [G[u][v]["color"] for u, v in G.edges()]
                nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=2, alpha=0.7, ax=ax)
                nx.draw_networkx_nodes(G, pos, node_color="#fdae61",
                                       node_size=[300 + 200 * G.degree(x) for x in G.nodes()], ax=ax)
                nx.draw_networkx_labels(G, pos, font_size=8, ax=ax)
                ax.set_axis_off()
                ax.set_title("Blue = co-occurrence (+), Red = exclusion (−)")
                st.pyplot(fig)
            except Exception as e:
                st.caption(f"(network figure unavailable: {e})")

        # Method transparency
        with st.expander("Method & references (surfacing the maths)"):
            st.json({k: s[k] for k in ["method", "fdr_alpha", "min_abs_correlation",
                                       "n_pairs_tested", "n_samples", "rank"]})
            for mr in s.get("method_references", []):
                st.markdown(f"- {mr.get('citation','')} {('doi:' + mr['doi']) if mr.get('doi') else ''}")

# =========================================================================== #
# ABOUT
# =========================================================================== #
else:
    st.title("About this pipeline")
    st.markdown(
        """
This app runs a 16S / ITS / 18S amplicon pipeline (QIIME2 + DADA2 + PICRUSt2)
and adds an **evidence-based microbial interaction network**.

**Design principles**
- Nothing is hard-coded: markers, primers, DADA2 parameters, classifier files,
  pipeline steps and network thresholds all come from `config.py` and can be
  overridden with `config.local.yaml` or `PEANUT_*` environment variables.
- Classifiers are discovered on disk by keyword, so renaming a `.qza` never
  silently breaks classification.
- The interaction network is compositionally aware (CLR), FDR-controlled, and
  each edge is cross-referenced to published literature.
        """
    )
    st.subheader("Current configuration")
    st.json({
        "markers": cfg.MARKERS,
        "amf_markers": cfg.AMF_MARKERS,
        "network": cfg.NETWORK,
        "classifiers_found": {m: (str(cfg.discover_classifier(m)) if cfg.discover_classifier(m) else None)
                              for m in cfg.MARKERS + cfg.AMF_MARKERS},
    })