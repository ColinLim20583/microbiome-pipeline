"""
Microbial interaction (co-occurrence) network analysis with an evidence base.

This module turns an ASV / feature abundance table into a statistically defensible
microbial interaction network and annotates each detected interaction with
published literature, so every edge carries BOTH the maths and the evidence.

Why it is built the way it is
-----------------------------
Amplicon abundance data are *compositional* (only relative fractions are
observed). Naively correlating raw counts produces spurious associations
(Friedman & Alm 2012, PLoS Comput Biol 8:e1002687). We therefore:

  1. Collapse ASVs to a chosen taxonomic rank (default Genus) - interactions are
     interpretable and matchable to literature at genus level.
  2. Filter by prevalence and abundance so rare taxa do not create noise edges.
  3. Apply a centred-log-ratio (CLR) transform (compositionally aware).
  4. Compute Spearman (rank) correlations between every taxon pair, with p-values.
  5. Control the false discovery rate across all pairs with Benjamini-Hochberg.
  6. Keep only edges that pass both an FDR threshold (q) and an effect-size
     threshold (|rho|).
  7. Cross-reference each surviving edge against a curated, literature-backed
     reference database (microbial_interactions.json) and label it CONSISTENT /
     DISCORDANT / NOVEL relative to what is published.

Nothing here is hard-coded to the peanut dataset: sample columns, taxa and
thresholds are all detected/configured dynamically. Thresholds come from
config.py (and can be overridden per call or from the Streamlit UI).

Usage
-----
CLI:
    python network_analysis.py                      # uses config defaults
    python network_analysis.py --input path.xlsx --rank Genus --sheet 16S

Importable:
    from network_analysis import run_network_analysis
    result = run_network_analysis(df=my_feature_table, rank="Genus")
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy import stats

# Make config importable whether run from repo root or from scripts/
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
try:
    import config as cfg
except Exception:  # pragma: no cover - allow standalone use without config.py
    cfg = None

TAX_RANKS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
NON_SAMPLE_COLUMNS = set(TAX_RANKS) | {"ASV ID", "ASV_ID", "Feature ID", "FeatureID", "feature-id", "OTU ID"}


# --------------------------------------------------------------------------- #
# Configuration container
# --------------------------------------------------------------------------- #
@dataclass
class NetworkParams:
    rank: str = "Genus"
    min_prevalence: float = 0.30
    min_abundance: float = 0.0
    clr_pseudocount: float = 1.0
    correlation_method: str = "spearman"   # "spearman" | "pearson"
    fdr_alpha: float = 0.05
    min_abs_correlation: float = 0.60
    reference_db: Path | None = None

    @classmethod
    def from_config(cls) -> "NetworkParams":
        if cfg is None:
            return cls()
        n = cfg.NETWORK
        ref = cfg.REFERENCE_DIR / n.get("reference_db", "microbial_interactions.json")
        return cls(
            rank=n.get("taxonomic_rank", "Genus"),
            min_prevalence=float(n.get("min_prevalence", 0.30)),
            min_abundance=float(n.get("min_abundance", 0.0)),
            clr_pseudocount=float(n.get("clr_pseudocount", 1.0)),
            correlation_method=n.get("correlation_method", "spearman"),
            fdr_alpha=float(n.get("fdr_alpha", 0.05)),
            min_abs_correlation=float(n.get("min_abs_correlation", 0.60)),
            reference_db=ref,
        )


@dataclass
class NetworkResult:
    params: NetworkParams
    taxa_table: pd.DataFrame                 # taxon x sample (filtered, raw counts)
    edges: pd.DataFrame                      # all tested pairs with stats
    significant_edges: pd.DataFrame          # edges passing FDR + effect size
    node_stats: pd.DataFrame                 # per-taxon degree / connectivity
    summary: dict[str, Any] = field(default_factory=dict)


# --------------------------------------------------------------------------- #
# 1. Load & detect
# --------------------------------------------------------------------------- #
def load_feature_table(path: str | Path, sheet: str | None = None) -> pd.DataFrame:
    """
    Load an abundance table from .xlsx (a specific sheet or the first taxonomic
    sheet) or a .tsv/.csv. Returns the raw DataFrame unchanged.
    """
    path = Path(path)
    if path.suffix.lower() in {".xlsx", ".xls"}:
        xls = pd.ExcelFile(path)
        if sheet is None:
            # pick the first sheet that is not a "_TOP_ORDERS" summary sheet
            candidates = [s for s in xls.sheet_names if "TOP" not in s.upper()]
            sheet = candidates[0] if candidates else xls.sheet_names[0]
        return xls.parse(sheet)
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    # feature tables exported from BIOM sometimes carry a leading comment row
    try:
        return pd.read_csv(path, sep=sep)
    except Exception:
        return pd.read_csv(path, sep=sep, skiprows=1)


def detect_sample_columns(df: pd.DataFrame, metadata_ids: Iterable[str] | None = None) -> list[str]:
    """
    Dynamically identify the sample (abundance) columns. No assumption that they
    start with 'sample'. A column is a sample column when it is:
      - not a taxonomy rank / feature-id column, AND
      - numeric (or coercible to numeric), AND
      - (optionally) present in the provided metadata sample-id list.
    """
    metadata_set = {str(s).strip() for s in metadata_ids} if metadata_ids else None
    sample_cols: list[str] = []
    for col in df.columns:
        if str(col) in NON_SAMPLE_COLUMNS:
            continue
        coerced = pd.to_numeric(df[col], errors="coerce")
        if coerced.notna().any():
            if metadata_set is None or str(col).strip() in metadata_set:
                sample_cols.append(col)
    return sample_cols


def load_metadata_sample_ids(metadata_path: str | Path | None) -> list[str] | None:
    if not metadata_path:
        return None
    metadata_path = Path(metadata_path)
    if not metadata_path.exists():
        return None
    md = pd.read_csv(metadata_path, sep="\t", dtype=str)
    id_col = md.columns[0]
    return [str(x).strip() for x in md[id_col].dropna().tolist()]


# --------------------------------------------------------------------------- #
# 2. Collapse to taxonomic rank
# --------------------------------------------------------------------------- #
def collapse_to_rank(df: pd.DataFrame, sample_cols: list[str], rank: str) -> pd.DataFrame:
    """
    Sum abundances by the chosen taxonomic rank -> DataFrame (taxon x sample).
    Unassigned / empty taxa are dropped so they cannot form meaningless edges.
    """
    if rank not in df.columns:
        raise ValueError(
            f"Rank '{rank}' not found. Available taxonomy columns: "
            f"{[c for c in df.columns if c in TAX_RANKS]}"
        )
    work = df[[rank] + sample_cols].copy()
    work[rank] = work[rank].fillna("Unassigned").astype(str).str.strip()
    for c in sample_cols:
        work[c] = pd.to_numeric(work[c], errors="coerce").fillna(0.0)
    grouped = work.groupby(rank)[sample_cols].sum()
    drop = {"Unassigned", "", "unassigned", "NA", "nan", "None", "__"}
    grouped = grouped.loc[~grouped.index.isin(drop)]
    grouped = grouped.loc[~grouped.index.str.fullmatch(r"[a-z]__?", case=False, na=False)]
    return grouped


# --------------------------------------------------------------------------- #
# 3. Prevalence / abundance filter
# --------------------------------------------------------------------------- #
def filter_taxa(table: pd.DataFrame, min_prevalence: float, min_abundance: float):
    """Keep taxa present in >= min_prevalence of samples and above min_abundance."""
    n_samples = table.shape[1]
    prevalence = (table > 0).sum(axis=1) / max(n_samples, 1)
    total = table.sum(axis=1)
    keep = (prevalence >= min_prevalence) & (total >= min_abundance)
    filtered = table.loc[keep].copy()
    prevalence = prevalence.loc[keep]
    return filtered, prevalence


# --------------------------------------------------------------------------- #
# 4. CLR transform (compositionally aware)
# --------------------------------------------------------------------------- #
def clr_transform(table: pd.DataFrame, pseudocount: float = 1.0) -> pd.DataFrame:
    """
    Centred-log-ratio transform. Input is taxon x sample; CLR is computed per
    sample (column) across taxa: clr = log(x) - mean(log(x)).
    Returns taxon x sample DataFrame of CLR values.
    """
    mat = table.to_numpy(dtype=float) + pseudocount
    log_mat = np.log(mat)
    gm = log_mat.mean(axis=0, keepdims=True)     # per-sample geometric mean (in log space)
    clr = log_mat - gm
    return pd.DataFrame(clr, index=table.index, columns=table.columns)


# --------------------------------------------------------------------------- #
# 5-6. Correlations, p-values, BH-FDR, edge selection
# --------------------------------------------------------------------------- #
def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR-adjusted q-values (1995)."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(1, n + 1))
    # enforce monotonicity from the largest p down
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty_like(q)
    out[order] = q
    return out


def compute_pairwise_correlations(clr_table: pd.DataFrame, raw_table: pd.DataFrame,
                                  method: str = "spearman") -> pd.DataFrame:
    """
    Compute correlations for every taxon pair across samples on the CLR matrix.
    Returns a long-format DataFrame: taxon_a, taxon_b, rho, p_value, n_samples,
    co_prevalence (fraction of samples where both taxa are present).
    """
    taxa = list(clr_table.index)
    n_taxa = len(taxa)
    n_samples = clr_table.shape[1]
    # orient samples x taxa for correlation across samples
    x = clr_table.to_numpy(dtype=float).T          # samples x taxa
    presence = (raw_table.to_numpy() > 0)          # taxa x samples

    if method == "pearson":
        # vectorised Pearson
        xc = x - x.mean(axis=0, keepdims=True)
        std = xc.std(axis=0, ddof=1, keepdims=True)
        std[std == 0] = np.nan
        corr = (xc.T @ xc) / (n_samples - 1) / (std.T @ std)
    else:
        # Spearman = Pearson on ranks; scipy returns full matrix + p-values,
        # but we recompute p from the t-approximation for consistency below.
        ranks = np.apply_along_axis(stats.rankdata, 0, x)
        rc = ranks - ranks.mean(axis=0, keepdims=True)
        std = rc.std(axis=0, ddof=1, keepdims=True)
        std[std == 0] = np.nan
        corr = (rc.T @ rc) / (n_samples - 1) / (std.T @ std)

    rows = []
    for i in range(n_taxa):
        for j in range(i + 1, n_taxa):
            rho = corr[i, j]
            if not np.isfinite(rho):
                continue
            rho = float(np.clip(rho, -0.999999, 0.999999))
            # two-sided p-value via t-approximation
            if n_samples > 2:
                t = rho * np.sqrt((n_samples - 2) / (1 - rho ** 2))
                p = 2 * stats.t.sf(abs(t), df=n_samples - 2)
            else:
                p = np.nan
            co_prev = float(np.mean(presence[i] & presence[j]))
            rows.append({
                "taxon_a": taxa[i],
                "taxon_b": taxa[j],
                "rho": rho,
                "p_value": float(p),
                "n_samples": int(n_samples),
                "co_prevalence": co_prev,
            })
    edges = pd.DataFrame(rows)
    if not edges.empty:
        edges["q_value"] = benjamini_hochberg(edges["p_value"].to_numpy())
        edges["direction"] = np.where(edges["rho"] >= 0, "co-occurrence (+)", "exclusion (-)")
    return edges


def select_significant(edges: pd.DataFrame, fdr_alpha: float, min_abs_corr: float) -> pd.DataFrame:
    if edges.empty:
        return edges
    keep = (edges["q_value"] <= fdr_alpha) & (edges["rho"].abs() >= min_abs_corr)
    return edges.loc[keep].sort_values("q_value").reset_index(drop=True)


# --------------------------------------------------------------------------- #
# 7. Evidence base cross-reference
# --------------------------------------------------------------------------- #
def load_reference_db(path: str | Path | None) -> dict:
    if not path:
        return {"interactions": [], "method_references": []}
    path = Path(path)
    if not path.exists():
        return {"interactions": [], "method_references": []}
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def _name_matches(name: str, token: str) -> bool:
    return token.lower() in str(name).lower()


def annotate_with_evidence(edges: pd.DataFrame, reference_db: dict) -> pd.DataFrame:
    """
    For each edge, look for a curated interaction whose two partners each match
    one of the edge's taxa (in either orientation). Adds columns:
      evidence_status  : SUPPORTED / NOVEL
      literature_verdict: CONSISTENT / DISCORDANT / CONTEXT / n.a.
      interaction_type, mechanism, references (joined), reference_ids
    """
    if edges.empty:
        return edges
    interactions = reference_db.get("interactions", [])
    out = edges.copy()
    status, verdict, itypes, mechs, refs, ids = [], [], [], [], [], []

    for _, row in out.iterrows():
        a, b = str(row["taxon_a"]), str(row["taxon_b"])
        observed_sign = "positive" if row["rho"] >= 0 else "negative"
        hit = None
        for entry in interactions:
            ta, tb = entry["taxon_a"], entry["taxon_b"]
            if (_name_matches(a, ta) and _name_matches(b, tb)) or \
               (_name_matches(a, tb) and _name_matches(b, ta)):
                hit = entry
                break
        if hit is None:
            status.append("NOVEL")
            verdict.append("n.a.")
            itypes.append("")
            mechs.append("")
            refs.append("")
            ids.append("")
        else:
            status.append("SUPPORTED")
            exp = hit.get("expected_sign", "context")
            if exp == "context":
                verdict.append("CONTEXT")
            elif exp == observed_sign:
                verdict.append("CONSISTENT")
            else:
                verdict.append("DISCORDANT")
            itypes.append(hit.get("interaction_type", ""))
            mechs.append(hit.get("mechanism", ""))
            refs.append(" | ".join(r.get("citation", "") for r in hit.get("references", [])))
            ids.append(hit.get("id", ""))

    out["evidence_status"] = status
    out["literature_verdict"] = verdict
    out["interaction_type"] = itypes
    out["mechanism"] = mechs
    out["references"] = refs
    out["reference_ids"] = ids
    return out


# --------------------------------------------------------------------------- #
# Graph metrics
# --------------------------------------------------------------------------- #
def compute_node_stats(sig_edges: pd.DataFrame) -> pd.DataFrame:
    if sig_edges.empty:
        return pd.DataFrame(columns=["taxon", "degree", "pos_degree", "neg_degree", "mean_abs_rho"])
    records: dict[str, dict[str, Any]] = {}
    for _, e in sig_edges.iterrows():
        for t in (e["taxon_a"], e["taxon_b"]):
            r = records.setdefault(t, {"degree": 0, "pos_degree": 0, "neg_degree": 0, "abs": []})
            r["degree"] += 1
            if e["rho"] >= 0:
                r["pos_degree"] += 1
            else:
                r["neg_degree"] += 1
            r["abs"].append(abs(e["rho"]))
    rows = [{
        "taxon": t,
        "degree": r["degree"],
        "pos_degree": r["pos_degree"],
        "neg_degree": r["neg_degree"],
        "mean_abs_rho": float(np.mean(r["abs"])),
    } for t, r in records.items()]
    return pd.DataFrame(rows).sort_values("degree", ascending=False).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Orchestrator
# --------------------------------------------------------------------------- #
def run_network_analysis(
    df: pd.DataFrame | None = None,
    input_path: str | Path | None = None,
    sheet: str | None = None,
    metadata_path: str | Path | None = None,
    params: NetworkParams | None = None,
    **overrides,
) -> NetworkResult:
    """
    End-to-end analysis. Provide either a DataFrame (df) or an input_path.
    Any NetworkParams field can be overridden via keyword args.
    """
    params = params or (NetworkParams.from_config() if cfg else NetworkParams())
    for k, v in overrides.items():
        if hasattr(params, k) and v is not None:
            setattr(params, k, v)

    if df is None:
        if input_path is None:
            raise ValueError("Provide either df or input_path.")
        df = load_feature_table(input_path, sheet=sheet)

    metadata_ids = load_metadata_sample_ids(metadata_path)
    sample_cols = detect_sample_columns(df, metadata_ids)
    if len(sample_cols) < 3:
        raise ValueError(
            f"Only {len(sample_cols)} sample column(s) detected; need >= 3 to "
            "estimate correlations. Check the input table / metadata."
        )

    taxa_all = collapse_to_rank(df, sample_cols, params.rank)
    taxa_table, prevalence = filter_taxa(taxa_all, params.min_prevalence, params.min_abundance)
    if taxa_table.shape[0] < 2:
        raise ValueError(
            f"Only {taxa_table.shape[0]} taxon passed the prevalence filter "
            f"(min_prevalence={params.min_prevalence}). Lower the threshold."
        )

    clr = clr_transform(taxa_table, params.clr_pseudocount)
    edges = compute_pairwise_correlations(clr, taxa_table, params.correlation_method)

    ref_db = load_reference_db(params.reference_db)
    sig = select_significant(edges, params.fdr_alpha, params.min_abs_correlation)
    sig = annotate_with_evidence(sig, ref_db)
    node_stats = compute_node_stats(sig)

    summary = {
        "rank": params.rank,
        "n_samples": len(sample_cols),
        "sample_columns": sample_cols,
        "n_taxa_total": int(taxa_all.shape[0]),
        "n_taxa_analyzed": int(taxa_table.shape[0]),
        "n_pairs_tested": int(edges.shape[0]),
        "n_significant_edges": int(sig.shape[0]),
        "n_supported_by_literature": int((sig["evidence_status"] == "SUPPORTED").sum()) if not sig.empty else 0,
        "n_consistent": int((sig["literature_verdict"] == "CONSISTENT").sum()) if not sig.empty else 0,
        "n_discordant": int((sig["literature_verdict"] == "DISCORDANT").sum()) if not sig.empty else 0,
        "method": params.correlation_method,
        "fdr_alpha": params.fdr_alpha,
        "min_abs_correlation": params.min_abs_correlation,
        "method_references": ref_db.get("method_references", []),
    }
    return NetworkResult(params, taxa_table, edges, sig, node_stats, summary)


def write_outputs(result: NetworkResult, out_dir: str | Path) -> dict[str, Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    p = out_dir / "interaction_edges_all.csv"
    result.edges.to_csv(p, index=False); paths["all_edges"] = p
    p = out_dir / "interaction_edges_significant.csv"
    result.significant_edges.to_csv(p, index=False); paths["significant_edges"] = p
    p = out_dir / "interaction_nodes.csv"
    result.node_stats.to_csv(p, index=False); paths["nodes"] = p
    p = out_dir / "interaction_summary.json"
    with p.open("w", encoding="utf-8") as fh:
        json.dump(result.summary, fh, indent=2)
    paths["summary"] = p
    return paths


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _default_input() -> Path | None:
    if cfg is None:
        return None
    candidate = cfg.EXPORTED_DIR / "asv_tables_combined.xlsx"
    return candidate if candidate.exists() else None


def main() -> None:
    ap = argparse.ArgumentParser(description="Microbial interaction network with evidence base.")
    ap.add_argument("--input", type=str, default=None, help="ASV table (.xlsx/.tsv/.csv).")
    ap.add_argument("--sheet", type=str, default=None, help="Excel sheet name (marker).")
    ap.add_argument("--metadata", type=str, default=None, help="metadata.tsv for sample-id matching.")
    ap.add_argument("--rank", type=str, default=None, help="Taxonomic rank (default from config).")
    ap.add_argument("--min-prevalence", type=float, default=None)
    ap.add_argument("--min-abs-correlation", type=float, default=None)
    ap.add_argument("--fdr-alpha", type=float, default=None)
    ap.add_argument("--out", type=str, default=None, help="Output directory.")
    args = ap.parse_args()

    input_path = args.input or _default_input()
    if input_path is None:
        print("No input table given and no default asv_tables_combined.xlsx found.")
        print("Run with:  python network_analysis.py --input path/to/table.xlsx")
        sys.exit(1)

    metadata_path = args.metadata
    if metadata_path is None and cfg is not None:
        md = cfg.DATA_DIR / "metadata.tsv"
        metadata_path = str(md) if md.exists() else None

    result = run_network_analysis(
        input_path=input_path, sheet=args.sheet, metadata_path=metadata_path,
        rank=args.rank, min_prevalence=args.min_prevalence,
        min_abs_correlation=args.min_abs_correlation, fdr_alpha=args.fdr_alpha,
    )

    out_dir = args.out or (str(cfg.RESULTS_DIR / "network") if cfg else "network_out")
    paths = write_outputs(result, out_dir)

    s = result.summary
    print("\n=== Microbial Interaction Network ===")
    print(f"Rank                 : {s['rank']}")
    print(f"Samples              : {s['n_samples']}")
    print(f"Taxa analysed        : {s['n_taxa_analyzed']} / {s['n_taxa_total']}")
    print(f"Pairs tested         : {s['n_pairs_tested']}")
    print(f"Significant edges    : {s['n_significant_edges']} "
          f"(FDR<= {s['fdr_alpha']}, |rho|>= {s['min_abs_correlation']})")
    print(f"Literature-supported : {s['n_supported_by_literature']} "
          f"(consistent={s['n_consistent']}, discordant={s['n_discordant']})")
    print("\nOutputs:")
    for k, v in paths.items():
        print(f"   {k:18s}: {v}")


if __name__ == "__main__":
    main()