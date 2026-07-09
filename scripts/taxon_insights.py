"""
Taxon Insights — explain WHY taxa are high/low, what it means, what to do, with evidence.

General-purpose (not peanut-specific): given any ASV/feature table with taxonomy,
this computes each taxon's mean relative abundance + prevalence, flags which taxa
are notably HIGH or LOW, and annotates them with a literature-backed knowledge
base (scripts/taxon_insights.json): ecological role, drivers (why high / why low),
agronomic/ecological implication, candidate interventions, and citations.

It reuses the dynamic sample-detection and rank-collapsing from network_analysis,
so nothing about sample naming or markers is hard-coded.

Usage
-----
CLI:
    python taxon_insights.py --input data/exported/asv_tables_combined.xlsx --sheet 16S

Importable:
    from taxon_insights import generate_taxon_insights
    result = generate_taxon_insights(df=my_table)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
try:
    import config as cfg
except Exception:
    cfg = None

# Reuse the tested helpers from the network module (dynamic sample detection etc.)
from network_analysis import (  # noqa: E402
    load_feature_table, detect_sample_columns, load_metadata_sample_ids, collapse_to_rank,
)

RANKS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


def _default_kb_path() -> Path:
    if cfg is not None:
        return cfg.REFERENCE_DIR / "taxon_insights.json"
    return Path(__file__).resolve().parent / "taxon_insights.json"


def load_knowledge_base(path: str | Path | None = None) -> dict:
    path = Path(path) if path else _default_kb_path()
    if not path.exists():
        return {"taxa": [], "factors": []}
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def mean_relative_abundance(df: pd.DataFrame, sample_cols: list[str], rank: str) -> pd.DataFrame:
    """Return per-taxon mean relative abundance (%) and prevalence at a rank."""
    table = collapse_to_rank(df, sample_cols, rank)          # taxon x sample (counts)
    col_sums = table.sum(axis=0).replace(0, np.nan)
    rel = table.divide(col_sums, axis=1).fillna(0.0)         # per-sample relative
    out = pd.DataFrame({
        "taxon": table.index,
        "mean_rel_abundance_pct": (rel.mean(axis=1) * 100).values,
        "prevalence": ((table > 0).sum(axis=1) / max(table.shape[1], 1)).values,
    })
    return out.sort_values("mean_rel_abundance_pct", ascending=False).reset_index(drop=True)


def _status(pct: float, series: pd.Series) -> str:
    """Label a taxon High/Medium/Low relative to the dataset distribution."""
    if pct <= 0:
        return "Absent"
    hi = series[series > 0].quantile(0.75)
    lo = series[series > 0].quantile(0.25)
    if pct >= hi:
        return "High"
    if pct <= lo:
        return "Low"
    return "Medium"


def _match(observed_name: str, kb_taxon: str) -> bool:
    o, k = str(observed_name).lower(), kb_taxon.lower()
    return o == k or k in o


def generate_taxon_insights(
    df: pd.DataFrame | None = None,
    input_path: str | Path | None = None,
    sheet: str | None = None,
    metadata_path: str | Path | None = None,
    kb_path: str | Path | None = None,
    top_n: int = 15,
) -> dict[str, Any]:
    """
    Returns a dict with:
      abundance   : full per-genus mean relative abundance table
      insights    : KB-annotated rows for taxa found (cause/meaning/solution/evidence)
      top_taxa    : the most abundant taxa overall (with KB note if available)
      factors     : general driver explanations from the KB
      summary     : counts
    """
    if df is None:
        if input_path is None:
            raise ValueError("Provide df or input_path.")
        df = load_feature_table(input_path, sheet=sheet)

    kb = load_knowledge_base(kb_path)
    metadata_ids = load_metadata_sample_ids(metadata_path)
    sample_cols = detect_sample_columns(df, metadata_ids)
    if not sample_cols:
        raise ValueError("No numeric sample columns detected in this table.")

    # Precompute abundance at every rank present, so both genus- and phylum-level
    # KB entries can be matched.
    available_ranks = [r for r in RANKS if r in df.columns]
    rank_tables: dict[str, pd.DataFrame] = {}
    for r in available_ranks:
        try:
            rank_tables[r] = mean_relative_abundance(df, sample_cols, r)
        except Exception:
            continue

    # Genus-level (or finest available) table for the headline "top taxa" view.
    headline_rank = "Genus" if "Genus" in rank_tables else (available_ranks[-1] if available_ranks else None)
    abundance = rank_tables.get(headline_rank, pd.DataFrame())

    # Annotate KB taxa found in the data.
    insight_rows = []
    for entry in kb.get("taxa", []):
        rank = entry.get("rank", "Genus")
        tbl = rank_tables.get(rank)
        if tbl is None or tbl.empty:
            continue
        hits = tbl[tbl["taxon"].apply(lambda x: _match(x, entry["taxon"]))]
        if hits.empty:
            status, pct, prev, observed = "Absent", 0.0, 0.0, ""
        else:
            row = hits.iloc[0]
            observed = row["taxon"]
            pct = float(row["mean_rel_abundance_pct"])
            prev = float(row["prevalence"])
            status = _status(pct, tbl["mean_rel_abundance_pct"])
        refs = " | ".join(r.get("citation", "") for r in entry.get("references", []))
        ref_urls = " | ".join(r.get("url", "") for r in entry.get("references", []))
        insight_rows.append({
            "taxon": entry["taxon"], "rank": rank, "group": entry.get("group", ""),
            "observed_as": observed, "status": status,
            "mean_rel_abundance_pct": round(pct, 3), "prevalence": round(prev, 3),
            "role": entry.get("role", ""),
            "why_high": entry.get("high_when", ""),
            "why_low": entry.get("low_when", ""),
            "implication": entry.get("implication", ""),
            "interventions": entry.get("interventions", ""),
            "references": refs, "reference_urls": ref_urls,
        })
    insights = pd.DataFrame(insight_rows)
    if not insights.empty:
        # Present found taxa first, ordered by abundance.
        insights["_found"] = insights["status"] != "Absent"
        insights = insights.sort_values(["_found", "mean_rel_abundance_pct"],
                                        ascending=[False, False]).drop(columns="_found").reset_index(drop=True)

    # Top taxa overall, with KB note where available.
    kb_names = {e["taxon"].lower(): e for e in kb.get("taxa", [])}
    top = abundance.head(top_n).copy() if not abundance.empty else pd.DataFrame()
    if not top.empty:
        def kb_note(name):
            for kname, e in kb_names.items():
                if _match(name, kname):
                    return e.get("group", "")
            return ""
        top["kb_group"] = top["taxon"].apply(kb_note)

    summary = {
        "n_samples": len(sample_cols),
        "headline_rank": headline_rank,
        "n_taxa_at_headline_rank": int(abundance.shape[0]) if not abundance.empty else 0,
        "n_kb_taxa_found": int((insights["status"] != "Absent").sum()) if not insights.empty else 0,
        "n_kb_taxa_total": len(kb.get("taxa", [])),
    }
    return {
        "abundance": abundance, "insights": insights, "top_taxa": top,
        "factors": kb.get("factors", []), "disclaimer": kb.get("disclaimer", ""),
        "summary": summary,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Taxon insights: why high/low + cause + solution + evidence.")
    ap.add_argument("--input", type=str, default=None)
    ap.add_argument("--sheet", type=str, default=None)
    ap.add_argument("--metadata", type=str, default=None)
    ap.add_argument("--top-n", type=int, default=15)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    input_path = args.input
    if input_path is None and cfg is not None:
        cand = cfg.EXPORTED_DIR / "asv_tables_combined.xlsx"
        input_path = str(cand) if cand.exists() else None
    if not input_path:
        print("Provide --input path to an ASV table (.xlsx/.tsv/.csv).")
        sys.exit(1)

    res = generate_taxon_insights(input_path=input_path, sheet=args.sheet,
                                  metadata_path=args.metadata, top_n=args.top_n)
    s = res["summary"]
    print("\n=== Taxon Insights ===")
    print(f"Samples: {s['n_samples']} | rank: {s['headline_rank']} | "
          f"KB taxa found: {s['n_kb_taxa_found']}/{s['n_kb_taxa_total']}")
    found = res["insights"][res["insights"]["status"] != "Absent"] if not res["insights"].empty else pd.DataFrame()
    if not found.empty:
        print("\nNotable taxa (with explanation):")
        for _, r in found.iterrows():
            print(f"\n• {r['taxon']} ({r['group']}) — {r['status']} "
                  f"[{r['mean_rel_abundance_pct']}% mean rel. abundance]")
            print(f"    role: {r['role']}")
            print(f"    why high: {r['why_high']}")
            print(f"    why low : {r['why_low']}")
            print(f"    meaning : {r['implication']}")
            print(f"    action  : {r['interventions']}")
            print(f"    evidence: {r['references']}")
    if args.out:
        out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
        res["abundance"].to_csv(out / "taxon_abundance.csv", index=False)
        res["insights"].to_csv(out / "taxon_insights.csv", index=False)
        print(f"\nSaved CSVs to {out}")


if __name__ == "__main__":
    main()
