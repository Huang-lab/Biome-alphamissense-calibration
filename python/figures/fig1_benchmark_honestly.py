"""Figure 1 — AM benchmarking honestly: easy ClinVar overstates accuracy.

Four panels per the rebuttal plan (Action 2 / R3):

  1A  AUROC stratified by ClinVar review-status star (1 / 2 / 3+).
      *Needs* a ClinVar VCF with the CLNREVSTAT field. If unavailable,
      panel emits a placeholder rectangle with the message and prints
      a warning — the figure file still renders so the rest of the
      pipeline doesn't break.

  1B  Per-gene FPR at 0.864 (global) vs at the Chen gene-specific
      threshold. Source = the BP4_* rows of the Chen subset; FPR =
      fraction of BP4 variants whose AM score (``vep_score``) crosses
      the threshold. The 21% headline number is the across-gene
      median at 0.864.

  1C  AM-on-VUS: what fraction of ClinVar VUSs (2+ stars, no conflicts)
      AM calls deleterious under 0.864 vs under gene-specific. Same
      ClinVar dependency as 1A.

  1D  Per-gene threshold map: 28 genes annotated with their Chen
      single-gene threshold, the domain-aggregate threshold they fall
      under (when single-gene is unavailable), or "uncovered". Source =
      ``refs/chen_summary_by_gene.tsv`` produced by step 00.

Inputs (paths):
    chen_subset           refs/chen_calibration.target_genes.tsv.gz
    chen_summary_by_gene  refs/chen_summary_by_gene.tsv
    clinvar_vcf           refs/clinvar.vcf.gz  (optional)
    out_dir               where to write fig1.{png,pdf}

Usage:
    python -m python.figures.fig1_benchmark_honestly \\
        --chen-subset refs/chen_calibration.target_genes.tsv.gz \\
        --chen-summary refs/chen_summary_by_gene.tsv \\
        --out-dir results/figures
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))  # so `import util` works
from util import LOG  # noqa: E402

from . import common


# ---------------------------------------------------------------------------
# Chen subset reader (BP4 rows for FPR computation)
# ---------------------------------------------------------------------------
def load_chen_bp4(chen_subset: str) -> "pd.DataFrame":
    """Return rows whose ``evidence`` starts with ``BP4`` (i.e. the
    anti-pathogenic / true-negative ground truth). Columns kept:
    gene_symbol, vep_score (AM score)."""
    import pandas as pd
    df = pd.read_csv(chen_subset, sep="\t", dtype=str, comment="#", header=None,
                     names=["chrom","pos","ref","alt","gene_symbol","evidence",
                            "points","calibration_approach","domain","vep_score"])
    df = df[df["evidence"].str.startswith("BP4", na=False)]
    df["vep_score"] = pd.to_numeric(df["vep_score"], errors="coerce")
    df = df.dropna(subset=["vep_score"])
    return df.reset_index(drop=True)


def per_gene_fpr(bp4: "pd.DataFrame", threshold: float) -> "pd.Series":
    """FPR per gene_symbol at ``threshold``: fraction of BP4 rows with
    ``vep_score >= threshold``. Genes with no BP4 rows are NaN."""
    if bp4.empty:
        import pandas as pd
        return pd.Series(dtype=float)
    return bp4.groupby("gene_symbol").apply(
        lambda g: (g["vep_score"] >= threshold).mean()
    )


# ---------------------------------------------------------------------------
# Per-gene calibration map (Chen single-gene threshold by gene)
# ---------------------------------------------------------------------------
def load_threshold_map(chen_summary_path: str, chen_subset_path: str
                       ) -> "pd.DataFrame":
    """Combine chen_summary_by_gene with the implied PP3_Moderate threshold per
    gene (the lowest vep_score in the subset whose evidence is PP3_Moderate
    or stronger). Output columns:
      gene, syndrome (filled by caller), threshold_PP3_Moderate, approach.
    """
    import pandas as pd
    summary = pd.read_csv(chen_summary_path, sep="\t")
    df = pd.read_csv(
        chen_subset_path, sep="\t", dtype=str, comment="#", header=None,
        names=["chrom","pos","ref","alt","gene_symbol","evidence","points",
               "calibration_approach","domain","vep_score"],
    )
    df["vep_score"] = pd.to_numeric(df["vep_score"], errors="coerce")
    positive = ["pp3_moderate", "pp3_moderate+", "pp3_strong", "pp3_strong+",
                "pp3_verystrong", "pp3_very_strong"]
    is_pos = df["evidence"].str.lower().str.replace(" ", "").str.replace("-", "_").isin(positive)
    pos_rows = df.loc[is_pos & df["vep_score"].notna()]
    thr = pos_rows.groupby("gene_symbol")["vep_score"].min().rename("threshold_PP3_Moderate")
    out = summary.merge(thr, left_on="gene", right_index=True, how="left")
    return out


# ---------------------------------------------------------------------------
# Panel functions
# ---------------------------------------------------------------------------
def _placeholder_panel(ax, msg: str, color: str = "#bbbbbb") -> None:
    ax.set_facecolor("#f5f5f5")
    ax.text(0.5, 0.5, msg, ha="center", va="center", fontsize=8,
            transform=ax.transAxes, color="#333333",
            wrap=True)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_edgecolor(color)


def panel_1A(ax, clinvar_path: Optional[str]) -> None:
    """AUROC by ClinVar review-status star. Falls back to a placeholder
    when no ClinVar VCF is provided (Q3 in the plan)."""
    if not clinvar_path or not os.path.isfile(clinvar_path):
        _placeholder_panel(
            ax,
            "Fig 1A pending\nClinVar VCF (CLNREVSTAT) required.\n"
            "Set references.clinvar_vcf_local in config.yaml\n"
            "to enable stratified AUROC by review-status star.",
        )
        ax.set_title("A  AUROC by ClinVar review-status star")
        return
    # Real implementation lands once ClinVar is wired into prepare_refs.
    _placeholder_panel(ax, "Fig 1A — ClinVar parsing TBD")
    ax.set_title("A  AUROC by ClinVar review-status star")


def panel_1B(ax, chen_subset: str, global_threshold: float = 0.864) -> None:
    """Per-gene FPR at 0.864 (global) vs at the gene-specific PP3_Moderate
    threshold. Bars are sorted by per-gene 0.864 FPR descending."""
    import numpy as np
    import pandas as pd

    bp4 = load_chen_bp4(chen_subset)
    if bp4.empty:
        _placeholder_panel(ax, "Fig 1B — no BP4 rows in Chen subset")
        ax.set_title("B  Per-gene FPR  (0.864 vs gene-specific)")
        return

    fpr_global = per_gene_fpr(bp4, global_threshold)
    # gene-specific FPR: per gene, count BP4 rows that cross THAT gene's
    # min PP3_Moderate threshold. Genes with no PP3_Moderate rows get NaN.
    thr_by_gene = (
        bp4  # dummy frame; we need PP3 rows for thresholds
    )
    # We don't have direct PP3 rows here because load_chen_bp4 already
    # filtered to BP4 — load the full subset just for the threshold lookup.
    full = pd.read_csv(
        chen_subset, sep="\t", dtype=str, comment="#", header=None,
        names=["chrom","pos","ref","alt","gene_symbol","evidence","points",
               "calibration_approach","domain","vep_score"],
    )
    full["vep_score"] = pd.to_numeric(full["vep_score"], errors="coerce")
    pp3 = full[full["evidence"].str.lower().isin(
        ["pp3_moderate","pp3_moderate+","pp3_strong","pp3_strong+",
         "pp3_verystrong","pp3_very_strong"])]
    gene_thr = pp3.groupby("gene_symbol")["vep_score"].min()

    fpr_gene = pd.Series(dtype=float)
    for g, thr in gene_thr.items():
        rows = bp4[bp4["gene_symbol"] == g]
        if rows.empty or pd.isna(thr):
            continue
        fpr_gene[g] = (rows["vep_score"] >= thr).mean()

    order = fpr_global.sort_values(ascending=False).index.tolist()
    x = np.arange(len(order))
    h_global = [fpr_global.get(g, np.nan) for g in order]
    h_gene   = [fpr_gene.get(g, np.nan)   for g in order]
    ax.bar(x - 0.2, h_global, width=0.4, color=common.CATEGORY_COLORS["AM_global_0864"],
           label="global 0.864")
    ax.bar(x + 0.2, h_gene,   width=0.4, color=common.CATEGORY_COLORS["AM_primary"],
           label="gene-specific (Chen)")
    ax.axhline(0.007, color="grey", lw=0.5, ls="--")
    ax.text(len(order) - 0.5, 0.012, "advertised 0.007",
            fontsize=6, ha="right", color="grey")
    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=90, fontsize=6)
    ax.set_ylabel("False-positive rate (BP4 rows)")
    ax.legend(loc="upper right")
    ax.set_title("B  Per-gene FPR  (0.864 vs gene-specific)")


def panel_1C(ax, clinvar_path: Optional[str]) -> None:
    """AM on ClinVar VUSs (2+ stars, no conflicts). Placeholder until ClinVar
    is wired in (same dependency as 1A)."""
    if not clinvar_path or not os.path.isfile(clinvar_path):
        _placeholder_panel(
            ax,
            "Fig 1C pending\nClinVar VCF needed for the VUS subset\n"
            "(2+ stars, no conflicts). Pairs with 1B + the\n"
            "AM-as-prioritizer-not-classifier framing.",
        )
        ax.set_title("C  AM on ClinVar VUSs  (2+ stars, no conflicts)")
        return
    _placeholder_panel(ax, "Fig 1C — ClinVar VUS subset TBD")
    ax.set_title("C  AM on ClinVar VUSs  (2+ stars, no conflicts)")


def panel_1D(ax, chen_summary_path: str, chen_subset_path: str) -> None:
    """Per-gene calibration map: 28 genes, color-coded by chen_approach,
    annotated with the gene-specific threshold when available."""
    import numpy as np
    import pandas as pd
    df = load_threshold_map(chen_summary_path, chen_subset_path)
    # Order genes by the panel order, fill in missing as 'not_in_chen'.
    df = df.set_index("gene").reindex(common.GENE_ORDER).reset_index()
    df["approach"] = df["calibration_approach_majority"].fillna("none")
    color_for = {
        "single_gene":      common.CATEGORY_COLORS["AM_primary"],
        "domain_aggregate": common.CATEGORY_COLORS["AM_only_non_PLP"],
        "none":             "#bbbbbb",
    }
    colors = [color_for.get(a, "#bbbbbb") for a in df["approach"]]
    y = np.arange(len(df))[::-1]
    ax.barh(y, [1.0] * len(df), color=colors, edgecolor="white", lw=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(df["gene"], fontsize=6)
    ax.set_xticks([])
    ax.set_xlim(0, 1.5)
    # Threshold label
    for yi, thr in zip(y, df["threshold_PP3_Moderate"]):
        if pd.notna(thr):
            ax.text(1.02, yi, f"thr={thr:.3f}", va="center", fontsize=6)
        else:
            ax.text(1.02, yi, "uncovered", va="center", fontsize=6,
                    color="#888888", style="italic")
    handles = [
        # Legend swatches
    ]
    from matplotlib.patches import Patch
    handles = [
        Patch(facecolor=color_for["single_gene"],      label="single_gene"),
        Patch(facecolor=color_for["domain_aggregate"], label="domain_aggregate"),
        Patch(facecolor=color_for["none"],             label="uncovered"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=6, ncol=3,
              bbox_to_anchor=(1.0, -0.08))
    ax.set_title("D  Per-gene calibration approach + threshold")


# ---------------------------------------------------------------------------
# make()
# ---------------------------------------------------------------------------
def make(out_dir: str, *,
         chen_subset: str,
         chen_summary: str,
         clinvar_vcf: Optional[str] = None) -> None:
    import matplotlib.pyplot as plt
    common.apply_rcparams()
    fig, axes = plt.subplots(2, 2, figsize=(common.FIG_WIDTH_DOUBLE[0],
                                            common.FIG_WIDTH_DOUBLE[0] * 1.0))
    panel_1A(axes[0, 0], clinvar_vcf)
    panel_1B(axes[0, 1], chen_subset)
    panel_1C(axes[1, 0], clinvar_vcf)
    panel_1D(axes[1, 1], chen_summary, chen_subset)
    fig.tight_layout()
    common.save_both(fig, Path(out_dir), "fig1")
    plt.close(fig)
    LOG.info("fig1: wrote fig1.png + fig1.pdf to %s", out_dir)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--chen-subset", required=True)
    ap.add_argument("--chen-summary", required=True)
    ap.add_argument("--clinvar-vcf", default=None)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args(argv)
    make(args.out_dir, chen_subset=args.chen_subset, chen_summary=args.chen_summary,
         clinvar_vcf=args.clinvar_vcf)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
