"""
Fig 5 — Carrier frequency by ancestry.

5A: Grouped bar chart — x = ancestry, y = carrier freq (%), 3 bars per group
    No error bars; N shown in x-axis labels; % at bar top.
5B: Gene × ancestry heatmap (seaborn style, one per variant class)
    color scale 0–2% white→red, annotate ≥1%.

Usage:
    python -m python.figures.fig5_ancestry --cohort cohortI
"""

import argparse
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests

from python.figures.common import (
    build_variant_table, load_metadata_groups,
    COLORS, VC_ACMG, VC_CLINVAR, VC_AM, VC_AM_ONLY,
    ALL_GENES_ORDERED, GENE_TO_SYNDROME,
    FIGW_DOUBLE,
    save_fig, save_table, cohort_label,
)

HEATMAP_VMAX = 3.0   # colorbar max 0–3%; annotate only cells >3%

# Ordered variant tracks shown in all fig5 panels.
VC_TRACKS = [VC_ACMG, VC_CLINVAR, VC_AM, VC_AM_ONLY]


def _vc_carrier_filter(vt, vc):
    """Per-track carrier subset of the variant table."""
    if vc == VC_ACMG:
        return vt[vt["in_acmg"] == True]
    if vc == VC_CLINVAR:
        return vt[vt["in_clinvar"] == True]
    if vc == VC_AM:
        return vt[vt["in_am"] == True]
    return vt[vt["variant_class"] == VC_AM_ONLY]


def make(cohort="cohortI"):
    print(f"[fig5_ancestry] building data for {cohort} …")
    vt   = build_variant_table(cohort)
    meta = load_metadata_groups(cohort)

    anc_map = (meta[["sample_id", "genetically_determined"]]
               .drop_duplicates(subset=["sample_id"])
               .set_index("sample_id")["genetically_determined"])
    vt["ancestry"] = vt["sample_id"].map(anc_map)

    anc_n = (meta.drop_duplicates(subset=["sample_id"])
             .groupby("genetically_determined")["sample_id"].nunique())
    valid_anc = sorted([a for a, n in anc_n.items() if n >= 10 and pd.notna(a)])
    print(f"  Ancestry groups (≥10 samples): {valid_anc}")

    # -----------------------------------------------------------------------
    # 5A — grouped bar chart (no error bars)
    # -----------------------------------------------------------------------
    rows = []
    for vc in VC_TRACKS:
        vc_filter = _vc_carrier_filter(vt, vc)
        for anc in valid_anc:
            anc_samples = set(
                meta[meta["genetically_determined"] == anc]
                .drop_duplicates(subset=["sample_id"])["sample_id"]
            )
            n_total    = len(anc_samples)
            n_carriers = vc_filter[vc_filter["ancestry"] == anc]["sample_id"].nunique()
            pct        = n_carriers / n_total * 100 if n_total > 0 else 0.0
            rows.append({
                "ancestry":      anc,
                "variant_class": vc,
                "n_carriers":    n_carriers,
                "n_total":       n_total,
                "pct":           pct,
            })

    df5a = pd.DataFrame(rows)
    save_table(df5a, "fig5A_ancestry_bar_data", cohort=cohort)

    n_vc    = len(VC_TRACKS)
    bar_w   = 0.8 / n_vc
    x       = np.arange(len(valid_anc))
    offsets = {vc: (i - (n_vc - 1) / 2.0) * bar_w for i, vc in enumerate(VC_TRACKS)}

    # x-axis labels with N
    x_labels = [f"{a}\n(N={int(anc_n.get(a, 0)):,})" for a in valid_anc]

    fig, ax = plt.subplots(figsize=(max(16, len(valid_anc) * 2.8 + 2), 7.0))

    for vc in VC_TRACKS:
        sub = df5a[df5a["variant_class"] == vc].set_index("ancestry").reindex(valid_anc)
        vals = sub["pct"].fillna(0).values
        xi   = x + offsets[vc]

        ax.bar(xi, vals, width=bar_w - 0.02,
               color=COLORS[vc], label=vc, alpha=0.95, edgecolor="none")

        # horizontal number label above bar (no "%" sign)
        for xi_val, v in zip(xi, vals):
            if v > 0:
                ax.text(xi_val, v + 0.03, f"{v:.1f}",
                        ha="center", va="bottom", fontsize=15, color="#222222")

    ymax_5a = df5a["pct"].max()
    ax.set_ylim(0, ymax_5a * 1.20)
    ax.set_xticks(x)
    ax.set_xticklabels(x_labels, fontsize=19)
    ax.tick_params(axis="y", labelsize=18)
    ax.set_xlabel("Ancestry", fontsize=24)
    ax.set_ylabel("Carrier frequency (%)", fontsize=24)
    ax.set_title(f"Carrier frequency by ancestry — {cohort_label(cohort)}",
                 fontsize=22, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.xaxis.grid(False)
    ax.yaxis.grid(False)

    legend_patches = [mpatches.Patch(color=COLORS[vc], label=vc) for vc in VC_TRACKS]
    ax.legend(handles=legend_patches, fontsize=17,
              bbox_to_anchor=(1.01, 1), loc="upper left",
              framealpha=0.9, edgecolor="#CCCCCC")

    save_fig(fig, "fig5A_ancestry_bar", cohort=cohort)

    # -----------------------------------------------------------------------
    # 5B — gene × ancestry heatmap (seaborn)
    # -----------------------------------------------------------------------
    genes = ALL_GENES_ORDERED

    for vc in VC_TRACKS:
        vc_filter = _vc_carrier_filter(vt, vc)
        matrix        = np.full((len(genes), len(valid_anc)), np.nan)
        n_carrier_mat = np.zeros((len(genes), len(valid_anc)), dtype=int)
        n_total_mat   = np.zeros((len(genes), len(valid_anc)), dtype=int)

        for gi, gene in enumerate(genes):
            gene_carriers = set(vc_filter[vc_filter["gene"] == gene]["sample_id"].unique())
            for ai, anc in enumerate(valid_anc):
                anc_samples = set(
                    meta[meta["genetically_determined"] == anc]
                    .drop_duplicates(subset=["sample_id"])["sample_id"]
                )
                n_total    = len(anc_samples)
                n_carriers = len(gene_carriers & anc_samples)
                n_total_mat[gi, ai]   = n_total
                n_carrier_mat[gi, ai] = n_carriers
                if n_total >= 10:
                    matrix[gi, ai] = n_carriers / n_total * 100

        # Annotate cells with carrier freq ≥ 0.05% (avoids "0.0" from tiny values)
        annot = np.full(matrix.shape, "", dtype=object)
        for gi in range(len(genes)):
            for ai in range(len(valid_anc)):
                if not np.isnan(matrix[gi, ai]) and matrix[gi, ai] >= 0.05:
                    annot[gi, ai] = f"{matrix[gi, ai]:.1f}"

        # Chi-squared test per gene (carrier vs non-carrier across ancestries)
        pvals = []
        for gi_raw in range(len(genes)):
            table = np.vstack([n_carrier_mat[gi_raw],
                               n_total_mat[gi_raw] - n_carrier_mat[gi_raw]])
            if table.sum() == 0 or (table[0] > 0).sum() < 2:
                pvals.append(1.0)
                continue
            try:
                _, p, _, _ = chi2_contingency(table)
            except Exception:
                p = 1.0
            pvals.append(p)
        _, qvals, _, _ = multipletests(pvals, method="fdr_bh")
        gene_qvals = dict(zip(genes, qvals))

        vc_slug = vc.replace(" ", "_").replace("/", "").replace("(", "").replace(")", "")
        rows_5b = []
        for gi, gene in enumerate(genes):
            for ai, anc in enumerate(valid_anc):
                rows_5b.append({
                    "gene":          gene,
                    "syndrome":      GENE_TO_SYNDROME.get(gene, ""),
                    "ancestry":      anc,
                    "variant_class": vc,
                    "n_carriers":    int(n_carrier_mat[gi, ai]),
                    "n_total":       int(n_total_mat[gi, ai]),
                    "pct":           matrix[gi, ai],
                })
        save_table(pd.DataFrame(rows_5b),
                   f"fig5B_heatmap_{vc_slug}_data", cohort=cohort)

        # Convert zero/near-zero to NaN so seaborn renders them as white
        matrix_plot = matrix.copy()
        matrix_plot[matrix_plot < 0.05] = np.nan

        df_pivot  = pd.DataFrame(matrix_plot, index=genes, columns=valid_anc)
        mat_clean = np.nan_to_num(matrix, nan=0)   # use original for clustering

        # Hierarchical clustering on both axes
        if mat_clean.shape[0] > 1:
            row_order = leaves_list(
                linkage(pdist(mat_clean, metric="euclidean"), method="average")
            )
            df_pivot = df_pivot.iloc[row_order]
            annot    = annot[row_order, :]

        if mat_clean.shape[1] > 1:
            col_order = leaves_list(
                linkage(pdist(mat_clean.T, metric="euclidean"), method="average")
            )
            df_pivot = df_pivot.iloc[:, col_order]
            annot    = annot[:, col_order]

        fig_w = max(FIGW_DOUBLE, len(valid_anc) * 1.1 + 3)
        fig_h = max(8, len(genes) * 0.38 + 2)
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))

        sns.heatmap(
            df_pivot, ax=ax,
            cmap="Reds", vmin=0, vmax=HEATMAP_VMAX,
            linewidths=0.5, linecolor="white",
            annot=annot, fmt="",
            annot_kws={"size": 10, "color": "black"},
            cbar_kws={
                "label": f"Carrier Frequency (%) — values >{HEATMAP_VMAX:.0f}% annotated",
                "shrink": 0.5,
            },
        )

        yticklabels = ax.get_yticklabels()
        for lbl in yticklabels:
            if gene_qvals.get(lbl.get_text(), 1.0) < 0.05:
                lbl.set_fontweight("bold")
        ax.set_yticklabels(yticklabels, rotation=0, fontsize=12)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize=13)
        ax.set_ylabel("Gene", fontsize=16)
        ax.set_xlabel("Ancestry", fontsize=16)
        ax.set_title(
            f"Gene × Ancestry carrier frequency ({vc}) — {cohort_label(cohort)}",
            fontsize=14, fontweight="bold", pad=10,
        )

        fname = {
            VC_ACMG:    "fig5B_heatmap_ACMGplp",
            VC_CLINVAR: "fig5B_heatmap_ClinVarPLP",
            VC_AM:      "fig5B_heatmap_AMcalibrated",
            VC_AM_ONLY: "fig5B_heatmap_AMcalibratedNotPLP",
        }[vc]
        save_fig(fig, fname, cohort=cohort)

    print(f"[fig5_ancestry] done — {cohort}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", default="cohortI", choices=["cohortI", "cohortII"])
    args = parser.parse_args()
    make(args.cohort)


if __name__ == "__main__":
    main()
