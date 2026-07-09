"""
Fig 4 — Carrier frequency by phenotype Group.

4A: Horizontal barplot — carrier freq per phenotype Group (≥45 samples)
    3 grouped bars per group: ACMG P/LP, AM_calibrated, AM_calibrated not P/LP
    No error bars; y-axis labels include (N=xxx); plain number at bar end.
    Sorted highest → lowest (highest at top); legend outside plot.

4B: Heatmap — gene × Group (seaborn style, one per variant class)
    white → dark red (0–5%), hierarchically clustered both axes.
    Annotate only cells > vmax (5%).

Usage:
    python -m python.figures.fig4_barplot_heatmap --cohort cohortI
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

from python.figures.common import (
    build_variant_table, load_metadata_groups,
    COLORS, VC_ACMG, VC_CLINVAR, VC_AM, VC_AM_ONLY,
    ALL_GENES_ORDERED, GENE_TO_SYNDROME,
    CONTROL_GROUPS_RAW,
    FIGW_DOUBLE, N_MIN_PHENOTYPE,
    save_fig, save_table, cohort_label,
)

# Ordered variant tracks shown in all fig4 panels.
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

HEATMAP_VMAX = 3.0   # colorbar max 0–3%; annotate only cells >3%


def make_barplot(vt, meta, case_groups, cohort):
    """4A: horizontal grouped barplot — highest group at top, plain number labels."""
    group_n = meta.groupby("Group")["sample_id"].nunique().to_dict()

    rows = []
    for vc in VC_TRACKS:
        vc_filter = _vc_carrier_filter(vt, vc)
        vc_carriers = set(vc_filter["sample_id"].unique())
        for grp in case_groups:
            grp_samples = set(meta[meta["Group"] == grp]["sample_id"].unique())
            n_total    = len(grp_samples)
            n_carriers = len(vc_carriers & grp_samples)
            pct        = n_carriers / n_total * 100 if n_total > 0 else 0.0
            rows.append({
                "Group":         grp,
                "variant_class": vc,
                "n_carriers":    n_carriers,
                "n_total":       n_total,
                "pct":           pct,
            })

    df = pd.DataFrame(rows)
    save_table(df, "fig4A_barplot_data", cohort=cohort)

    # Sort descending by AM_calibrated freq → highest group at top (with invert_yaxis)
    am_order = (df[df["variant_class"] == VC_AM]
                .sort_values("pct", ascending=False)["Group"].tolist())

    y_labels = [f"{grp} (N={group_n.get(grp, 0):,})" for grp in am_order]

    fig_h  = max(8, len(case_groups) * 0.42 + 2)
    fig, ax = plt.subplots(figsize=(18, fig_h))

    n_vc    = len(VC_TRACKS)
    bar_h   = 0.8 / n_vc   # total group height ~0.8, split evenly across tracks
    offsets = {vc: (i - (n_vc - 1) / 2.0) * bar_h for i, vc in enumerate(VC_TRACKS)}
    y_pos   = {grp: i for i, grp in enumerate(am_order)}

    for vc in VC_TRACKS:
        sub = df[df["variant_class"] == vc].copy()
        sub["y"] = sub["Group"].map(y_pos) + offsets[vc]
        sub = sub.dropna(subset=["y"])

        ax.barh(sub["y"], sub["pct"], height=bar_h - 0.02,
                color=COLORS[vc], label=vc, alpha=0.95, edgecolor="none")

        # Plain number at bar end (no % sign)
        for _, row in sub.iterrows():
            if row["pct"] > 0:
                ax.text(row["pct"] + 0.05, row["y"],
                        f"{row['pct']:.1f}",
                        va="center", ha="left", fontsize=13, color="#333333")

    # x headroom so bar-end labels don't run into the legend
    xmax = df["pct"].max()
    ax.set_xlim(0, xmax * 1.18)

    ax.set_yticks(list(range(len(am_order))))
    ax.set_yticklabels(y_labels, fontsize=16)
    ax.tick_params(axis="x", labelsize=15)
    ax.invert_yaxis()   # puts highest group (index 0) at the top

    ax.set_xlabel("Carrier frequency (%)", fontsize=20)
    ax.set_title(f"Carrier frequency by phenotype group — {cohort_label(cohort)}",
                 fontsize=20, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.xaxis.grid(False)
    ax.yaxis.grid(False)

    legend_patches = [mpatches.Patch(color=COLORS[vc], label=vc) for vc in VC_TRACKS]
    ax.legend(handles=legend_patches, fontsize=16,
              bbox_to_anchor=(1.01, 1), loc="upper left",
              framealpha=0.95, edgecolor="#CCCCCC")

    plt.tight_layout()
    save_fig(fig, "fig4A_barplot", cohort=cohort)


def make_heatmap(vt, meta, case_groups, cohort):
    """4B: gene × Group heatmap (seaborn, clustered, vmax=5%)."""
    genes = ALL_GENES_ORDERED

    for vc in VC_TRACKS:
        vc_filter = _vc_carrier_filter(vt, vc)
        matrix        = np.full((len(genes), len(case_groups)), np.nan)
        n_carrier_mat = np.zeros((len(genes), len(case_groups)), dtype=int)
        n_total_mat   = np.zeros((len(genes), len(case_groups)), dtype=int)

        for gi, gene in enumerate(genes):
            gene_carriers = set(vc_filter[vc_filter["gene"] == gene]["sample_id"].unique())
            for ci, grp in enumerate(case_groups):
                grp_samples = set(meta[meta["Group"] == grp]["sample_id"].unique())
                n_total    = len(grp_samples)
                n_carriers = len(gene_carriers & grp_samples)
                n_total_mat[gi, ci]   = n_total
                n_carrier_mat[gi, ci] = n_carriers
                if n_total >= 10:
                    matrix[gi, ci] = n_carriers / n_total * 100

        # Annotate only cells > HEATMAP_VMAX (saturated / clipped values)
        annot = np.full(matrix.shape, "", dtype=object)
        for gi in range(len(genes)):
            for ci in range(len(case_groups)):
                if not np.isnan(matrix[gi, ci]) and matrix[gi, ci] > HEATMAP_VMAX:
                    annot[gi, ci] = f"{matrix[gi, ci]:.1f}"

        # Save data
        vc_slug  = vc.replace(" ", "_").replace("/", "").replace("(", "").replace(")", "")
        rows_out = []
        for gi, gene in enumerate(genes):
            for ci, grp in enumerate(case_groups):
                rows_out.append({
                    "gene":          gene,
                    "syndrome":      GENE_TO_SYNDROME.get(gene, ""),
                    "Group":         grp,
                    "variant_class": vc,
                    "n_carriers":    int(n_carrier_mat[gi, ci]),
                    "n_total":       int(n_total_mat[gi, ci]),
                    "pct":           matrix[gi, ci],
                })
        save_table(pd.DataFrame(rows_out),
                   f"fig4B_heatmap_{vc_slug}_data", cohort=cohort)

        # Hierarchical clustering on both axes
        df_pivot    = pd.DataFrame(matrix, index=genes, columns=case_groups)
        mat_clean   = np.nan_to_num(matrix, nan=0)

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

        fig_w = max(FIGW_DOUBLE + 2, len(case_groups) * 0.42 + 3)
        fig_h = max(8, len(genes) * 0.38 + 2)
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))

        sns.heatmap(
            df_pivot, ax=ax,
            cmap="Reds", vmin=0, vmax=HEATMAP_VMAX,
            linewidths=0.5, linecolor="white",
            annot=annot, fmt="",
            annot_kws={"size": 12, "color": "black"},
            cbar_kws={
                "label": f"Carrier Frequency (%) — values >{HEATMAP_VMAX:.0f}% annotated",
                "shrink": 0.5,
            },
        )

        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=15)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=13)
        ax.set_ylabel("Gene", fontsize=19)
        ax.set_xlabel("Phenotype Group", fontsize=19)
        ax.set_title(
            f"Gene × Phenotype carrier frequency ({vc}) — {cohort_label(cohort)}",
            fontsize=17, fontweight="bold", pad=10,
        )
        cbar = ax.collections[0].colorbar
        cbar.ax.tick_params(labelsize=13)
        cbar.set_label(cbar.ax.get_ylabel(), fontsize=14)

        fname = {
            VC_ACMG:    "fig4B_heatmap_ACMGplp",
            VC_CLINVAR: "fig4B_heatmap_ClinVarPLP",
            VC_AM:      "fig4B_heatmap_AMcalibrated",
            VC_AM_ONLY: "fig4B_heatmap_AMcalibratedNotPLP",
        }[vc]
        save_fig(fig, fname, cohort=cohort)


def make(cohort="cohortI"):
    print(f"[fig4_barplot_heatmap] building data for {cohort} …")
    vt   = build_variant_table(cohort)
    meta = load_metadata_groups(cohort)

    group_n = meta.groupby("Group")["sample_id"].nunique()
    # Use N≥20 for barplot/heatmap to show more phenotype groups
    case_groups = sorted([
        g for g, n in group_n.items()
        if n >= 20 and g not in CONTROL_GROUPS_RAW
    ])
    print(f"  Case groups with ≥{N_MIN_PHENOTYPE} samples: {len(case_groups)}")

    make_barplot(vt, meta, case_groups, cohort)
    make_heatmap(vt, meta, case_groups, cohort)
    print(f"[fig4_barplot_heatmap] done — {cohort}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", default="cohortI", choices=["cohortI", "cohortII"])
    args = parser.parse_args()
    make(args.cohort)


if __name__ == "__main__":
    main()
