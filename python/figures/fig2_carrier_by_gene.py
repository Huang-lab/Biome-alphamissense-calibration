"""
Fig 2 — Carrier frequency by gene (vertical grouped bar chart).

x-axis: 28 genes grouped by syndrome
y-axis: carrier frequency (%)
3 bars per gene: ACMG P/LP, AM_calibrated, AM_calibrated not P/LP
Syndrome labels shown above each gene group.

Usage:
    python -m python.figures.fig2_carrier_by_gene --cohort cohortI
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

from python.figures.common import (
    build_variant_table, load_metadata,
    COLORS, VC_ACMG, VC_AM, VC_AM_ONLY,
    SYNDROME_ORDER, SYNDROME_LABELS, SYNDROMES,
    FIGW_DOUBLE,
    save_fig, save_table, cohort_label,
)


def make(cohort="cohortI"):
    print(f"[fig2_carrier_by_gene] building data for {cohort} …")
    vt   = build_variant_table(cohort)
    meta = load_metadata(cohort)
    id_col  = "SINAI_ID" if cohort == "cohortI" else "MASKED_MRN"
    total_n = meta[id_col].nunique()

    acmg_by_gene   = (vt[vt["in_acmg"] == True]
                      .groupby("gene")["sample_id"].nunique().rename("acmg"))
    am_by_gene     = (vt[vt["in_am"] == True]
                      .groupby("gene")["sample_id"].nunique().rename("am"))
    amonly_by_gene = (vt[vt["variant_class"] == VC_AM_ONLY]
                      .groupby("gene")["sample_id"].nunique().rename("amonly"))

    rows = []
    for syn_key in SYNDROME_ORDER:
        genes     = SYNDROMES[syn_key]["genes"]
        syn_label = SYNDROME_LABELS[syn_key]
        for gene in genes:
            rows.append({
                "gene":      gene,
                "syndrome":  syn_key,
                "syn_label": syn_label,
                "acmg":      acmg_by_gene.get(gene, 0),
                "am":        am_by_gene.get(gene, 0),
                "amonly":    amonly_by_gene.get(gene, 0),
            })

    df = pd.DataFrame(rows)
    df["pct_acmg"]   = df["acmg"]   / total_n * 100
    df["pct_am"]     = df["am"]     / total_n * 100
    df["pct_amonly"] = df["amonly"] / total_n * 100

    n_genes = len(df)
    bar_w   = 0.24
    group_w = 1.0   # x-spacing per gene group

    # -----------------------------------------------------------------------
    # Plot — vertical bars, x=gene groups
    # -----------------------------------------------------------------------
    fig_w = max(20, n_genes * group_w * 0.7 + 3)
    fig, ax = plt.subplots(figsize=(fig_w, 8))

    vc_configs = [
        ("pct_acmg",   VC_ACMG,    -bar_w),
        ("pct_am",     VC_AM,       0.0),
        ("pct_amonly", VC_AM_ONLY,  bar_w),
    ]

    x_centers = np.arange(n_genes) * group_w

    for col, vc, offset in vc_configs:
        xs   = x_centers + offset
        vals = df[col].values
        ax.bar(xs, vals, width=bar_w - 0.02,
               color=COLORS[vc], label=vc, alpha=0.9,
               edgecolor="white", linewidth=0.3)

    # Gene x-tick labels at 45°
    ax.set_xticks(x_centers)
    ax.set_xticklabels(df["gene"].tolist(), rotation=45, ha="right", fontsize=17)
    ax.tick_params(axis="y", labelsize=15)

    # Syndrome group separators + syndrome labels above bars
    prev_syn   = None
    syn_starts = {}   # syn_key → first gene index
    syn_ends   = {}   # syn_key → last gene index
    for i, row in df.iterrows():
        pos = list(df.index).index(i)
        s   = row["syndrome"]
        if s not in syn_starts:
            syn_starts[s] = pos
        syn_ends[s] = pos

    ymax = max(df[["pct_acmg", "pct_am", "pct_amonly"]].max()) * 1.05 or 1.0

    drawn_i = 0
    for syn_key in SYNDROME_ORDER:
        if syn_key not in syn_starts:
            continue
        s_start = syn_starts[syn_key]
        s_end   = syn_ends[syn_key]

        # Dashed vertical separator before each syndrome (except first)
        if s_start > 0:
            sep_x = (x_centers[s_start - 1] + x_centers[s_start]) / 2
            ax.axvline(sep_x, color="#CCCCCC", linewidth=0.8, linestyle="--")

        # Syndrome label centered above the group; stagger height so adjacent
        # single-gene syndromes (e.g. VHL, PTEN-HTS) don't overlap
        mid_x  = (x_centers[s_start] + x_centers[s_end]) / 2
        y_frac = 1.02 if drawn_i % 2 == 0 else 1.09
        ax.text(mid_x, ymax * y_frac, SYNDROME_LABELS[syn_key],
                ha="center", va="bottom", fontsize=15,
                color="#555555", style="italic")
        drawn_i += 1

    ax.set_ylabel("Carrier frequency (%)", fontsize=20)
    ax.set_title(f"Carrier frequency by gene — {cohort_label(cohort)}",
                 fontsize=20, fontweight="bold")
    ax.set_ylim(0, ymax * 1.22)
    ax.xaxis.grid(False)
    ax.yaxis.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    legend_patches = [
        mpatches.Patch(color=COLORS[VC_ACMG],    label=VC_ACMG),
        mpatches.Patch(color=COLORS[VC_AM],       label=VC_AM),
        mpatches.Patch(color=COLORS[VC_AM_ONLY],  label=VC_AM_ONLY),
    ]
    ax.legend(handles=legend_patches, fontsize=15,
              bbox_to_anchor=(1.01, 1), loc="upper left",
              framealpha=0.9, edgecolor="#CCCCCC")

    plt.tight_layout()
    save_fig(fig, "fig2_carrier_by_gene", cohort=cohort)

    # Save data table
    out_df = df[["gene", "syn_label", "acmg", "am", "amonly",
                 "pct_acmg", "pct_am", "pct_amonly"]].copy()
    out_df.columns = ["gene", "syndrome", "n_ACMG_PLP", "n_AM_calibrated",
                      "n_AM_not_PLP", "pct_ACMG_PLP", "pct_AM_calibrated", "pct_AM_not_PLP"]
    out_df["total_n"] = total_n
    save_table(out_df, "fig2_carrier_by_gene_data", cohort=cohort)
    print(f"[fig2_carrier_by_gene] done — {cohort}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", default="cohortI", choices=["cohortI", "cohortII"])
    args = parser.parse_args()
    make(args.cohort)


if __name__ == "__main__":
    main()
