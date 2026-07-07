"""
Fig 1 — Summary of variant counts and carrier frequencies.

Panels saved individually per cohort:
  fig1A_venn_{cohort}.png           – overlap of ACMG P/LP vs AM_calibrated (Venn)
  fig1B_variant_counts_{cohort}.png – unique variant counts per class
  fig1C_carrier_freq_{cohort}.png   – % of cohort with ≥1 carrier per class
  fig1D_carrier_by_group_{cohort}.png – carrier freq per phenotype group (3 subplots)

Usage:
    python -m python.figures.fig1_summary --cohort cohortI
"""

import argparse
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd

from python.figures.common import (
    build_variant_table, load_metadata, load_metadata_groups,
    COLORS, VARIANT_CLASSES, VC_ACMG, VC_AM, VC_AM_ONLY, VC_BOTH,
    FIGW_SINGLE, FIGW_DOUBLE,
    save_fig, save_table, cohort_label,
)

N_MIN_GROUP = 20   # minimum group size to show in fig1D barplot


def _count_unique_variants(vt, vc_list):
    sub = vt[vt["variant_class"].isin(vc_list)]
    return sub[["chrom", "pos", "ref", "alt"]].drop_duplicates().shape[0]


def make(cohort="cohortI"):
    print(f"[fig1_summary] building variant table for {cohort} …")
    vt   = build_variant_table(cohort)
    meta = load_metadata(cohort)
    id_col  = "SINAI_ID" if cohort == "cohortI" else "MASKED_MRN"
    total_n = meta[id_col].nunique()
    print(f"  Total samples: {total_n:,}")

    # -----------------------------------------------------------------------
    # Compute counts
    # -----------------------------------------------------------------------
    acmg_vars   = _count_unique_variants(vt, [VC_ACMG, VC_BOTH])
    am_vars     = _count_unique_variants(vt, [VC_AM_ONLY, VC_BOTH])
    amonly_vars = _count_unique_variants(vt, [VC_AM_ONLY])
    both_vars   = _count_unique_variants(vt, [VC_BOTH])

    acmg_carriers   = vt[vt["variant_class"].isin([VC_ACMG, VC_BOTH])]["sample_id"].nunique()
    am_carriers     = vt[vt["variant_class"].isin([VC_AM_ONLY, VC_BOTH])]["sample_id"].nunique()
    amonly_carriers = vt[vt["variant_class"] == VC_AM_ONLY]["sample_id"].nunique()

    labels    = [VC_ACMG, VC_AM, VC_AM_ONLY]
    colors    = [COLORS[l] for l in labels]
    var_counts  = [acmg_vars, am_vars, amonly_vars]
    carr_counts = [acmg_carriers, am_carriers, amonly_carriers]
    carr_pcts   = [c / total_n * 100 for c in carr_counts]

    # -----------------------------------------------------------------------
    # Panel A — Venn diagram
    # -----------------------------------------------------------------------
    fig_a, ax = plt.subplots(figsize=(6, 5))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7)
    ax.set_aspect("equal")
    ax.axis("off")

    c1 = patches.Circle((3.8, 3.5), 2.6,
                         facecolor=COLORS[VC_ACMG] + "55",
                         edgecolor=COLORS[VC_ACMG], linewidth=2)
    c2 = patches.Circle((6.2, 3.5), 2.6,
                         facecolor=COLORS[VC_AM] + "55",
                         edgecolor=COLORS[VC_AM], linewidth=2)
    ax.add_patch(c1)
    ax.add_patch(c2)

    acmg_only_vars = acmg_vars - both_vars
    ax.text(2.4, 3.5, f"{acmg_only_vars:,}\nvariants",
            ha="center", va="center", fontsize=14, fontweight="bold",
            color=COLORS[VC_ACMG])
    ax.text(5.0, 3.5, f"{both_vars:,}\nboth",
            ha="center", va="center", fontsize=13, color="#555555")
    ax.text(7.6, 3.5, f"{amonly_vars:,}\nvariants",
            ha="center", va="center", fontsize=14, fontweight="bold",
            color=COLORS[VC_AM])
    ax.text(2.4, 6.3, "ACMG P/LP", ha="center", va="center",
            fontsize=13, color=COLORS[VC_ACMG], fontweight="bold")
    ax.text(7.6, 6.3, "AM_calibrated", ha="center", va="center",
            fontsize=13, color=COLORS[VC_AM], fontweight="bold")
    ax.set_title(f"Variant overlap — {cohort_label(cohort)}", fontsize=14, fontweight="bold", pad=4)
    plt.tight_layout()
    save_fig(fig_a, "fig1A_venn", cohort=cohort)

    venn_data = pd.DataFrame({
        "category":          ["ACMG P/LP only", "AM_calibrated only", "both"],
        "n_unique_variants": [acmg_only_vars, amonly_vars, both_vars],
    })
    save_table(venn_data, "fig1A_venn_data", cohort=cohort)

    # -----------------------------------------------------------------------
    # Panel B — Variant counts bar chart
    # -----------------------------------------------------------------------
    fig_b, ax = plt.subplots(figsize=(7, 6))
    x    = np.arange(len(labels))
    bars = ax.bar(x, var_counts, color=colors, width=0.55, edgecolor="white", linewidth=0.5)
    for bar, val in zip(bars, var_counts):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max(var_counts) * 0.01,
                f"{val:,}", ha="center", va="bottom", fontsize=14, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=13, rotation=45, ha="right")
    ax.set_ylabel("Unique variants (n)", fontsize=16)
    ax.set_title(f"Variant counts — {cohort_label(cohort)}", fontsize=14, fontweight="bold")
    ax.yaxis.grid(True, linestyle=":", alpha=0.5, color="#CCCCCC")
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    save_fig(fig_b, "fig1B_variant_counts", cohort=cohort)

    # -----------------------------------------------------------------------
    # Panel C — Carrier frequency
    # -----------------------------------------------------------------------
    fig_c, ax = plt.subplots(figsize=(7, 6))
    bars = ax.bar(x, carr_pcts, color=colors, width=0.55, edgecolor="white", linewidth=0.5)
    for bar, pct, cnt in zip(bars, carr_pcts, carr_counts):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max(carr_pcts) * 0.02,
                f"{pct:.1f}%\n(n={cnt:,})",
                ha="center", va="bottom", fontsize=12)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=13, rotation=45, ha="right")
    ax.set_ylabel("Carrier frequency (%)", fontsize=16)
    ax.set_title(f"Carrier frequency — {cohort_label(cohort)}", fontsize=14, fontweight="bold",
                 pad=20)
    ax.yaxis.grid(True, linestyle=":", alpha=0.5, color="#CCCCCC")
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig_c.tight_layout(rect=[0, 0, 1, 0.92])
    save_fig(fig_c, "fig1C_carrier_freq", cohort=cohort)

    summary = pd.DataFrame({
        "variant_class":     labels,
        "n_unique_variants": var_counts,
        "n_carriers":        carr_counts,
        "pct_carriers":      [round(p, 4) for p in carr_pcts],
        "total_n":           total_n,
    })
    save_table(summary, "fig1_summary_data", cohort=cohort)

    # -----------------------------------------------------------------------
    # Panel D — Carrier frequency by phenotype group (3 variant classes)
    # -----------------------------------------------------------------------
    meta_grp  = load_metadata_groups(cohort)
    group_n   = meta_grp.groupby("Group")["sample_id"].nunique()
    all_groups = [g for g, n in group_n.items() if n > 50]

    rows_d = []
    for vc, vc_filter in [
        (VC_ACMG,    vt[vt["in_acmg"] == True]),
        (VC_AM,      vt[vt["in_am"]   == True]),
        (VC_AM_ONLY, vt[vt["variant_class"] == VC_AM_ONLY]),
    ]:
        vc_carriers_set = set(vc_filter["sample_id"].unique())
        overall_pct     = len(vc_carriers_set) / total_n * 100
        for grp in all_groups:
            grp_samples = set(
                meta_grp[meta_grp["Group"] == grp]["sample_id"].unique()
            )
            n_total    = len(grp_samples)
            n_carriers = len(vc_carriers_set & grp_samples)
            pct        = n_carriers / n_total * 100 if n_total > 0 else 0.0
            rows_d.append({
                "Group":         grp,
                "variant_class": vc,
                "n_total":       n_total,
                "n_carriers":    n_carriers,
                "pct":           pct,
                "overall_pct":   overall_pct,
            })

    df_d = pd.DataFrame(rows_d)
    save_table(df_d, "fig1D_carrier_by_group_data", cohort=cohort)

    n_groups  = len(all_groups)
    fig_h     = max(10, n_groups * 0.38 + 2)
    group_n_lookup = meta_grp.groupby("Group")["sample_id"].nunique().to_dict()

    vc_slugs_d = {
        VC_ACMG:    "fig1D_ACMGplp",
        VC_AM:      "fig1D_AMcalibrated",
        VC_AM_ONLY: "fig1D_AMcalibratedNotPLP",
    }

    for vc in [VC_ACMG, VC_AM, VC_AM_ONLY]:
        sub = df_d[df_d["variant_class"] == vc].copy()
        sub = sub.sort_values("pct", ascending=False).reset_index(drop=True)

        y_pos = list(range(len(sub)))

        fig_d, ax_d = plt.subplots(figsize=(14, fig_h))
        ax_d.barh(y_pos, sub["pct"], color=COLORS[vc], alpha=0.85, edgecolor="none")

        for i, row in sub.iterrows():
            ax_d.text(row["pct"] + max(sub["pct"]) * 0.01, list(sub.index).index(i),
                      f"{row['pct']:.1f}%",
                      va="center", ha="left", fontsize=10, color="#444444")

        ax_d.set_yticks(y_pos)
        y_labels_d = [f"{row['Group']}\n(N={group_n_lookup.get(row['Group'], row['n_total']):,})"
                      for _, row in sub.iterrows()]
        ax_d.set_yticklabels(y_labels_d, fontsize=12)
        ax_d.invert_yaxis()
        ax_d.set_xlabel("Carrier Frequency (%)", fontsize=16)
        ax_d.set_title(f"{vc} — {cohort_label(cohort)}", fontsize=13,
                       fontweight="bold", color=COLORS[vc])
        ax_d.spines["top"].set_visible(False)
        ax_d.spines["right"].set_visible(False)
        ax_d.xaxis.grid(False)
        ax_d.yaxis.grid(False)

        plt.tight_layout()
        save_fig(fig_d, vc_slugs_d[vc], cohort=cohort)

    print(f"[fig1_summary] done — {cohort}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", default="cohortI", choices=["cohortI", "cohortII"])
    args = parser.parse_args()
    make(args.cohort)


if __name__ == "__main__":
    main()
