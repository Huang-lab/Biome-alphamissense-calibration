"""
Fig 1 — Summary of variant counts and carrier frequencies.

Panels saved individually per cohort:
  fig1A_venn_{cohort}.png           – 3-set overlap of ACMG P/LP vs AM_calibrated vs ClinVar P/LP (Venn)
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
import numpy as np
import pandas as pd

from python.figures.common import (
    build_variant_table, load_metadata, load_metadata_groups,
    COLORS, VARIANT_CLASSES, VC_ACMG, VC_CLINVAR, VC_AM, VC_AM_ONLY, VC_BOTH,
    FIGW_SINGLE, FIGW_DOUBLE,
    save_fig, save_table, cohort_label,
)

N_MIN_GROUP = 20   # minimum group size to show in fig1D barplot


def _count_unique_variants(vt, vc_list):
    sub = vt[vt["variant_class"].isin(vc_list)]
    return sub[["chrom", "pos", "ref", "alt"]].drop_duplicates().shape[0]


def _venn_region_counts(vt):
    """Unique-variant counts for the 7 regions of the ACMG / AM_calibrated /
    ClinVar Venn. Dedupe by variant key and OR the membership flags."""
    key = ["chrom", "pos", "ref", "alt"]
    g = vt.groupby(key)[["in_acmg", "in_am", "in_clinvar"]].any()
    A, B, C = g["in_acmg"], g["in_am"], g["in_clinvar"]
    return {
        "Abc": int((A & ~B & ~C).sum()),   # ACMG only
        "aBc": int((~A & B & ~C).sum()),   # AM only
        "ABc": int((A & B & ~C).sum()),    # ACMG & AM
        "abC": int((~A & ~B & C).sum()),   # ClinVar only
        "AbC": int((A & ~B & C).sum()),    # ACMG & ClinVar
        "aBC": int((~A & B & C).sum()),    # AM & ClinVar
        "ABC": int((A & B & C).sum()),     # all three
    }


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
    reg = _venn_region_counts(vt)
    acmg_vars    = reg["Abc"] + reg["ABc"] + reg["AbC"] + reg["ABC"]   # any ACMG
    clinvar_vars = reg["abC"] + reg["AbC"] + reg["aBC"] + reg["ABC"]   # any ClinVar
    am_vars      = reg["aBc"] + reg["ABc"] + reg["aBC"] + reg["ABC"]   # any AM
    amonly_vars  = reg["aBc"]                                           # AM only (no clinical)

    acmg_carriers    = vt[vt["in_acmg"] == True]["sample_id"].nunique()
    clinvar_carriers = vt[vt["in_clinvar"] == True]["sample_id"].nunique()
    am_carriers      = vt[vt["in_am"] == True]["sample_id"].nunique()
    amonly_carriers  = vt[vt["variant_class"] == VC_AM_ONLY]["sample_id"].nunique()

    labels    = [VC_ACMG, VC_CLINVAR, VC_AM, VC_AM_ONLY]
    colors    = [COLORS[l] for l in labels]
    var_counts  = [acmg_vars, clinvar_vars, am_vars, amonly_vars]
    carr_counts = [acmg_carriers, clinvar_carriers, am_carriers, amonly_carriers]
    carr_pcts   = [c / total_n * 100 for c in carr_counts]

    # -----------------------------------------------------------------------
    # Panel A — 3-set Venn (ACMG P/LP vs AM_calibrated vs ClinVar P/LP)
    # -----------------------------------------------------------------------
    fig_a, ax = plt.subplots(figsize=(6.5, 5.5))
    # venn3 subset order: (Abc, aBc, ABc, abC, AbC, aBC, ABC) for sets (A,B,C).
    subsets = (reg["Abc"], reg["aBc"], reg["ABc"],
               reg["abC"], reg["AbC"], reg["aBC"], reg["ABC"])
    _drew_venn = False
    try:
        from matplotlib_venn import venn3, venn3_circles
        v = venn3(subsets=subsets,
                  set_labels=("ACMG P/LP", "AM_calibrated", "ClinVar P/LP"),
                  ax=ax)
        venn3_circles(subsets=subsets, ax=ax, linewidth=1.5)
        patch_colors = {"100": COLORS[VC_ACMG], "010": COLORS[VC_AM], "001": COLORS[VC_CLINVAR]}
        for pid, col in patch_colors.items():
            p = v.get_patch_by_id(pid)
            if p is not None:
                p.set_color(col)
                p.set_alpha(0.45)
        _drew_venn = True
    except Exception as e:  # noqa: BLE001 — fall back if matplotlib_venn absent
        print(f"  [warn] matplotlib_venn unavailable ({e}); drawing text fallback")
        ax.axis("off")
        ax.text(0.5, 0.5,
                "ACMG only: {Abc}\nAM only: {aBc}\nClinVar only: {abC}\n"
                "ACMG∩AM: {ABc}\nACMG∩ClinVar: {AbC}\nAM∩ClinVar: {aBC}\n"
                "all three: {ABC}".format(**reg),
                ha="center", va="center", fontsize=12, transform=ax.transAxes)
    ax.set_title(f"Variant overlap — {cohort_label(cohort)}",
                 fontsize=14, fontweight="bold", pad=4)
    plt.tight_layout()
    save_fig(fig_a, "fig1A_venn", cohort=cohort)

    venn_data = pd.DataFrame({
        "category": ["ACMG P/LP only", "AM_calibrated only", "ClinVar P/LP only",
                     "ACMG & AM", "ACMG & ClinVar", "AM & ClinVar", "all three"],
        "n_unique_variants": [reg["Abc"], reg["aBc"], reg["abC"],
                              reg["ABc"], reg["AbC"], reg["aBC"], reg["ABC"]],
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
        (VC_ACMG,    vt[vt["in_acmg"]    == True]),
        (VC_CLINVAR, vt[vt["in_clinvar"] == True]),
        (VC_AM,      vt[vt["in_am"]      == True]),
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
        VC_CLINVAR: "fig1D_ClinVarPLP",
        VC_AM:      "fig1D_AMcalibrated",
        VC_AM_ONLY: "fig1D_AMcalibratedNotPLP",
    }

    for vc in [VC_ACMG, VC_CLINVAR, VC_AM, VC_AM_ONLY]:
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
