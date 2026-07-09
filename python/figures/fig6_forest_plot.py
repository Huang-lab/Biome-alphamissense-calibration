"""
Fig 6 — Forest plot of logistic regression results.

One figure per variant class (3 total per cohort).
y-axis : phenotype label
x-axis : OR (log scale), reference line at 1.0
color  : syndrome (legend shown); alpha full = q_BH<0.05, faded = NS

Prerequisite: python python/run_regression.py --cohort cohortI

Usage:
    python -m python.figures.fig6_forest_plot --cohort cohortI
"""

import argparse
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import pandas as pd

from python.figures.common import (
    RESULTS_I, RESULTS_II,
    COLORS, VC_ACMG, VC_CLINVAR, VC_AM, VC_AM_ONLY,
    SYNDROME_ORDER, SYNDROME_LABELS, SYNDROME_COLORS, SYNDROMES,
    save_fig, save_table, cohort_label,
)

SIG_ALPHA   = 0.05
MARKER_SIZE = 11


def _load(cohort):
    out_dir = RESULTS_I if cohort == "cohortI" else RESULTS_II
    path    = os.path.join(out_dir, "regression_results.tsv")
    if not os.path.exists(path):
        raise FileNotFoundError(f"{path} not found — run run_regression.py first.")
    df = pd.read_csv(path, sep="\t")
    df = df[df["status"] == "ok"].copy()
    # Remove MUTYH (may still be in older regression output)
    df = df[df["syndrome_key"] != "MUTYH_Associated_Polyposis"].copy()
    return df


def _draw_panel(ax, sub, title):
    """Draw one forest panel: y=phenotype, color=syndrome. All rows are q<0.05."""
    if sub.empty:
        ax.set_visible(False)
        return

    # Sort ascending so highest ORs appear at top after invert
    sub = sub.sort_values("OR", ascending=True).reset_index(drop=True)

    y_pos = np.arange(len(sub))

    for i, row in sub.iterrows():
        syn_color = SYNDROME_COLORS.get(row["syndrome_key"], "#888888")

        # CI whisker
        ax.hlines(i, row["CI_low"], row["CI_high"],
                  color=syn_color, linewidth=2.0, alpha=1.0, zorder=2)
        # OR point
        ax.plot(row["OR"], i, "o",
                color=syn_color, markersize=MARKER_SIZE,
                markeredgewidth=1.0, markeredgecolor="black",
                alpha=1.0, zorder=3)

    # Reference line at OR = 1
    ax.axvline(1.0, color="#333333", linestyle="--", linewidth=1.0, alpha=0.7)

    # Keep rows from floating apart when there are only a few
    ax.set_ylim(-0.8, len(sub) - 0.2)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(sub["phenotype"].tolist(), fontsize=17)
    ax.set_xlabel("Odds Ratio (95% CI)", fontsize=18)
    ax.tick_params(axis="x", labelsize=15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.xaxis.grid(True, linestyle=":", alpha=0.4, color="#CCCCCC")
    ax.yaxis.grid(False)
    ax.set_axisbelow(True)

    # Linear x-axis starting at 0
    ax.set_xlim(left=0)

    return sub  # return sorted df so caller can build legend


def make(cohort="cohortI"):
    print(f"[fig6_forest_plot] loading regression results for {cohort} …")
    df = _load(cohort)
    print(f"  {len(df)} rows (status=ok, MUTYH excluded)")

    vc_slugs = {
        VC_ACMG:    "fig6_forest_ACMGplp",
        VC_CLINVAR: "fig6_forest_ClinVarPLP",
        VC_AM:      "fig6_forest_AMcalibrated",
        VC_AM_ONLY: "fig6_forest_AMcalibratedNotPLP",
    }

    for vc in [VC_ACMG, VC_CLINVAR, VC_AM, VC_AM_ONLY]:
        sub = df[(df["variant_class"] == vc) &
                 (~df["q_BH"].isna()) &
                 (df["q_BH"] < SIG_ALPHA)].copy()
        n_rows = len(sub)
        fig_h  = max(6, n_rows * 0.48 + 2.5)

        fig, ax = plt.subplots(figsize=(17, fig_h))
        sorted_sub = _draw_panel(ax, sub, vc)

        if sorted_sub is not None and not sorted_sub.empty:
            # Build syndrome legend — only syndromes present in this vc
            seen_syns = sorted_sub["syndrome_key"].unique()
            legend_handles = []
            for syn in [s for s in SYNDROME_ORDER if s in seen_syns]:
                gene_str = ", ".join(SYNDROMES[syn]["genes"])
                label    = f"{SYNDROME_LABELS[syn]} ({gene_str})"
                legend_handles.append(
                    mlines.Line2D([], [], color=SYNDROME_COLORS.get(syn, "#888"),
                                  marker="o", linestyle="None",
                                  markersize=MARKER_SIZE,
                                  label=label)
                )
            ax.legend(handles=legend_handles,
                      title="Syndrome (all q<0.05)",
                      title_fontsize=14,
                      bbox_to_anchor=(1.01, 1), loc="upper left",
                      fontsize=14, framealpha=0.9, edgecolor="#CCCCCC")

        fig.suptitle(
            f"Logistic regression — OR (95% CI)\n{vc} — {cohort_label(cohort)}",
            fontsize=17, fontweight="bold",
        )
        plt.tight_layout()
        save_fig(fig, vc_slugs[vc], cohort=cohort)

    save_table(df[["syndrome_key", "phenotype", "variant_class",
                   "OR", "CI_low", "CI_high", "p_value", "q_BH",
                   "n_cases_total", "n_controls_total",
                   "n_carriers_cases", "n_carriers_controls"]],
               "fig6_forest_plot_data", cohort=cohort)
    print(f"[fig6_forest_plot] done — {cohort}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", default="cohortI", choices=["cohortI", "cohortII"])
    args = parser.parse_args()
    make(args.cohort)


if __name__ == "__main__":
    main()
