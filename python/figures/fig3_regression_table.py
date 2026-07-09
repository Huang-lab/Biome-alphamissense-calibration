"""
Fig 3 — Logistic regression results table.

Reads results/cohortI/regression_results.tsv (produced by run_regression.py)
and renders a publication-ready styled matplotlib table.

Color coding:
  q_BH < 0.01 → light orange (#FDDBC7)
  q_BH < 0.05 → light yellow (#FFF9C4)
  otherwise   → white

Usage:
    python -m python.figures.fig3_regression_table --cohort cohortI
    # Prerequisite: python python/run_regression.py --cohort cohortI
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
    RESULTS_I, RESULTS_II,
    COLORS, VC_ACMG, VC_CLINVAR, VC_AM, VC_AM_ONLY,
    SYNDROMES,
    FIGW_DOUBLE,
    save_fig, save_table, cohort_label,
)

VC_COLOR = {
    VC_ACMG:    COLORS[VC_ACMG],
    VC_CLINVAR: COLORS[VC_CLINVAR],
    VC_AM:      COLORS[VC_AM],
    VC_AM_ONLY: COLORS[VC_AM_ONLY],
}


def _row_bg(q):
    if pd.isna(q):
        return "#FFFFFF"
    if q < 0.01:
        return "#FDDBC7"
    if q < 0.05:
        return "#FFF9C4"
    return "#FFFFFF"


def make(cohort="cohortI"):
    out_dir = RESULTS_I if cohort == "cohortI" else RESULTS_II
    rr_path = os.path.join(out_dir, "regression_results.tsv")

    if not os.path.exists(rr_path):
        print(f"[fig3] ERROR: {rr_path} not found. Run run_regression.py first.")
        return

    df = pd.read_csv(rr_path, sep="\t")
    print(f"[fig3_regression_table] {len(df)} regression rows for {cohort}")

    # Filter to meaningful results and sort; remove MUTYH (may be in older output)
    df = df[df["status"] == "ok"].copy()
    df = df[df["syndrome_key"] != "MUTYH_Associated_Polyposis"].copy()
    df = df.sort_values(["syndrome", "variant_class", "phenotype"]).reset_index(drop=True)

    # Format columns for display
    def fmt_or(row):
        if pd.isna(row["OR"]):
            return "—"
        return f"{row['OR']:.2f} ({row['CI_low']:.2f}–{row['CI_high']:.2f})"

    def fmt_p(val):
        if pd.isna(val):
            return "—"
        if val < 0.001:
            return "<0.001"
        return f"{val:.3f}"

    df["OR_fmt"] = df.apply(fmt_or, axis=1)
    df["p_fmt"]  = df["p_value"].apply(fmt_p)
    df["q_fmt"]  = df["q_BH"].apply(fmt_p)

    # Build genes column from syndrome_key
    df["genes"] = df["syndrome_key"].map(
        lambda k: ", ".join(SYNDROMES.get(k, {}).get("genes", []))
    )

    display_cols = ["syndrome", "genes", "phenotype", "variant_class",
                    "n_cases_total", "n_controls_total",
                    "n_carriers_cases", "n_carriers_controls",
                    "OR_fmt", "p_fmt", "q_fmt"]
    col_headers  = ["Syndrome", "Genes", "Phenotype", "Variant class",
                    "N cases", "N ctrl",
                    "Carriers\n(cases)", "Carriers\n(ctrl)",
                    "OR (95% CI)", "p", "q_BH"]

    tbl = df[display_cols].values.tolist()
    n_rows = len(tbl)

    if n_rows == 0:
        print("[fig3] No valid rows to display.")
        return

    # -----------------------------------------------------------------------
    # Render as matplotlib table (paginated — 40 rows per page)
    # -----------------------------------------------------------------------
    PAGE_SIZE = 30
    n_pages   = max(1, (n_rows + PAGE_SIZE - 1) // PAGE_SIZE)

    for page in range(n_pages):
        page_rows  = tbl[page * PAGE_SIZE : (page + 1) * PAGE_SIZE]
        page_q     = df["q_BH"].iloc[page * PAGE_SIZE : (page + 1) * PAGE_SIZE].values
        page_vc    = df["variant_class"].iloc[page * PAGE_SIZE : (page + 1) * PAGE_SIZE].values
        n_pr = len(page_rows)

        row_h  = 0.28
        col_h  = 0.55
        fig_h  = n_pr * row_h + col_h + 1.0
        fig_w  = FIGW_DOUBLE + 3.5

        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        ax.axis("off")

        all_rows = [col_headers] + page_rows
        cell_colors = [["#E8E8E8"] * len(col_headers)]  # header

        for i, (q, vc) in enumerate(zip(page_q, page_vc)):
            bg = _row_bg(q)
            row_colors = [bg] * len(col_headers)
            cell_colors.append(row_colors)

        tbl_artist = ax.table(
            cellText=all_rows,
            cellLoc="left",
            loc="upper left",
            cellColours=cell_colors,
        )
        tbl_artist.auto_set_font_size(False)
        tbl_artist.set_fontsize(10)
        tbl_artist.auto_set_column_width(list(range(len(col_headers))))

        # Bold header
        for j in range(len(col_headers)):
            tbl_artist[0, j].set_text_props(fontweight="bold", fontsize=11)

        # Colour variant-class text (col index 3: syndrome, genes, phenotype, variant_class)
        VC_COL_IDX = 3
        for i in range(1, n_pr + 1):
            vc = page_rows[i - 1][VC_COL_IDX]
            color = VC_COLOR.get(vc, "#000000")
            tbl_artist[i, VC_COL_IDX].set_text_props(color=color, fontweight="bold")

        suffix = f"_p{page + 1}" if n_pages > 1 else ""
        title  = (f"Logistic regression results — {cohort_label(cohort)}"
                  + (f" (page {page + 1}/{n_pages})" if n_pages > 1 else ""))
        ax.set_title(title, fontsize=14, fontweight="bold", pad=8, loc="left")

        legend_text = ("Row shading: light orange = q<0.01, light yellow = q<0.05")
        fig.text(0.01, 0.01, legend_text, fontsize=9, color="#555555")

        save_fig(fig, f"fig3_regression_table{suffix}", cohort=cohort)

    # Save data table (TSV with all numeric columns + genes column)
    save_table(df[["syndrome", "genes", "syndrome_key", "phenotype", "variant_class",
                   "n_cases_total", "n_controls_total",
                   "n_carriers_cases", "n_carriers_controls",
                   "OR", "CI_low", "CI_high", "p_value", "q_BH"]],
               "fig3_regression_data", cohort=cohort)
    print(f"[fig3_regression_table] done — {cohort}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", default="cohortI", choices=["cohortI", "cohortII"])
    args = parser.parse_args()
    make(args.cohort)


if __name__ == "__main__":
    main()
