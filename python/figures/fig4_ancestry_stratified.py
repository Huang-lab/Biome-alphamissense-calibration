"""Figure 4 — Ancestry-stratified case/control + PCA scatter (Action 6).

Three panels:
  4A  Case/control x ancestry x variant-category (3D bar / heatmap):
      per-ancestry carrier-frequency in cases vs controls under each
      of the three definitions. Direct answer to R1's ancestry note.
  4B  Per-syndrome OR by ancestry: rerun the run_stats logistic
      regression *within* each ancestry group; forest plot 15 x 5
      (ancestries x syndromes) with min-cell-size suppression.
  4C  PC1 x PC2 scatter colored by AM-only carrier status; if AM-only
      carriers cluster in one ancestry the calibration is artefactual,
      if they spread the signal is clinical.

This is a SCAFFOLD — the panels render a placeholder until 4B's
per-ancestry regression is wired up (it needs run_stats.py to grow a
``--stratify-by`` flag, planned for a later commit).

Usage:
    python -m python.figures.fig4_ancestry_stratified \\
        --c-table results/cohortI/C_regression_matrix.tsv \\
        --out-dir results/figures
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import List, Optional

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from util import LOG  # noqa: E402

from . import common


def _placeholder(ax, msg: str) -> None:
    ax.set_facecolor("#f5f5f5")
    ax.text(0.5, 0.5, msg, ha="center", va="center", fontsize=8,
            transform=ax.transAxes, color="#333333", wrap=True)
    ax.set_xticks([]); ax.set_yticks([])


def panel_4A(ax, df) -> None:
    """Ancestry x case/control x variant-category heatmap of carrier freq."""
    import numpy as np
    import pandas as pd
    is_ctrl = df["Group"].str.lower().str.contains("control", na=False)
    df = df.assign(_case=~is_ctrl, _ancestry=df["genetically_determined"].fillna("OTHER"))
    rows: List[dict] = []
    for cat in ["AMprimary", "AMonly", "ACMG"]:
        for anc in [a for a in common.ANCESTRY_ORDER if a in df["_ancestry"].unique()] \
                   + sorted(set(df["_ancestry"].unique()) - set(common.ANCESTRY_ORDER)):
            for is_case in [True, False]:
                sub = df[(df["_ancestry"] == anc) & (df["_case"] == is_case)]
                if sub.empty:
                    continue
                # any syndrome carrier = carrier_any_{cat}
                cols = [c for c in df.columns if c.startswith("carrier_") and c.endswith(f"_{cat}")]
                if not cols:
                    continue
                n_carrier = int((sub[cols].apply(lambda r: any(v.upper() == "TRUE" for v in r), axis=1)).sum())
                rows.append({
                    "category": cat, "ancestry": anc,
                    "case": "case" if is_case else "control",
                    "freq": n_carrier / max(len(sub), 1),
                    "n":    len(sub),
                })
    if not rows:
        _placeholder(ax, "no ancestry/category rows to plot")
        ax.set_title("A  Ancestry x case/control x variant-category")
        return
    pf = pd.DataFrame(rows)
    pf["row"] = pf["ancestry"] + " (" + pf["case"] + ")"
    mat = pf.pivot(index="row", columns="category", values="freq")
    mat = mat.reindex(columns=["ACMG", "AMprimary", "AMonly"])
    im = ax.imshow(mat.values, cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(mat.columns)))
    ax.set_xticklabels([common.CATEGORY_LABELS.get(
        {"ACMG":"ACMG_PLP","AMprimary":"AM_primary","AMonly":"AM_only_non_PLP"}[c], c)
        for c in mat.columns], rotation=45, ha="right", fontsize=6)
    ax.set_yticks(range(len(mat.index)))
    ax.set_yticklabels(mat.index, fontsize=6)
    ax.figure.colorbar(im, ax=ax, shrink=0.7, label="carrier freq")
    ax.set_title("A  Ancestry x case/control x variant-category")


def panel_4B(ax) -> None:
    _placeholder(
        ax,
        "Fig 4B pending\nPer-syndrome OR by ancestry forest plot.\n"
        "Needs run_stats.py --stratify-by genetically_determined\n"
        "(planned next commit); will read\n"
        "results/<cohort>/stats_syndrome_associations_by_ancestry.tsv.",
    )
    ax.set_title("B  Per-syndrome OR by ancestry")


def panel_4C(ax, df) -> None:
    """PC1 x PC2 scatter colored by AM-only carrier status."""
    import pandas as pd
    pf = df.assign(
        PC1=pd.to_numeric(df["PC1"], errors="coerce"),
        PC2=pd.to_numeric(df["PC2"], errors="coerce"),
    )
    am_only_cols = [c for c in df.columns if c.startswith("carrier_") and c.endswith("_AMonly")]
    pf["_am_only_any"] = pf[am_only_cols].apply(
        lambda r: any(v.upper() == "TRUE" for v in r), axis=1) if am_only_cols else False
    pf = pf.dropna(subset=["PC1", "PC2"])
    if pf.empty:
        _placeholder(ax, "Fig 4C — no samples with valid PC1/PC2\n"
                          "(check the id_bridge / PC join)")
        ax.set_title("C  PC1 x PC2  (AM-only carriers highlighted)")
        return
    ctrl = pf[~pf["_am_only_any"]]
    carr = pf[ pf["_am_only_any"]]
    ax.scatter(ctrl["PC1"], ctrl["PC2"], s=2, alpha=0.4, color="#bbbbbb", label="non-carrier")
    ax.scatter(carr["PC1"], carr["PC2"], s=8, alpha=0.9,
               color=common.CATEGORY_COLORS["AM_only_non_PLP"], label="AM-only carrier")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.legend(fontsize=6, loc="best")
    ax.set_title("C  PC1 x PC2  (AM-only carriers highlighted)")


def make(out_dir: str, *, c_table: str, cohort_name: Optional[str] = None) -> None:
    import matplotlib.pyplot as plt
    import pandas as pd
    common.apply_rcparams()
    df = pd.read_csv(c_table, sep="\t", dtype=str, na_filter=False)
    fig, axes = plt.subplots(1, 3, figsize=(common.FIG_WIDTH_DOUBLE[0],
                                            common.FIG_WIDTH_DOUBLE[0] * 0.4))
    panel_4A(axes[0], df)
    panel_4B(axes[1])
    panel_4C(axes[2], df)
    fig.tight_layout()
    basename = f"fig4.{cohort_name}" if cohort_name else "fig4"
    common.save_both(fig, Path(out_dir), basename)
    plt.close(fig)
    LOG.info("fig4: wrote %s.png + %s.pdf to %s", basename, basename, out_dir)


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--c-table", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--cohort-name", default=None,
                    help="if set, suffix output filenames with .<cohort>")
    args = ap.parse_args(argv)
    make(args.out_dir, c_table=args.c_table, cohort_name=args.cohort_name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
