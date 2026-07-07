"""Figure 5 — Cross-cohort replication: Regeneron (Cohort I) vs Sema4 (II).

Two panels comparing the two cohorts side-by-side on the same syndrome-
level OR matrix:

  5A  Cohort I vs Cohort II OR for canonical syndrome-phenotype cells:
      scatter with identity line. Replication = points on or near y=x.
  5B  Meta-analyzed OR (random-effects across the two cohorts) per
      canonical cell; forest plot. Builds the "replicates across
      independently sequenced cohorts" claim (Action 8).

Reads both ``results/<cohort>/stats_syndrome_associations.tsv``. Falls
back to a placeholder when one cohort's stats table is missing — useful
during incremental development.

Usage:
    python -m python.figures.fig5_cohort_replication \\
        --stats-cohortI  results/cohortI/stats_syndrome_associations.tsv \\
        --stats-cohortII results/cohortII/stats_syndrome_associations.tsv \\
        --out-dir results/figures
"""
from __future__ import annotations

import argparse
import math
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


def panel_5A(ax, c1, c2) -> None:
    """Scatter of OR(cohort I) vs OR(cohort II) for canonical cells."""
    import pandas as pd
    merged = c1.merge(c2, on=["syndrome", "phenotype", "variant_category"],
                      suffixes=("_I", "_II"))
    canon = merged[merged["is_canonical_I"].fillna(False)]
    canon = canon.dropna(subset=["OR_I", "OR_II"])
    if canon.empty:
        _placeholder(ax, "no canonical cells fittable in both cohorts")
        ax.set_title("A  Cohort I vs Cohort II OR  (canonical cells)")
        return
    for cat, sub in canon.groupby("variant_category"):
        ax.scatter(sub["OR_I"], sub["OR_II"], s=12, alpha=0.8,
                   color=common.CATEGORY_COLORS.get(cat, "#888"),
                   label=common.CATEGORY_LABELS.get(cat, cat),
                   edgecolor="black", linewidth=0.3)
    ax.axline((1, 1), slope=1, color="grey", lw=0.5, ls="--")
    ax.axhline(1.0, color="grey", lw=0.3, ls=":")
    ax.axvline(1.0, color="grey", lw=0.3, ls=":")
    ax.set_xlabel("OR (Cohort I)"); ax.set_ylabel("OR (Cohort II)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.legend(fontsize=6, loc="lower right")
    ax.set_title("A  Cohort I vs Cohort II OR  (canonical cells)")


def panel_5B(ax, c1, c2) -> None:
    """Per-canonical-cell meta-analyzed OR (random effects, DerSimonian-Laird).
    Implemented as a simple inverse-variance weighted average on log-OR
    with the Q-statistic-derived tau^2."""
    import math

    import numpy as np
    import pandas as pd

    merged = c1.merge(c2, on=["syndrome", "phenotype", "variant_category"],
                      suffixes=("_I", "_II"))
    canon = merged[merged["is_canonical_I"].fillna(False)]
    if canon.empty:
        _placeholder(ax, "no canonical cells available")
        ax.set_title("B  Random-effects meta OR per canonical cell")
        return
    rows = []
    for _, r in canon.iterrows():
        if any(pd.isna(r[c]) for c in ["OR_I","OR_95CI_lo_I","OR_95CI_hi_I",
                                       "OR_II","OR_95CI_lo_II","OR_95CI_hi_II"]):
            continue
        b1 = math.log(r["OR_I"]); se1 = (math.log(r["OR_95CI_hi_I"]) - math.log(r["OR_95CI_lo_I"])) / 3.92
        b2 = math.log(r["OR_II"]); se2 = (math.log(r["OR_95CI_hi_II"]) - math.log(r["OR_95CI_lo_II"])) / 3.92
        # fixed-effect first
        w1 = 1.0 / se1**2; w2 = 1.0 / se2**2
        b_fe = (b1*w1 + b2*w2) / (w1 + w2)
        Q = w1 * (b1 - b_fe)**2 + w2 * (b2 - b_fe)**2
        tau2 = max(0.0, (Q - 1) / (w1 + w2 - (w1*w1 + w2*w2) / (w1 + w2))) if (w1 + w2) else 0.0
        w1r = 1.0 / (se1**2 + tau2); w2r = 1.0 / (se2**2 + tau2)
        b_re = (b1*w1r + b2*w2r) / (w1r + w2r)
        se_re = math.sqrt(1.0 / (w1r + w2r))
        rows.append({
            "label": f"{common.SYNDROME_SHORT_LABEL.get(r['syndrome'], r['syndrome'])} -> {r['phenotype']}  ({r['variant_category']})",
            "cat":   r["variant_category"],
            "OR":    math.exp(b_re),
            "lo":    math.exp(b_re - 1.96 * se_re),
            "hi":    math.exp(b_re + 1.96 * se_re),
        })
    if not rows:
        _placeholder(ax, "no canonical cells fittable in both cohorts")
        ax.set_title("B  Random-effects meta OR per canonical cell")
        return
    pf = pd.DataFrame(rows).sort_values("OR", ascending=False).head(12)
    y = list(range(len(pf)))[::-1]
    for yi, (_, r) in zip(y, pf.iterrows()):
        ax.errorbar(r["OR"], yi,
                    xerr=[[r["OR"] - r["lo"]], [r["hi"] - r["OR"]]],
                    fmt="o", markersize=3,
                    color=common.CATEGORY_COLORS.get(r["cat"], "#888"),
                    ecolor=common.CATEGORY_COLORS.get(r["cat"], "#888"), lw=0.8)
    ax.axvline(1.0, color="grey", ls="--", lw=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(pf["label"], fontsize=6)
    ax.set_xscale("log")
    ax.set_xlabel("OR (95% CI, random-effects meta)")
    ax.set_title("B  Meta-analyzed canonical cells")


def make(out_dir: str, *, stats_cohortI: str, stats_cohortII: str) -> None:
    import matplotlib.pyplot as plt
    import pandas as pd

    common.apply_rcparams()
    fig, axes = plt.subplots(1, 2, figsize=(common.FIG_WIDTH_DOUBLE[0],
                                            common.FIG_WIDTH_DOUBLE[0] * 0.5))
    if not (os.path.isfile(stats_cohortI) and os.path.isfile(stats_cohortII)):
        for ax in axes:
            _placeholder(ax, "need stats from BOTH cohorts\n(rerun 05_run_stats on each)")
        common.save_both(fig, Path(out_dir), "fig5")
        plt.close(fig)
        LOG.warning("fig5: missing one stats table; emitted placeholder")
        return
    c1 = pd.read_csv(stats_cohortI, sep="\t")
    c2 = pd.read_csv(stats_cohortII, sep="\t")
    panel_5A(axes[0], c1, c2)
    panel_5B(axes[1], c1, c2)
    fig.tight_layout()
    common.save_both(fig, Path(out_dir), "fig5")
    plt.close(fig)
    LOG.info("fig5: wrote fig5.png + fig5.pdf to %s", out_dir)


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stats-cohortI", required=True)
    ap.add_argument("--stats-cohortII", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args(argv)
    make(args.out_dir, stats_cohortI=args.stats_cohortI,
         stats_cohortII=args.stats_cohortII)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
