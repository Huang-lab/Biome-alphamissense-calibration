"""Figure 3 — Case-control test: expected associations replicate before novel.

This is the highest-stakes figure of the rebuttal set (Action 4 is flagged
"highest-risk; run first" in the rebuttal letter). All panels read
``results/<cohort>/stats_syndrome_associations.tsv`` produced by
``python/run_stats.py``.

  3A  OR heatmap, 15 syndromes x N cancer phenotypes. Three sub-panels
      {ACMG/AMP P/LP, AM-primary, AM-only-non-P/LP}. Canonical cells
      (``is_canonical=True``) outlined in black. Color = log(OR), grey
      for cells with n_carriers < 5 (NaN OR).

  3B  Post-hoc power for the AM-only-non-P/LP cells: observed OR vs the
      minimum detectable OR at 80% power. Tells R1 whether a missing
      association is power-limited or a true null.

  3C  Forest plot for MUTYH-Polyposis + MEN + Lynch vs Colorectal/
      Thyroid/Endometrial (the three Action 11 cells most likely to
      reach significance in BioMe).

  3D  Top non-canonical hits: filter to is_canonical=False & q_BH<0.1,
      rank by |log OR|, show top 5-10 as a side forest plot. Framed as
      hypothesis-generating.

Usage:
    python -m python.figures.fig3_case_control_associations \\
        --stats results/cohortI/stats_syndrome_associations.tsv \\
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

VARIANT_CATEGORY_ORDER = ["ACMG_PLP", "AM_primary", "AM_only_non_PLP"]


def _placeholder_panel(ax, msg: str) -> None:
    ax.set_facecolor("#f5f5f5")
    ax.text(0.5, 0.5, msg, ha="center", va="center", fontsize=8,
            transform=ax.transAxes, color="#333333", wrap=True)
    ax.set_xticks([]); ax.set_yticks([])


def panel_3A(axes, df: "pd.DataFrame") -> None:
    """Three-panel heatmap: rows = 15 syndromes, cols = cancer phenotypes,
    cells outlined where is_canonical."""
    import numpy as np
    import pandas as pd
    from matplotlib.patches import Rectangle

    phenotypes = sorted(df["phenotype"].unique().tolist())
    syndromes = [s for s in common.SYNDROME_ORDER if s in df["syndrome"].unique()]
    if not phenotypes or not syndromes:
        for ax in axes:
            _placeholder_panel(ax, "no rows in stats table")
        return

    # symmetric log color scale; OR=1 -> 0
    log_ors = df["OR"].apply(lambda x: math.log(x) if (x and x > 0) else float("nan"))
    vmax = max(abs(log_ors.min(skipna=True) or 0), abs(log_ors.max(skipna=True) or 0), 1.0)

    for ax, cat in zip(axes, VARIANT_CATEGORY_ORDER):
        sub = df[df["variant_category"] == cat]
        if sub.empty:
            _placeholder_panel(ax, f"{cat}\n(no rows)")
            continue
        mat = sub.pivot(index="syndrome", columns="phenotype", values="OR")
        mat = mat.reindex(index=syndromes, columns=phenotypes)
        canonical = sub.pivot(index="syndrome", columns="phenotype",
                              values="is_canonical").reindex(
            index=syndromes, columns=phenotypes,
        ).fillna(False).astype(bool)

        # log-OR for color. .map replaces deprecated .applymap in pandas 3.0+.
        log_mat = mat.map(lambda x: math.log(x) if (pd.notna(x) and x > 0) else float("nan"))
        im = ax.imshow(log_mat.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                       aspect="auto")
        ax.set_xticks(np.arange(len(phenotypes)))
        ax.set_xticklabels(phenotypes, rotation=45, ha="right", fontsize=6)
        ax.set_yticks(np.arange(len(syndromes)))
        ax.set_yticklabels([common.SYNDROME_SHORT_LABEL.get(s, s) for s in syndromes],
                           fontsize=6)
        ax.set_title(common.CATEGORY_LABELS[cat], fontsize=8)
        # Annotate canonical cells with a thicker outline
        for i, syn in enumerate(syndromes):
            for j, phen in enumerate(phenotypes):
                if canonical.iat[i, j]:
                    ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1,
                                           edgecolor="black", facecolor="none", lw=1.0))
                # grey out NaN
                if pd.isna(log_mat.iat[i, j]):
                    ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1,
                                           edgecolor="none", facecolor="#e0e0e0", alpha=0.7))
        # Only first panel keeps y labels for compactness
        if cat != VARIANT_CATEGORY_ORDER[0]:
            ax.set_yticks([])
    # Shared colorbar on the rightmost panel
    cb = axes[-1].figure.colorbar(im, ax=axes[-1], shrink=0.7, pad=0.02)
    cb.set_label("log(OR)", fontsize=7)


def panel_3B(ax, df: "pd.DataFrame") -> None:
    """Observed OR vs min-detectable-OR for AM-only-non-P/LP canonical
    cells. y-axis: observed OR (with NaN as 'not testable'). x-axis: min
    detectable OR. Identity line = 'just-detectable'."""
    import numpy as np
    sub = df[(df["variant_category"] == "AM_only_non_PLP") & df["is_canonical"]]
    if sub.empty:
        _placeholder_panel(ax, "no canonical AM-only cells")
        return
    mdo = sub["min_detectable_OR_80pct_power"].astype(float)
    obs = sub["OR"].astype(float)
    mask = mdo.notna() & obs.notna()
    if not mask.any():
        _placeholder_panel(ax, "no canonical cells with sufficient cells\nto compute MDO + observed OR")
        return
    ax.scatter(mdo[mask], obs[mask], s=14,
               color=common.CATEGORY_COLORS["AM_only_non_PLP"], alpha=0.7,
               edgecolor="black", linewidth=0.3)
    lo, hi = 0.5, max(2.0, float(np.nanmax(np.concatenate([mdo.dropna(), obs.dropna()]))) + 0.5)
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.5)
    ax.axhline(1.0, color="grey", lw=0.4, ls=":")
    ax.set_xlabel("min detectable OR (80% power, α=0.05)")
    ax.set_ylabel("observed OR")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_title("B  AM-only power vs effect (canonical cells)")


def panel_3C(ax, df: "pd.DataFrame") -> None:
    """Forest plot for the three Action-11 syndromes against their canonical
    cancers, three series per row (ACMG / AM-primary / AM-only)."""
    target = [
        ("MUTYH_Associated_Polyposis", "colorectal"),
        ("Multiple_Endocrine_Neoplasia", "thyroid"),
        ("Lynch_Syndrome",               "colorectal"),
    ]
    rows = []
    for syndrome, phen_token in target:
        sub_syn = df[df["syndrome"] == syndrome]
        if sub_syn.empty:
            continue
        # any phenotype whose name (case-insensitive) contains the token
        sub = sub_syn[sub_syn["phenotype"].str.lower().str.contains(phen_token, na=False)]
        if sub.empty:
            continue
        for cat in VARIANT_CATEGORY_ORDER:
            r = sub[sub["variant_category"] == cat]
            if r.empty:
                continue
            rec = r.iloc[0]
            rows.append({
                "label": f"{common.SYNDROME_SHORT_LABEL.get(syndrome, syndrome)} -> {phen_token}",
                "cat":   cat,
                "OR":    rec["OR"],
                "lo":    rec["OR_95CI_lo"],
                "hi":    rec["OR_95CI_hi"],
            })
    if not rows:
        _placeholder_panel(ax, "no rows for MUTYH-AP / MEN / Lynch panels")
        return
    import pandas as pd
    pf = pd.DataFrame(rows)
    yspace = list(range(len(pf)))[::-1]
    for y, (_, r) in zip(yspace, pf.iterrows()):
        if pd.isna(r["OR"]):
            ax.text(1.0, y, "  n<5", va="center", fontsize=6, color="grey")
            continue
        ax.errorbar(r["OR"], y,
                    xerr=[[r["OR"] - r["lo"]], [r["hi"] - r["OR"]]],
                    fmt="o", markersize=3,
                    color=common.CATEGORY_COLORS[r["cat"]],
                    ecolor=common.CATEGORY_COLORS[r["cat"]], lw=0.8)
    ax.axvline(1.0, color="grey", ls="--", lw=0.5)
    ax.set_yticks(yspace)
    ax.set_yticklabels(pf["label"], fontsize=6)
    ax.set_xscale("log")
    ax.set_xlabel("OR (95% CI)")
    ax.set_title("C  MUTYH-AP / MEN / Lynch enrichment (Action 11)")


def panel_3D(ax, df: "pd.DataFrame", top_n: int = 8, q_cutoff: float = 0.1) -> None:
    """Top non-canonical hits with q_BH below the cutoff."""
    import pandas as pd
    sub = df[~df["is_canonical"] & df["q_BH"].notna() & (df["q_BH"] < q_cutoff)]
    if sub.empty:
        _placeholder_panel(
            ax,
            f"no non-canonical cells with q_BH < {q_cutoff}\n"
            f"(novel-association story not supported by data — Action 5\n"
            f"would need to soften from claim to hypothesis-generating)",
        )
        ax.set_title("D  Novel hits  (q < {:.2f})".format(q_cutoff))
        return
    sub = sub.assign(abs_logor=sub["OR"].apply(
        lambda x: abs(math.log(x)) if (pd.notna(x) and x > 0) else 0
    )).sort_values("abs_logor", ascending=False).head(top_n)
    yspace = list(range(len(sub)))[::-1]
    for y, (_, r) in zip(yspace, sub.iterrows()):
        ax.errorbar(r["OR"], y,
                    xerr=[[r["OR"] - r["OR_95CI_lo"]], [r["OR_95CI_hi"] - r["OR"]]],
                    fmt="o", markersize=3,
                    color=common.CATEGORY_COLORS[r["variant_category"]],
                    ecolor=common.CATEGORY_COLORS[r["variant_category"]], lw=0.8)
    ax.axvline(1.0, color="grey", ls="--", lw=0.5)
    ax.set_yticks(yspace)
    ax.set_yticklabels(
        [f"{common.SYNDROME_SHORT_LABEL.get(s, s)} -> {p}  ({c})"
         for s, p, c in zip(sub["syndrome"], sub["phenotype"], sub["variant_category"])],
        fontsize=6,
    )
    ax.set_xscale("log")
    ax.set_xlabel("OR (95% CI)")
    ax.set_title(f"D  Novel hits  (q_BH < {q_cutoff})")


def make(out_dir: str, *, stats: str) -> None:
    import matplotlib.pyplot as plt
    import pandas as pd

    common.apply_rcparams()
    df = pd.read_csv(stats, sep="\t")
    fig = plt.figure(figsize=(common.FIG_WIDTH_DOUBLE[0],
                              common.FIG_WIDTH_DOUBLE[0] * 1.1))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.6, 1.0, 1.0],
                          hspace=0.6, wspace=0.5)
    axesA = [fig.add_subplot(gs[0, j]) for j in range(3)]
    axB   = fig.add_subplot(gs[1, 0])
    axC   = fig.add_subplot(gs[1, 1:])
    axD   = fig.add_subplot(gs[2, :])
    panel_3A(axesA, df)
    panel_3B(axB,   df)
    panel_3C(axC,   df)
    panel_3D(axD,   df)
    common.save_both(fig, Path(out_dir), "fig3")
    plt.close(fig)
    LOG.info("fig3: wrote fig3.png + fig3.pdf to %s", out_dir)


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stats", required=True,
                    help="results/<cohort>/stats_syndrome_associations.tsv")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args(argv)
    make(args.out_dir, stats=args.stats)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
