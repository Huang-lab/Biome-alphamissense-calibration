"""Figure 3 — per-phenotype forest plots of significant case/control associations.

Reads ``results/<cohort>/stats_syndrome_associations.significant_q0.10.tsv``
(the significant-only sibling produced by ``python/run_stats.py``) and writes
ONE file per cancer phenotype that has at least one significant row:

    fig_logreg_forest_phenotype_<phenotype>.{png,pdf}

Each plotted row is one (syndrome, variant_category) pair with q_BH below the
threshold (default 0.10, configurable via ``--q-threshold``), showing OR with
95% CI on a log axis. Canonical syndrome × phenotype pairs are marked with an
asterisk in the y-label. Bars colored by variant_category using the shared
palette from ``common.CATEGORY_COLORS``. Phenotypes with zero significant
rows are skipped silently (logged at INFO).

If the explicit ``--significant-stats`` path is missing or empty, the script
falls back to reading the full ``--stats`` TSV and applying the q threshold
in-memory — that way callers don't need both files on disk.

Usage:
    python -m python.figures.fig3_case_control_associations \\
        --stats           results/cohortI/stats_syndrome_associations.tsv \\
        --significant     results/cohortI/stats_syndrome_associations.significant_q0.10.tsv \\
        --q-threshold     0.10 \\
        --out-dir         results/figures \\
        --cohort-name     cohortI
"""
from __future__ import annotations

import argparse
import math
import os
import re
import sys
from pathlib import Path
from typing import List, Optional

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from util import LOG  # noqa: E402

from . import common

VARIANT_CATEGORY_ORDER = [
    "ACMG_PLP",
    "AM_calibrated",
    "AM_calibrated_not_PLP",
    "AM_global_0.864",
]


def _safe_phenotype_slug(phenotype: str) -> str:
    """Lowercase, alnum-or-underscore for filenames. Keeps the original
    phenotype string usable as a figure title; the slug is filesystem-safe."""
    s = re.sub(r"[^A-Za-z0-9]+", "_", phenotype).strip("_").lower()
    return s or "unknown"


def _load_significant(*, significant_path: Optional[str],
                      stats_path: Optional[str],
                      q_threshold: float) -> "pd.DataFrame":
    """Return the significant-row subset as a DataFrame. Prefers the
    pre-filtered TSV; falls back to filtering the full stats TSV at
    ``q_threshold`` if the sibling isn't on disk yet."""
    import pandas as pd
    if significant_path and os.path.isfile(significant_path):
        df = pd.read_csv(significant_path, sep="\t")
        LOG.info("fig3: loaded %d significant rows from %s", len(df), significant_path)
        return df
    if not stats_path or not os.path.isfile(stats_path):
        raise FileNotFoundError(
            f"Neither significant-stats path ({significant_path!r}) nor "
            f"full stats path ({stats_path!r}) exists."
        )
    full = pd.read_csv(stats_path, sep="\t")
    sig = full[full["q_BH"].notna() & (full["q_BH"] < q_threshold)].copy()
    sig = sig.sort_values(["q_BH", "p_value"], ascending=[True, True])
    LOG.info("fig3: filtered %d/%d rows from %s with q_BH < %.3f",
             len(sig), len(full), stats_path, q_threshold)
    return sig


def _plot_forest_for_phenotype(df_p: "pd.DataFrame", phenotype: str,
                               *, q_threshold: float, cohort_name: Optional[str]):
    """Render one forest plot. Caller already filtered to a single phenotype.
    Returns the matplotlib figure (caller saves and closes it)."""
    import matplotlib.pyplot as plt
    import pandas as pd
    common.apply_rcparams()

    # Stable sort: significant rows already sorted by q_BH ascending in the
    # input; rank by OR descending for a visually intuitive forest plot.
    df_p = df_p.assign(
        _OR_for_sort=df_p["OR"].astype(float).fillna(0.0)
    ).sort_values("_OR_for_sort", ascending=True)  # ascending → largest at top

    n_rows = len(df_p)
    height = max(1.4, 0.28 * n_rows + 0.8)
    fig, ax = plt.subplots(figsize=(common.FIG_WIDTH_DOUBLE[0], height))

    y = list(range(n_rows))
    for yi, (_, r) in zip(y, df_p.iterrows()):
        cat = r["variant_category"]
        color = common.CATEGORY_COLORS.get(cat, "#444444")
        or_val = float(r["OR"]) if pd.notna(r["OR"]) else float("nan")
        lo = float(r["OR_95CI_lo"]) if pd.notna(r["OR_95CI_lo"]) else float("nan")
        hi = float(r["OR_95CI_hi"]) if pd.notna(r["OR_95CI_hi"]) else float("nan")
        if math.isnan(or_val):
            ax.text(1.0, yi, "  OR n/a", va="center", fontsize=6, color="grey")
            continue
        xerr_lo = max(or_val - lo, 0.0) if not math.isnan(lo) else 0.0
        xerr_hi = max(hi - or_val, 0.0) if not math.isnan(hi) else 0.0
        ax.errorbar(or_val, yi, xerr=[[xerr_lo], [xerr_hi]],
                    fmt="o", markersize=4,
                    color=color, ecolor=color, lw=1.0,
                    capsize=2.0)

    ax.axvline(1.0, color="grey", ls="--", lw=0.5)
    ax.set_yticks(y)
    ylabels = []
    for _, r in df_p.iterrows():
        syn = common.SYNDROME_SHORT_LABEL.get(r["syndrome"], r["syndrome"])
        cat_lbl = common.CATEGORY_LABELS.get(r["variant_category"], r["variant_category"])
        mark = "*" if bool(r["is_canonical"]) else " "
        q = float(r["q_BH"]) if pd.notna(r["q_BH"]) else float("nan")
        ylabels.append(f"{mark}{syn} · {cat_lbl}  (q={q:.2g})")
    ax.set_yticklabels(ylabels, fontsize=6)
    ax.set_xscale("log")
    ax.set_xlabel("OR (95% CI), log scale")
    title_bits = [f"phenotype: {phenotype}"]
    if cohort_name:
        title_bits.append(f"cohort: {cohort_name}")
    title_bits.append(f"q_BH < {q_threshold:g}, * = canonical syndrome-cancer pair")
    ax.set_title("  |  ".join(title_bits), fontsize=8)

    # Color-keyed legend (only the categories that actually appear)
    from matplotlib.patches import Patch
    cats_present = [c for c in VARIANT_CATEGORY_ORDER
                    if c in set(df_p["variant_category"])]
    handles = [Patch(facecolor=common.CATEGORY_COLORS[c],
                     label=common.CATEGORY_LABELS.get(c, c))
               for c in cats_present]
    if handles:
        ax.legend(handles=handles, loc="lower right", fontsize=6, frameon=False)

    fig.tight_layout()
    return fig


def make(out_dir: str, *, stats: Optional[str] = None,
         significant: Optional[str] = None,
         q_threshold: float = 0.10,
         cohort_name: Optional[str] = None) -> None:
    """Write per-phenotype forest plots into ``out_dir``. Skips phenotypes
    that have no rows below ``q_threshold`` (logs the skip at INFO)."""
    import matplotlib.pyplot as plt

    sig = _load_significant(
        significant_path=significant, stats_path=stats, q_threshold=q_threshold,
    )
    if sig.empty:
        LOG.warning("fig3: 0 significant rows total — no per-phenotype files written")
        return

    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    n_written = 0
    n_skipped = 0
    for phenotype, df_p in sig.groupby("phenotype"):
        if df_p.empty:
            n_skipped += 1
            LOG.info("fig3: skipping %r — 0 significant rows", phenotype)
            continue
        fig = _plot_forest_for_phenotype(
            df_p, phenotype,
            q_threshold=q_threshold, cohort_name=cohort_name,
        )
        slug = _safe_phenotype_slug(phenotype)
        base = f"fig_logreg_forest_phenotype_{slug}"
        if cohort_name:
            base = f"{base}.{cohort_name}"
        common.save_both(fig, out_path, base)
        plt.close(fig)
        n_written += 1
        LOG.info("fig3: wrote %s.png + .pdf (%d rows)", base, len(df_p))

    LOG.info("fig3: wrote %d per-phenotype figure(s); skipped %d phenotype(s) with 0 significant rows",
             n_written, n_skipped)


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stats", default=None,
                    help="results/<cohort>/stats_syndrome_associations.tsv "
                         "(fallback when --significant is missing)")
    ap.add_argument("--significant", default=None,
                    help="results/<cohort>/stats_syndrome_associations.significant_q0.10.tsv "
                         "— preferred input. Falls back to filtering --stats if absent.")
    ap.add_argument("--q-threshold", type=float, default=0.10,
                    help="q_BH cutoff used for the fallback filter. Has no effect "
                         "when --significant is provided AND that file exists.")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--cohort-name", default=None,
                    help="if set, suffix output filenames with .<cohort>")
    args = ap.parse_args(argv)
    if not args.stats and not args.significant:
        ap.error("at least one of --stats / --significant must be provided")
    make(args.out_dir, stats=args.stats, significant=args.significant,
         q_threshold=args.q_threshold, cohort_name=args.cohort_name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
