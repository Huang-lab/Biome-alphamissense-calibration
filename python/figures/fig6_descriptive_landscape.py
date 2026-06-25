"""Figure 6 — Descriptive carrier-frequency landscape (user-requested).

Five descriptive panels — no regression — that answer "where do the
carriers sit?" across phenotype, gene, and ancestry. Backed by Tables
5 / S10 / S11 / S12 (per the plan).

  6A  Per-phenotype carrier-frequency bars, sorted descending by
      AM-primary frequency. One bar group per cancer phenotype; three
      side-by-side bars per group for the three definitions.
  6B  Phenotype x gene carrier-frequency heatmap. Three sub-panels
      (one per carrier definition). Rows = phenotypes (cancer sites
      + control strata), cols = 28 genes ordered by syndrome group.
  6C  Ancestry-overall carrier-frequency bars; six ancestry groups
      with three series each. Shows AM-primary is not artefactually
      concentrated in one ancestry.
  6D  Ancestry x phenotype heatmap of AM-primary carrier frequency.
  6E  Ancestry x gene heatmap of AM-primary carrier frequency.

Reads ``C_regression_matrix.tsv`` (samples-by-syndrome) for 6A / 6C / 6D
and ``B_ACMG_vs_AM_comparison.tsv`` (carriers-by-variant with the gene
column) for 6B / 6E. Falls back to a placeholder when ``B_*.tsv`` is
absent — the script still emits 6A / 6C / 6D so the figure file is
valid mid-development.

Usage:
    python -m python.figures.fig6_descriptive_landscape \\
        --c-table results/cohortI/C_regression_matrix.tsv \\
        --b-table results/cohortI/B_ACMG_vs_AM_comparison.tsv \\
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


CATS = [("AMprimary", "AM_primary"),
        ("ACMG",      "ACMG_PLP"),
        ("AMonly",    "AM_only_non_PLP")]


def _any_carrier_per_sample(df, suffix: str) -> "pd.Series":
    """Per-sample bool: TRUE iff carrier_<any>_<suffix> is TRUE."""
    cols = [c for c in df.columns if c.startswith("carrier_") and c.endswith(f"_{suffix}")]
    if not cols:
        import pandas as pd
        return pd.Series([False] * len(df), index=df.index)
    return df[cols].apply(lambda row: any(v.upper() == "TRUE" for v in row), axis=1)


def panel_6A(ax, df) -> None:
    """Per-phenotype carrier-frequency bars, sorted by AM-primary freq desc."""
    import numpy as np
    import pandas as pd
    rows = []
    for suffix, label in CATS:
        carrier = _any_carrier_per_sample(df, suffix)
        grp = pd.DataFrame({"Group": df["Group"], "carrier": carrier})
        for phen, sub in grp.groupby("Group"):
            if not phen:
                continue
            rows.append({"phenotype": phen, "category": label,
                         "n_total": len(sub),
                         "n_carriers": int(sub["carrier"].sum()),
                         "freq": sub["carrier"].mean() if len(sub) else 0})
    if not rows:
        ax.text(0.5, 0.5, "no rows", ha="center", va="center",
                transform=ax.transAxes, color="grey")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title("A  Carrier freq per phenotype  (sorted by AM-primary)")
        return
    pf = pd.DataFrame(rows)
    ampr = pf[pf["category"] == "AM_primary"].set_index("phenotype")["freq"]
    order = ampr.sort_values(ascending=False).index.tolist()
    x = np.arange(len(order))
    w = 0.27
    for i, (_, label) in enumerate(CATS):
        h = [pf[(pf["phenotype"] == p) & (pf["category"] == label)]["freq"].iloc[0]
             if not pf[(pf["phenotype"] == p) & (pf["category"] == label)].empty else 0
             for p in order]
        ax.bar(x + (i - 1) * w, h, width=w,
               color=common.CATEGORY_COLORS[label],
               label=common.CATEGORY_LABELS[label])
    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=45, ha="right", fontsize=6)
    ax.set_ylabel("carrier frequency")
    ax.legend(fontsize=6, loc="upper right")
    ax.set_title("A  Carrier freq per phenotype  (sorted by AM-primary)")


def panel_6B(ax, b_df, c_df) -> None:
    """Phenotype x gene heatmap of AM-primary carrier frequency.

    Plotted only for the AM-primary category to keep the figure readable.
    The other two categories are computed but stored in companion TSVs
    (Tables S10/S12) — caller can subset further if needed.
    """
    import pandas as pd
    if b_df is None or "gene" not in (b_df.columns if b_df is not None else []):
        ax.text(0.5, 0.5, "B_ACMG_vs_AM_comparison.tsv not found",
                ha="center", va="center", transform=ax.transAxes, color="grey")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title("B  Phenotype x gene heatmap  (AM-primary)")
        return
    am_carriers = b_df[b_df["in_AM_primary"].str.lower() == "yes"]
    if am_carriers.empty:
        ax.text(0.5, 0.5, "no AM-primary carriers in B", ha="center", va="center",
                transform=ax.transAxes, color="grey")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title("B  Phenotype x gene heatmap  (AM-primary)")
        return
    pheno_by_sample = c_df.set_index("sample_id")["Group"].to_dict()
    am_carriers = am_carriers.assign(
        phenotype=am_carriers["sample_id"].map(pheno_by_sample),
    )
    grp = am_carriers.dropna(subset=["phenotype"]).groupby(["phenotype", "gene"]).size().rename("n").reset_index()
    pheno_total = c_df.groupby("Group").size().to_dict()
    grp["freq"] = grp.apply(lambda r: r["n"] / max(pheno_total.get(r["phenotype"], 1), 1), axis=1)
    mat = grp.pivot(index="phenotype", columns="gene", values="freq").fillna(0)
    gene_order = [g for g in common.GENE_ORDER if g in mat.columns]
    mat = mat.reindex(columns=gene_order)
    im = ax.imshow(mat.values, cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(mat.columns)))
    ax.set_xticklabels(mat.columns, rotation=90, fontsize=5)
    ax.set_yticks(range(len(mat.index)))
    ax.set_yticklabels(mat.index, fontsize=6)
    ax.figure.colorbar(im, ax=ax, shrink=0.7, label="freq")
    ax.set_title("B  Phenotype x gene  (AM-primary carrier freq)")


def panel_6C(ax, df) -> None:
    """Ancestry-overall carrier-frequency bars (three series)."""
    import numpy as np
    import pandas as pd
    rows = []
    ancestries = [a for a in common.ANCESTRY_ORDER if a in df["genetically_determined"].unique()]
    other = sorted(set(df["genetically_determined"].dropna().unique()) - set(ancestries))
    ancestries = ancestries + other
    for suffix, label in CATS:
        carrier = _any_carrier_per_sample(df, suffix)
        for a in ancestries:
            sub = df[df["genetically_determined"] == a]
            rows.append({"ancestry": a, "category": label,
                         "freq": carrier.loc[sub.index].mean() if len(sub) else 0})
    if not rows:
        ax.text(0.5, 0.5, "no ancestry rows", ha="center", va="center",
                transform=ax.transAxes, color="grey")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title("C  Carrier freq per ancestry  (any syndrome)")
        return
    pf = pd.DataFrame(rows)
    x = np.arange(len(ancestries))
    w = 0.27
    for i, (_, label) in enumerate(CATS):
        h = [pf[(pf["ancestry"] == a) & (pf["category"] == label)]["freq"].iloc[0]
             if not pf[(pf["ancestry"] == a) & (pf["category"] == label)].empty else 0
             for a in ancestries]
        ax.bar(x + (i - 1) * w, h, width=w,
               color=common.CATEGORY_COLORS[label], label=common.CATEGORY_LABELS[label])
    ax.set_xticks(x)
    ax.set_xticklabels(ancestries, fontsize=6)
    ax.set_ylabel("carrier frequency")
    ax.legend(fontsize=6, loc="upper right")
    ax.set_title("C  Carrier freq per ancestry  (any syndrome)")


def panel_6D(ax, df) -> None:
    """Ancestry x phenotype heatmap of AM-primary carrier frequency."""
    import pandas as pd
    carrier = _any_carrier_per_sample(df, "AMprimary")
    pf = pd.DataFrame({
        "ancestry": df["genetically_determined"],
        "phenotype": df["Group"],
        "carrier":  carrier.astype(int),
    })
    pf = pf[pf["ancestry"].astype(bool) & pf["phenotype"].astype(bool)]
    if pf.empty:
        ax.text(0.5, 0.5, "no rows", ha="center", va="center",
                transform=ax.transAxes, color="grey")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title("D  Ancestry x phenotype  (AM-primary freq)")
        return
    n = pf.groupby(["ancestry", "phenotype"]).size().rename("n")
    k = pf.groupby(["ancestry", "phenotype"])["carrier"].sum().rename("k")
    freq = (k / n).rename("freq").reset_index()
    mat = freq.pivot(index="ancestry", columns="phenotype", values="freq").fillna(0)
    ancestries = [a for a in common.ANCESTRY_ORDER if a in mat.index] + [a for a in mat.index if a not in common.ANCESTRY_ORDER]
    mat = mat.reindex(index=ancestries)
    im = ax.imshow(mat.values, cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(mat.columns)))
    ax.set_xticklabels(mat.columns, rotation=45, ha="right", fontsize=6)
    ax.set_yticks(range(len(mat.index)))
    ax.set_yticklabels(mat.index, fontsize=6)
    ax.figure.colorbar(im, ax=ax, shrink=0.7, label="AM-primary freq")
    ax.set_title("D  Ancestry x phenotype  (AM-primary freq)")


def panel_6E(ax, b_df, c_df) -> None:
    """Ancestry x gene heatmap of AM-primary carrier frequency."""
    import pandas as pd
    if b_df is None or "gene" not in b_df.columns:
        ax.text(0.5, 0.5, "B_ACMG_vs_AM_comparison.tsv not found",
                ha="center", va="center", transform=ax.transAxes, color="grey")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title("E  Ancestry x gene  (AM-primary freq)")
        return
    anc_by_sample = c_df.set_index("sample_id")["genetically_determined"].to_dict()
    n_by_anc = c_df.groupby("genetically_determined").size().to_dict()
    am = b_df[b_df["in_AM_primary"].str.lower() == "yes"].assign(
        ancestry=lambda d: d["sample_id"].map(anc_by_sample),
    )
    am = am.dropna(subset=["ancestry"])
    if am.empty:
        ax.text(0.5, 0.5, "no AM-primary carriers", ha="center", va="center",
                transform=ax.transAxes, color="grey")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title("E  Ancestry x gene  (AM-primary freq)")
        return
    k = am.groupby(["ancestry", "gene"]).size().rename("k").reset_index()
    k["freq"] = k.apply(lambda r: r["k"] / max(n_by_anc.get(r["ancestry"], 1), 1), axis=1)
    mat = k.pivot(index="ancestry", columns="gene", values="freq").fillna(0)
    ancestries = [a for a in common.ANCESTRY_ORDER if a in mat.index] + [a for a in mat.index if a not in common.ANCESTRY_ORDER]
    mat = mat.reindex(index=ancestries)
    gene_order = [g for g in common.GENE_ORDER if g in mat.columns]
    mat = mat.reindex(columns=gene_order)
    im = ax.imshow(mat.values, cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(mat.columns)))
    ax.set_xticklabels(mat.columns, rotation=90, fontsize=5)
    ax.set_yticks(range(len(mat.index)))
    ax.set_yticklabels(mat.index, fontsize=6)
    ax.figure.colorbar(im, ax=ax, shrink=0.7, label="AM-primary freq")
    ax.set_title("E  Ancestry x gene  (AM-primary freq)")


def make(out_dir: str, *, c_table: str, b_table: Optional[str] = None,
         cohort_name: Optional[str] = None) -> None:
    import matplotlib.pyplot as plt
    import pandas as pd
    common.apply_rcparams()
    c_df = pd.read_csv(c_table, sep="\t", dtype=str, na_filter=False)
    b_df = None
    if b_table and os.path.isfile(b_table):
        try:
            b_df = pd.read_csv(b_table, sep="\t", dtype=str, na_filter=False)
        except Exception as e:  # noqa: BLE001
            LOG.warning("fig6: cannot read B table %s: %s", b_table, e)

    fig = plt.figure(figsize=(common.FIG_WIDTH_DOUBLE[0],
                              common.FIG_WIDTH_DOUBLE[0] * 1.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.4, 1.4],
                          hspace=0.7, wspace=0.4)
    panel_6A(fig.add_subplot(gs[0, 0]), c_df)
    panel_6C(fig.add_subplot(gs[0, 1]), c_df)
    panel_6B(fig.add_subplot(gs[1, :]), b_df, c_df)
    panel_6D(fig.add_subplot(gs[2, 0]), c_df)
    panel_6E(fig.add_subplot(gs[2, 1]), b_df, c_df)
    basename = f"fig6.{cohort_name}" if cohort_name else "fig6"
    common.save_both(fig, Path(out_dir), basename)
    plt.close(fig)
    LOG.info("fig6: wrote %s.png + %s.pdf to %s", basename, basename, out_dir)


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--c-table", required=True)
    ap.add_argument("--b-table", default=None)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--cohort-name", default=None,
                    help="if set, suffix output filenames with .<cohort>")
    args = ap.parse_args(argv)
    make(args.out_dir, c_table=args.c_table, b_table=args.b_table,
         cohort_name=args.cohort_name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
