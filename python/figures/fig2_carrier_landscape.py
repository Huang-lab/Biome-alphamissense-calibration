"""Figure 2 — Carrier landscape: ACMG/AMP vs AM, calibrated honestly.

Three panels read straight from ``results/<cohort>/C_regression_matrix.tsv``
(samples-by-syndrome carrier flags), so no regression is needed for this
figure — it is the descriptive backbone for Actions 1, 2, 3 setup.

  2A  Per-syndrome carrier counts under three definitions
      {ACMG/AMP P/LP, AM-primary, AM-only-non-P/LP}.
      15 syndromes x 3 series side-by-side; the main descriptive table.

  2B  Stacked-bar share of {ACMG-only, AM-only-non-P/LP, both} per
      syndrome, normalized within syndrome. Shows AM-only is a
      meaningful but bounded slice of each syndrome's carriers.

  2C  Counterfactual: per-syndrome carrier counts under the old global
      0.864 threshold vs under the new Chen gene-specific calibration.
      Reads carrier_<group>_AM0864 (added by the pipeline) alongside
      carrier_<group>_AMprimary; lights up which genes/syndromes the
      old threshold over- vs under-claimed.

Usage:
    python -m python.figures.fig2_carrier_landscape \\
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


def _carrier_counts(df, suffix: str) -> "pd.Series":
    """Per-syndrome carrier counts for a given column suffix (ACMG /
    AMprimary / AMonly / AM0864)."""
    import pandas as pd
    out = {}
    for col in df.columns:
        if not col.startswith("carrier_") or not col.endswith(f"_{suffix}"):
            continue
        syndrome = col[len("carrier_"):-len(f"_{suffix}")]
        out[syndrome] = int((df[col].str.upper() == "TRUE").sum())
    return pd.Series(out)


def panel_2A(ax, df) -> None:
    """Per-syndrome counts under three definitions."""
    import numpy as np
    acmg = _carrier_counts(df, "ACMG")
    amp  = _carrier_counts(df, "AMprimary")
    amo  = _carrier_counts(df, "AMonly")
    syndromes = [s for s in common.SYNDROME_ORDER if s in amp.index]
    x = np.arange(len(syndromes))
    w = 0.27
    ax.bar(x - w, [acmg.get(s, 0) for s in syndromes], width=w,
           color=common.CATEGORY_COLORS["ACMG_PLP"], label=common.CATEGORY_LABELS["ACMG_PLP"])
    ax.bar(x,     [amp.get(s, 0)  for s in syndromes], width=w,
           color=common.CATEGORY_COLORS["AM_primary"], label=common.CATEGORY_LABELS["AM_primary"])
    ax.bar(x + w, [amo.get(s, 0)  for s in syndromes], width=w,
           color=common.CATEGORY_COLORS["AM_only_non_PLP"], label=common.CATEGORY_LABELS["AM_only_non_PLP"])
    ax.set_xticks(x)
    ax.set_xticklabels([common.SYNDROME_SHORT_LABEL.get(s, s) for s in syndromes],
                       rotation=45, ha="right", fontsize=6)
    ax.set_ylabel("# carriers")
    ax.set_title("A  Carrier counts per syndrome  (three definitions)")
    ax.legend(fontsize=6, loc="upper right", ncol=1)


def panel_2B(ax, df) -> None:
    """Stacked-bar share of {ACMG-only, AM-only-non-P/LP, both} per syndrome."""
    import numpy as np
    import pandas as pd
    rows = []
    for col in df.columns:
        if not col.startswith("carrier_") or not col.endswith("_ACMG"):
            continue
        syn = col[len("carrier_"):-len("_ACMG")]
        acmg = df[col].str.upper() == "TRUE"
        amp = df.get(f"carrier_{syn}_AMprimary", "").astype(str).str.upper() == "TRUE"
        both = (acmg & amp).sum()
        only_acmg = (acmg & ~amp).sum()
        only_am = (~acmg & amp).sum()
        total = both + only_acmg + only_am
        if total == 0:
            continue
        rows.append({"syndrome": syn, "ACMG_only": only_acmg / total,
                     "AM_only_non_PLP": only_am / total, "both": both / total})
    if not rows:
        ax.text(0.5, 0.5, "no carriers in cohort", ha="center", va="center",
                transform=ax.transAxes, fontsize=8, color="grey")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title("B  Stacked share of carrier provenance")
        return
    pf = pd.DataFrame(rows).set_index("syndrome")
    syndromes = [s for s in common.SYNDROME_ORDER if s in pf.index]
    pf = pf.reindex(syndromes)
    bottoms = np.zeros(len(pf))
    for key, color_key in [("ACMG_only", "ACMG_PLP"),
                           ("both",      "both"),
                           ("AM_only_non_PLP", "AM_only_non_PLP")]:
        ax.bar(range(len(pf)), pf[key], bottom=bottoms, width=0.7,
               color=common.CATEGORY_COLORS[color_key], label=key.replace("_", " "))
        bottoms = bottoms + pf[key].fillna(0).to_numpy()
    ax.set_xticks(range(len(pf)))
    ax.set_xticklabels([common.SYNDROME_SHORT_LABEL.get(s, s) for s in pf.index],
                       rotation=45, ha="right", fontsize=6)
    ax.set_ylabel("fraction of carriers")
    ax.set_ylim(0, 1)
    ax.set_title("B  Provenance share of each syndrome's carriers")
    ax.legend(fontsize=6, loc="upper right")


def panel_2C(ax, df) -> None:
    """Counterfactual: per-syndrome counts under AM>=0.864 (legacy) vs Chen
    gene-specific. Pairs with Fig 1B."""
    import numpy as np
    amp  = _carrier_counts(df, "AMprimary")
    am0864 = _carrier_counts(df, "AM0864")
    if amp.empty and am0864.empty:
        ax.text(0.5, 0.5, "AM0864 column missing\n(rerun pipeline batch 1)",
                ha="center", va="center", transform=ax.transAxes, color="grey")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title("C  Counterfactual: AM>=0.864 vs Chen gene-specific")
        return
    syndromes = [s for s in common.SYNDROME_ORDER if (s in amp.index or s in am0864.index)]
    x = np.arange(len(syndromes))
    w = 0.4
    ax.bar(x - w/2, [am0864.get(s, 0) for s in syndromes], width=w,
           color=common.CATEGORY_COLORS["AM_global_0864"],
           label=common.CATEGORY_LABELS["AM_global_0864"])
    ax.bar(x + w/2, [amp.get(s, 0)    for s in syndromes], width=w,
           color=common.CATEGORY_COLORS["AM_primary"],
           label=common.CATEGORY_LABELS["AM_primary"])
    ax.set_xticks(x)
    ax.set_xticklabels([common.SYNDROME_SHORT_LABEL.get(s, s) for s in syndromes],
                       rotation=45, ha="right", fontsize=6)
    ax.set_ylabel("# carriers")
    ax.set_title("C  Counterfactual: AM>=0.864 (legacy) vs Chen gene-specific")
    ax.legend(fontsize=6, loc="upper right")


def make(out_dir: str, *, c_table: str) -> None:
    import matplotlib.pyplot as plt
    import pandas as pd
    common.apply_rcparams()
    df = pd.read_csv(c_table, sep="\t", dtype=str, na_filter=False)
    fig, axes = plt.subplots(3, 1, figsize=(common.FIG_WIDTH_DOUBLE[0],
                                            common.FIG_WIDTH_DOUBLE[0] * 1.1))
    panel_2A(axes[0], df)
    panel_2B(axes[1], df)
    panel_2C(axes[2], df)
    fig.tight_layout()
    common.save_both(fig, Path(out_dir), "fig2")
    plt.close(fig)
    LOG.info("fig2: wrote fig2.png + fig2.pdf to %s", out_dir)


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--c-table", required=True, help="C_regression_matrix.tsv")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args(argv)
    make(args.out_dir, c_table=args.c_table)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
