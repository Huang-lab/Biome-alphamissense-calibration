"""Supplementary figures S1–S7 — one module so they share imports and styling.

Each ``sX(...)`` function renders one panel; the ``make(out_dir, ...)`` entry
point loops over the configured supplementaries. Stubs render a placeholder
panel + emit a logging warning when the underlying data is not yet wired in;
this keeps the figure pipeline runnable end-to-end as supplementaries land.

  S1  Age-cutoff sensitivity (Action 9, R1) — needs run_stats output for
      multiple ``--control-age-min`` settings.
  S2  AM false-negative variants vs the ANNOVAR/InterVar P/LP set
      (Action 13). Reads the variant-carrier TSV directly.
  S3  ClinVar VUS reclassification waterfall.
  S4  Per-cohort variant-level vs carrier-level concordance.
  S5  Threshold-policy sensitivity (Supporting vs Moderate vs Strong).
  S6  Pipeline replication audit (results vs prior run).
  S7  Sample inclusion / exclusion CONSORT.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from util import LOG  # noqa: E402

from . import common


SUPPL_LABELS: Dict[str, str] = {
    "S1": "Age-cutoff sensitivity",
    "S2": "AM false-negatives vs ANNOVAR/InterVar P/LP",
    "S3": "ClinVar VUS reclassification waterfall",
    "S4": "Variant-level vs carrier-level concordance",
    "S5": "Threshold-policy sensitivity",
    "S6": "Pipeline replication audit",
    "S7": "CONSORT (sample inclusion / exclusion)",
}


def _placeholder_fig(ax, sid: str, msg: str = "") -> None:
    ax.set_facecolor("#f5f5f5")
    label = SUPPL_LABELS.get(sid, sid)
    txt = f"Fig {sid} — {label}"
    if msg:
        txt += f"\n\n{msg}"
    else:
        txt += "\n\n(placeholder until data wired in)"
    ax.text(0.5, 0.5, txt, ha="center", va="center", fontsize=8,
            transform=ax.transAxes, color="#333333", wrap=True)
    ax.set_xticks([]); ax.set_yticks([])


def _save_single(out_dir: Path, sid: str, render) -> None:
    import matplotlib.pyplot as plt
    common.apply_rcparams()
    fig, ax = plt.subplots(figsize=common.FIG_WIDTH_SINGLE)
    render(ax)
    fig.tight_layout()
    common.save_both(fig, out_dir, f"fig{sid}")
    plt.close(fig)


def make(out_dir: str, **inputs) -> None:
    out_dir = Path(out_dir)
    for sid in SUPPL_LABELS:
        _save_single(out_dir, sid, lambda ax, _sid=sid: _placeholder_fig(ax, _sid))
    LOG.info("suppl: emitted %d placeholder figures to %s",
             len(SUPPL_LABELS), out_dir)


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args(argv)
    make(args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
