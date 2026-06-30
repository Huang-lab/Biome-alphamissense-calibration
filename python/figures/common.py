"""Shared figure styling, palette, and canonical orderings.

Import everything from here so the figure set reads as one paper:

    from python.figures import common
    common.apply_rcparams()
    fig, ax = plt.subplots(figsize=common.FIG_WIDTH_DOUBLE)

The constants here mirror the plan in
``/root/.claude/plans/you-are-generating-a-graceful-honey.md``. If you
add a new syndrome or carrier category to ``config/config.yaml``, mirror
it in :data:`SYNDROME_ORDER` / :data:`VARIANT_CATEGORIES` here so the
figures pick it up.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

# Figure dimensions in inches; matches Nature double-column width.
FIG_WIDTH_DOUBLE: Tuple[float, float] = (7.2, 4.0)
FIG_WIDTH_SINGLE: Tuple[float, float] = (3.5, 3.0)

# Qualitative palette for the four carrier definitions + their "both" overlap.
# Chosen for colorblind-friendliness (Wong 2011). Keys match the
# `variant_category` column values written by run_stats.py.
CATEGORY_COLORS: Dict[str, str] = {
    "ACMG_PLP":              "#0072B2",  # blue
    "AM_calibrated":         "#E69F00",  # amber — Chen per-variant lookup, evidence ≥ PP3_Moderate
    "AM_calibrated_not_PLP": "#009E73",  # green — AM_calibrated AND NOT ACMG_PLP (novel set)
    "both":                  "#56B4E9",  # sky blue (Fig 2B stacked share)
    "AM_global_0.864":       "#CC79A7",  # rose — legacy single-cutoff counterfactual
}

# Display labels — used in legends and axes. The 0.864 is in the label
# itself so reviewers don't have to look up which threshold "global" means.
CATEGORY_LABELS: Dict[str, str] = {
    "ACMG_PLP":              "ACMG/AMP P/LP",
    "AM_calibrated":         "AM-calibrated (Chen PP3≥Moderate)",
    "AM_calibrated_not_PLP": "AM-calibrated, not P/LP",
    "AM_global_0.864":       "AM ≥ 0.864 (legacy)",
}

# Canonical syndrome display order — top-to-bottom for the OR heatmap (Fig 3A).
# Driven by clinical importance + carrier yield in BioMe; mirrors the order in
# config.gene_groups so this list is the single source of truth for figures.
SYNDROME_ORDER: List[str] = [
    "Hereditary_Breast_and_Ovarian_Cancer_Syndrome",
    "Lynch_Syndrome",
    "Familial_Adenomatous_Polyposis",
    "MUTYH_Associated_Polyposis",
    "Multiple_Endocrine_Neoplasia",
    "Hereditary_Paraganglioma_Pheochromocytoma",
    "Li_Fraumeni_Syndrome",
    "PTEN_Hamartoma_Tumor_Syndrome",
    "Tuberous_Sclerosis_Complex",
    "Von_Hippel_Lindau_Syndrome",
    "Juvenile_Polyposis_Syndrome",
    "Peutz_Jeghers_Syndrome",
    "Neurofibromatosis_Type_2",
    "Hereditary_Retinoblastoma",
    "Wilms_Tumor_Syndrome",
]

# Compact label for the y-axis of Fig 3A (full name is too long).
SYNDROME_SHORT_LABEL: Dict[str, str] = {
    "Hereditary_Breast_and_Ovarian_Cancer_Syndrome": "HBOC",
    "Lynch_Syndrome":                                "Lynch",
    "Familial_Adenomatous_Polyposis":                "FAP",
    "MUTYH_Associated_Polyposis":                    "MUTYH-AP",
    "Multiple_Endocrine_Neoplasia":                  "MEN",
    "Hereditary_Paraganglioma_Pheochromocytoma":     "PGL/PHEO",
    "Li_Fraumeni_Syndrome":                          "Li-Fraumeni",
    "PTEN_Hamartoma_Tumor_Syndrome":                 "PHTS",
    "Tuberous_Sclerosis_Complex":                    "TSC",
    "Von_Hippel_Lindau_Syndrome":                    "VHL",
    "Juvenile_Polyposis_Syndrome":                   "JPS",
    "Peutz_Jeghers_Syndrome":                        "Peutz-Jeghers",
    "Neurofibromatosis_Type_2":                      "NF2",
    "Hereditary_Retinoblastoma":                     "Retinoblastoma",
    "Wilms_Tumor_Syndrome":                          "Wilms",
}

# Canonical 28-gene panel order (groups stay contiguous so the gene axis of
# Fig 6B / 6E reads visually grouped).
GENE_ORDER: List[str] = [
    # HBOC
    "BRCA1", "BRCA2", "PALB2",
    # Lynch
    "MLH1", "MSH2", "MSH6", "PMS2",
    # FAP / MUTYH-AP
    "APC", "MUTYH",
    # MEN
    "MEN1", "RET",
    # Paraganglioma-Pheo
    "SDHAF2", "SDHB", "SDHC", "SDHD", "TMEM127", "MAX",
    # Li-Fraumeni / PHTS / TSC / VHL
    "TP53", "PTEN", "TSC1", "TSC2", "VHL",
    # JPS / Peutz-Jeghers / NF2 / Retinoblastoma / Wilms
    "BMPR1A", "SMAD4", "STK11", "NF2", "RB1", "WT1",
]

# Ancestry display order (used by Fig 4 + Fig 6C/D/E). Matches the
# genetically_determined values BioMe uses.
ANCESTRY_ORDER: List[str] = ["EUR", "AFR", "HIS", "EAS", "SAS"]


def apply_rcparams() -> None:
    """Apply paper-grade matplotlib defaults. Idempotent; safe to call again."""
    import matplotlib as mpl

    mpl.rcParams.update({
        "figure.dpi":         150,
        "savefig.dpi":        300,
        "savefig.bbox":       "tight",
        "savefig.pad_inches": 0.05,
        "font.family":        "sans-serif",
        "font.sans-serif":    ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size":          8,
        "axes.titlesize":     9,
        "axes.labelsize":     8,
        "axes.spines.top":    False,
        "axes.spines.right":  False,
        "axes.linewidth":     0.6,
        "xtick.labelsize":    7,
        "ytick.labelsize":    7,
        "xtick.major.width":  0.6,
        "ytick.major.width":  0.6,
        "legend.fontsize":    7,
        "legend.frameon":     False,
        "lines.linewidth":    1.2,
        "lines.markersize":   3.5,
        "pdf.fonttype":       42,  # embed Type-3-free fonts so reviewers can edit text
        "ps.fonttype":        42,
    })


def save_both(fig, out_dir: Path, name: str) -> None:
    """Write ``out_dir/<name>.png`` and ``out_dir/<name>.pdf``."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{name}.png")
    fig.savefig(out_dir / f"{name}.pdf")
