"""
Shared utilities, style, and data loaders for BioMe AlphaMissense figures.

Variant class labels (used in all new figures):
  "ACMG P/LP"              – carriers in PLPorPTV file (P/LP or P/LP & PTV rows)
  "AM_calibrated"          – carriers with is_AM_carrier_primary==TRUE (may overlap ACMG)
  "AM_calibrated not P/LP" – AM_calibrated carriers NOT in ACMG P/LP list (novel only)
  "both"                   – in both lists; only used for Fig 1 Venn diagram
"""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT    = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RESULTS_I    = os.path.join(REPO_ROOT, "results", "cohortI")
RESULTS_II   = os.path.join(REPO_ROOT, "results", "cohortII")
RESULTS_DIR  = os.path.join(REPO_ROOT, "results")
STATS_DIR    = os.path.join(REPO_ROOT, "results", "stats")
METADATA_DIR = os.path.join(REPO_ROOT, "metadata")
REFS_DIR     = os.path.join(REPO_ROOT, "refs")
FIGS_DIR     = os.environ.get("BIOME_FIGS_DIR",
                               os.path.join(REPO_ROOT, "results", "figures_solely"))
TABLES_DIR   = os.environ.get("BIOME_TABLES_DIR",
                               os.path.join(REPO_ROOT, "results", "tables"))

METADATA_FILE = {
    "cohortI":  os.path.join(METADATA_DIR, "RegenWXS_HX_Newgroups.250109.tsv"),
    "cohortII": os.path.join(METADATA_DIR, "Sema4WXS_HX_Newgroups.250109.tsv"),
}

ACMG_FILE = {
    "cohortI":  os.path.join(METADATA_DIR, "Regen_VariantsInSamplesPLPorPTV.tsv"),
    "cohortII": os.path.join(METADATA_DIR, "Sema4_VariantsInSamplesPLPorPTV.tsv"),
}

PCS_FILE    = os.path.join(REFS_DIR, "GSA_GDA_PCA_V2.txt")
MRN_MAP_FILE = os.path.join(REFS_DIR, "Masked_mrn_map.txt")

for _d in (STATS_DIR, FIGS_DIR, TABLES_DIR):
    os.makedirs(_d, exist_ok=True)

# ---------------------------------------------------------------------------
# Publication style
# ---------------------------------------------------------------------------
FONT_FAMILY = "Arial"
plt.rcParams.update({
    "font.family":        FONT_FAMILY,
    "font.size":          14,
    "axes.labelsize":     16,
    "axes.titlesize":     14,
    "xtick.labelsize":    14,
    "ytick.labelsize":    14,
    "legend.fontsize":    13,
    "figure.dpi":         150,
    "savefig.dpi":        300,
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "axes.linewidth":     1.0,
    "xtick.major.width":  1.0,
    "ytick.major.width":  1.0,
    "lines.linewidth":    1.2,
    "pdf.fonttype":       42,
    "ps.fonttype":        42,
})

FIGW_DOUBLE = 7.2   # double-column journal width (inches)
FIGW_SINGLE = 3.5   # single-column

COHORT_META = {
    "cohortI":  {"label": "Cohort I",  "n": 28310},
    "cohortII": {"label": "Cohort II", "n": 13967},
}

# ---------------------------------------------------------------------------
# Variant class constants
# ---------------------------------------------------------------------------
VC_ACMG    = "ACMG P/LP"
VC_AM      = "AM_calibrated"
VC_AM_ONLY = "AM_calibrated not P/LP"
VC_BOTH    = "both"           # used only for Venn diagram in fig1

VARIANT_CLASSES = [VC_ACMG, VC_AM, VC_AM_ONLY]   # display order in all figures

COLORS = {
    VC_ACMG:    "#4D4D4D",   # dark gray
    VC_AM:      "#2166AC",   # blue
    VC_AM_ONLY: "#D6604D",   # red-orange
    VC_BOTH:    "#D1E5F0",   # light blue (Venn only)
    # legacy keys kept for archive scripts that may import this module
    "ACMG":         "#4D4D4D",
    "AMcalibrated": "#2166AC",
    "AMonly":       "#D6604D",
    "both":         "#D1E5F0",
    "old_0864":     "#B2ABD2",
    "single_gene":  "#2166AC",
    "domain_agg":   "#F4A582",
    "uncovered":    "#CCCCCC",
}

# ---------------------------------------------------------------------------
# Control / phenotype group definitions
# ---------------------------------------------------------------------------
CONTROL_GROUPS_RAW = {
    "Control (age≥50)", "Control (age<50)",
    "Family History (age≥50)", "Family History (age<50)",
}

# For regression: merge both age bins into one control and one FamHX group
CONTROL_LABEL  = "Control"
FAMHX_LABEL    = "Family History"

N_MIN_PHENOTYPE = 45   # minimum cases per Group to run regression
N_MIN_CARRIERS  = 3    # minimum carriers in both case and control

# ---------------------------------------------------------------------------
# Syndrome definitions
# ---------------------------------------------------------------------------
SYNDROMES = {
    "Hereditary_Breast_and_Ovarian_Cancer_Syndrome": {
        "label": "HBOC",
        "genes": ["BRCA1", "BRCA2", "PALB2"],
    },
    "Lynch_Syndrome": {
        "label": "Lynch",
        "genes": ["MLH1", "MSH2", "MSH6", "PMS2"],
    },
    "Familial_Adenomatous_Polyposis": {
        "label": "FAP",
        "genes": ["APC"],
    },
    "Multiple_Endocrine_Neoplasia": {
        "label": "MEN",
        "genes": ["MEN1", "RET"],
    },
    "Hereditary_Paraganglioma_Pheochromocytoma": {
        "label": "Para-Pheo",
        "genes": ["SDHAF2", "SDHB", "SDHC", "SDHD", "TMEM127", "MAX"],
    },
    "Li_Fraumeni_Syndrome": {
        "label": "LFS",
        "genes": ["TP53"],
    },
    "Tuberous_Sclerosis_Complex": {
        "label": "TSC",
        "genes": ["TSC1", "TSC2"],
    },
    "Von_Hippel_Lindau_Syndrome": {
        "label": "VHL",
        "genes": ["VHL"],
    },
    "PTEN_Hamartoma_Tumor_Syndrome": {
        "label": "PTEN-HTS",
        "genes": ["PTEN"],
    },
    "Juvenile_Polyposis_Syndrome": {
        "label": "JPS",
        "genes": ["BMPR1A", "SMAD4"],
    },
    "Hereditary_Retinoblastoma": {
        "label": "Rb",
        "genes": ["RB1"],
    },
    "Peutz_Jeghers_Syndrome": {
        "label": "PJS",
        "genes": ["STK11"],
    },
    "Neurofibromatosis_Type_2": {
        "label": "NF2",
        "genes": ["NF2"],
    },
    "Wilms_Tumor_Syndrome": {
        "label": "Wilms",
        "genes": ["WT1"],
    },
}

SYNDROME_ORDER  = list(SYNDROMES.keys())
SYNDROME_LABELS = {k: v["label"] for k, v in SYNDROMES.items()}

# One distinguishable color per syndrome (Set1 + Paired palette)
SYNDROME_COLORS = {
    "Hereditary_Breast_and_Ovarian_Cancer_Syndrome": "#E41A1C",
    "Lynch_Syndrome":                                "#377EB8",
    "Familial_Adenomatous_Polyposis":                "#4DAF4A",
    "Multiple_Endocrine_Neoplasia":                  "#984EA3",
    "Hereditary_Paraganglioma_Pheochromocytoma":     "#FF7F00",
    "Li_Fraumeni_Syndrome":                          "#A65628",
    "Tuberous_Sclerosis_Complex":                    "#F781BF",
    "Von_Hippel_Lindau_Syndrome":                    "#808080",
    "PTEN_Hamartoma_Tumor_Syndrome":                 "#1B9E77",
    "Juvenile_Polyposis_Syndrome":                   "#FC8D62",
    "Hereditary_Retinoblastoma":                     "#8DA0CB",
    "Peutz_Jeghers_Syndrome":                        "#E78AC3",
    "Neurofibromatosis_Type_2":                      "#A6D854",
    "Wilms_Tumor_Syndrome":                          "#E6AB02",
}

# Build gene → syndrome lookup
GENE_TO_SYNDROME = {}
for _skey, _sval in SYNDROMES.items():
    for _g in _sval["genes"]:
        GENE_TO_SYNDROME[_g] = _skey

ALL_GENES_ORDERED = []
for _skey in SYNDROME_ORDER:
    ALL_GENES_ORDERED.extend(SYNDROMES[_skey]["genes"])

# Canonical gene-disease pairs (kept for reference)
CANONICAL_PAIRS = [
    ("Hereditary_Breast_and_Ovarian_Cancer_Syndrome", "Breast"),
    ("Lynch_Syndrome",                                "Colon/Rectum"),
    ("Lynch_Syndrome",                                "Uterus/Endometrium"),
    ("Familial_Adenomatous_Polyposis",                "Colon/Rectum"),
    ("Multiple_Endocrine_Neoplasia",                  "Malignant neoplasm of thyroid gland"),
    ("Von_Hippel_Lindau_Syndrome",                    "Kidney"),
]

# ---------------------------------------------------------------------------
# Cohort helpers
# ---------------------------------------------------------------------------

def cohort_label(cohort):
    return COHORT_META.get(cohort, {"label": cohort})["label"]


def output_name(base, cohort=None):
    return f"{base}_{cohort}" if cohort else base


# ---------------------------------------------------------------------------
# Low-level I/O
# ---------------------------------------------------------------------------

def _tsv(path, **kw):
    return pd.read_csv(path, sep="\t", low_memory=False, **kw)


# ---------------------------------------------------------------------------
# Primary data loaders (new pipeline — no C_regression_matrix)
# ---------------------------------------------------------------------------

def load_metadata(cohort="cohortI"):
    """Full participant metadata for a cohort."""
    return _tsv(METADATA_FILE[cohort])


def load_metadata_groups(cohort="cohortI"):
    """Deduplicated DataFrame[sample_id, Group, Age_at_diagnosis, GENDER,
    genetically_determined] covering all recruited participants."""
    meta = load_metadata(cohort)
    if cohort == "cohortI":
        id_col = "SINAI_ID"
    else:
        id_col = "MASKED_MRN"
    cols = [id_col, "Group", "Age_at_diagnosis", "GENDER", "genetically_determined"]
    df = meta[cols].rename(columns={id_col: "sample_id"}).copy()
    df["sample_id"] = df["sample_id"].astype(str)
    return df.drop_duplicates(subset=["sample_id", "Group"]).reset_index(drop=True)


def load_acmg_carriers(cohort="cohortI"):
    """Load ACMG P/LP carrier variants (excludes pure PTV rows).

    Returns DataFrame with columns:
        sample_id, gene, chrom, pos, ref, alt
    """
    df = _tsv(ACMG_FILE[cohort])
    # Normalise column names (Regen and Sema4 differ slightly)
    df.columns = [c.lstrip("#") for c in df.columns]

    if cohort == "cohortI":
        # Regen: CHROM POS REF ALT Gene ExonicFunc_refGene SINAI_ID carrier annotation ...
        df = df.rename(columns={"CHROM": "chrom", "POS": "pos",
                                 "REF": "ref", "ALT": "alt",
                                 "Gene": "gene", "SINAI_ID": "sample_id"})
    else:
        # Sema4: CHROM POS REF ALT Gene PatientID carrier annotation SINAI_ID MASKED_MRN ...
        # cohortII sample_id matches A_variant_level_per_person.tsv: MASKED_MRN (numeric)
        df = df.rename(columns={"CHROM": "chrom", "POS": "pos",
                                 "REF": "ref", "ALT": "alt",
                                 "Gene": "gene", "MASKED_MRN": "sample_id"})
        df["sample_id"] = pd.to_numeric(df["sample_id"], errors="coerce").astype("Int64")

    # Keep only P/LP rows (drop pure PTV)
    ann_col = "annotation"
    df = df[df[ann_col].str.contains("P/LP", na=False)]
    return df[["sample_id", "gene", "chrom", "pos", "ref", "alt"]].drop_duplicates()


def load_am_carriers(cohort="cohortI"):
    """Load AM-calibrated carrier variants (is_AM_carrier_primary == TRUE).

    Returns DataFrame with columns:
        sample_id, gene, chrom (chr), pos, ref, alt
    """
    base = RESULTS_I if cohort == "cohortI" else RESULTS_II
    df = _tsv(os.path.join(base, "A_variant_level_per_person.tsv"))
    # Normalise boolean column
    col = "is_AM_carrier_primary"
    if df[col].dtype == object:
        df[col] = df[col].map({"TRUE": True, "FALSE": False,
                                "True": True, "False": False}).astype(bool)
    df = df[df[col] == True].copy()
    df = df.rename(columns={"chr": "chrom"})
    return df[["sample_id", "gene", "chrom", "pos", "ref", "alt"]].drop_duplicates()


def load_pcs(cohort="cohortI"):
    """Load PC1–PC10 for each sample.

    CohortI: PCs are keyed by MASKED_MRN; bridge via refs/Masked_mrn_map.txt
             RGNID → MASKED_MRN, then match to SINAI_ID via metadata.
    CohortII: PCs keyed by MASKED_MRN, which is already the sample_id.
    Returns DataFrame[sample_id, PC1, ..., PC10].
    """
    pc_cols = ["PC1", "PC2", "PC3", "PC4", "PC5",
               "PC6", "PC7", "PC8", "PC9", "PC10"]
    pcs = pd.read_csv(PCS_FILE, sep=r"\s+", engine="python")
    # ID column is ID2 (same as ID1 in this file = MASKED_MRN numeric)
    pcs = pcs[["ID2"] + [c for c in pc_cols if c in pcs.columns]].copy()
    pcs = pcs.rename(columns={"ID2": "masked_mrn"})
    pcs["masked_mrn"] = pcs["masked_mrn"].astype("Int64")

    if cohort == "cohortI":
        bridge = _tsv(MRN_MAP_FILE)
        # columns: PLATFORM, RGNID, MASKED_MRN
        bridge["masked_mrn"] = pd.to_numeric(bridge["MASKED_MRN"], errors="coerce").astype("Int64")
        bridge = bridge.rename(columns={"RGNID": "sample_id"})
        df = bridge[["sample_id", "masked_mrn"]].merge(
            pcs, on="masked_mrn", how="inner"
        )
    else:
        # cohortII sample_id IS the MASKED_MRN (numeric string)
        pcs["sample_id"] = pcs["masked_mrn"].astype(str)
        df = pcs.drop(columns=["masked_mrn"]).copy()

    return df[["sample_id"] + pc_cols].drop_duplicates(subset=["sample_id"])


def build_variant_table(cohort="cohortI"):
    """Build a per-variant-per-sample table with variant_class labels.

    variant_class values:
        "ACMG P/LP"              – in ACMG list only
        "AM_calibrated"          – in AM list only (or in both — AM_calibrated includes all AM)
        "AM_calibrated not P/LP" – in AM list but NOT in ACMG
        "both"                   – in both lists (used for Venn counts)

    Returns DataFrame[sample_id, gene, syndrome, chrom, pos, ref, alt,
                       in_acmg, in_am, variant_class]
    """
    acmg = load_acmg_carriers(cohort).copy()
    am   = load_am_carriers(cohort).copy()

    # Coerce both sample_id columns to the same string type before merging
    acmg["sample_id"] = acmg["sample_id"].astype(str)
    am["sample_id"]   = am["sample_id"].astype(str)

    acmg["in_acmg"] = True
    am["in_am"]     = True

    # Variant key for matching
    key_cols = ["sample_id", "gene", "pos", "ref", "alt"]

    merged = pd.merge(
        acmg[key_cols + ["chrom", "in_acmg"]],
        am[key_cols + ["in_am"]],
        on=key_cols, how="outer"
    )
    merged["in_acmg"] = merged["in_acmg"].notna() & (merged["in_acmg"] == True)
    merged["in_am"]   = merged["in_am"].notna()   & (merged["in_am"]   == True)

    def classify(row):
        if row["in_acmg"] and row["in_am"]:
            return VC_BOTH
        elif row["in_acmg"]:
            return VC_ACMG
        else:
            return VC_AM_ONLY

    merged["variant_class"] = merged.apply(classify, axis=1)

    # Add syndrome
    merged["syndrome"] = merged["gene"].map(GENE_TO_SYNDROME)

    # Fill missing chrom from AM table (ACMG-only rows didn't join chrom from am)
    am_chrom = am.set_index(key_cols)["chrom"].to_dict()
    merged["chrom"] = merged["chrom"].astype(object)
    missing_mask = merged["chrom"].isna()
    if missing_mask.any():
        keys = merged.loc[missing_mask, key_cols].apply(tuple, axis=1)
        merged.loc[missing_mask, "chrom"] = keys.map(lambda k: am_chrom.get(k))

    return merged.reset_index(drop=True)


def build_carrier_matrix(cohort="cohortI"):
    """Build a per-sample carrier matrix joined with metadata and PCs.

    For each sample, compute syndrome-level carrier flags for 3 variant classes:
        carrier_{syndrome}_{vc}  where vc in {ACMG_PLP, AM_calibrated, AM_not_PLP}

    Also merges control groups (age bins combined) and family history groups.

    Returns DataFrame with columns:
        sample_id, Group, Group_merged, Age_at_diagnosis, Age2, GENDER,
        genetically_determined, PC1..PC10,
        carrier_{syn}_{ACMG_PLP|AM_calibrated|AM_not_PLP}  (bool)
    """
    vt = build_variant_table(cohort)

    # --- ACMG P/LP carriers (ACMG-only + both) ---
    acmg_carriers = vt[vt["in_acmg"] == True][["sample_id", "gene", "syndrome"]].drop_duplicates()

    # --- AM_calibrated carriers (AM-only + both) ---
    am_carriers = vt[vt["in_am"] == True][["sample_id", "gene", "syndrome"]].drop_duplicates()

    # --- AM_calibrated not P/LP (AM-only) ---
    am_only_carriers = vt[vt["variant_class"] == VC_AM_ONLY][["sample_id", "gene", "syndrome"]].drop_duplicates()

    meta = load_metadata_groups(cohort)
    try:
        pcs = load_pcs(cohort)
        meta = meta.merge(pcs, on="sample_id", how="left")
    except Exception as e:
        print(f"  [warn] Could not load PCs for {cohort}: {e}")
        for pc in [f"PC{i}" for i in range(1, 11)]:
            meta[pc] = np.nan

    meta["Age2"] = meta["Age_at_diagnosis"] ** 2

    # Merge control groups
    def merge_group(g):
        if g in ("Control (age≥50)", "Control (age<50)"):
            return CONTROL_LABEL
        if g in ("Family History (age≥50)", "Family History (age<50)"):
            return FAMHX_LABEL
        return g

    meta["Group_merged"] = meta["Group"].apply(merge_group)

    # Build syndrome-level carrier columns
    all_sample_ids = meta["sample_id"].unique()
    carrier_df = pd.DataFrame({"sample_id": all_sample_ids})

    for syn_key in SYNDROME_ORDER:
        for vc_label, vc_df in [
            ("ACMG_PLP",      acmg_carriers),
            ("AM_calibrated", am_carriers),
            ("AM_not_PLP",    am_only_carriers),
        ]:
            syn_carriers = vc_df[vc_df["syndrome"] == syn_key]["sample_id"].unique()
            col = f"carrier_{syn_key}_{vc_label}"
            carrier_df[col] = carrier_df["sample_id"].isin(syn_carriers)

    result = meta.merge(carrier_df, on="sample_id", how="left")
    bool_cols = [c for c in result.columns if c.startswith("carrier_")]
    result[bool_cols] = result[bool_cols].fillna(False)
    return result


# ---------------------------------------------------------------------------
# Legacy loaders (kept for compatibility with archive scripts)
# ---------------------------------------------------------------------------

def load_table_A(cohort="cohortI"):
    base = RESULTS_I if cohort == "cohortI" else RESULTS_II
    return _tsv(os.path.join(base, "A_variant_level_per_person.tsv"))


def load_table_B(cohort="cohortI"):
    base = RESULTS_I if cohort == "cohortI" else RESULTS_II
    return _tsv(os.path.join(base, "B_ACMG_vs_AM_comparison.tsv"))


def load_table_B_summary(cohort="cohortI"):
    base = RESULTS_I if cohort == "cohortI" else RESULTS_II
    return _tsv(os.path.join(base, "B_summary_by_gene.tsv"))


# ---------------------------------------------------------------------------
# Figure save helpers
# ---------------------------------------------------------------------------

def save_fig(fig, name, tight=True, cohort=None):
    """Save figure as PNG (300 dpi) to results/figures/."""
    if tight:
        fig.tight_layout()
    name = output_name(name, cohort)
    png_path = os.path.join(FIGS_DIR, f"{name}.png")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    print(f"  Saved: {png_path}")
    plt.close(fig)


def save_panel(fig, panel_name, tight=True, cohort=None):
    """Alias for save_fig — saves a single panel."""
    save_fig(fig, panel_name, tight=tight, cohort=cohort)


def save_table(df, table_name, cohort=None):
    """Save a DataFrame as TSV to results/tables/."""
    table_name = output_name(table_name, cohort)
    out = os.path.join(TABLES_DIR, f"{table_name}.tsv")
    df.to_csv(out, sep="\t", index=False)
    print(f"  Saved table: {out}")


# ---------------------------------------------------------------------------
# Wilson 95% CI helper
# ---------------------------------------------------------------------------

def wilson_ci(n_pos, n_total, z=1.96):
    """Return (low, high) Wilson confidence interval for a proportion."""
    if n_total == 0:
        return (0.0, 0.0)
    p = n_pos / n_total
    denom = 1 + z**2 / n_total
    centre = (p + z**2 / (2 * n_total)) / denom
    half   = z * np.sqrt(p * (1 - p) / n_total + z**2 / (4 * n_total**2)) / denom
    return (max(0, centre - half), min(1, centre + half))


# ---------------------------------------------------------------------------
# Legacy forest-plot helper (kept for reference)
# ---------------------------------------------------------------------------

def forest_plot(ax, rows, y_labels, colors=None, dot_size=40,
                xlog=True, xlim=None, vline=1.0, title=None):
    if colors is None:
        colors = [COLORS["AM_calibrated"]] * len(rows)
    n = len(rows)
    ys = list(range(n - 1, -1, -1))
    for i, (row, y) in enumerate(zip(rows, ys)):
        color     = row.get("color", colors[i])
        up        = row.get("underpowered", False)
        sig       = row.get("significant", True)
        dot_color = color if sig and not up else "#AAAAAA"
        if not up:
            or_val = row["OR"]
            ax.errorbar(or_val, y,
                        xerr=[[or_val - row["CI_low"]], [row["CI_high"] - or_val]],
                        fmt="o", color=dot_color, ecolor=dot_color,
                        ms=5, capsize=2, lw=0.8, zorder=3)
    ax.axvline(vline, color="#888888", lw=0.8, ls="--", zorder=1)
    ax.set_yticks(ys)
    ax.set_yticklabels(y_labels, fontsize=7)
    if xlog:
        ax.set_xscale("log")
    if xlim:
        ax.set_xlim(xlim)
    if title:
        ax.set_title(title, fontsize=8)
    ax.set_xlabel("Odds Ratio (95% CI)", fontsize=7)
    ax.grid(axis="x", ls=":", lw=0.5, alpha=0.5)
