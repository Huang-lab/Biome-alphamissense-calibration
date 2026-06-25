"""Logistic-regression OR table at the syndrome-group level (rebuttal Table 4).

Reads ``results/<cohort>/C_regression_matrix.tsv`` and produces
``stats_syndrome_associations.tsv`` — one row per
(syndrome, phenotype, variant_category, cohort) combination. The model is

    logit(case) = beta_0 + beta_carrier * carrier_flag
                  + beta_age * age + beta_age2 * age^2
                  + beta_sex * sex + sum_k beta_PCk * PC_k

fit per cell with ``statsmodels.api.Logit``. We report exp(beta_carrier) as the
OR, the Wald 95% CI, and the p-value, plus an approximate minimum-detectable-OR
at 80% power for the observed sample/carrier configuration (used in Fig 3B to
answer R1's BRCA-paradox question).

Rows with too few carrier or non-carrier observations to estimate the model
(default ``--min-cell-size 5`` per cell of the 2x2) are listed in the table
with NaN statistics so the matrix stays rectangular for downstream heatmaps.

This is the *minimum viable* slice of run_stats.py. It deliberately only emits
the syndrome-level Table 4; per-cohort phenotype demographics (Table S9), the
ancestry × case/control matrix (Table 3), and the carrier-frequency landscapes
(Table 5 / S10 / S11 / S12) ship in later commits once the figure files that
consume them are in place.
"""
from __future__ import annotations

import argparse
import math
import os
import sys
from typing import Dict, List, Optional, Tuple

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from util import LOG, die, load_config, resolve, write_tsv  # noqa: E402


# ---------------------------------------------------------------------------
# Canonical syndrome -> cancer-site mapping for the `is_canonical` flag.
# Cells outside this mapping are still computed; they just don't get the
# canonical-outline highlight in Fig 3A. Keys MUST match config.gene_groups.
# Values are lowercase substrings — matched against the phenotype `Group`
# column case-insensitively, so "Breast" matches "Breast", "Breast Cancer",
# "breast invasive carcinoma", etc.
# ---------------------------------------------------------------------------
CANONICAL_SYNDROME_CANCERS: Dict[str, List[str]] = {
    "Hereditary_Breast_and_Ovarian_Cancer_Syndrome": ["breast", "ovarian", "prostate", "pancreatic"],
    "Lynch_Syndrome":                                 ["colorectal", "endometrial", "gastric", "ovarian", "urothelial"],
    "Familial_Adenomatous_Polyposis":                 ["colorectal"],
    "MUTYH_Associated_Polyposis":                     ["colorectal"],
    "Multiple_Endocrine_Neoplasia":                   ["thyroid", "parathyroid", "pancreatic"],
    "Hereditary_Paraganglioma_Pheochromocytoma":      ["paraganglioma", "pheochromocytoma"],
    "Li_Fraumeni_Syndrome":                           ["breast", "sarcoma", "brain", "adrenocortical"],
    "Tuberous_Sclerosis_Complex":                     ["kidney", "brain"],
    "Von_Hippel_Lindau_Syndrome":                     ["kidney", "pheochromocytoma", "hemangioblastoma"],
    "PTEN_Hamartoma_Tumor_Syndrome":                  ["breast", "thyroid", "endometrial"],
    "Juvenile_Polyposis_Syndrome":                    ["colorectal", "gastric"],
    "Hereditary_Retinoblastoma":                      ["retinoblastoma", "eye"],
    "Peutz_Jeghers_Syndrome":                         ["colorectal", "breast", "pancreatic"],
    "Neurofibromatosis_Type_2":                       ["brain", "schwannoma", "meningioma"],
    "Wilms_Tumor_Syndrome":                           ["wilms", "kidney"],
}

# Variant-category column suffix (in C_regression_matrix.tsv) -> human label.
VARIANT_CATEGORIES: List[Tuple[str, str]] = [
    ("ACMG",       "ACMG_PLP"),
    ("AMprimary",  "AM_primary"),
    ("AMonly",     "AM_only_non_PLP"),
]

CONTROL_GROUP_TOKENS: List[str] = ["control"]  # lowercase substring match


def _is_canonical(syndrome: str, phenotype: str) -> bool:
    expected = CANONICAL_SYNDROME_CANCERS.get(syndrome, [])
    pl = phenotype.lower()
    return any(tok in pl for tok in expected)


def _is_control_group(group: str) -> bool:
    gl = group.lower()
    return any(tok in gl for tok in CONTROL_GROUP_TOKENS)


def _bh_qvalues(pvals: List[Optional[float]]) -> List[Optional[float]]:
    """Benjamini-Hochberg adjustment over the non-None entries; the None
    entries pass through. Mirrors statsmodels.stats.multitest.multipletests
    with method='fdr_bh' but skips NaN cleanly."""
    n = sum(1 for p in pvals if p is not None and not math.isnan(p))
    if n == 0:
        return list(pvals)
    indexed = sorted(
        ((i, p) for i, p in enumerate(pvals) if p is not None and not math.isnan(p)),
        key=lambda ip: ip[1],
    )
    qs: List[float] = [float("nan")] * len(pvals)
    prev = 1.0
    for rank_from_top, (orig_i, p) in enumerate(reversed(indexed)):
        rank = n - rank_from_top  # 1..n
        q = min(prev, p * n / rank)
        qs[orig_i] = q
        prev = q
    return [qs[i] if pvals[i] is not None and not math.isnan(pvals[i]) else None
            for i in range(len(pvals))]


def _min_detectable_or(n_carrier_case: int, n_carrier_ctrl: int,
                       n_noncarrier_case: int, n_noncarrier_ctrl: int,
                       alpha: float = 0.05, power: float = 0.80) -> Optional[float]:
    """Approximate min-detectable two-sided OR for the observed 2x2 at the
    target power. Uses the standard log-OR variance approximation:

        SE(log OR) ~= sqrt(1/a + 1/b + 1/c + 1/d)

    and min |log OR| = (z_{alpha/2} + z_{1-beta}) * SE. Returns the OR
    (always >= 1.0) or None if any cell is zero (formula is undefined)."""
    from scipy.stats import norm
    cells = [n_carrier_case, n_carrier_ctrl, n_noncarrier_case, n_noncarrier_ctrl]
    if any(c <= 0 for c in cells):
        return None
    se = math.sqrt(sum(1.0 / c for c in cells))
    log_or = (norm.ppf(1 - alpha / 2) + norm.ppf(power)) * se
    return math.exp(log_or)


def _fit_logit(y, X) -> Optional[Tuple[float, float, float, float]]:
    """Return (beta_carrier, se, p, ok). The carrier flag is column 1 by
    convention. Returns None on convergence failure / perfect separation."""
    import numpy as np
    import statsmodels.api as sm
    try:
        # disp=0 silences convergence chatter; method='newton' falls back to
        # 'lbfgs' if it fails, mimicking R's glm() default behavior.
        res = sm.Logit(y, X).fit(disp=0, method="newton", maxiter=50)
    except (np.linalg.LinAlgError, ValueError, Exception):  # noqa: BLE001
        return None
    if not res.mle_retvals.get("converged", True):
        return None
    beta = float(res.params[1])
    se = float(res.bse[1])
    p = float(res.pvalues[1])
    return beta, se, p, 1.0


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=None)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--out-dir", default=None,
                    help="default: results/<cohort>/ from config")
    ap.add_argument("--min-cell-size", type=int, default=5,
                    help="skip cells where any of the 2x2 counts is below this; "
                         "row is still emitted with NaN stats")
    args = ap.parse_args(argv)

    import numpy as np
    import pandas as pd

    cfg = load_config(args.config)
    results_dir = resolve(cfg, cfg["paths"]["results_dir"])
    out_dir = args.out_dir or os.path.join(results_dir, args.cohort)
    c_path = os.path.join(results_dir, args.cohort, "C_regression_matrix.tsv")
    if not os.path.isfile(c_path):
        die(f"C_regression_matrix not found at {c_path} (step 04 must finish first)")

    pc_cols = list(cfg["pcs"]["pc_columns"])
    groups: Dict[str, List[str]] = cfg["gene_groups"]
    syndromes = list(groups.keys())

    df = pd.read_csv(c_path, sep="\t", dtype=str, na_filter=False)
    LOG.info("run_stats: loaded %d samples from %s", len(df), c_path)

    # numerify covariates ----------------------------------------------------
    df["_age"]  = pd.to_numeric(df["Age_at_diagnosis"], errors="coerce")
    df["_age2"] = pd.to_numeric(df["Age2"],            errors="coerce")
    df["_sex_M"] = (df["GENDER"].str.upper() == "M").astype(int)
    for pc in pc_cols:
        df[f"_{pc}"] = pd.to_numeric(df.get(pc, ""), errors="coerce")
    df["_is_control"] = df["Group"].apply(_is_control_group)

    # phenotype universe (cancer Groups only, controls excluded) -------------
    phenotypes = sorted(p for p in df["Group"].unique()
                        if p and not _is_control_group(p))
    LOG.info("run_stats: phenotype Groups detected: %s; controls: %d samples",
             phenotypes, int(df["_is_control"].sum()))

    covariate_cols = ["_age", "_age2", "_sex_M"] + [f"_{pc}" for pc in pc_cols]
    rows_out: List[List] = []

    for syndrome in syndromes:
        for phenotype in phenotypes:
            for col_suffix, label in VARIANT_CATEGORIES:
                carrier_col = f"carrier_{syndrome}_{col_suffix}"
                if carrier_col not in df.columns:
                    continue
                # Subset: cases of this phenotype + ALL controls.
                mask = (df["Group"] == phenotype) | df["_is_control"]
                sub = df.loc[mask].copy()
                sub["_carrier"] = (sub[carrier_col].str.upper() == "TRUE").astype(int)
                sub["_is_case"] = (sub["Group"] == phenotype).astype(int)
                # complete-case for covariates; lose samples with missing age/PCs
                cc_mask = sub[covariate_cols].notna().all(axis=1)
                sub_cc = sub.loc[cc_mask]
                a = int(((sub_cc["_carrier"] == 1) & (sub_cc["_is_case"] == 1)).sum())
                b = int(((sub_cc["_carrier"] == 1) & (sub_cc["_is_case"] == 0)).sum())
                c = int(((sub_cc["_carrier"] == 0) & (sub_cc["_is_case"] == 1)).sum())
                d = int(((sub_cc["_carrier"] == 0) & (sub_cc["_is_case"] == 0)).sum())
                n_total = a + b + c + d
                min_obs = min(a, b, c, d)

                row: List = [
                    args.cohort, syndrome, phenotype, label,
                    n_total, a + b, a, c, b, d,
                    _is_canonical(syndrome, phenotype),
                ]
                if min_obs < args.min_cell_size:
                    row += [float("nan")] * 5  # OR, lo, hi, p, mdo
                    row += ["age+age2+sex+" + "+".join(pc_cols)]
                    rows_out.append(row)
                    continue

                y = sub_cc["_is_case"].to_numpy(dtype=float)
                X = sub_cc[["_carrier"] + covariate_cols].to_numpy(dtype=float)
                # statsmodels wants the intercept explicit
                X = np.column_stack([np.ones(len(X)), X])  # so carrier is col 1
                fit = _fit_logit(y, X)
                if fit is None:
                    row += [float("nan")] * 4
                else:
                    beta, se, p, _ = fit
                    or_ = math.exp(beta)
                    lo = math.exp(beta - 1.96 * se)
                    hi = math.exp(beta + 1.96 * se)
                    row += [or_, lo, hi, p]
                mdo = _min_detectable_or(a, b, c, d)
                row.append(mdo if mdo is not None else float("nan"))
                row.append("age+age2+sex+" + "+".join(pc_cols))
                rows_out.append(row)

    # BH q across (cohort, variant_category) groups --------------------------
    df_out = pd.DataFrame(rows_out, columns=[
        "cohort", "syndrome", "phenotype", "variant_category",
        "n_total", "n_carriers", "n_cases_carriers", "n_cases_noncarriers",
        "n_ctrl_carriers", "n_ctrl_noncarriers", "is_canonical",
        "OR", "OR_95CI_lo", "OR_95CI_hi", "p_value",
        "min_detectable_OR_80pct_power", "covariates",
    ])
    df_out["q_BH"] = float("nan")
    for cat in df_out["variant_category"].unique():
        sel = df_out["variant_category"] == cat
        pv = df_out.loc[sel, "p_value"].tolist()
        pv_for_bh: List[Optional[float]] = [
            (None if (p is None or (isinstance(p, float) and math.isnan(p))) else float(p))
            for p in pv
        ]
        qv = _bh_qvalues(pv_for_bh)
        df_out.loc[sel, "q_BH"] = [float("nan") if q is None else q for q in qv]

    # output -----------------------------------------------------------------
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "stats_syndrome_associations.tsv")
    df_out = df_out[[
        "cohort", "syndrome", "phenotype", "variant_category",
        "n_total", "n_carriers",
        "n_cases_carriers", "n_cases_noncarriers",
        "n_ctrl_carriers", "n_ctrl_noncarriers",
        "OR", "OR_95CI_lo", "OR_95CI_hi", "p_value", "q_BH",
        "min_detectable_OR_80pct_power", "is_canonical", "covariates",
    ]]
    df_out.to_csv(out_path, sep="\t", index=False, na_rep="NaN")
    n_fit = int(df_out["OR"].notna().sum())
    n_can = int(df_out["is_canonical"].sum())
    LOG.info("run_stats: wrote %d rows to %s  (n_fit_ok=%d, canonical cells=%d)",
             len(df_out), out_path, n_fit, n_can)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
