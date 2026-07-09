"""
Logistic regression of syndrome carrier status vs phenotype groups.

Builds carrier matrix from primary data (ACMG PLPorPTV + A_variant_level),
runs logistic regression for each (syndrome, variant_class, phenotype) pair
with sufficient carriers (≥3 in cases AND ≥3 in controls), applies BH FDR,
and writes results/cohortI/regression_results.tsv.

Usage:
    python python/run_regression.py --cohort cohortI
    python python/run_regression.py --cohort cohortII
"""

import argparse
import sys
import os
import warnings

import numpy as np
import pandas as pd

# Allow running as a script from repo root
sys.path.insert(0, os.path.dirname(__file__))
from figures.common import (
    build_carrier_matrix, SYNDROME_ORDER, SYNDROME_LABELS,
    CONTROL_LABEL, FAMHX_LABEL, N_MIN_PHENOTYPE, N_MIN_CARRIERS,
    RESULTS_I, RESULTS_II,
)


def _bh_correction(pvals):
    """Benjamini-Hochberg FDR correction. Returns q-value array."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    if n == 0:
        return pvals.copy()
    order = np.argsort(pvals)
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n + 1)
    q = pvals * n / ranks
    # Enforce monotonicity (cumulative min from right)
    q_sorted = q[order]
    for i in range(n - 2, -1, -1):
        q_sorted[i] = min(q_sorted[i], q_sorted[i + 1])
    q[order] = q_sorted
    return np.minimum(q, 1.0)


def run_regression(cohort="cohortI"):
    print(f"[run_regression] loading carrier matrix for {cohort} …")
    mat = build_carrier_matrix(cohort)

    try:
        import statsmodels.formula.api as smf
    except ImportError:
        print("ERROR: statsmodels is required. Install with: pip install statsmodels")
        sys.exit(1)

    pc_cols = [f"PC{i}" for i in range(1, 11)]
    available_pcs = [c for c in pc_cols if c in mat.columns and mat[c].notna().any()]
    print(f"  PCs available: {available_pcs}")

    # Case groups: all Group_merged values that are not control or FamHX
    all_groups = mat["Group_merged"].dropna().unique()
    case_groups = sorted(set(all_groups) - {CONTROL_LABEL, FAMHX_LABEL})
    print(f"  Case groups: {len(case_groups)}")

    # Controls = samples with Group_merged == CONTROL_LABEL (all controls pooled)
    controls = mat[mat["Group_merged"] == CONTROL_LABEL].drop_duplicates(subset=["sample_id"])

    records = []

    for syn_key in SYNDROME_ORDER:
        syn_label = SYNDROME_LABELS[syn_key]
        for vc_label, col_suffix in [
            ("ACMG P/LP",              "ACMG_PLP"),
            ("ClinVar P/LP",           "ClinVar"),
            ("AM_calibrated",          "AM_calibrated"),
            ("AM_calibrated not P/LP", "AM_not_PLP"),
        ]:
            carrier_col = f"carrier_{syn_key}_{col_suffix}"
            if carrier_col not in mat.columns:
                continue

            for pheno in case_groups:
                # Cases: samples in this phenotype group
                cases = mat[mat["Group_merged"] == pheno].drop_duplicates(subset=["sample_id"])

                if len(cases) < N_MIN_PHENOTYPE:
                    continue

                n_carriers_cases = cases[carrier_col].sum()
                n_carriers_ctrl  = controls[carrier_col].sum()

                if n_carriers_cases < N_MIN_CARRIERS or n_carriers_ctrl < N_MIN_CARRIERS:
                    continue

                # Build regression DataFrame
                df_reg = pd.concat([
                    cases.assign(outcome=1),
                    controls.assign(outcome=0),
                ], ignore_index=True)

                df_reg = df_reg[["sample_id", "outcome", carrier_col,
                                  "Age_at_diagnosis", "Age2", "GENDER"]
                                 + available_pcs].dropna()
                df_reg = df_reg.rename(columns={carrier_col: "carrier"})
                df_reg["carrier"] = df_reg["carrier"].astype(int)

                pc_terms = " + ".join(available_pcs) if available_pcs else ""
                formula = ("outcome ~ carrier + Age_at_diagnosis + Age2 + C(GENDER)"
                           + (f" + {pc_terms}" if pc_terms else ""))

                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        model = smf.logit(formula, data=df_reg).fit(
                            cov_type="HC0", disp=False, maxiter=200
                        )
                    coef  = model.params["carrier"]
                    se    = model.bse["carrier"]
                    pval  = model.pvalues["carrier"]
                    OR    = float(np.exp(coef))
                    ci_lo = float(np.exp(coef - 1.96 * se))
                    ci_hi = float(np.exp(coef + 1.96 * se))
                    status = "ok"
                except Exception as exc:
                    OR, ci_lo, ci_hi, pval = np.nan, np.nan, np.nan, np.nan
                    status = f"failed: {exc}"

                records.append({
                    "syndrome":            syn_label,
                    "syndrome_key":        syn_key,
                    "phenotype":           pheno,
                    "variant_class":       vc_label,
                    "n_cases_total":       len(cases),
                    "n_controls_total":    len(controls),
                    "n_carriers_cases":    int(n_carriers_cases),
                    "n_carriers_controls": int(n_carriers_ctrl),
                    "OR":                  OR,
                    "CI_low":              ci_lo,
                    "CI_high":             ci_hi,
                    "p_value":             pval,
                    "status":              status,
                })

    results = pd.DataFrame(records)

    if results.empty:
        print("  No results produced — check carrier counts and data sources.")
        return results

    # BH correction across all valid tests
    valid = results["p_value"].notna()
    q = np.full(len(results), np.nan)
    q[valid] = _bh_correction(results.loc[valid, "p_value"].values)
    results["q_BH"] = q

    out_dir = RESULTS_I if cohort == "cohortI" else RESULTS_II
    out_path = os.path.join(out_dir, "regression_results.tsv")
    results.to_csv(out_path, sep="\t", index=False)
    print(f"  Saved: {out_path}  ({len(results)} rows)")
    return results


def main():
    parser = argparse.ArgumentParser(description="Run logistic regression for BioMe AM calibration.")
    parser.add_argument("--cohort", default="cohortI", choices=["cohortI", "cohortII"])
    args = parser.parse_args()
    run_regression(args.cohort)


if __name__ == "__main__":
    main()
