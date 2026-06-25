"""Call carriers using the Chen/Pejaver per-variant evidence label.

This is the **lookup-based** approach (Chen recommended): for each annotated
(sample, variant) row in the input, we look at the Chen `evidence` label
that step 02 attached, and:

- `is_AM_carrier_primary = TRUE` iff the label is in the configured
  positive set (default: PP3_Moderate, PP3_Moderate+, PP3_Strong, PP3_Strong+,
  PP3_VeryStrong).
- `is_AM_carrier_global_0864 = TRUE` iff am_pathogenicity >=
  calibration.global_threshold_legacy (default 0.864). This is the
  counterfactual the rebuttal letter relies on — what the original global
  threshold would have called, regardless of Chen calibration. Independent
  of `is_AM_carrier_primary`.
- `would_be_carrier_domain_aggregate = TRUE` is informational — TRUE iff
  `is_AM_carrier_primary` AND `chen_calibration_approach == domain_aggregate`.
- `threshold_source` is one of:
    single_gene             — variant in Chen with calibration_approach=single_gene
    domain_aggregate        — variant in Chen with calibration_approach=domain_aggregate
    not_in_chen_table       — variant not in Chen (primary=FALSE; AM score retained)

We do NOT derive per-gene thresholds anymore. Per Chen's framework, each
variant's evidence label already reflects the gene's (or its PFAM domain's)
calibrated posterior — re-deriving thresholds would duplicate work the table
encodes per-variant.

Input : intermediate/<cohort>/am_annotated.tsv
        (one row per (sample, variant); produced by python/annotate_am.py)
Output: intermediate/<cohort>/A_variant_level_per_person.tsv
        (super-set of the input + carrier-call columns; `passes_QC` is TRUE
         for every row because step 01 already masked failing GTs and step 02
         only emitted non-missing / non-homref samples).
"""
from __future__ import annotations

import argparse
import os
import sys
from typing import List, Optional, Set

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from util import LOG, die, load_config, read_tsv_dicts, resolve, write_tsv  # noqa: E402


def _normalize_evidence(label: str) -> str:
    if not label:
        return ""
    return label.strip().lower().replace(" ", "").replace("-", "_")


_PP3_AT_LEAST_MODERATE = {
    "pp3_moderate", "pp3_moderate+",
    "pp3_strong",   "pp3_strong+",
    "pp3_verystrong", "pp3_very_strong",
}
_PP3_AT_LEAST_SUPPORTING = _PP3_AT_LEAST_MODERATE | {"pp3_supporting", "pp3_supporting+"}
_PP3_AT_LEAST_STRONG = {"pp3_strong", "pp3_strong+", "pp3_verystrong", "pp3_very_strong"}
_PP3_AT_LEAST_VERYSTRONG = {"pp3_verystrong", "pp3_very_strong"}


def evidence_set_for(min_strength: str) -> Set[str]:
    s = (min_strength or "").strip().lower().replace(" ", "")
    if s in ("supporting", "pp3supporting", "pp3_supporting"):
        return _PP3_AT_LEAST_SUPPORTING
    if s in ("strong", "pp3strong", "pp3_strong"):
        return _PP3_AT_LEAST_STRONG
    if s in ("verystrong", "very_strong", "pp3verystrong", "pp3_verystrong"):
        return _PP3_AT_LEAST_VERYSTRONG
    return _PP3_AT_LEAST_MODERATE  # default


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--in", dest="inp", required=True, help="am_annotated.tsv")
    ap.add_argument("--out", required=True, help="A_variant_level_per_person.tsv")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    min_strength = cfg.get("calibration", {}).get("min_evidence_strength", "Moderate")
    positive_set = evidence_set_for(min_strength)
    global_legacy = float(cfg.get("calibration", {}).get("global_threshold_legacy", 0.864))

    # The locked spec says primary_threshold MUST be the published Chen label —
    # we honor `primary_threshold` if present for symmetry, but the only
    # legal value now is "gene_specific" (meaning: use Chen's published per-
    # variant label, derived from the gene's calibration where available and
    # the domain-aggregate calibration otherwise).
    policy = cfg.get("calibration", {}).get("primary_threshold", "gene_specific")
    if policy != "gene_specific":
        die(f"calibration.primary_threshold must be 'gene_specific' (got {policy!r}); "
            f"the spec only supports the Chen per-variant lookup.")

    id_col = cfg["cohorts"][args.cohort]["id_column"]

    _, rows = read_tsv_dicts(args.inp)
    out_header = [
        "sample_id", "id_column_used", "cohort", "chr", "pos", "ref", "alt",
        "gene", "transcript", "protein_variant",
        "am_pathogenicity",
        "chen_evidence", "chen_points", "chen_calibration_approach",
        "chen_domain", "chen_vep_score", "in_chen_table",
        "threshold_source", "genotype", "DP", "GQ", "passes_QC",
        "is_AM_carrier_primary", "is_AM_carrier_global_0864",
        "would_be_carrier_domain_aggregate",
    ]
    out_rows: List[List[str]] = []
    n_primary_single = 0
    n_primary_domain = 0
    n_subthreshold = 0
    n_not_in_chen = 0
    n_global_legacy = 0
    for r in rows:
        in_chen = (r.get("in_chen_table", "") or "").strip().upper() == "TRUE"
        ev = _normalize_evidence(r.get("chen_evidence", ""))
        approach = (r.get("chen_calibration_approach", "") or "").strip().lower()
        if not in_chen or not ev:
            source = "not_in_chen_table"; n_not_in_chen += 1
            is_primary = False
            would_da = False
        elif ev in positive_set:
            is_primary = True
            if approach.startswith("single") or approach == "gene_specific":
                source = "single_gene"; n_primary_single += 1
                would_da = False
            else:
                source = "domain_aggregate"; n_primary_domain += 1
                would_da = True
        else:
            # variant is in Chen but its label is below the min strength
            # (e.g. PP3_Supporting / BP4_* / NO_EVIDENCE)
            if approach.startswith("single") or approach == "gene_specific":
                source = "single_gene"
            elif approach:
                source = "domain_aggregate"
            else:
                source = "not_in_chen_table"
            n_subthreshold += 1
            is_primary = False
            would_da = False
        # Counterfactual under the original global threshold. Independent of
        # Chen membership — variants not in Chen still get this flag.
        amp_raw = (r.get("am_pathogenicity", "") or "").strip()
        try:
            is_global = float(amp_raw) >= global_legacy
        except ValueError:
            is_global = False
        if is_global:
            n_global_legacy += 1
        out_rows.append([
            r.get("sample_id", ""), id_col, args.cohort,
            r.get("chr", ""), r.get("pos", ""), r.get("ref", ""), r.get("alt", ""),
            r.get("gene", ""), r.get("transcript", ""), r.get("protein_variant", ""),
            r.get("am_pathogenicity", ""),
            r.get("chen_evidence", ""), r.get("chen_points", ""),
            r.get("chen_calibration_approach", ""),
            r.get("chen_domain", ""), r.get("chen_vep_score", ""),
            r.get("in_chen_table", "FALSE"),
            source,
            r.get("genotype", ""), r.get("DP", ""), r.get("GQ", ""),
            "TRUE",
            "TRUE" if is_primary else "FALSE",
            "TRUE" if is_global else "FALSE",
            "TRUE" if would_da else "FALSE",
        ])

    write_tsv(args.out, out_header, out_rows)
    LOG.info("%s call_carriers: rows=%d  primary_single=%d  primary_domain=%d  "
             "subthreshold=%d  not_in_chen=%d  global_legacy(>=%.3f)=%d  (min_strength=%s)",
             args.cohort, len(out_rows),
             n_primary_single, n_primary_domain, n_subthreshold, n_not_in_chen,
             global_legacy, n_global_legacy, min_strength)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
