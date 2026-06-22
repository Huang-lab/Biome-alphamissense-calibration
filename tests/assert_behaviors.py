"""Behavioral assertions for the synthetic test run.

8 cases (spec 7 + AB-filter syntax), updated for the Chen per-variant lookup:
  1. BRCA1 (S1@chr1:150) — Chen PP3_Moderate / single_gene → primary=TRUE,
     threshold_source=single_gene
  2. MEN1  (S1@chr2:550) — Chen PP3_Moderate+ / domain_aggregate → primary=TRUE
     (per Chen, the published label IS the carrier evidence regardless of
     calibration scope), threshold_source=domain_aggregate,
     would_be_carrier_domain_aggregate=TRUE
  3. WT1   (S1@chr3:950) — not in Chen table → threshold_source=not_in_chen_table,
     primary=FALSE, am_pathogenicity still recorded
  4. AM_only — S2 at BRCA1 chr1:150, category=AM_only in Table B
  5. Both    — S1 at BRCA1 chr1:150, category=both in Table B
  6. QC fail — S4 (DP=5 het) absent from Table A at chr1:150
  7. gVCF guard — checked in run_test.sh
  8. AB filter — chr1:200 S1 (AD=20,20) present; chr1:200 S3 (AD=19,1) absent
"""
from __future__ import annotations

import argparse
import csv
import os
import sys


def load_tsv(path: str):
    if not os.path.isfile(path):
        sys.exit(f"asserts: file not found: {path}")
    with open(path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        return list(rdr)


def find_row(rows, **kw):
    for r in rows:
        if all(str(r.get(k, "")).strip() == str(v).strip() for k, v in kw.items()):
            return r
    return None


def assert_case(label: str, cond: bool, extra: str = "") -> None:
    if cond:
        print(f"  [OK] {label}")
    else:
        print(f"  [FAIL] {label}  {extra}")
        sys.exit(1)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--config", required=True)
    args = ap.parse_args(argv)

    # Cohort I A and B tables
    a_path = os.path.join(args.run_dir, "results", "cohortI", "A_variant_level_per_person.tsv")
    b_path = os.path.join(args.run_dir, "results", "cohortI", "B_ACMG_vs_AM_comparison.tsv")
    c_path = os.path.join(args.run_dir, "results", "cohortI", "C_regression_matrix.tsv")
    a_rows = load_tsv(a_path)
    b_rows = load_tsv(b_path)
    c_rows = load_tsv(c_path)

    print(f"\nTable A rows: {len(a_rows)}")
    print(f"Table B rows: {len(b_rows)}")
    print(f"Table C rows: {len(c_rows)}\n")

    # Case 1: BRCA1 single_gene carrier (PP3_Moderate)
    r = find_row(a_rows, sample_id="S1", chr="chr1", pos="150", gene="BRCA1")
    assert_case("case 1: BRCA1 single_gene carrier (S1@chr1:150) present",
                r is not None, "no Table A row for S1@chr1:150/BRCA1")
    assert_case("case 1: threshold_source=single_gene",
                r is not None and r["threshold_source"] == "single_gene",
                f"got threshold_source={r and r.get('threshold_source')!r}")
    assert_case("case 1: is_AM_carrier_primary=TRUE",
                r is not None and r["is_AM_carrier_primary"].upper() == "TRUE")
    assert_case("case 1: chen_evidence is PP3_Moderate",
                r is not None and r.get("chen_evidence", "") == "PP3_Moderate",
                f"got chen_evidence={r and r.get('chen_evidence')!r}")

    # Case 2: MEN1 domain_aggregate carrier (PP3_Moderate+; primary=TRUE per
    # Chen — the published label is the carrier evidence regardless of scope)
    r = find_row(a_rows, sample_id="S1", chr="chr2", pos="550", gene="MEN1")
    assert_case("case 2: MEN1 domain_aggregate row (S1@chr2:550) present", r is not None)
    assert_case("case 2: threshold_source=domain_aggregate",
                r is not None and r["threshold_source"] == "domain_aggregate")
    assert_case("case 2: is_AM_carrier_primary=TRUE (Chen PP3_Moderate+ counts)",
                r is not None and r["is_AM_carrier_primary"].upper() == "TRUE",
                f"got is_AM_carrier_primary={r and r.get('is_AM_carrier_primary')!r}")
    assert_case("case 2: would_be_carrier_domain_aggregate=TRUE",
                r is not None and r["would_be_carrier_domain_aggregate"].upper() == "TRUE")

    # Case 3: WT1 not in Chen table → not_in_chen_table, primary=FALSE
    r = find_row(a_rows, sample_id="S1", chr="chr3", pos="950", gene="WT1")
    assert_case("case 3: WT1 not_in_chen row (S1@chr3:950) present", r is not None)
    assert_case("case 3: threshold_source=not_in_chen_table",
                r is not None and r["threshold_source"] == "not_in_chen_table",
                f"got threshold_source={r and r.get('threshold_source')!r}")
    assert_case("case 3: is_AM_carrier_primary=FALSE",
                r is not None and r["is_AM_carrier_primary"].upper() == "FALSE")
    am = r and r.get("am_pathogenicity", "")
    assert_case("case 3: am_pathogenicity still present (non-empty)",
                bool(am), f"am_pathogenicity={am!r}")

    # Case 4: AM_only — S2 at BRCA1
    r = find_row(b_rows, sample_id="S2", chr="chr1", pos="150", gene="BRCA1")
    assert_case("case 4: B row for S2@chr1:150/BRCA1 present", r is not None)
    assert_case("case 4: category=AM_only",
                r is not None and r["category"] == "AM_only",
                f"got category={r and r.get('category')!r}")

    # Case 5: both — S1 at BRCA1
    r = find_row(b_rows, sample_id="S1", chr="chr1", pos="150", gene="BRCA1")
    assert_case("case 5: B row for S1@chr1:150/BRCA1 present", r is not None)
    assert_case("case 5: category=both",
                r is not None and r["category"] == "both",
                f"got category={r and r.get('category')!r}")

    # Case 6: QC fail — S4 absent from Table A at chr1:150 (low DP=5)
    r = find_row(a_rows, sample_id="S4", chr="chr1", pos="150")
    assert_case("case 6: S4 (low DP=5 het) absent from Table A at chr1:150",
                r is None, "S4 should have been masked to ./. by QC and dropped")

    # Case 7: handled in run_test.sh (gVCF guard test before this script runs)

    # Case 8: AB filter — S1 (AD=20,20) kept; S3 (AD=19,1) dropped, both at chr1:200
    r_s1 = find_row(a_rows, sample_id="S1", chr="chr1", pos="200")
    r_s3 = find_row(a_rows, sample_id="S3", chr="chr1", pos="200")
    assert_case("case 8: balanced het S1@chr1:200 (AD=20,20) PRESENT in Table A",
                r_s1 is not None,
                "the AB filter expression may be syntactically wrong for this bcftools version — see lib/load_modules.sh banner")
    assert_case("case 8: imbalanced het S3@chr1:200 (AD=19,1) ABSENT from Table A",
                r_s3 is None,
                "the AB filter expression may be syntactically wrong for this bcftools version — see lib/load_modules.sh banner")

    # Bonus: Cohort II (Sema4) mode produced output via sliced_from_combined
    a2 = os.path.join(args.run_dir, "results", "cohortII", "A_variant_level_per_person.tsv")
    a2_rows = load_tsv(a2)
    r2 = find_row(a2_rows, sample_id="M1", chr="chr1", pos="150", gene="BRCA1")
    assert_case("bonus: Cohort II (sliced_from_combined) produced AM-primary row for M1@chr1:150/BRCA1",
                r2 is not None and r2["is_AM_carrier_primary"].upper() == "TRUE")

    print("\nALL ASSERTIONS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
