"""Behavioral assertions for the synthetic test run.

10 cases (spec 7 + AB-filter syntax + global-0864 counterfactual + Table C
carrier columns), updated for the Chen per-variant lookup:
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
  9. is_AM_carrier_global_0864 column present in Table A; TRUE iff
     am_pathogenicity >= 0.864 regardless of Chen membership.
 10. carrier_<group>_AM0864 column present in Table C; reflects the per-group
     OR of the global-legacy carrier flag across the group's genes.
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

    # Case 4: AM-only — S2 at BRCA1 (AM primary, not ACMG, not ClinVar)
    r = find_row(b_rows, sample_id="S2", chr="chr1", pos="150", gene="BRCA1")
    assert_case("case 4: B row for S2@chr1:150/BRCA1 present", r is not None)
    assert_case("case 4: category=AM (AM-only)",
                r is not None and r["category"] == "AM",
                f"got category={r and r.get('category')!r}")
    assert_case("case 4: in_ClinVar_PLP=no",
                r is not None and r["in_ClinVar_PLP"] == "no")

    # Case 5: ACMG+AM — S1 at BRCA1 chr1:150 (conflicting ClinVar excluded -> no)
    r = find_row(b_rows, sample_id="S1", chr="chr1", pos="150", gene="BRCA1")
    assert_case("case 5: B row for S1@chr1:150/BRCA1 present", r is not None)
    assert_case("case 5: category=ACMG_PLP+AM",
                r is not None and r["category"] == "ACMG_PLP+AM",
                f"got category={r and r.get('category')!r}")
    assert_case("case 5: in_ClinVar_PLP=no (chr1:150 was Conflicting -> dropped)",
                r is not None and r["in_ClinVar_PLP"] == "no",
                f"got {r and r.get('in_ClinVar_PLP')!r}")

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

    # Case 9: is_AM_carrier_global_0864 counterfactual column.
    # Fixture mini_AM has BRCA1@chr1:150 at am=0.900 (>=0.864) and BRCA1@chr1:200
    # at am=0.700 (<0.864). MEN1@chr2:550 at am=0.500 (<0.864) — so even
    # though it's a Chen primary carrier, the global-legacy flag is FALSE.
    r_brca = find_row(a_rows, sample_id="S1", chr="chr1", pos="150", gene="BRCA1")
    r_brca_low = find_row(a_rows, sample_id="S1", chr="chr1", pos="200", gene="BRCA1")
    r_men1 = find_row(a_rows, sample_id="S1", chr="chr2", pos="550", gene="MEN1")
    assert_case("case 9: is_AM_carrier_global_0864 column present in Table A",
                r_brca is not None and "is_AM_carrier_global_0864" in r_brca)
    assert_case("case 9: BRCA1 S1@chr1:150 (am=0.900) -> is_AM_carrier_global_0864=TRUE",
                r_brca is not None and r_brca["is_AM_carrier_global_0864"].upper() == "TRUE",
                f"got {r_brca and r_brca.get('is_AM_carrier_global_0864')!r}")
    assert_case("case 9: BRCA1 S1@chr1:200 (am=0.700) -> is_AM_carrier_global_0864=FALSE",
                r_brca_low is not None and r_brca_low["is_AM_carrier_global_0864"].upper() == "FALSE")
    assert_case("case 9: MEN1 S1@chr2:550 (am=0.500, primary=TRUE) -> global_0864=FALSE "
                "(primary and global are independent counterfactuals)",
                r_men1 is not None
                and r_men1["is_AM_carrier_primary"].upper() == "TRUE"
                and r_men1["is_AM_carrier_global_0864"].upper() == "FALSE")

    # Case 10: per-group AM0864 carrier columns in Table C.
    c_brca_row = find_row(c_rows, sample_id="S1")
    assert_case("case 10: Table C has carrier_<group>_AM0864 columns",
                c_brca_row is not None
                and "carrier_Hereditary_Breast_and_Ovarian_Cancer_Syndrome_AM0864" in c_brca_row)
    assert_case("case 10: S1 has carrier_HBOC_AM0864=TRUE (BRCA1 am=0.900)",
                c_brca_row is not None
                and c_brca_row["carrier_Hereditary_Breast_and_Ovarian_Cancer_Syndrome_AM0864"].upper() == "TRUE")
    assert_case("case 10: S1 has carrier_MEN_AM0864=FALSE (MEN1 am=0.500 < 0.864)",
                c_brca_row is not None
                and c_brca_row["carrier_Multiple_Endocrine_Neoplasia_AM0864"].upper() == "FALSE")

    # Case 11: ClinVar in the 3-way Table B carrier list. Fixture mini_clinvar:
    #   - MEN1  chr2:550 C>T  Pathogenic, 2*                        -> KEPT
    #   - WT1   chr3:950 G>A  Likely_pathogenic, 1*                 -> DROPPED (<2*)
    #   - BRCA1 chr1:300 A>G  Pathogenic, 3*, nonsense (PTV)        -> KEPT (exclude_ptv=false) but no carrier at chr1:300
    #   - BRCA1 chr1:150 A>G  Conflicting                          -> DROPPED (exclude_conflicting)
    #   - MLH1  chr1:180 A>AT Pathogenic, 2*, frameshift indel     -> KEPT (non-SNV)
    # S1 carries MEN1 chr2:550 (AM primary + ClinVar).
    assert_case("case 11: Table B has in_ClinVar_PLP column",
                len(b_rows) > 0 and "in_ClinVar_PLP" in b_rows[0])
    r = find_row(b_rows, sample_id="S1", chr="chr2", pos="550", gene="MEN1")
    assert_case("case 11: B row for S1@chr2:550/MEN1 present", r is not None)
    assert_case("case 11: MEN1 chr2:550 in_ClinVar_PLP=yes (Pathogenic, 2*)",
                r is not None and r["in_ClinVar_PLP"] == "yes",
                f"got {r and r.get('in_ClinVar_PLP')!r}")
    assert_case("case 11: MEN1 chr2:550 category=AM+ClinVar (AM primary + ClinVar, not ACMG)",
                r is not None and r["category"] == "AM+ClinVar",
                f"got category={r and r.get('category')!r}")
    assert_case("case 11: MEN1 chr2:550 clinvar_review_stars=2",
                r is not None and str(r.get("clinvar_review_stars", "")).strip() == "2")
    # WT1 chr3:950 was only 1* -> must NOT appear as a ClinVar carrier row
    r_wt1 = find_row(b_rows, sample_id="S1", chr="chr3", pos="950", gene="WT1")
    assert_case("case 11: WT1 chr3:950 (1*) not a ClinVar carrier in Table B",
                r_wt1 is None or r_wt1.get("in_ClinVar_PLP") != "yes",
                f"got {r_wt1 and r_wt1.get('in_ClinVar_PLP')!r}")

    # Case 12: non-SNV (indel) + PTV-included ClinVar capture via the ALL-VARIANT
    # QC pass. chr1:180 A>AT (frameshift) Pathogenic 2*, MLH1, carried by S1.
    # The SNV-only AM path drops it; 01c (all-variant) keeps it, so it appears in
    # Table B as a ClinVar carrier row. exclude_ptv=false keeps the frameshift.
    r_indel_b = find_row(b_rows, sample_id="S1", chr="chr1", pos="180", gene="MLH1")
    assert_case("case 12: indel chr1:180 A>AT present in Table B as ClinVar carrier (non-SNV kept)",
                r_indel_b is not None and r_indel_b["in_ClinVar_PLP"] == "yes",
                f"got row={r_indel_b!r}")
    assert_case("case 12: indel chr1:180 ref/alt preserved (A/AT)",
                r_indel_b is not None and r_indel_b.get("ref") == "A" and r_indel_b.get("alt") == "AT",
                f"got ref/alt={r_indel_b and (r_indel_b.get('ref'), r_indel_b.get('alt'))!r}")
    assert_case("case 12: indel chr1:180 in_AM_primary=no (SNV-only AM path never saw it)",
                r_indel_b is not None and r_indel_b["in_AM_primary"] == "no")
    # And confirm the indel really is absent from the SNV-only Table A (AM path).
    r_indel_a = find_row(a_rows, sample_id="S1", chr="chr1", pos="180")
    assert_case("case 12: indel chr1:180 ABSENT from SNV-only Table A (AM path unaffected)",
                r_indel_a is None, "SNV-only QC should have dropped the indel from the AM path")

    # Case 13: ACMG P/LP variants are now shown WITH coordinates (variant-level),
    # not gene-level blanks. S6 carries BRCA1 chr1:300_A_G (ACMG P/LP) but is not
    # a VCF/AM carrier, so it enters Table B via the ACMG variant column.
    r_acmg = find_row(b_rows, sample_id="S6", gene="BRCA1", pos="300")
    assert_case("case 13: ACMG-only S6 BRCA1 chr1:300 present in Table B with coords",
                r_acmg is not None, "ACMG variant coords should be parsed from the variant column")
    assert_case("case 13: S6 chr1:300 in_ACMG_PLP=yes, in_AM_primary=no, category=ACMG_PLP",
                r_acmg is not None and r_acmg["in_ACMG_PLP"] == "yes"
                and r_acmg["in_AM_primary"] == "no" and r_acmg["category"] == "ACMG_PLP",
                f"got {r_acmg!r}")
    assert_case("case 13: S6 chr1:300 ref/alt populated (A/G), not blank",
                r_acmg is not None and r_acmg.get("ref") == "A" and r_acmg.get("alt") == "G",
                f"got ref/alt={r_acmg and (r_acmg.get('ref'), r_acmg.get('alt'))!r}")

    print("\nALL ASSERTIONS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
