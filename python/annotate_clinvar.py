"""Standalone ClinVar P/LP (≥2★) carrier pass.

This does NOT touch the AM/Chen annotation path (annotate_am.py /
call_carriers.py). It reads the SAME step-01 QC'd per-chr VCF and, for each
carrier-like genotype, looks the variant up in the ClinVar P/LP ≥2★ subset
(refs/clinvar_plp_2star.target_genes.tsv.gz, built by prepare_refs.subset_clinvar).
A carrier row is emitted only when the variant is in that subset.

The resulting list is joined downstream by compare_tabulate.py exactly like the
ACMG P/LP list, giving a second clinical-truth comparator (`ClinVar_PLP`) and
redefining `AM_calibrated_not_PLP` = AM-calibrated AND NOT ACMG P/LP AND NOT
ClinVar P/LP.

Input : ALL-VARIANT QC'd per-chr VCF (intermediate/<cohort>/chr<N>.qc_allvar.vcf.gz)
        — same QC as step 01 but WITHOUT the SNV-only filter, so indels/MNVs are
        retained. Both this VCF and the ClinVar subset are left-normalized against
        the reference FASTA so non-SNV variants join reliably by (chrom,pos,ref,alt).
Output: intermediate/<cohort>/chr<N>.clinvar_carriers.tsv (per (sample, variant))
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from typing import Dict, List, Optional, Tuple

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from util import LOG, chrom_norm, die, load_config, min_rep, resolve, vcf_chrom_prefix, write_tsv  # noqa: E402
from annotate_am import stream_vcf_genotypes  # reuse the AM VCF reader  # noqa: E402

CvKey = Tuple[str, int, str, str]  # (chrom, pos, ref, alt)


def load_clinvar_for_region(clinvar_tbi: str, region: str) -> Dict[CvKey, dict]:
    """tabix the ClinVar subset for `region` (e.g. 'chr1'), index by
    (chrom,pos,ref,alt). Empty dict if the subset is missing (then no ClinVar
    carriers are called — loudly warned).

    ClinVar subset header: chrom pos ref alt gene_symbol clnsig review_stars mc
    """
    out: Dict[CvKey, dict] = {}
    if not os.path.isfile(clinvar_tbi):
        LOG.warning("ClinVar subset not found: %s (no ClinVar carriers will be called)", clinvar_tbi)
        return out
    if not os.path.isfile(clinvar_tbi + ".tbi"):
        LOG.warning("ClinVar subset tabix index missing: %s.tbi (no ClinVar carriers)", clinvar_tbi)
        return out
    cmd = ["tabix", clinvar_tbi, region]
    LOG.info("tabix %s %s", clinvar_tbi, region)
    p = subprocess.run(cmd, check=True, capture_output=True, text=True)
    for line in p.stdout.splitlines():
        f = line.rstrip("\n").split("\t")
        if len(f) < 7:
            continue
        chrom = chrom_norm(f[0])
        try:
            pos = int(f[1])
        except ValueError:
            continue
        ref, alt = f[2], f[3]
        # canonicalize indel padding so the key matches the VCF's representation
        pos, ref, alt = min_rep(pos, ref, alt)
        out[(chrom, pos, ref, alt)] = {
            "gene_symbol":  f[4],
            "clnsig":       f[5],
            "review_stars": f[6],
        }
    LOG.info("loaded %d ClinVar P/LP≥2★ rows for region %s", len(out), region)
    return out


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Flag ClinVar P/LP ≥2★ carriers from a QC'd VCF.")
    ap.add_argument("--config", default=None)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--chr", required=True, help="e.g. chr1")
    ap.add_argument("--vcf", required=True, help="QC'd per-chr VCF.gz")
    ap.add_argument("--out", required=True, help="output TSV (per-(sample,variant))")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    refs_dir = resolve(cfg, cfg["paths"]["refs_dir"])
    clinvar_subset = os.path.join(refs_dir, "clinvar_plp_2star.target_genes.tsv.gz")

    # The ClinVar subset is always chr-prefixed (chrom_norm output). The QC'd
    # VCF may use either convention; detect per-VCF so `bcftools query -r` matches.
    ref_region = args.chr if args.chr.startswith("chr") else f"chr{args.chr}"
    vcf_pfx = vcf_chrom_prefix(args.vcf)
    chr_num = args.chr[3:] if args.chr.startswith("chr") else args.chr
    vcf_region = f"{vcf_pfx}{chr_num}"
    LOG.info("region: ref=%s vcf=%s (vcf prefix=%s)", ref_region, vcf_region, vcf_pfx or "<none>")

    clinvar = load_clinvar_for_region(clinvar_subset, ref_region)

    header = ["sample_id", "cohort", "chr", "pos", "ref", "alt", "gene",
              "clinvar_clnsig", "clinvar_review_stars", "genotype", "DP", "GQ"]
    rows_out: List[List[str]] = []
    n_var_in = 0
    n_var_matched = 0
    n_carrier_rows = 0
    for chrom, pos, ref, alt, samples in stream_vcf_genotypes(args.vcf, region=vcf_region):
        n_var_in += 1
        # canonicalize the VCF variant the same way so indel padding differences
        # don't defeat the (chrom,pos,ref,alt) join with the ClinVar subset.
        mpos, mref, malt = min_rep(pos, ref, alt)
        rec = clinvar.get((chrom, mpos, mref, malt))
        if rec is None:
            continue
        n_var_matched += 1
        for sample, gt, dp, gq in samples:
            rows_out.append([
                sample, args.cohort, chrom, str(mpos), mref, malt,
                rec["gene_symbol"], rec["clnsig"], rec["review_stars"],
                gt, dp, gq,
            ])
            n_carrier_rows += 1

    write_tsv(args.out, header, rows_out)
    LOG.info("%s %s clinvar: variants_streamed=%d ClinVar_matched=%d carrier_rows=%d",
             args.cohort, ref_region, n_var_in, n_var_matched, n_carrier_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
