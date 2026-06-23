"""Annotate per-sample carrier-like genotypes against the AM subset AND the
Chen/Pejaver calibration table (variant-level lookup).

Input : QC'd per-chr VCF (failing GTs already masked to ./. by step 01)
Output: per (sample, variant) row for every AM-matched site where the sample
        has a non-missing, non-hom-ref genotype. Each row carries:
          - the AM score / class / transcript (from AlphaMissense_hg38)
          - the Chen calibration label fields (from chen_calibration.target_genes.tsv.gz):
              chen_evidence (e.g. PP3_Moderate, PP3_Moderate+, BP4_Supporting, ...),
              chen_points,
              chen_calibration_approach (single_gene | domain_aggregate),
              chen_domain (PFAM domain name, or no_pfam)
          - in_chen_table (TRUE/FALSE)

Joining is by (CHROM, POS, REF, ALT) for BOTH tables; we read the per-chr
slice via tabix and build small in-memory dicts (target-gene subsets fit
comfortably). We also assert that every AM-matched record is a single-base
SNV (Chen rows are also SNVs by construction).
"""
from __future__ import annotations

import argparse
import gzip
import os
import subprocess
import sys
from typing import Dict, Iterable, List, Optional, Tuple

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from util import LOG, chrom_norm, die, load_config, resolve, vcf_chrom_prefix, write_tsv  # noqa: E402


# ----------------------------------------------------------------------------
# AM lookup
# ----------------------------------------------------------------------------
AmKey = Tuple[str, int, str, str]  # (chrom, pos, ref, alt)


def load_am_for_region(am_tbi: str, region: str) -> Dict[AmKey, dict]:
    """tabix the AM subset for `region` (e.g. 'chr1') and index by (chrom,pos,ref,alt).

    AM subset header: chrom pos ref alt uniprot_id transcript_id protein_variant am_pathogenicity am_class gene_name
    """
    if not os.path.isfile(am_tbi):
        die(f"AM subset not found: {am_tbi}")
    if not os.path.isfile(am_tbi + ".tbi"):
        die(f"AM subset tabix index missing: {am_tbi}.tbi (run scripts/00_prepare_refs.sh)")
    out: Dict[AmKey, dict] = {}
    cmd = ["tabix", am_tbi, region]
    LOG.info("tabix %s %s", am_tbi, region)
    p = subprocess.run(cmd, check=True, capture_output=True, text=True)
    for line in p.stdout.splitlines():
        f = line.rstrip("\n").split("\t")
        if len(f) < 10:
            continue
        chrom = chrom_norm(f[0])
        try:
            pos = int(f[1])
        except ValueError:
            continue
        ref, alt = f[2], f[3]
        if len(ref) != 1 or len(alt) != 1:
            # AM-matched non-SNV → hard fail per locked decision C
            die(f"AM-matched non-SNV at {chrom}:{pos} {ref}/{alt}; AM should be missense-SNV only")
        out[(chrom, pos, ref, alt)] = {
            "uniprot_id": f[4],
            "transcript_id": f[5],
            "protein_variant": f[6],
            "am_pathogenicity": f[7],
            "am_class": f[8],
            "gene_name": f[9],
        }
    LOG.info("loaded %d AM rows for region %s", len(out), region)
    return out


def load_chen_for_region(chen_tbi: str, region: str) -> Dict[AmKey, dict]:
    """tabix the Chen calibration subset for `region` (e.g. 'chr1') and index
    by (chrom,pos,ref,alt). Returns an empty dict if the file doesn't exist —
    step 03 then sees `in_chen_table=FALSE` and never calls a carrier (loudly
    logged in REFERENCE_REPORT.md).

    Chen subset header:
      chrom pos ref alt gene_symbol evidence points calibration_approach domain vep_score
    """
    out: Dict[AmKey, dict] = {}
    if not os.path.isfile(chen_tbi):
        LOG.warning("Chen calibration subset not found: %s (no carriers will be called)", chen_tbi)
        return out
    if not os.path.isfile(chen_tbi + ".tbi"):
        LOG.warning("Chen calibration tabix index missing: %s.tbi (no carriers will be called)", chen_tbi)
        return out
    cmd = ["tabix", chen_tbi, region]
    LOG.info("tabix %s %s", chen_tbi, region)
    p = subprocess.run(cmd, check=True, capture_output=True, text=True)
    for line in p.stdout.splitlines():
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        chrom = chrom_norm(f[0])
        try:
            pos = int(f[1])
        except ValueError:
            continue
        ref, alt = f[2], f[3]
        key = (chrom, pos, ref, alt)
        # If the same (chrom,pos,ref,alt) appears twice across transcripts,
        # keep the strongest evidence label.
        rec = {
            "gene_symbol":          f[4],
            "evidence":             f[5],
            "points":               f[6],
            "calibration_approach": f[7],
            "domain":               f[8] if len(f) > 8 else "",
            "vep_score":            f[9] if len(f) > 9 else "",
        }
        prev = out.get(key)
        if prev is None or _evidence_rank(rec["evidence"]) > _evidence_rank(prev["evidence"]):
            out[key] = rec
    LOG.info("loaded %d Chen rows for region %s", len(out), region)
    return out


# Ranks used to break ties when the same variant appears twice in Chen (rare,
# but defensive). Mirrors prepare_refs._EVIDENCE_RANK; small inline copy avoids
# cross-module coupling.
_EV_RANK = {
    "pp3_supporting": 1, "pp3_moderate": 2, "pp3_strong": 3, "pp3_verystrong": 4,
}
def _evidence_rank(label: str) -> int:
    if not label:
        return 0
    k = label.strip().lower().replace(" ", "").replace("-", "_").rstrip("+")
    if k.startswith("bp4"):
        return 0
    return _EV_RANK.get(k, 0)


# ----------------------------------------------------------------------------
# VCF streaming via bcftools query
# ----------------------------------------------------------------------------
def stream_vcf_genotypes(vcf: str, region: Optional[str] = None) -> Iterable[Tuple[str, int, str, str, List[Tuple[str, str, str, str]]]]:
    """Yield (chrom, pos, ref, alt, [(sample, gt, dp, gq), ...]) per variant.

    Implementation: bcftools query with per-sample format expansion. We skip
    missing and hom-ref genotypes here (carrier-only emission). DP/GQ are
    carried through for the output table (genotypes already masked in step 01
    can still appear with low DP/GQ in legacy files — we DON'T re-filter on
    them since step 01 owns QC, but we record them).
    """
    fmt = '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT;%DP;%GQ]\n'
    cmd = ["bcftools", "query", "-f", fmt]
    if region:
        cmd += ["-r", region]
    cmd += [vcf]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    assert proc.stdout is not None
    try:
        for line in proc.stdout:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom = chrom_norm(parts[0])
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            ref, alt = parts[2], parts[3]
            samples: List[Tuple[str, str, str, str]] = []
            for cell in parts[4:]:
                # "S1=0/1;30;40"
                if "=" not in cell:
                    continue
                sample, rest = cell.split("=", 1)
                # rest = "GT;DP;GQ"
                bits = rest.split(";")
                gt = bits[0] if len(bits) > 0 else "."
                dp = bits[1] if len(bits) > 1 else "."
                gq = bits[2] if len(bits) > 2 else "."
                if gt in (".", "./.", ".|.", ""):
                    continue
                # drop hom-ref
                norm = gt.replace("|", "/")
                if norm in ("0/0",):
                    continue
                samples.append((sample, gt, dp, gq))
            if samples:
                yield chrom, pos, ref, alt, samples
    finally:
        rc = proc.wait()
        if rc != 0:
            die(f"bcftools query failed (rc={rc}) for {vcf}")


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Join per-chr filtered VCF to the AM subset.")
    ap.add_argument("--config", default=None)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--chr", required=True, help="e.g. chr1")
    ap.add_argument("--vcf", required=True, help="QC'd per-chr VCF.gz")
    ap.add_argument("--out", required=True, help="output TSV (per-(sample,variant))")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    refs_dir = resolve(cfg, cfg["paths"]["refs_dir"])
    am_subset = os.path.join(refs_dir, "AlphaMissense_hg38.subset_targets.tsv.gz")
    chen_subset = os.path.join(refs_dir, "chen_calibration.target_genes.tsv.gz")

    # Reference subsets (AM, Chen) are always chr-prefixed (gencode convention
    # from step 00). The QC'd VCF may use either convention depending on cohort;
    # step 02 detects it per-VCF so `bcftools query -r` matches the file.
    ref_region = args.chr if args.chr.startswith("chr") else f"chr{args.chr}"
    vcf_pfx = vcf_chrom_prefix(args.vcf)
    chr_num = args.chr[3:] if args.chr.startswith("chr") else args.chr
    vcf_region = f"{vcf_pfx}{chr_num}"
    LOG.info("region: ref=%s vcf=%s (vcf prefix=%s)", ref_region, vcf_region, vcf_pfx or "<none>")

    am = load_am_for_region(am_subset, ref_region)
    if not am:
        LOG.warning("no AM rows for %s; output will be empty header-only", ref_region)
    chen = load_chen_for_region(chen_subset, ref_region)

    header = ["sample_id", "cohort", "chr", "pos", "ref", "alt", "gene", "transcript",
              "protein_variant", "am_pathogenicity", "am_class",
              "chen_evidence", "chen_points", "chen_calibration_approach",
              "chen_domain", "chen_vep_score", "in_chen_table",
              "genotype", "DP", "GQ"]
    rows_out: List[List[str]] = []
    n_var_in = 0
    n_var_matched = 0
    n_carrier_rows = 0
    n_in_chen = 0
    for chrom, pos, ref, alt, samples in stream_vcf_genotypes(args.vcf, region=vcf_region):
        n_var_in += 1
        rec = am.get((chrom, pos, ref, alt))
        if rec is None:
            continue
        n_var_matched += 1
        crec = chen.get((chrom, pos, ref, alt))
        in_chen = "TRUE" if crec is not None else "FALSE"
        if crec is not None:
            n_in_chen += 1
        chen_evidence  = crec["evidence"]             if crec else ""
        chen_points    = crec["points"]               if crec else ""
        chen_approach  = crec["calibration_approach"] if crec else ""
        chen_domain    = crec["domain"]               if crec else ""
        chen_vep_score = crec["vep_score"]            if crec else ""
        for sample, gt, dp, gq in samples:
            rows_out.append([
                sample, args.cohort, chrom, str(pos), ref, alt,
                rec["gene_name"], rec["transcript_id"], rec["protein_variant"],
                rec["am_pathogenicity"], rec["am_class"],
                chen_evidence, chen_points, chen_approach, chen_domain, chen_vep_score, in_chen,
                gt, dp, gq,
            ])
            n_carrier_rows += 1

    write_tsv(args.out, header, rows_out)
    LOG.info("%s %s: variants_streamed=%d AM_matched_variants=%d in_Chen=%d carrier_rows=%d",
             args.cohort, ref_region, n_var_in, n_var_matched, n_in_chen, n_carrier_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
