#!/bin/bash
# =============================================================================
# lib/load_modules.sh — Minerva Lmod block. Source this file from every LSF
# script before calling bcftools/tabix/python.
# =============================================================================
# VERIFY: module names below are the typical Minerva names; if your site uses
# different names (e.g. bcftools/1.19 vs bcftools/1.17), edit them once here.
# After loading, this file prints the bcftools/tabix/python versions so any
# version-related issue (filter-expression syntax, tabix index format) is
# visible at the top of every job's stdout.
# =============================================================================

set +u  # Lmod scripts don't like nounset
if command -v module >/dev/null 2>&1; then
    module purge                  || true
    module load bcftools          || echo "VERIFY: module 'bcftools' not found"
    module load htslib            || echo "VERIFY: module 'htslib' not found"   # provides bgzip + tabix
    module load python            || echo "VERIFY: module 'python' not found"   # for pandas / pyyaml
else
    echo "VERIFY: 'module' command not found; assuming bcftools/tabix/python are on PATH (conda env?)" >&2
fi
set -u

# ---- version banner ---------------------------------------------------------
printf '\n[load_modules] tool versions (VERIFY these match the syntax used in build_qc_expr):\n'
printf '  bcftools : '; command -v bcftools >/dev/null && bcftools --version 2>&1 | head -n 1 || echo MISSING
printf '  tabix    : '; command -v tabix    >/dev/null && tabix    --version 2>&1 | head -n 1 || echo MISSING
printf '  bgzip    : '; command -v bgzip    >/dev/null && bgzip    --version 2>&1 | head -n 1 || echo MISSING
printf '  python3  : '; command -v python3  >/dev/null && python3  --version 2>&1                 || echo MISSING
printf '\n'

# ---- sanity-poke the filter expression syntax -------------------------------
# Construct a one-record stream and ensure the filter expression parses. This
# catches "old bcftools doesn't understand FMT/AD[0:1]" at job start rather
# than 6 hours into a chromosome.
if command -v bcftools >/dev/null 2>&1; then
    if ! echo -e '##fileformat=VCFv4.2\n##contig=<ID=chr1>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">\n##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\nchr1\t100\t.\tA\tG\t.\tPASS\t.\tGT:DP:GQ:AD\t0/1:30:40:15,15' \
        | bcftools view -e 'FMT/DP<10 | (GT="het" & (FMT/AD[0:1]/(FMT/AD[0:0]+FMT/AD[0:1])<0.2))' 2>/dev/null >/dev/null
    then
        echo "VERIFY: bcftools rejected the QC filter expression (likely an older version that indexes FMT/AD differently). Inspect build_qc_expr in lib/common.sh." >&2
    fi
fi
