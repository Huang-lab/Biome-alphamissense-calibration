#!/usr/bin/env bash
# Pre-commit hook: refuse to commit patient-identifiable BioMe data.
#
# Install:
#     ln -sf ../../scripts/pre-commit.sh .git/hooks/pre-commit
#     chmod +x .git/hooks/pre-commit
#
# Two checks run against staged content (not the working tree):
#   1. Path patterns matching real-data filename conventions are rejected.
#      Synthetic fixtures under tests/data/ are intentionally allowed.
#   2. Staged file *content* is scanned for SINAI_ID and MRN-map header
#      patterns. This catches leaks that ride in via a CSV/log/notebook
#      that doesn't match the filename rules.
#
# Bypass with --no-verify ONLY with explicit user OK. The hook prints
# which file failed which check so reviewers can see the exact match.

set -u

forbidden_paths=(
    'Regen.*VariantsInSamples.*\.tsv'
    'Sema4.*VariantsInSamples.*\.tsv'
    'WXS_HX_Newgroups.*\.tsv'
    'Masked_mrn_map.*\.txt'
    'GSA_GDA_PCA.*\.txt'
    'A_variant_master\.tsv$'
    'B_carriers_per_sample\.tsv$'
    'B_summary_by_gene\.tsv$'
    'C_regression_matrix\.tsv$'
    '\.vcf\.gz$'
)

# Content patterns. Synthetic fixture IDs (S1..S99) are NOT matched.
content_patterns=(
    'SINAI_[0-9]+_[A-Za-z0-9]+'
    '^PLATFORM[[:space:]]+RGNID[[:space:]]+MASKED_MRN'
)

staged=$(git diff --cached --name-only --diff-filter=ACM)
if [[ -z "${staged}" ]]; then
    exit 0
fi

fail=0

while IFS= read -r path; do
    # Allowlist: synthetic fixtures under tests/data/ may contain S1..S99 sample
    # IDs and synthetic VCFs.
    if [[ "${path}" == tests/data/* ]]; then
        continue
    fi
    for pat in "${forbidden_paths[@]}"; do
        if [[ "${path}" =~ ${pat} ]]; then
            printf 'REJECT path: %s (matches /%s/)\n' "${path}" "${pat}" >&2
            fail=1
        fi
    done
    if [[ ! -f "${path}" ]]; then
        continue
    fi
    # Skip binary blobs for the content scan.
    if ! git diff --cached --numstat -- "${path}" | awk '$1 != "-" {exit 0} {exit 1}'; then
        continue
    fi
    for pat in "${content_patterns[@]}"; do
        if git show ":${path}" 2>/dev/null | grep -nE "${pat}" >/dev/null; then
            hit=$(git show ":${path}" | grep -nE "${pat}" | head -3)
            printf 'REJECT content in %s (pattern /%s/):\n%s\n' "${path}" "${pat}" "${hit}" >&2
            fail=1
        fi
    done
done <<< "${staged}"

if [[ ${fail} -ne 0 ]]; then
    cat >&2 <<'MSG'

Commit blocked: staged content matches a sensitive-data rule.
If a match is a false positive, narrow the pattern in scripts/pre-commit.sh
(do not bypass with --no-verify without explicit project lead approval).
MSG
    exit 1
fi

exit 0
