#!/bin/bash
# =============================================================================
# lib/common.sh — shared bash helpers (logging, guards, gather/verify-22,
#                 QC expression builder). source this file from every script.
# =============================================================================
# All callers must use `set -euo pipefail`. This file registers an ERR trap
# that prints the failing command + stderr tail + a one-line cause hypothesis.
# =============================================================================

# ---- logging ----------------------------------------------------------------
_ts() { date "+%Y-%m-%dT%H:%M:%S%z"; }
log()  { printf '[%s] %s\n' "$(_ts)" "$*" >&2; }
warn() { printf '[%s] WARN: %s\n' "$(_ts)" "$*" >&2; }
die()  { printf '[%s] ERROR: %s\n' "$(_ts)" "$*" >&2; exit 1; }

# ---- error trap with cause hypothesis ---------------------------------------
# Heuristics inspect $? and the failing command to suggest the most likely cause.
__on_err() {
    local rc=$?
    local cmd="${BASH_COMMAND:-?}"
    local line="${BASH_LINENO[0]:-?}"
    local src="${BASH_SOURCE[1]:-?}"
    {
        printf '\n========== PIPELINE ERROR ==========\n'
        printf 'exit_code : %d\n' "$rc"
        printf 'at        : %s:%s\n' "$src" "$line"
        printf 'command   : %s\n' "$cmd"
        printf 'hypothesis: '
        case "$cmd" in
            *gvcf*|*NON_REF*|*"GVCFBlock"*)  echo "input appears to be a combined gVCF, not GenotypeGVCFs output" ;;
            *multiallelic*|*"awk"*ALT*)      echo "multiallelic record encountered; pipeline requires biallelic" ;;
            *tabix*)                         echo "tabix lookup failed; index may be missing or CHROM naming mismatch" ;;
            *bcftools*view*)                 echo "bcftools view failed; check region naming (chr<N>), index, and VCF integrity" ;;
            *bcftools*index*)                echo "bcftools index failed; output directory may be unwritable or input truncated" ;;
            *"transcript not found"*)        echo "transcript not in AlphaMissense file; check version-stripped fallback" ;;
            *"coordinates unresolved"*)      echo "could not resolve exon coordinates from gencode GTF for at least one gene" ;;
            *"ID column"*)                   echo "phenotype/PC ID column does not match VCF samples — see unmatched IDs above" ;;
            *"missing chr"*|*gather*)        echo "one or more chromosomes missing from gather inputs; an array task likely failed silently" ;;
            *requests*|*curl*|*wget*)        echo "network fetch failed; drop a local file in via config and re-run" ;;
            python*|*"py "*|*.py)            echo "python step failed; scroll up for the python traceback / explicit message" ;;
            *)                               echo "see stderr above and the failing command for context" ;;
        esac
        printf '====================================\n\n'
    } >&2
    exit "$rc"
}
trap __on_err ERR

# ---- repo-root discovery ----------------------------------------------------
# Resolve the repository root from the location of this file. Works whether the
# caller is in scripts/, tests/, or repo root.
repo_root() {
    local here
    here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    (cd "$here/.." && pwd)
}

# ---- yaml config reader -----------------------------------------------------
# Tiny yaml reader for the values this pipeline needs. We use python+pyyaml
# rather than yq so we have no extra system dependency.
# Usage: cfg_get <key.path>
cfg_get() {
    local key="$1"
    local cfg="${CONFIG_PATH:-$(repo_root)/config/config.yaml}"
    python3 - "$cfg" "$key" <<'PY'
import sys, yaml, json
cfg_path, key = sys.argv[1], sys.argv[2]
with open(cfg_path) as fh:
    doc = yaml.safe_load(fh)
node = doc
for part in key.split('.'):
    if isinstance(node, dict) and part in node:
        node = node[part]
    else:
        print('', end='')
        sys.exit(0)
if isinstance(node, (dict, list)):
    print(json.dumps(node))
else:
    print('' if node is None else node)
PY
}

# ---- guards -----------------------------------------------------------------

# guard_not_gvcf <vcf.gz>
# Fails with the exact gVCF abort message if the input looks like a gVCF.
guard_not_gvcf() {
    local vcf="$1"
    [[ -f "$vcf" ]] || die "guard_not_gvcf: file not found: $vcf"
    # Check header for ##GVCFBlock contigs.
    if bcftools view -h "$vcf" 2>/dev/null | grep -qE '##GVCFBlock|ID=MIN_DP|ID=RGQ'; then
        # Header alone isn't conclusive (Regeneron files reportedly have RGQ in
        # the header dictionary even after GenotypeGVCFs). Confirm with a body
        # peek: any <NON_REF> ALT or any reference-only block (ALT == '.') is
        # disqualifying.
        if bcftools view -H "$vcf" 2>/dev/null | head -n 200 | awk -F'\t' '
            $5 == "<NON_REF>" || $5 == "." { found=1; exit }
            $5 ~ /<NON_REF>/                 { found=1; exit }
            END { exit (found ? 0 : 1) }
        '; then
            die "Input appears to be a combined gVCF, not genotyped calls; point at the GenotypeGVCFs output. [$vcf]"
        fi
    fi
    # Even without GVCFBlock header markers, body check for <NON_REF> catches it.
    if bcftools view -H "$vcf" 2>/dev/null | head -n 500 | awk -F'\t' '
        $5 == "<NON_REF>" || $5 ~ /,<NON_REF>/ { found=1; exit }
        END { exit (found ? 0 : 1) }
    '; then
        die "Input appears to be a combined gVCF, not genotyped calls; point at the GenotypeGVCFs output. [$vcf]"
    fi
}

# guard_no_multiallelic <vcf.gz>
# Fails if any record has comma in ALT.
guard_no_multiallelic() {
    local vcf="$1"
    [[ -f "$vcf" ]] || die "guard_no_multiallelic: file not found: $vcf"
    local found
    found="$(bcftools view -H "$vcf" 2>/dev/null | awk -F'\t' '$5 ~ /,/ {print; exit}' || true)"
    if [[ -n "$found" ]]; then
        die "Multiallelic record encountered; pipeline requires biallelic. Offender: $(echo "$found" | cut -f1-5)"
    fi
}

# guard_snv_only <tsv>
# Fails if any (REF,ALT) pair in the joined AM table is not a single-base SNV.
# Expects columns 'ref' and 'alt' to be discoverable from header.
guard_snv_only() {
    local tsv="$1"
    [[ -f "$tsv" ]] || die "guard_snv_only: file not found: $tsv"
    python3 - "$tsv" <<'PY'
import sys, csv
path = sys.argv[1]
opener = open
if path.endswith('.gz'):
    import gzip
    opener = lambda p: gzip.open(p, 'rt')
with opener(path) as fh:
    rdr = csv.reader(fh, delimiter='\t')
    header = next(rdr)
    try:
        ri = header.index('ref'); ai = header.index('alt')
    except ValueError:
        # try common alternatives
        lc = [h.lower() for h in header]
        ri = lc.index('ref'); ai = lc.index('alt')
    bad = 0
    for row in rdr:
        if len(row[ri]) != 1 or len(row[ai]) != 1:
            bad += 1
            if bad <= 3:
                print('non-SNV:', row[:5], file=sys.stderr)
    if bad:
        sys.exit(f'AM-matched non-SNV rows: {bad}')
PY
}

# gather_verify_22 <dir> <pattern_with_{chr}>
# Lists missing chr1..chr22 files and exits non-zero if any are missing.
gather_verify_22() {
    local dir="$1" pat="$2"
    local missing=()
    for i in $(seq 1 22); do
        local p="${dir}/${pat//\{chr\}/chr${i}}"
        if [[ ! -s "$p" ]]; then
            missing+=("chr${i}: $p")
        fi
    done
    if (( ${#missing[@]} > 0 )); then
        printf 'missing chr in gather:\n' >&2
        printf '  %s\n' "${missing[@]}" >&2
        die "gather_verify_22: ${#missing[@]} chromosome(s) missing; aborting"
    fi
    log "gather_verify_22: all 22 outputs present under $dir"
}

# build_qc_expr <dp_min> <gq_min> <ab_lo> <ab_hi>
# Returns a bcftools `+setGT -i ...` FAILURE expression. Genotypes matching
# the expression are converted to missing (./.) by setGT -n .. The expression
# is centralized HERE so the exact syntax is in ONE place.
#
# VERIFY: this AB syntax was tested against bcftools 1.17 and 1.19. Older
# bcftools (<1.10) indexed FMT/AD differently (FMT/AD[0] vs FMT/AD[*:0]); the
# loaded module's `bcftools --version` is logged by lib/load_modules.sh so any
# divergence is visible at the top of every job's stdout.
build_qc_expr() {
    local dp="$1" gq="$2" ab_lo="$3" ab_hi="$4"
    # NB: in setGT context, the expression is evaluated per (variant, sample).
    # FMT/AD[0] = REF count of the current sample, FMT/AD[1] = ALT count.
    # We mask if: DP too low, OR GQ too low, OR (het AND AB outside band).
    # We also mask hom-ref and missing here so step 02 sees only "carrier-like"
    # genotypes for the AM-matched variants.
    printf '%s' "FMT/DP<${dp} | FMT/GQ<${gq} | GT=\"RR\" | GT=\"mis\" | (GT=\"het\" & (FMT/AD[0:1]/(FMT/AD[0:0]+FMT/AD[0:1])<${ab_lo} | FMT/AD[0:1]/(FMT/AD[0:0]+FMT/AD[0:1])>${ab_hi}))"
}

# count_records <vcf.gz>
count_records() {
    bcftools view -H "$1" 2>/dev/null | wc -l | awk '{print $1}'
}
