#!/bin/bash
# =============================================================================
# submit_clinvar.sh — ClinVar-only rerun for one cohort:
#   01c (all-variant QC, no SNV filter) -> 02c annotate -> 02c_gather -> 04.
# Use this after adding the ClinVar category. It does NOT touch the AM/Chen path
# (01 SNV QC, 02, 03); Table A is reused as-is. Only step 01 is re-run in an
# ALL-VARIANT mode (separate chr<N>.qc_allvar.vcf.gz output) so ClinVar can see
# non-SNV P/LP variants.
#
# Mirrors submit_cohort*.sh: it bakes COHORT / CONFIG_PATH / BIOAM_REPO_ROOT
# (and QC_VARIANT_MODE for 01c) into an on-the-fly wrapper so they survive LSF's
# stdin copy (raw `bsub <` does NOT propagate these on the Minerva LSF install).
#
# Usage: bash scripts/submit_clinvar.sh cohortI     (then cohortII)
# Requires: bsub on PATH; scripts/00_prepare_refs.sh already ran (ClinVar subset
# present); the cohort's raw input VCFs available (01c re-reads them).
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/.." && pwd)"
source "$REPO_ROOT/lib/common.sh"

CONFIG_PATH="${CONFIG_PATH:-$REPO_ROOT/config/config.yaml}"
export CONFIG_PATH
COHORT="${1:-}"
[[ "$COHORT" == "cohortI" || "$COHORT" == "cohortII" ]] \
    || die "usage: bash scripts/submit_clinvar.sh cohortI|cohortII"

command -v bsub >/dev/null 2>&1 || die "bsub not on PATH; submit from a Minerva login node"
[[ -s "$REPO_ROOT/refs/clinvar_plp_2star.target_genes.tsv.gz" ]] \
    || die "refs/clinvar_plp_2star.target_genes.tsv.gz missing; run scripts/00_prepare_refs.sh first"
INTERMED="$REPO_ROOT/intermediate/$COHORT"
mkdir -p "$REPO_ROOT/logs" "$INTERMED"

REF_FASTA="$(cfg_get references.reference_fasta 2>/dev/null || echo '')"
[[ -n "$REF_FASTA" ]] \
    || warn "references.reference_fasta unset — indel matching will be best-effort. Set it in config.local.yaml for reliable indel normalization."
# Cohort II (sliced_from_combined) needs its indexed local symlink first.
if [[ "$COHORT" == "cohortII" ]]; then
    LOCAL_DIR="$(cfg_get cohorts.cohortII.local_source_dir)"
    case "$LOCAL_DIR" in /*) ;; *) LOCAL_DIR="$REPO_ROOT/$LOCAL_DIR";; esac
    SRC="$(cfg_get cohorts.cohortII.combined_vcf)"
    [[ -s "${LOCAL_DIR}/$(basename "$SRC").csi" ]] \
        || warn "cohortII local index missing under $LOCAL_DIR; submit scripts/index_cohortII_source.lsf first if 01c fails"
fi

PROJECT="$(cfg_get lsf.project)"
QUEUE="$(cfg_get lsf.queue)"

# emit_wrapper STEP_KEY ARRAY_SPEC WAIT_CLAUSE LSF_SCRIPT_NAME [EXTRA_ENV]
#   EXTRA_ENV: optional extra "export FOO=bar" line baked into the wrapper.
emit_wrapper() {
    local step_key="$1" array_spec="$2" wait_clause="$3" lsf_name="$4" extra_env="${5:-}"
    local job_name="biome_am_${step_key#step_}_${COHORT}"
    local n W mem
    n="$(cfg_get "lsf.${step_key}.n")"
    W="$(cfg_get "lsf.${step_key}.W")"
    mem="$(cfg_get "lsf.${step_key}.mem")"
    local jname="${job_name}${array_spec}"
    cat <<WRAPPER
#!/bin/bash
#BSUB -J ${jname}
#BSUB -P ${PROJECT}
#BSUB -q ${QUEUE}
#BSUB -n ${n}
#BSUB -W ${W}
#BSUB -R rusage[mem=${mem}]
#BSUB -R span[hosts=1]
#BSUB -o ${REPO_ROOT}/logs/${step_key#step_}.%J.%I.stdout
#BSUB -eo ${REPO_ROOT}/logs/${step_key#step_}.%J.%I.stderr
#BSUB -L /bin/bash
$( [[ -n "$wait_clause" ]] && echo "#BSUB -w \"${wait_clause}\"" )

# >>> values baked in by submit_clinvar.sh at submission time <<<
export BIOAM_REPO_ROOT=${REPO_ROOT}
export COHORT=${COHORT}
export CONFIG_PATH=${CONFIG_PATH}
${extra_env}

cd "\$BIOAM_REPO_ROOT"
exec bash "\$BIOAM_REPO_ROOT/scripts/${lsf_name}"
WRAPPER
}

J01C="biome_am_01c_qc_clinvar_${COHORT}"
J02C="biome_am_02c_annotate_clinvar_${COHORT}"
J02CG="biome_am_02c_gather_${COHORT}"

log "$COHORT: submitting 01c all-variant QC array [1-22] (QC_VARIANT_MODE=all)"
emit_wrapper step_01c_qc_clinvar   "[1-22]" ""               01_qc_missense.lsf       "export QC_VARIANT_MODE=all" | bsub

log "$COHORT: submitting 02c ClinVar annotate array [1-22] gated on 01c"
emit_wrapper step_02c_annotate_clinvar "[1-22]" "done(${J01C})"   02c_annotate_clinvar.lsf    | bsub

log "$COHORT: submitting 02c_gather gated on 02c"
emit_wrapper step_02c_gather           ""       "done(${J02C})"   02c_gather_clinvar.lsf      | bsub

log "$COHORT: submitting 04 (re-tabulate with ClinVar) gated on 02c_gather"
emit_wrapper step_04_compare_tab       ""       "done(${J02CG})"  04_compare_and_tabulate.lsf | bsub

log "$COHORT ClinVar chain submitted. Monitor with: bjobs -J 'biome_am_*_${COHORT}'"
