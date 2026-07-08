#!/bin/bash
# =============================================================================
# submit_clinvar.sh — ClinVar-only rerun: chain 02c -> 02c_gather -> 04 for one
# cohort. Use this after adding the ClinVar category when steps 01/02/03 (QC'd
# VCF + Table A) are ALREADY DONE from a prior run and do not need rerunning.
#
# Mirrors submit_cohort*.sh: it bakes COHORT / CONFIG_PATH / BIOAM_REPO_ROOT
# into an on-the-fly wrapper so they survive LSF's stdin copy (raw `bsub <` does
# NOT propagate these on the Minerva LSF install).
#
# Usage: bash scripts/submit_clinvar.sh cohortI     (then cohortII)
# Requires: bsub on PATH; scripts/00_prepare_refs.sh already ran (ClinVar subset
# present); intermediate/<cohort>/chr<N>.qc.vcf.gz present from step 01.
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
[[ -s "$INTERMED/chr1.qc.vcf.gz" ]] \
    || die "$INTERMED/chr1.qc.vcf.gz missing; step 01 (QC) must have run for $COHORT"
mkdir -p "$REPO_ROOT/logs" "$INTERMED"

PROJECT="$(cfg_get lsf.project)"
QUEUE="$(cfg_get lsf.queue)"

# emit_wrapper STEP_KEY ARRAY_SPEC WAIT_CLAUSE LSF_SCRIPT_NAME  (see submit_cohortI.sh)
emit_wrapper() {
    local step_key="$1" array_spec="$2" wait_clause="$3" lsf_name="$4"
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

cd "\$BIOAM_REPO_ROOT"
exec bash "\$BIOAM_REPO_ROOT/scripts/${lsf_name}"
WRAPPER
}

J02C="biome_am_02c_annotate_clinvar_${COHORT}"
J02CG="biome_am_02c_gather_${COHORT}"

log "$COHORT: submitting 02c ClinVar annotate array [1-22]"
emit_wrapper step_02c_annotate_clinvar "[1-22]" ""                02c_annotate_clinvar.lsf    | bsub

log "$COHORT: submitting 02c_gather gated on 02c"
emit_wrapper step_02c_gather           ""       "done(${J02C})"   02c_gather_clinvar.lsf      | bsub

log "$COHORT: submitting 04 (re-tabulate with ClinVar) gated on 02c_gather"
emit_wrapper step_04_compare_tab       ""       "done(${J02CG})"  04_compare_and_tabulate.lsf | bsub

log "$COHORT ClinVar chain submitted. Monitor with: bjobs -J 'biome_am_*_${COHORT}'"
