#!/bin/bash
# =============================================================================
# submit_cohortI.sh — chain 01 -> 02 -> 02_gather -> 03 -> 04 for Cohort I.
#
# Instead of relying on `bsub -env "all, COHORT=..., CONFIG_PATH=..."` (which
# does not propagate user vars under the Minerva LSF install we hit), we
# generate a tiny WRAPPER script per step on the fly. The wrapper:
#   1. carries the #BSUB directives (resources from config.yaml)
#   2. exports COHORT / CONFIG_PATH / BIOAM_REPO_ROOT as plain bash statements
#      so they survive LSF's stdin-copy to /local/JOBS/<id>/shell
#   3. exec's bash on the real .lsf work script
#
# Usage: bash scripts/submit_cohortI.sh
# Requires: bsub on PATH. On the login node, also confirm that
# 00_prepare_refs.sh ran and refs/REFERENCE_REPORT.md was reviewed.
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/.." && pwd)"
source "$REPO_ROOT/lib/common.sh"

CONFIG_PATH="${CONFIG_PATH:-$REPO_ROOT/config/config.yaml}"
export CONFIG_PATH
COHORT=cohortI

command -v bsub >/dev/null 2>&1 || die "bsub not on PATH; submit from a Minerva login node"
[[ -s "$REPO_ROOT/refs/REFERENCE_REPORT.md" ]] || die "refs/REFERENCE_REPORT.md missing; run scripts/00_prepare_refs.sh first"
mkdir -p "$REPO_ROOT/logs" "$REPO_ROOT/intermediate/$COHORT"

# LSF resources (read once from config)
PROJECT="$(cfg_get lsf.project)"
QUEUE="$(cfg_get lsf.queue)"

# emit_wrapper STEP_KEY ARRAY_SPEC WAIT_CLAUSE LSF_SCRIPT_NAME
#   STEP_KEY:        cfg path under `lsf.`, e.g. step_01_qc_missense
#   ARRAY_SPEC:      "[1-22]" or ""  (empty for non-array jobs)
#   WAIT_CLAUSE:     bsub -w expression or ""  (e.g. "done(biome_am_01_cohortI)")
#   LSF_SCRIPT_NAME: name under scripts/, e.g. 01_qc_missense.lsf
#
# Prints the full wrapper to stdout. Caller pipes it to `bsub`.
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

# >>> values baked in by submit_cohortI.sh at submission time <<<
export BIOAM_REPO_ROOT=${REPO_ROOT}
export COHORT=${COHORT}
export CONFIG_PATH=${CONFIG_PATH}

cd "\$BIOAM_REPO_ROOT"
exec bash "\$BIOAM_REPO_ROOT/scripts/${lsf_name}"
WRAPPER
}

J01="biome_am_01_qc_missense_${COHORT}"
J02="biome_am_02_annotate_am_${COHORT}"
J02G="biome_am_02_gather_${COHORT}"
J03="biome_am_03_call_carriers_${COHORT}"
J04="biome_am_04_compare_tab_${COHORT}"

log "Cohort I: submitting 01 array [1-22]"
emit_wrapper step_01_qc_missense   "[1-22]" ""               01_qc_missense.lsf       | bsub

log "Cohort I: submitting 02 array [1-22] gated on 01"
emit_wrapper step_02_annotate_am   "[1-22]" "done(${J01})"   02_annotate_am.lsf       | bsub

log "Cohort I: submitting 02_gather gated on 02"
emit_wrapper step_02_gather        ""       "done(${J02})"   02_gather.lsf            | bsub

log "Cohort I: submitting 03 gated on 02_gather"
emit_wrapper step_03_call_carriers ""       "done(${J02G})"  03_call_carriers.lsf     | bsub

log "Cohort I: submitting 04 gated on 03"
emit_wrapper step_04_compare_tab   ""       "done(${J03})"   04_compare_and_tabulate.lsf | bsub

log "Cohort I chain submitted. Monitor with: bjobs -J 'biome_am_*_${COHORT}'"
