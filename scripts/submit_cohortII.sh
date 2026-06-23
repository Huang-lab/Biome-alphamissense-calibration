#!/bin/bash
# =============================================================================
# submit_cohortII.sh — chain index -> 01 -> 02 -> 02_gather -> 03 -> 04 for Cohort II.
#
# The discrete index_cohortII_source job runs FIRST so all 22 array tasks see
# the .csi (no race). The source VCF stays in /sc/private/ untouched; the .csi
# is written under intermediate/cohortII/.
#
# Same wrapper-generation approach as submit_cohortI.sh — see that file for
# the rationale (Minerva LSF drops -env user vars; we bake them into the
# wrapper script body as plain `export` statements).
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/.." && pwd)"
source "$REPO_ROOT/lib/common.sh"

CONFIG_PATH="${CONFIG_PATH:-$REPO_ROOT/config/config.yaml}"
export CONFIG_PATH
COHORT=cohortII

command -v bsub >/dev/null 2>&1 || die "bsub not on PATH; submit from a Minerva login node"
[[ -s "$REPO_ROOT/refs/REFERENCE_REPORT.md" ]] || die "refs/REFERENCE_REPORT.md missing; run scripts/00_prepare_refs.sh first"
mkdir -p "$REPO_ROOT/logs" "$REPO_ROOT/intermediate/$COHORT"

PROJECT="$(cfg_get lsf.project)"
QUEUE="$(cfg_get lsf.queue)"

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

# >>> values baked in by submit_cohortII.sh at submission time <<<
export BIOAM_REPO_ROOT=${REPO_ROOT}
export COHORT=${COHORT}
export CONFIG_PATH=${CONFIG_PATH}

cd "\$BIOAM_REPO_ROOT"
exec bash "\$BIOAM_REPO_ROOT/scripts/${lsf_name}"
WRAPPER
}

JIDX="biome_am_index_cohortII_${COHORT}"
J01="biome_am_01_qc_missense_${COHORT}"
J02="biome_am_02_annotate_am_${COHORT}"
J02G="biome_am_02_gather_${COHORT}"
J03="biome_am_03_call_carriers_${COHORT}"

log "Cohort II: submitting source-index pre-step"
emit_wrapper step_index_cohortII   ""       ""               index_cohortII_source.lsf | bsub

log "Cohort II: submitting 01 array [1-22] gated on index"
emit_wrapper step_01_qc_missense   "[1-22]" "done(${JIDX})"  01_qc_missense.lsf       | bsub

log "Cohort II: submitting 02 array [1-22] gated on 01"
emit_wrapper step_02_annotate_am   "[1-22]" "done(${J01})"   02_annotate_am.lsf       | bsub

log "Cohort II: submitting 02_gather gated on 02"
emit_wrapper step_02_gather        ""       "done(${J02})"   02_gather.lsf            | bsub

log "Cohort II: submitting 03 gated on 02_gather"
emit_wrapper step_03_call_carriers ""       "done(${J02G})"  03_call_carriers.lsf     | bsub

log "Cohort II: submitting 04 gated on 03"
emit_wrapper step_04_compare_tab   ""       "done(${J03})"   04_compare_and_tabulate.lsf | bsub

log "Cohort II chain submitted. Monitor with: bjobs -J 'biome_am_*_${COHORT}'"
