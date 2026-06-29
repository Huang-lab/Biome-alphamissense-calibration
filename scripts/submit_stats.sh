#!/bin/bash
# =============================================================================
# submit_stats.sh — submit 05_run_stats for one or both cohorts via LSF.
#
# Like submit_cohortI.sh / submit_cohortII.sh, this generates a wrapper script
# per submission with COHORT / CONFIG_PATH / BIOAM_REPO_ROOT baked in as plain
# `export` statements — bsub's stdin-copy carries the wrapper to the execute
# host, so the env vars survive (the Minerva LSF install does NOT propagate
# user vars passed via `bsub -env "..."`).
#
# Usage:
#   bash scripts/submit_stats.sh                # both cohorts
#   bash scripts/submit_stats.sh cohortI        # one cohort
#   bash scripts/submit_stats.sh cohortII
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/.." && pwd)"
source "$REPO_ROOT/lib/common.sh"

CONFIG_PATH="${CONFIG_PATH:-$REPO_ROOT/config/config.yaml}"
export CONFIG_PATH

command -v bsub >/dev/null 2>&1 || die "bsub not on PATH; submit from a Minerva login node"
mkdir -p "$REPO_ROOT/logs"

PROJECT="$(cfg_get lsf.project)"
QUEUE="$(cfg_get lsf.queue)"
N="$(cfg_get lsf.step_05_run_stats.n)"
W="$(cfg_get lsf.step_05_run_stats.W)"
MEM="$(cfg_get lsf.step_05_run_stats.mem)"

emit_wrapper_and_submit() {
    local cohort="$1"
    local jname="biome_am_05_run_stats_${cohort}"
    log "submitting ${jname}"
    cat <<WRAPPER | bsub
#!/bin/bash
#BSUB -J ${jname}
#BSUB -P ${PROJECT}
#BSUB -q ${QUEUE}
#BSUB -n ${N}
#BSUB -W ${W}
#BSUB -R rusage[mem=${MEM}]
#BSUB -R span[hosts=1]
#BSUB -o ${REPO_ROOT}/logs/05_run_stats.%J.${cohort}.stdout
#BSUB -eo ${REPO_ROOT}/logs/05_run_stats.%J.${cohort}.stderr
#BSUB -L /bin/bash

# >>> values baked in by submit_stats.sh at submission time <<<
export BIOAM_REPO_ROOT=${REPO_ROOT}
export COHORT=${cohort}
export CONFIG_PATH=${CONFIG_PATH}

cd "\$BIOAM_REPO_ROOT"
exec bash "\$BIOAM_REPO_ROOT/scripts/05_run_stats.lsf"
WRAPPER
}

if [[ $# -eq 0 ]]; then
    emit_wrapper_and_submit cohortI
    emit_wrapper_and_submit cohortII
else
    for c in "$@"; do
        case "$c" in
            cohortI|cohortII) emit_wrapper_and_submit "$c" ;;
            *) die "unknown cohort: $c (expected cohortI or cohortII)" ;;
        esac
    done
fi

log "Monitor with: bjobs -J 'biome_am_05_run_stats_*'"
