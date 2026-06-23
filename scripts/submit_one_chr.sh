#!/bin/bash
# =============================================================================
# submit_one_chr.sh — bsub a single chromosome of one step. Useful for smoke
# testing or rerunning a failed array task without re-submitting the whole
# pipeline.
#
# Usage:
#   bash scripts/submit_one_chr.sh <cohort> <step> <chr_index>
#
# Examples:
#   bash scripts/submit_one_chr.sh cohortI 01_qc_missense 17
#   bash scripts/submit_one_chr.sh cohortI 02_annotate_am 17
#   bash scripts/submit_one_chr.sh cohortII 01_qc_missense 7
#
# This generates the same kind of wrapper that submit_cohortI.sh does, but
# for a single array index (no `[1-22]` and no `-w done(...)`).
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/.." && pwd)"
source "$REPO_ROOT/lib/common.sh"

COHORT="${1:-}"; STEP="${2:-}"; CHR_IDX="${3:-}"
[[ -n "$COHORT" && -n "$STEP" && -n "$CHR_IDX" ]] || \
    die "usage: bash scripts/submit_one_chr.sh <cohort: cohortI|cohortII> <step: 01_qc_missense|02_annotate_am> <chr_index 1..22>"
[[ "$COHORT" == "cohortI" || "$COHORT" == "cohortII" ]] || die "cohort must be cohortI or cohortII (got $COHORT)"
[[ "$CHR_IDX" =~ ^[0-9]+$ ]] || die "chr_index must be a positive integer (got $CHR_IDX)"

CONFIG_PATH="${CONFIG_PATH:-$REPO_ROOT/config/config.yaml}"
export CONFIG_PATH

# Map STEP -> lsf-script-name and config key
case "$STEP" in
    01_qc_missense)  LSF_NAME="01_qc_missense.lsf";  STEP_KEY="step_01_qc_missense";;
    02_annotate_am)  LSF_NAME="02_annotate_am.lsf";  STEP_KEY="step_02_annotate_am";;
    *) die "unsupported step: $STEP (only 01_qc_missense and 02_annotate_am are array steps)";;
esac

PROJECT="$(cfg_get lsf.project)"
QUEUE="$(cfg_get lsf.queue)"
N="$(cfg_get "lsf.${STEP_KEY}.n")"
W="$(cfg_get "lsf.${STEP_KEY}.W")"
MEM="$(cfg_get "lsf.${STEP_KEY}.mem")"

JNAME="biome_am_${STEP}_${COHORT}_one_chr${CHR_IDX}[${CHR_IDX}]"

log "submitting one-chr job: cohort=$COHORT step=$STEP chr=chr$CHR_IDX (job name: $JNAME)"

cat <<WRAPPER | bsub
#!/bin/bash
#BSUB -J ${JNAME}
#BSUB -P ${PROJECT}
#BSUB -q ${QUEUE}
#BSUB -n ${N}
#BSUB -W ${W}
#BSUB -R rusage[mem=${MEM}]
#BSUB -R span[hosts=1]
#BSUB -o ${REPO_ROOT}/logs/${STEP}.%J.%I.stdout
#BSUB -eo ${REPO_ROOT}/logs/${STEP}.%J.%I.stderr
#BSUB -L /bin/bash

# >>> values baked in by submit_one_chr.sh at submission time <<<
export BIOAM_REPO_ROOT=${REPO_ROOT}
export COHORT=${COHORT}
export CONFIG_PATH=${CONFIG_PATH}

cd "\$BIOAM_REPO_ROOT"
exec bash "\$BIOAM_REPO_ROOT/scripts/${LSF_NAME}"
WRAPPER

log "done. Monitor with: bjobs -J 'biome_am_${STEP}_${COHORT}_one_chr${CHR_IDX}*'"
