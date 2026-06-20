#!/bin/bash
# =============================================================================
# submit_cohortII.sh — chain index -> 01 -> 02 -> 02_gather -> 03 -> 04 for Cohort II.
#
# The discrete index_cohortII_source job runs FIRST so all 22 array tasks see
# the .csi (no race). The source VCF stays in /sc/private/ untouched; the .csi
# is written under intermediate/cohortII/.
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

JIDX="biome_am_index_${COHORT}"
J01="biome_am_01_${COHORT}"
J02="biome_am_02_${COHORT}"
J02G="biome_am_02gather_${COHORT}"
J03="biome_am_03_${COHORT}"
J04="biome_am_04_${COHORT}"

log "Cohort II: submitting source-index pre-step"
bsub -J "$JIDX" -env "all, CONFIG_PATH=$CONFIG_PATH" < "$REPO_ROOT/scripts/index_cohortII_source.lsf"

log "Cohort II: submitting 01 array [1-22] with -w done($JIDX)"
bsub -J "${J01}[1-22]" -w "done(${JIDX})" -env "all, COHORT=$COHORT, CONFIG_PATH=$CONFIG_PATH" < "$REPO_ROOT/scripts/01_qc_missense.lsf"

log "Cohort II: submitting 02 array [1-22] with -w done($J01)"
bsub -J "${J02}[1-22]" -w "done(${J01})" -env "all, COHORT=$COHORT, CONFIG_PATH=$CONFIG_PATH" < "$REPO_ROOT/scripts/02_annotate_am.lsf"

log "Cohort II: submitting 02_gather with -w done($J02)"
bsub -J "$J02G" -w "done(${J02})" -env "all, COHORT=$COHORT, CONFIG_PATH=$CONFIG_PATH" < "$REPO_ROOT/scripts/02_gather.lsf"

log "Cohort II: submitting 03 with -w done($J02G)"
bsub -J "$J03" -w "done(${J02G})" -env "all, COHORT=$COHORT, CONFIG_PATH=$CONFIG_PATH" < "$REPO_ROOT/scripts/03_call_carriers.lsf"

log "Cohort II: submitting 04 with -w done($J03)"
bsub -J "$J04" -w "done(${J03})" -env "all, COHORT=$COHORT, CONFIG_PATH=$CONFIG_PATH" < "$REPO_ROOT/scripts/04_compare_and_tabulate.lsf"

log "Cohort II chain submitted. Monitor with: bjobs -J 'biome_am_*_${COHORT}*'"
