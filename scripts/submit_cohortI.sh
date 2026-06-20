#!/bin/bash
# =============================================================================
# submit_cohortI.sh — chain 01 -> 02 -> 02_gather -> 03 -> 04 for Cohort I.
#
# Usage: bash scripts/submit_cohortI.sh
# Requires: bsub on PATH (Minerva). On the login node, also confirm that
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

J01="biome_am_01_${COHORT}"
J02="biome_am_02_${COHORT}"
J02G="biome_am_02gather_${COHORT}"
J03="biome_am_03_${COHORT}"
J04="biome_am_04_${COHORT}"

log "Cohort I: submitting 01 array [1-22]"
bsub -J "${J01}[1-22]" -env "all, COHORT=$COHORT, CONFIG_PATH=$CONFIG_PATH" < "$REPO_ROOT/scripts/01_qc_missense.lsf"

log "Cohort I: submitting 02 array [1-22] with -w done($J01)"
bsub -J "${J02}[1-22]" -w "done(${J01})" -env "all, COHORT=$COHORT, CONFIG_PATH=$CONFIG_PATH" < "$REPO_ROOT/scripts/02_annotate_am.lsf"

log "Cohort I: submitting 02_gather with -w done($J02)"
bsub -J "$J02G" -w "done(${J02})" -env "all, COHORT=$COHORT, CONFIG_PATH=$CONFIG_PATH" < "$REPO_ROOT/scripts/02_gather.lsf"

log "Cohort I: submitting 03 with -w done($J02G)"
bsub -J "$J03" -w "done(${J02G})" -env "all, COHORT=$COHORT, CONFIG_PATH=$CONFIG_PATH" < "$REPO_ROOT/scripts/03_call_carriers.lsf"

log "Cohort I: submitting 04 with -w done($J03)"
bsub -J "$J04" -w "done(${J03})" -env "all, COHORT=$COHORT, CONFIG_PATH=$CONFIG_PATH" < "$REPO_ROOT/scripts/04_compare_and_tabulate.lsf"

log "Cohort I chain submitted. Monitor with: bjobs -J 'biome_am_*_${COHORT}*'"
