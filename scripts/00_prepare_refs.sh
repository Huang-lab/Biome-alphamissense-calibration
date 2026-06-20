#!/bin/bash
# =============================================================================
# 00_prepare_refs.sh — LOGIN NODE script (internet required), no LSF.
#
# Builds: refs/calibration_thresholds.tsv, refs/gene_transcript_map.tsv,
#         refs/target_genes.exons.bed, refs/AlphaMissense_hg38.subset_targets.tsv.gz (+ .tbi),
#         refs/REFERENCE_REPORT.md
# Hard-gates: prints coverage counts and tells you to review REFERENCE_REPORT.md
#             before submitting 01+.
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/.." && pwd)"
# shellcheck source=../lib/common.sh
source "$REPO_ROOT/lib/common.sh"
# shellcheck source=../lib/load_modules.sh
source "$REPO_ROOT/lib/load_modules.sh"

CONFIG_PATH="${CONFIG_PATH:-$REPO_ROOT/config/config.yaml}"
export CONFIG_PATH

log "00_prepare_refs: starting on $(hostname) at $(_ts)"
log "config: $CONFIG_PATH"

cd "$REPO_ROOT"

# requirements check
for bin in python3 bgzip tabix; do
    command -v "$bin" >/dev/null 2>&1 || die "required binary not on PATH: $bin"
done

python3 "$REPO_ROOT/python/prepare_refs.py" --config "$CONFIG_PATH" "$@"

log "00_prepare_refs: done. Review refs/REFERENCE_REPORT.md before submitting 01+."
