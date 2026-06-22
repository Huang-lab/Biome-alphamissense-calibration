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

# --- pre-fetch Chen/Pejaver variant-level calibration table via wget ---------
# The python fetcher fails on Minerva login nodes whose outbound HTTPS allow-list
# excludes zenodo.org. wget is a simpler dependency and the file's a direct
# download (no API call). Skip if already extracted.
ZENODO_URL="$(python3 - <<'PY'
import os, yaml
cfg = os.environ['CONFIG_PATH']
with open(cfg) as fh: d = yaml.safe_load(fh)
print((d.get('calibration', {}) or {}).get('zenodo_file_url', '') or '')
PY
)"
if [[ -n "$ZENODO_URL" ]]; then
    if command -v wget >/dev/null 2>&1; then
        ZEN_DIR="$REPO_ROOT/refs/zenodo"
        mkdir -p "$ZEN_DIR"
        TGZ="$ZEN_DIR/calibration_download.tar.gz"
        # already-extracted CSV(s)?
        if ! ls "$ZEN_DIR"/*.csv >/dev/null 2>&1; then
            if [[ ! -s "$TGZ" ]]; then
                log "wget calibration tarball: $ZENODO_URL"
                wget --tries=3 --timeout=120 -O "$TGZ" "$ZENODO_URL" \
                    || die "wget failed for $ZENODO_URL; cluster may block outbound HTTPS — drop file in via calibration.local_file"
            fi
            log "extracting $TGZ -> $ZEN_DIR"
            tar -C "$ZEN_DIR" -xzf "$TGZ"
        fi
        CSV="$(find "$ZEN_DIR" -maxdepth 3 -type f -name '*.csv' | head -n 1)"
        [[ -s "$CSV" ]] || die "no .csv found under $ZEN_DIR after extraction"
        export BIOAM_CALIBRATION_CSV="$CSV"
        log "calibration CSV ready: $CSV"
    else
        warn "wget not on PATH; skipping Zenodo pre-fetch (will fall back to python requests or placeholder)"
    fi
fi

python3 "$REPO_ROOT/python/prepare_refs.py" --config "$CONFIG_PATH" "$@"

log "00_prepare_refs: done. Review refs/REFERENCE_REPORT.md before submitting 01+."
