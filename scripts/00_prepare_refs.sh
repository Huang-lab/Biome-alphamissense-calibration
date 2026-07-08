#!/bin/bash
# =============================================================================
# 00_prepare_refs.sh — LOGIN NODE script (internet required), no LSF.
#
# Builds:
#   refs/chen_calibration.target_genes.tsv.gz (+ .tbi)   <- per-variant Chen labels
#   refs/chen_summary_by_gene.tsv                        <- per-gene Chen coverage
#   refs/gene_transcript_map.tsv
#   refs/target_genes.exons.bed
#   refs/AlphaMissense_hg38.subset_targets.tsv.gz (+ .tbi)
#   refs/REFERENCE_REPORT.md
# Hard-gates: prints Chen coverage counts and tells you to review
# REFERENCE_REPORT.md before submitting 01+.
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

# --only-clinvar fast path: (re)build ONLY the ClinVar subset, skipping the slow
# Chen/AM/gencode ref rebuild. Detect it so we can also skip the Zenodo fetch.
ONLY_CLINVAR=0
for _a in "$@"; do [[ "$_a" == "--only-clinvar" ]] && ONLY_CLINVAR=1; done

# requirements check
for bin in python3 bgzip tabix; do
    command -v "$bin" >/dev/null 2>&1 || die "required binary not on PATH: $bin"
done

# --- pre-fetch Chen/Pejaver variant-level calibration table via wget ---------
# The python fetcher fails on Minerva login nodes whose outbound HTTPS allow-list
# excludes zenodo.org. wget is a simpler dependency and the file's a direct
# download (no API call). Skip if already extracted, or entirely in --only-clinvar.
if [[ "$ONLY_CLINVAR" -eq 1 ]]; then
    log "--only-clinvar: skipping Chen/Zenodo fetch and AM/gencode ref rebuild"
fi
if [[ "$ONLY_CLINVAR" -eq 0 ]]; then
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
fi  # end: skip Chen/Zenodo fetch under --only-clinvar

# --- pre-fetch ClinVar VCF (GRCh38) via wget --------------------------------
# prepare_refs.subset_clinvar() reads it from references.clinvar_vcf_local, or
# from $BIOAM_CLINVAR_VCF which we export here after downloading. Skip if a
# local path is already configured, or if we've downloaded it before.
CLINVAR_LOCAL="$(cfg_get references.clinvar_vcf_local 2>/dev/null || echo '')"
if [[ -z "$CLINVAR_LOCAL" ]]; then
    CLINVAR_URL="$(cfg_get references.clinvar_vcf_url 2>/dev/null || echo '')"
    if [[ -n "$CLINVAR_URL" ]] && command -v wget >/dev/null 2>&1; then
        CV_DIR="$REPO_ROOT/refs/clinvar"
        mkdir -p "$CV_DIR"
        CV_VCF="$CV_DIR/clinvar.vcf.gz"
        if [[ ! -s "$CV_VCF" ]]; then
            log "wget ClinVar VCF: $CLINVAR_URL"
            wget --tries=3 --timeout=180 -O "$CV_VCF" "$CLINVAR_URL" \
                || die "wget failed for $CLINVAR_URL; drop the file in via references.clinvar_vcf_local"
        fi
        export BIOAM_CLINVAR_VCF="$CV_VCF"
        log "ClinVar VCF ready: $CV_VCF"
    else
        warn "no ClinVar URL configured or wget unavailable; ClinVar_PLP category will be empty"
    fi
fi

# --- left-normalize ClinVar against the reference FASTA (for indel matching) --
# subset_clinvar reads BIOAM_CLINVAR_VCF; if we normalize, repoint it at the
# normalized file so the target-gene subset (and thus indel coordinates) match
# the normalized BioMe all-variant QC VCF. Requires bcftools + reference_fasta.
REF_FASTA="$(cfg_get references.reference_fasta 2>/dev/null || echo '')"
SRC_CLINVAR="${BIOAM_CLINVAR_VCF:-$CLINVAR_LOCAL}"
if [[ -n "$REF_FASTA" && -s "$REF_FASTA" && -n "$SRC_CLINVAR" && -s "$SRC_CLINVAR" ]]; then
    if command -v bcftools >/dev/null 2>&1; then
        CV_NORM="$REPO_ROOT/refs/clinvar/clinvar.norm.vcf.gz"
        mkdir -p "$(dirname "$CV_NORM")"
        log "left-normalizing ClinVar against $REF_FASTA -> $CV_NORM"
        bcftools norm -f "$REF_FASTA" -m -any -Oz -o "$CV_NORM" "$SRC_CLINVAR" \
            || die "bcftools norm failed on ClinVar VCF ($SRC_CLINVAR)"
        export BIOAM_CLINVAR_VCF="$CV_NORM"
    else
        warn "bcftools not on PATH; cannot normalize ClinVar — indel matching will be best-effort"
    fi
elif [[ -z "$REF_FASTA" ]]; then
    warn "references.reference_fasta unset; ClinVar indels matched best-effort (SNVs exact). Set it in config.local.yaml for reliable indel matching."
fi

python3 "$REPO_ROOT/python/prepare_refs.py" --config "$CONFIG_PATH" "$@"

log "00_prepare_refs: done. Review refs/REFERENCE_REPORT.md before submitting 01+."
