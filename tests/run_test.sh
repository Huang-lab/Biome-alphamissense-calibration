#!/bin/bash
# =============================================================================
# tests/run_test.sh — local dry-run of stages 01→04 without LSF, using
# synthetic mini data. Verifies 8 behavioral cases (the 7 from the spec plus
# an AB-filter-syntax check) and the gVCF guard.
#
# Requires: bcftools, bgzip, tabix, python3 (with pyyaml, requests).
# Usage:    bash tests/run_test.sh
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/.." && pwd)"
RUN="$REPO_ROOT/tests/run"

source "$REPO_ROOT/lib/common.sh"

# tools check
for bin in bcftools bgzip tabix python3; do
    command -v "$bin" >/dev/null 2>&1 || die "missing required tool: $bin (install via conda env file environment.yml)"
done

rm -rf "$RUN"
mkdir -p "$RUN"/{refs,intermediate,results,logs,inputs}
mkdir -p "$RUN/inputs/cohortI" "$RUN/inputs/cohortII"

DATA="$REPO_ROOT/tests/data"

# ---- bgzip + index per-chr Cohort I VCFs -----------------------------------
log "preparing Cohort I per-chr VCFs"
for chr in chr1 chr2 chr3; do
    cp "$DATA/cohortI/${chr}.vcf" "$RUN/inputs/cohortI/${chr}.vcf"
    bgzip -f "$RUN/inputs/cohortI/${chr}.vcf"
    bcftools index --csi -f "$RUN/inputs/cohortI/${chr}.vcf.gz"
done
# empty per-chr files for chr4..chr22 (so gather_verify_22 doesn't fail)
for i in $(seq 4 22); do
    chr="chr${i}"
    cp "$DATA/cohortI/empty.vcf" "$RUN/inputs/cohortI/${chr}.vcf"
    bgzip -f "$RUN/inputs/cohortI/${chr}.vcf"
    bcftools index --csi -f "$RUN/inputs/cohortI/${chr}.vcf.gz"
done

# ---- bgzip Cohort II combined VCF ------------------------------------------
log "preparing Cohort II combined VCF"
cp "$DATA/cohortII/combined.vcf" "$RUN/inputs/cohortII/combined.vcf"
bgzip -f "$RUN/inputs/cohortII/combined.vcf"
# NOTE: we deliberately do NOT index the combined VCF here. The
# index_cohortII_source.lsf pre-step is responsible (testing that pre-step).

# ---- write test config -----------------------------------------------------
log "writing test config"
CONFIG="$RUN/test_config.yaml"
cat > "$CONFIG" <<EOF
paths:
  repo_root: "."
  refs_dir: "$RUN/refs"
  intermediate_dir: "$RUN/intermediate"
  results_dir: "$RUN/results"
  logs_dir: "$RUN/logs"

references:
  alphamissense_tsv: "$RUN/inputs/mini_AM.tsv.gz"
  gencode_transcript_map: "$DATA/mini_gencode_map.tsv"
  gencode_gtf_url: ""
  gencode_gtf_local: "$DATA/mini_gencode.gtf"

calibration:
  zenodo_record: "TEST"
  zenodo_api_url: ""
  local_file: "$DATA/mini_calibration.tsv"
  primary_threshold: "gene_specific"

target_genes:
  - APC
  - BMPR1A
  - BRCA1
  - BRCA2
  - MAX
  - MEN1
  - MLH1
  - MSH2
  - MSH6
  - NF2
  - PALB2
  - PMS2
  - PTEN
  - RB1
  - RET
  - SDHAF2
  - SDHB
  - SDHC
  - SDHD
  - SMAD4
  - STK11
  - TMEM127
  - TP53
  - TSC1
  - TSC2
  - VHL
  - WT1
  - MUTYH

gene_groups:
  Hereditary_Breast_and_Ovarian_Cancer_Syndrome: [BRCA1, BRCA2, PALB2]
  Lynch_Syndrome: [MLH1, MSH2, MSH6, PMS2]
  Familial_Adenomatous_Polyposis: [APC]
  MUTYH_Associated_Polyposis: [MUTYH]
  Multiple_Endocrine_Neoplasia: [MEN1, RET]
  Hereditary_Paraganglioma_Pheochromocytoma: [SDHAF2, SDHB, SDHC, SDHD, TMEM127, MAX]
  Li_Fraumeni_Syndrome: [TP53]
  Tuberous_Sclerosis_Complex: [TSC1, TSC2]
  Von_Hippel_Lindau_Syndrome: [VHL]
  PTEN_Hamartoma_Tumor_Syndrome: [PTEN]
  Juvenile_Polyposis_Syndrome: [BMPR1A, SMAD4]
  Hereditary_Retinoblastoma: [RB1]
  Peutz_Jeghers_Syndrome: [STK11]
  Neurofibromatosis_Type_2: [NF2]
  Wilms_Tumor_Syndrome: [WT1]

qc:
  filter_pass_only: true
  min_dp: 10
  min_gq: 20
  ab_lo: 0.2
  ab_hi: 0.8
  drop_homref: true
  drop_missing: true

cohorts:
  cohortI:
    name: "cohortI"
    mode: "per_chr_input"
    vcf_dir: "$RUN/inputs/cohortI"
    vcf_pattern: "{chr}.vcf.gz"
    phenotype_tsv: "$DATA/mini_phenotype.tsv"
    id_column: "SINAI_ID"
    acmg_tsv: "$DATA/mini_acmg_cohortI.tsv"
  cohortII:
    name: "cohortII"
    mode: "sliced_from_combined"
    combined_vcf: "$RUN/inputs/cohortII/combined.vcf.gz"
    local_source_dir: "$RUN/intermediate/cohortII"
    phenotype_tsv: "$DATA/mini_phenotype.tsv"
    id_column: "MASKED_MRN"
    acmg_tsv: "$DATA/mini_acmg_cohortII.tsv"

pcs:
  pcs_tsv: "$DATA/mini_pcs.tsv"
  id_column: ""
  pc_columns: ["PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"]

acmg:
  annotation_filter_substr: "P/LP"
  annotation_column: "annotation"

lsf:
  project: "acc_test"
  queue: "test"
  step_00_prepare_refs:  { n: 2, W: "00:30", mem: "4000" }
  step_01_qc_missense:   { n: 2, W: "00:30", mem: "4000" }
  step_02_annotate_am:   { n: 2, W: "00:30", mem: "4000" }
  step_02_gather:        { n: 2, W: "00:30", mem: "4000" }
  step_03_call_carriers: { n: 2, W: "00:30", mem: "4000" }
  step_04_compare_tab:   { n: 2, W: "00:30", mem: "4000" }
  step_index_cohortII:   { n: 2, W: "00:30", mem: "4000" }
EOF
export CONFIG_PATH="$CONFIG"

# ---- bgzip the mini AM (the config points at .gz; do this BEFORE 00) -------
log "preparing mini AM (bgzip+tabix happens inside 00)"
cp "$DATA/mini_AM.tsv" "$RUN/inputs/mini_AM.tsv"
bgzip -f "$RUN/inputs/mini_AM.tsv"
# NOTE: prepare_refs.subset_am will write a fresh subset+tabix; the raw mini
# AM .gz doesn't need a tabix index.

# ---- 00 prepare refs (skip network; uses local_file and local GTF) ---------
log "==== 00_prepare_refs ===="
bash "$REPO_ROOT/scripts/00_prepare_refs.sh"

# ---- 01 + 02 per cohort, manual array loop --------------------------------
for COHORT in cohortI cohortII; do
    export COHORT
    if [[ "$COHORT" == "cohortII" ]]; then
        log "==== index_cohortII_source ===="
        bash "$REPO_ROOT/scripts/index_cohortII_source.lsf"
    fi
    for i in $(seq 1 22); do
        export TEST_CHR_IDX=$i
        log "==== 01_qc_missense  cohort=$COHORT chr$i ===="
        bash "$REPO_ROOT/scripts/01_qc_missense.lsf"
    done
    for i in $(seq 1 22); do
        export TEST_CHR_IDX=$i
        log "==== 02_annotate_am  cohort=$COHORT chr$i ===="
        bash "$REPO_ROOT/scripts/02_annotate_am.lsf"
    done
    unset TEST_CHR_IDX
    log "==== 02_gather  cohort=$COHORT ===="
    bash "$REPO_ROOT/scripts/02_gather.lsf"
    log "==== 03_call_carriers  cohort=$COHORT ===="
    bash "$REPO_ROOT/scripts/03_call_carriers.lsf"
    log "==== 04_compare_and_tabulate  cohort=$COHORT ===="
    bash "$REPO_ROOT/scripts/04_compare_and_tabulate.lsf"
done

# ---- gVCF guard test (separate; expects 01 to ABORT) -----------------------
log "==== gVCF guard test (expects 01 to abort) ===="
cp "$DATA/gvcf_sample.vcf" "$RUN/inputs/gvcf_sample.vcf"
bgzip -f "$RUN/inputs/gvcf_sample.vcf"
bcftools index --csi -f "$RUN/inputs/gvcf_sample.vcf.gz"

GVCFCFG="$RUN/gvcf_test_config.yaml"
sed "s|$RUN/inputs/cohortI|$RUN/inputs/gvcf_only|g" "$CONFIG" > "$GVCFCFG"
mkdir -p "$RUN/inputs/gvcf_only"
cp "$RUN/inputs/gvcf_sample.vcf.gz" "$RUN/inputs/gvcf_only/chr1.vcf.gz"
cp "$RUN/inputs/gvcf_sample.vcf.gz.csi" "$RUN/inputs/gvcf_only/chr1.vcf.gz.csi"

# `set +e` does NOT suppress the inherited ERR trap from lib/common.sh, so we
# use `|| rc=$?` (which is excluded from the ERR-trap conditions) to capture
# the expected non-zero exit without aborting this script.
GVCF_RC=0
CONFIG_PATH="$GVCFCFG" COHORT=cohortI TEST_CHR_IDX=1 \
    bash "$REPO_ROOT/scripts/01_qc_missense.lsf" 2>"$RUN/gvcf_stderr.log" \
    || GVCF_RC=$?
if [[ $GVCF_RC -eq 0 ]]; then
    cat "$RUN/gvcf_stderr.log" >&2
    die "gVCF guard FAILED to abort step 01 (expected non-zero exit)"
fi
if ! grep -q "Input appears to be a combined gVCF" "$RUN/gvcf_stderr.log"; then
    cat "$RUN/gvcf_stderr.log" >&2
    die "gVCF guard aborted, but stderr did not contain the prescribed message"
fi
log "gVCF guard: OK (exited $GVCF_RC with the prescribed message)"

# ---- behavioral assertions -------------------------------------------------
log "==== behavioral assertions ===="
python3 "$REPO_ROOT/tests/assert_behaviors.py" --run-dir "$RUN" --config "$CONFIG"

log "ALL TESTS PASSED"
