# BioMe AlphaMissense Calibration Pipeline

Code-only LSF pipeline (manual `bsub`) for applying gene-specific calibrated AlphaMissense thresholds (Chen/Pejaver 2026) to 28 cancer-predisposition genes in BioMe Cohort I (Regeneron) and Cohort II (Sema4), comparing against an existing ACMG P/LP set, and producing regression-ready tables.

## Run order (manual on Minerva)

1. **Prepare references — login node, internet required.**
   ```bash
   bash scripts/00_prepare_refs.sh
   ```
   Builds the calibration table, gene→transcript map, exon BED, AM subset (bgzip+tabix), and writes `refs/REFERENCE_REPORT.md`. This is a **hard gate** — review the 28-gene table and coverage counts before submitting anything else. If the Zenodo fetch failed, drop the published calibration TSV at the path you set in `config.calibration.local_file` and rerun. Likewise for the gencode GTF if the exon-BED source falls back to `AM_minmax_fallback`.

2. **Local behavioral tests — no LSF.**
   ```bash
   bash tests/run_test.sh
   ```
   Runs steps 01→04 in-process on a synthetic mini dataset, asserts 8 behavioral cases (the 7 from the spec plus an AB-filter-syntax check), and confirms the gVCF guard aborts with the prescribed message.

3. **Cohort I (Regeneron, per-chr VCFs) — submit on a login node.**
   ```bash
   bash scripts/submit_cohortI.sh
   ```
   Chains `01[1-22]` → `02[1-22]` → `02_gather` → `03` → `04` via `bsub -w "done(...)"`.

4. **Cohort II (Sema4, single combined VCF in `/sc/private/`) — submit on a login node.**
   ```bash
   bash scripts/submit_cohortII.sh
   ```
   Chains `index_cohortII_source` → `01[1-22]` → `02[1-22]` → `02_gather` → `03` → `04`. The pre-step writes a `.csi` for the read-only source under `intermediate/cohortII/`; the 22 array tasks then `bcftools view -r chr<N>` against `source##idx##local.csi`.

5. **Collect deliverables.**
   - `results/<cohort>/A_variant_level_per_person.tsv`
   - `results/<cohort>/B_ACMG_vs_AM_comparison.tsv` + `B_summary_by_gene.tsv`
   - `results/<cohort>/C_regression_matrix.tsv`
   - `results/SUMMARY.md` (per-cohort match rates, carrier counts per gene group, gene-coverage table)

## How the pieces fit together

| stage | script | input | output |
|------|--------|-------|--------|
| 00   | `scripts/00_prepare_refs.sh`         | AM TSV, gencode map, Zenodo calibration, gencode v32 GTF | `refs/` artifacts + `REFERENCE_REPORT.md` |
| 01   | `scripts/01_qc_missense.lsf` (array) | per-chr VCFs *(Cohort I)* OR sliced combined VCF *(Cohort II)* | `intermediate/<cohort>/chr<N>.qc.vcf.gz` |
| 02   | `scripts/02_annotate_am.lsf` (array) | per-chr filtered VCF, AM subset | `intermediate/<cohort>/chr<N>.am_annotated.tsv` |
| 02g  | `scripts/02_gather.lsf` (one-shot)   | the 22 per-chr TSVs              | `intermediate/<cohort>/am_annotated.tsv` (single header) |
| 03   | `scripts/03_call_carriers.lsf`       | gathered annotated TSV + thresholds | `intermediate/<cohort>/A_variant_level_per_person.tsv` |
| 04   | `scripts/04_compare_and_tabulate.lsf`| Table A, phenotype, PCs, ACMG P/LP | `results/<cohort>/{A,B,B_summary,C}.tsv`, `results/SUMMARY.md` |
| pre  | `scripts/index_cohortII_source.lsf`  | `/sc/private/.../Sema4.vcf.gz`   | `intermediate/cohortII/<basename>.csi` |

## Locked design decisions (enforced in code)

- **Threshold policy.** Primary call uses **gene-specific calibrated thresholds only**. Domain-aggregate thresholds are RECORDED in the output but NEVER used to promote a carrier. Genes with neither calibration are retained in Table A with their `am_pathogenicity` and `threshold_source=uncovered`, but excluded from the AM columns of Table C. **No global 0.864 fallback anywhere.** Change `calibration.primary_threshold` in `config/config.yaml` if (and only if) the policy needs to change — there is no other knob.
- **Missense scope = AM-join.** A variant is in-scope iff it matches an AM row for one of the 28 target genes' transcripts. No separate consequence annotator (VEP/SnpEff/...) is invoked.
- **SNV-only.** Multiallelic records cause a hard abort; AM-matched variants must be single-base SNVs (also enforced).
- **No AF filter.** The VCFs have no site AF field and the design explicitly forbids one.
- **QC.** `FILTER==PASS` site filter; per-genotype `DP≥10`, `GQ≥20`, het AB in `[0.2, 0.8]`; hom-ref and missing dropped. All knobs in `config.yaml > qc.*`. The filter expression itself is centralized in `lib/common.sh::build_qc_expr` (single source of truth).
- **gVCF guard.** Aborts with the exact message
  `Input appears to be a combined gVCF, not genotyped calls; point at the GenotypeGVCFs output.`
- **gather-and-verify gate.** Step 03 independently re-verifies all 22 per-chr outputs even when submitted with `bsub -w "done(...)"` — array tasks can silently fail.
- **Read-only Cohort II source.** Nothing is ever written under `/sc/private/`. The source is indexed via a **discrete pre-step** (NOT inside the array) so tasks 2–22 don't race the index.
- **ACMG comparison.** Default filter: `annotation` column contains `"P/LP"`. PTV rows are excluded with a warning that PTV is not comparable to AM-missense (`config.acmg.annotation_filter_substr`).
- **Per-cohort ID columns.** Cohort I links phenotype on `SINAI_ID`; Cohort II on `MASKED_MRN`. Configurable; match rate is reported and low rates trigger a warning that prints up to 5 unmatched IDs from each side.

## Where to change things

- **Paths / config**: `config/config.yaml` — everything is there. Do not hard-code paths in scripts.
- **QC parameters**: `config.yaml > qc:`. Reflected automatically into `lib/common.sh::build_qc_expr`.
- **Module names** (Minerva-specific): `lib/load_modules.sh`. Lines are marked `VERIFY` because Minerva module names are site-specific.
- **bcftools filter syntax** (per-version differences): `lib/common.sh::build_qc_expr`. The script's banner prints `bcftools --version` and pokes the expression at job start so a mismatch fails fast.
- **Gene panel / Table S5 groups**: `config.yaml > target_genes` and `gene_groups`.

## Hardware / environment

- **CPU only.** GPU not used: this is `bcftools` streaming + pandas joins; no model inference, no GPU-amenable workload here. AM scores are precomputed in the public table.
- **Modules expected on Minerva (VERIFY against your site names):** `bcftools`, `htslib` (provides `bgzip`+`tabix`), `python` (≥3.9 with `pandas`, `pyyaml`, `requests`).
- **Conda fallback** when modules aren't available: `conda env create -f environment.yml && conda activate biome-am`.

## Assumptions and limits (called out in code with `VERIFY`/`WARN`)

- VCFs are GenotypeGVCFs output (biallelic, real GTs). The gVCF guard catches the most common slip.
- AM transcript IDs match gencode v32 transcript IDs by exact versioned ID (with a version-stripped fallback). Any gene whose AM transcript can't be located is listed in `REFERENCE_REPORT.md`.
- PCs are joined on a column auto-discovered by maximum overlap with phenotype IDs; can be overridden via `config.yaml > pcs.id_column`.
- Calibration thresholds are fetched/parsed heuristically; if the Zenodo file format changes, `python/prepare_refs.py::_try_parse_calibration_table` is the only spot to update.

## Feedback / iterating

The pipeline is intentionally stage-wise and idempotent. To rerun a single chromosome: delete its `intermediate/<cohort>/chr<N>.qc.vcf.gz` (and the matching `.am_annotated.tsv`) and `bsub` only that array index (e.g. `bsub -J biome_am_01_cohortI[5]` after editing the array spec, or run the body directly with `LSB_JOBINDEX=5 COHORT=cohortI bash scripts/01_qc_missense.lsf`). To rerun a whole cohort, delete `intermediate/<cohort>/` and rerun `submit_<cohort>.sh`.
