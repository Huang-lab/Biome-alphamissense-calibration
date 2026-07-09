# BioMe AlphaMissense Calibration Pipeline

Code-only LSF pipeline (manual `bsub`) for applying gene-specific calibrated AlphaMissense thresholds [Chen/Pejaver 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13174790/#S2) to 28 cancer-predisposition genes in BioMe Cohort I (Regeneron) and Cohort II (Sema4), comparing against ACMG P/LP and ClinVar P/LP truth sets, and producing regression-ready tables.

## Variant categories (how each is defined)

Every carrier list, comparison table, and figure is built from these four
categories. All are matched at the **variant level** (`chrom, pos, ref, alt`,
with reference-free minimal-representation canonicalization so indel padding
differences still join):

| Category | Definition |
|---|---|
| **ACMG P/LP** | Carrier of a variant classified **Pathogenic / Likely-pathogenic by the ACMG interpretation** (ANNOVAR/InterVar), taken from the per-cohort `*_VariantsInSamplesPLPorPTV.tsv`: rows whose `annotation` contains `P/LP`. PTV-only rows are excluded (not comparable to AM missense). |
| **ClinVar P/LP** | Carrier of a variant that is **ClinVar Pathogenic / Likely-pathogenic with ≥ 2 gold stars** (review status `criteria_provided,_multiple_submitters,_no_conflicts` / `reviewed_by_expert_panel` / `practice_guideline`), excluding `Conflicting`. Covers **all variant types** (missense + indels/other), from a dedicated all-variant QC pass. |
| **AM_calibrated** | Carrier of a variant whose **Chen/Pejaver per-variant `evidence` label is ≥ PP3_Moderate** — i.e. in `{PP3_Moderate, PP3_Moderate+, PP3_Strong, PP3_Strong+, PP3_VeryStrong}`. This is a PP3 **evidence-label lookup**, *not* a numeric AlphaMissense-score cutoff: Chen already encodes the gene-specific (or domain-aggregate) calibration in that label. `BP4_*`, `PP3_Supporting`, and variants absent from the table are **not** carriers. Threshold set by `calibration.min_evidence_strength` (default `Moderate`). Missense SNVs only (AlphaMissense scores missense). |
| **AM_calibrated not P/LP** | An **AM_calibrated** carrier that is **neither ACMG P/LP nor ClinVar P/LP** — i.e. the novel AM calls no clinical source already flags. |

## Where to start

**The `refs/` directory is empty by design** — it is populated at runtime by
`scripts/00_prepare_refs.sh` and is never committed. Choose your path:

| I want to… | What to do |
|---|---|
| Verify the pipeline logic without any patient data | [Synthetic smoke test](#reproduce-without-minerva-no-phi-required) — runs on any laptop in ~30 s |
| Run the real BioMe cohorts | [Minerva full run](#run-order-manual-on-minerva) — start with step 00 below |

### Reference inputs (public download URLs)

The pipeline needs five external references. **Step 00 auto-downloads three of
them** (Chen calibration table, ClinVar VCF, and the Gencode GTF when needed);
the other **two you must pre-stage** and point at via `config/config.yaml`
(`references.alphamissense_tsv`, `references.gencode_transcript_map`). All are
public — you do **not** need Minerva access to obtain them.

| Reference | Public source | Auto-fetched by step 00? |
|---|---|---|
| AlphaMissense hg38 scores (`AlphaMissense_hg38.tsv.gz`) | [Zenodo — *Predictions for AlphaMissense* (record 8360242)](https://zenodo.org/records/8360242) | **No** — pre-stage; set `references.alphamissense_tsv` |
| Gencode v32 transcript→gene map (2 cols: `transcript_id`, `gene_name`) | derive from the [Gencode v32 GTF](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz) (one-liner below) | **No** — pre-stage; set `references.gencode_transcript_map` |
| Chen/Pejaver variant-level calibration table | [Zenodo record 18668684](https://zenodo.org/records/18668684) (`variant_level_calibration_table_AlphaMissense.csv.tar.gz`) | Yes (`wget`); fallback: `calibration.local_file` |
| Gencode v32 GTF (for the exon BED) | [EBI Gencode v32](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz) | Yes (`wget`); fallback: `references.gencode_gtf_local` |
| ClinVar VCF, GRCh38 (`clinvar.vcf.gz`) | [NCBI ClinVar](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz) | Yes (`wget`); fallback: `references.clinvar_vcf_local` |

> Huang-lab Minerva copies of the two pre-staged files live under
> `/sc/arion/projects/rg_huangk06/variants_PLP_BioMe/AlphaMissense/data/`.

Build the transcript→gene map from the Gencode GTF (needs `gawk`):

```bash
zcat gencode.v32.primary_assembly.annotation.gtf.gz \
  | gawk -F'\t' '$3=="transcript"{ \
      match($9,/transcript_id "([^"]+)"/,t); match($9,/gene_name "([^"]+)"/,g); \
      print t[1]"\t"g[1] }' \
  > gencode.v32.transcriptID_genename.tsv
```

Once the two pre-staged files are in place, **run step 00 first** (login node,
internet required):

```bash
bash scripts/00_prepare_refs.sh
```

This downloads the Chen/Pejaver calibration table and the ClinVar VCF from the
sources above, subsets the AM table to the 28 target genes, and writes everything
under `refs/`. It prints a `REFERENCE_REPORT.md` that you must review before
submitting steps 01+. Do not skip it — downstream steps will silently produce no
carriers if `refs/` is empty.

If the Minerva login node blocks outbound HTTPS (`zenodo.org`, `ncbi.nlm.nih.gov`),
download the files manually from the URLs above and set the corresponding
`*_local` / `local_file` config keys before running step 00.

> **ClinVar-only rerun:** if the AM/Chen references already exist and you only
> need to (re)build the ClinVar subset, use `bash scripts/00_prepare_refs.sh
> --only-clinvar` (skips the slow AM-TSV rescan).

## Reproduce without Minerva (no PHI required)

The real cohorts run only on Minerva (patient-identifiable inputs cannot be
shared). To validate the pipeline end-to-end on the bundled **synthetic**
fixtures — no Minerva, no PHI — clone the repo and run:

```bash
conda env create -f environment.yml && conda activate biome-am
bash tests/run_test.sh          # runs steps 01→04 in-process + asserts behaviors
```

This exercises QC, AM annotation, carrier calling, and the ACMG-vs-AM
comparison against `tests/data/` and prints `ALL TESTS PASSED` on success.
Everything below is the full Minerva workflow for the real cohorts.

## Cohort input files (Minerva paths)

All paths below are also encoded in `config/config.yaml` — that file is the single source of truth. This table is a quick reference.

### Cohort I (Regeneron WXS)
| Purpose | Path | Notes |
|---|---|---|
| Per-chr VCFs | `/sc/arion/projects/rg_huangk06/variants_PLP_BioMe/Regeneron/cancerChr/SINAI_Freeze_Two.GL.pVCF.PASS.onTarget.biallelic.{chr}.cancerGenes.vcf.gz` | `{chr}` = `chr1`..`chr22` (also `chrX` exists but no panel genes are on it, so we only iterate autosomes) |
| Phenotype TSV | `/sc/arion/projects/rg_huangk06/variants_PLP_BioMe/Regeneron/metadata/RegenWXS_HX_Newgroups.250109.tsv` | ID column: **`SINAI_ID`** |
| ACMG P/LP set | `/sc/arion/projects/rg_huangk06/variants_PLP_BioMe/AlphaMissense/2cohort/CarrFreq/Regen_VariantsInSamplesPLPorPTV.tsv` | PTV rows excluded (not comparable to AM-missense) |

### Cohort II (Sema4 WXS)
| Purpose | Path | Notes |
|---|---|---|
| Combined VCF | `/sc/private/regen/data/Sema4/BioMe_Sema4_WES.vcf.gz` | **READ-ONLY**; never write beside it |
| Local `.csi` for the combined VCF | `intermediate/cohortII/BioMe_Sema4_WES.vcf.gz.csi` | Written by `index_cohortII_source.lsf` (the discrete pre-step). 01 reads via `source##idx##local.csi` |
| Phenotype TSV | `/sc/arion/projects/rg_huangk06/variants_PLP_BioMe/metadata/Sema4WXS_HX_Newgroups.250109.tsv` | ID column: **`MASKED_MRN`** |
| ACMG P/LP set | `/sc/arion/projects/rg_huangk06/variants_PLP_BioMe/AlphaMissense/2cohort/CarrFreq/Sema4_VariantsInSamplesPLPorPTV.tsv` | Same schema as Cohort I |

### Shared inputs
| Purpose | Path | Notes |
|---|---|---|
| AlphaMissense scores | `/sc/arion/projects/rg_huangk06/variants_PLP_BioMe/AlphaMissense/data/AlphaMissense_hg38.tsv.gz` | Step 00 subsets to the 28 target genes |
| Gencode v32 transcript→gene map | `/sc/arion/projects/rg_huangk06/variants_PLP_BioMe/AlphaMissense/data/gencode.v32.transcriptID_genename.tsv` | Step 00 picks the canonical transcript per gene |
| Ancestry PCs | `/sc/arion/projects/rg_huangk06/variants_PLP_BioMe/GSA_GDA_PCA_V2.txt` | PC id column auto-discovered (max overlap with phenotype IDs) |
| Chen/Pejaver variant-level calibration table | `refs/zenodo/newAM_calibration_table_20260206.csv` | Step 00 wgets from Zenodo record 18668684 and untars |
| Repo location (this code) | `/sc/arion/projects/rg_huangk06/variants_PLP_BioMe/Biome-alphamissense-calibration` | Submit scripts compute this at runtime |

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
   Chains `01[1-22]` → `02[1-22]` → `02_gather` → `03` → `04` via `bsub -w "done(...)"`. Each step is submitted via a generated wrapper script that hardcodes `BIOAM_REPO_ROOT`, `COHORT`, and `CONFIG_PATH` as `export` statements inside the wrapper body — this is necessary because Minerva LSF does not propagate user env vars through `bsub -env` reliably.

4. **Cohort II (Sema4, single combined VCF in `/sc/private/`) — submit on a login node.**
   ```bash
   bash scripts/submit_cohortII.sh
   ```
   Chains `index_cohortII_source` → `01[1-22]` → `02[1-22]` → `02_gather` → `03` → `04`. The pre-step writes a `.csi` for the read-only source under `intermediate/cohortII/`; the 22 array tasks then `bcftools view -r chr<N>` against `source##idx##local.csi`.

### Smoke-test a single chromosome
For debugging or re-running a failed array task without re-submitting the whole pipeline:
```bash
bash scripts/submit_one_chr.sh cohortI 01_qc_missense 17
bash scripts/submit_one_chr.sh cohortI 02_annotate_am 17
```
Both `01_qc_missense` and `02_annotate_am` are array steps and can be re-run per chromosome. Other steps (`02_gather`, `03`, `04`) are whole-cohort and have no array index.

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

## Sensitive data + git hygiene

The phenotype TSVs, the masked-MRN bridge, the PCA file, and any Tables A/B/C
derived from real BioMe samples are **patient-identifiable** and must never
land in git. Two layers of defense:

1. `.gitignore` blocks the filename patterns directly (real-data file names plus
   `intermediate/` and `results/`). Synthetic fixtures under `tests/data/` are
   intentionally allowed.
2. **Install the pre-commit hook** in each clone before committing:
   ```bash
   ln -sf ../../scripts/pre-commit.sh .git/hooks/pre-commit
   chmod +x .git/hooks/pre-commit
   ```
   The hook rejects staged paths matching the real-data patterns and scans
   staged content for `SINAI_*` IDs and the MRN-map header. Run
   `bash scripts/pre-commit.sh` directly to test.

If a sensitive file ever slips in, halt and rewrite history before pushing.

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
