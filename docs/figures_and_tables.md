# Figures and Tables — auditor reference

Each row maps a published panel or table to (a) the script that builds it,
(b) the data file it consumes, and (c) the reviewer Action it answers.
This is the index a future auditor uses to reproduce any figure from raw
data, and the index the rebuttal-letter author uses to cite the right
panel next to each reviewer quote.

For the reviewer-quote ↔ panel mapping, see
`python/figures/reviewer_map.yaml`.

## Pipeline files (regenerate everything)

| Stage | Script | Outputs |
|---|---|---|
| 00 | `scripts/00_prepare_refs.sh` | `refs/REFERENCE_REPORT.md`, `refs/chen_calibration.target_genes.tsv.gz` (BP4 + PP3 + every label), `refs/chen_summary_by_gene.tsv` |
| 01–04 | `scripts/0{1..4}_*.lsf` | `intermediate/<cohort>/...`, `results/<cohort>/{A,B,B_summary,C}.tsv` |
| 05 | `scripts/05_run_stats.lsf` (per cohort) | `results/<cohort>/stats_syndrome_associations.tsv` |
| 06 | `scripts/06_make_figures.lsf` | `results/figures/fig{1..6,S1..S7}.{png,pdf}` |

Wrapper: `scripts/submit_figures.sh` chains 05[cohortI] + 05[cohortII] (parallel) → 06.

## Figures

| Panel | Module | Reads | Reviewer Action |
|---|---|---|---|
| 1A | `python/figures/fig1_benchmark_honestly.py` `panel_1A` | ClinVar VCF (optional; Q3 pending) | Action 2 / R3 |
| 1B | same `panel_1B` | `refs/chen_calibration.target_genes.tsv.gz` BP4 rows | Action 2 / R2,R3 |
| 1C | same `panel_1C` | ClinVar VCF VUS subset (Q3 pending) | Action 1,2 / R3 |
| 1D | same `panel_1D` | `refs/chen_summary_by_gene.tsv` | Action 2 / R2 |
| 2A | `python/figures/fig2_carrier_landscape.py` `panel_2A` | `results/cohortI/C_regression_matrix.tsv` | Actions 1,2,3 setup / R1,R2 |
| 2B | same `panel_2B` | same | Actions 1,2 |
| 2C | same `panel_2C` | same (uses `carrier_<group>_AM0864`) | Action 2 / R2 |
| 3A | `python/figures/fig3_case_control_associations.py` `panel_3A` | `results/cohortI/stats_syndrome_associations.tsv` | Action 4 / R1,R2 (highest stakes) |
| 3B | same `panel_3B` | same (AM-only canonical cells) | Action 4 / R1 |
| 3C | same `panel_3C` | same (MUTYH-AP + MEN + Lynch) | Action 11 / R1 |
| 3D | same `panel_3D` | same (non-canonical, q_BH<0.1) | Action 5 |
| 4A | `python/figures/fig4_ancestry_stratified.py` `panel_4A` | `results/cohortI/C_regression_matrix.tsv` | Action 6 / R1 |
| 4B | same `panel_4B` (stub) | stratified run_stats output (not yet wired) | Action 6 / R1 |
| 4C | same `panel_4C` | C_regression_matrix.tsv PC columns | Actions 6,18 |
| 5A | `python/figures/fig5_cohort_replication.py` `panel_5A` | both cohorts' stats tables | Action 8 |
| 5B | same `panel_5B` | same | Action 8 |
| 6A | `python/figures/fig6_descriptive_landscape.py` `panel_6A` | C_regression_matrix.tsv | user-requested |
| 6B | same `panel_6B` | `B_ACMG_vs_AM_comparison.tsv` + C | user-requested |
| 6C | same `panel_6C` | C_regression_matrix.tsv | user-requested + Action 6 |
| 6D | same `panel_6D` | C | user-requested + Action 6 |
| 6E | same `panel_6E` | B + C | user-requested + Action 6 |
| S1–S7 | `python/figures/suppl.py` | (placeholders; see module docstring for each) | Various |

## Tables

| Table | Script | Status |
|---|---|---|
| Table 1 — per-gene Chen calibration map | `python/run_stats.py::build_table1` | pending — needs ClinVar stratified AUROC inputs |
| Table 2 — per-gene × per-cohort carrier counts | `python/run_stats.py::build_table2` | pending — descriptive, easy follow-up |
| Table 3 — case/control × variant-category × FH × age | `python/run_stats.py::build_table3` | pending — Action 3 |
| **Table 4 — syndrome-level OR (rebuttal cornerstone)** | `python/run_stats.py` (`stats_syndrome_associations.tsv`) | **built**; one row per (cohort, syndrome, phenotype, variant_category) |
| Table 5 — per-phenotype carrier-frequency rank | `python/run_stats.py::build_table5` | pending; backs Fig 6A |
| Table S8 / S10 — phenotype × gene long | `python/run_stats.py::build_tableS10` | pending; backs Fig 6B |
| Table S9 — phenotype demographics | `python/run_stats.py::build_tableS9` | pending; Action 10 |
| Table S11 — ancestry × phenotype long | `python/run_stats.py::build_tableS11` | pending; backs Figs 6C/D |
| Table S12 — ancestry × gene long | `python/run_stats.py::build_tableS12` | pending; backs Fig 6E |

## Sensitive-data note

Every figure here aggregates over ≥10 samples per cell when possible (HIPAA
safe-harbor floor). Figures and tables are written to `results/`, which is
gitignored end-to-end — only the generating *code* lives in the repo.
The pre-commit hook (`scripts/pre-commit.sh`) blocks any accidental
commit of patient-identifiable inputs.
