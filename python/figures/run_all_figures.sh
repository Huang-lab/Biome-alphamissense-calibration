#!/usr/bin/env bash
# Run all figure scripts in the correct order.
# Must be run from the repo root directory.
set -e

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO_ROOT"

export MPLCONFIGDIR="${MPLCONFIGDIR:-/private/tmp/biome_mplconfig}"
mkdir -p "$MPLCONFIGDIR"

echo "=== BioMe AlphaMissense Rebuttal Figures ==="
echo "Repo root: $REPO_ROOT"
echo ""

# Step 1: Run regressions (skip per-gene for speed; run separately if needed)
echo "--- Step 1: run_stats.py (cohortI) ---"
python3 python/figures/run_stats.py --cohort cohortI --skip-gene

echo "--- Step 1b: run_stats.py (cohortII) ---"
python3 python/figures/run_stats.py --cohort cohortII --skip-gene

# Step 2: Build existing figures (each now saves individual panels + data tables)
echo ""
echo "--- Step 2: Building figures (individual panels + composites) ---"

echo "  Figure 3 (associations — credibility anchor) ..."
python3 python/figures/fig3_associations.py

echo "  Figure 1 (threshold motivation) ..."
python3 python/figures/fig1_threshold.py --cohort cohortI
python3 python/figures/fig1_threshold.py --cohort cohortII

echo "  Figure 2 (carrier landscape) ..."
python3 python/figures/fig2_carrier_landscape.py

echo "  Figure 4 (case-control) ..."
python3 python/figures/fig4_casecontrol.py --cohort cohortI
python3 python/figures/fig4_casecontrol.py --cohort cohortII

echo "  Figure 5 (replication + ancestry) ..."
python3 python/figures/fig5_replication.py

echo "  Supplementary figures ..."
python3 python/figures/suppl.py

# Step 3: New figures
echo ""
echo "--- Step 3: New gene × phenotype heatmaps ---"
python3 python/figures/fig_gene_pheno_heatmap.py --cohort cohortI
python3 python/figures/fig_gene_pheno_heatmap.py --cohort cohortII

echo "--- Step 3b: Phenotype carrier frequency barplot ---"
python3 python/figures/fig_pheno_barplot.py --cohort cohortI
python3 python/figures/fig_pheno_barplot.py --cohort cohortII

# Step 4: Figure guide
echo ""
echo "--- Step 4: Figure guide document ---"
python3 python/figures/write_figure_guide.py

echo ""
echo "=== All figures complete ==="
echo "Individual panels:  ${BIOME_FIGS_DIR:-$REPO_ROOT/results/figures}"
echo "Data tables:        ${BIOME_TABLES_DIR:-$REPO_ROOT/results/tables}"
echo "Figure guide:       ${BIOME_FIGS_DIR:-$REPO_ROOT/results/figures}/figure_guide.md"
echo ""
find "${BIOME_FIGS_DIR:-$REPO_ROOT/results/figures}" -maxdepth 1 -type f -name '*.png' -print | sort
