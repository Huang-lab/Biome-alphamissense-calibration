#!/bin/bash
# Submit the rebuttal-figures chain on Minerva:
#   05_run_stats[cohortI]  + 05_run_stats[cohortII]   (parallel)
#       -> 06_make_figures   (cross-cohort, blocks on both)
#
# Run on a login node from the repo root. No arguments.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/.." && pwd)"
cd "$REPO_ROOT"

# Reuse the env-export trick from submit_cohortI/II.sh: hardcode the values
# inside a wrapper script so LSF picks them up reliably (bsub -env propagation
# is flaky on Minerva).
PROJ="$(python3 -c "import yaml; print(yaml.safe_load(open('config/config.yaml'))['lsf']['project'])")"
QUEUE="$(python3 -c "import yaml; print(yaml.safe_load(open('config/config.yaml'))['lsf']['queue'])")"

submit_step() {
    local cohort="$1"
    local step_script="$2"
    local jobname="$3"
    local wrapper="$REPO_ROOT/intermediate/wrappers/${jobname}.sh"
    mkdir -p "$(dirname "$wrapper")"
    cat > "$wrapper" <<EOF
#!/bin/bash
export BIOAM_REPO_ROOT="$REPO_ROOT"
export CONFIG_PATH="$REPO_ROOT/config/config.yaml"
export COHORT="$cohort"
exec bash "$REPO_ROOT/scripts/${step_script}.lsf"
EOF
    chmod +x "$wrapper"
    bsub -J "$jobname" \
         -P "$PROJ" -q "$QUEUE" \
         -n 4 -W 02:00 -R rusage[mem=16000] -R span[hosts=1] \
         -o "logs/${jobname}.%J.stdout" -eo "logs/${jobname}.%J.stderr" \
         -L /bin/bash \
         < "$wrapper"
}

submit_step cohortI  05_run_stats biome_am_05_cohortI
submit_step cohortII 05_run_stats biome_am_05_cohortII

# Cross-cohort: depend on both 05 jobs.
WRAPPER="$REPO_ROOT/intermediate/wrappers/biome_am_06_figures.sh"
mkdir -p "$(dirname "$WRAPPER")"
cat > "$WRAPPER" <<EOF
#!/bin/bash
export BIOAM_REPO_ROOT="$REPO_ROOT"
export CONFIG_PATH="$REPO_ROOT/config/config.yaml"
exec bash "$REPO_ROOT/scripts/06_make_figures.lsf"
EOF
chmod +x "$WRAPPER"
bsub -J biome_am_06_figures \
     -w "done(biome_am_05_cohortI) && done(biome_am_05_cohortII)" \
     -P "$PROJ" -q "$QUEUE" \
     -n 4 -W 02:00 -R rusage[mem=16000] -R span[hosts=1] \
     -o "logs/biome_am_06_figures.%J.stdout" \
     -eo "logs/biome_am_06_figures.%J.stderr" \
     -L /bin/bash \
     < "$WRAPPER"

echo "submitted: 05[cohortI] + 05[cohortII] (parallel) -> 06_figures"
echo "monitor: bjobs -J 'biome_am_*'"
