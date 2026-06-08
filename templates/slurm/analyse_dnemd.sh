#!/bin/bash
#SBATCH --job-name=dnemd_analyse
#SBATCH --output=logs/analyse_dnemd_%a.out
#SBATCH --error=logs/analyse_dnemd_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --partition=cpu
#SBATCH --array=0-4          # adjust: number of time_points_ps entries - 1
                              # e.g. [0, 10, 100, 1000, 5000] → array=0-4

# ── environment ──────────────────────────────────────────────────────────────
source $(conda info --base)/etc/profile.d/conda.sh
conda activate pyDNEMD

# ── run ──────────────────────────────────────────────────────────────────────
# task ID maps to the index in time_points_ps list in config.yaml
mkdir -p logs

dnemd-analyse-dnemd --config config.yaml --task-id $SLURM_ARRAY_TASK_ID