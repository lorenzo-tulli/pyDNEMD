#!/bin/bash
#SBATCH --job-name=dnemd_np
#SBATCH --output=logs/np_%a.out
#SBATCH --error=logs/np_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --partition=cpu
#SBATCH --array=0-199        # adjust: n_runs * n_timepoints - 1

# ── environment ──────────────────────────────────────────────────────────────
source $(conda info --base)/etc/profile.d/conda.sh
conda activate pyDNEMD

module load GROMACS

# ── run ──────────────────────────────────────────────────────────────────────
# task ID is mapped to (run_id, time_ns) automatically using
# extract_start_ps / extract_end_ps / extract_frequency_ps from config.yaml
mkdir -p logs

dnemd-run-np --config config.yaml --task-id $SLURM_ARRAY_TASK_ID