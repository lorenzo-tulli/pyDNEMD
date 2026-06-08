#!/bin/bash
#SBATCH --job-name=dnemd_extract
#SBATCH --output=logs/extract_%a.out
#SBATCH --error=logs/extract_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --partition=cpu
#SBATCH --array=0-39         # adjust: number of ns windows - 1
                              # ns windows = (extract_ns_end - extract_ns_start)
                              #              / extract_ns_interval

# ── environment ──────────────────────────────────────────────────────────────
source $(conda info --base)/etc/profile.d/conda.sh
conda activate pyDNEMD

module load GROMACS

# ── run ──────────────────────────────────────────────────────────────────────
# task ID maps to the corresponding ns window defined by
# extract_ns_start / extract_ns_end / extract_ns_interval in config.yaml
mkdir -p logs

dnemd-extract --config config.yaml --task-id $SLURM_ARRAY_TASK_ID