#!/bin/bash
#SBATCH --job-name=dnemd_eq
#SBATCH --output=logs/eq_%a.out
#SBATCH --error=logs/eq_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --array=1-5          # one task per replicate — adjust to n_runs in config.yaml

# ── environment ──────────────────────────────────────────────────────────────
source $(conda info --base)/etc/profile.d/conda.sh
conda activate pyDNEMD

module load GROMACS             # adjust to your HPC module name

# ── run ──────────────────────────────────────────────────────────────────────
mkdir -p logs

dnemd-equilibrium --config config.yaml --run $SLURM_ARRAY_TASK_ID