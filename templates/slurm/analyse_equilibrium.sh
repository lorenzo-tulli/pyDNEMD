#!/bin/bash
#SBATCH --job-name=dnemd_analyse_eq
#SBATCH --output=logs/analyse_eq.out
#SBATCH --error=logs/analyse_eq.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --partition=cpu

# ── environment ──────────────────────────────────────────────────────────────
source $(conda info --base)/etc/profile.d/conda.sh
conda activate pyDNEMD

module load GROMACS             # needed for gmx rms / rmsf calls

# ── run ──────────────────────────────────────────────────────────────────────
mkdir -p logs

dnemd-analyse-equilibrium --config config.yaml