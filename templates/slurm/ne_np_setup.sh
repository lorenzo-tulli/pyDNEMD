#!/bin/bash
#SBATCH --job-name=dnemd_ne_np_setup
#SBATCH --output=logs/ne_np_setup.out
#SBATCH --error=logs/ne_np_setup.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --partition=cpu

# ── environment ──────────────────────────────────────────────────────────────
source $(conda info --base)/etc/profile.d/conda.sh
conda activate pyDNEMD

module load GROMACS

# ── run ──────────────────────────────────────────────────────────────────────
# --start and --frequency override extract_start_ps and extract_frequency_ps
# from config.yaml if needed; otherwise remove them to use config defaults.
mkdir -p logs

dnemd-create-ne-np --config config.yaml --start 50000 --frequency 5000