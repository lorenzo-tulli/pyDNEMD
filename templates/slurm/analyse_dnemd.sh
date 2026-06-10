#!/bin/bash
#SBATCH --job-name=dnemd_analyse
#SBATCH --output=logs/analyse_dnemd_%a.out
#SBATCH --error=logs/analyse_dnemd_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --account=chem021482
#SBATCH --time=10:00:00
#SBATCH --partition=compute

module load openmpi/5.0.3
module load gromacs/2024.2-netlib-lapack

# task ID maps to the index in time_points_ps list in config.yaml
mkdir -p logs

dnemd-analyse-dnemd --config config_test.yaml --task-id $SLURM_ARRAY_TASK_ID
