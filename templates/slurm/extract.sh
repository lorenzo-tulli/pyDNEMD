#!/bin/bash
#SBATCH --job-name=dnemd_extract
#SBATCH --output=logs/extract_%a.out
#SBATCH --error=logs/extract_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --account=chem021482
#SBATCH --time=10:00:00
#SBATCH --partition=compute

module load openmpi/5.0.3
module load gromacs/2024.2-netlib-lapack

# task ID maps to the corresponding ns window defined by
# extract_ns_start / extract_ns_end / extract_ns_interval in config.yaml
mkdir -p logs

dnemd-extract --config config_test.yaml --task-id $SLURM_ARRAY_TASK_ID
