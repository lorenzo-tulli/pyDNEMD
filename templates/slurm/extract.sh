#!/bin/bash
#SBATCH --job-name=dnemd_extract
#SBATCH --output=logs/05_extract_%a.out
#SBATCH --error=logs/05_extract_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --account=XYZ
#SBATCH --time=10:00:00
#SBATCH --partition=compute

module load openmpi/5.0.3
module load gromacs/2024.2-netlib-lapack

########################################################################
## ---Resources requested to obtain the results in examples/output ---##
########################################################################

mkdir -p logs

dnemd-extract --config config_test.yaml --task-id $SLURM_ARRAY_TASK_ID
