#!/bin/bash
#SBATCH --job-name=dnemd_ne_np_setup
#SBATCH --output=logs/ne_np_setup_%a.out
#SBATCH --error=logs/ne_np_setup_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --account=chem021482
#SBATCH --time=1:00:00
#SBATCH --partition=test
#SBATCH --array=1-2          # adjust to n_runs in config

module load openmpi/5.0.3
module load gromacs/2024.2-netlib-lapack

# Each array task processes one replicate.
# --skip-topology-test avoids race conditions when tasks run simultaneously.
mkdir -p logs

dnemd-create-ne-np --config config_test.yaml --run $SLURM_ARRAY_TASK_ID --skip-topology-test
