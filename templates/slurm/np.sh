#!/bin/bash
#SBATCH --job-name=dnemd_np
#SBATCH --output=logs/np_%a.out
#SBATCH --error=logs/np_%a.err
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH --partition=gpu,mwvdk
#SBATCH --array=0-4
#SBATCH --account=chem021482

module load openmpi/5.0.3
module load gromacs/2024.2-netlib-lapack

# task ID is mapped to (run_id, time_ns) automatically using
# extract_start_ps / extract_end_ps / extract_frequency_ps from config.yaml
mkdir -p logs

dnemd-run-np --config config_test.yaml --task-id $SLURM_ARRAY_TASK_ID
