#!/bin/bash
#SBATCH --job-name=dnemd_eq
#SBATCH --output=logs/eq_%a.out
#SBATCH --error=logs/eq_%a.err
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH --partition=gpu,mwvdk
#SBATCH --array=1-2
#SBATCH --account=chem021482

module load openmpi/5.0.3
module load gromacs/2024.2-netlib-lapack

mkdir -p logs

dnemd-equilibrium --config config_test.yaml --run $SLURM_ARRAY_TASK_ID
