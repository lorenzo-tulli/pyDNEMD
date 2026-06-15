#!/bin/bash
#SBATCH --job-name=dnemd_eq
#SBATCH --output=logs/01_eq_%a.out
#SBATCH --error=logs/01_eq_%a.err
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH --partition=gpu,mwvdk
#SBATCH --array=1-2
#SBATCH --account=XYZ

module load openmpi/5.0.3
module load gromacs/2024.2-netlib-lapack

########################################################################
## ---Resources requested to obtain the results in examples/output ---##
########################################################################

mkdir -p logs

dnemd-equilibrium --config config_test.yaml --run $SLURM_ARRAY_TASK_ID
