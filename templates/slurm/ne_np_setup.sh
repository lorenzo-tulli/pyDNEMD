#!/bin/bash
#SBATCH --job-name=dnemd_ne_np_setup
#SBATCH --output=logs/03_ne_np_setup_%a.out
#SBATCH --error=logs/03_ne_np_setup_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --account=XYZ
#SBATCH --time=1:00:00
#SBATCH --partition=test
#SBATCH --array=1-2

module load openmpi/5.0.3
module load gromacs/2024.2-netlib-lapack

########################################################################
## ---Resources requested to obtain the results in examples/output ---##
########################################################################

# --skip-topology-test avoids race conditions when tasks run simultaneously.
mkdir -p logs

dnemd-create-ne-np --config config_test.yaml --run $SLURM_ARRAY_TASK_ID --skip-topology-test
