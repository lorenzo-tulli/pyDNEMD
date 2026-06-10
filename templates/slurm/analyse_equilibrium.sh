#!/bin/bash
#SBATCH --job-name=dnemd_analyse_eq
#SBATCH --output=logs/analyse_eq.out
#SBATCH --error=logs/analyse_eq.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --account=chem021482
#SBATCH --time=10:00:00
#SBATCH --partition=compute

module load openmpi/5.0.3
module load gromacs/2024.2-netlib-lapack

mkdir -p logs

dnemd-analyse-equilibrium --config config_test.yaml
