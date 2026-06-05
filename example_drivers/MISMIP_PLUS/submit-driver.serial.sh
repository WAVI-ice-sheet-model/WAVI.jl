#!/bin/bash -l
#SBATCH --job-name=MISMIP-serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --partition=rocky
#SBATCH --account=rocky
#SBATCH --mem=64G

#SBATCH --output=logs/%A_%a.%N.out
#SBATCH --error=logs/%A_%a.%N.err

echo "START: `date +%s`"
export JULIA_DEBUG=WAVI

RUN_NAME="${SLURM_JOB_NAME}"
echo "RUN_NAME = $RUN_NAME"

julia --project=../.. MISMIP_PLUS.jl

echo "FINISH: `date +%s`"
