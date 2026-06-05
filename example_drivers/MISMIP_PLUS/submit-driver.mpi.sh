#!/bin/bash -l
#SBATCH --job-name=MISMIP-mpi
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --partition=rocky
#SBATCH --account=rocky
#SBATCH --mem=64G

#SBATCH --output=logs/%A_%a.%N.out
#SBATCH --error=logs/%A_%a.%N.err

module load mpi/mpich-x86_64

echo "START: $(date +%s)"
export JULIA_DEBUG=WAVI

RUN_NAME="${SLURM_JOB_NAME}"
echo "RUN_NAME = $RUN_NAME"

time mpiexecjl --project=../.. -n $SLURM_NTASKS julia MISMIP_PLUS.jl "$RUN_NAME"

echo "FINISH: $(date +%s)"
