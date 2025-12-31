#!/bin/bash
set -e

# Script location
S_DIR="$(cd "$(dirname "$0")" && pwd)"
# Project root
ROOT_DIR="$(cd "$S_DIR/../.." && pwd)"

# Base directory for all test runs
BASE_DIR="$S_DIR/outputs"
mkdir -p "$BASE_DIR"

echo "--------------------------------------------------"
echo "Running Iceberg - BasicSpec (Serial)"
echo "--------------------------------------------------"
mkdir -p "$BASE_DIR/basic"
cd "$BASE_DIR/basic"
julia --project="$ROOT_DIR" "$S_DIR/Iceberg.jl"
mv outputs/* .
rmdir outputs

echo "--------------------------------------------------"
echo "Running Iceberg - ThreadedSpec"
echo "--------------------------------------------------"
for t in 2; do
    echo "  Threads: $t"
    mkdir -p "$BASE_DIR/threaded_t$t"
    cd "$BASE_DIR/threaded_t$t"
    julia --project="$ROOT_DIR" -t "$t" "$S_DIR/Iceberg.jl"
    mv outputs/* .
    rmdir outputs
done

echo "--------------------------------------------------"
echo "Running Iceberg - MPISpec"
echo "--------------------------------------------------"
for n in 2; do
    echo "  MPI Processes: $n"
    mkdir -p "$BASE_DIR/mpi_n$n"
    cd "$BASE_DIR/mpi_n$n"
    mpiexecjl -n "$n" julia --project="$ROOT_DIR" "$S_DIR/Iceberg.jl"
    mv outputs/* .
    rmdir outputs
done

echo "--------------------------------------------------"
echo "All Iceberg tests completed. Outputs in $BASE_DIR"
echo "--------------------------------------------------"
