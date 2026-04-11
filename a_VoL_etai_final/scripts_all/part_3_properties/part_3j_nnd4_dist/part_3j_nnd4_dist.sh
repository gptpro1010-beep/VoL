#!/bin/bash
#SBATCH --account=ant_colonies
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --time=0:50:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G

module reset

g++ -O3 -std=c++17 -fopenmp part_3j_nnd4_dist.cpp -o part_3j_nnd4_dist_exec

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

./part_3j_nnd4_dist_exec
