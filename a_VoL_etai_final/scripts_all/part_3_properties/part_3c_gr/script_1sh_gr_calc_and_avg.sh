#!/bin/bash
#SBATCH --account=ant_colonies
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=8G

module reset

g++ -O3 -std=c++17 -fopenmp script_1_gr_calc_and_avg.cpp -o script_1_gr_calc_and_avg_exec

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

./script_1_gr_calc_and_avg_exec
