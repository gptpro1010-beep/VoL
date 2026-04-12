#!/bin/bash
#SBATCH --account=ant_colonies
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --time=0:40:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G

module reset

g++ -O3 -std=c++17 script_1_lcc_indices_calc.cpp -o script_1_lcc_indices_calc_exec

./script_1_lcc_indices_calc_exec
