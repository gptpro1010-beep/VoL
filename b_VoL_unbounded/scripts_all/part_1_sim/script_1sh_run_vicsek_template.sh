#!/bin/bash
#SBATCH --account=ant_colonies
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --time=0:10:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=4G
#SBATCH --job-name=vicsek_unbounded_rho_{RHO}_speed_{SPEED}

module reset

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

noise=0.01
rho={RHO}
speed={SPEED}

echo "----------------------------------"
echo "Running unbounded Vicsek simulation"
echo "rho   = $rho"
echo "speed = $speed"
echo "noise = $noise"
echo "----------------------------------"

if [ ! -f vicsek_exec ]; then
    echo "Compiling simulation code..."
    g++ -O3 -std=c++17 -fopenmp script_1_Vicsek_sim.cpp -o vicsek_exec
fi

echo "Starting simulation..."
./vicsek_exec

echo "Simulation finished."
