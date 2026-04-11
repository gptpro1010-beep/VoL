# a_VoL_etai_final

This folder contains scripts used in the VoL workflow.

## What is here

- `scripts_all/part_1_sim/`: simulation template files and run commands.
- `scripts_all/part_2_anim/`: MATLAB scripts for movies and snapshots.
- `scripts_all/part_3_properties/`: C++ and shell scripts for property calculations.
- `scripts_all/part_4_post/`: MATLAB scripts for post-processing.
- `params.txt`: notes about key parameters.
- `file_tree.txt`: the detailed file and folder list.

## Simulation parameter values (part_1_sim)

These values are defined in `scripts_all/part_1_sim/commands.txt` and `script_1_Vicsek_sim_template.cpp`.

- **Densities (`rho`)**: `1.00, 2.00, 4.00, 8.00, 16.00`
- **Speeds (`v`)**: `0.01, 0.02, 0.04, 0.08, 0.16`
- **Noise (`eta`)**: `0.01`
- **Number of particles (`N`)**: `100`
- **Initial lattice side (`n_side`)**: `10`
- **Number of realizations (`n_sims`)**: `20` (`sim_01` ... `sim_20`)
- **Interaction radius (`r`)**: `1.0`
- **Box size (`L`)**: `sqrt(N / rho)`
- **Characteristic time (`T`)**: `round(L / v)`
- **Max integration steps (`T_max`)**: `5000 * T`
- **Snapshot interval in steps (`dt_out`)**: `5 * T`
- **Order-parameter write interval in steps**: every `T` steps
- **Scheduler walltime request**: `#SBATCH --time=0:10:00`

## Expected output file structure (simulation output)

Expected structure for one parameter combination:
`sim_data/vicsek_density_<rho>_speed_<speed>_noise_0.01/`

Below is the expected tree for **`sim_01` only** (as requested):

```text
sim_data/vicsek_density_<rho>_speed_<speed>_noise_0.01/
├── script_1_Vicsek_sim.cpp
├── run_vicsek.sh
├── vicsek_exec
└── sim_01/
    ├── order_param_noise_0.01_sim_01_timeseries.txt
    └── sim_data/
        ├── sim_data_noise_0.01_sim_01_snap_00000.txt
        ├── sim_data_noise_0.01_sim_01_snap_00001.txt
        ├── ...
        └── sim_data_noise_0.01_sim_01_snap_<last>.txt
```

Notes:
- Output for `sim_02` to `sim_20` follows the same structure and naming pattern.
- Figure outputs are intentionally omitted here.

## Note

- `part_5_verification/` is listed in `file_tree.txt` for consistency, but it is not in the current `scripts_all/` folder.
