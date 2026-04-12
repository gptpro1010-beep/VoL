# a_VoL_etai_final

This project simulates particles that start from a lattice and evolve using Vicsek-type alignment dynamics with periodic boundaries. The workflow sweeps density (`rho`) and speed (`s`) at fixed noise (`eta`), then saves order-parameter time series and particle snapshots for analysis. Downstream scripts compute spatial/temporal statistics and generate post-processed summaries. See `file_tree.txt` for the project file structure and folder-level details.

## Simulation parameter values (part_1_sim)

- **Densities (`rho`)**: `1.00, 2.00, 4.00, 8.00, 16.00`
- **Speeds (`s`)**: `0.01, 0.02, 0.04, 0.08, 0.16`
- **Noise (`eta`)**: `0.01`
- **Number of particles (`N`)**: `100`
- **Initial lattice side (`n_side`)**: `10`
- **Number of realizations (`n_sims`)**: `20` (`sim_01` ... `sim_20`)
- **Interaction radius (`r`)**: `1.0`
- **Box size (`L`)**: `sqrt(N / rho)`
- **Characteristic time (`T`)**: `round(L / s)`
- **Max integration steps (`T_max`)**: `5000 * T`
- **Snapshot interval in steps (`dt_out`)**: `5 * T`
- **Order-parameter write interval**: every `T` steps
- **Scheduler walltime request**: `#SBATCH --time=0:10:00`

## Expected output file structure (simulation output)

Expected structure for one parameter combination:
`sim_data/vicsek_density_<rho>_speed_<s>_noise_0.01/`

Expected tree for **`sim_01` only**:

```text
sim_data/vicsek_density_<rho>_speed_<s>_noise_0.01/
└── sim_01/
    ├── order_param_noise_0.01_sim_01_timeseries.txt
    └── sim_data/
        ├── sim_data_noise_0.01_sim_01_snap_00000.txt
        ├── sim_data_noise_0.01_sim_01_snap_00001.txt
        ├── ...
        └── sim_data_noise_0.01_sim_01_snap_<last>.txt
```

- Output for `sim_02` to `sim_20` follows the same structure and naming pattern.
- Figure outputs are intentionally omitted here.
