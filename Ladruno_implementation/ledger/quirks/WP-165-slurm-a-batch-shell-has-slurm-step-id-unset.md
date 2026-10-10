---
wp: WP-165
title: "Slurm: a batch shell has SLURM_STEP_ID UNSET; each srun rank has a real step id (verified on Esmeralda, WP-165)"
date: 2026-10-03
---
### Slurm: a batch shell has `SLURM_STEP_ID` UNSET; each `srun` rank has a real step id (verified on Esmeralda, WP-165)
- **Observed (Slurm on Esmeralda, 2026-10-03):** inside `sbatch -n 2` the batch shell has `SLURM_NTASKS=2`,
  `SLURM_PROCID=0`, `SLURM_STEP_ID` unset; under `srun --mpi=pmix_v3` each rank has `SLURM_PROCID=<rank>`,
  `SLURM_STEP_ID=0` (first step). This is what `Ladruno_LaunchEnv.h` relies on (SLURM pair trusted only with a
  real step id; the batch/extern pseudo-steps use 0xFFFFFFFE/0xFFFFFFFD where a Slurm version does set them).
- **Also:** Esmeralda compute nodes have no `cmake`; OpenSees builds there run on the login host (niced), and
  the serial `OpenSees` binary needs `LD_LIBRARY_PATH=/mnt/nfshare/lib` on compute nodes (MKL lives in `/lib`
  only on the login host).
