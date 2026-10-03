---
wp: LEGACY
title: "mpirun works inside a 1-node SLURM allocation and dies instantly across 2 nodes"
legacy_seq: 226
---
### `mpirun` works inside a 1-node SLURM allocation and dies instantly across 2 nodes
- **Bites:** the identical sweep script runs fine on `--nodes=1` and, on `--nodes=2`, every launch fails immediately with `[[...]] FORCE-TERMINATE AT (null):1 - error plm_slurm_module.c(471)` / `An internal error has occurred in ORTE`. OpenMPI's ORTE SLURM launcher, not the model.
- **Workaround/status:** use `/opt/slurm/bin/srun --cpu-bind=cores --mpi=pmix_v3 <wrapper> deck.tcl` for anything multi-node — the launcher `02_esmeralda_linux_build_guide.md` §7 already documents. `srun --mpi=list` confirms `pmix_v3` is available. *2026-07-26 (ADR-75 P2h).*
