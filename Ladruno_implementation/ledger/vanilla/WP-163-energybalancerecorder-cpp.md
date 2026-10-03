---
wp: WP-163
title: "WP-163 vanilla rows"
pr: "#917"
date: 2026-10-03
files: ["`SRC/recorder/EnergyBalanceRecorder.cpp`"]
table: "main"
---
| `SRC/recorder/EnergyBalanceRecorder.cpp` | `// Ladruno WP-163`: (M1) skip SP ShadowSubdomains in the max-DOF / hourglass-probe sizing loop (the shared kernel skips them in the energy sweep — a Subdomain has `getNodePtrs()==0`); (M4/MP-3) the per-rank `.part-<rank>` filename now comes from the shared `Ladruno_LaunchEnv.h` probe — SLURM only inside an srun step, SIZE > 1 without a valid RANK refused (was: an inline copy that wrote `part-0` for a sequential run inside `sbatch`). | [#917](https://github.com/nmorabowen/OpenSees/pull/917) |
