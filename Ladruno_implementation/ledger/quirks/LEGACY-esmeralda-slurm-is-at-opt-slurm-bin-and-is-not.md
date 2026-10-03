---
wp: LEGACY
title: "esmeralda: SLURM is at /opt/slurm/bin and is NOT on PATH — this has now cost two different things"
legacy_seq: 225
---
### esmeralda: SLURM is at `/opt/slurm/bin` and is NOT on `PATH` — this has now cost two different things
- **Bites (1), the big one:** `which sbatch` / `which sinfo` return nothing, which reads as "there is no scheduler / the cluster is down". It is not. `esmeralda` had **33 days of uptime** while an entire ADR-75 session recorded the cluster as down and deferred both cluster-gated items on that basis. `/opt/slurm/bin/sinfo` answers instantly: 18 nodes, 32 cores / 60 GB each.
- **Bites (2):** a `#SBATCH` script that calls bare `srun` dies with **`rc=127` (command not found)** and an **empty log**. Two consecutive sweep submissions (jobs 144449, 144451) produced zero results this way. It is easy to misread as an MPI/launcher problem, because the *previous* failure on the same script genuinely was one (see the next row).
- **Workaround/status:** call `/opt/slurm/bin/{sbatch,srun,squeue,sinfo,sacct}` by absolute path from scripts, or prepend it to `PATH` at the top of every batch script. *2026-07-26 (ADR-75 P2h).*
