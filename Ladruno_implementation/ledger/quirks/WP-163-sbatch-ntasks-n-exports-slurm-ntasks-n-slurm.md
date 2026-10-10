---
wp: WP-163
title: "sbatch --ntasks=N exports SLURM_NTASKS=N / SLURM_PROCID=0 into the batch shell itself — trust the SLURM pair only inside an srun step (WP-163)"
date: 2026-10-03
---
### `sbatch --ntasks=N` exports `SLURM_NTASKS=N` / `SLURM_PROCID=0` into the batch shell itself — trust the SLURM pair only inside an srun step (WP-163)
- **Bites:** the recorders pick their partition file from the launcher env. A plain SEQUENTIAL run inside a batch
  script (no `srun`) saw `SLURM_NTASKS=4` and wrote `<stem>.part-0.ladruno` with `NUM_PARTITIONS=4`; `-envelope`
  of reactions was then refused as "partitioned". A launcher outside the probe table inheriting the batch vars gave
  EVERY rank `SLURM_PROCID=0` → all ranks truncated the same `part-0`.
- **Workaround/status:** `SRC/recorder/Ladruno_LaunchEnv.h` (WP-163 M4/MP-3): the SLURM pair counts only when
  `SLURM_STEP_ID` is a real step number (srun sets it; the batch/extern pseudo-steps use 0xFFFFFFFE/0xFFFFFFFD);
  SIZE > 1 with a missing/out-of-range RANK is an error. One probe for the ladruno, EnergyBalance and Monitor
  recorders — keep it the only copy.
