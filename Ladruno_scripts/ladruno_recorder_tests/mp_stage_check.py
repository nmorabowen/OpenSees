"""WP-163/165 multi-rank gate -- checker for mp_stage_model.tcl (+ m4_sequential_model.tcl).

    python3 mp_stage_check.py <out_dir> [np]

Asserts, for an np-rank (default 2) srun of mp_stage_model.tcl:
  * stage.part-<r>.ladruno for every rank, PARTITIONED=1, NUM_PARTITIONS=np,
    PARTITION_ID=r; one shared RUN_ID with RUN_ID_SCOPE="launcher"      (M4, MP-9)
  * exactly ONE MODEL_STAGE per file despite the mid-run SP pattern; DISPLACEMENT
    and energyBalance hold all 4 recorded steps                         (R6, R3)
  * region.part-1.ladruno: the stage is EMPTY_PARTITION=1 with zero-length
    NODES/ID; region.part-0.ladruno holds rank 0's two region nodes      (MP-8)
  * mon.part-<r>.h5 for every rank and no shared mon.h5                 (MP-6)
and, if present, m4_seq.ladruno (sequential run inside sbatch): PARTITIONED=0 (M4).
"""
from __future__ import annotations

import os
import sys

import h5py
import numpy as np

OUT = sys.argv[1]
NP = int(sys.argv[2]) if len(sys.argv) > 2 else 2
fails: list[str] = []


def check(cond, msg):
    print(("  ok  " if cond else " FAIL ") + msg)
    if not cond:
        fails.append(msg)


def a(group, name):
    v = group.attrs[name]
    v = v.flat[0] if hasattr(v, "flat") else v
    return v.decode() if isinstance(v, bytes) else v


def stages(f):
    return sorted(k for k in f if k.startswith("MODEL_STAGE"))


run_ids = set()
for r in range(NP):
    p = os.path.join(OUT, f"stage.part-{r}.ladruno")
    if not os.path.exists(p):
        check(False, f"{p} exists")
        continue
    with h5py.File(p, "r") as f:
        info = f["INFO"]
        check(int(a(info, "PARTITIONED")) == 1, f"stage.part-{r}: PARTITIONED == 1")
        check(int(a(info, "NUM_PARTITIONS")) == NP, f"stage.part-{r}: NUM_PARTITIONS == {NP}")
        check(int(a(info, "PARTITION_ID")) == r, f"stage.part-{r}: PARTITION_ID == {r}")
        check(str(a(info, "RUN_ID_SCOPE")) == "launcher",
              f"stage.part-{r}: RUN_ID_SCOPE == launcher (got {a(info, 'RUN_ID_SCOPE')})")
        run_ids.add(str(a(info, "RUN_ID")))
        st = stages(f)
        check(len(st) == 1, f"stage.part-{r}: one MODEL_STAGE despite the SP pattern (got {st})")
        if st:
            res = f[f"{st[0]}/RESULTS"]
            check(res["ON_NODES/DISPLACEMENT/DATA"].shape[0] == 4,
                  f"stage.part-{r}: DISPLACEMENT has 4 rows")
            check(res["ON_DOMAIN/energyBalance/DATA"].shape[0] == 4,
                  f"stage.part-{r}: energyBalance has 4 rows")
check(len(run_ids) == 1, f"one RUN_ID across ranks (got {sorted(run_ids)})")

for r in range(NP):
    p = os.path.join(OUT, f"region.part-{r}.ladruno")
    if not os.path.exists(p):
        check(False, f"{p} exists")
        continue
    with h5py.File(p, "r") as f:
        st = stages(f)
        last = f[st[-1]] if st else None
        if last is None:
            check(False, f"region.part-{r}: has a MODEL_STAGE")
            continue
        n_ids = last["MODEL/NODES/ID"].shape[0]
        empty = int(np.asarray(last.attrs.get("EMPTY_PARTITION", 0)).flat[0])
        if r == 0:
            check(n_ids == 2 and empty == 0, f"region.part-0: 2 region nodes, not empty ({n_ids}, {empty})")
        else:
            check(n_ids == 0 and empty == 1,
                  f"region.part-{r}: EMPTY_PARTITION with zero-length NODES ({n_ids}, {empty})")

for r in range(NP):
    check(os.path.exists(os.path.join(OUT, f"mon.part-{r}.h5")), f"mon.part-{r}.h5 exists")
check(not os.path.exists(os.path.join(OUT, "mon.h5")), "no shared mon.h5")

seq = os.path.join(OUT, "m4_seq.ladruno")
if os.path.exists(seq) or os.path.exists(os.path.join(OUT, "m4_seq.part-0.ladruno")):
    check(os.path.exists(seq), "sequential run in sbatch kept its filename (m4_seq.ladruno)")
    if os.path.exists(seq):
        with h5py.File(seq, "r") as f:
            check(int(a(f["INFO"], "PARTITIONED")) == 0, "m4_seq: PARTITIONED == 0")

if fails:
    print(f"\nMP_STAGE CHECK: {len(fails)} failures")
    sys.exit(1)
print("\nMP_STAGE CHECK: all assertions passed")
