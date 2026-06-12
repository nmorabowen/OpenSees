"""Rank-env probe gate -- checker.

Opens each file rank_env_model.py produced and asserts the partition
manifest (INFO attrs PARTITION_ID / PARTITIONED / NUM_PARTITIONS) matches
the faked launcher environment:

  * plain      -> verbatim filename, PARTITIONED=0, NUM_PARTITIONS=1
  * pmi        -> .part-2, NUM_PARTITIONS=4   (Intel/MS-MPI names)
  * ompi       -> .part-5, NUM_PARTITIONS=8   (OpenMPI names)
  * slurm      -> .part-7, NUM_PARTITIONS=16  (srun names, pmix included)
  * precedence -> .part-2, NUM_PARTITIONS=4   (PMI pair outranks ambient SLURM)

Run with the venv python (has h5py):
    python rank_env_check.py [out_dir]
"""

from __future__ import annotations

import os
import sys

import h5py
import numpy as np

OUT = sys.argv[1] if len(sys.argv) > 1 else "."


def attr_int(group, name) -> int:
    """Manifest attrs are written as 1-element arrays; read either shape."""
    return int(np.asarray(group.attrs[name]).flat[0])

# (filename, expected PARTITIONED, PARTITION_ID, NUM_PARTITIONS)
EXPECT = [
    ("rank_env_plain.ladruno", 0, 0, 1),
    ("rank_env_pmi.part-2.ladruno", 1, 2, 4),
    ("rank_env_ompi.part-5.ladruno", 1, 5, 8),
    ("rank_env_slurm.part-7.ladruno", 1, 7, 16),
    ("rank_env_precedence.part-2.ladruno", 1, 2, 4),
]

fails: list[str] = []


def check(cond: bool, msg: str) -> None:
    print(("  ok  " if cond else " FAIL ") + msg)
    if not cond:
        fails.append(msg)


for fname, partitioned, part_id, num_parts in EXPECT:
    path = os.path.join(OUT, fname)
    if not os.path.exists(path):
        check(False, f"{fname}: file missing")
        continue
    with h5py.File(path, "r") as f:
        info = f["INFO"]
        check(attr_int(info, "PARTITIONED") == partitioned,
              f"{fname}: PARTITIONED == {partitioned}")
        check(attr_int(info, "PARTITION_ID") == part_id,
              f"{fname}: PARTITION_ID == {part_id}")
        check(attr_int(info, "NUM_PARTITIONS") == num_parts,
              f"{fname}: NUM_PARTITIONS == {num_parts}")

if fails:
    print(f"\nRANK_ENV CHECK: {len(fails)} failures")
    sys.exit(1)
print("\nRANK_ENV CHECK: all assertions passed")
