"""Parallel MPCO_Ladruno output gate -- checker (venv python, h5py).

Given the expected partition count N and the output dir, assert:
  1. exactly N files "mp_out.part-0.ladruno" .. "mp_out.part-(N-1).ladruno" exist
     (0-based, contiguous from 0 -- the apeGmsh stitch contract) and NO bare
     "mp_out.ladruno" (the pre-fix collision target);
  2. each part file passes the schema-v1 validator;
  3. each file's INFO says PARTITIONED=1, NUM_PARTITIONS=N, PARTITION_ID=its index;
  4. each file holds a nonempty MODEL/NODES (each rank wrote its own partition).

    python mp_parallel_check.py <out_dir> <N>
"""

from __future__ import annotations

import os
import re
import sys

import h5py
import numpy as np

import ladruno_format as lf


def _attr(g, k):
    v = g.attrs[k]
    a = np.atleast_1d(v)
    v = a[0] if a.size == 1 else a
    return v.decode() if isinstance(v, bytes) else v


def main() -> int:
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    problems = 0

    if os.path.exists(os.path.join(out, "mp_out.ladruno")):
        print("  [FAIL] bare mp_out.ladruno exists -> ranks collided on one file")
        problems += 1

    found = sorted(
        int(m.group(1))
        for f in os.listdir(out)
        for m in [re.match(r"^mp_out\.part-(\d+)\.ladruno$", f)]
        if m
    )
    print(f"partition files found: {found} (expected 0..{n - 1})")
    if found != list(range(n)):
        print(f"  [FAIL] partition set {found} != contiguous 0..{n - 1}")
        problems += 1

    for idx in found:
        path = os.path.join(out, f"mp_out.part-{idx}.ladruno")
        probs = lf.validate(path)
        if probs:
            print(f"  [FAIL] part-{idx} schema violations: {probs}")
            problems += 1
        with h5py.File(path, "r") as f:
            info = f["INFO"]
            part_id = int(_attr(info, "PARTITION_ID"))
            nparts = int(_attr(info, "NUM_PARTITIONS"))
            partitioned = int(_attr(info, "PARTITIONED"))
            stage = [k for k in f if k.startswith("MODEL_STAGE")][0]
            nnodes = f[f"{stage}/MODEL/NODES/ID"].shape[0]
            ok = (part_id == idx and nparts == n and partitioned == 1 and nnodes > 0)
            tag = "OK" if ok else "FAIL"
            if not ok:
                problems += 1
            print(f"  [{tag}] part-{idx}: PARTITION_ID={part_id} NUM_PARTITIONS={nparts} "
                  f"PARTITIONED={partitioned} nodes={nnodes}")

    print("\n" + "=" * 56)
    if problems == 0:
        print("MP_PARALLEL_CHECK: ALL PASS")
        return 0
    print(f"MP_PARALLEL_CHECK: {problems} PROBLEM(S)")
    return 1


if __name__ == "__main__":
    sys.exit(main())
