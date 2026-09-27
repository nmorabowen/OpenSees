"""Partitioned-reaction gate (WP-126) -- checker (h5py only, no apeGmsh).

Pairs with mp_reaction_model.py (the WP-126 reproduction: a fixed support node SHARED by
two partitions). Asserts the recorder half of the stitching contract:

  1. each part's REACTION_FORCE carries PARTITION_REDUCTION=SUM, DISPLACEMENT carries NONE;
  2. stitching per the contract (SUM -> add the parts' rows for a shared node; NONE -> any
     copy) reproduces the SERIAL reaction exactly, at every step;
  3. the partitioned -envelope files contain NO REACTION_FORCE envelope (refused), and
     the serial -envelope file does (one partition: exact).

    mpiexec -n 1 <py3.12> -S mp_reaction_model.py mp <openseesmp_dir> <out>
    mpiexec -n 2 <py3.12> -S mp_reaction_model.py mp <openseesmp_dir> <out>
    python mp_reaction_check.py <out>            (any python with h5py + numpy)
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import h5py
import numpy as np

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
NODE = 1


def _attr(g, k):
    """The recorder writes attributes as rank-1 arrays (same writer as COMPONENTS)."""
    v = g.attrs.get(k)
    if v is None:
        return None
    v = np.atleast_1d(v).flat[0]
    return v.decode() if isinstance(v, bytes) else str(v)


def _stage(f):
    return [k for k in f if k.startswith("MODEL_STAGE")][0]


def rows(path, result):
    """{node_id: DATA[T x nComp]} and the result's PARTITION_REDUCTION."""
    with h5py.File(path, "r") as f:
        g = f[f"{_stage(f)}/RESULTS/ON_NODES/{result}"]
        ids = g["ID"][:, 0]
        data = g["DATA"][...]                     # [T x nIds x nComp]
        return {int(n): data[:, i, :] for i, n in enumerate(ids)}, _attr(g, "PARTITION_REDUCTION")


def stitch(parts, result):
    merged, kinds = {}, set()
    for p in parts:
        r, kind = rows(p, result)
        kinds.add(kind)
        for n, v in r.items():
            if n in merged and kind == "SUM":
                merged[n] = merged[n] + v
            elif n not in merged:
                merged[n] = v.copy()
    return merged, kinds


def envelope_names(path):
    with h5py.File(path, "r") as f:
        base = f"{_stage(f)}/RESULTS/ENVELOPES/ON_NODES"
        return set(f[base]) if base in f else set()


def main() -> int:
    out = Path(sys.argv[1])
    problems = 0

    def check(ok, msg):
        nonlocal problems
        print(f"  [{'OK' if ok else 'FAIL'}] {msg}")
        problems += 0 if ok else 1

    serial = out / "mp_np1_stream.ladruno"
    parts = sorted(out.glob("mp_np2_stream.part-*.ladruno"))
    check(serial.exists() and len(parts) == 2, f"serial file + 2 part files present ({len(parts)} parts)")
    if problems:
        return 1

    for p in parts:
        _, kr = rows(p, "REACTION_FORCE")
        _, kd = rows(p, "DISPLACEMENT")
        check(kr == "SUM" and kd == "NONE", f"{p.name}: REACTION_FORCE={kr}, DISPLACEMENT={kd}")

    ref, _ = rows(serial, "REACTION_FORCE")
    got, _ = stitch(parts, "REACTION_FORCE")
    ok = NODE in got and np.array_equal(got[NODE], ref[NODE])
    check(ok, f"stitched reaction at shared node {NODE} == serial: "
              f"{got.get(NODE, np.array([]))[-1].tolist()} vs {ref[NODE][-1].tolist()}")
    first_wins = rows(parts[0], "REACTION_FORCE")[0][NODE][-1].tolist()
    print(f"       (first-partition-wins would give {first_wins} -- the pre-WP-126 apeGmsh answer)")

    dref, _ = rows(serial, "DISPLACEMENT")
    dgot, _ = stitch(parts, "DISPLACEMENT")
    check(all(np.array_equal(dgot[n], dref[n]) for n in dref if n in dgot),
          "stitched displacement (NONE, first copy) == serial on every shared node")

    for p in sorted(out.glob("mp_np2_env.part-*.ladruno")):
        names = envelope_names(p)
        check("REACTION_FORCE" not in names, f"{p.name}: reaction envelope refused ({sorted(names)})")
    names = envelope_names(out / "mp_np1_env.ladruno")
    check("REACTION_FORCE" in names, f"mp_np1_env.ladruno: serial reaction envelope kept ({sorted(names)})")

    print("\nMP_REACTION_CHECK:", "ALL PASS" if problems == 0 else f"{problems} PROBLEM(S)")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
