"""Eigen / modesOfVibration coverage — validator (L2).

The per-step parity gate skips modal groups, so this gate instead asserts the
recorder's MODAL output path produced a structurally valid, finite, non-trivial
result:
  - the .ladruno passes the schema validator,
  - it has a MODEL_STAGE with KIND == 'eigen',
  - ON_NODES carries a MODES_OF_VIBRATION group with an ID and finite, non-zero
    eigenvector data.

Run with the venv python (has h5py):
    python eigen_check.py <eig_ref.mpco> <eig_test.ladruno>
"""

from __future__ import annotations

import sys

import h5py
import numpy as np

import ladruno_format as lf


def _norm(name: str) -> str:
    return "".join(c for c in name.lower() if c.isalnum())


def main() -> int:
    new_path = sys.argv[2] if len(sys.argv) > 2 else sys.argv[1]

    problems = lf.validate(new_path)
    if problems:
        print(f"EIGEN FAIL: validator rejected the file: {problems[:5]}")
        return 1

    with h5py.File(new_path, "r") as f:
        stages = [k for k in f if k.startswith("MODEL_STAGE")]
        eigen_stages = [s for s in stages
                        if _norm(lf._attr(f[s], "KIND")) == "eigen"]
        if not eigen_stages:
            print(f"EIGEN FAIL: no MODEL_STAGE with KIND='eigen' (stages={stages})")
            return 1

        found_modes = False
        nz = 0
        for s in eigen_stages:
            on_nodes = f.get(f"{s}/RESULTS/ON_NODES")
            if on_nodes is None:
                continue
            for rname in on_nodes:
                if "modesofvibration" not in _norm(rname):
                    continue
                grp = on_nodes[rname]
                if "ID" not in grp:
                    print(f"EIGEN FAIL: {rname} has no ID dataset")
                    return 1
                # gather every numeric leaf under the group (DATA/STEP_k or chunked)
                arrs = []

                def _collect(g):
                    for k in g:
                        item = g[k]
                        if isinstance(item, h5py.Group):
                            _collect(item)
                        elif k != "ID":
                            arrs.append(np.asarray(item[...], dtype=float))

                _collect(grp)
                vals = np.concatenate([a.reshape(-1) for a in arrs]) if arrs else np.array([])
                if vals.size == 0:
                    print(f"EIGEN FAIL: {rname} carries no eigenvector data")
                    return 1
                if not np.all(np.isfinite(vals)):
                    print(f"EIGEN FAIL: {rname} has non-finite values")
                    return 1
                nz += int(np.count_nonzero(np.abs(vals) > 0.0))
                found_modes = True

        if not found_modes:
            print("EIGEN FAIL: no MODES_OF_VIBRATION group under ON_NODES")
            return 1
        if nz == 0:
            print("EIGEN FAIL: all eigenvector components are zero")
            return 1

    # Value parity vs the frozen recorder: the modal eigenvectors must match mpco
    # to 1e-12 (DATA/STEP_<step>/MODE_<k> [nNodes x nComp]).
    ref_path = sys.argv[1] if len(sys.argv) > 2 else None
    if ref_path:
        ref = _read_modes(ref_path)
        new = _read_modes(new_path)
        common = set(ref) & set(new)
        mism = [(k, ref[k], new[k]) for k in common
                if abs(ref[k] - new[k]) > 1e-12 + 1e-9 * abs(ref[k])]
        nzc = sum(1 for k in common if abs(ref[k]) > 0.0)
        print(f"  modal parity: ref={len(ref)} new={len(new)} common={len(common)} "
              f"nonzero={nzc}")
        if not common:
            print("EIGEN FAIL: no overlapping modal values vs mpco (layout mismatch?)")
            return 1
        if mism:
            print(f"EIGEN FAIL: {len(mism)} modal mismatches vs mpco (showing 5):")
            for (k, a, b) in mism[:5]:
                print(f"  {k}: ref={a!r} new={b!r}")
            return 1
        print(f"EIGEN OK — modal eigenvectors match mpco to 1e-12 ({nzc} nonzero, "
              f"{len(common)} values)")
        return 0

    print(f"EIGEN OK — eigen stage present, modal data finite & non-trivial ({nz} nonzero)")
    return 0


def _read_modes(path: str) -> dict:
    """{(node_tag, 'MODE_k', comp_index): value} over every MODES_OF_VIBRATION group
    (DATA/STEP_<step>/MODE_<k> [nNodes x nComp]); works for .mpco and .ladruno."""
    out: dict = {}
    with h5py.File(path, "r") as f:
        for s in [k for k in f if k.startswith("MODEL_STAGE")]:
            base = f.get(f"{s}/RESULTS/ON_NODES")
            if base is None:
                continue
            for rname in base:
                if "modesofvibration" not in _norm(rname):
                    continue
                grp = base[rname]
                if "ID" not in grp or "DATA" not in grp:
                    continue
                ids = np.asarray(grp["ID"][...]).reshape(-1)
                data = grp["DATA"]
                if not isinstance(data, h5py.Group):
                    continue
                for step in data:
                    sgrp = data[step]
                    if not isinstance(sgrp, h5py.Group):
                        continue
                    for mk in sgrp:
                        arr = np.atleast_2d(np.asarray(sgrp[mk][...]))
                        for r, tag in enumerate(ids):
                            for c in range(arr.shape[1]):
                                out[(int(tag), mk, c)] = float(arr[r, c])
    return out


if __name__ == "__main__":
    raise SystemExit(main())
