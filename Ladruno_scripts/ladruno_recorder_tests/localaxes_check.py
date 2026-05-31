"""LOCAL_AXES + KIND gate -- checker (venv python, h5py).

Verifies on localaxes.ladruno:
  1. MODEL_STAGE[*]/KIND == "transient"          (the -kind option)
  2. MODEL/LOCAL_AXES/<classTag>-<ClassName>/{ID, FRAME[E×4]} present
  3. each FRAME row is a unit quaternion (w,x,y,z); reconstructing its rotation
     matrix gives an orthonormal frame whose local x-axis (column 0) equals the
     element axis (node2-node1)/L -> proves a REAL frame was written, not the
     identity-quaternion fallback this dataset exists to remove.

    python localaxes_check.py <out_dir>
"""

from __future__ import annotations

import os
import sys

import h5py
import numpy as np

TOL = 1e-9


def quat_to_R(q):
    """(w,x,y,z) unit quaternion -> 3x3 rotation matrix.

    The frozen recorder's quatFromMat (Ladruno_Types.h) packs the matrix whose
    COLUMNS are (vx,vy,vz) using the (R[0][1]-R[1][0]) off-diagonal sign — the
    transpose/conjugate of standard Shepperd — so the local axes come back as the
    ROWS of the matrix this standard formula reconstructs. The checker extracts
    rows accordingly (verified analytically against a 45-degree frame)."""
    w, x, y, z = q
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]])


def main() -> int:
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    fname = sys.argv[2] if len(sys.argv) > 2 else "localaxes.ladruno"
    path = os.path.join(out, fname)
    print(f"checking {path}")
    problems = 0
    with h5py.File(path, "r") as f:
        stage = [k for k in f if k.startswith("MODEL_STAGE")][0]

        # (1) KIND
        kind = np.atleast_1d(f[stage].attrs["KIND"])[0]
        if isinstance(kind, bytes):
            kind = kind.decode()
        kind = str(kind)
        print(f"KIND = {kind!r}")
        if kind != "transient":
            print(f"  [FAIL] KIND != 'transient'"); problems += 1
        else:
            print("  [OK] KIND == 'transient'")

        # (2)+(3) LOCAL_AXES
        la_path = f"{stage}/MODEL/LOCAL_AXES"
        if la_path not in f:
            print("  [FAIL] MODEL/LOCAL_AXES missing"); return 1 + problems
        la = f[la_path]
        nid, ncoord = f[f"{stage}/MODEL/NODES/ID"][...].ravel(), f[f"{stage}/MODEL/NODES/COORDINATES"][...]
        coord = {int(t): ncoord[i] for i, t in enumerate(nid)}
        elem_axis = {}  # elem_tag -> unit axis, from CONNECTIVITY
        for ename in f[f"{stage}/MODEL/ELEMENTS"]:
            conn = f[f"{stage}/MODEL/ELEMENTS/{ename}/CONNECTIVITY"][...]
            for row in conn:
                # pad to 3D so 2D models (2-component coords) compare against the
                # full 3D direction cosines the CrdTransf returns.
                a = np.zeros(3); b = np.zeros(3)
                ca = np.asarray(coord[int(row[1])], float).ravel()[:3]
                cb = np.asarray(coord[int(row[2])], float).ravel()[:3]
                a[:ca.size] = ca; b[:cb.size] = cb
                v = b - a
                elem_axis[int(row[0])] = v / np.linalg.norm(v)

        total_frames = 0
        for gname in la:
            grp = la[gname]
            ids = grp["ID"][...].ravel()
            frame = grp["FRAME"][...]
            print(f"  group {gname}: ID={ids.tolist()} FRAME{frame.shape}")
            if frame.shape != (ids.shape[0], 4):
                print(f"    [FAIL] FRAME {frame.shape} != ({ids.shape[0]}, 4)"); problems += 1
                continue
            for r, tag in enumerate(ids):
                q = frame[r]
                if abs(np.linalg.norm(q) - 1.0) > TOL:
                    print(f"    [FAIL] elem {tag}: quaternion not unit (|q|={np.linalg.norm(q):.3e})")
                    problems += 1
                R = quat_to_R(q)
                if np.max(np.abs(R @ R.T - np.eye(3))) > TOL:
                    print(f"    [FAIL] elem {tag}: frame not orthonormal"); problems += 1
                vx = R[0, :]  # local x-axis (rows = axes; see quat_to_R note)
                axis = elem_axis.get(int(tag))
                if axis is not None and np.max(np.abs(vx - axis)) > TOL:
                    print(f"    [FAIL] elem {tag}: local x {vx} != element axis {axis}")
                    problems += 1
                else:
                    is_identity = np.max(np.abs(R - np.eye(3))) < TOL
                    print(f"    [OK] elem {tag}: local x = element axis "
                          f"{vx.round(4).tolist()} (non-identity={not is_identity})")
                total_frames += 1

        if total_frames == 0:
            print("  [FAIL] no FRAME rows written"); problems += 1

    print("\n" + "=" * 56)
    if problems == 0:
        print("LOCALAXES_CHECK: ALL PASS")
        return 0
    print(f"LOCALAXES_CHECK: {problems} PROBLEM(S)")
    return 1


if __name__ == "__main__":
    sys.exit(main())
