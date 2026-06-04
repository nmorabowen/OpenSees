"""Shell per-LAYER material-stress emission — check (L1/L2).

Asserts the fix that re-enabled the per-layer read through the natural
``material.fiber.<resp>`` shell verb:

  1. the material-verb file contains a real ON_ELEMENTS per-layer bucket
     (COLUMN_MAP present; ncols == nGP * nLayers * nStressComp), NOT empty;
  2. it carries nonzero stress data;
  3. it is byte-identical (per (elem, column, step)) to the section-verb file,
     i.e. material.fiber.stress == section.fiber.stress.

Pure h5py+numpy (runs under the venv check python).
    python shell_layer_check.py <material.ladruno> <section.ladruno>
"""
from __future__ import annotations

import os
import sys

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import ladruno_format as lf  # noqa: E402

EXPECT_NGP = 4        # ShellMITC4
EXPECT_NLAYERS = 3
EXPECT_NCOMP = 5      # plate-fiber stress components
EXPECT_NCOLS = EXPECT_NGP * EXPECT_NLAYERS * EXPECT_NCOMP   # 60
EXPECT_NELEM = 4


def _strip_rname(nm):
    out = {}
    for (t, gp, sec, fib, rc, step), v in nm.items():
        comp = rc.split(":", 1)[1]
        out[(t, gp, sec, fib, comp, step)] = v
    return out


def _bucket_cols(reader, stage):
    base_path = f"{stage}/RESULTS/ON_ELEMENTS"
    if base_path not in reader.f:
        return 0, 0
    base = reader.f[base_path]
    ncols = nids = 0
    for rname in base:
        for bname in base[rname]:
            b = base[rname][bname]
            if "COLUMN_MAP" not in b:
                continue
            ncols += len(lf.decode_columns(b["COLUMN_MAP"]))
            nids = max(nids, int(b["ID"].shape[0]))
    return ncols, nids


def main(mat_path: str, sec_path: str) -> int:
    fails = []
    for p in (mat_path, sec_path):
        if not os.path.exists(p):
            print(f"[FAIL] missing file: {p}")
            return 1

    with lf.LadrunoReader(mat_path) as rm, lf.LadrunoReader(sec_path) as rs:
        stm, sts = rm.stages()[0], rs.stages()[0]

        m_cols, m_ids = _bucket_cols(rm, stm)
        s_cols, s_ids = _bucket_cols(rs, sts)

        # (1) material verb must now emit a real per-layer bucket
        if m_cols == 0:
            fails.append("material.fiber.stress emitted NO per-layer bucket (the bug)")
        if m_cols != EXPECT_NCOLS:
            fails.append(f"material ncols {m_cols} != expected {EXPECT_NCOLS} "
                         f"(nGP*nLayers*nComp = {EXPECT_NGP}*{EXPECT_NLAYERS}*{EXPECT_NCOMP})")
        if m_ids != EXPECT_NELEM:
            fails.append(f"material nIDs {m_ids} != {EXPECT_NELEM}")
        # section verb is the reference that already worked
        if s_cols != EXPECT_NCOLS:
            fails.append(f"section ncols {s_cols} != expected {EXPECT_NCOLS}")

        nm = _strip_rname(lf.normalize_element(rm, stm))
        ns = _strip_rname(lf.normalize_element(rs, sts))

        # (2) nonzero data
        nonzero = sum(1 for v in nm.values() if abs(v) > 1e-12)
        if nonzero == 0:
            fails.append("material per-layer data is all zero")

        # (3) byte-identical to the section verb
        shared = set(nm) & set(ns)
        if not shared:
            fails.append("material/section share no comparable (elem,col,step) keys")
        else:
            maxd = max(abs(nm[k] - ns[k]) for k in shared)
            if len(nm) != len(ns) or len(shared) != len(nm):
                fails.append(f"key mismatch: |mat|={len(nm)} |sec|={len(ns)} shared={len(shared)}")
            if maxd > 1e-12:
                fails.append(f"material vs section maxdiff {maxd:.3e} > 1e-12")

    if fails:
        for f in fails:
            print(f"[FAIL] SHELL LAYER STRESS: {f}")
        return 1
    print(f"SHELL_LAYER_CHECK: ALL PASS (ncols={EXPECT_NCOLS}, nelem={EXPECT_NELEM}, "
          f"material.fiber.stress == section.fiber.stress)")
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("usage: shell_layer_check.py <material.ladruno> <section.ladruno>")
        sys.exit(2)
    sys.exit(main(sys.argv[1], sys.argv[2]))
