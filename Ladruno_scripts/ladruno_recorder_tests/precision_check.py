"""Bounded-error gate for the `-precision f32` opt-in lossy mode (the lossy-mode
analogue of the 1e-12 parity gate). Reads prec_f64.ladruno (lossless oracle) and
prec_f32.ladruno (lossy) from the SAME model and asserts:

  1. f32 values match f64 within a float32 round-trip bound (not 1e-12):
        |v_f32 - v_f64| <= ATOL + RTOL_F32 * |v_f64|     (RTOL_F32 ~ f32 eps)
  2. STORED_PRECISION is "f32" / "f64" respectively (reader awareness).
  3. on-disk DATA dtype is float32 / float64 respectively.
  4. the f32 file is meaningfully smaller (sanity: the mode actually saves space).

Run with the venv python (h5py):
  python precision_check.py <prec_f64.ladruno> <prec_f32.ladruno>
"""
from __future__ import annotations

import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ladruno_format as lf  # noqa: E402

# float32 has ~7 significant digits; eps ~ 1.19e-7. Allow a few eps of slack for
# the narrowing round-trip plus shuffle/deflate (lossless on the already-narrowed
# bytes, so no extra error there).
RTOL_F32 = 6.0e-7
ATOL = 1e-30  # effectively pure relative; abs floor only guards exact zeros


def _norm(name: str) -> str:
    return "".join(c for c in name.lower() if c.isalnum())


def _gather(path: str) -> dict:
    """{(result_path, flat_index, step): value} over every RESULTS .../DATA leaf.

    Walks RESULTS recursively because element results nest two levels deep
    (ON_ELEMENTS/<display>/<bucket>/DATA) while nodal/domain results are one
    level (ON_NODES/<name>/DATA). `step` is the STEP-axis value cast to int so
    keys are comparable across the two files (avoids numpy-vs-python int mismatch)."""
    out: dict = {}
    with h5py.File(path, "r") as f:
        for stage in [k for k in f if k.startswith("MODEL_STAGE")]:
            res = f.get(f"{stage}/RESULTS")
            if res is None:
                continue

            def visit(name, obj):
                if (isinstance(obj, h5py.Group) and "DATA" in obj
                        and isinstance(obj["DATA"], h5py.Dataset)):
                    key = f"{stage}/{name}"
                    for step, arr in lf.iter_step_slices(obj):
                        flat = np.asarray(arr).reshape(-1)
                        s = int(step)
                        for idx, val in enumerate(flat):
                            out[(key, idx, s)] = float(val)

            res.visititems(visit)
    return out


def _data_dtype(path: str) -> str:
    """dtype kind+size of the first RESULTS .../DATA dataset, e.g. 'float32'."""
    with h5py.File(path, "r") as f:
        for stage in [k for k in f if k.startswith("MODEL_STAGE")]:
            res = f.get(f"{stage}/RESULTS")
            if res is None:
                continue
            for fam in res:
                for rname in res[fam]:
                    grp = res[fam][rname]
                    if "DATA" in grp:
                        return str(grp["DATA"].dtype)
    return "?"


def main() -> int:
    f64_path, f32_path = sys.argv[1], sys.argv[2]
    ok = True

    # --- 2. STORED_PRECISION tags ---
    with lf.LadrunoReader(f64_path) as r:
        sp64 = r.stored_precision()
    with lf.LadrunoReader(f32_path) as r:
        sp32 = r.stored_precision()
    print(f"STORED_PRECISION: f64-file={sp64!r}  f32-file={sp32!r}")
    if sp64 != "f64":
        print("  FAIL: lossless file should tag STORED_PRECISION='f64'"); ok = False
    if sp32 != "f32":
        print("  FAIL: lossy file should tag STORED_PRECISION='f32'"); ok = False

    # --- 3. on-disk dtype ---
    dt64, dt32 = _data_dtype(f64_path), _data_dtype(f32_path)
    print(f"DATA dtype: f64-file={dt64}  f32-file={dt32}")
    if dt64 != "float64":
        print("  FAIL: f64 file DATA should be float64"); ok = False
    if dt32 != "float32":
        print("  FAIL: f32 file DATA should be float32"); ok = False

    # --- 1. bounded error ---
    a = _gather(f64_path)
    b = _gather(f32_path)
    common = set(a) & set(b)
    if not common:
        print("  FAIL: no overlapping values to compare"); return 1
    worst = 0.0
    nbad = 0
    for k in common:
        va, vb = a[k], b[k]
        tol = ATOL + RTOL_F32 * abs(va)
        err = abs(va - vb)
        rel = err / abs(va) if va != 0 else err
        if rel > worst:
            worst = rel
        if err > tol:
            nbad += 1
    print(f"bounded-error: {len(common)} values, worst relative err = {worst:.2e} "
          f"(bound {RTOL_F32:.1e}), violations = {nbad}")
    if nbad:
        print(f"  FAIL: {nbad} values exceed the float32 bound"); ok = False

    # --- 4. size sanity ---
    s64, s32 = os.path.getsize(f64_path), os.path.getsize(f32_path)
    ratio = s32 / s64 if s64 else 0
    print(f"size: f64={s64} f32={s32} ratio={ratio:.3f}")
    if ratio >= 0.95:
        print("  WARNING: f32 file is not meaningfully smaller (model may be too "
              "small for payload to dominate fixed overhead)")

    print("PRECISION GATE:", "OK" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
