"""Element-result ENVELOPE gate -- checker (venv python, h5py).

For each block (quad stress / fiber-beam section.fiber.stress):
  1. schema validate the envelope file (checks ENVELOPES shapes,
     ABSMAX==max(|MIN|,|MAX|), and the element COLUMN_MAP consistency);
  2. reduce the time series to componentwise MIN/MAX/ABSMAX/ARG_STEP over steps
     and match the ENVELOPES datasets (1e-9);
  3. verify the envelope COLUMN_MAP round-trips the time-series COLUMN_MAP -- the
     per-column (gauss_id, section_tag, fiber_id, comp) structure must be identical
     (the whole point of item (1): gauss/section/fiber structure on envelopes).

    python element_envelope_check.py <out_dir>
"""

from __future__ import annotations

import os
import sys

import numpy as np

import ladruno_format as lf

TOL = 1e-9

BLOCKS = [("quad stress", "ee_quad"), ("fiber-beam section.fiber.stress", "ee_beam")]


def reduce_elements(path):
    """name ("<display>/<bucket>") -> (ids, MIN, MAX, ABSMAX, ARG_STEP, columns)
    reduced over all ON_ELEMENTS time-series steps."""
    out = {}
    with lf.LadrunoReader(path) as r:
        stage = r.stages()[0]
        base = r.f[f"{stage}/RESULTS/ON_ELEMENTS"]
        for rname in base:
            rgrp = base[rname]
            for bname in rgrp:
                bucket = rgrp[bname]
                if "COLUMN_MAP" not in bucket:
                    continue
                ids = np.asarray(bucket["ID"][...]).reshape(-1)
                cols = lf.decode_columns(bucket["COLUMN_MAP"])
                mn = mx = am = arg = None
                for k, arr in lf.iter_step_slices(bucket):
                    a = np.atleast_2d(np.asarray(arr, dtype=float))
                    if mn is None:
                        mn, mx = a.copy(), a.copy()
                        am = np.abs(a)
                        arg = np.full(a.shape, k, dtype=int)
                    else:
                        mn = np.minimum(mn, a)
                        mx = np.maximum(mx, a)
                        na = np.abs(a)
                        upd = na > am          # strict: first occurrence wins (matches sink)
                        am = np.where(upd, na, am)
                        arg = np.where(upd, k, arg)
                out[f"{rname}/{bname}"] = (ids, mn, mx, am, arg, cols)
    return out


def _col_keys(cols):
    return [(c["gauss_id"], c["section_tag"], c["fiber_id"], c["comp"]) for c in cols]


def check_block(label, stem, out):
    ts_path = os.path.join(out, f"{stem}_ts.ladruno")
    env_path = os.path.join(out, f"{stem}_env.ladruno")
    problems = 0

    probs = lf.validate(env_path)
    if probs:
        print(f"  [FAIL] {label}: envelope file schema violations: {probs}")
        problems += 1
    else:
        print(f"  [OK] {label}: envelope file conformant")

    ts = reduce_elements(ts_path)
    with lf.LadrunoReader(env_path) as r:
        stage = r.stages()[0]
        env = r.envelopes(stage).get("ON_ELEMENTS", {})

    if set(env) != set(ts):
        print(f"  [FAIL] {label}: envelope buckets {set(env)} != time-series {set(ts)}")
        return problems + 1, 0

    total = 0
    for name, (ids, mn, mx, am, arg, cols_ts) in ts.items():
        e = env[name]
        ok = True
        if not np.array_equal(np.asarray(e["id"]).reshape(-1), ids):
            print(f"  [FAIL] {label}/{name}: ID order differs ts vs env")
            ok = False
            problems += 1
        checks = {
            "MIN": (e["min"], mn), "MAX": (e["max"], mx),
            "ABSMAX": (e["absmax"], am), "ARG_STEP": (e["arg_step"], arg),
        }
        for key, (got, want) in checks.items():
            got = np.atleast_2d(np.asarray(got))
            want = np.atleast_2d(want)
            if got.shape != want.shape or not np.allclose(got, want, atol=TOL):
                print(f"  [FAIL] {label}/{name}/{key}: envelope != time-series reduction")
                ok = False
                problems += 1
        # COLUMN_MAP round-trip: the per-column gauss/section/fiber/comp structure
        # on the envelope must equal the time-series COLUMN_MAP (item (1)).
        cols_env = e.get("columns")
        if cols_env is None:
            print(f"  [FAIL] {label}/{name}: envelope has no COLUMN_MAP")
            ok = False
            problems += 1
        elif _col_keys(cols_env) != _col_keys(cols_ts):
            print(f"  [FAIL] {label}/{name}: COLUMN_MAP differs ts vs env")
            ok = False
            problems += 1
        if ok:
            n = mn.size
            ncol = mn.shape[1] if mn.ndim == 2 else 1
            print(f"  [OK] {label}/{name}: MIN/MAX/ABSMAX/ARG_STEP match reduction + "
                  f"COLUMN_MAP {ncol} cols ({n} values, "
                  f"range [{mn.min():.3e}, {mx.max():.3e}])")
            total += n
    return problems, total


def main() -> int:
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    problems = 0
    total = 0
    for label, stem in BLOCKS:
        p, t = check_block(label, stem, out)
        problems += p
        total += t

    print("\n" + "=" * 56)
    if problems == 0:
        print(f"ELEMENT_ENVELOPE_CHECK: ALL PASS ({total} values vs reduction, "
              f"COLUMN_MAP round-trips)")
        return 0
    print(f"ELEMENT_ENVELOPE_CHECK: {problems} PROBLEM(S)")
    return 1


if __name__ == "__main__":
    sys.exit(main())
