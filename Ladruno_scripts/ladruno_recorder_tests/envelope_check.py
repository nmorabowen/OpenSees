"""Envelope (ADR D7) gate -- checker (venv python, h5py).

The envelope must equal the manual componentwise reduction of the time series it
summarizes. Reduce env_ts.ladruno (MIN/MAX/ABSMAX/ARG_STEP over steps) and compare
to the ENVELOPES datasets in env_env.ladruno (1e-9). Also runs the schema
validator (which checks ENVELOPES shapes + ABSMAX==max(|MIN|,|MAX|)).

    python envelope_check.py <out_dir>
"""

from __future__ import annotations

import os
import sys

import numpy as np

import ladruno_format as lf

TOL = 1e-9


def reduce_timeseries(path):
    """rname -> (ids, MIN, MAX, ABSMAX, ARG_STEP) over all ON_NODES steps."""
    out = {}
    with lf.LadrunoReader(path) as r:
        stage = r.stages()[0]
        base = r.f[f"{stage}/RESULTS/ON_NODES"]
        for rname in base:
            grp = base[rname]
            ids = np.asarray(grp["ID"][...]).reshape(-1)
            mn = mx = am = arg = None
            for k, arr in lf.iter_step_slices(grp):
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
            out[rname] = (ids, mn, mx, am, arg)
    return out


def main() -> int:
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    ts_path = os.path.join(out, "env_ts.ladruno")
    env_path = os.path.join(out, "env_env.ladruno")
    problems = 0

    probs = lf.validate(env_path)
    if probs:
        print(f"  [FAIL] envelope file schema violations: {probs}")
        problems += 1
    else:
        print("  [OK] envelope file conformant")

    ts = reduce_timeseries(ts_path)
    # time-series component names per result (for the COMPONENTS round-trip check)
    ts_comps = {}
    with lf.LadrunoReader(r_path := ts_path) as r:
        stage = r.stages()[0]
        base = r.f[f"{stage}/RESULTS/ON_NODES"]
        for rname in base:
            ts_comps[rname] = [c for c in lf._attr(base[rname], "COMPONENTS").split(",") if c]

    with lf.LadrunoReader(env_path) as r:
        stage = r.stages()[0]
        env = r.envelopes(stage).get("ON_NODES", {})

    if set(env) != set(ts):
        print(f"  [FAIL] envelope result names {set(env)} != time-series {set(ts)}")
        return problems + 1

    total = 0
    for rname, (ids, mn, mx, am, arg) in ts.items():
        e = env[rname]
        checks = {
            "MIN": (e["min"], mn), "MAX": (e["max"], mx),
            "ABSMAX": (e["absmax"], am), "ARG_STEP": (e["arg_step"], arg),
        }
        ok = True
        for key, (got, want) in checks.items():
            got = np.atleast_2d(np.asarray(got))
            want = np.atleast_2d(want)
            if got.shape != want.shape or not np.allclose(got, want, atol=TOL):
                print(f"  [FAIL] {rname}/{key}: envelope != time-series reduction")
                ok = False
                problems += 1
        # COMPONENTS round-trip (finding (h)): the envelope group must carry the
        # same authoritative component names as the time-series result.
        env_comps = e.get("components")
        if env_comps != ts_comps.get(rname):
            print(f"  [FAIL] {rname}/COMPONENTS: envelope {env_comps} != "
                  f"time-series {ts_comps.get(rname)}")
            ok = False
            problems += 1
        if ok:
            n = mn.size
            print(f"  [OK] {rname}: MIN/MAX/ABSMAX/ARG_STEP match reduction + "
                  f"COMPONENTS={env_comps} ({n} values, range [{mn.min():.3e}, {mx.max():.3e}])")
            total += n

    print("\n" + "=" * 56)
    if problems == 0:
        print(f"ENVELOPE_CHECK: ALL PASS ({total} values vs time-series reduction)")
        return 0
    print(f"ENVELOPE_CHECK: {problems} PROBLEM(S)")
    return 1


if __name__ == "__main__":
    sys.exit(main())
