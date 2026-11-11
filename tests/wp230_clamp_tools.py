"""WP-230 helpers (not collected: no `test_` prefix).

FIXTURE   tests/data/wp230_oracle_paths.json, the oracle's clamp paths
          (Ladruno_files/testbed/sanisand_lowp_clamp/clamp_oracle.py --fixture).
BASELINE  tests/data/wp230_cutoff_byteid_baseline.json: WP-152's three element legs with
          -sasTensionCutoff ("cutoff") and with R1 alone ("none"), recorded on a build
          WITHOUT WP-230 by

              py -3.12 -S Ladruno_scripts/ci_run_pytest.py -- wp230_clamp_tools.py::record -q

          (or `python wp230_clamp_tools.py --record` with the engine on the path).
          A row per step: rc, then float.hex of stress(6), alpha(6), alpha_in(6) of Gauss
          point 1, and the first 44 sasStats columns (the 44 that exist without WP-230; this
          WP appends two) and the first 12 sasOptions values.
"""
from __future__ import annotations

import json
import os

from _testbed import ops
import wp152_cutoff_tools as T

HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURE = os.path.join(HERE, "data", "wp230_oracle_paths.json")
BASELINE = os.path.join(HERE, "data", "wp230_cutoff_byteid_baseline.json")
LEGS = {"cutoff": T.BASE + T.R1 + T.CUTOFF, "none": T.BASE + T.R1}


def _hex(v):
    return [float(x).hex() for x in v]


def byteid_rows(leg):
    """name -> rows, for WP-152's iso / te / cyc paths. Without the cutoff ("none") the path
    stops at its first refusal; the refusal itself is part of the row."""
    out = {}
    for name, incs in T.paths().items():
        rows = []

        def per_step(k, rec):
            rows.append([rec["rc"]] + _hex(T.mresp("stress")) + _hex(T.mresp("alpha"))
                        + _hex(T.mresp("alpha_in")) + _hex(rec["sas"][:44])
                        + _hex(T.mresp("sasOptions")[:12]))
        T.run(incs, LEGS[leg], per_step=per_step)
        out[name] = rows
    return out


def run_aniso(incs, opts, e_stage0):
    """As T.run, but stage 0 (elastic) ends at the ANISOTROPIC strain e_stage0 = (ex, ey, ez),
    compression positive, so the post-flip state carries a deviator (a K0-like state)."""
    orig = T._series

    def series(tag0, incs_, ev0=None):
        hist = [[-e_stage0[i] * (k + 1) / T.N0 for i in range(3)] for k in range(T.N0)]
        cur = list(hist[-1])
        for de in incs_:
            cur = [cur[i] - de[i] for i in range(3)]
            hist.append(list(cur))
        for i in range(3):
            ops.timeSeries("Path", tag0 + i, "-dt", 1.0, "-values", 0.0, *[h[i] for h in hist], hist[-1][i])
        return hist
    T._series = series
    try:
        return T.run(incs, opts)
    finally:
        T._series = orig


def record(path=BASELINE):
    data = dict(note="WP-230 byte-identity baseline: WP-152's element legs on a build WITHOUT WP-230 "
                     "(wp230_clamp_tools.record)", build=ops.ladrunoBuild().strip().splitlines()[0],
                rows={leg: byteid_rows(leg) for leg in LEGS})
    ops.wipe()
    json.dump(data, open(path, "w"))
    return {leg: {n: len(r) for n, r in rows.items()} for leg, rows in data["rows"].items()}


if __name__ == "__main__":
    import sys
    if "--record" in sys.argv:
        print(record())
