"""WP-232k R row (reported, not gated): the granular-column collapse (Lube et al. 2004, 2005; the
ADR §8 validation row) with soil particles, mode F on bspline2, against DruckerPrager.

Plane strain, half column (symmetry fix on xmin): width a = 1 m, height 2a (aspect ratio 2), on
a no-slip base (fix ymin xy; `--base wall` gives M0b's Coulomb wall, mu 0.4). Grid h = 0.25 m,
x in [0, 8] m, ppc 2 (128 particles),
`bspline2`. Materials: the k0_ref_column decks of tests/mpm_soil_cases.txt (kN, t, m, s, kPa;
rho 2.0, g 9.81) -- SANISAND (tag 3, stage 0 -> 1), HS (4, 0 -> 1), the gated DP (2, 0 -> 2) and
dp_mildcone (6, 0 -> 2). Timeline: stage 0 (linear elastic in every soil model) under gravity at
once with local damping 0.8 for SETTLE s, the flip, then the collapse UNDAMPED (local damping is
not Galilean invariant) for COLLAPSE s; dt = 0.25 dt_c, from `mpm criticalTimeStep` before
every 0.05 s chunk (the estimate moves with the flow; a fixed dt was refused by the CFL check).
DP on f235a001f (Windows, 128 particles): runout 3.41 (dp) / 3.44 (dp_mildcone), 9-12 s.

Recorded: runout (L_inf - a) / a from the front particle centre, the final height H_inf / 2a,
a history (every 0.05 s: front, top, kinetic energy), floor hits (particle-samples at p <= the
material's p floor, K5 / D23, from the flip on), flip resets (particles at p below the floor at
the flip, which SANISAND's Elastic2Plastic resets to -Pmin I), SANISAND sasStats when its scheme
is 129, and a refusal abort (W17) with its time -- an acceptable outcome, recorded.

    python tests/mpm_soil_collapse.py --leg sanisand --wall-cap 7200 --out <json>
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ladruno_mpm_soil_harness as K  # noqa: E402
import ladruno_mpm_soil_s_harness as S  # noqa: E402

A, HC = 1.0, 2.0                  # half width, height
H = 0.25
CELLS = (32, 10)                  # x in [0, 8], y in [0, 2.5]
MU = 0.4
SETTLE, COLLAPSE, SAMPLE = 1.0, 4.0, 0.05
DT_FACTOR = 0.25
LEGS = {"sanisand": (3, 0, 1), "hs": (4, 0, 1), "dp": (2, 0, 2), "dp_mildcone": (6, 0, 2)}


def _p(sig, szz, k):
    return -(sig[3 * k] + sig[3 * k + 1] + szz[k]) / 3.0


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--leg", choices=sorted(LEGS), required=True)
    ap.add_argument("--wall-cap", type=float, default=7200.0)
    ap.add_argument("--settle", type=float, default=SETTLE)
    ap.add_argument("--collapse", type=float, default=COLLAPSE)
    # the default is the no-slip base: with M0b's Coulomb wall (mu 0.4) a base particle's centre
    # reaches y = 0 and the next step is refused as leaving the grid (A7) -- DP at t = 2.08 s,
    # f235a001f; kept as an option
    ap.add_argument("--base", choices=["wall", "rough"], default="rough",
                    help="rough: fix ymin xy (no slip); wall: Coulomb wall mu 0.4 on ymin")
    ap.add_argument("--out")
    a = ap.parse_args(argv)
    from _testbed import ops
    tag, s0, s1 = LEGS[a.leg]
    kind, toks = K.reference_materials()[tag]
    w0 = time.time()
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial(*K.material_command(tag, kind, toks)[1:])
    ops.mpm("grid", K.HOST, "-origin", 0.0, 0.0, "-h", H, "-cells", *CELLS, "-basis", "bspline2")
    ops.mpm("mode", "explicit")
    ops.mpm("body", 1, "box", 0.0, 0.0, A, HC, "-ppc", 2, "-mat", tag, "-rho", K.RHO)
    ops.mpm("fix", "xmin", "x")
    if a.base == "rough":
        ops.mpm("fix", "ymin", "xy")
    else:
        ops.mpm("wall", "ymin", MU)
    ops.mpm("gravity", 0.0, -K.G)
    ops.mpm("damping", "local", 0.8)
    ops.mpm("flip", 1.0)
    ops.mpm("stage", tag, s0)
    pfloor = K.P_FLOOR.get(a.leg)
    out = {"leg": a.leg, "material": [kind, toks], "build": str(ops.ladrunoBuild()),
           "geometry": dict(a=A, height=HC, h=H, cells=CELLS, ppc=2, basis="bspline2",
                            base=(f"wall ymin mu {MU}" if a.base == "wall" else "fix ymin xy"),
                            symmetry="xmin x"),
           "timeline": dict(settle=a.settle, collapse=a.collapse, sample=SAMPLE,
                            dt_factor=DT_FACTOR),
           "p_floor": pfloor, "history": [], "stopped": None, "floor_hits": 0}
    dt = DT_FACTOR * float(K._as_list(ops.mpm("criticalTimeStep"))[0])
    out["dt_settle"] = dt

    def advance(T, dt0, phase):
        t_end = float(ops.mpm("info", "time")) + T
        while True:
            # the CFL estimate moves with the flow: dt from the current dt_c every chunk
            dt = min(dt0, DT_FACTOR * float(K._as_list(ops.mpm("criticalTimeStep"))[0]))
            out["dt_min"] = min(out.get("dt_min", dt), dt)
            if float(ops.mpm("info", "time")) >= t_end - 0.5 * dt:
                break
            n_chunk = max(1, int(round(SAMPLE / dt)))
            if time.time() - w0 > a.wall_cap:
                out["stopped"] = f"wall cap {a.wall_cap:.0f} s at t = {ops.mpm('info', 'time')}"
                return False
            try:
                ops.mpm("analyze", n_chunk, dt)
            except Exception as exc:
                aborted = int(K._as_list(ops.mpm("debug", "aborted"))[0])
                out["stopped"] = (f"{'abort (refusal at commit, W17)' if aborted else 'failed step'}"
                                  f" in {phase} near t = {ops.mpm('info', 'time')}: {exc}")
                return False
            x = K._as_list(ops.mpm("particleData", "x"))
            N = len(x) // 2
            sig = K._as_list(ops.mpm("particleData", "sig"))
            szz = K._as_list(ops.mpm("particleData", "sigZZ"))
            if pfloor is not None and phase == "collapse":
                out["floor_hits"] += sum(1 for k in range(N) if _p(sig, szz, k) <= pfloor)
            s = K._as_list(ops.mpm("info", "summary"))
            out["history"].append([float(ops.mpm("info", "time")), phase,
                                   max(x[2 * k] for k in range(N)),
                                   max(x[2 * k + 1] for k in range(N)), s[-1]])
        return True

    ok = advance(a.settle, dt, "settle")
    if ok:
        sig = K._as_list(ops.mpm("particleData", "sig"))
        szz = K._as_list(ops.mpm("particleData", "sigZZ"))
        N = len(szz)
        out["flip_resets"] = (sum(1 for k in range(N) if _p(sig, szz, k) < pfloor)
                              if pfloor is not None else None)
        ops.mpm("stage", tag, s1)
        ops.mpm("damping", "local", 0.0)
        dtc = float(K._as_list(ops.mpm("criticalTimeStep"))[0])
        dt = min(dt, DT_FACTOR * dtc)
        out["dt_collapse"] = dt
        advance(a.collapse, dt, "collapse")
    x = K._as_list(ops.mpm("particleData", "x"))
    N = len(x) // 2
    front = max(x[2 * k] for k in range(N))
    top = max(x[2 * k + 1] for k in range(N) if x[2 * k] < H)
    out.update(particles=N, t_end=float(ops.mpm("info", "time")), front=front,
               runout=(front - A) / A, height_final=top / HC, wall_s=time.time() - w0,
               reported="runout against DP; not gated")
    if a.leg == "sanisand":
        try:
            out["sas_counts"] = S.sas_counts(ops)
        except Exception as exc:                     # scheme 1 decks carry no SAS-ME stats
            out["sas_counts"] = f"n/a: {exc}"
    print(json.dumps({k: v for k, v in out.items() if k != "history"}, indent=1))
    if a.out:
        with open(a.out, "w", encoding="utf-8") as f:
            json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
