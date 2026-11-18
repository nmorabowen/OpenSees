"""WP-232k R row (reported, not gated): the granular-column collapse (Lube et al. 2004, 2005; the
ADR §8 validation row) with soil particles, mode F on bspline2, against DruckerPrager.

Plane strain, half column (symmetry fix on xmin): width a = 1 m, height 2a (aspect ratio 2), on
a no-slip base (fix ymin xy; `--base wall` gives M0b's Coulomb wall, mu 0.4). Grid h = 0.25 m,
x in [0, 8] m, ppc 2 (128 particles), `bspline2`. Materials: the k0_ref_column decks of
tests/mpm_soil_cases.txt (kN, t, m, s, kPa; rho 2.0, g 9.81) -- SANISAND (tag 3, stage 0 -> 1),
HS (4, 0 -> 1), the gated DP (2, 0 -> 2) and dp_mildcone (6, 0 -> 2).

Timeline (the GATE staging, default): gravity and a smooth gate on the column's free face, both
as per-particle body forces (the gate: a lateral pressure K_GATE rho g (2a - y) on the outer
particle column, per unit mass); stage 0 (linear elastic in every soil model) with local damping
0.8 for SETTLE s; the flip under the gate; GATE_HOLD s more, damped; then the gate is released
(gravity only) and the collapse runs UNDAMPED (local damping is not Galilean invariant) for
COLLAPSE s. dt = 0.25 dt_c, re-read before every 0.05 s chunk (the estimate moves with the flow;
a fixed dt was refused by the CFL check). `--no-gate` is the first staging (flip on the
unconfined elastic column), kept for the record: on it SANISAND and HS do not collapse --
SANISAND's stage switch raises M_c at every point whose stress ratio lies outside the bounding
surface (vanilla ManzariDafalias, "Outside Bounding!"; face q/p up to 2.8 against M_c 1.33), and
HS seeds its hardening from the committed stress and carries c = 10 kPa (a vertical cut stands
to 4 c / gamma tan(45 + phi/2) = 4.2 m > 2 m) -- while DP's flip returns the stress to its yield
surface and the column falls.

Recorded: runout (L_inf - L_0) / a with L the front particle edge (centre + half the particle
spacing; 0 for a column that does not move), the final height H_inf / 2a (top particle edge,
inner column), a history (every 0.05 s: time, phase, front, top, kinetic energy, max |v|), volume
checks (sum V_p, det F range) at the flip, the release and the end, floor hits (particle-samples
at p <= the material's p floor, K5 / D23, from the release on), flip resets (particles at p below
the floor at the flip, which SANISAND's Elastic2Plastic resets to -Pmin I), the largest face
stress ratio q/p at the flip, SANISAND sasStats when its scheme is 129, and a refusal abort (W17)
with its time -- an acceptable outcome, recorded.

    python tests/mpm_soil_collapse.py --leg sanisand --wall-cap 7200 --out <json>
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ladruno_mpm_soil_harness as K  # noqa: E402
import ladruno_mpm_soil_s_harness as S  # noqa: E402

A, HC = 1.0, 2.0                  # half width, height
H = 0.25
SP = H / 2.0                      # particle spacing at ppc 2
CELLS = (32, 10)                  # x in [0, 8], y in [0, 2.5]
MU = 0.4
SETTLE, GATE_HOLD, COLLAPSE, SAMPLE = 1.0, 0.5, 4.0, 0.05
K_GATE = 0.5
DT_FACTOR = 0.25
LEGS = {"sanisand": (3, 0, 1), "hs": (4, 0, 1), "dp": (2, 0, 2), "dp_mildcone": (6, 0, 2)}


def _p(sig, szz, k):
    return -(sig[3 * k] + sig[3 * k + 1] + szz[k]) / 3.0


def _q(sig, szz, k):
    sx, sy, sz, txy = sig[3 * k], sig[3 * k + 1], szz[k], sig[3 * k + 2]
    return math.sqrt(0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2) + 3.0 * txy * txy)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--leg", choices=sorted(LEGS), required=True)
    ap.add_argument("--wall-cap", type=float, default=7200.0)
    ap.add_argument("--settle", type=float, default=SETTLE)
    ap.add_argument("--gate-hold", type=float, default=GATE_HOLD)
    ap.add_argument("--collapse", type=float, default=COLLAPSE)
    ap.add_argument("--no-gate", action="store_true",
                    help="flip on the unconfined elastic column (the first staging)")
    # the default is the no-slip base: with M0b's Coulomb wall (mu 0.4) a base particle's centre
    # reaches y = 0 and the next step is refused as leaving the grid (A7) -- DP at t = 2.08 s,
    # f235a001f; kept as an option
    ap.add_argument("--base", choices=["wall", "rough"], default="rough",
                    help="rough: fix ymin xy (no slip); wall: Coulomb wall mu 0.4 on ymin")
    ap.add_argument("--out")
    a = ap.parse_args(argv)
    gate = not a.no_gate
    from _testbed import ops
    L = K._as_list
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
    x0 = L(ops.mpm("particleData", "x"))
    m0 = L(ops.mpm("particleData", "m"))
    N = len(m0)
    front0 = max(x0[0::2]) + 0.5 * SP

    def body_force(with_gate):
        """gravity on every particle (+ the gate's lateral pressure on the outer column)."""
        bf = []
        for k in range(N):
            fx = 0.0
            if with_gate and x0[2 * k] > A - SP:
                fx = -K_GATE * K.RHO * K.G * (HC - x0[2 * k + 1]) * SP / m0[k]
            bf += [fx, -K.G]
        ops.mpm("particleBodyForce", *bf)

    if gate:
        body_force(True)
    else:
        ops.mpm("gravity", 0.0, -K.G)
    ops.mpm("damping", "local", 0.8)
    ops.mpm("flip", 1.0)
    ops.mpm("stage", tag, s0)
    pfloor = K.P_FLOOR.get(a.leg)
    out = {"leg": a.leg, "material": [kind, toks], "build": str(ops.ladrunoBuild()),
           "geometry": dict(a=A, height=HC, h=H, cells=CELLS, ppc=2, basis="bspline2",
                            base=(f"wall ymin mu {MU}" if a.base == "wall" else "fix ymin xy"),
                            symmetry="xmin x"),
           "staging": (dict(kind="gate", k_gate=K_GATE, settle=a.settle, gate_hold=a.gate_hold)
                       if gate else dict(kind="no gate (flip on the unconfined column)",
                                         settle=a.settle)),
           "timeline": dict(collapse=a.collapse, sample=SAMPLE, dt_factor=DT_FACTOR),
           "p_floor": pfloor, "history": [], "stopped": None, "floor_hits": 0, "checks": {}}
    dt = DT_FACTOR * float(L(ops.mpm("criticalTimeStep"))[0])
    out["dt_settle"] = dt

    def check(tag_):
        V = L(ops.mpm("particleData", "V"))
        F = L(ops.mpm("particleData", "F"))
        J = [F[4 * k] * F[4 * k + 3] - F[4 * k + 1] * F[4 * k + 2] for k in range(N)]
        out["checks"][tag_] = dict(t=float(ops.mpm("info", "time")), sum_V=sum(V),
                                   sum_V0=sum(L(ops.mpm("particleData", "V0"))),
                                   J_min=min(J), J_max=max(J))

    def advance(T, dt0, phase):
        t_end = float(ops.mpm("info", "time")) + T
        while True:
            # the CFL estimate moves with the flow: dt from the current dt_c every chunk
            dt = min(dt0, DT_FACTOR * float(L(ops.mpm("criticalTimeStep"))[0]))
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
                aborted = int(L(ops.mpm("debug", "aborted"))[0])
                out["stopped"] = (f"{'abort (refusal at commit, W17)' if aborted else 'failed step'}"
                                  f" in {phase} near t = {ops.mpm('info', 'time')}: {exc}")
                return False
            x = L(ops.mpm("particleData", "x"))
            v = L(ops.mpm("particleData", "v"))
            if pfloor is not None and phase == "collapse":
                sig = L(ops.mpm("particleData", "sig"))
                szz = L(ops.mpm("particleData", "sigZZ"))
                out["floor_hits"] += sum(1 for k in range(N) if _p(sig, szz, k) <= pfloor)
            s = L(ops.mpm("info", "summary"))
            out["history"].append([float(ops.mpm("info", "time")), phase,
                                   max(x[0::2]) + 0.5 * SP, max(x[1::2]) + 0.5 * SP, s[-1],
                                   max(math.hypot(v[2 * k], v[2 * k + 1]) for k in range(N))])
        return True

    ok = advance(a.settle, dt, "settle")
    if ok:
        sig = L(ops.mpm("particleData", "sig"))
        szz = L(ops.mpm("particleData", "sigZZ"))
        x = L(ops.mpm("particleData", "x"))
        out["flip_resets"] = (sum(1 for k in range(N) if _p(sig, szz, k) < pfloor)
                              if pfloor is not None else None)
        out["face_q_over_p_max_at_flip"] = max(
            _q(sig, szz, k) / max(_p(sig, szz, k), 1e-9) for k in range(N) if x[2 * k] > A - SP)
        check("flip")
        st_flip = L(ops.mpm("particleResponse", "state")) if a.leg == "sanisand" else None
        ops.mpm("stage", tag, s1)
        if gate:
            ok = advance(a.gate_hold, dt, "gate_hold")
            body_force(False)                         # the gate is released
        if ok:
            check("release")
            st_rel = L(ops.mpm("particleResponse", "state")) if a.leg == "sanisand" else None
            ops.mpm("damping", "local", 0.0)
            out["dt_collapse"] = min(dt, DT_FACTOR * float(L(ops.mpm("criticalTimeStep"))[0]))
            advance(a.collapse, out["dt_collapse"], "collapse")
            if a.leg == "sanisand":
                st = L(ops.mpm("particleResponse", "state"))
                n = len(st) // N
                out["alpha_change_max"] = {
                    "flip_to_release": max(abs(st_rel[n * k + i] - st_flip[n * k + i])
                                           for k in range(N) for i in range(6, 12)),
                    "release_to_end": max(abs(st[n * k + i] - st_rel[n * k + i])
                                          for k in range(N) for i in range(6, 12))}
    check("end")
    x = L(ops.mpm("particleData", "x"))
    front = max(x[0::2]) + 0.5 * SP
    top = max(x[2 * k + 1] for k in range(N) if x[2 * k] < H) + 0.5 * SP
    out.update(particles=N, t_end=float(ops.mpm("info", "time")), front0=front0, front=front,
               runout=(front - front0) / A, height_final=top / HC, wall_s=time.time() - w0,
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
