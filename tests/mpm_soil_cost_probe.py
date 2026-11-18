"""WP-232k R row: cost per particle-step at particles, with the strain increment particles see.

Runs the K0p column in mode F on bspline2 (capped, as the K0p records) for the elastic, DP,
SANISAND and HS legs through the K0 harness's probe hook (called after every 100-step chunk):
- wall time per chunk, by phase -> seconds per particle-step (soil + floor-block + cap particles);
- every 20th chunk of the pre-shear and the translation, ONE extra explicit step bracketed by
  particle F reads: the per-step strain increment |sym(F_{n+1} F_n^-1) - I| per soil particle
  (median and 90th percentile). The extra steps change the run by < 0.1 % of its steps; this is
  a probe, not a meter.
The material cost of a leg is its time per particle-step minus the elastic leg's (same engine
work); the SANISAND / DP ratio of those is the number M0b's 341 (at 2.56e-4) and 647 (at 9.44e-4)
are per plastic call.

    python tests/mpm_soil_cost_probe.py --out <json> [--legs elastic dp sanisand hs]
"""
from __future__ import annotations

import argparse
import json
import math
import os
import statistics
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ladruno_mpm_soil_harness as K  # noqa: E402

PHASES = ("preshear", "translate")


def _inc(F0, F1, N):
    """|sym(F1 F0^-1) - I| (Frobenius) per particle, 2D."""
    out = []
    for k in range(N):
        a, b, c, d = F0[4 * k:4 * k + 4]
        det = a * d - b * c
        inv = (d / det, -b / det, -c / det, a / det)
        e, f, g, h = F1[4 * k:4 * k + 4]
        m = (e * inv[0] + f * inv[2], e * inv[1] + f * inv[3],
             g * inv[0] + h * inv[2], g * inv[1] + h * inv[3])
        exx, eyy, exy = m[0] - 1.0, m[3] - 1.0, 0.5 * (m[1] + m[2])
        out.append(math.sqrt(exx * exx + eyy * eyy + 2.0 * exy * exy))
    return out


def run_leg(ops, leg):
    dt = K.TIMELINE["F"]["translate"][1]
    st = {"t": None, "n": 0, "chunks": {}, "inc": {p: [] for p in PHASES}, "N": None,
          "Ntot": None}

    def probe(ops_, phase):
        now = time.perf_counter()
        if st["N"] is None:
            body = K._as_list(ops_.mpm("particleData", "body"))
            st["Ntot"] = len(body)
            st["N"] = sum(1 for b in body if b == 1.0)
        if st["t"] is not None:
            c = st["chunks"].setdefault(phase, [0.0, 0])
            c[0] += now - st["t"]
            c[1] += 1
        st["n"] += 1
        if phase in PHASES and st["n"] % 20 == 0:
            N = st["N"]
            F0 = K._as_list(ops_.mpm("particleData", "F"))
            ops_.mpm("analyze", 1, dt)
            F1 = K._as_list(ops_.mpm("particleData", "F"))
            st["inc"][phase] += _inc(F0, F1, N)
        st["t"] = time.perf_counter()

    w0 = time.time()
    rec = K.run(ops, leg, "F", "bspline2", translate=True, capped=True, k0p=True, probe=probe)
    out = {"leg": leg, "particles_soil": st["N"], "particles_total": st["Ntot"],
           "wall_s": time.time() - w0, "meter": rec.get("meter"), "phases": {}}
    for ph, (sec, nch) in st["chunks"].items():
        out["phases"][ph] = {"chunks": nch, "s_per_particle_step": sec / (nch * 100 * st["Ntot"])}
    for ph in PHASES:
        v = sorted(st["inc"][ph])
        if v:
            out["phases"].setdefault(ph, {}).update(
                inc_median=statistics.median(v), inc_p90=v[int(0.9 * (len(v) - 1))],
                inc_samples=len(v))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--legs", nargs="+", default=["elastic", "dp", "sanisand", "hs"])
    ap.add_argument("--out")
    a = ap.parse_args(argv)
    from _testbed import ops
    res = {"build": ops.ladrunoBuild(), "case": "K0p, mode F, bspline2, capped", "legs": {}}
    for leg in a.legs:
        res["legs"][leg] = run_leg(ops, leg)
        print(leg, json.dumps(res["legs"][leg]["phases"]), flush=True)
    E = res["legs"].get("elastic")
    if E:
        res["material_cost"] = {}
        for ph in PHASES:
            e = E["phases"].get(ph, {}).get("s_per_particle_step")
            row = {}
            for leg, L in res["legs"].items():
                s = L["phases"].get(ph, {}).get("s_per_particle_step")
                if s is not None and e is not None:
                    row[leg] = s - e
            if row.get("dp") and row.get("sanisand") is not None:
                row["sanisand_over_dp"] = row["sanisand"] / row["dp"] if row["dp"] > 0 else None
            if row.get("dp") and row.get("hs") is not None:
                row["hs_over_dp"] = row["hs"] / row["dp"] if row["dp"] > 0 else None
            res["material_cost"][ph] = row
        print(json.dumps(res["material_cost"], indent=1))
    if a.out:
        with open(a.out, "w", encoding="utf-8") as f:
            json.dump(res, f, indent=1)


if __name__ == "__main__":
    main()
