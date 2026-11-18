"""WP-232k harness: the K0 SANISAND cell-crossing meter (interface tests/mpm_soil_interface.md
D17-D23, gate K0 of wp-briefs.md "232k" and its freeze).

Written from the frozen interface and the cases file's reference blocks only
(tests/mpm_soil_cases.txt `k0_ref_column`); it never reads tests/mpm_soil_reference.py. The
reference is the JSON record those blocks name: legs.<mat>.after_hold.particles, one entry per
particle row (Y, sxx, syy, szz, sxy), read per particle by its INITIAL Y.

The case (D19): kN, t, m, s, kPa; rho 2.0, g 9.81, l0 10, column two cells wide x in [0, 1],
y in [0, 10], h 0.5, ppc 2 (160 particles); grid origin (0, -3.5), cells 2 x 28; `fix xmin x`,
`fix xmax x`, no `fix ymin`; the floor is a `gridvel ... profile ... translate y` (D17) carrying
the column down D = 2.5 m (5 cells). The twin runs the identical timeline with v = 0.

Mode A (D20): quasi-static LoadControl; gravity ramp over 10 steps of 0.1; the stage flip; a
2-step hold; the translation at v = 1, dt = h/8 = 0.0625 (40 steps, Tr = 0); a 2-step hold.
Mode F (D21): local damping 0.8, FLIP 1, dt 5e-4; gravity ramp 1.0 s, settle 0.5 s, the flip,
hold 0.25 s, the translation at v = 0.1 with Tr = 2 s ramps (t1 = t0 + 27 s), settle 0.5 s:
58,500 steps. Its inertia ratio |dP_y/dt| / W (100-step differences of `mpm info summary`) is
gated <= 1 % at the measurement; its maximum over the run is reported.

The meter: m = RMS over particles of |sigma_p - sigma_ref(Y_p)| / (rho g l0) over xx, yy, xy and
zz, after the full translation (mode F: after the ramp-down and settle).

Usage (the orchestrator launches the long legs):
    python tests/ladruno_mpm_soil_harness.py k0 --mode A --basis gimp [--mat sanisand ...]
writes tests/data/mpm/soil/k0_meter_<mode>_<basis>[_<mat>].json.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CASES = REPO / "tests" / "mpm_soil_cases.txt"
REF_JSON = REPO / "tests" / "data" / "mpm" / "soil" / "k0_column.json"
OUT_DIR = REPO / "tests" / "data" / "mpm" / "soil"
HOST = 232
TS_TAG = 9001

RHO, G, L0, W = 2.0, 9.81, 10.0, 1.0
RHO_G_L0 = RHO * G * L0                          # 196.2 kPa
H = 0.5
ORIGIN = (0.0, -3.5)
CELLS = (2, 28)
D_TRANSLATE = 2.5

# leg -> (material tag in k0_ref_column, stage before gravity, stage at the flip)
# (orchestrator ruling: `mpm stage` passes s through; DP is 0 / 2, SANISAND and HS 0 / 1)
LEGS = {
    "elastic": (1, None, None),
    "dp": (2, 0, 2),
    "sanisand": (3, 0, 1),
    "hs": (4, 0, 1),
    "dp_mildcone": (6, 0, 2),
}
GATED = ("elastic", "dp", "sanisand", "hs")
# the material floor in p (kPa) a particle-step is counted at (K5 / D23): SANISAND -Pmin
# default 1e-3 P_atm, HS p_min = max(1e-6, 1e-3 p_ref); the elastic and DP legs have none
P_FLOOR = {"sanisand": 1e-3 * 101.0, "hs": 1e-3 * 100.0}

TIMELINE = {
    "A": dict(gravity=(10, 0.1), hold1=(2, 0.1), translate=(40, 0.0625), hold2=(2, 0.1),
              v=1.0, Tr=0.0, Tg=1.0),
    "F": dict(gravity=(3000, 5e-4), hold1=(500, 5e-4), translate=(54000, 5e-4), hold2=(1000, 5e-4),
              v=0.1, Tr=2.0, Tg=1.0),
}


# --------------------------------------------------------------------------- #
#  The reference blocks                                                        #
# --------------------------------------------------------------------------- #
def reference_materials(case="k0_ref_column"):
    """{tag: (kind, [tokens])} from the cases file's reference block."""
    mats, inside = {}, False
    for raw in CASES.read_text(encoding="utf-8").splitlines():
        t = raw.split("#", 1)[0].split()
        if not t:
            continue
        if t[0] == "case":
            inside = (t[1] == case)
            continue
        if inside and t[0] == "end":
            break
        if inside and t[0] == "material":
            mats[int(t[1])] = (t[2], t[3:])
    if not mats:
        raise RuntimeError(f"no material lines in {case} ({CASES})")
    return mats


def material_command(tag, kind, tokens):
    def num(x):
        try:
            return int(x) if x.lstrip("-").isdigit() else float(x)
        except ValueError:
            return x
    args = [num(x) for x in tokens]
    name = {"elastic": "ElasticIsotropic", "druckerprager": "DruckerPrager",
            "sanisand": "LadrunoSANISAND", "hs": "LadrunoHS"}[kind]
    return ["nDMaterial", name, tag, *args]


def reference_rows(leg, path=REF_JSON):
    rec = json.loads(Path(path).read_text(encoding="utf-8"))
    return rec["legs"][leg]["after_hold"]["particles"], rec.get("build")


def _ref_at(rows, Y):
    best = min(rows, key=lambda r: abs(r["Y"] - Y))
    if abs(best["Y"] - Y) > 1e-6:
        raise ValueError(f"no reference row at Y = {Y}")
    return best


# --------------------------------------------------------------------------- #
#  One run                                                                     #
# --------------------------------------------------------------------------- #
def _as_list(v):
    if isinstance(v, (list, tuple)):
        return [float(x) for x in v]
    return [float(v)]


def build(ops, leg, mode, basis, translate=True):
    tag, s0, _ = LEGS[leg]
    tl = TIMELINE[mode]
    mats = reference_materials()
    kind, toks = mats[tag]
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial(*material_command(tag, kind, toks)[1:])
    ops.mpm("grid", HOST, "-origin", *ORIGIN, "-h", H, "-cells", *CELLS, "-basis", basis)
    ops.mpm("mode", "quasiStatic" if mode == "A" else "explicit")
    ops.mpm("body", 1, "box", 0.0, 0.0, W, L0, "-ppc", 2, "-mat", tag, "-rho", RHO)
    ops.mpm("fix", "xmin", "x")
    ops.mpm("fix", "xmax", "x")
    t0 = sum(n * dt for n, dt in (tl["gravity"], tl["hold1"]))
    n, dt = tl["translate"]
    t1 = t0 + n * dt
    v = -tl["v"] if translate else 0.0
    # the floor: the column base row and the row below; it moves with the column (D17)
    ops.mpm("gridvel", "floor", "box", 0.0, -1.0, W, 0.0, "v", "free", v,
            "profile", t0, tl["Tr"], t1, "translate", "y")
    ops.mpm("gravity", 0.0, -G)
    ops.timeSeries("Path", TS_TAG, "-time", 0.0, tl["Tg"], 1.0e9, "-values", 0.0, 1.0, 1.0)
    ops.mpm("timeSeries", TS_TAG)
    if mode == "F":
        ops.mpm("damping", "local", 0.8)
        ops.mpm("flip", 1.0)
    else:
        ops.constraints("Transformation")
        ops.numberer("RCM")
        ops.system("UmfPack")
        ops.test("NormDispIncr", 1.0e-10, 25)
        ops.algorithm("Newton")
        ops.integrator("LoadControl", tl["gravity"][1])
        ops.analysis("Static")
    if s0 is not None:
        ops.mpm("stage", tag, s0)
    return t0, t1


def _state(ops, leg):
    if leg != "sanisand":
        return None
    return _as_list(ops.mpm("particleResponse", "state"))


def _floor_hits(sig, szz, pfloor):
    if pfloor is None:
        return 0
    hits = 0
    for k in range(len(szz)):
        p = -(sig[3 * k] + sig[3 * k + 1] + szz[k]) / 3.0
        if p <= pfloor:
            hits += 1
    return hits


def run(ops, leg, mode, basis, translate=True, chunk=None):
    """One leg; returns the record (meter, crossings, inertia, floor hits, state)."""
    tag, _, s1 = LEGS[leg]
    tl = TIMELINE[mode]
    wall0 = time.time()
    build(ops, leg, mode, basis, translate)
    x0 = _as_list(ops.mpm("particleData", "x"))
    N = len(x0) // 2
    Y = [x0[2 * k + 1] for k in range(N)]
    cell = [math.floor((Y[k] - ORIGIN[1]) / H) for k in range(N)]
    crossings = [0] * N
    floor_hits = 0
    pfloor = P_FLOOR.get(leg)
    W_col = RHO * G * L0 * W
    inertia, inertia_max, mom_prev, t_prev = None, 0.0, None, None
    chunk = chunk or (1 if mode == "A" else 100)

    def sample(phase):
        nonlocal floor_hits, mom_prev, t_prev, inertia, inertia_max
        x = _as_list(ops.mpm("particleData", "x"))
        for k in range(N):
            c = math.floor((x[2 * k + 1] - ORIGIN[1]) / H)
            if c != cell[k]:
                crossings[k] += abs(c - cell[k])
                cell[k] = c
        if pfloor is not None:
            floor_hits += _floor_hits(_as_list(ops.mpm("particleData", "sig")),
                                      _as_list(ops.mpm("particleData", "sigZZ")), pfloor)
        if mode == "F":
            s = _as_list(ops.mpm("info", "summary"))
            t = float(ops.mpm("info", "time"))
            if mom_prev is not None and t > t_prev:
                inertia = abs((s[2] - mom_prev) / (t - t_prev)) / W_col
                if phase in ("translate", "hold2"):
                    inertia_max = max(inertia_max, inertia)
            mom_prev, t_prev = s[2], t

    def advance(n, dt, phase):
        done = 0
        while done < n:
            k = min(chunk, n - done)
            if ops.mpm("analyze", k, dt) != 0:
                raise RuntimeError(f"{leg} {mode} {basis}: mpm analyze failed in {phase}")
            done += k
            sample(phase)

    advance(*tl["gravity"], "gravity")
    if s1 is not None:
        ops.mpm("stage", tag, s1)
    advance(*tl["hold1"], "hold1")
    state_before = _state(ops, leg)
    advance(*tl["translate"], "translate")
    advance(*tl["hold2"], "hold2")
    sig = _as_list(ops.mpm("particleData", "sig"))
    szz = _as_list(ops.mpm("particleData", "sigZZ"))
    x = _as_list(ops.mpm("particleData", "x"))
    rows, ref_build = reference_rows(leg)
    sq = 0.0
    for k in range(N):
        r = _ref_at(rows, Y[k])
        d = [sig[3 * k] - r["sxx"], sig[3 * k + 1] - r["syy"], sig[3 * k + 2] - r["sxy"],
             szz[k] - r["szz"]]
        sq += sum(e * e for e in d)
    meter = math.sqrt(sq / N) / RHO_G_L0
    drop = sum(x[2 * k + 1] - Y[k] for k in range(N)) / N
    rec = dict(leg=leg, mode=mode, basis=basis, translate=translate, particles=N,
               meter=meter, mean_displacement_y=drop,
               crossings=dict(min=min(crossings), max=max(crossings),
                              mean=sum(crossings) / N, per_particle=crossings),
               floor_hits=floor_hits, p_floor=pfloor,
               state_before=state_before, state_after=_state(ops, leg),
               ref_build=ref_build, wall_s=time.time() - wall0)
    if mode == "F":
        rec["inertia_ratio"] = inertia
        rec["inertia_ratio_max"] = inertia_max
    return rec


def state_drift(rt, rz):
    """SANISAND: per-particle |e_t - e_twin| and |alpha_t - alpha_twin| (state: alpha at 6-11,
    e at 24, 26 values per particle), max and RMS."""
    a, b = rt.get("state_after"), rz.get("state_after")
    if not a or not b:
        return None
    n = len(a) // 26
    de = [abs(a[26 * k + 24] - b[26 * k + 24]) for k in range(n)]
    da = [math.sqrt(sum((a[26 * k + i] - b[26 * k + i]) ** 2 for i in range(6, 12)))
          for k in range(n)]
    return dict(de_max=max(de), de_rms=math.sqrt(sum(v * v for v in de) / n),
                dalpha_max=max(da), dalpha_rms=math.sqrt(sum(v * v for v in da) / n))


def k0(ops, mode, basis, legs=None):
    legs = legs or list(LEGS)
    out = dict(case="K0", mode=mode, basis=basis, build=str(ops.ladrunoBuild()),
               rho_g_l0=RHO_G_L0, legs={})
    for leg in legs:
        rt = run(ops, leg, mode, basis, translate=True)
        rz = run(ops, leg, mode, basis, translate=False)
        out["legs"][leg] = dict(translated=rt, twin=rz, state_drift=state_drift(rt, rz))
    L = out["legs"]
    if "sanisand" in L and "dp" in L:
        ms, md = L["sanisand"]["translated"]["meter"], L["dp"]["translated"]["meter"]
        out["A_sanisand"] = ms / md if md > 0 else float("inf")
        out["gate"] = dict(rule="A <= 2 or m(SANISAND) <= 0.01",
                           passed=bool(out["A_sanisand"] <= 2.0 or ms <= 0.01))
    for leg in ("elastic", "hs", "dp_mildcone"):
        if leg in L and "dp" in L:
            md = L["dp"]["translated"]["meter"]
            out[f"A_{leg}"] = L[leg]["translated"]["meter"] / md if md > 0 else float("inf")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["k0"])
    ap.add_argument("--mode", choices=["A", "F"], required=True)
    ap.add_argument("--basis", choices=["gimp", "bspline2"], required=True)
    ap.add_argument("--mat", action="append", choices=list(LEGS))
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from _testbed import ops   # noqa: E402
    res = k0(ops, a.mode, a.basis, a.mat)
    suffix = "" if not a.mat else "_" + "_".join(a.mat)
    path = Path(a.out) if a.out else OUT_DIR / f"k0_meter_{a.mode}_{a.basis}{suffix}.json"
    path.write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(json.dumps({k: v for k, v in res.items() if k != "legs"}, indent=1))
    for leg, r in res["legs"].items():
        print(leg, "meter", r["translated"]["meter"], "twin", r["twin"]["meter"],
              "crossings", r["translated"]["crossings"], "floor_hits",
              r["translated"]["floor_hits"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
