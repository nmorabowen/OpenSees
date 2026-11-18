"""WP-232k harness: the K0 SANISAND cell-crossing meter (interface tests/mpm_soil_interface.md
D17-D23, gate K0 of wp-briefs.md "232k" and its freeze).

Written from the frozen interface and the cases file's reference blocks only
(tests/mpm_soil_cases.txt `k0_ref_column`); it never reads tests/mpm_soil_reference.py. The
reference is the JSON record those blocks name: legs.<mat>.after_hold.particles, one entry per
particle row (Y, sxx, syy, szz, sxy), read per particle by its INITIAL Y.

The case (D19): kN, t, m, s, kPa; rho 2.0, g 9.81, l0 10, column two cells wide x in [0, 1],
y in [0, 10], h 0.5, ppc 2 (160 particles); `fix xmin x`, `fix xmax x`, no `fix ymin`; the
column moves down D = 2.5 m (5 cells). The twin runs the identical timeline with v = 0.

The floor (amendment D19a, proposed by 232k-C before any soil mode-F number): a bonded elastic
block y in [-1.5, 0] (the elastic leg's E, nu, rho; body 2, excluded from the meter) under the
column; grid origin (0, -5), cells 2 x 31. The frozen D19 floor (a gridvel box directly under
the column) put the soil's bottom row on GIMP tails each time a node row left the moving set, and
released the full column reaction at that row: mode F gimp elastic 0.188 (twin 0.015), bottom-row
syy error -0.64 rho g l0, inertia max 0.208, drop -2.522 (the record of the diagnosis).
  Mode A ("block"): the `gridvel ... profile ... translate y` (D17) prescribes rows 2-4 cells
  below the soil, inside the block, so the set's entries and exits never touch a soil support.
  Mode F ("float"): no grid condition; the block carries the whole weight as an upward body
  force (`mpm particleBodyForce`, sum of external forces 0) and the translation is the uniform
  body acceleration -v phi'(t) on every particle, with local damping off during it (FLAC local
  damping is not Galilean invariant and would brake the acceleration itself; it is 0.8 in every
  other phase). Gravity is applied at once in stage 0 (elastic in every leg) and settles.

Mode A (D20): quasi-static LoadControl; gravity ramp over 10 steps of 0.1; the stage flip; a
2-step hold; the translation at v = 1, dt = h/8 = 0.0625 (40 steps, Tr = 0); a 2-step hold.
Mode F (D21): local damping 0.8, FLIP 1, dt 5e-4; gravity (at once) and settle 1.5 s, the flip,
hold 0.25 s, the translation at v = 0.1 with Tr = 2 s ramps (t1 = t0 + 27 s), settle 0.5 s:
58,500 steps. Its inertia ratio |dP_y/dt| / (M g) (100-step differences of `mpm info summary`)
is gated <= 1 % at the measurement; its maximum over the run is reported.

Floor counts (K5 / D23): particle-steps with p <= the material's floor, sampled after every
chunk, counted from the flip on (`floor_hits`); the gravity stage's count is reported apart
(`floor_hits_stage0`) because stage 0 of both soil models is pressure independent.

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
REF_JSON_CAPPED = REPO / "tests" / "data" / "mpm" / "soil" / "k0_column_capped.json"
CAP_TAG = 98                                     # D23: the bonded elastic cap (body 3)
OUT_DIR = REPO / "tests" / "data" / "mpm" / "soil"
HOST = 232
TS_TAG = 9001

RHO, G, L0, W = 2.0, 9.81, 10.0, 1.0
RHO_G_L0 = RHO * G * L0                          # 196.2 kPa
H = 0.5
D_TRANSLATE = 2.5
# The floor (amendment D19a, 232k-C round 2, before any soil mode-F number). "block": a bonded
# elastic floor block y in [-1.5, 0] (3 cells, the elastic leg's E, nu, rho) under the column;
# the translating gridvel prescribes the nodes 2 to 4 cells BELOW the soil, inside the block, so
# the prescribed set's node entries and exits never touch a soil particle's support (<= 1.5 h).
# "nodes": the frozen D19 floor -- the gridvel box itself under the column (y in [-1, 0]); its
# node exits put the soil's bottom row on GIMP tails and lock strain there (diagnosis record).
# "float" (mode F): the same block, no grid condition at all -- the block carries the whole
# column's weight as an upward body force (sum of external forces = 0), and the translation is a
# uniform body acceleration a(t) = -v phi'(t) on every particle, so the set moves rigidly without
# any boundary crossing the grid (diagnosis: a prescribed row that leaves the moving set releases
# the full column reaction at once, a kick of up to 13 % of W per cell on bspline2).
FLOOR = "block"                                  # mode A (mode F: "float")
FLOOR_TOP = 0.0                                  # "nodes" only: the box top above the base
BLOCK_T = 1.5                                    # "block": thickness
BLOCK_TAG = 99
FLOAT_DAMP_TRANSLATE = 0.0                       # "float": local damping during the translation
GRID = {"nodes": ((0.0, -3.5), (2, 28)), "block": ((0.0, -5.0), (2, 31)),
        "float": ((0.0, -5.0), (2, 31))}
ORIGIN, CELLS = GRID[FLOOR]

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


def cap_spec():
    """D23's cap from the cases block k0_ref_column_capped: (kind, tokens, y0, y1, rho)."""
    mats = reference_materials("k0_ref_column_capped")
    inside, body = False, None
    for raw in CASES.read_text(encoding="utf-8").splitlines():
        t = raw.split("#", 1)[0].split()
        if not t:
            continue
        if t[0] == "case":
            inside = (t[1] == "k0_ref_column_capped")
            continue
        if inside and t[0] == "end":
            break
        if inside and t[0] == "body" and t[1] == "cap":
            body = t
    # body cap y <y0> <y1> mat <tag> rho <rho>
    y0, y1 = float(body[3]), float(body[4])
    tag = int(body[body.index("mat") + 1])
    rho = float(body[body.index("rho") + 1])
    kind, toks = mats[tag]
    return kind, toks, y0, y1, rho


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


def build(ops, leg, mode, basis, translate=True, floor_top=None, floor=None, capped=False):
    tag, s0, _ = LEGS[leg]
    tl = TIMELINE[mode]
    mats = reference_materials()
    kind, toks = mats[tag]
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    floor = floor or ("float" if mode == "F" else FLOOR)
    origin, cells = GRID[floor]
    if capped:                                   # D23: one more cell row above the cap
        cells = (cells[0], cells[1] + 1)
    ops.nDMaterial(*material_command(tag, kind, toks)[1:])
    ops.mpm("grid", HOST, "-origin", *origin, "-h", H, "-cells", *cells, "-basis", basis)
    ops.mpm("mode", "quasiStatic" if mode == "A" else "explicit")
    ops.mpm("body", 1, "box", 0.0, 0.0, W, L0, "-ppc", 2, "-mat", tag, "-rho", RHO)
    if floor in ("block", "float"):
        ek, et = mats[1]                         # the elastic leg's operator
        ops.nDMaterial(*material_command(BLOCK_TAG, ek, et)[1:])
        ops.mpm("body", 2, "box", 0.0, -BLOCK_T, W, 0.0, "-ppc", 2, "-mat", BLOCK_TAG, "-rho", RHO)
    if capped:                                   # D23: the bonded elastic cap, body 3
        ck, ct, y0, y1, crho = cap_spec()
        ops.nDMaterial(*material_command(CAP_TAG, ck, ct)[1:])
        ops.mpm("body", 3, "box", 0.0, y0, W, y1, "-ppc", 2, "-mat", CAP_TAG, "-rho", crho)
    ops.mpm("fix", "xmin", "x")
    ops.mpm("fix", "xmax", "x")
    t0 = sum(n * dt for n, dt in (tl["gravity"], tl["hold1"]))
    n, dt = tl["translate"]
    t1 = t0 + n * dt
    v = -tl["v"] if translate else 0.0
    # the floor; it moves with the column (D17)
    if floor == "float":
        lo = hi = None                                # no grid condition (body forces, run())
    elif floor == "block":
        lo, hi = -BLOCK_T - H, -BLOCK_T + H          # 2 to 4 cells below the soil
    else:
        lo, hi = -1.0, FLOOR_TOP if floor_top is None else floor_top
    if floor != "float":
        ops.mpm("gridvel", "floor", "box", 0.0, lo, W, hi, "v", "free", v,
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


def run(ops, leg, mode, basis, translate=True, chunk=None, floor_top=None, history=False,
        floor=None, capped=False, probe=None):
    """One leg; returns the record (meter, crossings, inertia, floor hits, state)."""
    tag, _, s1 = LEGS[leg]
    tl = TIMELINE[mode]
    wall0 = time.time()
    floor = floor or ("float" if mode == "F" else FLOOR)
    origin = GRID[floor][0]
    build(ops, leg, mode, basis, translate, floor_top, floor, capped)
    x0 = _as_list(ops.mpm("particleData", "x"))
    body = _as_list(ops.mpm("particleData", "body"))
    N = sum(1 for b in body if b == 1.0)          # the soil (body 1, seeded first)
    assert all(b == 1.0 for b in body[:N])
    Y = [x0[2 * k + 1] for k in range(N)]
    cell = [math.floor((Y[k] - origin[1]) / H) for k in range(N)]
    crossings = [0] * N
    floor_hits = 0
    floor_stage0 = [0]
    pfloor = P_FLOOR.get(leg)
    W_all = float(_as_list(ops.mpm("info", "summary"))[0]) * G     # every particle's weight
    inertia, inertia_max, mom_prev, t_prev = None, 0.0, None, None
    hist = []
    chunk = chunk or (1 if mode == "A" else 100)

    def sample(phase):
        nonlocal floor_hits, mom_prev, t_prev, inertia, inertia_max
        if probe is not None:                    # diagnosis hook: probe(ops, phase)
            probe(ops, phase)
        x = _as_list(ops.mpm("particleData", "x"))
        for k in range(N):
            c = math.floor((x[2 * k + 1] - origin[1]) / H)
            if c != cell[k]:
                crossings[k] += abs(c - cell[k])
                cell[k] = c
        # the floors act only after the flip: stage 0 of both soil models is pressure
        # independent (ManzariDafalias.h:546-550; LadrunoHS.cpp stage-0 branch); counted from
        # the flip on, and the gravity stage's count is reported separately
        if pfloor is not None:
            nh = _floor_hits(_as_list(ops.mpm("particleData", "sig"))[:3 * N],
                             _as_list(ops.mpm("particleData", "sigZZ"))[:N], pfloor)
            if phase == "gravity":
                floor_stage0[0] += nh
            else:
                floor_hits += nh
        if mode == "F":
            s = _as_list(ops.mpm("info", "summary"))
            t = float(ops.mpm("info", "time"))
            if mom_prev is not None and t > t_prev:
                inertia = abs((s[2] - mom_prev) / (t - t_prev)) / W_all
                if phase in ("translate", "hold2"):
                    inertia_max = max(inertia_max, inertia)
            mom_prev, t_prev = s[2], t
            if history:
                hist.append((t, inertia, s[-1], phase))

    # "float": per-particle body force (per unit mass, not scaled by the series): gravity on
    # every particle, the block's upward support (the whole weight over the block's mass), and
    # the translation's uniform acceleration acc. Gravity is applied at once (stage 0 is elastic
    # in every leg) and settles under the damping.
    m_all = _as_list(ops.mpm("particleData", "m"))
    M_tot = sum(m_all)
    M_blk = sum(mm for mm, b in zip(m_all, body) if b == 2.0)

    def set_body_force(acc):
        if floor != "float":
            return
        up = G * M_tot / M_blk if M_blk > 0 else 0.0
        bf = []
        for b in body:
            bf += [0.0, -G + acc + (up if b == 2.0 else 0.0)]
        ops.mpm("particleBodyForce", *bf)

    v_tr = -tl["v"] if translate else 0.0
    Tr = tl["Tr"]

    def advance(n, dt, phase):
        done = 0
        if floor == "float" and phase == "translate":
            # the trapezoid's acceleration: -v/Tr over [t0, t0+Tr], 0, +v/Tr over [t1-Tr, t1]
            nr = int(round(Tr / dt)) if Tr > 0 else 0
            pieces = [(nr, v_tr / Tr if Tr > 0 else 0.0), (n - 2 * nr, 0.0),
                      (nr, -v_tr / Tr if Tr > 0 else 0.0)]
            if Tr == 0.0:
                raise ValueError("float needs Tr > 0 (a velocity jump has no finite acceleration)")
            # local damping is not Galilean invariant (interface §6.5): it would brake the
            # rigid acceleration itself, so the translation runs undamped (FLOAT_DAMP_TRANSLATE)
            ops.mpm("damping", "local", FLOAT_DAMP_TRANSLATE)
            for m_, acc in pieces:
                set_body_force(acc)
                adv(m_, dt, phase)
            set_body_force(0.0)
            ops.mpm("damping", "local", 0.8)
            return
        set_body_force(0.0)
        adv(n, dt, phase)

    def adv(n, dt, phase):
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
    rows, ref_build = reference_rows(leg, REF_JSON_CAPPED if capped else REF_JSON)
    sq = 0.0
    for k in range(N):
        r = _ref_at(rows, Y[k])
        d = [sig[3 * k] - r["sxx"], sig[3 * k + 1] - r["syy"], sig[3 * k + 2] - r["sxy"],
             szz[k] - r["szz"]]
        sq += sum(e * e for e in d)
    meter = math.sqrt(sq / N) / RHO_G_L0
    drop = sum(x[2 * k + 1] - Y[k] for k in range(N)) / N
    rec = dict(leg=leg, mode=mode, basis=basis, translate=translate, particles=N, floor=floor,
               capped=capped,
               meter=meter, mean_displacement_y=drop,
               crossings=dict(min=min(crossings), max=max(crossings),
                              mean=sum(crossings) / N, per_particle=crossings),
               floor_hits=floor_hits, floor_hits_stage0=floor_stage0[0], p_floor=pfloor,
               state_before=state_before, state_after=_state(ops, leg),
               ref_build=ref_build, wall_s=time.time() - wall0)
    rec["y_error"] = [(Y[k], (sig[3 * k + 1] - _ref_at(rows, Y[k])["syy"]) / RHO_G_L0,
                       (sig[3 * k] - _ref_at(rows, Y[k])["sxx"]) / RHO_G_L0) for k in range(N)]
    if history:
        rec["history"] = hist
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


def k0(ops, mode, basis, legs=None, capped=False):
    legs = legs or list(LEGS)
    out = dict(case="K0", mode=mode, basis=basis, build=str(ops.ladrunoBuild()),
               rho_g_l0=RHO_G_L0, capped=capped, legs={})
    for leg in legs:
        rt = run(ops, leg, mode, basis, translate=True, capped=capped)
        rz = run(ops, leg, mode, basis, translate=False, capped=capped)
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
    ap.add_argument("--capped", action="store_true", help="D23: the bonded 20 kPa elastic cap")
    a = ap.parse_args(argv)
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from _testbed import ops   # noqa: E402
    res = k0(ops, a.mode, a.basis, a.mat, capped=a.capped)
    suffix = ("_capped" if a.capped else "") + ("" if not a.mat else "_" + "_".join(a.mat))
    path = Path(a.out) if a.out else OUT_DIR / f"k0_meter_{a.mode}_{a.basis}{suffix}.json"
    path.write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(json.dumps({k: v for k, v in res.items() if k != "legs"}, indent=1))
    for leg, r in res["legs"].items():
        print(leg, "meter", r["translated"]["meter"], "twin", r["twin"]["meter"],
              "crossings", {k: v for k, v in r["translated"]["crossings"].items() if k != "per_particle"},
              "drop", r["translated"]["mean_displacement_y"], "floor_hits",
              r["translated"]["floor_hits"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
