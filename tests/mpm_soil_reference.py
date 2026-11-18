"""WP-232k (232k-V): the FE references behind tests/mpm_soil_cases.txt -- never MPM.

Never reads the MPM soil C++ (SRC/domain/pattern/ladrunoMPM/). Every number comes from the fork's
own FE (LadrunoBrick slabs and one-element drivers) or a closed form. Units kN, t, m, s, kPa
(interface D19: rho = 2.0 t/m^3, g = 9.81, rho g l0 = 196.2 kPa).

    python3.12 tests/mpm_soil_reference.py k0            # K0 FE columns, 4 materials (+ capped, D23)
    python3.12 tests/mpm_soil_reference.py onepoint      # S1 one-point drivers + Delta F paths
    python3.12 tests/mpm_soil_reference.py s2            # S2 column + surcharge
    python3.12 tests/mpm_soil_reference.py s4            # S4 low-g one-point reference
    python3.12 tests/mpm_soil_reference.py footing <hs|sanisand> <n_across_b> [variant]

FE legs need a built fork binary: LADRUNO_BIN=<dir with opensees.pyd>, else this worktree's
dist/bin. Each leg writes tests/data/mpm/soil/<leg>.json with the build hash (ops.ladrunoBuild()),
the deck, the date and the numbers; tests/test_mpm_soil_reference.py checks those files (no
binary needed).

STAGE ROUTES, recorded in every JSON (`flip_route`):
* SANISAND and HS flip with `updateMaterialStage -material 1 -stage 1`, as their decks do
  (tests/test_ladruno_sanisand.py:614/:620). `MaterialStageParameter::setDomain` stops at the first
  answering element, so that command reaches ONE element's Gauss points; the others flip lazily
  (SANISAND: ladrunoRunStageFlipOnce, without Elastic2Plastic; HS: syncStageFromRegistry, the full
  transition). For SANISAND that difference is inert where every point has p > Pmin and
  eta < M_c at the flip; each record carries min p and max eta to show it.
* DruckerPrager is staged by `parameter ... element e materialState` over every element (all 8
  points each): materialState 0 before gravity, 1 at the flip.
"""
from __future__ import annotations

import datetime as _dt
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
DATA = HERE / "data" / "mpm" / "soil"

# ----------------------------------------------------------------------------------------------
# Units and the decks (quoted; see tests/mpm_soil_cases.txt for the same lines)
# ----------------------------------------------------------------------------------------------
G_GRAV = 9.81
RHO = 2.0

# SANISAND: Gorini's set, tests/test_ladruno_sanisand.py:245-263 (kPa), scheme 1 defaults
# (`_build_confined`, :590, passes opts=()): -maxSubsteps 0 (uncapped), honorTolR 0 (TolE 1e-4).
SAN_PARAMS = [264.32, 0.3129, 0.6944, 1.33090, 0.71, 0.027, 0.83, 0.45, 101.0, 0.005, 1.3, 0.968,
              3.5, 0.05, 5.75, 12.5, 1100.0, RHO]
SAN_PMIN = 1.0e-3 * SAN_PARAMS[8]          # the class default, LadrunoSANISAND.cpp:245
SAN_MC = SAN_PARAMS[3]


def san_g0() -> float:
    """Stage-0 shear modulus: pressure-independent, G0 P_atm (2.97 - e)^2 / (1 + e)
    (ManzariDafalias.cpp GetElasticModuli, mElastFlag == 0 branch)."""
    G0, e, Pa = SAN_PARAMS[0], SAN_PARAMS[2], SAN_PARAMS[8]
    return G0 * Pa * (2.97 - e) ** 2 / (1.0 + e)


# HS: the strip-footing deck of tests/test_ladruno_hs_material.py:704-719, without -sig0 (K3/D6).
HS_ARGS = [35000.0, 35000.0, 105000.0, 0.5, 10.0, 38.0, 8.0,
           "-nu", 0.2, "-pref", 100.0, "-Rf", 0.9, "-OCR", 1.5]
HS_K0NC = 1.0 - math.sin(math.radians(38.0))   # the default K0nc = 1 - sin(phi)

# ElasticIsotropic: SANISAND's stage-0 operator (G, nu), so the elastic leg is SANISAND's stage 0.
EL_NU = SAN_PARAMS[1]
EL_E = 2.0 * san_g0() * (1.0 + EL_NU)

# DruckerPrager (UW): the program's admitted deck, tests/mpm_implicit_cases.txt:286-301 (the
# 232f/232j mild cone): K 1.5e5, G 15517.24, sigma_y 0.2, rho_dp 0.210128, rho_bar 0 (non-assoc.).
DP_ARGS = [1.5e5, 15517.241379310342, 0.2, 0.21012821217029673, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
           0.0, 0.0]
DP_NU = 0.45

MATERIALS = ("elastic", "dp", "sanisand", "hs")

# The K0 column, interface D19: x in [0, 1.0], y in [0, 10], h 0.5, ppc 2 (160 particles).
COL = dict(W=1.0, L0=10.0, h=0.5, ppc=2)
# D23: a bonded elastic cap one cell thick whose weight is 20 kPa.
CAP = dict(q=20.0, t=0.5)


def cap_rho() -> float:
    return CAP["q"] / (G_GRAV * CAP["t"])


def column_particle_y(L0=COL["L0"], h=COL["h"], ppc=COL["ppc"]):
    """Particle initial Y of the lattice rule (cell-centred sub-cells): (k + 1/2) h / ppc."""
    n = int(round(L0 / (h / ppc)))
    return [(k + 0.5) * h / ppc for k in range(n)]


# ----------------------------------------------------------------------------------------------
# The fork FE engine
# ----------------------------------------------------------------------------------------------
_OPS = None


def fe_engine():
    global _OPS
    if _OPS is not None:
        return _OPS
    b = os.environ.get("LADRUNO_BIN") or str(REPO / "dist" / "bin")
    if not (Path(b) / "opensees.pyd").exists() and not list(Path(b).glob("opensees*.so")):
        raise SystemExit(f"no fork binary at {b} (set LADRUNO_BIN); FE legs need one")
    os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")
    os.environ["PATH"] = b + os.pathsep + os.environ.get("PATH", "")
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(b)
    sys.path.insert(0, b)
    import opensees as ops  # noqa: E402
    if os.path.normcase(os.path.dirname(os.path.abspath(ops.__file__))) != os.path.normcase(b):
        raise SystemExit(f"WRONG ENGINE: imported {ops.__file__}, expected {b}")
    _OPS = ops
    return ops


def _provenance(leg: str, deck: str) -> dict:
    ops = fe_engine()
    return dict(leg=leg, build=ops.ladrunoBuild(), engine=os.path.abspath(ops.__file__),
                deck=f"tests/mpm_soil_reference.py::{deck}",
                date=_dt.date.today().isoformat(),
                units="kN, t, m, s, kPa; stress tension-positive (OpenSees sign)")


def _save(name: str, payload: dict) -> Path:
    DATA.mkdir(parents=True, exist_ok=True)
    p = DATA / f"{name}.json"
    payload["output"] = str(p.relative_to(REPO)).replace("\\", "/")
    p.write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    return p


def decks() -> dict:
    return dict(
        elastic=dict(cmd="nDMaterial ElasticIsotropic 1 E nu", E=EL_E, nu=EL_NU,
                     note="SANISAND's stage-0 operator: G = G0 P_atm (2.97-e)^2/(1+e)"),
        dp=dict(cmd="nDMaterial DruckerPrager 1 " + " ".join(repr(x) for x in DP_ARGS),
                note="tests/mpm_implicit_cases.txt:286-301 (232f/232j mild cone); staged "
                     "materialState 0 first"),
        sanisand=dict(cmd="nDMaterial LadrunoSANISAND 1 " + " ".join(repr(x) for x in SAN_PARAMS),
                      note="tests/test_ladruno_sanisand.py:245-263, scheme 1 defaults "
                           "(-maxSubsteps 0, honorTolR 0)"),
        hs=dict(cmd="nDMaterial LadrunoHS 1 " + " ".join(str(x) for x in HS_ARGS) + " -stage 0",
                note="tests/test_ladruno_hs_material.py:704-719 without -sig0 (K3/D6)"),
        units=dict(rho=RHO, g=G_GRAV))


def _material(ops, mat: str, tag: int = 1, san_opts=()) -> None:
    if mat == "elastic":
        ops.nDMaterial("ElasticIsotropic", tag, EL_E, EL_NU)
    elif mat == "dp":
        ops.nDMaterial("DruckerPrager", tag, *DP_ARGS)
    elif mat == "sanisand":
        ops.nDMaterial("LadrunoSANISAND", tag, *SAN_PARAMS, *san_opts)
    elif mat == "hs":
        ops.nDMaterial("LadrunoHS", tag, *HS_ARGS, "-stage", 0)
    else:
        raise ValueError(mat)


def _stage(ops, mat: str, stage: int, elems) -> str:
    """Stage every soil point the way its deck does; return the route recorded in the JSON."""
    if mat == "dp":
        ops.parameter(900 + stage)
        for e in elems:
            ops.addToParameter(900 + stage, "element", e, "materialState")
        ops.updateParameter(900 + stage, float(stage))
        return "parameter element e materialState over every element (all points)"
    if mat in ("sanisand", "hs"):
        ops.updateMaterialStage("-material", 1, "-stage", stage)
        return ("updateMaterialStage -material 1: reaches ONE element's points; the rest flip "
                "lazily (SANISAND without Elastic2Plastic; HS full transition)")
    return "none (elastic)"


# ----------------------------------------------------------------------------------------------
# K0: the FE column (no translation) on the D19 domain
# ----------------------------------------------------------------------------------------------
_GP = 1.0 / math.sqrt(3.0)


def _column_model(ops, mat, g=G_GRAV, capped=False, L0=COL["L0"], h=COL["h"], W=COL["W"]):
    """Plane-strain slab of LadrunoBrick B-bar, nx = W/h, ny = L0/h (+1 cap row), thickness h,
    z fixed everywhere; sides x = 0 and x = W rollers (fix x), base fixed. Soil is material 1,
    the cap (D23) ElasticIsotropic tag 2 with rho_cap g t = 20 kPa."""
    nx, ny = int(round(W / h)), int(round(L0 / h))
    nyt = ny + (1 if capped else 0)
    t = h

    def nid(i, j, k):
        return 1 + i + (nx + 1) * (j + (nyt + 1) * k)

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for k in (0, 1):
        for j in range(nyt + 1):
            for i in range(nx + 1):
                n = nid(i, j, k)
                ops.node(n, i * h, j * h, k * t)
                if j == 0:
                    ops.fix(n, 1, 1, 1)
                else:
                    ops.fix(n, 1 if i in (0, nx) else 0, 0, 1)
    _material(ops, mat)
    if capped:
        ops.nDMaterial("ElasticIsotropic", 2, EL_E, EL_NU)
    soil, cap, e = [], [], 0
    for j in range(nyt):
        for i in range(nx):
            e += 1
            isc = capped and j == ny
            rho = cap_rho() if isc else RHO
            ops.element("LadrunoBrick", e, nid(i, j, 0), nid(i + 1, j, 0), nid(i + 1, j + 1, 0),
                        nid(i, j + 1, 0), nid(i, j, 1), nid(i + 1, j, 1), nid(i + 1, j + 1, 1),
                        nid(i, j + 1, 1), 2 if isc else 1, "-formulation", "bbar",
                        "-b", 0.0, -rho * g, 0.0)
            (cap if isc else soil).append((e, j))
    return dict(nx=nx, ny=ny, nyt=nyt, soil=soil, cap=cap, nid=nid, ne=e)


def _analysis(ops, tol=1e-8, iters=40):
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormUnbalance", tol, iters, 0)
    ops.algorithm("Newton")


def _row_table(ops, m, h=COL["h"], which="soil"):
    """Per element row: (Y_centre, sxx, syy, szz, sxy), the mean over the row's Gauss points.

    In this 1D column the B-bar hex has u_y linear in y per element, so its stress is constant
    per element (GP spread recorded) and equals the exact field at the element centre only; an
    interpolation of FE sigma_yy between GPs would be off by up to rho g h / 2. Hence `_ratios`."""
    byrow, spread, rows = {}, 0.0, []
    for e, j in m[which]:
        byrow.setdefault(j, []).extend(
            list(ops.eleResponse(e, "material", gp + 1, "stress"))[:4] for gp in range(8))
    for j in sorted(byrow):
        a = np.array(byrow[j], dtype=float)
        spread = max(spread, float(np.max(np.ptp(a, axis=0))))
        rows.append([(j + 0.5) * h, *[float(v) for v in a.mean(axis=0)]])
    return rows, spread


def _ratios(table, ys, syy_exact):
    """sigma_ref at each Y: syy from equilibrium (closed form); sxx, szz, sxy from the FE stress
    RATIOS to sigma_yy (the constitutive answer, which is what the FE contributes), linear in Y
    between element centres, linear extrapolation from the end pair."""
    T = np.array(table)
    Yc = T[:, 0]
    rx, rz, rs = T[:, 1] / T[:, 2], T[:, 3] / T[:, 2], T[:, 4] / T[:, 2]
    out = []
    for y in ys:
        i = int(np.clip(np.searchsorted(Yc, y) - 1, 0, len(Yc) - 2))
        f = (y - Yc[i]) / (Yc[i + 1] - Yc[i])

        def lin(a):
            return float(a[i] + f * (a[i + 1] - a[i]))
        sy = syy_exact(y)
        out.append(dict(Y=float(y), sxx=lin(rx) * sy, syy=sy, szz=lin(rz) * sy, sxy=lin(rs) * sy,
                        ratio_xx=lin(rx), ratio_zz=lin(rz)))
    return out


def _flip_inertness(table):
    """SANISAND's Elastic2Plastic is a no-op where p > Pmin (tr sigma >= 3 Pmin) and eta <= M_c
    (ManzariDafalias.cpp:5787-5810). p, q from the 3D stress (compression positive)."""
    pmin, etamax = math.inf, 0.0
    for _, sxx, syy, szz, sxy in table:
        sc = np.array([[-sxx, -sxy, 0.0], [-sxy, -syy, 0.0], [0.0, 0.0, -szz]])
        p = np.trace(sc) / 3.0
        s = sc - p * np.eye(3)
        q = math.sqrt(1.5 * float(np.sum(s * s)))
        pmin = min(pmin, p)
        etamax = max(etamax, q / p if p > 0 else math.inf)
    return dict(p_min_at_flip=pmin, eta_max_at_flip=etamax, Pmin=SAN_PMIN, Mc=SAN_MC,
                elastic2plastic_inert=bool(pmin > SAN_PMIN and etamax <= SAN_MC))


def fe_column(mat, capped=False, g=G_GRAV, n_grav=10, n_hold=2, surcharge=None, n_sur=10):
    """K0 / S2 column: stage 0 gravity ramp (n_grav LoadControl steps) -> flip -> n_hold hold steps
    (LoadControl 0) [-> surcharge on the top nodes, tributary, in n_sur steps]. Reports, at each
    checkpoint, the GP-station table and the values at the MPM particles' initial Y (soil only),
    with sigma_yy also in closed form."""
    ops = fe_engine()
    t0 = time.time()
    m = _column_model(ops, mat, g=g, capped=capped)
    soil_elems = [e for e, _ in m["soil"]]
    route0 = _stage(ops, mat, 0, soil_elems)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", *range(1, m["ne"] + 1), "-type", "-brickW")
    _analysis(ops)
    ops.integrator("LoadControl", 1.0 / n_grav)
    ops.analysis("Static")
    for st in range(n_grav):
        if ops.analyze(1) != 0:
            raise RuntimeError(f"{mat} column gravity step {st + 1} failed")
    ys = column_particle_y()
    q_top = cap_rho() * g * CAP["t"] if capped else 0.0

    def syy_exact(y, extra=0.0):
        return -(RHO * g * (COL["L0"] - y) + q_top + extra)

    def checkpoint(extra=0.0):
        tab, spread = _row_table(ops, m)
        cp = dict(rows=tab, gp_spread=spread,
                  centre_syy_error=max(abs(r[2] - syy_exact(r[0], extra)) for r in tab),
                  particles=_ratios(tab, ys, lambda y: syy_exact(y, extra)))
        if capped:
            cp["cap_rows"] = _row_table(ops, m, which="cap")[0]
        return cp

    out = dict(material=mat, capped=capped, g=g, rho=RHO, L0=COL["L0"], W=COL["W"], h=COL["h"],
               hexes=m["ne"], element="LadrunoBrick -formulation bbar -geom linear, plane-strain "
               "slab (z fixed), sides fix x, base fixed", rho_g_l0=RHO * g * COL["L0"],
               q_top=q_top, stage0_route=route0, n_grav=n_grav, n_hold=n_hold)
    if capped:
        out["cap"] = dict(material="ElasticIsotropic", E=EL_E, nu=EL_NU, rho=cap_rho(),
                          thickness=CAP["t"], weight_kPa=CAP["q"])
    out["end_stage0"] = checkpoint()
    ops.loadConst("-time", 0.0)
    out["flip_route"] = _stage(ops, mat, 1, soil_elems)
    if mat == "sanisand":
        out["flip_check"] = _flip_inertness(out["end_stage0"]["rows"])
    out["after_flip"] = checkpoint()
    ops.integrator("LoadControl", 0.0)
    for st in range(n_hold):
        if ops.analyze(1) != 0:
            raise RuntimeError(f"{mat} column hold step {st + 1} failed")
    out["after_hold"] = checkpoint()
    if surcharge is not None:
        top = [m["nid"](i, m["nyt"], k) for k in (0, 1) for i in range(m["nx"] + 1)]
        h = COL["h"]
        ops.timeSeries("Linear", 2)
        ops.pattern("Plain", 2, 2)
        for n in top:
            i = (n - 1) % (m["nx"] + 1)
            ops.load(n, 0.0, -surcharge * (h * h / 4.0) * (1 if i in (0, m["nx"]) else 2), 0.0)
        ops.integrator("LoadControl", 1.0 / n_sur)
        for st in range(n_sur):
            if ops.analyze(1) != 0:
                raise RuntimeError(f"{mat} column surcharge step {st + 1} failed")
        out["surcharge"] = surcharge
        out["after_surcharge"] = checkpoint(extra=surcharge)
    out["wall_s"] = round(time.time() - t0, 2)
    r = out["after_hold"]["particles"]
    print(f"COL {mat} capped={capped} g={g}: ratio_xx {r[0]['ratio_xx']:.6f} .. "
          f"{r[-1]['ratio_xx']:.6f}; centre syy err {out['after_hold']['centre_syy_error']:.2e}"
          f" ({out['wall_s']} s)", flush=True)
    return out


def leg_k0() -> dict:
    res = {"decks": decks()}
    for capped in (False, True):
        legs = {mat: fe_column(mat, capped=capped) for mat in MATERIALS}
        name = "k0_column_capped" if capped else "k0_column"
        payload = dict(_provenance(name, "fe_column"), **res, legs=legs,
                       closed_form=dict(
                           K0_elastic=EL_NU / (1 - EL_NU), K0_sanisand_stage0=EL_NU / (1 - EL_NU),
                           K0_hs_stage0=HS_K0NC, K0_dp_stage0=DP_NU / (1 - DP_NU),
                           note="laterally confined (oedometric) isotropic elastic stage 0 with a "
                                "constant nu: sigma_xx = sigma_zz = nu/(1-nu) sigma_yy; HS stage 0 "
                                "uses nu_s0 = max(nu_ur, K0nc/(1+K0nc)), so the ratio is K0nc "
                                "(LadrunoHS.cpp integrate3D stage 0)"),
                       particle_y=column_particle_y(),
                       meter_reference=(
                           "sigma_ref(Y) = (sxx, syy, szz) of `particles` at the particle's initial"
                           " Y, checkpoint after_hold (the stress the translation must keep): syy "
                           "closed form, sxx and szz = the FE ratio to syy times syy"))
        _save(name, payload)
    return res


# ----------------------------------------------------------------------------------------------
def main(argv=None) -> None:
    argv = list(sys.argv[1:] if argv is None else argv)
    what = argv[0] if argv else "k0"
    if what == "k0":
        leg_k0()
    else:
        raise SystemExit(f"unknown leg {what!r}")


if __name__ == "__main__":
    main()
