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

GRAVITY: every body force goes through `eleLoad -type -BrickW` (constant BRICKW). Issue #1111:
the Python interpreter matches that type case-sensitively and silently ignores `-brickW`, after
which LadrunoBrick applies `-b` unconditionally. The K0, S2 and S3 HS records were re-run with
the real ramp on 2026-11-18 (fa68005ac); every stage-0 law here is linear, and they agree with
the earlier records to 1e-11 kPa.

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

# The body-force eleLoad type, spelled as the Python interpreter matches it, case-sensitively
# (issue #1111): OpenSeesPatternCommands.cpp:542 `else if (strcmp(type,"-BrickW") == 0)`, and an
# unknown type falls through to `return 0` (:1072-1074), silently. Measured 2026-11-18
# (b94c9d75d): `-type -brickW` attaches NOTHING, and LadrunoBrick then applies its `-b`
# unconditionally (LadrunoBrick.cpp:1516-1519, `if (applyLoad == 0) bodyForce -= ... b[p]`), so
# a "ramped" gravity is in fact full gravity from the first step. Records before this fix (K0,
# S2, S3 HS at 42c11333a) were taken that way; every stage-0 law here is linear elastic, so their
# end-of-stage-0 states are the same, and they were re-run with the ramp.
BRICKW = "-BrickW"

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

# DruckerPrager (UW). GATED leg `dp` (orchestrator's ruling, 2026-11-17): SANISAND's stage-0 K
# and G (so A = m(SANISAND)/m(DP) shares the elastic operator) with the mild cone's strength and
# flow: sigma_y 0.2, rho_dp 0.210128, rho_bar 0 (non-associated), the rest 0.
DP_K = EL_E / (3.0 * (1.0 - 2.0 * EL_NU))
DP_G = EL_E / (2.0 * (1.0 + EL_NU))
DP_ARGS = [DP_K, DP_G, 0.2, 0.21012821217029673, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
DP_NU = EL_NU
# REPORTED leg `dp_mildcone`: the program's admitted deck, tests/mpm_implicit_cases.txt:286-301
# (232f/232j): K 1.5e5, G 15517.24 (nu 0.45), the same cone.
DPM_ARGS = [1.5e5, 15517.241379310342, 0.2, 0.21012821217029673, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0]
DPM_NU = 0.45

MATERIALS = ("elastic", "dp", "sanisand", "hs", "dp_mildcone")

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
                note="GATED: SANISAND's stage-0 K, G with the mild cone's sigma_y, rho, rho_bar; "
                     "staged materialState 0 first"),
        dp_mildcone=dict(cmd="nDMaterial DruckerPrager 1 " + " ".join(repr(x) for x in DPM_ARGS),
                         note="REPORTED: tests/mpm_implicit_cases.txt:286-301 (232f/232j mild "
                              "cone); staged materialState 0 first"),
        sanisand=dict(cmd="nDMaterial LadrunoSANISAND 1 " + " ".join(repr(x) for x in SAN_PARAMS),
                      note="tests/test_ladruno_sanisand.py:245-263, scheme 1 defaults "
                           "(-maxSubsteps 0, honorTolR 0)"),
        hs=dict(cmd="nDMaterial LadrunoHS 1 " + " ".join(str(x) for x in HS_ARGS) + " -stage 0",
                note="tests/test_ladruno_hs_material.py:704-719 without -sig0 (K3/D6)"),
        units=dict(rho=RHO, g=G_GRAV))


def _material(ops, mat: str, tag: int = 1, san_opts=()) -> None:
    """`san_opts`: extra tokens after the deck -- SANISAND's positional/flag options, or HS's
    flags (placed before HS's `-stage 0`)."""
    if mat == "elastic":
        ops.nDMaterial("ElasticIsotropic", tag, EL_E, EL_NU)
    elif mat == "dp":
        ops.nDMaterial("DruckerPrager", tag, *DP_ARGS)
    elif mat == "dp_mildcone":
        ops.nDMaterial("DruckerPrager", tag, *DPM_ARGS)
    elif mat == "sanisand":
        ops.nDMaterial("LadrunoSANISAND", tag, *SAN_PARAMS, *san_opts)
    elif mat == "hs":
        ops.nDMaterial("LadrunoHS", tag, *HS_ARGS, *san_opts, "-stage", 0)
    else:
        raise ValueError(mat)


def _stage(ops, mat: str, stage: int, elems) -> str:
    """Stage every soil point the way its deck does; return the route recorded in the JSON."""
    if mat in ("dp", "dp_mildcone"):
        # UW DruckerPrager's materialState: 0 elastic; 1 STILL ELASTIC, with K and G rescaled once
        # by sqrt(1 + sigma_m / P_atm) (updateElasticParam, DruckerPrager.cpp:823-838, and the
        # `|| mElastFlag < 2` elastic exit at :446); 2 elastoplastic (the constructor default,
        # :160). So DP's plastic flip is materialState 2, not 1.
        val = 0 if stage == 0 else 2
        ops.parameter(900 + stage)
        for e in elems:
            ops.addToParameter(900 + stage, "element", e, "materialState")
        ops.updateParameter(900 + stage, float(val))
        return (f"parameter element e materialState over every element (all points), value {val}"
                " (DP: 0 elastic, 2 elastoplastic; 1 is elastic with rescaled moduli)")
    if mat in ("sanisand", "hs"):
        ops.updateMaterialStage("-material", 1, "-stage", stage)
        return ("updateMaterialStage -material 1: reaches ONE element's points; the rest flip "
                "lazily (SANISAND without Elastic2Plastic; HS full transition)")
    return "none (elastic)"


# ----------------------------------------------------------------------------------------------
# K0: the FE column (no translation) on the D19 domain
# ----------------------------------------------------------------------------------------------
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
    ops.eleLoad("-ele", *range(1, m["ne"] + 1), "-type", BRICKW)
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
                           K0_dp_mildcone_stage0=DPM_NU / (1 - DPM_NU),
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
# K0p: the plastic cell-crossing meter's reference (gate K0p, pre-registered 2026-11-18). K0's
# capped column, stage 0 -> flip -> hold, then a stage-1 pre-shear: nu -> 0.05 (SANISAND
# `poissonRatio`; DP `bulkModulus` + `shearModulus` for nu 0.05 at the same E), gravity +30 % in
# 10 steps (every body: soil and cap), then a 2-step hold.
# ----------------------------------------------------------------------------------------------
K0P = dict(nu=0.05, dg=0.30, n_pre=10, n_hold=2)


def _flip_every_point(ops, mat, elems) -> str:
    """Flip EVERY soil point (not updateMaterialStage's one element)."""
    if mat in ("dp", "dp_mildcone"):
        return _stage(ops, mat, 1, elems)
    if mat == "elastic":
        return "none (elastic)"
    ops.parameter(960)
    for e in elems:
        # SANISAND: ManzariDafalias::setParameter wants ('materialState', tag) -> id 5, which runs
        # Elastic2Plastic on s = 1 (ManzariDafalias.cpp:975-981); HS: 'materialState' -> id 13,
        # the full applyStageTransition (LadrunoHS.cpp:1189-1192, :1225-1232).
        ops.addToParameter(960, "element", e, "materialState", 1)
    ops.updateParameter(960, 1.0)
    return ("parameter element e materialState 1 over every soil element (all 8 points), "
            "updateParameter 1.0" + (" (id 5: Elastic2Plastic at every point)" if mat == "sanisand"
                                     else " (id 13: applyStageTransition at every point)"))


def _preshear_params(ops, mat, elems) -> str:
    nu = K0P["nu"]
    if mat == "sanisand":
        # ManzariDafalias.cpp:955-956 `else if (strcmp(argv[0],"poissonRatio") == 0) { ...
        # return param.addObject(7, this);` and :995-997 `else if (responseID == 7) { m_nu =
        # info.theDouble; }`; LadrunoSANISAND::setParameter falls through to the base (:5193).
        ops.parameter(961)
        for e in elems:
            ops.addToParameter(961, "element", e, "poissonRatio", 1)
        ops.updateParameter(961, nu)
        return "SANISAND poissonRatio (id 7) -> 0.05 at every soil point"
    if mat in ("dp", "dp_mildcone"):
        # DruckerPrager.cpp:1056-1061 `shearModulus` -> addObject(10), `bulkModulus` ->
        # addObject(11); :1097-1107 `mG = info.theDouble; mCe = mK*mIIvol + 2*mG*mIIdev;` and
        # `mK = info.theDouble; mCe = ...`. E is the leg's own: E = 9 K G / (3 K + G).
        K0_, G0_ = (DP_K, DP_G) if mat == "dp" else (DPM_ARGS[0], DPM_ARGS[1])
        E = 9.0 * K0_ * G0_ / (3.0 * K0_ + G0_)
        K, G = E / (3.0 * (1.0 - 2.0 * nu)), E / (2.0 * (1.0 + nu))
        for pid, name, val in ((962, "shearModulus", G), (963, "bulkModulus", K)):
            ops.parameter(pid)
            for e in elems:
                ops.addToParameter(pid, "element", e, name)
            ops.updateParameter(pid, val)
        return f"DP shearModulus (id 10) -> {G!r}, bulkModulus (id 11) -> {K!r} (nu 0.05 at E {E!r})"
    if mat == "hs":
        return ("none: LadrunoHS answers no Poisson parameter (setParameter names: E50, Eoed, Eur, "
                "c, phi, psi, G0, gamma07, OCR, pop, rho, stage; LadrunoHS.cpp:1172-1250); "
                "gravity +30 % only")
    return "none (elastic): gravity +30 % only"


def _row_alpha(ops, m, h=COL["h"]):
    """SANISAND alpha (state 6-11) per element row, mean over the row's points, and per point."""
    rows, pts = {}, []
    for e, j in m["soil"]:
        for gp in range(8):
            a = list(ops.eleResponse(e, "material", gp + 1, "state"))[6:12]
            rows.setdefault(j, []).append(a)
            pts.append(a)
    return [[(j + 0.5) * h, *np.mean(np.array(rows[j]), axis=0).tolist()] for j in sorted(rows)], pts


def _ladder_advance(ops, t_from, t_to, dt0, floor=5.0e-4):
    """LoadControl from t_from to t_to with nominal step dt0 on the K9 ladder, halving on a full
    ladder failure down to floor * dt0, regrowing x2 after a success. Returns the counts."""
    dt, nfail, nrelax, nhalv, steps = dt0, 0, 0, 0, 0
    while ops.getTime() < t_to - 1e-12:
        dt = min(dt, t_to - ops.getTime())
        ops.integrator("LoadControl", dt)
        ok = False
        for i, (alg, tl, it) in enumerate(K9):
            ops.test("NormUnbalance", tl, it, 0)
            ops.algorithm(alg)
            if ops.analyze(1) == 0:
                ok, nrelax = True, nrelax + (i == 2)
                break
            nfail += 1
        if not ok:
            dt *= 0.5
            nhalv += 1
            if dt < floor * dt0:
                raise RuntimeError(f"ladder floor at t = {ops.getTime()}")
            continue
        steps += 1
        dt = min(2.0 * dt, dt0)
    return dict(steps=steps, failed_attempts=nfail, relaxed=nrelax, halvings=nhalv,
                ladder=[list(x) for x in K9])


def fe_column_k0p(mat, g=G_GRAV):
    ops = fe_engine()
    t0 = time.time()
    m = _column_model(ops, mat, g=g, capped=True)
    soil = [e for e, _ in m["soil"]]
    route0 = _stage(ops, mat, 0, soil)
    # ONE gravity pattern on a Path series: 0 -> 1 over t in [0, 1] (stage 0), 1 -> 1.3 over
    # [1, 2] (the pre-shear), held to t = 3 (PathSeries returns 0 AT its last time point).
    ops.timeSeries("Path", 1, "-time", 0.0, 1.0, 2.0, 3.0,
                   "-values", 0.0, 1.0, 1.0 + K0P["dg"], 1.0 + K0P["dg"])
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", *range(1, m["ne"] + 1), "-type", BRICKW)
    _analysis(ops)
    ops.integrator("LoadControl", 0.1)
    ops.analysis("Static")
    for st in range(10):
        if ops.analyze(1) != 0:
            raise RuntimeError(f"K0p {mat} gravity step {st + 1} failed")
    ys = column_particle_y()
    q_top = cap_rho() * g * CAP["t"]

    def checkpoint(factor):
        tab, spread = _row_table(ops, m)
        sy = lambda y: -factor * (RHO * g * (COL["L0"] - y) + q_top)  # noqa: E731
        cp = dict(rows=tab, gp_spread=spread,
                  centre_syy_error=max(abs(r[2] - sy(r[0])) for r in tab),
                  particles=_ratios(tab, ys, sy), cap_rows=_row_table(ops, m, which="cap")[0])
        if mat == "sanisand":
            cp["alpha_rows"], cp["_alpha_pts"] = _row_alpha(ops, m)
        return cp

    out = dict(material=mat, capped=True, g=g, rho=RHO, preshear=K0P, stage0_route=route0)
    out["end_stage0"] = checkpoint(1.0)
    out["flip_route"] = _flip_every_point(ops, mat, soil)
    if mat == "sanisand":
        out["flip_check"] = _flip_inertness(out["end_stage0"]["rows"])
    ops.integrator("LoadControl", 0.0)
    for st in range(2):
        if ops.analyze(1) != 0:
            raise RuntimeError(f"K0p {mat} hold step {st + 1} failed")
    out["after_hold"] = checkpoint(1.0)
    out["preshear_params"] = _preshear_params(ops, mat, soil)
    # The pre-shear: n_pre steps of dt = 1/n_pre on the K9 ladder (KrylovNewton 1e-6/40 ->
    # NewtonLineSearch 1e-6/60 -> KrylovNewton 1e-5/80), halving to 5e-4 of the step. Newton
    # 1e-8 (stage 0) stalls at step 3 for SANISAND: its scheme-1 TolE 1e-4 sets a residual floor.
    lad = _ladder_advance(ops, 1.0, 2.0, 1.0 / K0P["n_pre"])
    out["after_preshear"] = checkpoint(1.0 + K0P["dg"])
    out["preshear_ladder"] = lad
    ops.integrator("LoadControl", 0.0)
    for st in range(K0P["n_hold"]):
        ok = -1
        for alg, tl, it in K9:
            ops.test("NormUnbalance", tl, it, 0)
            ops.algorithm(alg)
            if ops.analyze(1) == 0:
                ok = 0
                break
        if ok != 0:
            raise RuntimeError(f"K0p {mat} post-pre-shear hold step {st + 1} failed")
    out["after_hold2"] = checkpoint(1.0 + K0P["dg"])
    if mat == "sanisand":
        a0 = np.array(out["after_hold"].pop("_alpha_pts"))
        a1 = np.array(out["after_preshear"].pop("_alpha_pts"))
        a2 = np.array(out["after_hold2"].pop("_alpha_pts"))
        d = np.linalg.norm(a1 - a0, axis=1)
        out["alpha_moved"] = dict(
            tol=1e-10, points=int(a0.shape[0]), moved_by_preshear=int(np.sum(d > 1e-10)),
            fraction=float(np.mean(d > 1e-10)), max_norm_dalpha=float(d.max()),
            min_norm_dalpha=float(d.min()),
            hold2_max_norm_dalpha=float(np.linalg.norm(a2 - a1, axis=1).max()),
            note="FE Gauss points, ||alpha(after_preshear) - alpha(after_hold)||; the MPM validity "
                 "condition (>= 50 % of particles) is the harness's, on particles")
        out["end_stage0"].pop("_alpha_pts", None)
    out["wall_s"] = round(time.time() - t0, 2)
    r = out["after_hold2"]["particles"]
    print(f"K0P {mat}: ratio_xx after pre-shear+hold {r[0]['ratio_xx']:.6f} .. {r[-1]['ratio_xx']:.6f}"
          + (f"; alpha moved {out['alpha_moved']['fraction']:.3f}" if mat == "sanisand" else "")
          + f" ({out['wall_s']} s)", flush=True)
    return out


def leg_k0p() -> dict:
    legs = {mat: fe_column_k0p(mat) for mat in MATERIALS}
    payload = dict(_provenance("k0p_column_capped", "fe_column_k0p"), decks=decks(), legs=legs,
                   particle_y=column_particle_y(), preshear=K0P,
                   meter_reference=("sigma_ref(Y) = legs.<mat>.after_hold2.particles at the "
                                    "particle's initial Y: syy = -1.3 (rho g (l0 - Y) + 20) closed "
                                    "form (gravity +30 % acts on the cap too), sxx, szz = FE ratio "
                                    "x syy"))
    _save("k0p_column_capped", payload)
    return payload


# ----------------------------------------------------------------------------------------------
# S1: one-point drivers on LadrunoBrick -geom hypo -formulation std
# ----------------------------------------------------------------------------------------------
# Every path is a sequence of deformation gradients F_0 = I, F_1, ..., F_N (row-major 3x3; plane
# strain has F[2][2] = 1 and no z coupling). The FE driver is one unit-cube LadrunoBrick with
# every node prescribed u_a = (F_n - I) X_a (sp + Path series, zero free equations), so its
# hypoIncrement sees exactly F_n -> F_n+1 at every Gauss point. The MPM harness imposes
# dF_n = F_n+1 F_n^-1 (L_n = (dF_n - I) / dt). Stage 0 runs the confinement leg; the flip (one
# element: updateMaterialStage reaches all 8 of its points, Elastic2Plastic included) comes after
# it, before the path's first step.
S1 = dict(p0=100.0, n_conf=10, eps_ax=5.0e-3, n_ax=40, lat_triax=0.25, gamma_a=2.0e-3,
          n_quarter=10, cycles=3)


def _hs_k_stage0() -> float:
    """HS stage-0 bulk modulus (LadrunoHS.cpp integrate3D stage 0): nu_s0 = max(nu_ur,
    min(0.49, K0nc/(1+K0nc))), Gur_ref = Eur_ref / (2 (1 + nu_ur)), K = Gur 2(1+nu_s0)/(3(1-2nu_s0))."""
    nu_ur, eur = 0.2, HS_ARGS[2]
    nu_s0 = max(nu_ur, min(0.49, HS_K0NC / (1.0 + HS_K0NC)))
    gur = eur / (2.0 * (1.0 + nu_ur))
    return gur * 2.0 * (1.0 + nu_s0) / (3.0 * (1.0 - 2.0 * nu_s0))


def _k_stage0(mat: str) -> float:
    # SANISAND stage 0: K = 2(1+nu)/(3(1-2nu)) G = E/(3(1-2nu)) with E = 2G(1+nu), i.e. DP_K
    return {"sanisand": DP_K, "hs": _hs_k_stage0()}[mat]


def s1_path(mat: str, dim: str, path: str) -> dict:
    """The F sequence: isotropic confinement to p0 (stage 0, n_conf steps), then the path."""
    nd = 3 if dim == "3d" else 2
    a = S1["p0"] / (nd * _k_stage0(mat))          # per-direction compressive strain, stage 0
    Fs = [np.eye(3)]
    for k in range(1, S1["n_conf"] + 1):
        F = np.eye(3)
        for i in range(nd):
            F[i, i] = 1.0 - a * k / S1["n_conf"]
        Fs.append(F)
    Fc = Fs[-1].copy()
    steps = []
    if path in ("triax", "oedo"):
        lat = S1["lat_triax"] if path == "triax" else 0.0
        for k in range(1, S1["n_ax"] + 1):
            e = S1["eps_ax"] * k / S1["n_ax"]
            D = np.eye(3)
            D[1, 1] = 1.0 - e
            D[0, 0] = 1.0 + lat * e
            if nd == 3:
                D[2, 2] = 1.0 + lat * e
            steps.append(D @ Fc)
    elif path == "tension":
        # lateral (x) stretch only, 40 steps to F_xx = 1 + 5e-3: drives sigma_xx through zero into
        # tension, the only path here that reaches HS's -sigmat cut-off
        for k in range(1, S1["n_ax"] + 1):
            D = np.eye(3)
            D[0, 0] = 1.0 + S1["eps_ax"] * k / S1["n_ax"]
            steps.append(D @ Fc)
    elif path == "shear":
        q, ga = S1["n_quarter"], S1["gamma_a"]
        g = []    # per cycle: 0 -> +ga (q steps), +ga -> -ga (2q), -ga -> 0 (q)
        for _ in range(S1["cycles"]):
            g += [ga * k / q for k in range(1, q + 1)]
            g += [ga - 2 * ga * k / (2 * q) for k in range(1, 2 * q + 1)]
            g += [-ga + ga * k / q for k in range(1, q + 1)]
        for gk in g:
            S = np.eye(3)
            S[0, 1] = gk                                   # isochoric simple shear (det = 1)
            steps.append(S @ Fc)
    else:
        raise ValueError(path)
    return dict(F=Fs + steps, n_stage0=S1["n_conf"], conf_strain=a)


def _onepoint_model(ops, mat, dim, Fs, san_opts=()):
    X = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for a, x in enumerate(X, start=1):
        ops.node(a, *map(float, x))
        if dim == "ps":
            ops.fix(a, 0, 0, 1)
    _material(ops, mat, san_opts=san_opts)
    ops.element("LadrunoBrick", 1, *range(1, 9), 1, "-geom", "hypo")
    n = len(Fs)
    tags = 0
    for a, x in enumerate(X, start=1):
        for i in range(3 if dim == "3d" else 2):
            vals = [float(((F - np.eye(3)) @ np.array(x, float))[i]) for F in Fs]
            tags += 1
            ops.timeSeries("Path", tags, "-time", *[float(t) for t in range(n)], "-values", *vals)
            ops.pattern("Plain", tags, tags)
            ops.sp(a, i + 1, 1.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-13, 25, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


def _point_record(ops, mat):
    sig = [list(ops.eleResponse(1, "material", gp + 1, "stress")) for gp in range(8)]
    spread = float(np.max(np.ptp(np.array(sig), axis=0)))
    rec = dict(sig=sig[0], gp_spread=spread, state=list(ops.eleResponse(1, "material", 1, "state")))
    if mat == "sanisand":
        rec["substeps"] = list(ops.eleResponse(1, "material", 1, "substeps"))
        rec["substeps_gp"] = [int(ops.eleResponse(1, "material", gp + 1, "substeps")[0])
                              for gp in range(8)]
    rec["gp_rel_spread"] = spread / max(1.0, max(abs(v) for v in sig[0]))
    return rec


def fe_onepoint(mat, dim, path, san_opts=(), flip=True):
    ops = fe_engine()
    t0 = time.time()
    P = s1_path(mat, dim, path)
    Fs = P["F"]
    _onepoint_model(ops, mat, dim, Fs, san_opts)
    route = _stage(ops, mat, 0, [1])
    recs, flip_at = [], P["n_stage0"]
    for n in range(1, len(Fs)):
        if flip and n == flip_at + 1:
            route = _stage(ops, mat, 1, [1])
        ok = ops.analyze(1)
        if ok != 0:
            recs.append(dict(step=n, failed=ok))
            break
        r = _point_record(ops, mat)
        r["step"] = n
        recs.append(r)
    dF = [(Fs[n + 1] @ np.linalg.inv(Fs[n])).tolist() for n in range(len(Fs) - 1)]
    ok_steps = [r for r in recs if "failed" not in r]
    first_spread = next((r["step"] for r in ok_steps if r["gp_rel_spread"] > 1e-8), None)
    first_sub = next((r["step"] for r in ok_steps
                      if "substeps_gp" in r and len(set(r["substeps_gp"])) > 1), None)
    out = dict(material=mat, dim=dim, path=path, n_stage0=P["n_stage0"], flip_before_step=flip_at + 1,
               fe_self_consistency=dict(
                   first_step_gp_rel_spread_gt_1e8=first_spread,
                   first_step_gp_substeps_differ=first_sub,
                   max_gp_rel_spread=max(r["gp_rel_spread"] for r in ok_steps),
                   note="the 8 Gauss points of the homogeneous brick see the same F to round-off "
                        "(different dN/dX); where they part, the FE itself cannot hold rel 1e-8"),
               conf_strain=P["conf_strain"], flip_route=route + " (one element: all 8 points)",
               F=[F.tolist() for F in Fs], dF=dF, steps=recs, wall_s=round(time.time() - t0, 2),
               element="LadrunoBrick -geom hypo -formulation std, unit cube, every node "
                       "prescribed (zero free equations)")
    last = recs[-1]
    print(f"S1 {mat} {dim} {path}: {len(recs)} steps, sig_last {np.round(last.get('sig', []), 4)}"
          f" spread {last.get('gp_spread')} ({out['wall_s']} s)", flush=True)
    return out


def leg_onepoint() -> dict:
    legs = {}
    for mat in ("sanisand", "hs"):
        for dim in ("3d", "ps"):
            for path in ("triax", "oedo", "shear"):
                legs[f"{mat}_{dim}_{path}"] = fe_onepoint(mat, dim, path)
    payload = dict(_provenance("s1_onepoint", "fe_onepoint"), decks=decks(), paths=S1, legs=legs,
                   conventions=dict(
                       F="row-major 3x3 per step, F[0] = I; dF[n] = F[n+1] F[n]^-1; with "
                         "quasi-static dt, L[n] = (dF[n] - I) / dt",
                       steps="steps[k] is the committed state after F[k+1] (k = 0-based)",
                       sig="GP 1 of the brick, 3D Voigt (xx, yy, zz, xy, yz, zx), tension positive; "
                           "all 8 GPs agree to gp_spread",
                       state="SANISAND: 26 (eps_e 0-5, alpha 6-11, fabric 12-17, alpha_in 18-23, "
                             "e 24, dGamma 25); HS: its 7-vector 'state'",
                       substeps="SANISAND 'substeps' response [ME substeps, capHit] of the last "
                                "(converged) update, GP 1; substeps_gp: ME substeps at all 8 GPs"))
    _save("s1_onepoint", payload)
    return legs


# S1 option legs (critic AMEND, 2026-11-18): the admitted options, one-point FE, same paths.
S1_OPTIONS = {
    # SANISAND positional IntScheme (LadrunoSANISAND.cpp:64-66, oData[0], defaults :229-237)
    "sanisand_is2": ("sanisand", (2,), ("triax", "oedo", "shear")),
    "sanisand_is129": ("sanisand", (129,), ("triax", "oedo", "shear")),
    # -implex: companion scheme 1, which REQUIRES a cap -- "-implex on IntScheme 1 REQUIRES
    # -maxSubsteps > 0" (LadrunoSANISAND.cpp:3132) -- so -maxSubsteps 1000 (#792's cap, as the S3
    # v4 variant); -implexDt pseudo (ops_Dt = the LoadControl increment, 1 here), -implexFactor
    # fixed, -implexControl OFF, all defaults (LadrunoSANISAND.cpp:66-71)
    "sanisand_implex": ("sanisand", ("-implex", "-maxSubsteps", 1000), ("triax", "oedo", "shear")),
    "hs_implex": ("hs", ("-implex",), ("triax", "oedo", "shear")),
    # HSSmall through -smallStrain G0_ref gamma_0.7 (the deck of test_ladruno_hs_material.py:778-787)
    "hs_smallstrain": ("hs", ("-smallStrain", 140000.0, 1.5e-4), ("shear",)),
    # tension cut-off sigma_t (the deck has c = 10 kPa; 5 kPa < c cot(phi) = 12.8 kPa). Measured:
    # on triax / oedo / shear it is bit-identical to the deck (no principal stress goes tensile),
    # so the `tension` path is added, with the deck's own tension leg beside it.
    "hs_sigmat": ("hs", ("-sigmat", 5.0), ("triax", "oedo", "shear", "tension")),
    "hs_deck": ("hs", (), ("tension",)),
}


def leg_onepoint_options() -> dict:
    legs = {}
    for key, (mat, opts, paths) in S1_OPTIONS.items():
        for dim in ("3d", "ps"):
            for path in paths:
                r = fe_onepoint(mat, dim, path, san_opts=opts)
                r["options"] = list(opts)
                legs[f"{key}_{dim}_{path}"] = r
    payload = dict(_provenance("s1_onepoint_options", "fe_onepoint"), decks=decks(), paths=S1,
                   options={k: dict(material=v[0], opts=list(v[1]), paths=list(v[2]))
                            for k, v in S1_OPTIONS.items()},
                   legs=legs,
                   conventions="as s1_onepoint.json (F, dF, steps, fe_self_consistency); the "
                               "material is the deck plus `options`")
    _save("s1_onepoint_options", payload)
    return legs


# ----------------------------------------------------------------------------------------------
# S2: the K0 column (D19 domain, no translation) + a surcharge rho g l0 on the top, mode A QS
# ----------------------------------------------------------------------------------------------
S2 = dict(n_sur=(10, 20))


def leg_s2() -> dict:
    legs = {}
    q = RHO * G_GRAV * COL["L0"]
    for mat in ("sanisand", "hs"):
        for n in S2["n_sur"]:
            legs[f"{mat}_nsur{n}"] = fe_column(mat, surcharge=q, n_sur=n)
    # step-count sensitivity of the reference itself (the MPM may use another step count)
    sens = {}
    for mat in ("sanisand", "hs"):
        a = legs[f"{mat}_nsur{S2['n_sur'][0]}"]["after_surcharge"]["particles"]
        b = legs[f"{mat}_nsur{S2['n_sur'][1]}"]["after_surcharge"]["particles"]
        sens[mat] = max(abs(x["sxx"] - y["sxx"]) for x, y in zip(a, b)) / q
    payload = dict(_provenance("s2_column", "fe_column(surcharge)"), decks=decks(), legs=legs,
                   surcharge=q, nsur_sensitivity_over_rho_g_l0=sens,
                   note="surcharge as nodal top loads (tributary), ramped in n_sur LoadControl "
                        "steps after the flip and the 2-step hold; syy = -(rho g (l0 - Y) + q) "
                        "closed form; sxx, szz from the FE ratio, as K0")
    _save("s2_column", payload)
    print("S2 n_sur sensitivity (max |dsxx| / rho g l0):", sens, flush=True)
    return legs


# ----------------------------------------------------------------------------------------------
# S4: the low-g twin -- the one-point FE at each particle depth, flipped by a direct
# updateParameter (materialState, id 5), so Elastic2Plastic runs at that point
# ----------------------------------------------------------------------------------------------
S4 = dict(g=0.5)


def s4_g_threshold(rows: int = 1) -> float:
    """The largest g with `rows` particle rows (from the top, depth (k + 1/2) h/ppc) at p < Pmin
    at the flip: p = (1 + 2 K0)/3 rho g d with SANISAND's stage-0 K0 = nu / (1 - nu)."""
    k0 = EL_NU / (1.0 - EL_NU)
    d = (rows - 0.5) * COL["h"] / COL["ppc"]
    return SAN_PMIN / ((1.0 + 2.0 * k0) / 3.0 * RHO * d)


def fe_flip_point(syy: float, flip: str = "parameter") -> dict:
    """One SANISAND point: stage-0 oedometric compression (eps_xx = eps_zz = 0) to sigma_yy =
    syy in 10 steps (-geom hypo), then the flip, then one hold step."""
    ops = fe_engine()
    M = DP_K + 4.0 * DP_G / 3.0
    e = syy / M
    Fs = [np.eye(3)]
    for k in range(1, 11):
        F = np.eye(3)
        F[1, 1] = 1.0 + e * k / 10.0
        Fs.append(F)
    Fs.append(Fs[-1].copy())                     # the hold
    _onepoint_model(ops, "sanisand", "3d", Fs)
    _stage(ops, "sanisand", 0, [1])
    for _ in range(10):
        if ops.analyze(1) != 0:
            raise RuntimeError("S4 point stage 0 failed")
    pre = _point_record(ops, "sanisand")
    if flip == "parameter":
        ops.parameter(950)
        # ManzariDafalias::setParameter wants (name, matTag): argc < 2 returns -1 SILENTLY
        # (ManzariDafalias.cpp:931-937), so `... materialState` without the tag flips nothing.
        ops.addToParameter(950, "element", 1, "materialState", 1)
        ops.updateParameter(950, 1.0)
    else:
        ops.updateMaterialStage("-material", 1, "-stage", 1)
    post = _point_record(ops, "sanisand")
    ok = ops.analyze(1)
    hold = _point_record(ops, "sanisand") if ok == 0 else dict(failed=ok)
    p_pre = -sum(pre["sig"][:3]) / 3.0
    return dict(syy_target=syy, p_pre=p_pre, below_Pmin=bool(p_pre < SAN_PMIN), pre_flip=pre,
                after_flip=post, after_hold=hold)


def leg_s4(g: float = S4["g"]) -> dict:
    rows = []
    depths = [(k + 0.5) * COL["h"] / COL["ppc"] for k in range(int(COL["L0"] / (COL["h"] / COL["ppc"])))]
    for d in depths:
        r = fe_flip_point(-RHO * g * d)
        r["depth"] = d
        rows.append(r)
    nres = sum(r["below_Pmin"] for r in rows)
    payload = dict(_provenance("s4_lowg", "fe_flip_point"), decks=decks(), g=g, rows=rows,
                   rows_below_Pmin=nres, Pmin=SAN_PMIN,
                   g_threshold=dict(one_row=s4_g_threshold(1), two_rows=s4_g_threshold(2)),
                   flip_route="parameter element 1 materialState + updateParameter 1.0 "
                              "(ManzariDafalias id 5: Elastic2Plastic at every point of the element)",
                   rule="Elastic2Plastic (ManzariDafalias.cpp:5787-5795): tr sigma < 3 Pmin (p < "
                        "Pmin) -> sigma = sigma_n = Pmin I, alpha = alpha_n = 0; else sigma kept "
                        "(and eta <= Mc here, no Outside Bounding). Rows by depth d from the top: "
                        "d = (k + 1/2) 0.25, the same for the S2 twin (l0 10) and a 2 x 4-cell "
                        "column (l0 2, k < 8).")
    _save("s4_lowg", payload)
    print(f"S4 g={g}: {nres} row(s) below Pmin; thresholds {payload['g_threshold']}", flush=True)
    for r in rows[:3]:
        print("  d", r["depth"], "p_pre", round(r["p_pre"], 5), "post sig",
              [round(x, 6) for x in r["after_flip"]["sig"][:3]], "alpha",
              [round(x, 6) for x in r["after_flip"]["state"][6:9]], flush=True)
    return payload


# ----------------------------------------------------------------------------------------------
# D12: does SANISAND -implex -maxSubsteps 1 refuse at commit on one stage-1 shear step?
# ----------------------------------------------------------------------------------------------
def leg_d12(dgamma: float = 1.0e-3) -> dict:
    ops = fe_engine()
    P = s1_path("sanisand", "3d", "oedo")
    Fs = P["F"][:P["n_stage0"] + 1]
    S = np.eye(3)
    S[0, 1] = dgamma
    Fs.append(S @ Fs[-1])
    opts = ("-implex", "-maxSubsteps", 1)
    _onepoint_model(ops, "sanisand", "3d", Fs, san_opts=opts)
    _stage(ops, "sanisand", 0, [1])
    for _ in range(P["n_stage0"]):
        if ops.analyze(1) != 0:
            raise RuntimeError("D12 stage 0 failed")
    _stage(ops, "sanisand", 1, [1])
    ok = ops.analyze(1)
    stats = list(ops.eleResponse(1, "material", 1, "commitStats"))
    out = dict(_provenance("d12_trigger", "leg_d12"), opts=list(opts), dgamma=dgamma,
               analyze_return=ok, commitStats=stats,
               fires=bool(ok != 0),
               note="one LadrunoBrick -geom hypo, isotropic stage 0 to p0 = 100 kPa, "
                    "updateMaterialStage 1, one simple-shear step dgamma")
    _save("d12_trigger", out)
    print(f"D12: analyze returned {ok}; commitStats {stats}", flush=True)
    return out


# ----------------------------------------------------------------------------------------------
# S3: the strip footing (232j L4 domain), gravity + 20 kPa surcharge, staged, K9 ladder
# ----------------------------------------------------------------------------------------------
FOOT = dict(B=2.0, X=12.0, D=6.0, q0=20.0, ds=0.002, s_max=0.1, floor=5.0e-4,
            s_over_b=(0.01, 0.02, 0.05))
# The freeze's fixed probe order for the S3 SANISAND deck (first that reaches s = 0.05 B at r1
# AND r2 on the K9 ladder wins); options appended to SAN_PARAMS.
SAN_VARIANTS = {
    "deck": (),
    # v0 (ruling 2026-11-17, ahead of v1-v4): the TIMs model-of-record option set on our 18
    # constants, as tims-workbench/Tries/2d-model/model/strip.py emits it -- positional
    # "IntScheme": 129, "TanType": 0, ... "TolR": 1.0e-4 (:466, over TOYOURA_DM04's
    # "JacoType=1, TolF=1.0e-7" at :189), then (:1736-1758) "-Presidual" 0.0, "-honorTolR" 0,
    # "-flipAlphaIn" "init", "-Pmin" 0.0101, "-maxSubsteps" 2000,
    # "-sasHFloor", 1.0, "-sasReseatHyst", 1, "-sasSoftCap", 0.5 and "-sasTensionCutoff", 0.5, 1.0.
    "v0_tims_r1sep": (129, 0, 1, 1.0e-7, 1.0e-4,
                      "-Presidual", 0.0, "-honorTolR", 0, "-flipAlphaIn", "init",
                      "-Pmin", 0.0101, "-maxSubsteps", 2000,
                      "-sasHFloor", 1.0, "-sasReseatHyst", 1, "-sasSoftCap", 0.5,
                      "-sasTensionCutoff", 0.5, 1.0),
    "v1_me_cap1000": ("-maxSubsteps", 1000),
    "v2_sasme": (129,),
    "v3_cppm_refuse": (2, "-cppmOnFail", "refuse"),
    "v4_implex_cap1000": ("-implex", "-maxSubsteps", 1000),
}
K9 = [("KrylovNewton", 1.0e-6, 40), ("NewtonLineSearch", 1.0e-6, 60), ("KrylovNewton", 1.0e-5, 80)]


# LadrunoSANISAND `sasStats` columns (LadrunoSANISAND.cpp:5616-5640, per instance, cumulative,
# counted at commit -- :4969): 15 sas_refusals; WP-152 separation 36 sas_sepEntriesTension,
# 37 sas_sepEntriesLowP, 38 sas_sepExits, 39 sas_sepActive (1 while separated).
_SAS_COLS = dict(refusals=15, sep_entries_tension=36, sep_entries_lowp=37, sep_exits=38)


def _sas_counts(ops, elems) -> dict:
    """Separation events summed over every Gauss point of the model (8 per hex), from the
    material's `sasStats` response; sep_active_points = points separated at the end."""
    tot = {k: 0.0 for k in _SAS_COLS}
    active = 0
    for e in elems:
        for gp in range(8):
            v = list(ops.eleResponse(e, "material", gp + 1, "sasStats"))
            if len(v) <= 39:
                return dict(available=False, note=f"sasStats has {len(v)} columns")
            for k, c in _SAS_COLS.items():
                tot[k] += v[c]
            active += int(v[39] > 0.0)
    return dict(available=True, response="sasStats (LadrunoSANISAND.cpp:5616-5640)",
                separation_events=int(tot["sep_entries_tension"] + tot["sep_entries_lowp"]),
                **{k: int(v) for k, v in tot.items()}, sep_active_points=active)


def fe_footing_soil(mat, n_across_b, variant="deck", max_wall=None, curve=True):
    """Half model x in [0, 12], y in [-6, 0], h = B / n_across_b, LadrunoBrick B-bar plane-strain
    slab. Stage 0: gravity (-BrickW, rho 2.0, g 9.81) and the 20 kPa surcharge (tributary nodal
    loads on every top node) ramp together in 10 LoadControl steps; loadConst; the flip; 2 hold
    steps (LoadControl 0). Push: the rough footing (top nodes x <= B/2: u_y prescribed, u_x held)
    descends 2 mm per step to s = 0.05 B = 0.1 m on the K9 ladder (KrylovNewton 1e-6/40 ->
    NewtonLineSearch 1e-6/60 -> KrylovNewton 1e-5/80 relaxed -> halving, floor 5e-4 of ds).

    q(s) (the GATED observable, `q_at`, curve column 2) is 232f/232j's convention, as
    tests/mpm_implicit_reference.py fe_footing: (F + r_corr) / (B/2 t), F = -sum R_y over the
    footing nodes after `reactions`, r_corr = the surcharge lumped on those nodes (their
    tributary reaches h/2 past B/2, so q(0) = q0 (1 + h/B) + F0/(B/2 t), not q0) -- the MPM
    case's -cond_force(footing, y) / (B/2), which also carries the surcharge particles under the
    footing's nodes. Curve column 3: q0 + (F - F0) / (B/2 t), F0 = F at the zero-displacement
    start (a LoadControl 0 step with the footing SPs on), reported."""
    ops = fe_engine()
    f = FOOT
    h = f["B"] / n_across_b
    nx, ny = int(round(f["X"] / h)), int(round(f["D"] / h))
    xs = np.linspace(0.0, f["X"], nx + 1)
    ys = np.linspace(-f["D"], 0.0, ny + 1)
    t = h

    def nid(i, j, k):
        return 1 + i + (nx + 1) * (j + (ny + 1) * k)

    t0 = time.time()
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for k in (0, 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                n = nid(i, j, k)
                ops.node(n, float(xs[i]), float(ys[j]), k * t)
                if j == 0:
                    ops.fix(n, 1, 1, 1)
                else:
                    ops.fix(n, 1 if i in (0, nx) else 0, 0, 1)
    san_opts = SAN_VARIANTS[variant] if mat == "sanisand" else ()
    _material(ops, mat, san_opts=san_opts)
    ne = 0
    for j in range(ny):
        for i in range(nx):
            ne += 1
            ops.element("LadrunoBrick", ne, nid(i, j, 0), nid(i + 1, j, 0), nid(i + 1, j + 1, 0),
                        nid(i, j + 1, 0), nid(i, j, 1), nid(i + 1, j, 1), nid(i + 1, j + 1, 1),
                        nid(i, j + 1, 1), 1, "-formulation", "bbar", "-b", 0.0, -RHO * G_GRAV, 0.0)
    elems = list(range(1, ne + 1))
    route0 = _stage(ops, mat, 0, elems)
    top = [(i, k) for k in (0, 1) for i in range(nx + 1)]
    trib = {nid(i, ny, k): (h * t / 4.0) * (1 if i in (0, nx) else 2) for i, k in top}
    foot = [nid(i, ny, k) for i, k in top if xs[i] <= f["B"] / 2 + 1e-9]
    r_corr = f["q0"] * sum(trib[n] for n in foot)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", *elems, "-type", BRICKW)
    for n, A in trib.items():
        ops.load(n, 0.0, -f["q0"] * A, 0.0)
    _analysis(ops, tol=1.0e-6)
    ops.integrator("LoadControl", 0.1)
    ops.analysis("Static")
    for st in range(10):
        if ops.analyze(1) != 0:
            raise RuntimeError(f"footing stage-0 step {st + 1} failed")
    ops.loadConst("-time", 0.0)
    route1 = _stage(ops, mat, 1, elems)
    ops.integrator("LoadControl", 0.0)
    for st in range(2):
        if ops.analyze(1) != 0:
            raise RuntimeError(f"footing hold step {st + 1} failed")
    flip_check = None
    if mat == "sanisand":
        pmin, etamax = math.inf, 0.0
        for e in elems:
            s = ops.eleResponse(e, "material", 1, "stress")
            flip_check = _flip_inertness([[0.0, s[0], s[1], s[2], s[3]]])
            pmin = min(pmin, flip_check["p_min_at_flip"])
            etamax = max(etamax, flip_check["eta_max_at_flip"])
        flip_check = dict(p_min_after_hold=pmin, eta_max_after_hold=etamax, Pmin=SAN_PMIN, Mc=SAN_MC,
                          note="GP 1 of every element after the flip and hold")
    s_max = f["s_max"]
    uy0 = {n: ops.nodeDisp(n, 2) for n in foot}
    ux0 = {n: ops.nodeDisp(n, 1) for n in foot}
    for kk, n in enumerate(foot):
        ops.timeSeries("Path", 100 + kk, "-time", 0.0, 1.0, "-values", uy0[n], uy0[n] - s_max,
                       "-useLast")
        ops.pattern("Plain", 100 + kk, 100 + kk)
        ops.sp(n, 2, 1.0)
        ops.sp(n, 1, ux0[n])
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormUnbalance", 1.0e-6, 40, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.0)
    ops.analysis("Static")
    if ops.analyze(1) != 0:
        raise RuntimeError("footing zero step failed")
    ops.reactions()
    F0 = -sum(ops.nodeReaction(n, 2) for n in foot)
    width = f["B"] / 2 * t
    rows, dsc, nfail, nrelax, nhalv, stopped = [], f["ds"], 0, 0, 0, None
    tp = time.time()
    while True:
        s_now = ops.getTime() * s_max
        if s_now >= s_max - 1e-12:
            break
        if max_wall is not None and time.time() - t0 > max_wall:
            stopped = f"wall cap {max_wall} s"
            break
        dsc = min(dsc, s_max - s_now)
        ops.integrator("LoadControl", dsc / s_max)
        ok = False
        for i, (alg, tl, it) in enumerate(K9):
            ops.test("NormUnbalance", tl, it, 0)
            ops.algorithm(alg)
            if ops.analyze(1) == 0:
                ok, nrelax = True, nrelax + (i == 2)
                break
            nfail += 1
        if not ok:
            dsc *= 0.5
            nhalv += 1
            if dsc < f["floor"] * f["ds"]:
                stopped = f"halving floor at s = {s_now:.6f}"
                break
            continue
        dsc = min(2 * dsc, f["ds"])
        ops.reactions()
        F = -sum(ops.nodeReaction(n, 2) for n in foot)
        rows.append((ops.getTime() * s_max, (F + r_corr) / width, f["q0"] + (F - F0) / width))
    sas = _sas_counts(ops, elems) if mat == "sanisand" else None
    s = np.array([r[0] for r in rows]) if rows else np.zeros(1)
    q = np.array([r[1] for r in rows]) if rows else np.zeros(1)
    s_end = float(s[-1]) if rows else 0.0
    out = dict(material=mat, variant=variant, san_opts=list(san_opts), n_across_b=n_across_b, h=h,
               hexes=ne, element="LadrunoBrick -formulation bbar (small strain), plane-strain slab",
               stage0_route=route0, flip_route=route1, flip_check=flip_check, F0=F0,
               sas_counts=sas,
               q0=f["q0"], s_end=s_end, reached=bool(s_end >= f["s_max"] - 1e-9), stopped=stopped,
               q_at={f"{sb:.2f}": float(np.interp(sb * f["B"], s, q))
                     for sb in f["s_over_b"] if rows and sb * f["B"] <= s_end + 1e-12},
               q_max=float(q.max()) if rows else None,
               steps=len(rows), failed_attempts=nfail, relaxed=nrelax, halvings=nhalv,
               wall_s=round(time.time() - t0, 1), push_wall_s=round(time.time() - tp, 1),
               ladder=[list(x) for x in K9] + [["halving", "floor", f["floor"]]])
    if curve:
        out["curve"] = [[float(a), float(b), float(c)] for a, b, c in rows]
    print(f"FOOT {mat} {variant} n={n_across_b}: s_end {s_end:.4f} reached {out['reached']} "
          f"q_at {out['q_at']} fails {nfail} relaxed {nrelax} halvings {nhalv} sas {sas} "
          f"({out['wall_s']} s) {stopped or ''}", flush=True)
    return out


def leg_footing(mat, n_across_b, variant="deck", max_wall=None) -> dict:
    res = fe_footing_soil(mat, n_across_b, variant, max_wall=max_wall)
    name = f"s3_fe_{mat}_{variant}_n{n_across_b}" if mat == "sanisand" else f"s3_fe_{mat}_n{n_across_b}"
    payload = dict(_provenance(name, "fe_footing_soil"), decks=decks(), **res)
    _save(name, payload)
    return payload


# ----------------------------------------------------------------------------------------------
def main(argv=None) -> None:
    argv = list(sys.argv[1:] if argv is None else argv)
    what = argv[0] if argv else "k0"
    if what == "k0":
        leg_k0()
    elif what == "k0p":
        leg_k0p()
    elif what == "onepoint":
        leg_onepoint()
    elif what == "onepoint_options":
        leg_onepoint_options()
    elif what == "s2":
        leg_s2()
    elif what == "s4":
        leg_s4(float(argv[1]) if len(argv) > 1 else S4["g"])
    elif what == "d12":
        leg_d12()
    elif what == "footing":
        mw = None
        if "--max-wall" in argv:
            mw = float(argv[argv.index("--max-wall") + 1])
        leg_footing(argv[1], int(argv[2]), argv[3] if len(argv) > 3 and not argv[3].startswith("--")
                    else "deck", max_wall=mw)
    else:
        raise SystemExit(f"unknown leg {what!r}")


if __name__ == "__main__":
    main()
