"""WP-232k (232k-V): checks the REFERENCES of tests/mpm_soil_cases.txt -- never MPM.

* the stored fork FE records (tests/data/mpm/soil/*.json) carry their provenance and were taken
  on the decks this file names;
* their closed-form parts reproduce: sigma_yy by equilibrium, the stage-0 K0 ratios
  (nu / (1 - nu), HS's K0nc), SANISAND's flip rule (p < Pmin resets to Pmin I, alpha = 0), the
  S4 g thresholds;
* the S1 F paths are consistent (F[n+1] = dF[n] F[n], the shear is isochoric), every Gauss
  point of the homogeneous brick agrees, and the S2 reference does not depend on its step count;
* the reference blocks of the cases file point at records that exist.

numpy only. No fork binary needed: the FE legs are re-run by
`python3.12 tests/mpm_soil_reference.py <leg>`.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mpm_soil_reference as R  # noqa: E402

pytestmark = pytest.mark.zone_a

DATA = R.DATA
K0_EL = R.EL_NU / (1.0 - R.EL_NU)


def _json(name):
    p = DATA / f"{name}.json"
    if not p.exists():
        pytest.skip(f"{p.name} not recorded yet")
    return json.loads(p.read_text(encoding="utf-8"))


def _provenance_ok(d):
    assert len(d["build"]) == 40 and all(c in "0123456789abcdef" for c in d["build"])
    assert d["deck"].startswith("tests/mpm_soil_reference.py::")
    assert d["units"].startswith("kN, t, m, s, kPa")


# ----------------------------------------------------------------------------------------------
# Decks and closed forms (no record needed)
# ----------------------------------------------------------------------------------------------


def test_decks_are_the_quoted_ones():
    assert R.SAN_PARAMS[:9] == [264.32, 0.3129, 0.6944, 1.33090, 0.71, 0.027, 0.83, 0.45, 101.0]
    assert R.SAN_PARAMS[-1] == R.RHO == 2.0
    assert math.isclose(R.san_g0(), 81588.19155115394, rel_tol=1e-14)
    assert math.isclose(R.EL_E, 214234.27337502, rel_tol=1e-13)
    # the gated DP shares SANISAND's stage-0 operator (orchestrator's ruling)
    assert math.isclose(R.DP_K, R.EL_E / (3 * (1 - 2 * R.EL_NU)), rel_tol=1e-15)
    assert math.isclose(R.DP_G, R.san_g0(), rel_tol=1e-14)
    assert R.DP_ARGS[2:4] == R.DPM_ARGS[2:4] == [0.2, 0.21012821217029673]
    assert R.HS_ARGS[:7] == [35000.0, 35000.0, 105000.0, 0.5, 10.0, 38.0, 8.0]
    assert "-sig0" not in R.HS_ARGS and "-sig03" not in R.HS_ARGS     # K3 / D6
    assert R.G_GRAV == 9.81 and math.isclose(R.RHO * R.G_GRAV * R.COL["L0"], 196.2)


def test_particle_lattice_and_cap():
    ys = R.column_particle_y()
    assert len(ys) == 40 and ys[0] == 0.125 and ys[-1] == 9.875
    assert math.isclose(R.cap_rho() * R.G_GRAV * R.CAP["t"], 20.0, rel_tol=1e-14)


def test_s4_thresholds_closed_form():
    k0 = K0_EL
    for rows in (1, 2):
        g = R.s4_g_threshold(rows)
        d = (rows - 0.5) * 0.25
        assert math.isclose((1 + 2 * k0) / 3 * R.RHO * g * d, R.SAN_PMIN, rel_tol=1e-13)
    assert R.s4_g_threshold(2) < R.S4["g"] < R.s4_g_threshold(1)     # exactly one row at g = 0.5


# ----------------------------------------------------------------------------------------------
# K0 columns
# ----------------------------------------------------------------------------------------------
_RATIO = {"elastic": K0_EL, "sanisand": K0_EL, "hs": R.HS_K0NC,
          "dp_mildcone": R.DPM_NU / (1 - R.DPM_NU)}


@pytest.mark.parametrize("name", ["k0_column", "k0_column_capped"])
def test_k0_records(name):
    d = _json(name)
    _provenance_ok(d)
    assert set(d["legs"]) == set(R.MATERIALS)
    q_top = 20.0 if name.endswith("capped") else 0.0
    for mat, L in d["legs"].items():
        assert L["g"] == 9.81 and L["rho"] == 2.0
        for cp in ("end_stage0", "after_flip", "after_hold"):
            P = L[cp]["particles"]
            assert [p["Y"] for p in P] == R.column_particle_y()
            for p in P:   # equilibrium, closed form
                assert math.isclose(p["syy"], -(2.0 * 9.81 * (10.0 - p["Y"]) + q_top),
                                    rel_tol=1e-13)
            assert L[cp]["centre_syy_error"] < 1e-9 * 196.2      # FE centre = closed form
            assert L[cp]["gp_spread"] < 1e-9 * 196.2             # 1D: constant per element
        if mat in _RATIO:   # elastic K0 states are kept through the flip and the hold
            for cp in ("end_stage0", "after_flip", "after_hold"):
                for p in L[cp]["particles"]:
                    assert abs(p["ratio_xx"] - _RATIO[mat]) < 1e-9
                    assert abs(p["ratio_zz"] - _RATIO[mat]) < 1e-9
        if mat == "dp":     # outside the cone at eta 0.855: it returns at the first hold step
            assert all(abs(p["ratio_xx"] - K0_EL) < 1e-9 for p in L["end_stage0"]["particles"])
            r = [p["ratio_xx"] for p in L["after_hold"]["particles"]]
            assert min(r) > K0_EL + 1e-3 and max(r) < 0.5
            assert "value 2" in L["flip_route"]
        if mat == "sanisand":
            fc = L["flip_check"]
            assert fc["elastic2plastic_inert"] and fc["p_min_at_flip"] > fc["Pmin"]
            assert fc["eta_max_at_flip"] < fc["Mc"]
            assert "ONE element" in L["flip_route"]


# ----------------------------------------------------------------------------------------------
# S1
# ----------------------------------------------------------------------------------------------


def test_s1_paths_and_records():
    d = _json("s1_onepoint")
    _provenance_ok(d)
    assert len(d["legs"]) == 12
    for key, L in d["legs"].items():
        F = [np.array(x) for x in L["F"]]
        dF = [np.array(x) for x in L["dF"]]
        assert np.array_equal(F[0], np.eye(3))
        for n in range(len(dF)):
            assert np.allclose(dF[n] @ F[n], F[n + 1], rtol=0, atol=1e-14)
        if L["dim"] == "ps":
            for f in F:
                assert f[2, 2] == 1.0 and f[0, 2] == f[1, 2] == f[2, 0] == f[2, 1] == 0.0
        if L["path"] == "shear":
            fc = F[L["n_stage0"]]
            for f in F[L["n_stage0"]:]:
                assert math.isclose(np.linalg.det(f), np.linalg.det(fc), rel_tol=1e-14)
        assert len(L["steps"]) == len(F) - 1 and all("failed" not in s for s in L["steps"])
        assert L["flip_before_step"] == L["n_stage0"] + 1
        sc = L["fe_self_consistency"]
        if L["material"] == "hs":       # HS: every Gauss point agrees to round-off
            assert sc["max_gp_rel_spread"] < 1e-12 and sc["first_step_gp_rel_spread_gt_1e8"] is None
        else:
            # SANISAND scheme 1 (TolE 1e-4): round-off in dN/dX parts the 8 points after a few
            # plastic steps, sometimes BEFORE their substep counts differ; the FE reference itself
            # then holds only ~1e-4, inside 10 x TolE = 1e-3.
            assert sc["max_gp_rel_spread"] < 1e-3
            for s in L["steps"]:
                assert len(s["substeps"]) == 2 and len(s["substeps_gp"]) == 8
            for s in L["steps"][:L["n_stage0"] + 2]:   # stage 0 and the first two plastic steps
                assert s["gp_rel_spread"] < 1e-8
        # stage 0 reaches the confinement p0 = 100 kPa (3D) within the hypo/small-strain gap
        sig0 = L["steps"][L["n_stage0"] - 1]["sig"]
        if L["dim"] == "3d":
            assert abs(-sum(sig0[:3]) / 3 - 100.0) < 0.1


# ----------------------------------------------------------------------------------------------
# S2, S4, D12
# ----------------------------------------------------------------------------------------------


def test_s2_record():
    d = _json("s2_column")
    _provenance_ok(d)
    q = d["surcharge"]
    assert math.isclose(q, 196.2)
    for key, L in d["legs"].items():
        for p in L["after_surcharge"]["particles"]:
            assert math.isclose(p["syy"], -(2.0 * 9.81 * (10.0 - p["Y"]) + q), rel_tol=1e-13)
        if key.startswith("sanisand"):   # K0 loading stays on the cone axis: elastic, exact ratio
            assert all(abs(p["ratio_xx"] - K0_EL) < 1e-9 for p in L["after_surcharge"]["particles"])
    assert d["nsur_sensitivity_over_rho_g_l0"]["sanisand"] < 1e-12
    assert d["nsur_sensitivity_over_rho_g_l0"]["hs"] < 1e-3      # far below the 0.02 gate


def test_s4_record():
    d = _json("s4_lowg")
    _provenance_ok(d)
    assert d["g"] == 0.5 and d["rows_below_Pmin"] == 1
    for r in d["rows"]:
        post = r["after_flip"]
        if r["below_Pmin"]:      # Elastic2Plastic: sigma = Pmin I, alpha = 0
            assert post["sig"][:3] == [-R.SAN_PMIN] * 3 and post["sig"][3:] == [0.0] * 3
            assert post["state"][6:12] == [0.0] * 6
        else:                    # sigma kept, eta <= Mc
            assert post["sig"] == r["pre_flip"]["sig"]
        assert "failed" not in r["after_hold"]


def test_d12_trigger_fires():
    d = _json("d12_trigger")
    _provenance_ok(d)
    assert d["fires"] and d["analyze_return"] != 0
    assert d["opts"] == ["-implex", "-maxSubsteps", 1]


@pytest.mark.parametrize("n", [4, 8, 16])
def test_s3_hs_records(n):
    d = _json(f"s3_fe_hs_n{n}")
    _provenance_ok(d)
    assert d["reached"] and d["stopped"] is None
    assert set(d["q_at"]) == {"0.01", "0.02", "0.05"}
    q = [d["q_at"][k] for k in ("0.01", "0.02", "0.05")]
    assert q[0] < q[1] < q[2]


# ----------------------------------------------------------------------------------------------
# The cases file's reference blocks
# ----------------------------------------------------------------------------------------------


def _blocks():
    out, cur, name = {}, None, None
    for raw in (HERE / "mpm_soil_cases.txt").read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        tok = line.split()
        if tok[0] == "case":
            name, cur = tok[1], {}
        elif tok[0] == "end":
            out[name], cur = cur, None
        elif cur is not None:
            cur.setdefault(tok[0], []).append(tok[1:])
    return out


def test_reference_blocks_point_at_records():
    refs = {k: v for k, v in _blocks().items() if v.get("kind") == [["reference"]]}
    assert {"k0_ref_column", "k0_ref_column_capped"} <= set(refs)
    for name, b in refs.items():
        for (path,) in b["data"]:
            assert (R.REPO / path).exists(), f"{name}: {path} missing"
    mats = {int(t[0]): t[1:] for t in refs["k0_ref_column"]["material"]}
    assert [float(x) for x in mats[2][1:]] == R.DP_ARGS
    assert [float(x) for x in mats[3][1:]] == R.SAN_PARAMS
    assert [float(x) for x in mats[1][1:]] == [R.EL_E, R.EL_NU]


def test_k0p_record():
    d = _json("k0p_column_capped")
    _provenance_ok(d)
    assert set(d["legs"]) == set(R.MATERIALS)
    for mat, L in d["legs"].items():
        for cp in ("after_preshear", "after_hold2"):
            for p in L[cp]["particles"]:
                assert math.isclose(p["syy"], -1.3 * (2.0 * 9.81 * (10.0 - p["Y"]) + 20.0),
                                    rel_tol=1e-13)
            assert L[cp]["centre_syy_error"] < 1e-6 * 196.2
        if mat != "elastic":
            assert "every" in L["flip_route"]
    am = d["legs"]["sanisand"]["alpha_moved"]
    assert am["fraction"] >= 0.5 and am["min_norm_dalpha"] > 1e-3      # the pre-shear is plastic
    # nu -> 0.05 moved SANISAND off its stage-0 ratio; elastic keeps it
    r = [p["ratio_xx"] for p in d["legs"]["sanisand"]["after_hold2"]["particles"]]
    assert max(r) < K0_EL - 1e-2
    assert all(abs(p["ratio_xx"] - K0_EL) < 1e-9 for p in d["legs"]["elastic"]["after_hold2"]["particles"])
    assert "poissonRatio" in d["legs"]["sanisand"]["preshear_params"]
    assert "shearModulus" in d["legs"]["dp"]["preshear_params"]


def test_s1_options_record():
    d = _json("s1_onepoint_options")
    _provenance_ok(d)
    base = _json("s1_onepoint")["legs"]
    legs = d["legs"]
    assert len(legs) == 6 * 3 + 6 + 2 + 8 + 2
    for key, L in legs.items():
        assert all("failed" not in s for s in L["steps"]), key
        F = [np.array(x) for x in L["F"]]
        for n, dF in enumerate(L["dF"]):
            assert np.allclose(np.array(dF) @ F[n], F[n + 1], rtol=0, atol=1e-14)
    # IntScheme 129's points never part; HS options never part
    for key, L in legs.items():
        if key.startswith(("sanisand_is129", "hs_")):
            assert L["fe_self_consistency"]["max_gp_rel_spread"] < 1e-12, key
    # vacuity facts recorded in the cases block
    for dim in ("3d", "ps"):
        for p in ("triax", "oedo", "shear"):
            assert legs[f"sanisand_implex_{dim}_{p}"]["steps"][-1]["sig"] == \
                base[f"sanisand_{dim}_{p}"]["steps"][-1]["sig"]
            assert legs[f"hs_sigmat_{dim}_{p}"]["steps"][-1]["sig"] == \
                base[f"hs_{dim}_{p}"]["steps"][-1]["sig"]
        st = legs[f"hs_sigmat_{dim}_tension"]["steps"][-1]["sig"]
        sd = legs[f"hs_deck_{dim}_tension"]["steps"][-1]["sig"]
        assert abs(sum(st[:3]) / 3 - 5.0) < 1e-6       # mean stress held at sigma_t
        assert abs(sum(sd[:3]) / 3) < 1e-6              # the deck: sigma_t = 0
