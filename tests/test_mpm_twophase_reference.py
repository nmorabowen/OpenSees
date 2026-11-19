"""WP-232l (232l-V): checks the REFERENCES of the two-phase gates -- never MPM.

Round 1: the closed forms (Terzaghi copy against its source, the ramp Duhamel, Mandel's series
against the independent Laplace solution) and the provenance and shape of the stored LadrunoUP
records (tests/data/mpm/twophase/*.json). numpy + scipy only; no fork binary needed. The FE legs
are re-run by `python3.12 tests/mpm_twophase_reference.py <leg>`.
"""
from __future__ import annotations

import ast
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mpm_twophase_reference as R  # noqa: E402

pytestmark = pytest.mark.zone_a


def _json(name):
    p = R.DATA / f"{name}.json"
    if not p.exists():
        pytest.skip(f"{p.name} not recorded yet")
    return json.loads(p.read_text(encoding="utf-8"))


def _provenance_ok(d):
    assert len(d["build"]) == 40 and all(c in "0123456789abcdef" for c in d["build"])
    assert d["deck"].startswith("tests/mpm_twophase_reference.py::")
    assert d["units"].startswith("kN, t, m, s, kPa")


def _source_functions(path, names):
    """Exec the named top-level functions of a test module without importing it (its import
    binds an engine, test_ladruno_up_element_analytic.py:74-79)."""
    tree = ast.parse(Path(path).read_text(encoding="utf-8"))
    ns = {"np": np}
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in names:
            exec(compile(ast.Module([node], []), str(path), "exec"), ns)
    return ns


def test_terzaghi_copy_matches_its_source():
    ns = _source_functions(HERE / "test_ladruno_up_element_analytic.py",
                           {"terzaghi_series", "terzaghi_U"})
    for Z in (0.0, 0.3, 0.9):
        for Tv in (0.01, 0.1, 0.5):
            assert ns["terzaghi_series"](Z, Tv) == R.terzaghi_series(Z, Tv)
    for Tv in (0.01, 0.1, 0.5):
        assert ns["terzaghi_U"](Tv) == R.terzaghi_U(Tv)


def test_terzaghi_case_is_the_fe_tests():
    c = R.TERZ
    assert (c["E"], c["nu"], c["rho"], c["Kf"], c["poro"], c["rhoF"], c["kbar"], c["L"],
            c["q"], c["nely"], c["dt"], c["gamma"], c["beta"]) == (
        1.0e4, 0.2, 2.0, 4444.0, 0.4, 1.0, 1.0e-4, 10.0, 10.0, 10, 0.02, 0.6, 0.3025)


def test_terzaghi_ramp_limits():
    # ramp -> step; long after the ramp the ramp response lags the step by ~ Tvr/2
    assert abs(R.terzaghi_ramp(0.4, 0.2, 1e-6) - R.terzaghi_series(0.4, 0.2)) < 1e-5
    a = R.terzaghi_ramp(0.5, 0.3, 0.02)
    b = R.terzaghi_series(0.5, 0.3 - 0.01)
    assert abs(a - b) < 2e-4
    k = R.terzaghi_constants(1e4, 0.2, 4444.0, 0.4, 1e-4, 10.0, 10.0)
    assert math.isclose(R.terzaghi_settle_ramp(50.0, 0.02, k["s0"], k["sinf"]), k["sinf"],
                        rel_tol=1e-12)


def test_mandel_closed_form_selfcheck():
    out = R.closed_selfcheck()
    assert out["mandel_c_vs_cd_form"] < 1e-14
    assert out["mandel_series_vs_laplace_max_over_p0"] < 1e-8
    assert out["mandel_p_0plus_uniform_err"] < 1e-3


def test_mandel_peak_meets_the_kickoff_floor():
    mc = R.mandel_constants(**R.MANDEL_SOIL)
    pk, tv = R.mandel_peak(mc, nt=1500)
    assert pk >= 1.05 and 0.05 < tv < 0.09


def test_cases_file_reference_blocks_point_at_records():
    text = (HERE / "mpm_twophase_cases.txt").read_text(encoding="utf-8")
    for line in text.splitlines():
        if line.startswith("data ") and "<dt>" not in line:
            assert (R.REPO / line.split()[1]).exists(), line


def test_records_t2_terzaghi():
    d = _json("t2_terzaghi_fe")
    _provenance_ok(d)
    for run in d["runs"]:
        assert [r["Tv"] for r in run["rows"]] == list(R.TV_T2)
        for r in run["rows"][1:]:
            assert r["p_err_gimp_points"] < 0.02 and r["U_err"] < 0.02


def test_records_t3_mandel():
    d = _json("t3_mandel_fe")
    _provenance_ok(d)
    assert "UNVERIFIED" in d["closed_form_note"]
    assert d["selfcheck"]["mandel_series_vs_laplace_max_over_p0"] < 1e-8
    pk = d["closed_form_peak_ramp"]["p_over_p0"]
    assert pk >= 1.05 and math.isclose(d["peak_gate"], 1 + 0.5 * (pk - 1), rel_tol=1e-15)
    assert len(d["gated"]) == 2
    for run in d["gated"]:
        assert run["t_ramp"] == R.MANDEL["t_ramp"] and run["dt"] == 0.02
        assert [r["tv"] for r in run["rows"]] == list(R.MANDEL_TVS)
        assert all(b == 0.02 for b in run["band"])          # FE error < 0.01 everywhere


def test_records_t4_t5():
    h = _json("t4_hydro_fe")
    _provenance_ok(h)
    for run in h["gated"] + h["reported"]:
        assert run["p_err_over_gwH"] < 1e-6
        # sigma' is element-constant: the point error is the h/4 (GIMP) and h/(2 sqrt 3) offsets
        hL = 1.0 / run["nely"]
        assert math.isclose(run["sig_err_gimp_over_rhopgH"], hL / 4, rel_tol=1e-5)
        assert math.isclose(run["sig_err_gauss_over_rhopgH"], hL / (2 * math.sqrt(3)), rel_tol=1e-5)
    b = _json("t5_b4_fe")
    _provenance_ok(b)
    cb = {json.dumps(r["stab"]): r["CB_lap"] for r in b["runs"]}
    assert cb['"off"'] >= 0.5 and cb['["auto", 0.25]'] <= 0.3


# ----------------------------------------------------------------------------------------------
# T0 oracle and T1 records
# ----------------------------------------------------------------------------------------------
def test_t1_records_and_oracle_reproduce_ladrunoup():
    """The oracle's particle element (interface §1), summed over the Gauss-point particles under
    Newmark's predictor from the P5 state, reproduces LadrunoUP's A and b on every T1 leg."""
    d = _json("t1_fe")
    _provenance_ok(d)
    assert len(d["legs"]) == 12
    for name, leg in d["legs"].items():
        ea, eb = R.oracle_vs_fe_t1(leg)
        assert ea < 1e-13 and eb < 1e-13, f"{name}: A {ea:.2e}, b {eb:.2e}"


def test_t1_dp_leg_stays_elastic():
    d = _json("t1_fe")
    for m in R.T1_MESHES:
        for st in ("off", "auto"):
            a = np.array(d["legs"][f"{m}_elastic_{st}"]["dU"])
            b = np.array(d["legs"][f"{m}_dp_{st}"]["dU"])
            assert np.max(np.abs(a - b)) <= 1e-12 * np.max(np.abs(a))


def test_t0_blocks_have_the_route_e_pattern():
    """[K, -Q; 0, H], [C_Ray, 0; Q^T, S*], [M, 0; 0, 0]; the Q in K is the Q^T in C."""
    for b in R.t0_particle_dump():
        K, C, M = (np.array(b[k]) for k in ("K", "C", "M"))
        nen = len(b["nodes"])
        pI = [3 * k + 2 for k in range(nen)]
        uI = [3 * k + i for k in range(nen) for i in range(2)]
        assert np.all(K[np.ix_(pI, uI)] == 0.0)
        assert np.all(C[np.ix_(uI, pI)] == 0.0)
        assert np.all(M[pI, :] == 0.0) and np.all(M[:, pI] == 0.0)
        assert np.max(np.abs(-K[np.ix_(uI, pI)] - C[np.ix_(pI, uI)].T)) < 1e-15
        assert np.max(np.abs(np.sum(np.array(b["dS"]), axis=0))) < 1e-14     # partition of unity
        assert abs(sum(b["S"]) - 1.0) < 1e-14


def test_t0_record_is_the_oracle():
    d = _json("t0_oracle")
    assert d["source"].startswith("tests/mpm_twophase_interface.md")
    fresh = {"t0_i": R.t0_particle_dump(R.T0A), "t0_ii": R.t0_modeA(R.T0A),
             "t0_iii": R.t0_routeS(R.T0S)}
    assert json.loads(json.dumps(fresh)) == {k: d[k] for k in fresh}


def test_t0_drained_nodes_and_p_activity():
    """Mode A: the pfix nodes hold p = 0 after each step. Route S: under ruling 3d3699a50 the
    p-active set is the 25 real nodes at every step (the light ghost row is dropped)."""
    out = R.t0_modeA(R.T0A)
    for h in out["history"][1:]:
        g = h["grid"]
        for k, idx in enumerate(g["node_idx"]):
            if idx[1] == 2:
                assert g["U"][3 * k + 2] == 0.0
    s = R.t0_routeS(R.T0S)
    for h in s["history"][1:]:
        assert len(h["grid"]["nodes"]) == 25
        assert min(i[1] for i in h["grid"]["node_idx"]) == 0


# ----------------------------------------------------------------------------------------------
# §7b: P5a, fix 2
# ----------------------------------------------------------------------------------------------
def test_p5a_projection_is_exact_on_grid_fields():
    """M^-1 P^T V P = I: a grid field interpolated to the particles projects back to itself."""
    grid = R.OGrid(R.T0A["x0"], R.T0A["h"], R.T0A["cells"])
    X, V0, m, lpv = R.lattice(grid, *R.T0A["box"], R.T0A["ppc"], R.T0A["rho"])
    sups, _ = R.renormalize([grid.support(x, [lpv, lpv]) for x in X], [m] * len(X), 0.05)
    nodes = sorted({i for s in sups for i, _, _ in s}, key=lambda t: t[::-1])
    rng = np.random.default_rng(7)
    g = {n: rng.normal() for n in nodes}
    vals = [[sum(S * g[i] for i, S, _ in s)] for s in sups]
    Xh, eig = R.consistent_l2([[(i, S) for i, S, _ in s] for s in sups], [V0] * len(X), vals, nodes)
    assert eig > 0
    assert max(abs(Xh[a, 0] - g[n]) for a, n in enumerate(nodes)) < 1e-12


def test_t0_and_t1_use_p5a():
    d = _json("t0_oracle")
    for h in d["t0_ii"]["history"][1:] + d["t0_iii"]["history"][1:]:
        assert h["grid"]["M_eigmin"] > 1e-3          # SPD with margin at every step
    t = _json("t1_fe")
    for leg in t["legs"].values():
        assert leg["projection"].startswith("P5a") and leg["M_eigmin"] > 1e-3


def test_sealed_column_invariant_bit_for_bit():
    """§7b: motionless, sealed, kbar = 0, pd = pdd = 0 -- particle p unchanged bit for bit, for a
    field the grid cannot represent (its invisible residual is O(1))."""
    worst, resid, _ = R.sealed_column_check(steps=3)
    assert worst == 0.0
    assert resid > 1.0


def test_fix2_affine_image_reproduces_hydrostatics():
    """§7b fix 2: with the affine even image on the impervious base, the hydrostatic field is
    reproduced (values, gradients: zero Darcy flux) and survives the P5a projection; P9's plain
    even image does not (negative control)."""
    e = R.fix2_hydrostatic_check()
    assert max(e) < 1e-13
    bad = R.fix2_hydrostatic_check(affine=False)
    assert bad[0] > 1e-3 and bad[1] > 0.1


def test_t0_iii_section_7c():
    """§7c: the light ghost row is dropped whole (u side), the deck's dt is admitted by
    s x the drained pencil."""
    d = _json("t0_oracle")
    hist = d["t0_iii"]["history"][1:]
    assert hist[0]["grid"]["u_dropped"] == []
    # corrected 7c (1): only ghost nodes (i = -1 or j = -1) are ever dropped
    assert all(i[0] == -1 or i[1] == -1 for h in hist for i in h["grid"]["u_dropped"])
    assert all(len(h["grid"]["nodes"]) == 25 for h in hist)          # p-active: the real nodes
    assert d["case_routeS"]["dt"] <= d["t0_iii"]["safety"] * min(d["t0_iii"]["dt_drained_pencil"])


def test_r_rows_route_e_record():
    """Reported rows (never gated): provenance and shape only."""
    d = _json("r_rows_route_e")
    _provenance_ok(d)
    assert d["reported_never_gated"] is True
    assert len(d["r1_first_settlement"]) == 36
    assert {(r["basis"], r["pressureBlend"]) for r in d["r2_projection_diffusion"]} == {
        ("gimp", 1.0), ("gimp", 0.0), ("bspline2", 1.0), ("bspline2", 0.0)}
    assert len(d["r4_cost"]) == 6
    assert {(r["basis"], r["translate"]) for r in d["r3_cell_crossing"]} == {
        ("gimp", True), ("gimp", False), ("bspline2", True), ("bspline2", False)}
    assert d["r5_undrained_simple_shear"]["status"].startswith("not run")


def test_r_rows_route_s_record():
    """Reported rows (never gated): provenance, the u(0+) definition, shape."""
    d = _json("r_rows_route_s")
    _provenance_ok(d)
    assert d["reported_never_gated"] is True
    assert "20 undrained wave transits" in d["definition_u0"]
    assert len(d["r1_first_settlement"]) == 30          # 2 Kf x 3 stab x (FE + 2 routes x 2 bases)
    s_rows = [r for r in d["r2_cost_dt"] if r["route"] == "S"]
    assert all(r["dt_run"] <= r["dt_admitted_s0.5"] for r in s_rows)
    assert "cost_and_dt_T2column" in d["r3_route_comparison"]


def test_t13_free_translation_is_exact():
    """§7c (1) corrected: a free gimp two-phase block in rigid translation keeps v, sigma', p."""
    r = R.t13_translation()
    assert r["v_dev"] < 1e-12 and r["sig_dev"] < 1e-9 and r["p_dev"] < 1e-9 and r["x_dev"] < 1e-12
