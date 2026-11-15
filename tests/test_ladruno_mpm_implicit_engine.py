"""WP-232f: MPM mode A (particles as elements, the user's implicit analysis) through the
shipped binaries, against 232f-V's cases (tests/mpm_implicit_cases.txt; ADR 232 §8 P1 rows,
the pre-registered gate table of .claude/skills/ladruno-mpm/references/wp-briefs.md "232f").

The light rows run here, in Zone-A: the self-weight column (implicit quasi-static), the 1D
bar under Newmark (the W1 re-seed leg), the failed-step retry, the no-cell-crossing re-seed
(normal / forced / skipped), recorder synchrony, restart (database) both ways in both shipped
engines, every refusal, and engine parity. The fine-resolution gates (cantilever and footing
r2 / r3, the plate with a hole) are tests/test_ladruno_mpm_implicit_gates.py.

An expected value still `TBD` (232f-V filling it) is skipped loudly, never passed.

This file joins the MPM mutation family's suite (Ladruno_scripts/mutation_gate.py): ZERO
drops the particle element's internal force, SCALE the batched G2P's position increment.

Bootstrap: dist/bin/opensees.pyd on PYTHONPATH; OpenSees.exe from dist/bin or LADRUNO_TCL_EXE.
Written from tests/mpm_implicit_interface.md and the cases file only.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_mpm_implicit_harness as H  # noqa: E402

pytestmark = pytest.mark.zone_a

REPO = Path(__file__).resolve().parents[1]

LIGHT_RUN = ["column_selfweight_small_qs", "bar_implicit_2d"]
REPORT_RUN = ["cantilever_qs_r1", "cantilever_dyn_r1", "footing_mildcone_r1"]
PROTOCOLS = ["retry_bar", "reseed_bar_forced", "reseed_bar_skip", "recorder_sync_cantilever"]
REFUSALS = [c["name"] for c in H.parse_cases() if c["kind"] == "refusal"]


def _tcl_exe():
    exe = os.environ.get("LADRUNO_TCL_EXE")
    if exe:
        return exe if Path(exe).exists() else None
    cand = REPO / "dist" / "bin" / ("OpenSees.exe" if os.name == "nt" else "OpenSees")
    return str(cand) if cand.exists() else None


TCL_EXE = _tcl_exe()


def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    if not hasattr(ops, "mpm"):
        pytest.skip("this opensees build has no `mpm` command")
    return ops


_CACHE: dict = {}


def run(engine, name):
    key = (engine, name)
    if key not in _CACHE:
        c = H.case_by_name(name)
        if engine == "python":
            _CACHE[key] = H.run_python(c, _ops())
        else:
            if TCL_EXE is None:
                pytest.skip("classic OpenSees exe not found (build it, or set LADRUNO_TCL_EXE)")
            _CACHE[key] = H.run_tcl(c, TCL_EXE)
    return _CACHE[key]


def check(engine, name, allow_report=True):
    c = H.case_by_name(name)
    r = run(engine, name)
    res = H.evaluate(c, r, runs=lambda n: run(engine, n))
    failed = [f"{lab}: {det}" for lab, v, det in res if v == "fail"]
    tbd = [lab for lab, v, _ in res if v == "tbd"]
    for lab, v, det in res:
        print(f"[{engine}] {v.upper():6s} {lab}: {det}")
    assert not failed, "\n".join(failed)
    return res, tbd


def _skip_if_only_tbd(res, tbd):
    if tbd and all(v in ("tbd", "report") for _, v, _ in res):
        pytest.skip("expected values TBD (232f-V): " + "; ".join(tbd))


@pytest.mark.parametrize("engine", ["python", "tcl"])
@pytest.mark.parametrize("name", LIGHT_RUN)
def test_mode_a_case(engine, name):
    res, tbd = check(engine, name)
    _skip_if_only_tbd(res, tbd)


@pytest.mark.xfail(strict=True, reason="interface A4 / quirk WP-232b-ugimp-fixed-domains: the 20% "
                   "compression self-weight column misses the frozen 0.03 in the implicit QS leg too")
def test_mode_a_selfweight_column_20pct():
    check("python", "column_selfweight_qs")


def _band(c, r, obs):
    e = [x for x in c["expects"] if x["obs"] == obs][0]
    got = H.observable(c, r, e)[0]
    at, rt = e["tol"]
    return got, e["vals"][0], abs(got - e["vals"][0]) <= at + rt * abs(e["vals"][0])


@pytest.mark.parametrize("name, obs", [("cantilever_qs_r1", "tip_disp_y"),
                                       ("cantilever_dyn_r1", "first_peak_tip_disp_y"),
                                       ("cantilever_dyn_r1", "first_period_tip")])
def test_mode_a_r1_band(name, obs):
    """Review item 4: the coarse cantilever rows are `report` in the cases file, but they sit
    inside their pre-registered 5 % bands -- asserted here, so Zone-A enforces a physical P1
    band (and SCALE has a second mode-A position oracle)."""
    c = H.case_by_name(name)
    r = run("python", name)
    got, ref, ok = _band(c, r, obs)
    assert r["stopped"] is None and ok, f"{name} {obs}: {got} against {ref} (5 %)"


@pytest.mark.parametrize("name", REPORT_RUN)
def test_mode_a_report_case(name):
    """Coarse rows: status `report` -- computed and printed, failed only if the run stops."""
    r = run("python", name)
    assert r["stopped"] is None, f"{name}: the run stopped at step {r['stopped']}"
    res, tbd = check("python", name)
    _skip_if_only_tbd(res, tbd)


@pytest.mark.parametrize("name", PROTOCOLS)
def test_mode_a_protocol(name):
    res, tbd = check("python", name)
    _skip_if_only_tbd(res, tbd)


def test_retry_counts_and_schedule():
    """Freeze decision 7 + gate row 'failed-step retry': the forced failure is retried with a
    cut dt (no regrowth), G2P runs once per accepted step, and the accepted dt sequence is the
    reference's schedule."""
    r = run("python", "retry_bar")
    assert r["info_retries"] == 1, r["info_retries"]
    assert r["info_g2p"] == r["info_acceptedSteps"] == 4 + 392
    assert r["dtHistory"][:4] == [0.02] * 4
    assert all(abs(d - 0.01) < 1e-15 for d in r["dtHistory"][4:]), r["dtHistory"][4:8]


def test_reseed_counters():
    """The shipped run raises domainChange() on every step; the skip switch really skipped
    some (the bar's supports change only where a GIMP edge h + lp is crossed, so most steps
    of reseed_bar_skip have no support change)."""
    rn, rs = run("python", "reseed_bar_normal"), run("python", "reseed_bar_skip")
    assert rn["info_domainChanges"] == rn["info_acceptedSteps"]
    assert 0 < rs["info_domainChanges"] < rs["info_acceptedSteps"], rs["info_domainChanges"]


REFUSAL_RULE = {
    "refuse_integrator_centraldifference": "transient integrator class tag",
    "refuse_ladrunohht": "LadrunoHHT (33013)",
    "refuse_ladrunogeneralizedalpha": "LadrunoGeneralizedAlpha (33014)",
    "refuse_variabletransient": "VariableTransient is refused",
    "refuse_arclength": "arc-length is out of scope",
    "refuse_dispcontrol_gridnode": "DisplacementControl is refused",
}


@pytest.mark.parametrize("engine", ["python", "tcl"])
@pytest.mark.parametrize("name", REFUSALS)
def test_mode_a_refusal(engine, name):
    """Each refusal is loud, changes nothing, and is refused for ITS rule (review item 12:
    any exception used to pass)."""
    check(engine, name)
    text = run(engine, name).get("err_text", "")
    assert REFUSAL_RULE[name] in text, f"{name}: rule text missing from:\n{text[-1500:]}"


@pytest.mark.parametrize("engine", ["python", "tcl"])
@pytest.mark.parametrize("name", ["restart_bar", "restart_footing"])
def test_restart_wipe_restore(engine, name):
    """Save at n, restore (Python: `wipe` + restore in this process; Tcl: a second OpenSees.exe), then
    continue: particle state equal to the uninterrupted run (rel 1e-10)."""
    check(engine, name)


@pytest.mark.parametrize("name", ["restart_bar", "restart_footing"])
def test_restart_same_process_tcl(name):
    """Freeze decision 4: classic Tcl also restores after `wipe` in the same OpenSees.exe."""
    if TCL_EXE is None:
        pytest.skip("classic OpenSees exe not found (build it, or set LADRUNO_TCL_EXE)")
    c = H.case_by_name(name)
    r = H.run_tcl(c, TCL_EXE, inproc=True)
    rr = run("tcl", c["compares"][0][2])
    worst = max(H.field_relerr(r, rr, f) for f in H.PROTOCOL_FIELDS)
    assert worst <= 1e-10, f"{name}: Tcl wipe+restore differs (worst rel {worst:.3g})"


@pytest.mark.parametrize("name", ["restart_bar", "restart_footing"])
def test_restart_fresh_process_python(name):
    """Freeze decision 4: restore in a fresh second Python process as well."""
    c = H.case_by_name(name)
    r2 = H.run_python_two_process(c, _ops())
    ref = c["compares"][0][2]
    rr = run("python", ref)
    worst = max(H.field_relerr(r2, rr, f) for f in H.PROTOCOL_FIELDS)
    assert worst <= 1e-10, f"{name}: fresh-process restore differs (worst rel {worst:.3g})"


def test_engine_parity_bar():
    """Classic Tcl and Python give the same particle output (ADR 232 §8)."""
    rp, rt = run("python", "reseed_bar_normal"), run("tcl", "reseed_bar_normal")
    for f in ("x", "v", "sig"):
        assert H.field_relerr(rt, rp, f) <= 1e-10, f


def test_cantilever_qs_data():
    """Freeze decision 9 and the D2.3 study, DATA not gates: the cantilever QS tip under the
    mode-A default (renormalize, c_m 5e-2), the old default (exclude 1e-6, which diverges at
    step 7), exclude 1e-3 and ghost 1e-2 at c_m 1e-2; Newton iterations per step (the hypo
    tangent is not the exact linearisation). Printed for the handoff; only the default must
    converge."""
    ops = _ops()
    c = H.case_by_name("cantilever_qs_r1")
    rows = []
    for label, mc, en in (("renormalize, c_m 5e-2", None, None),
                          ("exclude, c_m 1e-6", 1e-6, ("exclude",)),
                          ("exclude, c_m 1e-3", 1e-3, ("exclude",)),
                          ("ghost 1e-2, c_m 1e-2", 1e-2, ("ghost", 1e-2))):
        r = H.run_python(c, ops, mass_cutoff=mc, empty_nodes=en)
        rows.append((label, r["hist_tip"][9.875][-1], r["stopped"]))
    print("\ncantilever_qs_r1 tip_disp_y (reference -2.3723):")
    for lab, tip, st in rows:
        print(f"  {lab:24s} {tip:.6f}  stopped={st}")
    for name, args in H.model_commands(c):
        getattr(ops, name)(*H._py_args(args))
    for name, args in H.analysis_commands(c):
        getattr(ops, name)(*H._py_args(args))
    its = []
    for _ in range(c["steps"]):
        ops.mpm("analyze", 1, c["dt"])
        its.append(int(ops.testIter()))
    print("  Newton iterations per step (default):", its)
    assert rows[0][2] is None, "the mode-A default must converge on the cantilever QS r1"


def test_condforce_is_the_reaction():
    """Footing diagnosis, step 2: on a laterally confined block compressed by a full-width
    gridvel, `mpm condForce` equals the particle stress times the loaded width (one QS step),
    for ElasticIsotropic and DruckerPrager. The footing observable measures what it says."""
    ops = _ops()
    for mat in ("elastic", "dp"):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        if mat == "elastic":
            ops.nDMaterial("ElasticIsotropic", 1, 1.0e5, 0.3)
        else:
            ops.nDMaterial("DruckerPrager", 1, 1.5e5, 15517.241379310342, 0.2, 0.21012821217029673,
                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 0.25, "-cells", 8, 10)
        ops.mpm("mode", "quasiStatic")
        ops.mpm("body", 1, "box", 0.0, 0.0, 2.0, 2.0, "-ppc", 2, "-mat", 1, "-rho", 1.0)
        ops.mpm("fix", "xmin", "x")
        ops.mpm("fix", "xmax", "x")
        ops.mpm("fix", "ymin", "y")
        ops.mpm("gridvel", "top", "box", -1.0, 2.0, 3.0, 1.0e9, "v", "free", -0.001)
        ops.constraints("Transformation")
        ops.numberer("RCM")
        ops.system("UmfPack")
        ops.test("NormDispIncr", 1e-12, 30)
        ops.algorithm("Newton")
        ops.integrator("LoadControl", 1.0)
        ops.analysis("Static")
        ops.mpm("analyze", 1, 1.0)
        cf = list(ops.mpm("condForce", "top"))
        syy = list(ops.mpm("particleData", "sig"))[1::3]
        load = 2.0 * sum(syy) / len(syy)
        assert max(syy) - min(syy) <= 1e-9 * abs(load), (mat, min(syy), max(syy))
        assert abs(cf[1] - load) <= 1e-9 * abs(load), (mat, cf[1], load)


def test_particle_tangent_against_finite_differences():
    """Footing diagnosis, step 1, kept as a regression: the particle element's tangent against
    central differences of its resisting force (`mpm debug tangentCheck`), ElasticIsotropic,
    plane strain, on the self-weight column -- undeformed (consistent to 3e-6: the hypo
    kinematics) and after five quasi-static steps (about 1e-3: the tangent is LadrunoBrick's
    hypo assembly, not the exact linearisation of the Hughes-Winget increment). The FD has
    converged when eps 1e-8 and 1e-9 agree."""
    ops = _ops()
    c = H.case_by_name("column_selfweight_small_qs")
    for name, args in H.model_commands(c):
        getattr(ops, name)(*H._py_args(args))
    for name, args in H.analysis_commands(c):
        getattr(ops, name)(*H._py_args(args))
    a = list(ops.mpm("debug", "tangentCheck", 1e-6, 1e-8))
    b = list(ops.mpm("debug", "tangentCheck", 1e-6, 1e-9))
    assert a[4] == 80 and a[0] < 1e-5 and abs(a[0] - b[0]) < 0.1 * a[0] + 1e-9, (a, b)
    for _ in range(5):
        ops.mpm("analyze", 1, c["dt"])
    a = list(ops.mpm("debug", "tangentCheck", 1e-6, 1e-8))
    b = list(ops.mpm("debug", "tangentCheck", 1e-6, 1e-9))
    assert a[0] < 5e-3 and abs(a[0] - b[0]) < 0.01 * a[0], (a, b)


def test_mode_a_rigid_translation():
    """Closed form, the G2P position increment: every grid node carries the same prescribed
    velocity (a gridvel over the whole grid), so every particle moves by exactly v * t with no
    strain and no stress -- quasi-static mode A over five steps. This is the mode-A oracle the
    MPM mutation family's SCALE mode (the batched G2P position increment x1.5) turns red; the
    other Zone-A mode-A rows read stresses, velocities or a second run of the same binary."""
    ops = _ops()
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e5, 0.3)
    ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 0.25, "-cells", 12, 12)
    ops.mpm("mode", "quasiStatic")
    ops.mpm("body", 1, "box", 1.0, 1.0, 2.0, 2.0, "-ppc", 2, "-mat", 1, "-rho", 1.0)
    v = (0.01, -0.02)
    ops.mpm("gridvel", "all", "box", -1.0, -1.0, 4.0, 4.0, "v", *v)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-12, 10)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    X = list(ops.mpm("particleData", "x"))
    ops.mpm("analyze", 5, 1.0)
    x = list(ops.mpm("particleData", "x"))
    sig = list(ops.mpm("particleData", "sig"))
    err = max(abs(x[k] - (X[k] + 5.0 * v[k % 2])) for k in range(len(x)))
    assert err < 1e-12, f"rigid translation off by {err}"
    assert max(abs(s) for s in sig) < 1e-6, max(abs(s) for s in sig)



# --------------------------------------------------------------------------- #
#  Review round 1 (PR #1083)                                                   #
# --------------------------------------------------------------------------- #
def _small_block(ops, mode="quasiStatic", mat="elastic", ndm=2, cells=(12, 12)):
    ops.wipe()
    ops.model("basic", "-ndm", ndm, "-ndf", ndm)
    if mat == "elastic":
        ops.nDMaterial("ElasticIsotropic", 1, 1.0e5, 0.3)
    else:
        ops.nDMaterial("DruckerPrager", 1, 1.5e5, 15517.241379310342, 0.2, 0.21012821217029673,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    ops.mpm("grid", 232, "-origin", *([0.0] * ndm), "-h", 0.25, "-cells", *cells[:ndm])
    ops.mpm("mode", mode)


def _analysis(ops, mode, handler="Transformation", integ=("Newmark", 0.5, 0.25), test=1e-12, iters=20):
    ops.constraints(handler, *((1.0e12, 1.0e12) if handler == "Penalty" else ()))
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", test, iters)
    ops.algorithm("Newton")
    if mode == "quasiStatic":
        ops.integrator("LoadControl", 1.0)
        ops.analysis("Static")
    else:
        ops.integrator(*integ)
        ops.analysis("Transient")


def test_plain_analyze_after_failed_mpm_analyze_is_refused():
    """Review item 1: after an `mpm analyze` that exhausts its cuts the particle elements are
    still armed; a plain `analyze` must not run them (no G2P, particles unchanged)."""
    ops = _ops()
    c = H.case_by_name("column_selfweight_small_qs")
    for name, args in H.model_commands(c):
        getattr(ops, name)(*H._py_args(args))
    for name, args in H.analysis_commands(c):
        getattr(ops, name)(*H._py_args(args))
    ops.mpm("analyze", 2, c["dt"])
    ops.mpm("cut", 2.0, 1)
    ops.mpm("debug", "failAttempt", 3, 5)
    with pytest.raises(Exception):
        ops.mpm("analyze", 1, c["dt"])
    x0 = list(ops.mpm("particleData", "x"))
    g0 = int(ops.mpm("info", "g2p"))
    rc = ops.analyze(1)
    assert rc != 0, "a plain analyze ran the MPM particle elements"
    assert list(ops.mpm("particleData", "x")) == x0 and int(ops.mpm("info", "g2p")) == g0


@pytest.mark.parametrize("handler", ["Transformation", "Penalty"])
@pytest.mark.parametrize("flip", [1.0, 0.0])
def test_transient_gridvel_translation(handler, flip):
    """Review item 2: a whole-grid gridvel in TRANSIENT mode A from rest. The prescribed nodes
    must start each step at v (not the projected particle velocity): every particle then
    moves at exactly v, FLIP and PIC alike, under the Transformation and the Penalty handler."""
    ops = _ops()
    _small_block(ops, mode="implicit")
    ops.mpm("body", 1, "box", 1.0, 1.0, 2.0, 2.0, "-ppc", 2, "-mat", 1, "-rho", 1.0)
    ops.mpm("flip", flip)
    v = (0.01, -0.02)
    ops.mpm("gridvel", "all", "box", -1.0, -1.0, 4.0, 4.0, "v", *v)
    _analysis(ops, "implicit", handler=handler)
    X = list(ops.mpm("particleData", "x"))
    tol = 1e-10 if handler == "Transformation" else 1e-6
    for k in range(1, 6):
        ops.mpm("analyze", 1, 0.1)
        x = list(ops.mpm("particleData", "x"))
        vp = list(ops.mpm("particleData", "v"))
        ex = max(abs(x[i] - (X[i] + 0.1 * k * v[i % 2])) for i in range(len(x)))
        ev = max(abs(vp[i] - v[i % 2]) for i in range(len(vp)))
        assert ex < tol and ev < tol * 100, (k, ex, ev)


def _peaks(t, y):
    out = []
    for k in range(1, len(y) - 1):
        if y[k] > y[k - 1] and y[k] >= y[k + 1] and y[k] > 0:
            mins = H._minima(t[k - 1:k + 2], [-a for a in y[k - 1:k + 2]])
            out.append(-mins[0][1] if mins else y[k])
    return out


def test_dynamic_rayleigh_betak():
    """The element guide's mandatory dynamic Rayleigh regression (betaK != 0, transient), on
    the Newmark bar: stiffness-proportional damping betaK = 0.1 at omega = 2 pi / 10 gives
    zeta = betaK omega / 2 = 0.0314, so successive positive velocity peaks of the first mode
    decay by exp(-2 pi zeta / sqrt(1 - zeta^2)) = 0.8209 (2 %)."""
    ops = _ops()
    c = H.case_by_name("bar_implicit_2d")
    for name, args in H.model_commands(c):
        getattr(ops, name)(*H._py_args(args))
    ops.mpm("rayleigh", 0.0, 0.1)
    for name, args in H.analysis_commands(c):
        getattr(ops, name)(*H._py_args(args))
    t, cv = [0.0], []
    s0 = list(ops.mpm("info", "summary"))
    cv.append(s0[1] / s0[0])
    for k in range(2500):
        ops.mpm("analyze", 1, 0.01)
        sm = list(ops.mpm("info", "summary"))
        t.append((k + 1) * 0.01)
        cv.append(sm[1] / sm[0])
    pk = _peaks(t, cv)
    assert len(pk) >= 2, pk
    zeta = 0.1 * (2 * 3.141592653589793 / 10.0) / 2
    want = __import__("math").exp(-2 * 3.141592653589793 * zeta / (1 - zeta * zeta) ** 0.5)
    assert abs(pk[1] / pk[0] - want) <= 0.02 * want, (pk[:3], want)


@pytest.mark.parametrize("integ", [("HHT", 0.9), ("GeneralizedAlpha", 1.0, 0.9)])
def test_bar_period_other_integrators(integ):
    """Review item 5: the whitelisted HHT and GeneralizedAlpha re-seed through their own
    domainChanged(); the Newmark bar's period (10, 1 %) holds under both."""
    ops = _ops()
    c = H.case_by_name("bar_implicit_2d")
    c["integrator_command"] = [str(x) for x in integ]
    r = H.run_python(c, ops)
    e = [x for x in c["expects"] if x["obs"] == "period_comx"][0]
    got = H.observable(c, r, e)[0]
    assert r["stopped"] is None and abs(got - 10.0) <= 0.1, got


def test_3d_rigid_translation():
    """Review item 5: mode A in 3D -- the whole-grid gridvel translation, closed form."""
    ops = _ops()
    _small_block(ops, ndm=3, cells=(6, 6, 6))
    ops.mpm("body", 1, "box", 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, "-ppc", 2, "-mat", 1, "-rho", 1.0)
    v = (0.01, -0.02, 0.015)
    ops.mpm("gridvel", "all", "box", -1.0, -1.0, -1.0, 3.0, 3.0, 3.0, "v", *v)
    _analysis(ops, "quasiStatic")
    X = list(ops.mpm("particleData", "x"))
    ops.mpm("analyze", 4, 1.0)
    x = list(ops.mpm("particleData", "x"))
    err = max(abs(x[k] - (X[k] + 4.0 * v[k % 3])) for k in range(len(x)))
    assert err < 1e-12, err
    assert max(abs(s) for s in ops.mpm("particleData", "sig")) < 1e-6


def test_3d_bar_period():
    """Review item 5: the vibrating bar in 3D (rollers on y and z faces), Newmark average
    acceleration, dt 0.02: period 10 within 1 %."""
    ops = _ops()
    c = H.case_by_name("bar_implicit_2d")
    c["ndm"] = 3
    c["origin"] = [0.0, 0.0, 0.0]
    c["cells"] = [30, 1, 1]
    c["bodies"] = [["1", "box", "0.0", "0.0", "0.0", "25.0", "1.0", "1.0", "ppc", "2", "mat", "1",
                    "rho", "1.0"]]
    c["fixes"] = c["fixes"] + [("zmin", "z"), ("zmax", "z")]
    c["dt"], c["steps"] = 0.02, 1000
    r = H.run_python(c, ops)
    e = [x for x in c["expects"] if x["obs"] == "period_comx"][0]
    got = H.observable(c, r, e)[0]
    assert r["stopped"] is None and abs(got - 10.0) <= 0.1, got


def test_mpm_reset_after_a_removed_particle_element():
    """Review item 6: `remove element` on a particle element, then `mpm reset`: no crash (the
    host had a null material slot), and the next `mpm analyze` is refused loudly."""
    ops = _ops()
    _small_block(ops)
    ops.mpm("body", 1, "box", 1.0, 1.0, 2.0, 2.0, "-ppc", 2, "-mat", 1, "-rho", 1.0)
    ops.mpm("gridvel", "all", "box", -1.0, -1.0, 4.0, 4.0, "v", 0.0, -0.001)
    _analysis(ops, "quasiStatic")
    ops.mpm("analyze", 1, 1.0)
    tag = int(min(ops.getEleTags()))
    ops.remove("element", tag)
    ops.mpm("reset")
    with pytest.raises(Exception):
        ops.mpm("analyze", 1, 1.0)


def test_plain_reset_keeps_particle_plastic_state():
    """Review item 7: a plain `reset` must not revert particle materials behind the host's back
    (it used to drop DruckerPrager's plastic history while x / F / strain stayed deformed). A
    block compressed by a gridvel: 3 steps, `reset`, 1 step equals 4 steps straight."""
    ops = _ops()
    out = []
    for do_reset in (False, True):
        _small_block(ops, mat="dp")
        ops.mpm("body", 1, "box", 0.0, 0.0, 2.0, 2.0, "-ppc", 2, "-mat", 1, "-rho", 1.0)
        ops.mpm("fix", "xmin", "x")
        ops.mpm("fix", "ymin", "y")
        ops.mpm("gridvel", "top", "box", -1.0, 2.0, 4.0, 1.0e9, "v", "free", -0.004)
        _analysis(ops, "quasiStatic")
        ops.mpm("analyze", 3, 1.0)
        if do_reset:
            ops.reset()
        ops.mpm("analyze", 1, 1.0)
        out.append(list(ops.mpm("particleData", "sig")))
    a, b = out
    assert max(abs(x - y) for x, y in zip(a, b)) <= 1e-9 * max(abs(x) for x in a), "plain reset reverted particles"


def test_retry_after_trial_on_a_plastic_state():
    """Review item 13: the forced failure now happens AFTER the trial (materials carry a trial
    state), on DruckerPrager in compression; the retried step equals the cut-dt run (1e-12).
    The compression rate (1 mm per step) and Newton (1e-10, 60 iterations) are chosen so the
    run itself never needs a natural cut: the one retry is the forced one."""
    ops = _ops()
    runs = []
    for variant in ("retry", "reference"):
        _small_block(ops, mat="dp")
        ops.mpm("body", 1, "box", 0.0, 0.0, 2.0, 2.0, "-ppc", 2, "-mat", 1, "-rho", 1.0)
        ops.mpm("fix", "xmin", "x")
        ops.mpm("fix", "ymin", "y")
        ops.mpm("gridvel", "top", "box", -1.0, 2.0, 4.0, 1.0e9, "v", "free", -0.001)
        _analysis(ops, "quasiStatic", test=1e-10, iters=60)
        if variant == "retry":
            ops.mpm("debug", "failAttempt", 3, 1, "afterTrial")
            ops.mpm("analyze", 6, 1.0)
            assert int(ops.mpm("info", "retries")) == 1
        else:
            ops.mpm("analyze", 2, 1.0)
            ops.mpm("analyze", 8, 0.5)
        runs.append({f: list(ops.mpm("particleData", f)) for f in ("x", "sig", "F")})
    for f in ("x", "sig", "F"):
        a, b = runs[0][f], runs[1][f]
        den = max(abs(v) for v in b) or 1.0
        assert max(abs(x - y) for x, y in zip(a, b)) / den <= 1e-12, f
