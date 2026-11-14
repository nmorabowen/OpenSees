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


@pytest.mark.parametrize("engine", ["python", "tcl"])
@pytest.mark.parametrize("name", REFUSALS)
def test_mode_a_refusal(engine, name):
    check(engine, name)


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
