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


def test_reseed_bar_has_no_crossing():
    """The re-seed protocols are only meaningful if no support changes in the run."""
    r = run("python", "reseed_bar_normal")
    assert r["info_connectivityChanges"] == 1, r["info_connectivityChanges"]   # the first step builds it


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
    """Freeze decisions 2 and 9, DATA not gates: the cantilever QS tip at massCutoff 1e-6
    (mode-A default) and 1e-10, exclude against ghost (gamma 1e-6), and Newton iterations per
    step (the hypo tangent is not the exact linearisation). Printed for the handoff."""
    ops = _ops()
    c = H.case_by_name("cantilever_qs_r1")
    rows = []
    for label, mc, en in (("exclude, c_m 1e-6", None, None), ("exclude, c_m 1e-10", 1e-10, None),
                          ("ghost 1e-6, c_m 1e-6", None, ("ghost", 1e-6))):
        r = H.run_python(c, ops, mass_cutoff=mc, empty_nodes=en)
        rows.append((label, r["hist_tip"][9.875][-1], r["stopped"]))
    print("\ncantilever_qs_r1 tip_disp_y:")
    for lab, tip, st in rows:
        print(f"  {lab:24s} {tip:.6f}  stopped={st}")
    # Newton iterations per step
    for name, args in H.model_commands(c):
        getattr(ops, name)(*H._py_args(args))
    for name, args in H.analysis_commands(c):
        getattr(ops, name)(*H._py_args(args))
    its = []
    for _ in range(c["steps"]):
        ops.mpm("analyze", 1, c["dt"])
        its.append(int(ops.testIter()))
    print("  Newton iterations per step:", its)
    assert all(st is None for _, _, st in rows)
