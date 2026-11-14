"""WP-232c: the explicit MPM engine through the shipped binaries, against 232b's
frozen cases (tests/mpm_cases.txt; ADR 232 §8, the S1 rows the cases cover).

Every case runs in BOTH shipped engines — Python (opensees.pyd) and classic Tcl
(OpenSees.exe) — through the `mpm` command, with nDMaterial ElasticIsotropic,
which must reproduce the interface's hypoelastic law (§9). The expectations are
the case file's own: closed forms (translation, bar period, incline, K0 quiet
start, PIC < FLIP energy), the oracle's particle fields after 1, 3 and 5 steps,
the unit-level hooks (shape, critical dt, active nodes: interface Q4), and the
critical-dt refusal. A last test holds the two engines to each other.

This is the MPM family's mutation suite (Ladruno_scripts/mutation_gate.py):
ZERO drops the P2G internal force (the bar, the particle fields and the K0 start
notice), SCALE scales the G2P position increment (translation, incline,
particle fields notice).

Bootstrap: dist/bin/opensees.pyd on PYTHONPATH (Ladruno_internal/BUILD_GOTCHAS.md
§4); OpenSees.exe from dist/bin or LADRUNO_TCL_EXE.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_mpm_harness as H  # noqa: E402

pytestmark = pytest.mark.zone_a

REPO = Path(__file__).resolve().parents[1]

# Interface A4: the 20%-compression self-weight column misses the frozen 0.03 in
# the oracle too (0.072); strict xfail until the orchestrator rules.
KNOWN_XFAIL = {("column_selfweight_2d_gimp", "selfweight_error")}

CASES = H.parse_cases()
NAMES = [c["name"] for c in CASES]


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


def _check(engine, name):
    c = H.case_by_name(name)
    r = run(engine, name)
    results = H.evaluate(c, r)
    assert results, f"{name}: no evaluable expectation"
    failed, xfailed = [], []
    for label, ok, detail in results:
        obs = label.split()[1]
        if (name, obs) in KNOWN_XFAIL:
            xfailed.append((label, ok, detail))
            continue
        if not ok:
            failed.append(f"{label}: {detail}")
    assert not failed, "\n".join(failed)
    return xfailed


@pytest.mark.parametrize("engine", ["python", "tcl"])
@pytest.mark.parametrize("name", NAMES)
def test_mpm_case(engine, name):
    _check(engine, name)


@pytest.mark.xfail(strict=True, reason="interface A4: uGIMP self-weight column at 20% "
                   "compression misses the frozen 0.03 (oracle 0.072)")
@pytest.mark.parametrize("engine", ["python", "tcl"])
def test_mpm_selfweight_column_gate(engine):
    xf = _check(engine, "column_selfweight_2d_gimp")
    assert xf and all(ok for _, ok, _ in xf), xf


@pytest.mark.parametrize("engine", ["python", "tcl"])
def test_mpm_case_relations(engine):
    """The case file's `compare` lines (PIC dissipates more than FLIP; the linear
    hat's cell-crossing meter exceeds uGIMP's)."""
    for c in CASES:
        for cp in c["compares"]:
            obs, args, op, other = cp[1], cp[2:-2], cp[-2], cp[-1]
            a = H.scalar(c, run(engine, c["name"]), obs, args)
            oc = H.case_by_name(other)
            b = H.scalar(oc, run(engine, other), obs, args)
            assert (a < b) if op == "lt" else (a > b), f"{c['name']} {obs} {a} {op} {other} {b}"


@pytest.mark.parametrize("name", ["onestep_2d_gimp", "fivestep_2d_gimp", "threestep_2d_linear",
                                  "onestep_3d_gimp", "fivestep_3d_gimp", "bar_flip_2d",
                                  "impact_discs_2d"])
def test_mpm_engine_parity(name):
    """Classic Tcl and Python give the same particle output (ADR 232 §8)."""
    rp, rt = run("python", name), run("tcl", name)
    for f in H.FIELDS_END:
        a, b = rp["end_" + f], rt["end_" + f]
        assert len(a) == len(b), f
        worst = max((abs(x - y) / (1e-15 + 1e-12 * abs(x)) for x, y in zip(a, b)), default=0.0)
        assert worst <= 1.0, f"{name} {f}: Python and Tcl differ (worst {worst:.3g} x tol)"
