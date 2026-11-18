"""WP-232k S1, S9 (the Zone-A soil oracle) and S10 (engines): the one-point identity of soil
particles against 232k-V's one-point FE drivers (tests/data/mpm/soil/s1_onepoint.json; the
reference blocks of tests/mpm_soil_cases.txt). Every particle of a homogeneous block sees the
record's per-step dF through an affine grid motion (`mpm gridvel ... affine`); sigma and the
material state are compared per step with GP 1. Tolerance (freeze): rel 1e-8 before the FE's
own Gauss-point spread exceeds 1e-8, rel 1e-3 from there.

S9: the plane-strain triaxial legs are in Zone-A; they go red if the stage flip or the material
update is skipped at particles (the stage-0 operator differs from the stage-1 one from step 11).
The full 12 legs x 2 modes are marked slow (run by the warrant).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_mpm_soil_s_harness as S  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
LEGS = [f"{m}_{d}_{p}" for m in ("sanisand", "hs") for d in ("3d", "ps")
        for p in ("triax", "oedo", "shear")]


def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    if not hasattr(ops, "mpm"):
        pytest.skip("this opensees build has no `mpm` command")
    return ops


def _tcl_exe():
    exe = os.environ.get("LADRUNO_TCL_EXE")
    if exe:
        return exe if Path(exe).exists() else None
    cand = REPO / "dist" / "bin" / ("OpenSees.exe" if os.name == "nt" else "OpenSees")
    return str(cand) if cand.exists() else None


def _check(out):
    bad = [r for r in out["steps"] if not (r["err_sig"] <= r["tol"] and r["err_state"] <= r["tol"])]
    assert not bad, (f"{out['leg']} mode {out['mode']}: first failing step {bad[0]}"
                     f" (max err {out['max_err']:.3g}; first MPM/FE substep difference at"
                     f" step {out['first_substep_diff']})")


@pytest.mark.zone_a
@pytest.mark.parametrize("mode", ["A", "F"])
@pytest.mark.parametrize("leg", ["hs_ps_triax", "sanisand_ps_triax"])
def test_s9_onepoint_triax(leg, mode):
    _check(S.s1_run(_ops(), leg, mode))


@pytest.mark.slow
@pytest.mark.parametrize("mode", ["A", "F"])
@pytest.mark.parametrize("leg", LEGS)
def test_s1_onepoint_all(leg, mode):
    _check(S.s1_run(_ops(), leg, mode))


@pytest.mark.zone_a
def test_s10_tcl_equals_python():
    """S10: classic Tcl (OpenSees.exe) and Python give the same S1 triaxial answer, bit for bit
    in the printed doubles (the S1 leg hs_ps_triax, mode A, 20 steps across the flip)."""
    exe = _tcl_exe()
    if exe is None:
        pytest.skip("no OpenSees.exe")
    ops = _ops()
    rec = S.load("s1_onepoint.json")
    S.s1_run(ops, "hs_ps_triax", "A", rec=rec, nsteps=20)
    py = S._l(ops.mpm("particleData", "sig"))
    tcl = S.run_tcl(S.s1_tcl_deck(rec, "hs_ps_triax", nsteps=20), exe)
    assert len(py) == len(tcl)
    scale = max(abs(v) for v in py)
    assert max(abs(a - b) for a, b in zip(py, tcl)) <= 1e-12 * scale
