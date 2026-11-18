"""WP-232k S1, S9 (the Zone-A soil oracle) and S10 (engines): the one-point identity of soil
particles against 232k-V's one-point FE drivers (tests/data/mpm/soil/s1_onepoint.json; the
reference blocks of tests/mpm_soil_cases.txt). Every particle of a homogeneous block sees the
record's per-step dF through an affine grid motion (`mpm gridvel ... affine`); sigma and the
material state are compared per step with GP 1. Tolerance (freeze): rel 1e-8 before the FE's
own Gauss-point spread exceeds 1e-8, rel 1e-3 from there.

S9: the plane-strain triaxial legs are in Zone-A; they go red if the stage flip or the material
update is skipped at particles (the stage-0 operator differs from the stage-1 one from step 11).
The full 12 legs x 2 modes are marked slow (measured: about 20 s in all on the Windows box,
build 7c932089c; the SANISAND 3D shear legs are the longest at 7 s).
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


S1_3D_A_MISS = pytest.mark.xfail(strict=True, reason=(
    "owner-accepted 2026-11-18 recorded miss (ruling at 926d27153; build 7c932089c): SANISAND 3D in"
    " mode A at step 12, the"
    " second plastic step, one step before the FE's own Gauss-point spread passes 1e-8 (step 13):"
    " triax rel 3.2e-6, oedo 1.06e-8, with the FE's substep counts (5, 4). Mode A's kinematics differ"
    " from mode F's by 1.3e-12 in epsHat (the SP displacement path; mode F 7e-17 between"
    " particles) and ModifiedEuler's adaptive substep sizes amplify it ~1e6; mode F holds 3.4e-10"
    " and 1.1e-12 there. Not a kinematic error: the band presumes round-off-free strain input."))


def _s1_param(leg, mode):
    marks = [S1_3D_A_MISS] if (mode == "A" and leg in ("sanisand_3d_triax", "sanisand_3d_oedo")) else []
    return pytest.param(leg, mode, marks=marks, id=f"{leg}-{mode}")


@pytest.mark.slow
@pytest.mark.parametrize("leg,mode", [_s1_param(lg, m) for lg in LEGS for m in ("A", "F")])
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
