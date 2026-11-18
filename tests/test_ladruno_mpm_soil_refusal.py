"""WP-232k S5 (the ADR §8 P5 row; interface D12-D14, amended at the freeze).

SANISAND, the W17 abort: `-implex -maxSubsteps 1`, isotropic stage 0 to ~100 kPa, the flip, one
simple-shear step dgamma = 1e-3. The IMPL-EX companion hits the substep cap in commitState and
declares the refusal (LadrunoSANISAND.cpp:4623-4683; 232k-V confirmed the trigger on the one-point
FE, tests/data/mpm/soil/d12_trigger.json): the run aborts, the particles keep the last accepted
step bit for bit, the message names the material and its tag, and a later `mpm analyze` is
refused until `mpm reset`. Both modes, classic Tcl and Python.

HS has no reachable commit refusal at particles (its trial refusal is forwarded): one simple-
shear step gamma = 10 trips |deps|^2 > 25 (LadrunoHSKernel.h:1721). (gamma = 6, D14's number, does
not: the Hughes-Winget increment de-rotated by R_mid has |deps|^2 = 20.7 < 25 -- measured on the
mode-F leg, which then integrated the step.) The step fails loudly, the
message names the material and tag, the particles keep the last accepted step, and there is NO
abort: a later step is admitted.
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_mpm_soil_s_harness as S  # noqa: E402

pytestmark = pytest.mark.zone_a
REPO = Path(__file__).resolve().parents[1]

SAN = ("LadrunoSANISAND 1 264.32 0.3129 0.6944 1.3309 0.71 0.027 0.83 0.45 101.0 0.005 1.3 0.968 "
       "3.5 0.05 5.75 12.5 1100.0 2.0 -implex -maxSubsteps 1")
HS = ("LadrunoHS 1 35000.0 35000.0 105000.0 0.5 10.0 38.0 8.0 -nu 0.2 -pref 100.0 -Rf 0.9 "
      "-OCR 1.5 -stage 0")


def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    if not hasattr(ops, "mpm"):
        pytest.skip("this opensees build has no `mpm` command")
    return ops


def _commands(mat, mode, dt):
    """A 2 x 2-cell plane-strain block in a wide grid (gamma = 10 moves particles +-3.75 m),
    every node with mass on one affine gridvel."""
    t = mat.split()
    args = [t[0], 1] + [S._num(x) for x in t[2:]]
    cmds = [("wipe", []), ("model", ["basic", "-ndm", 2, "-ndf", 2]), ("nDMaterial", args),
            ("mpm", ["grid", 232, "-origin", -4.0, -1.0, "-h", 0.5, "-cells", 18, 6]),
            ("mpm", ["mode", "quasiStatic" if mode == "A" else "explicit"]),
            ("mpm", ["body", 1, "box", 0.0, 0.0, 1.0, 1.0, "-ppc", 2, "-mat", 1, "-rho", 2.0]),
            ("mpm", ["massCutoff", 0.0]),
            ("mpm", ["gridvel", "all", "box", -4.6, -1.6, 5.6, 2.6, "affine", 0.0, 0.0, 0.0, 0.0,
                     "ref", 0.5, 0.5])]
    if mode == "A":
        cmds += [("constraints", ["Transformation"]), ("numberer", ["RCM"]),
                 ("system", ["UmfPack"]), ("test", ["NormDispIncr", 1e-12, 10]),
                 ("algorithm", ["Newton"]), ("integrator", ["LoadControl", dt]),
                 ("analysis", ["Static"]), ("mpm", ["cut", 2.0, 0])]
    cmds.append(("mpm", ["stage", 1, 0]))
    return cmds


def _confine_steps(dt, a=1.5e-4, n=10):
    """Isotropic plane-strain compression in n steps of strain a/n each, as L = d/dt."""
    d = -a / n
    return [[d / dt, 0.0, 0.0, d / dt]] * n


def _drive(ops, mat, mode):
    dt = 1.0 if mode == "A" else 1e-4
    for name, args in _commands(mat, mode, dt):
        getattr(ops, name)(*args)
    for L in _confine_steps(dt):
        ops.mpm("gridvel", "all", "setAffine", *L)
        ops.mpm("analyze", 1, dt)
    ops.mpm("stage", 1, 1)
    return dt


def _state(ops):
    return {f: S._l(ops.mpm("particleData", f)) for f in ("x", "F", "sig", "epsHat")}


@pytest.mark.parametrize("mode", ["A", "F"])
def test_s5_sanisand_commit_refusal_aborts(mode, capfd):
    ops = _ops()
    dt = _drive(ops, SAN, mode)
    before = _state(ops)
    ops.mpm("gridvel", "all", "setAffine", 0.0, 1e-3 / dt, 0.0, 0.0)
    capfd.readouterr()
    with pytest.raises(Exception):
        ops.mpm("analyze", 1, dt)
    cap = capfd.readouterr()
    err = cap.err + cap.out
    assert int(ops.mpm("debug", "aborted")) == 1
    assert "LadrunoSANISAND" in err and "tag 1" in err, err[-2000:]
    assert _state(ops) == before                      # unmoved, bit for bit
    with pytest.raises(Exception):
        ops.mpm("analyze", 1, dt)                      # refused until `mpm reset`


@pytest.mark.parametrize("mode", ["A", "F"])
def test_s5_hs_trial_refusal_fails_loudly_without_abort(mode, capfd):
    ops = _ops()
    dt = _drive(ops, HS, mode)
    before = _state(ops)
    ops.mpm("gridvel", "all", "setAffine", 0.0, 10.0 / dt, 0.0, 0.0)
    capfd.readouterr()
    with pytest.raises(Exception):
        ops.mpm("analyze", 1, dt)
    cap = capfd.readouterr()
    err = cap.err + cap.out
    assert "LadrunoHS" in err and "tag 1" in err, err[-2000:]
    assert int(ops.mpm("debug", "aborted")) == 0
    assert _state(ops) == before
    ops.mpm("gridvel", "all", "setAffine", 0.0, 1e-4 / dt, 0.0, 0.0)
    assert int(ops.mpm("analyze", 1, dt)) == 0          # no abort: the next step runs


def test_s5_sanisand_abort_classic_tcl():
    exe = os.environ.get("LADRUNO_TCL_EXE") or str(REPO / "dist" / "bin" / "OpenSees.exe")
    if not Path(exe).exists():
        pytest.skip("no OpenSees.exe")
    lines = []
    for name, args in _commands(SAN, "A", 1.0):
        lines.append(" ".join([name] + [repr(a) if isinstance(a, float) else str(a) for a in args]))
    for L in _confine_steps(1.0):
        lines.append("mpm gridvel all setAffine " + " ".join(repr(v) for v in L))
        lines.append("mpm analyze 1 1.0")
    lines.append("mpm stage 1 1")
    lines.append("mpm gridvel all setAffine 0.0 0.001 0.0 0.0")
    lines.append('if {[catch {mpm analyze 1 1.0}]} {puts "REFUSED"}')
    lines.append('puts "ABORTED [mpm debug aborted]"')
    with tempfile.TemporaryDirectory() as td:
        f = Path(td) / "s5.tcl"
        f.write_text("\n".join(lines) + "\n", encoding="utf-8")
        r = subprocess.run([exe, str(f)], capture_output=True, text=True, cwd=td, timeout=300,
                           env=dict(os.environ, LADRUNO_OPENSEES_QUIET="1"))
    out = r.stdout + r.stderr
    assert "REFUSED" in out and "ABORTED 1" in out, out[-2000:]
    assert "LadrunoSANISAND" in out, out[-2000:]
