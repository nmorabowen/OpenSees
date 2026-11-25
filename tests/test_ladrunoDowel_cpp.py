"""Ladruno WP-257: the bar dowel kernel (SRC/material/nD/LadrunoDowelKernel.h) in g++, without OpenSees.

tests/_testbed/dowel_gpp.cpp compiles the header alone:
- self-checks D1-D5 (capacity and envelope closed forms, tangent vs central differences, the
  origin-oriented path, the axial interaction);
- a cyclic (gamma, sigma) path with commits against an independent numpy transcription of the law
  (dowel_ref below), to round-off.

Skipped without a C++ compiler.
"""
import math
import os
import shutil
import subprocess
import tempfile

import numpy as np
import pytest
from _cxx import cxx_runtime_flags  # WP-199: static GCC runtime on Windows

pytestmark = [pytest.mark.zone_a, pytest.mark.t1]

_HERE = os.path.dirname(os.path.abspath(__file__))
_ND = os.path.normpath(os.path.join(_HERE, "..", "SRC", "material", "nD"))
_GATE = os.path.join(_HERE, "_testbed", "dowel_gpp.cpp")


def dowel_ref(fc, fy, db, ls, kd, su, rows, cd=0.0):
    """tau per row of (gamma, sigma, commit): Dulacska / Vintzeleou-Tassios capacity, exponential
    envelope to 95 % at su, secant to the origin inside the largest committed slip; with cd > 0 the
    cyclic decay 1/(1 + cd D), D the cumulative re-slip below smax in su (implicit: the trial adds
    its own re-slip since the last commit)."""
    su = su if su > 0 else 0.05 * db
    s0 = su / 3.0
    tu = kd * 4.0 / math.pi * math.sqrt(fc * fy)
    smax, sc, D, out = 0.0, 0.0, 0.0, []
    for g, sig, c in rows:
        a = abs(sig) / fy
        cap = tu * math.sqrt(1.0 - a * a) if a < 1.0 else 0.0
        s = abs(g) * ls
        h = 1.0 - math.exp(-s / s0) if s >= smax else (1.0 - math.exp(-smax / s0)) * s / smax
        up = min(s, smax) - sc
        Dt = D + up / su if up > 0 else D
        b = 1.0 / (1.0 + cd * Dt) if cd > 0 else 1.0
        out.append(math.copysign(1.0, g) * cap * h * b if g != 0 else 0.0)
        if c:
            D, sc, smax = Dt, s, max(smax, s)
    return np.array(out)


@pytest.fixture(scope="module")
def exe():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    td = tempfile.mkdtemp()
    out = os.path.join(td, "dowel.exe")
    cc = subprocess.run([gpp, "-O2", "-std=c++17", "-I", _ND, _GATE, "-o", out, *cxx_runtime_flags()],
                        stdin=subprocess.DEVNULL, capture_output=True, text=True, check=False)
    assert cc.returncode == 0, f"g++ failed to compile the dowel gate:\n{cc.stderr}"
    yield out
    shutil.rmtree(td, ignore_errors=True)


def test_self_checks(exe):
    r = subprocess.run([exe], stdin=subprocess.DEVNULL, capture_output=True, text=True, check=False)
    assert r.returncode == 0, f"g++ dowel gate FAILED:\n{r.stdout}\n{r.stderr}"
    for head in ("D1 capacity", "D2 envelope", "D3 tangent", "D4 origin-oriented path", "D5 axial interaction", "D6 cyclic decay",
                 "GPP DOWEL GATE: PASS"):
        assert head in r.stdout, r.stdout


@pytest.mark.parametrize("su, cd", [(-1.0, 0.0), (2.0, 0.0), (-1.0, 0.3), (2.0, 2.0)],
                         ids=["su_default", "su_2mm", "decay", "decay_su_2mm"])
def test_cyclic_path_vs_reference(exe, su, cd):
    par = (33.0, 331.0, 25.4, 150.0, 0.4, su, cd)
    g = np.concatenate([np.linspace(0, 4e-3, 9), np.linspace(4e-3, -2e-3, 13), np.linspace(-2e-3, 6e-3, 17)])
    sig = 300.0 * np.sin(np.arange(g.size) / 5.0)             # the bar's axial stress wanders, |sigma| < fy
    sig[-3:] = [331.0, -400.0, 0.0]                          # yielded both ways, then elastic again
    rows = [(float(a), float(b), int(i % 3 != 1)) for i, (a, b) in enumerate(zip(g, sig))]   # some trials uncommitted
    text = " ".join(repr(float(x)) for x in par) + f"\n{len(rows)}\n"
    text += "".join(f"{a!r} {b!r} {c}\n" for a, b, c in rows)
    r = subprocess.run([exe, "path"], input=text, capture_output=True, text=True, check=False)
    assert r.returncode == 0, r.stdout + r.stderr
    got = np.array([[float(x) for x in ln.split()] for ln in r.stdout.strip().splitlines()])
    exp = dowel_ref(*par[:6], rows, cd=cd)
    assert np.abs(got[:, 0] - exp).max() < 1e-12 * np.abs(exp).max()
    assert got[-3, 0] == 0.0 and got[-2, 0] == 0.0           # yielded: no dowel
