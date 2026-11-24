"""Ladruno WP-254: the Model C MCFT strip kernel against the numpy reference, without OpenSees.

tests/_testbed/mcft_strip_gpp.cpp compiles SRC/material/nD/LadrunoMCFTKernel.h (and the WP-245 shear
map) alone. Three gates:

- the g++ self-checks K1-K6 (lazy root scan = full scan bit for bit, steel law, path independence,
  transverse equilibrium, tangent, root rule);
- strip paths: monotonic and cyclic (eps_x, gamma_xy) paths with commits, a core strip (light ties,
  rho fyt < fcr, so the cracking fold is crossed) and a cover strip (no ties), against
  tests/_testbed/mcft_strip_ref.py, to round-off;
- section paths: a Model C column section (40 strips, 8 bars, parabolic gamma, V = sum tau A) through
  the shear map, on a cracking cyclic (eps0, kappa, gamma) path with commits, against the reference
  Section.forces, to round-off.

Skipped without a C++ compiler.
"""
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pytest
from _cxx import cxx_runtime_flags  # WP-199: static GCC runtime on Windows

pytestmark = [pytest.mark.zone_a, pytest.mark.t1]

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.normpath(os.path.join(_HERE, "..", "SRC", "material"))
_GATE = os.path.join(_HERE, "_testbed", "mcft_strip_gpp.cpp")
sys.path.insert(0, os.path.join(_HERE, "_testbed"))
import mcft_strip_ref as ref

# round-off: stresses are compared to the strip's own scale (f'c), forces to the section's
TOL_STRESS = 1e-9
TOL_EY = 1e-12


@pytest.fixture(scope="module")
def exe():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    td = tempfile.mkdtemp()
    out = os.path.join(td, "mcft_strip.exe")
    cc = subprocess.run([gpp, "-O2", "-std=c++17", "-I", os.path.join(_SRC, "nD"), "-I", os.path.join(_SRC, "section"),
                         _GATE, "-o", out, *cxx_runtime_flags()],
                        stdin=subprocess.DEVNULL, capture_output=True, text=True, check=False)
    assert cc.returncode == 0, f"g++ failed to compile the MCFT strip gate:\n{cc.stderr}"
    yield out
    shutil.rmtree(td, ignore_errors=True)


def _f(x):
    return repr(float(x))


def _run(exe, mode, text):
    r = subprocess.run([exe, mode], input=text, capture_output=True, text=True, check=False)
    assert r.returncode == 0, f"gate exe failed ({mode}):\n{r.stdout}\n{r.stderr}"
    return np.array([[float(x) for x in ln.split()] for ln in r.stdout.strip().splitlines()])


def test_self_checks(exe):
    r = subprocess.run([exe], stdin=subprocess.DEVNULL, capture_output=True, text=True, check=False)
    assert r.returncode == 0, f"g++ MCFT strip gate FAILED:\n{r.stdout}\n{r.stderr}"
    for head in ("K1 lazy scan = full scan", "K2 steel law", "K3 path independence", "K4 transverse equilibrium",
                 "K5 tangent", "K6 root rule", "GPP MCFT STRIP GATE: PASS"):
        assert head in r.stdout, r.stdout


# ------------------------------------------------------------------------------------------- strips
SEC, DATA = ref.column_section()


def _strip_params(i):
    c = SEC.conc
    return [c.fc, DATA["fp"][i], DATA["ep"][i], DATA["e20"][i], c.smx, c.smy, c.agg, c.k_beta, c.c_vci,
            SEC.rho[i], SEC.tie_steel.fy, SEC.tie_steel.fu, SEC.tie_steel.Es, SEC.tie_steel.esh, SEC.tie_steel.esu]


def _ref_strip(i):
    c = SEC.conc
    conc = ref.Concrete(c.fc, np.array([DATA["fp"][i]]), np.array([DATA["ep"][i]]), np.array([DATA["e20"][i]]),
                        c.smx, c.smy, c.agg, c.k_beta, c.c_vci)
    return ref.Strips(conc, np.array([SEC.rho[i]]), SEC.tie_steel)


def _paths():
    """(name, rows of (ex, g, commit)). Uncommitted rows are trial probes from the same committed state."""
    P = {}
    n = 120
    P["mono_shear_compression"] = [(-2e-4, 1.2e-4 * k, 1) for k in range(n)]
    P["mono_shear_tension"] = [(1e-5 * k, 8e-5 * k, 1) for k in range(n)]
    P["crush_and_unload"] = ([(-1.5e-4 * k, 2e-5 * k, 1) for k in range(80)]
                             + [(-0.012 + 1.5e-4 * k, 1.6e-3, 1) for k in range(80)])
    amp = [0.002, -0.002, 0.004, -0.004, 0.008, -0.008, 0.012, -0.012]
    rows, g = [], 0.0
    for a in amp:
        steps = 30
        g0 = g
        for k in range(1, steps + 1):
            g = g0 + (a - g0) * k / steps
            rows.append((-3e-4 + 2e-4 * np.sin(len(rows) / 9.0), g, 1))
            if len(rows) % 7 == 0:
                rows.append((rows[-1][0] + 3e-4, rows[-1][1] * 1.3, 0))      # an uncommitted probe
    P["cyclic_shear"] = rows
    rng = np.random.default_rng(254)
    ex, g, rows = 0.0, 0.0, []
    for _ in range(200):
        ex = float(np.clip(ex + rng.normal(0, 3e-4), -0.008, 0.01)); g = float(np.clip(g + rng.normal(0, 6e-4), -0.02, 0.02))
        rows.append((ex, g, int(rng.random() < 0.8)))
    P["random_walk"] = rows
    return P


CORE = int(np.nonzero(DATA["core"])[0][5])
COVER = 0


@pytest.mark.parametrize("which", ["core", "cover"])
@pytest.mark.parametrize("path", list(_paths()))
def test_strip_path_vs_reference(exe, which, path):
    i = CORE if which == "core" else COVER
    rows = _paths()[path]
    text = "STRIP " + " ".join(_f(v) for v in _strip_params(i)) + f"\n{len(rows)}\n"
    text += "".join(f"{_f(ex)} {_f(g)} {c}\n" for ex, g, c in rows)
    got = _run(exe, "strip", text)
    st_ = _ref_strip(i)
    st = st_.new_state()
    exp = []
    for ex, g, c in rows:
        o = st_.stresses(np.array([ex]), np.array([g]), st)
        exp.append([o["sx"][0], o["txy"][0], o["ey"][0]])
        if c:
            st_.commit(np.array([ex]), np.array([g]), st)
    exp = np.array(exp)
    scale = SEC.conc.fc
    ds = np.abs(got[:, :2] - exp[:, :2]).max() / scale
    dy = np.abs(got[:, 2] - exp[:, 2]).max()
    assert ds < TOL_STRESS and dy < TOL_EY, f"{which} {path}: max |ds|/fc = {ds:.3e}, max |d eps_y| = {dy:.3e}"
    if which == "core" and path in ("mono_shear_compression", "cyclic_shear"):
        # the path must actually crack and use the scan (otherwise the gate proves little)
        assert got[:, 3].max() > 10 * SEC.conc.ecr and got[:, 7].sum() > len(rows) // 2


# ------------------------------------------------------------------------------------------ section
def _section_text(rows):
    lines = [f"SECTION {_f(SEC.h)} {_f(SEC.b)} {SEC.nstrip} {len(SEC.bars_y)}"]
    for i in range(SEC.nstrip):
        lines.append(f"{_f(SEC.y[i])} 0.0 {_f(SEC.Ac[i])} " + " ".join(_f(v) for v in _strip_params(i)))
    bs = SEC.bar_steel
    for yb, A in zip(SEC.bars_y, SEC.bars_A):
        lines.append(" ".join(_f(v) for v in (yb, 0.0, A, bs.fy, bs.fu, bs.Es, bs.esh, bs.esu)))
    lines.append(str(len(rows)))
    lines += [f"{_f(e0)} {_f(k)} {_f(g)} {c}" for e0, k, g, c in rows]
    return "\n".join(lines) + "\n"


def _section_path():
    rows = []
    amp = [(4e-6, 0.002), (-4e-6, -0.002), (1e-5, 0.005), (-1e-5, -0.005), (2e-5, 0.01), (-2e-5, -0.01)]
    k0, g0 = 0.0, 0.0
    for ka, ga in amp:
        for j in range(1, 41):
            t = j / 40
            k, g = k0 + (ka - k0) * t, g0 + (ga - g0) * t
            rows.append((-3e-4 + 0.4 * abs(k) * 228.6, k, g, 1))
            if j % 9 == 0:
                rows.append((rows[-1][0], 1.1 * k, 1.1 * g, 0))
        k0, g0 = ka, ga
    return rows


def test_section_path_vs_reference(exe):
    rows = _section_path()
    got = _run(exe, "section", _section_text(rows))
    st = SEC.new_state()
    exp = []
    for e0, k, g, c in rows:
        exp.append(SEC.forces((e0, k, g), st))
        if c:
            SEC.commit((e0, k, g), st)
    exp = np.array(exp)
    scale = np.array([SEC.conc.fc * SEC.b * SEC.h, SEC.conc.fc * SEC.b * SEC.h ** 2, SEC.conc.fc * SEC.b * SEC.h])
    d = (np.abs(got - exp) / scale).max(axis=0)
    assert d.max() < TOL_STRESS, f"section N, M, V max rel. diff {d}"
    assert np.abs(exp[:, 2]).max() > 0.3 * np.abs(exp[:, 2]).max(initial=0) and np.abs(exp[:, 2]).max() > 1e5
