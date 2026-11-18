"""WP-232k K4 (D7-D11): the routes into the one stage walk, observed through what only the walk
does -- SANISAND's Elastic2Plastic at EVERY particle. On a low-g column (l0 2, g 0.25) the
particles near the top sit at p < Pmin at the flip (mode A gimp: the top cell's four; mode F
bspline2: the top row's two); a route that reached one particle alone (the
MaterialStageParameter first-answer stop, MaterialStageParameter.cpp:79) would leave the others
to the lazy flip, which skips the reset. Every route must reset every one of them to
sigma = -Pmin I.

Routes: `mpm stage` (both modes); `parameter ... loadPattern ... materialState|updateMaterialStage`
+ `updateParameter` (mode F); mode A's `updateMaterialStage`; and the same in classic Tcl.
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
PMIN = 0.101


def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    if not hasattr(ops, "mpm"):
        pytest.skip("this opensees build has no `mpm` command")
    return ops


def _staged(ops, mode):
    cmd = S.load("s4_lowg.json")["decks"]["sanisand"]["cmd"]
    S.column_build(ops, cmd, "bspline2" if mode == "F" else "gimp", mode, l0=2.0, g=0.25)
    if mode == "A":
        ops.mpm("analyze", 10, 0.1)
    else:
        ops.mpm("analyze", 3000, 5e-4)


def _low(sig, szz):
    """the particles with p < Pmin (the ones Elastic2Plastic resets)"""
    return [k for k in range(len(szz)) if -(sig[3 * k] + sig[3 * k + 1] + szz[k]) / 3.0 < PMIN]


def _count(low, sig, szz):
    return sum(1 for k in low if sig[3 * k] == -PMIN and sig[3 * k + 1] == -PMIN and szz[k] == -PMIN)


@pytest.mark.parametrize("route", ["mpm_stage_A", "mpm_stage_F", "parameter_F",
                                   "parameter_updateMaterialStage_F", "updateMaterialStage_A"])
def test_route_resets_every_particle(route):
    ops = _ops()
    mode = route[-1]
    _staged(ops, mode)
    low = _low(S._l(ops.mpm("particleData", "sig")), S._l(ops.mpm("particleData", "sigZZ")))
    if route.startswith("mpm_stage"):
        ops.mpm("stage", 1, 1)
    elif route == "parameter_F":
        ops.parameter(7001, "loadPattern", 232, "materialState", 1)
        ops.updateParameter(7001, 1.0)
    elif route == "parameter_updateMaterialStage_F":
        ops.parameter(7002, "loadPattern", 232, "updateMaterialStage", 1)
        ops.updateParameter(7002, 1.0)       # a DOUBLE: the host reads theDouble, not theInt
    else:
        ops.updateMaterialStage("-material", 1, "-stage", 1)
    n = _count(low, S._l(ops.mpm("particleData", "sig")), S._l(ops.mpm("particleData", "sigZZ")))
    assert len(low) >= 2 and n == len(low), (route, len(low), n)


def _tcl_exe():
    exe = os.environ.get("LADRUNO_TCL_EXE") or str(REPO / "dist" / "bin" / "OpenSees.exe")
    return exe if Path(exe).exists() else None


@pytest.mark.parametrize("route", ["mpm_stage", "updateMaterialStage", "parameter"])
def test_tcl_routes(route):
    """The routes in classic Tcl (mode A, gimp): every top-cell particle reset."""
    exe = _tcl_exe()
    if exe is None:
        pytest.skip("no OpenSees.exe")
    cmd = S.load("s4_lowg.json")["decks"]["sanisand"]["cmd"]
    lines = ["model basic -ndm 2 -ndf 2", cmd,
             "mpm grid 232 -origin 0.0 0.0 -h 0.5 -cells 2 5 -basis gimp",
             "mpm mode quasiStatic",
             "mpm body 1 box 0.0 0.0 1.0 2.0 -ppc 2 -mat 1 -rho 2.0",
             "mpm fix xmin x", "mpm fix xmax x", "mpm fix ymin y", "mpm gravity 0.0 -0.25",
             "timeSeries Path 9001 -time 0.0 1.0 1.0e9 -values 0.0 1.0 1.0",
             "mpm timeSeries 9001",
             "constraints Plain", "numberer RCM", "system UmfPack",
             "test NormDispIncr 1.0e-10 25", "algorithm Newton", "integrator LoadControl 0.1",
             "analysis Static", "mpm stage 1 0", "mpm analyze 10 0.1",
             'puts "PSIG [mpm particleData sig]"', 'puts "PSZZ [mpm particleData sigZZ]"']
    lines += {"mpm_stage": ["mpm stage 1 1"],
              "updateMaterialStage": ["updateMaterialStage -material 1 -stage 1"],
              "parameter": ["parameter 7001 loadPattern 232 materialState 1",
                            "updateParameter 7001 1.0"]}[route]
    lines += ['puts "SIG [mpm particleData sig]"', 'puts "SZZ [mpm particleData sigZZ]"',
              'puts "X [mpm particleData x]"']
    with tempfile.TemporaryDirectory() as td:
        f = Path(td) / "stage.tcl"
        f.write_text("\n".join(lines) + "\n", encoding="utf-8")
        r = subprocess.run([exe, str(f)], stdin=subprocess.DEVNULL, capture_output=True, text=True,
                           cwd=td, timeout=300, env=dict(os.environ, LADRUNO_OPENSEES_QUIET="1"))
    vals = {}
    for line in (r.stdout + r.stderr).splitlines():
        t = line.split()
        if t and t[0] in ("SIG", "SZZ", "X", "PSIG", "PSZZ"):
            vals[t[0]] = [float(v) for v in t[1:]]
    assert set(vals) == {"SIG", "SZZ", "X", "PSIG", "PSZZ"}, (r.stdout + r.stderr)[-2000:]
    low = _low(vals["PSIG"], vals["PSZZ"])
    n = _count(low, vals["SIG"], vals["SZZ"])
    assert len(low) >= 2 and n == len(low), (route, len(low), n)
