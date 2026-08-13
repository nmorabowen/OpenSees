"""Dense SOEs killed the whole process on a model with ZERO free equations. CI gate.

A model whose every DOF is either fixed or sp-constrained (under `constraints
Transformation`) numbers ZERO equations. `system UmfPack` handled that (rc=0);
the classic dense SOEs did not:

  * FullGenLinSOE / BandGenLinSOE / BandSPDLinSOE / ProfileSPDLinSOE /
    SProfileSPDLinSOE all create their `vectX`/`vectB` (FullGen also `matA`)
    Vector wrappers inside `if (size != oldSize)` in setSize(). With zero
    equations `size == oldSize == 0`, the block never runs, the wrappers stay
    null from the default constructor, and the first `getB()`/`getX()` hits
    the `FATAL ... exit(-1)` path. DiagonalSOE excluded `size == 0` from that
    block *explicitly*, with the same consequence.
  * The two ProfileSPD variants additionally read `iDiagLoc[size-1]` ==
    `iDiagLoc[-1]` (UB) when an existing SOE is RESIZED down to zero equations
    with a stale `iDiagLoc`, and left `profileSize` stale.

SYMPTOM THAT MADE THIS EXPENSIVE TO BISECT: `exit(-1)` == whole-process death
with exit code 255 (Windows) and NO Python traceback; the FATAL message goes to
opserr, which the Python module redirects, so nothing prints. faulthandler is
silent because it is a clean exit(), not a signal. It looks exactly like a bug
in the test file. See LEDGER_quirks.md ("zero-free-equation models").

Each system therefore runs in a SUBPROCESS here: on a regression the child dies
with a nonzero exit code and the parent reports it as an ordinary test failure
instead of the whole pytest run being killed. Two Windows-specific traps in that
subprocess bootstrap are commented at their call sites (an oversized PYTHONPATH
and an inherited stdin handle) -- both produce failures that mimic the very
crash this file detects, and both appear only under the full suite.

Two scenarios per system, both on the 8-node unit stdBrick cube with
1/8-symmetry fixes and every remaining DOF driven by an sp (uniform-strain
patch field, so the exact solution is known):

  fresh : zero equations from the first setSize(). This is the case that CRASHES
          pre-fix -- 6 of the 7 systems die (UmfPack is the control).
  shrink: analyze with node 7 free (3 equations), then sp node 7's DOFs in a new
          pattern and analyze again -> the same SOE is resized 3 -> 0.

STATE PLAINLY WHAT THE SHRINK HALF DOES AND DOES NOT CATCH: it does **not**
reproduce the crash, and it passes on the pre-fix binary. A 3 -> 0 resize has
`size != oldSize`, so the wrapper block DOES run and the wrappers are rebuilt at
zero length -- the null-wrapper path is reachable only on the FIRST setSize.
What shrink guards is the second, quieter defect on the same route: the two
ProfileSPD variants' `profileSize = iDiagLoc[size-1]` reads `iDiagLoc[-1]` when
size is 0, which is UB that happens not to fault here. It is kept as a
non-crashing regression guard on a path with no other coverage, not as a
reproducer. Pre-fix score, measured: 6 of 13 fail, all of them `fresh`.

Upstream elements/material only (stdBrick + ElasticIsotropic) => zone_a.
"""
import json
import os
import subprocess
import sys

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# The child must load the SAME binary the parent resolved -- hand it that one
# directory explicitly. Do NOT forward the parent's whole `sys.path` as
# PYTHONPATH: under the full suite that string grows past Windows' 32 KB
# environment-block limit and every `CreateProcess` dies with
# `OSError: [WinError 50] The request is not supported` -- which looks exactly
# like the crash this file is meant to detect. Standalone runs pass, full-suite
# runs fail, on a change that has nothing to do with either.
ENGINE_DIR = os.path.dirname(os.path.abspath(ops.__file__))

EPS_X = 1.0e-3      # prescribed ux on the x=1 face
EPS_Y = -2.5e-4     # prescribed uy on the y=1 face
EPS_Z = 5.0e-4      # prescribed uz on the z=1 face

# UmfPack always worked -- kept as the control that the scenario itself is sane.
SYSTEMS = ["FullGeneral", "BandGeneral", "BandSPD", "ProfileSPD", "SProfileSPD",
           "UmfPack"]

# `Diagonal` has the same null-wrapper hole and is fixed with it, but it only
# ever solves the DIAGONAL of the tangent -- a Newton static on a stdBrick with
# real free equations is then Jacobi iteration and would not converge in the
# shrink scenario's phase A. It gates the zero-equation path only, where there
# is nothing to solve and the distinction is moot.
SYSTEMS_FRESH = SYSTEMS + ["Diagonal"]

CHILD = r'''
import json, os, sys

_D = %(ENGINE_DIR)r
if os.path.isdir(_D):
    os.environ["PATH"] = _D + os.pathsep + os.environ.get("PATH", "")
    try:
        os.add_dll_directory(_D)
    except (OSError, FileNotFoundError):
        pass
    sys.path.insert(0, _D)

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

SYS = sys.argv[1]
SCEN = sys.argv[2]
EPS_X, EPS_Y, EPS_Z = %(EPS_X)r, %(EPS_Y)r, %(EPS_Z)r

# unit cube; fixes are the 1/8-symmetry set (planes x=0, y=0, z=0)
NODES = {
    1: ((0.0, 0.0, 0.0), (1, 1, 1)),
    2: ((1.0, 0.0, 0.0), (0, 1, 1)),
    3: ((1.0, 1.0, 0.0), (0, 0, 1)),
    4: ((0.0, 1.0, 0.0), (1, 0, 1)),
    5: ((0.0, 0.0, 1.0), (1, 1, 0)),
    6: ((1.0, 0.0, 1.0), (0, 1, 0)),
    7: ((1.0, 1.0, 1.0), (0, 0, 0)),
    8: ((0.0, 1.0, 1.0), (1, 0, 0)),
}


def build(skip_node7_sps):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (xyz, fx) in NODES.items():
        ops.node(tag, *xyz)
        if any(fx):
            ops.fix(tag, *fx)
    ops.nDMaterial("ElasticIsotropic", 1, 200000.0, 0.25)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)

    # hold the factor at 1.0 past t=1 -- Path returns 0 outside its time range
    ops.timeSeries("Path", 1, "-time", 0.0, 1.0, 2.0, "-values", 0.0, 1.0, 1.0)
    ops.timeSeries("Path", 2, "-time", 0.0, 1.0, 2.0, "-values", 0.0, 1.0, 1.0)

    ops.pattern("Plain", 1, 1)
    for tag, (xyz, fx) in NODES.items():
        if skip_node7_sps and tag == 7:
            continue
        if xyz[0] == 1.0 and fx[0] == 0:
            ops.sp(tag, 1, EPS_X)

    ops.pattern("Plain", 2, 2)
    for tag, (xyz, fx) in NODES.items():
        if skip_node7_sps and tag == 7:
            continue
        if xyz[1] == 1.0 and fx[1] == 0:
            ops.sp(tag, 2, EPS_Y)
        if xyz[2] == 1.0 and fx[2] == 0:
            ops.sp(tag, 3, EPS_Z)

    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system(SYS)
    ops.test("NormDispIncr", 1.0e-10, 10, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.5)
    ops.analysis("Static")


out = {"system": SYS, "scenario": SCEN}

if SCEN == "fresh":
    # zero free equations from the very first setSize()
    build(skip_node7_sps=False)
    out["rc1"] = ops.analyze(2)                       # t: 0 -> 1
    out["u7"] = [ops.nodeDisp(7, d) for d in (1, 2, 3)]
    out["u1"] = [ops.nodeDisp(1, d) for d in (1, 2, 3)]
else:
    # node 7 free first (3 equations), then resize the SAME SOE down to zero
    build(skip_node7_sps=True)
    out["rc1"] = ops.analyze(2)                       # t: 0 -> 1
    out["u7_free"] = [ops.nodeDisp(7, d) for d in (1, 2, 3)]
    ops.timeSeries("Path", 3, "-time", 1.0, 2.0, "-values", 0.0, 1.0)
    ops.pattern("Plain", 3, 3)
    ops.sp(7, 1, EPS_X)
    ops.sp(7, 2, EPS_Y)
    ops.sp(7, 3, EPS_Z)
    out["rc2"] = ops.analyze(2)                       # t: 1 -> 2, size 3 -> 0
    out["u7"] = [ops.nodeDisp(7, d) for d in (1, 2, 3)]

out["engine"] = os.path.dirname(os.path.abspath(ops.__file__))
print("RESULT " + json.dumps(out))
'''


def _run_child(system, scenario):
    script = CHILD % {"EPS_X": EPS_X, "EPS_Y": EPS_Y, "EPS_Z": EPS_Z,
                      "ENGINE_DIR": ENGINE_DIR}
    proc = subprocess.run(
        [sys.executable, "-c", script, system, scenario],
        # stdin=DEVNULL is load-bearing on Windows: with the inherited stdin,
        # `Popen` tries to `DuplicateHandle` whatever pytest's capture left
        # there and raises `OSError: [WinError 50] The request is not
        # supported` -- but only once some OTHER module in the suite has run,
        # so the file passes standalone and every case fails under `pytest
        # tests/`. That failure is indistinguishable here from the crash this
        # file exists to detect.
        stdin=subprocess.DEVNULL,
        capture_output=True, text=True, timeout=300,
        encoding="utf-8", errors="replace",
    )
    assert proc.returncode == 0, (
        f"{system}/{scenario}: child process died (exit {proc.returncode}) -- "
        f"the zero-free-equation crash is back.\n"
        f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
    )
    for line in proc.stdout.splitlines():
        if line.startswith("RESULT "):
            out = json.loads(line[len("RESULT "):])
            # a child that quietly loaded a DIFFERENT binary would report on
            # something nobody built -- the classic stale-engine false pass
            assert os.path.normcase(out["engine"]) \
                == os.path.normcase(ENGINE_DIR), (
                f"{system}/{scenario}: child loaded {out['engine']}, parent is "
                f"{ENGINE_DIR}"
            )
            return out
    raise AssertionError(
        f"{system}/{scenario}: no RESULT line.\n"
        f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
    )


@pytest.mark.parametrize("system", SYSTEMS_FRESH)
def test_zero_free_equations_fresh(system):
    """All DOFs fixed or sp'd from the start: analyze must succeed, sps enforced."""
    out = _run_child(system, "fresh")
    assert out["rc1"] == 0, out
    ux, uy, uz = out["u7"]
    assert abs(ux - EPS_X) < 1e-12 and abs(uy - EPS_Y) < 1e-12 \
        and abs(uz - EPS_Z) < 1e-12, out
    assert max(abs(v) for v in out["u1"]) < 1e-15, out


@pytest.mark.parametrize("system", SYSTEMS)
def test_zero_free_equations_shrink(system):
    """Same SOE resized 3 -> 0 equations mid-run (stale iDiagLoc path)."""
    out = _run_child(system, "shrink")
    assert out["rc1"] == 0 and out["rc2"] == 0, out
    # phase A only needs to have produced a sane solve (the free corner is NOT
    # on the uniform field -- constant stress leaves a nonzero nodal force at a
    # free corner, so it relaxes); the gate is the resize-to-zero that follows
    assert all(abs(v) < 1.0 for v in out["u7_free"]), out
    # after the shrink the corner is prescribed exactly
    ux, uy, uz = out["u7"]
    assert abs(ux - EPS_X) < 1e-12 and abs(uy - EPS_Y) < 1e-12 \
        and abs(uz - EPS_Z) < 1e-12, out
