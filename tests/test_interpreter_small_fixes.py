"""Small interpreter fixes.

  eigen_*        repeated openseespy eigen calls, with and without an analysis

Python cases run in a child interpreter so that one failing case cannot take
down the run.
"""
import json
import os
import subprocess
import sys

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


ENGINE = os.path.abspath(ops.__file__)
REL_TOL = 1.0e-8

CHILD = r'''
import json, os, sys
try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

SCEN = sys.argv[1]
out = {"engine": os.path.abspath(ops.__file__)}


def chain():
    # truss chain with lateral springs of distinct stiffness: 20 free DOFs,
    # no repeated eigenvalues
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.uniaxialMaterial("Elastic", 1, 3000.0)
    ops.node(0, 0.0, 0.0)
    ops.fix(0, 1, 1)
    for i in range(1, 11):
        ops.node(i, 10.0 * i, 0.0)
        ops.mass(i, 5.0, 5.0)
        ops.element("Truss", i, i - 1, i, 10.0, 1)
        ops.node(100 + i, 10.0 * i, 10.0)
        ops.fix(100 + i, 1, 1)
        ops.element("Truss", 100 + i, 100 + i, i, 4.0 + 0.5 * i, 1)


def transient():
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("BandGen")
    ops.test("NormDispIncr", 1.0e-8, 10)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")


def eig(*args):
    try:
        return list(ops.eigen(*args))
    except Exception as ex:
        return "raised: " + repr(ex)


if SCEN == "eigen_same_type":
    chain()
    out["runs"] = [eig(3), eig(3), eig(3)]
elif SCEN == "eigen_explicit_type":
    chain()
    out["runs"] = [eig("-genBandArpack", 3), eig("-genBandArpack", 3)]
elif SCEN == "eigen_switched_type":
    chain()
    out["runs"] = [eig("-genBandArpack", 3), eig("-fullGenLapack", 3),
                   eig("-genBandArpack", 3)]
elif SCEN == "eigen_with_analysis":
    chain()
    transient()
    out["runs"] = [eig("-fullGenLapack", 3), eig("-fullGenLapack", 3)]
elif SCEN == "eigen_then_analysis":
    chain()
    first = eig(3)
    transient()
    out["runs"] = [first, eig(3), eig(3)]

print("RESULT " + json.dumps(out))
'''


def _run(scenario):
    proc = subprocess.run(
        [sys.executable, "-c", CHILD, scenario],
        stdin=subprocess.DEVNULL, capture_output=True, text=True,
        encoding="utf-8", errors="replace", timeout=300,
    )
    assert proc.returncode == 0, (
        f"{scenario}: child exited with {proc.returncode}\n"
        f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}")
    lines = [l for l in proc.stdout.splitlines() if l.startswith("RESULT ")]
    assert lines, f"{scenario}: no result\n{proc.stdout}\n{proc.stderr}"
    out = json.loads(lines[-1][len("RESULT "):])
    assert os.path.normcase(out["engine"]) == os.path.normcase(ENGINE), out
    return out


def _assert_same_modes(runs, tol=REL_TOL):
    for r in runs:
        assert isinstance(r, list) and len(r) == 3, runs
        assert all(x > 0.0 for x in r), runs
    for r in runs[1:]:
        for x, y in zip(runs[0], r):
            assert abs(x - y) <= tol * abs(x), runs


@pytest.mark.parametrize("scenario", ["eigen_same_type", "eigen_explicit_type",
                                      "eigen_with_analysis", "eigen_then_analysis"])
def test_repeated_eigen(scenario):
    _assert_same_modes(_run(scenario)["runs"])


def test_repeated_eigen_switched_type():
    runs = _run("eigen_switched_type")["runs"]
    _assert_same_modes(runs, tol=1.0e-6)
    _assert_same_modes([runs[0], runs[2]])
