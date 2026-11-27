"""WP-261: the whole LadrunoIMKBeam battery through the safeguarded hinge solve.

The plain local Newton of ``ladrunoIMKSolveAxis`` / ``ladrunoIMKSolveBiaxial``
solves almost every hinge state, so the fallback that takes over when it fails
(``ladrunoIMKSafeguardedSolve`` with the elements' own residual and Jacobian
lambdas) is reached only by hard states: unloading from a plastic branch. The
dedicated tests reach it on the strong-axis, two-hinge path only.

``LADRUNO_IMK_HINGE_FALLBACK_ONLY=1`` (read once per process) skips the plain
loop, so EVERY hinge solve goes through the fallback. This test reruns the IMK
test files in a child process with the switch set. Their oracles must still
hold: single hinges, released ends, series shear hinges, the coupled biaxial
NDMaterial path, weak-axis hinges, end zones, P-M scaling, IMPL-EX, Rayleigh,
hinge-status responses. A sign or index error in any of the fallback's 3x3 or
6x6 assemblies fails here.

Excluded: ``test_ladruno_imk_beam_3d_simultaneous_biaxial_and_shear_hinge_softening``
(fails with or without WP-261, issue #1169).
"""
import os
import subprocess

import pytest

from _testbed import ops  # noqa: F401  (pins the engine the child must use)
from _testbed.subprocess_run import pinned_child

pytestmark = [pytest.mark.zone_a]

_HERE = os.path.dirname(os.path.abspath(__file__))

_FILES = [
    "test_ladrunoIMKBeam2d_element.py",
    "test_ladrunoIMKBeam_element.py",
    "test_ladrunoIMKBeam_unload.py",
    "test_ladruno_imk_endzone_loads.py",
    "test_imk_betakc.py",
    "test_rayleigh_inertia_bezier_imk.py",
    "test_ladruno_hinge_status.py",
    "test_ladruno_beam_consistency.py",
    "test_ladruno_response_tokens.py",
    "test_partial_trapezoidal_beam_loads.py",
    "test_ladruno_force_beam_column.py",
]
_DESELECT = [
    "test_ladruno_force_beam_column.py::"
    "test_ladruno_imk_beam_3d_simultaneous_biaxial_and_shear_hinge_softening",   # #1169
]


def test_imk_battery_passes_with_fallback_only():
    env = dict(os.environ)
    env["LADRUNO_IMK_HINGE_FALLBACK_ONLY"] = "1"
    args = ["-m", "pytest", "-q", "-p", "no:cacheprovider", "-x"]
    for d in _DESELECT:
        args += ["--deselect", d]
    args += _FILES
    argv, env = pinned_child(args, env)
    p = subprocess.run(argv, env=env, cwd=_HERE, stdin=subprocess.DEVNULL,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=1200)
    tail = "\n".join((p.stdout + p.stderr).splitlines()[-40:])
    assert p.returncode == 0, tail
    assert " passed" in p.stdout, tail


def test_fallback_only_switch_routes_every_solve_through_the_fallback():
    """The child really ran the fallback: in switch mode a hinged element's
    first nonzero state is solved by it (hingeSolveStats[0] > 0)."""
    script = r'''
import sys
sys.path.insert(0, r"%s")
from _testbed import ops
ops.wipe(); ops.model("basic", "-ndm", 2, "-ndf", 3)
ops.node(1, 0.0, 0.0); ops.node(2, 3.0, 0.0)
ops.fix(1, 1, 1, 1); ops.fix(2, 1, 1, 0)
ops.geomTransf("Linear", 1)
ops.uniaxialMaterial("Steel01", 10, 1.0e5, 8.0e8, 1.0e-4)
ops.element("LadrunoIMKBeam", 1, 1, 2, 1.0e-2, 2.0e11, 2.0e-5, 1, "-hinge", "both", "-matZ", 10)
ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(2, 0.0, 0.0, 1.0)
ops.system("BandGeneral"); ops.numberer("Plain"); ops.constraints("Transformation")
ops.test("NormDispIncr", 1.0e-12, 25); ops.algorithm("Newton")
ops.integrator("DisplacementControl", 2, 3, 1.0e-3); ops.analysis("Static")
assert ops.analyze(3) == 0
print("STATS", *ops.eleResponse(1, "hingeSolveStats"))
''' % _HERE
    env = dict(os.environ)
    env["LADRUNO_IMK_HINGE_FALLBACK_ONLY"] = "1"
    argv, env = pinned_child(["-u", "-c", script], env)
    p = subprocess.run(argv, env=env, stdin=subprocess.DEVNULL, capture_output=True,
                       text=True, encoding="utf-8", errors="replace", timeout=300)
    assert p.returncode == 0, p.stdout + p.stderr
    line = [ln for ln in p.stdout.splitlines() if ln.startswith("STATS")][0]
    n_fb, n_ev, n_bad = (float(t) for t in line.split()[1:])
    assert n_fb >= 3 and n_ev >= n_fb and n_bad == 0, line
