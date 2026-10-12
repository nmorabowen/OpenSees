"""WP concrete3d-hang-diagnosis review M3: hosts that DROPPED setTrialStrain's return code.

LadrunoSolidShell (inside its enhanced-strain Newton) and LogStrainNDMaterial (around the inner small-strain material) called
setTrialStrain and ignored the result. A material that refuses a wild trial therefore never cut the step at the TRIAL: only
the WP-99 commit latch fired, Domain::commit() then aborted with the model INCONSISTENT (nodes and the sibling integration
points had already committed -- the analysis reports failure but the committed model is not the last converged state).

The discriminator used here is exactly that inconsistency. A refused step must (1) fail (rc < 0), (2) leave the committed
NODAL state untouched (nodal displacements equal the pre-step values -- the latch abort leaves the wild displacements
committed), and (3) be recoverable: the same model, stepped again with a SMALLER increment, converges. That is "the analysis
cuts the step" as opposed to "aborts".

Wild strain: the review's probe_absorb.cpp direction after a realistic 30 x 2e-5 crack phase (shared with
test_ladrunoConcrete3D_wild_trial_refusal.py; shear doubled to engineering). At full magnitude the trial is refused; at 5 %
it converges.
"""
import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_NODES = {
    1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0), 3: (1.0, 1.0, 0.0), 4: (0.0, 1.0, 0.0),
    5: (0.0, 0.0, 1.0), 6: (1.0, 0.0, 1.0), 7: (1.0, 1.0, 1.0), 8: (0.0, 1.0, 1.0),
}
_CONN = [1, 2, 3, 4, 5, 6, 7, 8]
_CRACK_DIR = (1.0, -0.2, -0.2, 0.6, 0.0, 0.0)
_CRACK_STEP, _CRACK_STEPS = 2.0e-5, 30
_WILD_RAW = (0.174, -0.502, 0.956, 0.964, -0.732, -1.448)


def _affine_disp(x, y, z, e6):
    exx, eyy, ezz, gxy, gyz, gxz = e6
    return (exx * x + 0.5 * gxy * y + 0.5 * gxz * z,
            0.5 * gxy * x + eyy * y + 0.5 * gyz * z,
            0.5 * gxz * x + 0.5 * gyz * y + ezz * z)


def _affine_disp_hencky(x, y, z, e6):
    """Displacement of node X for the deformation gradient F = exp(E), E the symmetric tensor of the ENGINEERING strain e6:
    the LogStrain host then sees the Hencky strain E exactly (b = F F^T = exp(2E)), i.e. the same trial the small-strain
    hosts see."""
    exx, eyy, ezz, gxy, gyz, gxz = e6
    E = np.array([[exx, 0.5 * gxy, 0.5 * gxz], [0.5 * gxy, eyy, 0.5 * gyz], [0.5 * gxz, 0.5 * gyz, ezz]])
    w, V = np.linalg.eigh(E)
    F = V @ np.diag(np.exp(w)) @ V.T
    u = (F - np.eye(3)) @ np.array([x, y, z])
    return float(u[0]), float(u[1]), float(u[2])


def _concrete(tag):
    ops.nDMaterial("LadrunoConcrete3D", tag, 30000.0, 0.2, 30.0, 3.0, 0.1, 30.0,
                   "-Df", 0.85, "-hardening", 0.3, 0.01,
                   "-tensionLaw", "bilinear", "-epsFc", 1.0e-3,
                   "-flowPotential", "cdpm2", "-compressionDrive", "cdpm2", "-tcTemper", "proj")


def _build(kind):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _NODES.items():
        ops.node(tag, x, y, z)
    _concrete(1)
    ops.fix(1, 1, 1, 1)
    if kind == "solidshell":
        ops.element("LadrunoSolidShell", 1, *_CONN, 1)
    elif kind == "logstrain":
        ops.nDMaterial("LogStrain", 2, 1)
        ops.element("LadrunoBrick", 1, *_CONN, 2, "-formulation", "std", "-geom", "finite")
    else:                                                    # control: an element that already forwards the sentinel
        ops.element("LadrunoBrick", 1, *_CONN, 1, "-formulation", "bbar")
    # pattern 1: the realistic cracking field, ramped by the Linear series over t in [0, 1]
    e6_crack = tuple(c * _CRACK_STEP * _CRACK_STEPS for c in _CRACK_DIR)
    disp = _affine_disp_hencky if kind == "logstrain" else _affine_disp
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n, (x, y, z) in _NODES.items():
        if n != 1:
            ux, uy, uz = disp(x, y, z, e6_crack)
            ops.sp(n, 1, ux); ops.sp(n, 2, uy); ops.sp(n, 3, uz)
    # pattern 2: the wild field, factor 0 until t = 1 and then a Linear ramp to 1 at t = 2 (held afterwards)
    ops.timeSeries("Path", 2, "-time", 1.0, 2.0, "-values", 0.0, 1.0, "-useLast")
    ops.pattern("Plain", 2, 2)
    e6_wild = _WILD_RAW
    for n, (x, y, z) in _NODES.items():
        if n != 1:
            ux, uy, uz = disp(x, y, z, e6_wild)
            ops.sp(n, 1, ux); ops.sp(n, 2, uy); ops.sp(n, 3, uz)
    ops.system("UmfPack")                                    # the CDPM2 tangent is non-symmetric
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-6, 20)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / _CRACK_STEPS)
    ops.analysis("Static")


def _disp():
    return [tuple(ops.nodeDisp(n)) for n in _NODES]


@pytest.mark.parametrize("kind", ["brick_control", "solidshell", "logstrain"])
def test_refused_trial_cuts_the_step_and_the_model_stays_consistent(kind):
    _build(kind)
    for _ in range(_CRACK_STEPS):
        assert ops.analyze(1) == 0, "the realistic cracking steps must converge"
    d0 = _disp()
    t0 = ops.getTime()

    ops.integrator("LoadControl", 1.0)                       # the full wild increment in one step: refused
    rc = ops.analyze(1)
    assert rc < 0, f"[{kind}] the wild step must fail (rc = {rc})"
    d1 = _disp()
    worst = max(abs(a - b) for u0, u1 in zip(d0, d1) for a, b in zip(u0, u1))
    assert worst < 1.0e-12, (
        f"[{kind}] the refused step left the committed NODAL state changed by {worst:.3e}: the step was ABORTED at commit "
        f"(model inconsistent), not cut at the trial")
    assert ops.getTime() == pytest.approx(t0, abs=1e-12), f"[{kind}] pseudo-time not reverted"

    # the cut: what a step-halving ladder does -- retry the SAME step with a smaller increment until it converges. Every
    # failed attempt must leave the committed nodal state untouched (a refusal is a cut, never an abort).
    dlam, rc = 0.05, -1
    for _ in range(8):
        ops.integrator("LoadControl", dlam)
        rc = ops.analyze(1)
        if rc == 0:
            break
        dl = max(abs(a - b) for u0, u1 in zip(d0, _disp()) for a, b in zip(u0, u1))
        assert dl < 1.0e-12, f"[{kind}] a refused retry (dlam {dlam}) changed the committed nodal state by {dl:.3e}"
        dlam *= 0.5
    assert rc == 0, f"[{kind}] no step-halving retry converged (last dlam {dlam})"
    assert all(abs(c) < 1.0 for u in _disp() for c in u)
