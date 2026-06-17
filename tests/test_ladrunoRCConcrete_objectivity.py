"""LadrunoRCConcrete — rigid-rotation (frame-indifference) objectivity gate (Phase 2b.2c).

The blocking design decision of the fixed-crack cyclic model (ADR 19, D1/D5) is that the
stored crack NORMAL + the aggregate-interlock direction are *directional internal variables*.
Scalar dt/dc damage is trivially frame-indifferent; a stored crack frame is NOT, unless the
constitutive response is genuinely objective — i.e. rotating the strain frame by Q must rotate
the stress frame by exactly the same Q. The supported large-rotation route is the corotational
ELEMENT (ASDShellQ4 -corotational), which feeds the small-strain material the *de-rotated*
strain Q^T E Q; so the property that makes that route objective is precisely

        sigma( Q E Q^T )  ==  Q sigma(E) Q^T            for every step of an arbitrary
                                                        (cracked, cyclic, X-cracked) history.

This test drives a homogeneous unit cube through a tension-then-cyclic-shear path that forms a
fixed crack and exercises the interlock, in the reference frame AND in a frame rigidly rotated
about z by a large angle, and asserts the two stress trajectories coincide after the Q-transform
to machine-ish precision. A non-objective directional state (a crack normal or interlock
projector that did not co-rotate) would break the identity wherever the crack is engaged.

Driver: the same all-DOF-prescribed (Penalty) homogeneous-strain stdBrick probe as the material
battery, but prescribing a FULL symmetric strain tensor u = E.X (so a rotation about z mixes the
in-plane components), built from two load patterns (tension scale, shear scale) carrying the
ROTATED basis tensors Q A Q^T, Q B Q^T.

Plan: Ladruno_implementation/19_ladruno_rc_shell_adr.md (Phase 2b.2c — rigid-rotation objectivity).
"""
import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}

E, NU, KC = 30000.0, 0.2, 2.0 / 3.0
CE = [0.0, 0.0007, 0.0020, 0.0100]
CS = [0.0, 24.0,   30.0,   5.0]
CD = [0.0, 0.0,    0.25,   1.0 - 5.0 / 45.0]
TE = [0.0, 0.0001, 0.0010]
TS = [0.0, 3.0,    0.5]
TD = [0.0, 0.0,    1.0 - 0.5 / 5.0]


def _mat(tag, **kw):
    args = ["LadrunoRCConcrete", tag, E, NU,
            "-Ce", *CE, "-Cs", *CS, "-Cd", *CD,
            "-Te", *TE, "-Ts", *TS, "-Td", *TD, "-Kc", KC,
            "-interlock", "-agg", 16.0, "-cyclic", "-crackSpacing", 50.0]
    for k, v in kw.items():
        args += ["-" + k] + ([] if v is True else [v])
    ops.nDMaterial(*args)


def _rot_z(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


# engineering -> tensor strain basis (symmetric): tension exx=1, shear gxy=1 (=> eps_xy=0.5)
_A = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
_B = np.array([[0.0, 0.5, 0.0], [0.5, 0.0, 0.0], [0.0, 0.0, 0.0]])


def _voigt_to_mat(s6):
    return np.array([[s6[0], s6[3], s6[5]],
                     [s6[3], s6[1], s6[4]],
                     [s6[5], s6[4], s6[2]]])


def _run_rotated_path(mat_fn, exx, gamma, Q, ts1, ts2, nper=40):
    """Homogeneous stdBrick driven through the strain history
        E(t) = ts1(t)*exx*(Q A Q^T) + ts2(t)*gamma*(Q B Q^T),
    prescribing u = E.X at every node/DOF via Penalty. Returns the list of GP stress 6-vectors
    over the recorded stages (stage 1 onward). With Q=I this is the standard tension-then-shear
    cyclic crack-shear path; with Q=R(theta) it is the same physical deformation in a rotated
    frame."""
    At = Q @ _A @ Q.T
    Bt = Q @ _B @ Q.T
    nst = len(ts1) - 1
    times = [float(i) for i in range(nst + 2)]
    ts1 = list(ts1) + [ts1[-1]]
    ts2 = list(ts2) + [ts2[-1]]
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, cc in _CUBE.items():
        ops.node(t, *cc)
    mat_fn(1)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Path", 1, "-time", *times, "-values", *ts1)
    ops.timeSeries("Path", 2, "-time", *times, "-values", *ts2)
    # pattern 1 = tension basis (scaled exx); pattern 2 = shear basis (scaled gamma)
    ops.pattern("Plain", 1, 1)
    for t, cc in _CUBE.items():
        X = np.array(cc, float)
        u = exx * (At @ X)
        for d in range(3):
            ops.sp(t, d + 1, float(u[d]))
    ops.pattern("Plain", 2, 2)
    for t, cc in _CUBE.items():
        X = np.array(cc, float)
        u = gamma * (Bt @ X)
        for d in range(3):
            ops.sp(t, d + 1, float(u[d]))
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Penalty", 1.0e14, 1.0e14)
    ops.test("NormDispIncr", 1.0e-9, 300, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nper)
    ops.analysis("Static")
    out = []
    for stage in range(nst):
        for _ in range(nper):
            assert ops.analyze(1) == 0, f"objectivity path analyze failed (stage {stage}, Q-rot)"
            if stage >= 1:
                ops.eleResponse(1, "forces")
                out.append(np.array(list(ops.eleResponse(1, "stresses"))[0:6]))
    return out


@pytest.mark.t1
@pytest.mark.parametrize("deg", [30.0, 90.0, 127.0])
def test_cyclic_cracked_state_is_frame_indifferent(deg):
    """sigma(Q E Q^T) == Q sigma(E) Q^T over a full reversing cracked path, for a large rigid
    frame rotation Q about z. The fixed crack normal, interlock projectors and friction state
    must all co-rotate so the stress trajectory is identical after the frame transform."""
    exx, gamma = 5.0e-4, 2.0e-3
    ts1 = [0, 1, 1, 1, 1]
    ts2 = [0, 0, 1, -1, 1]
    ref = _run_rotated_path(lambda t: _mat(t), exx, gamma, np.eye(3), ts1, ts2)
    Q = _rot_z(np.deg2rad(deg))
    rot = _run_rotated_path(lambda t: _mat(t), exx, gamma, Q, ts1, ts2)
    assert len(ref) == len(rot) and len(ref) > 0
    speak = max(np.max(np.abs(s)) for s in ref) + 1.0e-12
    worst = 0.0
    for sref, srot in zip(ref, rot):
        expected = Q @ _voigt_to_mat(sref) @ Q.T     # the reference stress in the rotated frame
        got = _voigt_to_mat(srot)
        worst = max(worst, float(np.max(np.abs(got - expected))))
    assert worst < 1.0e-5 * speak, (
        f"frame-indifference broken at {deg} deg: max |sigma(QEQ^T) - Q sigma(E) Q^T| = "
        f"{worst} (peak |sigma| = {speak})")


@pytest.mark.t1
def test_xcrack_state_is_frame_indifferent():
    """Same objectivity identity with the full Phase-2b.2b X-cracking + cumulative wear engaged
    (a biaxial path that forms BOTH orthogonal cracks): the second crack normal + the governing-
    cap selection + the slip-driven wear must all be frame-indifferent too."""
    # biaxial tension (exx on both in-plane axes via the rotated basis is awkward; instead form
    # crack 2 by adding a transverse tension component through a 2nd tension-like basis). Here a
    # simple route: tension + reversing shear with -xcrack; crack 2 forms once the rotated shear
    # opens the orthogonal plane. The identity must hold regardless of which cracks engage.
    exx, gamma = 6.0e-4, 3.0e-3
    ts1 = [0, 1, 1, 1, 1, 1]
    ts2 = [0, 0, 1, -1, 1, -1]
    ref = _run_rotated_path(lambda t: _mat(t, xcrack=True, degKappa=0.5, degSlipRef=0.02),
                            exx, gamma, np.eye(3), ts1, ts2, nper=40)
    Q = _rot_z(np.deg2rad(63.0))
    rot = _run_rotated_path(lambda t: _mat(t, xcrack=True, degKappa=0.5, degSlipRef=0.02),
                            exx, gamma, Q, ts1, ts2, nper=40)
    speak = max(np.max(np.abs(s)) for s in ref) + 1.0e-12
    worst = max(float(np.max(np.abs(_voigt_to_mat(srot) - Q @ _voigt_to_mat(sref) @ Q.T)))
                for sref, srot in zip(ref, rot))
    assert worst < 1.0e-5 * speak, f"xcrack frame-indifference broken: worst={worst} peak={speak}"
