"""Remaining-beam ``localAxes`` recorder response (Ladruno) — Zone-A.

The Ladruno recorder records ``MODEL/LOCAL_AXES`` for any element that answers the
``"localAxes"`` element response (id 30 -> 9 packed direction cosines from the
element's ``CrdTransf``). Originally only the Elastic/Force/Disp beam families were
wired; this battery covers the *remaining* standard-transform beams now wired:

  * MixedBeamColumn3d / MixedBeamColumn2d
  * GradientInelasticBeamColumn3d / GradientInelasticBeamColumn2d

What it proves (and would catch the failure mode we deliberately avoided on
DispBeamColumn2dInt, whose transform's ``getLocalAxes`` zeros the frame):

  * the response exists and returns 9 finite numbers,
  * the frame is NON-degenerate (not all-zero) and ORTHONORMAL,
  * the local x-axis equals the normalized element axis (so a wrong/identity
    frame, or zeros, fails).

Construction note: every element takes ``(tag, iNode, jNode, transfTag,
integrationTag[, lc])`` and consumes a ``beamIntegration`` rule + section(s).
GradientInelastic requires identical internal-section tags, so we use a single
Elastic section for all integration points. localAxes depends only on the
CrdTransf, so an Elastic section is sufficient and keeps the test fast.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

E = 2.0e5
G = 8.0e4
A = 100.0
IZ = 1.0e4
IY = 8.0e3
J = 1.5e4
NP = 5          # Lobatto points (>=3 for the GI internal-point rule)
LC = 0.1        # GI characteristic length (0 < lc < L)


def _norm(v):
    return math.sqrt(sum(c * c for c in v))


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def _assert_orthonormal_frame(la, axis_expected, *, twod):
    assert len(la) == 9, f"localAxes must be 9 packed cosines, got {len(la)}"
    assert all(math.isfinite(c) for c in la), f"non-finite localAxes: {la}"
    vx, vy, vz = la[0:3], la[3:6], la[6:9]

    # NON-degenerate: the deliberate failure mode (base CrdTransf::getLocalAxes
    # zeros the frame) would trip exactly here.
    for name, v in (("vx", vx), ("vy", vy), ("vz", vz)):
        assert _norm(v) > 1e-9, f"{name} is a zero/degenerate axis: {v}"

    # unit length
    for name, v in (("vx", vx), ("vy", vy), ("vz", vz)):
        assert abs(_norm(v) - 1.0) <= 1e-9, f"{name} not unit: |{name}|={_norm(v)}"

    # mutually orthogonal
    assert abs(_dot(vx, vy)) <= 1e-9, f"vx.vy={_dot(vx, vy)}"
    assert abs(_dot(vx, vz)) <= 1e-9, f"vx.vz={_dot(vx, vz)}"
    assert abs(_dot(vy, vz)) <= 1e-9, f"vy.vz={_dot(vy, vz)}"

    # local x == normalized element axis
    n = _norm(axis_expected)
    ax = [c / n for c in axis_expected]
    for i in range(3):
        assert abs(vx[i] - ax[i]) <= 1e-9, f"vx{tuple(vx)} != element axis {tuple(ax)}"

    if twod:
        # 2D transforms pack the out-of-plane normal as vz=(0,0,1)
        assert abs(vz[0]) <= 1e-9 and abs(vz[1]) <= 1e-9 and abs(vz[2] - 1.0) <= 1e-9, \
            f"2D vz should be (0,0,1), got {tuple(vz)}"


def _build_3d(elem_kind):
    """elem_kind in {'mixed','gi'}. Returns (eleTag, expected_axis)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    iN, jN = (0.0, 0.0, 0.0), (1.0, 2.0, 2.0)   # length 3, deliberately oblique
    ops.node(1, *iN)
    ops.node(2, *jN)
    ops.section("Elastic", 1, E, A, IZ, IY, G, J)
    ops.beamIntegration("Lobatto", 1, 1, NP)
    # vecxz not parallel to the element axis
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    if elem_kind == "mixed":
        ops.element("mixedBeamColumn", 10, 1, 2, 1, 1)
    else:
        ops.element("gradientInelasticBeamColumn", 10, 1, 2, 1, 1, LC)
    return 10, [jN[i] - iN[i] for i in range(3)]


def _build_2d(elem_kind):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    iN, jN = (0.0, 0.0), (1.0, 2.0)             # length sqrt(5), oblique
    ops.node(1, *iN)
    ops.node(2, *jN)
    ops.section("Elastic", 1, E, A, IZ)
    ops.beamIntegration("Lobatto", 1, 1, NP)
    ops.geomTransf("Linear", 1)
    if elem_kind == "mixed":
        ops.element("mixedBeamColumn", 10, 1, 2, 1, 1)
    else:
        ops.element("gradientInelasticBeamColumn", 10, 1, 2, 1, 1, LC)
    # pack expected axis as a 3-vector with 0 out-of-plane component
    return 10, [jN[0] - iN[0], jN[1] - iN[1], 0.0]


def test_mixedBeamColumn3d_localAxes():
    tag, axis = _build_3d("mixed")
    la = ops.eleResponse(tag, "localAxes")
    _assert_orthonormal_frame(la, axis, twod=False)


def test_mixedBeamColumn2d_localAxes():
    tag, axis = _build_2d("mixed")
    la = ops.eleResponse(tag, "localAxes")
    _assert_orthonormal_frame(la, axis, twod=True)


def test_gradientInelasticBeamColumn3d_localAxes():
    tag, axis = _build_3d("gi")
    la = ops.eleResponse(tag, "localAxes")
    _assert_orthonormal_frame(la, axis, twod=False)


def test_gradientInelasticBeamColumn2d_localAxes():
    tag, axis = _build_2d("gi")
    la = ops.eleResponse(tag, "localAxes")
    _assert_orthonormal_frame(la, axis, twod=True)
