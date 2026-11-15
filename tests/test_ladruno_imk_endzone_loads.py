"""LadrunoIMKBeam(2d): element loads under -endZones (the end-zone load fix).

Before this fix, an element load on a LadrunoIMKBeam with `-endZones dI dJ`
(1) left out the face reactions' lever arm: the node moments were short by
p0 * dI (p0 * dJ), so a free node got a spurious unbalanced moment and the
support moments were wrong; (2) dropped any load lying on a rigid zone; and
(3) read partial/point positions as fractions of the clear span Lc instead of
the full length L, so the loads LadrunoMembrane pushes to its -frames shrank
by Lc / L. `hingeMoment` also reported the node basic moments, not the face
moments the hinges carry.

Oracle: the same member built from three `elasticBeamColumn` pieces, the two
end pieces 1e6 times stiffer in bending and torsion (the rigid zones; the IMK
beam's axial stiffness is EA/L over the full L, so the stubs keep A), carrying
the same load. With
elastic ends the IMK beam must match its node displacements and support
reactions. Loads: uniform, a partial load crossing both zone boundaries, and
point loads on the i zone, the clear span and the j zone; 3-D and 2-D. On a
cantilever the free tip is a node whose equilibrium the bug broke.
"""
import math

import numpy as np
import pytest
from _testbed import ops
from _testbed.roundtrip import database_roundtrip

pytestmark = [pytest.mark.zone_a]

L, DI, DJ = 5.0, 0.4, 0.7
E, G, A, J, IY, IZ = 30.0e6, 12.5e6, 0.18, 0.004, 0.003, 0.005
STIFF = 1.0e6
XS = (0.0, DI, L - DJ, L)            # reference node abscissae


def _solve():
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-14, 20)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    ops.reactions()


# --- loads: (kind, params) on the full member, positions as fractions of L ---
LOADS_3D = {
    "uniform": ("uniform", (-12.0, 7.0, 3.0)),
    "partial": ("partial", (-10.0, 4.0, 2.0, 0.02, 0.97, -25.0, -6.0, 5.0)),
    # eleLoad input order (3-D): Py Pz xL Px
    "point_zone_i": ("point", (-30.0, 12.0, 0.05, 4.0)),
    "point_clear": ("point", (-30.0, 12.0, 0.45, 4.0)),
    "point_zone_j": ("point", (-30.0, 12.0, 0.93, 4.0)),
}
LOADS_2D = {
    "uniform": ("uniform", (-12.0, 3.0)),
    # eleLoad input order (2-D): Wya Wxa aOverL bOverL Wyb Wxb
    "partial": ("partial", (-10.0, 2.0, 0.02, 0.97, -25.0, 5.0)),
    # eleLoad input order (2-D): Py xL Px
    "point_zone_i": ("point", (-30.0, 0.05, 4.0)),
    "point_clear": ("point", (-30.0, 0.45, 4.0)),
    "point_zone_j": ("point", (-30.0, 0.93, 4.0)),
}


def _lerp(w1, w2, x1, x2, x):
    return w1 + (w2 - w1) * (x - x1) / (x2 - x1)


def _apply_on_piece(ndm, ele, x0, x1, kind, prm):
    """Apply the member load `kind/prm` restricted to the piece [x0, x1]."""
    ln = x1 - x0
    if kind == "uniform":
        ops.eleLoad("-ele", ele, "-type", "-beamUniform", *prm)
        return
    if kind == "point":
        k = 2 if ndm == 3 else 1          # index of xL in the input
        x = prm[k] * L
        if x0 <= x <= x1 and not (x == x1 and x1 < L):
            local = list(prm)
            local[k] = (x - x0) / ln
            ops.eleLoad("-ele", ele, "-type", "-beamPoint", *local)
        return
    # partial, linear between xa and xb
    if ndm == 3:
        wya, wza, wxa, fa, fb, wyb, wzb, wxb = prm
        wa, wb = (wya, wza, wxa), (wyb, wzb, wxb)
    else:
        wta, waa, fa, fb, wtb, wab = prm
        wa, wb = (wta, waa), (wtb, wab)
    xa, xb = fa * L, fb * L
    c1, c2 = max(x0, xa), min(x1, xb)
    if c2 <= c1:
        return
    w1 = [_lerp(a, b, xa, xb, c1) for a, b in zip(wa, wb)]
    w2 = [_lerp(a, b, xa, xb, c2) for a, b in zip(wa, wb)]
    r1, r2 = (c1 - x0) / ln, (c2 - x0) / ln
    if ndm == 3:
        ops.eleLoad("-ele", ele, "-type", "-beamUniform", *w1, r1, r2, *w2)
    else:
        ops.eleLoad("-ele", ele, "-type", "-beamUniform", w1[0], w1[1], r1, r2, w2[0], w2[1])


def _cantilever(ndm, model, kind, prm):
    """Node 1 fixed at x = 0, node 2 free at x = L. Returns (tip disp, base reaction)."""
    ops.wipe()
    ndf = 6 if ndm == 3 else 3
    ops.model("basic", "-ndm", ndm, "-ndf", ndf)
    coords = (lambda x: (x, 0.0, 0.0)) if ndm == 3 else (lambda x: (x, 0.0))
    ops.node(1, *coords(0.0))
    ops.node(2, *coords(L))
    ops.fix(1, *([1] * ndf))
    if ndm == 3:
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    else:
        ops.geomTransf("Linear", 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    if model == "imk":
        if ndm == 3:
            ops.element("LadrunoIMKBeam", 1, 1, 2, A, E, G, J, IY, IZ, 1, "-endZones", DI, DJ)
        else:
            ops.element("LadrunoIMKBeam", 1, 1, 2, A, E, IZ, 1, "-endZones", DI, DJ)
        if kind == "uniform":
            ops.eleLoad("-ele", 1, "-type", "-beamUniform", *prm)
        elif kind == "point":
            ops.eleLoad("-ele", 1, "-type", "-beamPoint", *prm)
        else:
            ops.eleLoad("-ele", 1, "-type", "-beamUniform", *prm)
    else:
        ops.node(3, *coords(XS[1]))
        ops.node(4, *coords(XS[2]))
        pieces = [(11, 1, 3, STIFF), (12, 3, 4, 1.0), (13, 4, 2, STIFF)]
        for (tag, ni, nj, k), x0, x1 in zip(pieces, XS[:-1], XS[1:]):
            if ndm == 3:
                ops.element("elasticBeamColumn", tag, ni, nj, A, E, G, J * k, IY * k, IZ * k, 1)
            else:
                ops.element("elasticBeamColumn", tag, ni, nj, A, E, IZ * k, 1)
            _apply_on_piece(ndm, tag, x0, x1, kind, prm)
    _solve()
    return np.array(ops.nodeDisp(2)), np.array(ops.nodeReaction(1))


@pytest.mark.parametrize("ndm", [3, 2])
@pytest.mark.parametrize("case", list(LOADS_3D))
def test_endzone_loads_match_rigid_stub_reference(ndm, case):
    kind, prm = (LOADS_3D if ndm == 3 else LOADS_2D)[case]
    u_imk, r_imk = _cantilever(ndm, "imk", kind, prm)
    u_ref, r_ref = _cantilever(ndm, "ref", kind, prm)
    np.testing.assert_allclose(r_imk, r_ref, rtol=1e-5, atol=1e-6 * np.abs(r_ref).max())
    # atol floor 1e-9: the 1e6-stiff stubs are not rigid, so a load on a zone of the
    # reference still moves the tip ~1e-11 where the IMK tip correctly stays at 0.
    np.testing.assert_allclose(u_imk, u_ref, rtol=1e-4,
                               atol=max(1e-5 * np.abs(u_ref).max(), 1e-9))


@pytest.mark.parametrize("ndm", [3, 2])
def test_endzone_uniform_load_is_fully_equilibrated(ndm):
    """The base reaction carries the whole load w L and its moment w L^2 / 2,
    the rigid-zone parts included (they were dropped before the fix)."""
    w = -12.0
    prm = (w, 0.0, 0.0) if ndm == 3 else (w, 0.0)
    _, r = _cantilever(ndm, "imk", "uniform", prm)
    fy, mz = (r[1], r[5]) if ndm == 3 else (r[1], r[2])
    assert math.isclose(fy, -w * L, rel_tol=1e-12)
    assert math.isclose(mz, -w * L * L / 2.0, rel_tol=1e-12)


def _fixed_fixed_3d(w, my=None, steps=1, lam=1.0, dz=0.5):
    """One element, both ends fixed, -endZones dz dz, uniform local-y load."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.node(3, L / 2, 0.0, 0.0)  # unused, keeps tags apart from the cantilever
    ops.fix(1, *[1] * 6)
    ops.fix(2, *[1] * 6)
    ops.fix(3, *[1] * 6)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    hinge = []
    if my is not None:
        kh = 1.0e4 * 4.0 * E * IZ / (L - 2 * dz)
        ops.uniaxialMaterial("ElasticPP", 10, kh, my / kh)
        hinge = ["-matZ", 10]
    ops.element("LadrunoIMKBeam", 1, 1, 2, A, E, G, J, IY, IZ, 1, *hinge, "-endZones", dz, dz)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 1, "-type", "-beamUniform", w, 0.0, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-14, 50)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", lam / steps)
    ops.analysis("Static")


def test_hinge_moment_is_the_face_moment():
    """Elastic fixed-fixed clear span: hingeMoment = w Lc^2 / 12 at both faces;
    the node moment adds the face shear's and the zone load's lever arms."""
    w, dz = -20.0, 0.5
    lc = L - 2 * dz
    _fixed_fixed_3d(w, dz=dz)
    assert ops.analyze(1) == 0
    ops.reactions()
    hm = ops.eleResponse(1, "hingeMoment")
    face = abs(w) * lc * lc / 12.0
    assert math.isclose(abs(hm[0]), face, rel_tol=1e-10)
    assert math.isclose(abs(hm[1]), face, rel_tol=1e-10)
    m_node = face + (abs(w) * lc / 2.0) * dz + abs(w) * dz * dz / 2.0
    assert math.isclose(abs(ops.nodeReaction(1)[5]), m_node, rel_tol=1e-10)


def test_face_hinges_yield_at_the_predicted_load_factor():
    """My = 0.8 x the elastic face moment: the face hinges yield at lambda = 0.8
    and then hold My. With end springs kh the symmetric elastic face moment is
    w Lc^2 / 12 / (1 + 2 EI / (kh Lc))."""
    w, dz = -20.0, 0.5
    lc = L - 2 * dz
    kh = 1.0e4 * 4.0 * E * IZ / lc          # as _fixed_fixed_3d builds it
    face_el = abs(w) * lc * lc / 12.0 / (1.0 + 2.0 * E * IZ / (kh * lc))
    my = 0.8 * face_el
    _fixed_fixed_3d(w, my=my, steps=50, lam=1.0, dz=dz)
    onset = None
    for _ in range(50):
        assert ops.analyze(1) == 0
        hm = ops.eleResponse(1, "hingeMoment")
        if onset is None and abs(hm[0]) >= my * (1.0 - 1e-6):  # the hinge solve's own tolerance
            onset = ops.getLoadFactor(1)
    assert onset is not None and abs(onset - 0.8) < 1e-9, onset
    assert math.isclose(abs(hm[0]), my, rel_tol=1e-9)
    assert math.isclose(abs(hm[1]), my, rel_tol=1e-9)


def test_endzone_load_state_round_trips():
    """dM0 is serialized: localForce of the restored element is bit-identical."""
    kind, prm = LOADS_3D["partial"]

    def build():
        _cantilever(3, "imk", kind, prm)

    database_roundtrip(build, probe_nodes=[2], ndf=6,
                       probe_fn=lambda: list(ops.eleResponse(1, "localForce")),
                       rounds=2)
