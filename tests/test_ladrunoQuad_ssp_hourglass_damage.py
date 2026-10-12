"""LadrunoQuad ``-formulation ssp`` — Tier-A hourglass degradation under softening (C1).

The ssp hourglass stiffness ``Kstab`` is built from the INITIAL tangent and scaled by
``s = max(floor, 1 - omega)``. Before C1, omega was read at the centroid only and the floor
was 0.01. On the Kormeling-Reinhardt notched beam that left a 0.28-0.36 kN post-peak plateau
(std: 0.07-0.12 kN). Instrumenting the stabilization work (``eleResponse(e, 'hourglass')``,
``sum d.f_stab / u_load``) showed that 0.15-0.20 kN of the 0.27-0.32 kN tail was stabilization
work, ~96 % of it in the three UNDAMAGED top-row elements under the load: the last ligament
bends in the element's hourglass mode, the centroid sits on the neutral axis and never cracks
(omega = 0), so the elastic hourglass spring is an uncrackable hinge. Lowering the floor alone
(0.01 -> 1e-4) only moved the tail from 0.32 to 0.27 kN.

C1: omega is also sampled at 4 shadow material points driven by the full 2x2-Gauss bilinear
strain (committed, one step lag), plus the secant energy loss at those points, and the floor
is 1e-4. ``-hgDamage centroid`` / ``-hgLegacy`` / ``-hourglassFloor 1`` keep the pre-C1 and the
stock-SSPquad behaviours reachable.

This file: the closed-form single-element hinge (T1), the escapes, the work identity of the
``hourglass`` diagnostic, and the parser guards. The structural gate is the K&R benchmark in
the validation repo (numbers in LadrunoPlaneElements_guide.md section 1.5).
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

E, NU, FC, FT, GF, GC = 20.0e9, 0.2, 24.0e6, 2.4e6, 100.0, 25000.0
L = 0.01           # 10 mm element
THK = 0.1
EPS_CR = FT / E


def _model(quad_flags, material="ladruno"):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in {1: (0, 0), 2: (L, 0), 3: (L, L), 4: (0, L)}.items():
        ops.node(tag, float(x), float(y))
    if material == "ladruno":
        ops.nDMaterial("LadrunoConcrete3D", 1, E, NU, FC, FT, GF, GC, "-lch", L)
    else:
        ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    ops.element("LadrunoQuad", 1, 1, 2, 3, 4, 1, "-formulation", "ssp",
                "-type", "PlaneStrain", "-thick", THK, *quad_flags)
    # one dummy free dof (a spring to ground) so the all-prescribed quad has an equation
    ops.node(5, 2 * L, 0.0)
    ops.node(6, 3 * L, 0.0)
    ops.fix(6, 1, 1)
    ops.fix(5, 0, 1)
    ops.uniaxialMaterial("Elastic", 9, 1.0e6)
    ops.element("truss", 2, 5, 6, 1.0, 9)


def _drive_hourglass(quad_flags, amp_over_cr=60.0, nsteps=120, material="ladruno"):
    """Impose the pure x-hourglass (bending) mode u_x = a*xi*eta, u_y = 0 on every node.
    Centroid strain is exactly zero for this square; the Gauss points see
    eps_xx = +/- a/(sqrt(3) L/2) -- half of them in tension."""
    _model(quad_flags, material)
    a_max = amp_over_cr * EPS_CR * (L / 2.0) * math.sqrt(3.0)   # GP strain = amp_over_cr * eps_cr
    h = {1: 1.0, 2: -1.0, 3: 1.0, 4: -1.0}                      # xi*eta at the nodes
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (1, 2, 3, 4):
        ops.sp(n, 1, h[n] * a_max)
        ops.sp(n, 2, 0.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 20)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    out = []
    for _ in range(nsteps):
        assert ops.analyze(1) == 0
        f = ops.eleResponse(1, "force")
        hg = ops.eleResponse(1, "hourglass")
        out.append((ops.getTime(), math.sqrt(sum(v * v for v in f)), hg))
    return out


def test_centroid_sampling_is_blind_to_the_hourglass_crack():
    """Pre-C1 (-hgDamage centroid): the centroid never cracks, s stays 1 and the hourglass
    force grows elastically to 60x the cracking strain -- the uncrackable hinge."""
    r = _drive_hourglass(["-hgDamage", "centroid"])
    assert all(row[2][0] == 1.0 for row in r)
    assert r[-1][1] > 0.95 * r[-1][0] / r[0][0] * r[0][1]    # linear in the amplitude


def test_gauss_sampling_cracks_the_hinge():
    """C1 default: the shadow Gauss points crack, s falls to the floor band and the
    element's force at 60x the cracking strain is < 2 % of the pre-C1 value."""
    old = _drive_hourglass(["-hgDamage", "centroid"])
    new = _drive_hourglass([])
    s_end = new[-1][2][0]
    assert s_end < 0.02, s_end
    assert new[-1][2][1] > 0.98                               # omega seen (shadow)
    assert new[-1][1] < 0.02 * old[-1][1], (new[-1][1], old[-1][1])
    # before cracking the two are identical (shadows are elastic -> omega = 0)
    assert new[0][1] == pytest.approx(old[0][1], rel=1e-12)


def test_hglegacy_is_centroid_plus_old_floor():
    leg = _drive_hourglass(["-hgLegacy"], nsteps=40)
    cen = _drive_hourglass(["-hgDamage", "centroid", "-hourglassFloor", 0.01], nsteps=40)
    for a, b in zip(leg, cen):
        assert a[1] == b[1]


def test_elastic_material_unchanged_and_no_shadow_cost():
    """No damage channel => s = 1 for every option (the ssp vs SSPquad parity tests in
    test_ladrunoPlane_element.py stay valid); the stiffness is bit-identical."""
    Ks = []
    for flags in ([], ["-hgLegacy"], ["-hourglassFloor", 1.0], ["-hgDamage", "centroid"]):
        _model(flags, material="elastic")
        Ks.append(ops.eleResponse(1, "stiff"))
    for K in Ks[1:]:
        assert K == Ks[0]


def test_hourglass_response_work_identity():
    """d.f_stab + d.f_phys == d.f_int (no body force / pressure) at every step."""
    _model([])
    h = {1: 1.0, 2: -1.0, 3: 1.0, 4: -1.0}
    a = 3.0 * EPS_CR * L
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (1, 2, 3, 4):
        ops.sp(n, 1, h[n] * a + 0.5 * a * (ops.nodeCoord(n)[0] / L))   # hourglass + stretch
        ops.sp(n, 2, 0.2 * a * (ops.nodeCoord(n)[1] / L))
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 20)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.05)
    ops.analysis("Static")
    for _ in range(20):
        assert ops.analyze(1) == 0
        f = ops.eleResponse(1, "force")
        d = [c for n in (1, 2, 3, 4) for c in ops.nodeDisp(n)]
        s, om, dfs, dfp = ops.eleResponse(1, "hourglass")
        work = sum(x * y for x, y in zip(d, f))
        assert dfs + dfp == pytest.approx(work, rel=1e-9, abs=1e-12 * abs(work) + 1e-20)
        assert 1.0e-4 <= s <= 1.0


@pytest.mark.parametrize("flags", [["-hourglassFloor", 1.5], ["-hourglassFloor", -0.1],
                                   ["-hgDamage", "bogus"]])
def test_parser_refuses_bad_values(flags):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in {1: (0, 0), 2: (L, 0), 3: (L, L), 4: (0, L)}.items():
        ops.node(tag, float(x), float(y))
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    with pytest.raises(Exception):
        ops.element("LadrunoQuad", 1, 1, 2, 3, 4, 1, "-formulation", "ssp", *flags)
