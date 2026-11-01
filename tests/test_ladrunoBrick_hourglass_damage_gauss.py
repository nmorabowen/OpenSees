"""LadrunoBrick ``ssp`` / ``uri -hourglass stiffness`` — Tier-A degradation sampled at the
Gauss points (Ladruno C3a, the C1 LadrunoQuad treatment).

The single-point stabilization (``Kstab`` for ssp, FB ``kappa`` for uri-stiffness) is scaled
by ``s = max(floor, 1 - omega)``. Before C3a, omega was the CENTROID damage and the floor 0.01:
an element deformed in a bending/hourglass mode has zero centroid strain, never damages at
its centroid, and its elastic stabilization is an uncrackable hinge. C3a adds the committed
damage (and secant loss) of 8 shadow material copies driven with the full trilinear strain at
the 2x2x2 Gauss points, and floor 1e-4. Escapes: ``-hgDamage centroid``, ``-hourglassFloor f``
(1 = no degradation), ``-hgLegacy`` (= centroid + 0.01, the pre-C3a element).
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

E, NU, FT = 30000.0, 0.2, 3.0
ET0 = FT / E
TE = [ET0, 1.0e-3, 8.0e-3]
TS = [FT, 1.0, 0.05]
TD = [0.0, 1.0 - 1.0 / (E * 1.0e-3), 1.0 - 0.05 / (E * 8.0e-3)]
CE, CS, CD = [-ET0, -1.0e-3], [-FT, -1.0], [0.0, 0.0]
L = 1.0
CUBE = {1: (0, 0, 0), 2: (L, 0, 0), 3: (L, L, 0), 4: (0, L, 0),
        5: (0, 0, L), 6: (L, 0, L), 7: (L, L, L), 8: (0, L, L)}
FORMS = [["-formulation", "ssp"], ["-formulation", "uri", "-hourglass", "stiffness"]]


def _mat(kind):
    if kind == "asd":
        ops.nDMaterial("ASDConcrete3D", 1, E, NU, "-rho", 0.0,
                       "-Te", *TE, "-Ts", *TS, "-Td", *TD,
                       "-Ce", *CE, "-Cs", *CS, "-Cd", *CD, "-autoRegularization", 1.0)
    elif kind == "ladruno":
        ops.nDMaterial("LadrunoConcrete3D", 1, E, NU, 30.0, FT, 0.1, 25.0, "-lch", L)
    elif kind == "j2":
        ops.nDMaterial("J2Plasticity", 1, E / (3 * (1 - 2 * NU)), E / (2 * (1 + NU)),
                       20.0, 30.0, 10.0, 100.0)
    else:
        ops.nDMaterial("ElasticIsotropic", 1, E, NU)


def _model(kind, flags):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in CUBE.items():
        ops.node(t, *map(float, c))
    _mat(kind)
    ops.element("LadrunoBrick", 1, *CUBE.keys(), 1, *flags)


def _drive_bending_hourglass(kind, flags, amp_over_cr=60.0, nsteps=120):
    """u_x = a * xi * eta on every node (pure x-bending hourglass): the centroid strain is
    exactly zero; the Gauss points see eps_xx = +/- a*eta_g/(L/2), half in tension."""
    _model(kind, flags)
    a = amp_over_cr * ET0 * (L / 2.0) * math.sqrt(3.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in CUBE.items():
        xi, eta = 2 * x / L - 1, 2 * y / L - 1
        ops.sp(t, 1, a * xi * eta)
        ops.sp(t, 2, 0.0)
        ops.sp(t, 3, 0.0)
    ops.constraints("Penalty", 1e14, 1e14)
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-12, 20)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    out = []
    for _ in range(nsteps):
        assert ops.analyze(1) == 0
        f = ops.eleResponse(1, "forces")
        out.append((math.sqrt(sum(v * v for v in f)), list(ops.eleResponse(1, "hourglass"))))
    return out


@pytest.mark.parametrize("kind", ["asd", "ladruno"])
@pytest.mark.parametrize("form", FORMS)
def test_centroid_blind_vs_gauss_cracked_hinge(kind, form):
    """-hgDamage centroid: s stays 1 through 60x the cracking strain (the hinge). Default
    (gauss): the shadows crack, s falls to <= 2 %, and the element force at the end is
    < 5 % of the centroid-sampled one. Before cracking both are identical."""
    old = _drive_bending_hourglass(kind, form + ["-hgDamage", "centroid"])
    new = _drive_bending_hourglass(kind, form)
    assert all(r[1][0] == 1.0 for r in old)
    assert new[-1][1][0] < 0.02, new[-1][1]
    assert new[-1][0] < 0.05 * old[-1][0], (new[-1][0], old[-1][0])
    assert new[0][0] == pytest.approx(old[0][0], rel=1e-12)


@pytest.mark.parametrize("form", FORMS)
def test_hglegacy_equals_centroid_plus_old_floor(form):
    a = _drive_bending_hourglass("asd", form + ["-hgLegacy"], nsteps=40)
    b = _drive_bending_hourglass("asd", form + ["-hgDamage", "centroid", "-hourglassFloor", 0.01],
                                 nsteps=40)
    for x, y in zip(a, b):
        assert x[0] == y[0]


@pytest.mark.parametrize("kind", ["elastic", "j2"])
@pytest.mark.parametrize("form", FORMS)
def test_no_damage_channel_stiffness_bit_identical(kind, form):
    """Materials without a damage channel: s = 1 and no shadows under every option."""
    Ks = []
    for flags in ([], ["-hgLegacy"], ["-hourglassFloor", 1.0], ["-hgDamage", "centroid"]):
        _model(kind, form + flags)
        Ks.append(list(ops.eleResponse(1, "stiff")))
    for K in Ks[1:]:
        assert K == Ks[0]


@pytest.mark.parametrize("form", FORMS)
def test_hourglass_response_work_identity(form):
    """dfStab == 2*hourglassEnergy and dfStab + dfPhys == d.F_int, every step."""
    _model("asd", form)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    a = 20 * ET0 * L
    for t, (x, y, z) in CUBE.items():
        xi, eta = 2 * x / L - 1, 2 * y / L - 1
        ops.sp(t, 1, a * xi * eta + 0.5 * a * x)
        ops.sp(t, 2, 0.1 * a * y)
        ops.sp(t, 3, 0.0)
    ops.constraints("Penalty", 1e14, 1e14)
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 0.05)
    ops.analysis("Static")
    for _ in range(20):
        assert ops.analyze(1) == 0
        f = ops.eleResponse(1, "forces")
        d = [c for n in CUBE for c in ops.nodeDisp(n)]
        s, om, dfs, dfp = ops.eleResponse(1, "hourglass")
        ehg = ops.eleResponse(1, "hourglassEnergy")[0]
        W = sum(x * y for x, y in zip(d, f))
        assert dfs == pytest.approx(2.0 * ehg, rel=1e-12, abs=1e-30)
        assert dfs + dfp == pytest.approx(W, rel=1e-9, abs=1e-18)
        assert 1e-4 <= s <= 1.0


@pytest.mark.parametrize("flags", [["-hourglassFloor", 2.0], ["-hourglassFloor", -1.0],
                                   ["-hgDamage", "nope"]])
def test_parser_refuses(flags):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in CUBE.items():
        ops.node(t, *map(float, c))
    _mat("elastic")
    with pytest.raises(Exception):
        ops.element("LadrunoBrick", 1, *CUBE.keys(), 1, "-formulation", "ssp", *flags)
