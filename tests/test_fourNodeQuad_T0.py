"""T0 worked example (feature class: element) — the single-element battery on
the existing FourNodeQuad, so a new continuum element author has a concrete,
runnable template for rank + constant-strain patch.

  * Rank / zero-energy: a free quad must have EXACTLY 3 zero modes (2D rigid
    body). 4 would be a spurious hourglass mode.
  * Constant-strain patch: a distorted 4-element patch with ONE free interior
    node, boundary nodes prescribed to a linear field. Every Gauss point of
    every element must report the same closed-form constant stress to 1e-9.

The strong (interior-node-free, distorted) form is deliberate — a rectangular
single element passes trivially and hides B-bar / incompatible-mode bugs.
"""
import pytest

from _testbed import ops
from _testbed.fem_checks import assert_zero_energy, check_constant_stress

pytestmark = [pytest.mark.zone_a, pytest.mark.t0]

E = 1000.0
NU = 0.25
THICK = 1.0


# ---------------------------------------------------------------- rank / modes
def test_rank_zero_energy():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    # one free, mildly distorted quad
    ops.node(1, 0.0, 0.0)
    ops.node(2, 1.0, 0.0)
    ops.node(3, 1.2, 1.1)
    ops.node(4, 0.1, 0.9)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    ops.element("quad", 1, 1, 2, 3, 4, THICK, "PlaneStress", 1)
    # exactly 3 rigid-body modes, no spurious hourglass
    assert_zero_energy([1, 2, 3, 4], ndf=2, ndm=2)


# ---------------------------------------------------------------- patch test
def _plane_stress_D(E, nu):
    f = E / (1.0 - nu * nu)
    g = E / (2.0 * (1.0 + nu))
    return f, g


def test_constant_strain_patch():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)

    # 3x3 distorted grid; node 5 is the free interior node.
    coords = {
        1: (0.0, 0.0), 2: (1.0, 0.0), 3: (2.0, 0.0),
        4: (0.0, 1.0), 5: (1.15, 1.05), 6: (2.0, 1.0),
        7: (0.0, 2.0), 8: (1.0, 2.0), 9: (2.0, 2.0),
    }
    for n, (x, y) in coords.items():
        ops.node(n, x, y)

    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    quads = {1: (1, 2, 5, 4), 2: (2, 3, 6, 5), 3: (4, 5, 8, 7), 4: (5, 6, 9, 8)}
    for e, nodes in quads.items():
        ops.element("quad", e, *nodes, THICK, "PlaneStress", 1)

    # linear displacement field -> constant strain
    ax, ay, bx, by = 1e-3, 4e-4, 2e-4, -3e-4
    eps_xx, eps_yy, gam_xy = ax, by, ay + bx

    def field(x, y):
        return ax * x + ay * y, bx * x + by * y

    boundary = [1, 2, 3, 4, 6, 7, 8, 9]  # everything except interior node 5
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in boundary:
        ux, uy = field(*coords[n])
        ops.sp(n, 1, ux)
        ops.sp(n, 2, uy)

    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static", "-noWarnings")
    ops.analyze(1)

    # interior node must land exactly on the field
    ux5, uy5 = field(*coords[5])
    assert abs(ops.nodeDisp(5, 1) - ux5) <= 1e-9 * abs(ux5) + 1e-12
    assert abs(ops.nodeDisp(5, 2) - uy5) <= 1e-9 * abs(uy5) + 1e-12

    # closed-form constant stress (plane stress)
    f, g = _plane_stress_D(E, NU)
    target = [f * (eps_xx + NU * eps_yy), f * (NU * eps_xx + eps_yy), g * gam_xy]

    strain_ref = max(abs(eps_xx), abs(eps_yy), abs(gam_xy))
    for e in quads:
        check_constant_stress(e, n_gauss=4, target_stress=target,
                              E=E, strain_ref=strain_ref)
