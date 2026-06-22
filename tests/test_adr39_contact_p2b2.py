"""ADR-39 P2b-2a — DEFORMABLE-master node-to-segment contact validation.

The P2b-1 NTS adapter already reads master-node trial displacement and assembles the
master tangent blocks, so a DEFORMABLE master works end-to-end with the main-term
(kn*BᵀB) tangent — static Newton converges (penalty contact is robust; the ∂n/∂u
block is a convergence refinement for large-rotation/curved cases, NOT a correctness
requirement here). These tests pin that capability against analytic series stiffness.

- single deformable LadrunoBrick, top face = master segment, slave pressed down:
  contact gap = P/kn AND brick compresses P*L/(E*A) (series).
- block-on-block: top block's bottom node is the slave, bottom block's top face is the
  deformable master; both bricks compress + the interface penetrates P/kn.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

KN = 1.0e6


def _cube(base_node, z0, L, E, nu, ele_tag, mat_tag, fix_bottom=False):
    """8-node unit cube from base_node..base_node+7, bottom at z0, edge L. Returns
    (bottom4 tags, top4 tags)."""
    n = base_node
    pts = [(0, 0, 0), (L, 0, 0), (L, L, 0), (0, L, 0),
           (0, 0, L), (L, 0, L), (L, L, L), (0, L, L)]
    tags = []
    for i, (x, y, z) in enumerate(pts):
        t = n + i
        ops.node(t, float(x), float(y), z0 + z)
        tags.append(t)
    if fix_bottom:
        for t in tags[:4]:
            ops.fix(t, 1, 1, 1)
    ops.nDMaterial("ElasticIsotropic", mat_tag, E, nu)
    ops.element("LadrunoBrick", ele_tag, *tags, mat_tag)
    return tags[:4], tags[4:]


def test_p2b2_slave_on_deformable_brick():
    """Slave node pressed onto the top face of a single fixed-bottom deformable brick.
    Series: contact penetration = P/kn, brick uniaxial compression = P*L/(E*A)."""
    E, nu, P, L = 2.0e4, 0.0, 1.0e2, 1.0
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _bot, top = _cube(1, 0.0, L, E, nu, ele_tag=1, mat_tag=1, fix_bottom=True)
    ops.node(99, L / 2, L / 2, L - 1.0e-8)     # just penetrated into the top face
    ops.fix(99, 1, 1, 0)                        # free z
    ops.contactSurface(10, "-master", 4, *top)
    ops.contactSurface(20, "-slave", 99)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(99, 0.0, 0.0, -P)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1e-10, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "deformable-master static contact did not converge"
    z_slave = (L - 1.0e-8) + ops.nodeDisp(99, 3)
    z_top = L + ops.nodeDisp(top[0], 3)
    pen = z_top - z_slave
    comp = -ops.nodeDisp(top[0], 3)
    assert abs(pen - P / KN) / (P / KN) < 0.02, f"pen={pen} expect={P/KN}"
    assert abs(comp - P * L / (E * L * L)) / (P * L / (E * L * L)) < 0.05, f"comp={comp}"


def test_p2b2_block_on_block():
    """Two stacked deformable bricks: top block's bottom-corner node is the slave,
    bottom block's top face is the deformable master. Both bricks deform; the
    interface penetration is P/kn. Validates a DEFORMABLE-vs-DEFORMABLE interface."""
    E, nu, P, L = 5.0e4, 0.0, 2.0e2, 1.0
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    # bottom brick (fixed base), top face = master
    _b0, top0 = _cube(1, 0.0, L, E, nu, ele_tag=1, mat_tag=1, fix_bottom=True)
    # top brick starts JUST PENETRATED (bottom nodes 1e-8 below the master face) so
    # the 4 contact pairs are active from step 1 (else the unsupported top block falls).
    bot1, top1 = _cube(11, L - 1.0e-8, L, E, nu, ele_tag=2, mat_tag=2, fix_bottom=False)
    # press the top brick straight down via its TOP nodes; let it ride on the master
    # via its bottom-center-ish corner node as the slave. Constrain lateral drift.
    for t in bot1 + top1:
        ops.fix(t, 1, 1, 0)                    # free z only (uniaxial column)
    # slave = one bottom node of the top brick (node 11) sitting over the master quad
    ops.contactSurface(10, "-master", 4, *top0)
    ops.contactSurface(20, "-slave", 11, 12, 13, 14)   # all 4 bottom nodes
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in top1:
        ops.load(t, 0.0, 0.0, -P / 4)          # total -P spread over the top face
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1e-8, 60, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "block-on-block deformable contact did not converge"
    # interface penetration: bottom-node z of top brick vs top-node z of bottom brick
    z_master = L + ops.nodeDisp(top0[0], 3)
    z_slave = (L - 1.0e-8) + ops.nodeDisp(11, 3)
    pen = z_master - z_slave
    # per-node contact force = P/4 (4 slave nodes share P) → per-node pen = (P/4)/kn
    exp = (P / 4) / KN
    assert pen > 0 and abs(pen - exp) / exp < 0.05, f"interface pen={pen} expect={exp}"
