"""ADR-41 C4.1 — mortar MESH-TYING: penalty tie force + tangent (a permanent bond).

Mesh-tying welds a slave surface to a non-matching master surface so they deform as one
continuous body (the zero-gap limit of contact). The active set is frozen ON for every slave
node and the FULL 3-vector weighted relative DISPLACEMENT r_I = Σ_J D_IJ u_s,J − Σ_K M_IK u_m,K
is driven to zero (normal AND tangential, no KKT, no clamp, no friction). C4.1 is the penalty
tie (λ_tie≡0): the force t_I = epsTie·(r_I/a_I) via D/−M + the SPD tangent epsTie·B̃ᵀB̃⊗I₃.

Gates (capstone mesh-tying row, penalty phase):
  - a split structure tied across a NON-MATCHING interface == the monolithic solution to a
    penalty tol (the headline patch test — a uniform stress field transmits across the tie);
  - the tie holds in TENSION (an equality bond — contact would separate; the tie does not);
  - static Newton converges (the tie tangent is SPD — no singular DOF; no -consistanttan);
  - SOE equation count constant (no new DOFs);
  - -tie ⊥ friction (-mu/-cohesion/-tauMax refused);
  - SHARED-NODE order-independence (LEDGER_quirks MAJOR-1 resolution): a slave node shared by
    ≥2 tie facets gives an order-INDEPENDENT bond (the global-accumulator treatment λ_N uses).

The mechanics were validated oracle-first in contact_prototypes/proto_c4_mortar_tie.py.
"""
import os
import sys

import pytest

from _testbed import ops

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "Ladruno_scripts"))
from analyze_augmented import analyze_augmented  # noqa: E402

pytestmark = [pytest.mark.zone_a]

E, NU = 2.0e4, 0.0          # uniaxial, nu=0 ⇒ constant-strain bricks reproduce the bar solution
PTOT = 20.0                 # total end load ⇒ sigma = PTOT/A = 20, eps = 1e-3


def _static_setup():
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 60, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


def _split_column(epsTie, swap_facets=False, tension=True, fix_xy=True):
    """A column [0,1]²×[0,2] split at z=1 into a bottom block (1 LadrunoBrick — the MASTER face)
    and a top block (TWO bricks split in x at 0.5 — the SLAVE face, NON-MATCHING + a SHARED edge
    node at x=0.5). Bottom z=0 fixed; top z=2 loaded uniaxially. The tie welds the interface.

    swap_facets: list the two slave facets in reversed order (to probe facet-eval-order
    independence). tension: +z end load (the tie must HOLD — contact would separate). Returns the
    interface + tip z-displacements and the SOE size.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)

    # bottom block (master): nodes 1..4 (z=0, fixed), 5..8 (z=1, master facet)
    bot = [(1, 0, 0, 0), (2, 1, 0, 0), (3, 1, 1, 0), (4, 0, 1, 0),
           (5, 0, 0, 1), (6, 1, 0, 1), (7, 1, 1, 1), (8, 0, 1, 1)]
    for t, x, y, z in bot:
        ops.node(t, float(x), float(y), float(z))
    for t in (1, 2, 3, 4):
        ops.fix(t, 1, 1, 1)
    ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)

    # top block (slave): interface z=1 nodes 9..14 (x=0,0.5,1 × y=0,1), top z=2 nodes 15..20.
    top = [(9, 0, 0, 1), (10, 0.5, 0, 1), (11, 1, 0, 1), (12, 0, 1, 1), (13, 0.5, 1, 1), (14, 1, 1, 1),
           (15, 0, 0, 2), (16, 0.5, 0, 2), (17, 1, 0, 2), (18, 0, 1, 2), (19, 0.5, 1, 2), (20, 1, 1, 2)]
    for t, x, y, z in top:
        ops.node(t, float(x), float(y), float(z))
    # brick A [0,0.5]×[0,1]×[1,2]; brick B [0.5,1]×[0,1]×[1,2]
    ops.element("LadrunoBrick", 2, 9, 10, 13, 12, 15, 16, 19, 18, 1)
    ops.element("LadrunoBrick", 3, 10, 11, 14, 13, 16, 17, 20, 19, 1)

    if fix_xy:                                   # uniaxial: x,y fixed on every free node
        for t in range(5, 21):
            ops.fix(t, 1, 1, 0)

    # MASTER facet (bottom block top face); SLAVE = the two top-block bottom facets (z=1, CCW).
    ops.contactSurface(1, "-master", 4, 5, 6, 7, 8)
    facetA = [9, 10, 13, 12]
    facetB = [10, 11, 14, 13]
    slave = (facetB + facetA) if swap_facets else (facetA + facetB)
    ops.contactSurface(2, "-slave-segments", 4, *slave)
    # -tie: a permanent bond. -outward +z (the slave's allowed half-space; the interface is
    # coincident at z=1 so the auto orientDir degenerates — pass it explicitly, as C2/C3 do).
    ops.contact(1, 1, 2, "-mortar", "-tie", "-epsTie", epsTie, "-outward", 0.0, 0.0, 1.0)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    sgn = 1.0 if tension else -1.0
    # consistent nodal loads for a uniform end traction t=PTOT over the 2-quad top face:
    # corners (15,17,18,20) get t/8, the shared mid-edge nodes (16,19) get t/4. Σ = t.
    for t in (15, 17, 18, 20):
        ops.load(t, 0.0, 0.0, sgn * PTOT / 8.0)
    for t in (16, 19):
        ops.load(t, 0.0, 0.0, sgn * PTOT / 4.0)

    _static_setup()
    ops.integrator("LoadControl", 1.0)
    rc = ops.analyze(1)
    sysN = ops.systemSize()
    if rc != 0:
        return None
    # interface z-disp (a slave node and a master node — the bond should make them equal) + tip.
    return dict(rc=rc, sysN=sysN,
                u_iface_slave=ops.nodeDisp(9, 3), u_iface_master=ops.nodeDisp(5, 3),
                u_tip=ops.nodeDisp(15, 3),
                tie_res=ops.ladrunoMortarTieResidual())


def test_c4_1_patch_test_split_matches_monolithic():
    """THE HEADLINE: a column split across a NON-MATCHING tied interface reproduces the monolithic
    uniaxial bar solution (constant stress σ=P/A transmits exactly) to a penalty tolerance. The tie
    is loaded in TENSION — an equality bond holds it (contact would separate). u(z) = (P/EA)·z ⇒
    interface (z=1) → P/E, tip (z=2) → 2P/E."""
    r = _split_column(epsTie=1.0e8, tension=True)
    assert r is not None, "tied split column did not converge"
    u_iface_exact = PTOT / E          # 1.0e-3
    u_tip_exact = 2.0 * PTOT / E      # 2.0e-3
    # the bond closes: slave and master interface z-disp agree (penalty leaves O(1/epsTie)).
    assert abs(r["u_iface_slave"] - r["u_iface_master"]) < 1.0e-6, (
        f"interface not tied: slave {r['u_iface_slave']:.4e} vs master {r['u_iface_master']:.4e}")
    # the transmitted field == the monolithic bar solution (constant-stress patch test) to a
    # PENALTY tolerance: the penalty tie leaves a residual bond r = O(load/epsTie) ≈ 2e-7 here
    # (epsTie=1e8) — the O(1/epsTie) gap that C4.2 (ALM) then drives to machine precision.
    assert abs(r["u_iface_slave"] - u_iface_exact) < 1.0e-6, (
        f"interface disp {r['u_iface_slave']:.5e} != monolithic {u_iface_exact:.5e} (penalty tol)")
    assert abs(r["u_tip"] - u_tip_exact) < 1.0e-6, (
        f"tip disp {r['u_tip']:.5e} != monolithic {u_tip_exact:.5e} (penalty tol)")


def test_c4_1_spd_newton_converges_no_new_dofs():
    """The tie tangent is SPD (full identity I₃, no active-set, no Csl) ⇒ static Newton converges
    with a symmetric system and no singular DOF; and λ_tie is Domain-side state, so the SOE size is
    just the free structural DOFs (the tie adds NO equations)."""
    r = _split_column(epsTie=1.0e7, tension=True)
    assert r is not None and r["rc"] == 0, "SPD tie Newton did not converge"
    # free z-DOFs: nodes 5..20 (16 nodes) are z-free; x,y fixed. 16 equations, none added by the tie.
    assert r["sysN"] == 16, f"unexpected SOE size {r['sysN']} (expected 16 free z-DOFs; tie adds none)"


def test_c4_1_tie_refuses_friction():
    """A mesh-tie is an EQUALITY bond — it has no friction cone. -tie with any of
    -mu/-cohesion/-tauMax must be refused (returns < 0 / raises), at the command surface."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(1 + i, float(x), float(y), 0.0)
        ops.node(5 + i, float(x), float(y), 0.0)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    for bad in (("-mu", 0.3), ("-cohesion", 1.0), ("-tauMax", 2.0)):
        with pytest.raises(Exception):
            ops.contact(1, 1, 2, "-mortar", "-tie", "-epsTie", 1.0e5, *bad)
    # -tie WITHOUT -mortar is also refused (mesh-tying is a mortar formulation).
    with pytest.raises(Exception):
        ops.contact(2, 1, 2, "-tie", "-epsTie", 1.0e5)


def test_c4_1_shared_node_order_independent():
    """SHARED-NODE regression (LEDGER_quirks MAJOR-1 resolution): the x=0.5 interface nodes are
    shared by BOTH slave tie facets. The tie state uses the order-INDEPENDENT global accumulator
    λ_N already uses (r_I is a linear accumulation, not a return-map output), so listing the two
    slave facets in either order must give a BIT-IDENTICAL converged solution."""
    fwd = _split_column(epsTie=1.0e7, swap_facets=False, tension=True)
    rev = _split_column(epsTie=1.0e7, swap_facets=True, tension=True)
    assert fwd is not None and rev is not None
    # facet-eval order must not change the result (the force is a per-facet sum; λ_tie uses the
    # global accumulator) — agreement to ~machine precision, NOT merely a penalty tolerance.
    assert abs(fwd["u_tip"] - rev["u_tip"]) < 1.0e-12, (
        f"shared-node tie is facet-order-dependent: fwd {fwd['u_tip']:.6e} vs rev {rev['u_tip']:.6e}")
    assert abs(fwd["u_iface_slave"] - rev["u_iface_slave"]) < 1.0e-12


def test_c4_1_compression_also_ties():
    """The bond is an EQUALITY (no clamp): it transmits COMPRESSION too (and identically to
    tension, up to sign) — the symmetric counterpart of the tension headline test."""
    rt = _split_column(epsTie=1.0e8, tension=True)
    rc = _split_column(epsTie=1.0e8, tension=False)
    assert rt is not None and rc is not None
    # compression mirrors tension: u_tip(compression) == −u_tip(tension).
    assert abs(rc["u_tip"] + rt["u_tip"]) < 1.0e-9, (
        f"compression not the mirror of tension: {rc['u_tip']:.5e} vs {-rt['u_tip']:.5e}")
