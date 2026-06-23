"""ADR-41 C2.2 — frictionless commit-cycle Uzawa ALM for mortar contact.

The headline accuracy win. Penalty enforcement (C2.1) leaves a residual penetration
g~ = O(1/epsN): accuracy is bought only with an ever-stiffer (ill-conditioned) epsN.
Augmented Lagrange drives g~ -> an epsN-INDEPENDENT tolerance at a FINITE epsN by carrying
a per-GLOBAL-slave-node multiplier lambda_N, updated once per Domain::commit():

    lambda_I <- min(0, lambda_I + epsN * gbar_I),   gbar_I = g~_I^global / a_I^global

lambda_N + the global weighted gap live on LadrunoContactDomain keyed (contactTag, slaveNodeTag);
the per-facet adapters report g~_I^facet/a_I^facet into the node slot and read the global gap
back (oracle proto_c2_alm.py T7 pins that the per-facet sum with one shared global pressure
reproduces the global assembly). No new DOFs, no saddle-point solver: the SOE equation count is
constant across augmentations. The held-load `analyze_augmented` proc (Ladruno_scripts) re-commits
at a zero load increment to drive ||g~||_inf < augTol within a step.

Gates (capstone C2 / ADR-41 P2):
  - penetration g~ -> an epsN-INDEPENDENT tol within maxAug (the headline penalty can't reach);
  - release/reopen -> contact force = 0 and lambda -> 0 (the active set empties);
  - SOE equation count CONSTANT across augmentations (no new DOFs);
  - ALM on a REAL elastic mesh (a LadrunoBrick block): the interface gap closes to ~0 and the
    converged block state is epsN-independent (the elastic coupling the C2 oracle could not model).

The mechanics were validated oracle-first in contact_prototypes/proto_c2_alm.py (25/25).
"""
import os
import sys

import pytest

from _testbed import ops

# the analyze_augmented held-load Uzawa proc ships in Ladruno_scripts/
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "Ladruno_scripts"))
from analyze_augmented import analyze_augmented  # noqa: E402

pytestmark = [pytest.mark.zone_a]

PTOT = 400.0    # total downward load


def _unit_square(tags, z):
    c = [(tags[0], 0.0, 0.0, z), (tags[1], 1.0, 0.0, z),
         (tags[2], 1.0, 1.0, z), (tags[3], 0.0, 1.0, z)]
    for t, x, y, zz in c:
        ops.node(t, x, y, zz)
    return tags


def _static_setup():
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-11, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


def _rigid_master_slave_facet(epsN, augTol=1.0e-9, maxAug=40):
    """One 4-node slave facet (seeded penetrating) on a coincident FIXED master quad, under a
    uniform downward pressure. Returns (penalty_pen, alm_pen, n_aug, alm_z, sysN_before, sysN_after).
    penalty_pen = the C2.1 penalty depth after the load step; alm_* = after held-load augmentation."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _unit_square([1, 2, 3, 4], 0.0)              # master, fixed
    for t in (1, 2, 3, 4):
        ops.fix(t, 1, 1, 1)
    _unit_square([5, 6, 7, 8], -1.0e-4)          # slave, seeded penetrating; x,y fixed, z free
    for t in (5, 6, 7, 8):
        ops.fix(t, 1, 1, 0)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    ops.contact(1, 1, 2, "-mortar", "-epsN", epsN, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):                        # uniform pressure: corner a=1/4 each (A=1)
        ops.load(t, 0.0, 0.0, -PTOT / 4.0)
    _static_setup()
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0, "penalty load step did not converge"
    sysN_before = ops.systemSize()
    penalty_pen = ops.ladrunoMortarPenetration()
    penalty_z = ops.nodeDisp(5, 3)
    status, n_aug, alm_pen = analyze_augmented(ops, maxAug=maxAug, augTol=augTol)
    sysN_after = ops.systemSize()
    assert status == 0, f"augmentation did not converge (status {status}, pen {alm_pen:.2e})"
    return (penalty_pen, penalty_z, alm_pen, n_aug,
            -1.0e-4 + ops.nodeDisp(5, 3), sysN_before, sysN_after)


def test_c2_2_penetration_epsN_independent():
    """THE HEADLINE: penalty penetration is O(1/epsN) (shrinks with stiffer epsN — accuracy
    bought only with conditioning), but Uzawa augmentation drives it to a FIXED augTol at ANY
    epsN, and the converged interface position is epsN-INDEPENDENT (penalty alone cannot reach)."""
    augTol = 1.0e-9
    lo = _rigid_master_slave_facet(epsN=1.0e4, augTol=augTol)
    hi = _rigid_master_slave_facet(epsN=1.0e6, augTol=augTol)

    # penalty alone: penetration ~ P/epsN, so 100x stiffer epsN => ~100x less penetration
    # (it is epsN-DEPENDENT — the very thing ALM removes).
    assert lo[0] > 50.0 * hi[0], (
        f"penalty penetration not O(1/epsN): epsN=1e4 -> {lo[0]:.2e}, 1e6 -> {hi[0]:.2e}")

    # ALM: penetration driven below augTol at BOTH penalties (the win penalty can't reach).
    assert lo[2] < augTol, f"epsN=1e4 ALM penetration {lo[2]:.2e} !< {augTol}"
    assert hi[2] < augTol, f"epsN=1e6 ALM penetration {hi[2]:.2e} !< {augTol}"

    # and the converged interface (slave z) is epsN-INDEPENDENT: both settle to gap ~0 (z~0),
    # within a tol set by augTol — NOT by epsN. The penalty depths differed by ~100x; the ALM
    # depths agree.
    assert abs(lo[4] - hi[4]) < 1.0e-7, (
        f"converged position epsN-dependent: 1e4 -> {lo[4]:.3e}, 1e6 -> {hi[4]:.3e}")
    assert abs(lo[4]) < 1.0e-7, f"ALM did not close the gap (slave z {lo[4]:.3e} !~ 0)"


def test_c2_2_eqn_count_constant_across_augmentations():
    """No new DOFs: lambda_N is Domain-side state, not an equation. The SOE size is identical
    before and after the whole augmentation sweep (the active set is frozen within the sweep,
    so no domainChanged()/re-number is triggered between re-solves)."""
    res = _rigid_master_slave_facet(epsN=1.0e5)
    sysN_before, sysN_after = res[5], res[6]
    assert sysN_before == sysN_after, (
        f"equation count changed across augmentations: {sysN_before} -> {sysN_after}")
    # a frictionless mortar pair adds NO equations beyond the 4 free slave z-DOFs.
    assert sysN_before == 4, f"unexpected system size {sysN_before} (expected 4 free z-DOFs)"


def test_c2_2_release_reopen_zero_force():
    """Release/reopen: a slave in ALM contact that is then unloaded (pulled away) drops out of
    the active set -> contact force = 0 and lambda -> 0 (the min(0,.) clamp opens the node).
    A residual compressive lambda would spuriously pull the slave back; we check it does not."""
    epsN = 1.0e5
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _unit_square([1, 2, 3, 4], 0.0)
    for t in (1, 2, 3, 4):
        ops.fix(t, 1, 1, 1)
    # an elastic spring per slave node (zeroLength to ground) gives z-stiffness so the slave is
    # well-posed when contact is INACTIVE (a contact-only slave is singular once it releases).
    ks = 1.0e3
    _unit_square([5, 6, 7, 8], -1.0e-4)
    ops.uniaxialMaterial("Elastic", 1, ks)
    for i, t in enumerate((5, 6, 7, 8)):
        ops.node(100 + t, *[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0),
                            (1.0, 1.0, 0.0), (0.0, 1.0, 0.0)][i], )
        ops.fix(100 + t, 1, 1, 1)
        ops.fix(t, 1, 1, 0)
        ops.element("zeroLength", 200 + t, 100 + t, t, "-mat", 1, "-dir", 3)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    ops.contact(1, 1, 2, "-mortar", "-epsN", epsN, "-outward", 0.0, 0.0, 1.0)
    # CONSTANT series: the applied load is time-independent, so the second (reversing) step
    # does not re-scale the first pattern — the net load is exactly p1+p2 at every step.
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    # phase 1: push DOWN into the master, augment to close the gap.
    for t in (5, 6, 7, 8):
        ops.load(t, 0.0, 0.0, -PTOT / 4.0)
    _static_setup()
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    status, _, _ = analyze_augmented(ops, maxAug=40, augTol=1.0e-9)
    assert status == 0, "phase-1 augmentation failed"
    assert ops.ladrunoMortarPenetration() < 1.0e-9

    # phase 2: REVERSE the load to pull the slave UP and away (a second pattern that cancels
    # then reverses). After re-solving + a commit, the slave separates: penetration query = 0,
    # and the slave rides UP with the spring (no residual contact pull-back => lambda -> 0).
    ops.pattern("Plain", 2, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, 0.0, 0.0, +2.0 * PTOT / 4.0)   # net load now +P (tension, away from master)
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0, "release step did not converge"
    # a few held-load commits to settle the multiplier release
    ops.integrator("LoadControl", 0.0)
    for _ in range(5):
        assert ops.analyze(1) == 0

    assert ops.ladrunoMortarPenetration() == 0.0, "node still flagged in contact after release"
    # with NO contact force, each slave node rides its spring alone: net load +P/4 per node
    # against stiffness ks => displacement +P/(4*ks). A residual compressive lambda (failed
    # release) would pull the node back down, shrinking this displacement.
    u_free = (PTOT / 4.0) / ks
    u_slave = ops.nodeDisp(5, 3)
    assert abs(u_slave - u_free) < 1.0e-6, (
        f"contact still acting after release: slave u {u_slave:.4e} != free {u_free:.4e}")


def _elastic_block(L, E, nu, z0):
    """8-node LadrunoBrick on [0,L]^2 x [z0, z0+L]. Bottom face (z0) = slave facet nodes 5..8,
    top face loaded. x,y fixed (uniaxial); bottom + top z free. Returns (bottom4, top4)."""
    pts = [(0, 0, z0), (L, 0, z0), (L, L, z0), (0, L, z0),
           (0, 0, z0 + L), (L, 0, z0 + L), (L, L, z0 + L), (0, L, z0 + L)]
    tags = list(range(5, 13))
    for t, (x, y, z) in zip(tags, pts):
        ops.node(t, float(x), float(y), float(z))
        ops.fix(t, 1, 1, 0)                       # uniaxial: x,y fixed, z free
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    ops.element("LadrunoBrick", 1, *tags, 1)
    return tags[:4], tags[4:]


def _elastic_block_alm(epsN, augTol=1.0e-9):
    """A REAL elastic mesh: a LadrunoBrick block, bottom face = mortar slave on a fixed master
    quad at z=0, loaded DOWN on the top face. ALM closes the interface gap (bottom -> z=0); the
    converged state must be epsN-independent (the elastic coupling the C2 oracle couldn't model).
    Returns (penalty_bottom_z, alm_bottom_z, top_disp)."""
    L, E, nu = 1.0, 2.0e4, 0.0
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _unit_square([1, 2, 3, 4], 0.0)              # fixed master quad at z=0
    for t in (1, 2, 3, 4):
        ops.fix(t, 1, 1, 1)
    bottom, top = _elastic_block(L, E, nu, z0=-1.0e-4)   # block bottom seeded just penetrating
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, *bottom)
    ops.contact(1, 1, 2, "-mortar", "-epsN", epsN, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    P = 1.0e-3 * E * L * L                        # 0.1%-modulus contact pressure
    for t in top:
        ops.load(t, 0.0, 0.0, -P / 4.0)
    _static_setup()
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0, "elastic-block penalty step did not converge"
    penalty_bottom = -1.0e-4 + ops.nodeDisp(bottom[0], 3)
    status, _, _ = analyze_augmented(ops, maxAug=60, augTol=augTol)
    assert status == 0, "elastic-block augmentation did not converge"
    alm_bottom = -1.0e-4 + ops.nodeDisp(bottom[0], 3)
    top_disp = ops.nodeDisp(top[0], 3)
    return penalty_bottom, alm_bottom, top_disp


def test_c2_2_real_elastic_mesh_alm():
    """ALM on a deformable LadrunoBrick (the elastic coupling the mortar-mechanics oracle could
    not model — where the deferred Hertz gate lands). The interface gap closes to ~0 under
    augmentation, and the converged interface position is epsN-INDEPENDENT — penalty alone leaves
    an epsN-dependent residual penetration at the contact face."""
    lo = _elastic_block_alm(epsN=1.0e4)
    hi = _elastic_block_alm(epsN=1.0e6)

    # penalty: the bottom face penetrates the master by ~P/(epsN*A), epsN-dependent (z < 0).
    assert lo[0] < -1.0e-6, f"no penalty penetration at epsN=1e4 (bottom z {lo[0]:.2e})"
    assert abs(lo[0]) > 10.0 * abs(hi[0]), (
        f"penalty penetration not O(1/epsN): 1e4 -> {lo[0]:.2e}, 1e6 -> {hi[0]:.2e}")

    # ALM: the interface gap closes to ~0 (bottom -> z=0) regardless of epsN, and the two agree.
    assert abs(lo[1]) < 1.0e-7, f"ALM did not close the gap at epsN=1e4 (bottom z {lo[1]:.3e})"
    assert abs(hi[1]) < 1.0e-7, f"ALM did not close the gap at epsN=1e6 (bottom z {hi[1]:.3e})"
    assert abs(lo[1] - hi[1]) < 1.0e-7, (
        f"converged interface epsN-dependent: 1e4 -> {lo[1]:.3e}, 1e6 -> {hi[1]:.3e}")
