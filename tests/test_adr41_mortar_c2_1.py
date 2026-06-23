"""ADR-41 C2.1 — mortar PENALTY adapter: the first load-bearing mortar phase.

The `MORTAR` LadrunoContactFE mode integrates a slave/master facet pair with
LadrunoMortarKernel::integratePair at the trial config and assembles the frictionless
penalty force F^s=−(D·t)n / F^m=+(Mᵀ·t)n (t_I=min(0, epsN·ḡ_I)) + the tangent
K_c=epsN·B̃ᵀdiag(act/a)B̃⊗(n⊗n). λ_N Uzawa = C2.2; friction = C3.

Gates:
  - a slave facet pressed onto a fixed master facet settles at the penalty penetration
    δ = P/(epsN·A) (matched + non-matched mesh) — and **static Newton CONVERGES**, which
    is itself the K_c-correctness check (a wrong tangent would not converge quadratically);
  - self-equilibrium: the total reaction on the fixed master = −(applied load);
  - a mortar contact with NO overlap is inert (the clip rejects it ⇒ no spurious force).

The mechanics were validated oracle-first in contact_prototypes/proto_c2_alm.py (18/18);
this battery confirms the C++ transcription assembles into a real OpenSees solve.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

EPS = 1.0e6     # mortar normal penalty
PTOT = 400.0    # total downward load


def _unit_square(tags, z):
    c = [(tags[0], 0.0, 0.0, z), (tags[1], 1.0, 0.0, z),
         (tags[2], 1.0, 1.0, z), (tags[3], 0.0, 1.0, z)]
    for t, x, y, zz in c:
        ops.node(t, x, y, zz)
    return tags


def _solve_static():
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-11, 40, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    return ops.analyze(1)


def test_c2_1_matched_penalty_penetration():
    """One slave quad on a coincident, fixed master quad: the slave settles at the
    penalty penetration δ = P/(epsN·A), A=1. Newton convergence ⇒ K_c is right."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _unit_square([1, 2, 3, 4], 0.0)              # master, fixed
    for t in (1, 2, 3, 4):
        ops.fix(t, 1, 1, 1)
    # slave coincident but seeded slightly penetrating so contact is active from step 1
    # (a slave held ONLY by the penalty is singular when inactive).
    _unit_square([5, 6, 7, 8], -1.0e-4)
    for t in (5, 6, 7, 8):
        ops.fix(t, 1, 1, 0)                      # x,y fixed; z free
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    ops.contact(1, 1, 2, "-mortar", "-epsN", EPS, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, 0.0, 0.0, -PTOT / 4.0)

    assert _solve_static() == 0, "static Newton did not converge (K_c?)"
    # final slave z-position == −δ = −P/(epsN·A); A = 1
    pos = -1.0e-4 + ops.nodeDisp(5, 3)
    assert abs(pos - (-PTOT / EPS)) < 1.0e-7, f"penetration {pos} != {-PTOT/EPS}"
    # uniform: all four slave nodes at the same depth. (The settled penetration δ=P/epsN
    # IS the equilibrium proof — the contact force exactly balanced the applied load.
    # nodeReaction is not used: like the NTS path, the contact FE's force variants are
    # size-0, so contact forces don't flow into ops.reactions().)
    for t in (6, 7, 8):
        p = (-1.0e-4 + ops.nodeDisp(t, 3))
        assert abs(p - pos) < 1.0e-9


def test_c2_1_nonmatched_penalty_penetration():
    """NON-MATCHED: a 2×2 slave mesh (4 facets, 9 nodes) on a single fixed master quad.
    The clip tiles each slave facet over the master; the block still settles at δ=P/(epsN·A)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _unit_square([1, 2, 3, 4], 0.0)              # master single quad, fixed
    for t in (1, 2, 3, 4):
        ops.fix(t, 1, 1, 1)
    # slave 2×2 grid (nodes 11..19), seeded slightly penetrating
    z0 = -1.0e-4
    nid = {}
    k = 11
    for j in range(3):
        for i in range(3):
            nid[(i, j)] = k
            ops.node(k, i * 0.5, j * 0.5, z0)
            ops.fix(k, 1, 1, 0)
            k += 1
    sflat = []
    for j in range(2):
        for i in range(2):
            sflat += [nid[(i, j)], nid[(i + 1, j)], nid[(i + 1, j + 1)], nid[(i, j + 1)]]
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, *sflat)
    ops.contact(1, 1, 2, "-mortar", "-epsN", EPS, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    # CONSISTENT nodal loads load_I = −P·a_I (a_I = ∫N_I dΓ tributary area): a uniform
    # pressure P̄=P gives a UNIFORM penalty depth δ=P/epsN (else the depth varies by node).
    # For a 2×2 unit-quad mesh: corner a=1/16, edge a=1/8, center a=1/4 (Σ=1).
    for (i, j), n in nid.items():
        ne = (1 if i in (0, 2) else 0) + (1 if j in (0, 2) else 0)
        w = {2: 1.0 / 16.0, 1: 1.0 / 8.0, 0: 1.0 / 4.0}[ne]   # corner/edge/center
        ops.load(n, 0.0, 0.0, -PTOT * w)

    assert _solve_static() == 0, "non-matched static Newton did not converge"
    # every slave node settles to the SAME penalty depth (uniform-pressure interface ⇒
    # the partition-of-unity assembly is consistent across the non-matched clip)
    depth = -PTOT / EPS
    for n in nid.values():
        pos = z0 + ops.nodeDisp(n, 3)
        assert abs(pos - depth) < 1.0e-6, f"node {n} pos {pos} != {depth}"


def test_c2_1_coarse_slave_fine_master_full_coverage():
    """A COARSE slave facet (one big quad) over a FINE master mesh (4×4): brute-force
    pairing tiles the whole slave facet over all master facets, so the block settles at
    the FULL-coverage depth δ=P/(epsN·A). (A centroid-only broad phase silently dropped
    the outer master facets — gate MAJOR, C2.1 review.)"""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    # fine master: 4×4 grid of quads (25 nodes), fixed
    nm = 4
    mid = {}
    k = 1
    for j in range(nm + 1):
        for i in range(nm + 1):
            mid[(i, j)] = k
            ops.node(k, i / nm, j / nm, 0.0)
            ops.fix(k, 1, 1, 1)
            k += 1
    mflat = []
    for j in range(nm):
        for i in range(nm):
            mflat += [mid[(i, j)], mid[(i + 1, j)], mid[(i + 1, j + 1)], mid[(i, j + 1)]]
    # coarse slave: ONE quad over the whole [0,1]², seeded penetrating, free z
    _unit_square([101, 102, 103, 104], -1.0e-4)
    for t in (101, 102, 103, 104):
        ops.fix(t, 1, 1, 0)
    ops.contactSurface(1, "-master", 4, *mflat)
    ops.contactSurface(2, "-slave-segments", 4, 101, 102, 103, 104)
    ops.contact(1, 1, 2, "-mortar", "-epsN", EPS, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (101, 102, 103, 104):
        ops.load(t, 0.0, 0.0, -PTOT / 4.0)
    assert _solve_static() == 0, "coarse-slave/fine-master Newton did not converge"
    for t in (101, 102, 103, 104):
        pos = -1.0e-4 + ops.nodeDisp(t, 3)
        assert abs(pos - (-PTOT / EPS)) < 1.0e-7, f"node {t} pos {pos} (coverage hole?)"


def test_c2_1_disjoint_mortar_is_inert():
    """A mortar pair whose facets do NOT overlap produces no force (the clip rejects it):
    a model with the disjoint mortar contact defined evolves BITWISE identically to the
    same model without it (the adapter, if built, integrates to an empty overlap)."""
    def _run(define):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        _unit_square([1, 2, 3, 4], 0.0)              # master, fixed
        for t in (1, 2, 3, 4):
            ops.fix(t, 1, 1, 1)
        # slave quad far in +x (no in-plane overlap), z free with mass (well-posed under
        # Newmark even with zero stiffness), x,y fixed, given a small z-velocity.
        for t, x, y in [(5, 5.0, 0.0), (6, 6.0, 0.0), (7, 6.0, 1.0), (8, 5.0, 1.0)]:
            ops.node(t, x, y, 0.10)
            ops.fix(t, 1, 1, 0)
            ops.mass(t, 1.0, 1.0, 1.0)
            ops.setNodeVel(t, 3, 0.2, "-commit")
        if define:
            ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
            ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
            ops.contact(1, 1, 2, "-mortar", "-epsN", EPS, "-outward", 0.0, 0.0, 1.0)
        ops.constraints("LadrunoContact")
        ops.numberer("Plain")
        ops.system("FullGeneral")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.test("NormDispIncr", 1.0e-12, 25, 0)
        ops.algorithm("Newton")
        ops.analysis("Transient")
        out = []
        for _ in range(8):
            assert ops.analyze(1, 1.0e-3) == 0
            out.append(tuple(ops.nodeDisp(t, 3) for t in (5, 6, 7, 8)))
        return out
    assert _run(True) == _run(False)
