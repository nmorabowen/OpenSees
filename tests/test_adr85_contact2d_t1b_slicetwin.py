"""ADR-85 T1b -- the SLICE TWIN and the plane-stress headline (gate G-T1b c, d).

WHAT G-T1b(c) ACTUALLY REQUIRES, VERBATIM
    "slice-twin gate, honestly stated: first measure the element-formulation
     floor (2D plane vs one-layer brick with contact replaced by a tied/merged
     mesh -- the floor is a recorded artifact), then require contact
     OBSERVABLES (summed normal force, per-node gaps) to agree within
     floor + 2%."

    The word doing the work is "honestly".  A naive twin test compares a 2D
    contact deck against a one-element-thick 3D contact deck and attributes the
    whole difference to the contact lane.  That is wrong twice over: part of the
    difference is the ELEMENTS (a FourNodeQuad and a stdBrick are not obliged to
    agree to the last bit), and quoting a single number hides which part is
    which.  So this file measures the element difference first:

      PHASE 1 (_measure_floor)  NO contact -- the interface is merged into one
          continuous body.  Whatever the 2D and 3D answers differ by here is the
          floor, and it is PRINTED, not just asserted.
      PHASE 2 (test_slice_twin_contact_observables)  the interface opened, with
          an NTS pair across it.  The contact observables must agree within
          floor + 2%.

    The 5% ceiling on the floor itself is a sanity bound, not the gate: the two
    elements are exact siblings under this deck (see below), so a floor anywhere
    near 5% means the twin is not a twin and phase 2's budget is meaningless.

WHY PHASE 1 AND PHASE 2 DO NOT SHARE AN INTERFACE TOPOLOGY
    Phase 2's master surface must OVERHANG its slaves (the full argument is in
    the PHASE-2 section below: a slave on a master boundary makes the 3D leg
    refuse and the 2D leg accept, so the twin ends up measuring end policy).
    Merging is the one topology that can never have an overhang -- merged means
    the slave nodes ARE master nodes, i.e. exactly on the boundary.  The two
    phases therefore cannot use the same interface, and they do not.

    That is not a loss, because the floor is a property of the ELEMENT PAIR --
    quad versus one-layer brick under the same material, BC pattern and load --
    and not of the interface, which phase 1 does not have at all.  Phase 2's
    budget is floor + 2% and the floor is round-off (see below), so the budget
    is 2% either way.

WHY THE TWO ELEMENTS ARE EXACT SIBLINGS HERE
    An 8-node trilinear brick, one layer thick with u_z fixed on every node and
    a z-symmetric problem, has a z-INDEPENDENT solution.  Its strain energy then
    integrates to (z-extent) x the plane-strain energy of the bilinear quad
    built on the same in-plane nodes, so in exact arithmetic the in-plane
    answers are IDENTICAL, not merely close.  Verified numerically on the paired
    stiffness at BOTH z-extents this file uses (1 and 4): the residual is 1e-15
    relative.  That the identity scales with the z-extent is what lets phase 2's
    slab overhang in z while its 2D twin stays a single quad of thickness 4.

    Both elements are chosen upstream and plain -- `quad` (FourNodeQuad, 2x2
    Gauss) against `stdBrick` (Brick, 2x2x2), rather than the fork's
    LadrunoQuad/LadrunoBrick -- precisely so the measured floor is round-off and
    nothing else.  A large floor is therefore a finding, and it is printed
    either way.

THE kn SCALING IS A CONVENTION, NOT A FUDGE
    kn is a penalty on a NODAL gap, force/length, and carries no area and no
    thickness (ADR-85 How/7).  The 2D interface has 2 slave nodes; its 3D twin
    has 4 (two z-layers of the same two in-plane positions).  Four springs in
    parallel across the same interface must each be half as stiff for the same
    interface compliance, so the twin uses kn3 = kn2/2 and the per-node gaps
    become directly comparable.  Anything else would compare a stiffer interface
    with a softer one and call the difference a 2D bug.  (`auto` on both sides
    would scale itself, but that would make this gate depend on the auto-kn
    sibling, which has its own test in the NTS file.)

    The master overhang does NOT disturb this.  Only the SLAVE side sets the
    spring count, and the slave body is still one element thick in z, so it is
    still 2 nodes against 4 -- the factor stays 2 whatever the master footprint
    is.  The same goes for gap comparability: per-node force halves and kn
    halves, so the gap is unchanged.

WHY PLANE STRESS CANNOT BE A TWIN TEST
    It is the ADR's headline reason for the whole lane: the one-element-thick
    slice is KINEMATICALLY plane strain (u_z = 0 => eps_zz = 0), and plane
    stress lives in the 2D element family's condensation.  There is no 3D twin
    to compare against, by construction.  So test_2d_plane_stress_contact is
    2D-versus-ANALYTIC, and it asserts both halves: the measurement matches the
    plane-stress closed form, and it is separated from the plane-strain one by
    the factors (1+nu) and (1-nu^2) -- i.e. the slice provably could not have
    produced this number.

NOT MARKED slow
    Four decks of at most four elements each, one static step apiece.  The
    conftest slow tier is for minutes-scale physics and requires a MEASURED wall
    time in the docstring and the ledger row; this file has neither because it
    needs neither.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# ---- shared deck constants
E = 1.0e4
NU = 0.25          # nonzero: the twin must exercise the coupled plane-strain
                   # modulus, not a decoupled nu = 0 special case
THICK = 1.0        # 2D thickness of the SLAVE block (its 3D twin is one TZ layer)
TZ = 1.0           # 3D out-of-plane slice thickness of the SLAVE block; also the
                   # slave spacing in z, which is what the master overhang is
                   # measured in
W = 1.0            # slave-block width == the deck's in-plane length scale, and
                   # the slave spacing in x
H = 1.0            # block height
KN = 1.0e6         # 2D interface penalty (3D twin uses KN/2 -- see docstring)
SEED = 1.0e-8      # seeded penetration: the pair is active from step 1
P = 1.0e2          # total applied load on the top face

# Deliberately NON-UNIFORM so u_y varies along x and the elements are exercised
# in shear as well as compression.  A uniform load would make both twins
# reproduce a linear field exactly and the floor would gate nothing.  (It is
# also what tilts the master, which is what exposed the boundary-slave problem
# the PHASE-2 overhang fixes.)
LOADW = (1.0 / 3.0, 2.0 / 3.0)      # left, right; sums to 1

XS = (0.0, W)                        # the SLAVE block's x-stations (one element
                                     # across).  In phase 2 the master slab spans
                                     # much more than this -- see MX0/MX1 below.

# y levels: 0 = base, 1 = master row, 2 = slave row (seeded), 3 = top row.
# Used by the PHASE-1 merged decks only, which share level 1 between the blocks.
YL = {0: 0.0, 1: H, 2: H - SEED, 3: 2.0 * H}


def _static(tol=1.0e-12, iters=40, what="slice twin"):
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", tol, iters, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0, f"{what}: did not converge"


# --------------------------------------------------------------------------
# PHASE-1 (MERGED) decks -- the element-formulation calibration.
#
# These carry NO contact, so the master-overhang rule below does not apply to
# them; merging means the interface nodes ARE shared, which is the one topology
# an overhang can never have.  See the module docstring ("WHY PHASE 1 AND PHASE
# 2 DO NOT SHARE AN INTERFACE TOPOLOGY").
# --------------------------------------------------------------------------
def _t2(ix, lvl):
    return 100 + 10 * lvl + ix


def _t3(ix, lvl, iz):
    return 100 + 100 * iz + 10 * lvl + ix


def _build_2d_merged():
    """Two stacked PlaneStrain quads sharing their interface row -- one
    continuous body.  x is fixed on every node (the twin compares the in-plane
    response of the two element families, and a laterally confined column
    removes every rigid-body mode without needing asymmetric single-node fixes
    that would have to be mirrored exactly in 3D)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    for lvl in (0, 1, 3):
        for ix, x in enumerate(XS):
            t = _t2(ix, lvl)
            ops.node(t, x, YL[lvl])
            ops.fix(t, 1, 1 if lvl == 0 else 0)
    ops.element("quad", 1, _t2(0, 0), _t2(1, 0), _t2(1, 1), _t2(0, 1),
                THICK, "PlaneStrain", 1)
    ops.element("quad", 2, _t2(0, 1), _t2(1, 1), _t2(1, 3), _t2(0, 3),
                THICK, "PlaneStrain", 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for ix, wgt in enumerate(LOADW):
        ops.load(_t2(ix, 3), 0.0, -P * wgt)


def _build_3d_merged():
    """The same body as one layer of stdBrick, thickness TZ in z, u_z fixed on
    every node -- kinematically plane strain.  Each 2D node becomes a z-PAIR, so
    every nodal load is halved."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    for lvl in (0, 1, 3):
        for iz, z in enumerate((0.0, TZ)):
            for ix, x in enumerate(XS):
                t = _t3(ix, lvl, iz)
                ops.node(t, x, YL[lvl], z)
                ops.fix(t, 1, 1 if lvl == 0 else 0, 1)   # x, z always fixed
    for tag, (lo, hi) in ((1, (0, 1)), (2, (1, 3))):
        # bottom z-face CCW, then top z-face CCW (the shipped _cube ordering)
        ops.element("stdBrick", tag,
                    _t3(0, lo, 0), _t3(1, lo, 0), _t3(1, hi, 0), _t3(0, hi, 0),
                    _t3(0, lo, 1), _t3(1, lo, 1), _t3(1, hi, 1), _t3(0, hi, 1),
                    1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for iz in (0, 1):
        for ix, wgt in enumerate(LOADW):
            ops.load(_t3(ix, 3, iz), 0.0, -P * wgt / 2.0, 0.0)


# --------------------------------------------------------------------------
# PHASE-2 (CONTACT) decks -- master OVERHANGS the slave footprint.
#
# WHY THE OVERHANG (the fix for the first gate run's only failure).  The
# original deck made the master surface exactly co-extensive with the slave
# footprint, so every slave sat on a master boundary.  Under the deliberately
# non-uniform load the master tilts; the closest-point projection of a boundary
# slave then slides marginally past the end (offset ~ gap * slope), and the
# SHIPPED 3D NTS lane HARD-REFUSES an out-of-bounds projection at a facet edge.
# That is its contract, not a bug -- it is the reason SOFT=2 exists.  The 2D leg
# meanwhile converged, because T1b added a terminal acceptance window at open
# chain ends.  So the two legs had DIFFERENT end semantics and the deck sat
# exactly on the difference: the twin was measuring end policy, not the contact
# lane.  Overhanging the master on every side is standard NTS deck practice and
# removes end semantics from this gate entirely.  (End behaviour is covered
# separately -- the 2D window is exercised by the patch deck in
# tests/test_adr85_contact2d_t1b_nts.py, whose outermost slaves DO sit on the
# master chain's ends.)
#
# GEOMETRY.  Slave spacing is W in x and TZ in z; the master overhangs by 1.5x
# that on all four sides, so no slave is anywhere near a boundary:
#
#   slab   x in [-1.5W, 2.5W], y in [0,H], z in [-1.5TZ, 2.5TZ]  -- ONE element
#   block  x in [0, W],        y in [H,2H], z in [0, TZ]         -- ONE element
#   slaves land at (xi, eta) = (0.375, 0.375) .. (0.625, 0.625)
#
# ONE master element, not several: a subdivided master would put a slave on a
# shared segment vertex / facet edge, whose ownership is arbitrated in 2D (the
# T1a rule) and NOT in 3D -- the twins would then disagree for a reason that is
# not a 2D bug.
#
# THE TWIN SURVIVES THE CHANGE.  The slab is still ONE element through z, so the
# paired-brick identity still holds -- it scales with the z-extent, verified
# numerically at both tz = 1 and tz = 4 (rel 1e-15).  The 2D twin of a slab with
# a 4*TZ z-extent is simply a quad of thickness 4*TZ, which is why THICK_SLAB
# below is not 1.  kn3 = kn2/2 and per-node gap comparability are untouched: the
# top block is still one element thick, so the interface still has 2 slave nodes
# in 2D against 4 in 3D.
# --------------------------------------------------------------------------
MX0, MX1 = -1.5 * W, 2.5 * W         # master slab x-span (overhang 1.5 W each side)
MZ0, MZ1 = -1.5 * TZ, 2.5 * TZ       # master slab z-span (overhang 1.5 TZ each side)
THICK_SLAB = MZ1 - MZ0               # the 2D twin of that z-extent

# 2D tags: 1,2 slab base | 3,4 slab top (MASTER) | 5,6 slave row | 7,8 top row
C2_BASE, C2_MAST, C2_SLAVE, C2_TOP = (1, 2), (3, 4), (5, 6), (7, 8)


def _build_2d_contact(kn=KN):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    for t, x in zip(C2_BASE, (MX0, MX1)):
        ops.node(t, x, 0.0)
        ops.fix(t, 1, 1)
    for t, x in zip(C2_MAST, (MX0, MX1)):
        ops.node(t, x, H)
        ops.fix(t, 1, 0)
    ops.element("quad", 1, C2_BASE[0], C2_BASE[1], C2_MAST[1], C2_MAST[0],
                THICK_SLAB, "PlaneStrain", 1)
    for t, x in zip(C2_SLAVE, XS):
        ops.node(t, x, H - SEED)
        ops.fix(t, 1, 0)
    for t, x in zip(C2_TOP, XS):
        ops.node(t, x, 2.0 * H)
        ops.fix(t, 1, 0)
    ops.element("quad", 2, C2_SLAVE[0], C2_SLAVE[1], C2_TOP[1], C2_TOP[0],
                THICK, "PlaneStrain", 1)
    ops.contactSurface(10, "-master", 2, *C2_MAST)
    ops.contactSurface(20, "-slave", *C2_SLAVE)
    ops.contact(1, 10, 20, kn, 0.0, 0.0, "-outward", 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t, wgt in zip(C2_TOP, LOADW):
        ops.load(t, 0.0, -P * wgt)


# 3D tags: slab 11..18 (z=MZ0 face then z=MZ1 face), block 21..28
C3_SLAB = tuple(range(11, 19))
C3_BLOCK = tuple(range(21, 29))
C3_MAST = (C3_SLAB[3], C3_SLAB[2], C3_SLAB[6], C3_SLAB[7])      # the y=H face
C3_SLAVE = ((C3_BLOCK[0], C3_BLOCK[4]), (C3_BLOCK[1], C3_BLOCK[5]))  # per x, both z
C3_TOP = ((C3_BLOCK[3], C3_BLOCK[7]), (C3_BLOCK[2], C3_BLOCK[6]))    # per x, both z


def _build_3d_contact(kn=KN / 2.0):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    slab_pts = [(MX0, 0.0, MZ0), (MX1, 0.0, MZ0), (MX1, H, MZ0), (MX0, H, MZ0),
                (MX0, 0.0, MZ1), (MX1, 0.0, MZ1), (MX1, H, MZ1), (MX0, H, MZ1)]
    for t, (x, y, z) in zip(C3_SLAB, slab_pts):
        ops.node(t, x, y, z)
        ops.fix(t, 1, 1 if y == 0.0 else 0, 1)      # x, z always fixed
    ops.element("stdBrick", 1, *C3_SLAB, 1)
    y0, y1 = H - SEED, 2.0 * H
    blk_pts = [(XS[0], y0, 0.0), (XS[1], y0, 0.0), (XS[1], y1, 0.0), (XS[0], y1, 0.0),
               (XS[0], y0, TZ), (XS[1], y0, TZ), (XS[1], y1, TZ), (XS[0], y1, TZ)]
    for t, (x, y, z) in zip(C3_BLOCK, blk_pts):
        ops.node(t, x, y, z)
        ops.fix(t, 1, 0, 1)
    ops.element("stdBrick", 2, *C3_BLOCK, 1)
    ops.contactSurface(10, "-master", 4, *C3_MAST)
    ops.contactSurface(20, "-slave", *[t for pair in C3_SLAVE for t in pair])
    ops.contact(1, 10, 20, kn, 0.0, 0.0, "-outward", 0.0, 1.0, 0.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for pair, wgt in zip(C3_TOP, LOADW):
        for t in pair:
            ops.load(t, 0.0, -P * wgt / 2.0, 0.0)


# --------------------------------------------------------------------------
# master-surface interpolation -- the slaves are no longer over master NODES,
# so the gap has to be read off the deformed master SURFACE at the slave's
# parametric station (this is what the contact kernel itself does).
# --------------------------------------------------------------------------
def _xi(x):
    return (x - MX0) / (MX1 - MX0)


def _eta(z):
    return (z - MZ0) / (MZ1 - MZ0)


def _gap_2d(slave_tag, x):
    s = _xi(x)
    um = (1.0 - s) * ops.nodeDisp(C2_MAST[0], 2) + s * ops.nodeDisp(C2_MAST[1], 2)
    return ops.nodeDisp(slave_tag, 2) - um - SEED


def _gap_3d(slave_tag, x, z):
    s, t = _xi(x), _eta(z)
    # bilinear over the facet corners in their DECLARED order:
    # (MX0,MZ0), (MX1,MZ0), (MX1,MZ1), (MX0,MZ1)
    wts = ((1.0 - s) * (1.0 - t), s * (1.0 - t), s * t, (1.0 - s) * t)
    um = sum(wi * ops.nodeDisp(ti, 2) for wi, ti in zip(wts, C3_MAST))
    return ops.nodeDisp(slave_tag, 2) - um - SEED


# --------------------------------------------------------------------------
# PHASE 1 -- the element-formulation floor
# --------------------------------------------------------------------------
def _measure_floor():
    """Merged (contact-free) 2D vs 3D.  Returns (floor, detail rows).

    Only u_y is compared: x (and z) are fixed on every node, so u_x is
    identically zero and a "relative difference" on it would be 0/0.  The
    normalizer is the deck's largest |u_y|, which is the honest scale for a
    field that legitimately passes through zero at the fixed base.
    """
    _build_2d_merged()
    _static(what="2D merged twin")
    u2 = {(ix, lvl): ops.nodeDisp(_t2(ix, lvl), 2)
          for lvl in (0, 1, 3) for ix in range(len(XS))}

    _build_3d_merged()
    _static(what="3D merged twin")
    u3 = {(ix, lvl, iz): ops.nodeDisp(_t3(ix, lvl, iz), 2)
          for lvl in (0, 1, 3) for iz in (0, 1) for ix in range(len(XS))}

    scale = max(abs(v) for v in u2.values())
    assert scale > 0.0, "the merged twin did not deform -- the deck is inert"
    rows = []
    for (ix, lvl), a in sorted(u2.items()):
        for iz in (0, 1):
            rows.append((ix, lvl, iz, a, u3[(ix, lvl, iz)],
                         abs(u3[(ix, lvl, iz)] - a) / scale))
    return max(r[5] for r in rows), rows


def test_slice_twin_floor():
    """PHASE 1: measure and RECORD the element-formulation floor.

    This is the artifact G-T1b(c) asks for.  It is not a pass/fail claim about
    contact -- there is no contact in either deck -- it is the number phase 2
    is allowed to spend before its own budget starts.

    The 5% ceiling is a sanity bound on the twin itself.  Under this deck the
    two elements are exact siblings in exact arithmetic (module docstring), so
    the expected floor is round-off; a percent-scale floor means one of the two
    elements is not the plain isoparametric element it is assumed to be, and
    phase 2 would then be measuring that instead of the contact lane.
    """
    floor, rows = _measure_floor()
    print(f"\n[slice-twin FLOOR] max relative u_y difference = {floor:.3e}")
    for ix, lvl, iz, a, b, rel in rows:
        print(f"  ix={ix} lvl={lvl} iz={iz}  2D {a: .12e}  3D {b: .12e}  rel {rel:.3e}")
    assert floor < 0.05, (
        f"element-formulation floor {floor:.3e} exceeds the 5% sanity ceiling: "
        "a one-layer stdBrick with u_z fixed should reproduce the PlaneStrain "
        "FourNodeQuad exactly, so the twin decks are not the same problem")


def test_slice_twin_contact_observables():
    """PHASE 2: the CONTACT observables agree within floor + 2% (G-T1b(c)).

    Observables, exactly the two the gate names:

      SUMMED INTERFACE NORMAL FORCE.  Both twins must transmit the applied P
        (equilibrium) and must transmit the same P as each other.  The equality
        is trivial once both are right, which is why it is not the only check.
      PER-NODE GAPS.  This is the discriminating one: it compares the interface
        COMPLIANCE and the way the load splits between the two in-plane
        positions.  With the non-uniform top load the two nodes carry different
        forces, so a 2D lane that (say) distributed the slave forces through the
        wrong shape functions would land the same total and the wrong split.
        The 3D twin's two z-layers are summed per in-plane position before
        comparison, which is the same pairing that makes kn3 = kn2/2 the right
        scaling.

    The gaps are geometric, read off the deformed configuration rather than
    inferred from the force.  Since the master now OVERHANGS the slaves (see the
    PHASE-2 section above) the slaves no longer sit over master NODES, so the
    master surface is interpolated at each slave's parametric station -- linear
    in 2D, bilinear in 3D -- exactly as the contact kernel does it.

    PREDICTED VALUES (independent reduced-DOF linear solve of this deck: x fixed
    everywhere makes the y-only reduction exact, with the confined modulus
    lam+2mu on u,y and mu on u,x, plus the NTS penalty and the seeded gap):

        contact forces   47.649169, 52.350831   (total exactly 100 = P)
        gaps            -4.76491693e-05, -5.23508307e-05
        no liftoff      min force 47.65, nowhere near zero

    Both stations carry clearly different force AND gap, which is what makes the
    per-station comparison a real check rather than a symmetric tautology.
    """
    floor, _rows = _measure_floor()
    tol = floor + 0.02

    _build_2d_contact()
    _static(what="2D contact twin")
    f2 = [ops.ladrunoContactForce(t) for t in C2_SLAVE]
    g2 = [_gap_2d(t, x) for t, x in zip(C2_SLAVE, XS)]

    _build_3d_contact()
    _static(what="3D contact twin")
    f3 = [sum(ops.ladrunoContactForce(t) for t in pair) for pair in C3_SLAVE]
    g3 = [[_gap_3d(t, x, z) for t, z in zip(pair, (0.0, TZ))]
          for pair, x in zip(C3_SLAVE, XS)]

    print(f"\n[slice-twin CONTACT] floor {floor:.3e}, budget {tol:.3e}")
    print(f"  forces 2D {f2}\n  forces 3D {f3}")
    print(f"  gaps   2D {g2}\n  gaps   3D {g3}")

    assert sum(f2) > 0.0 and sum(f3) > 0.0, (
        f"an interface transmitted nothing (2D {sum(f2)}, 3D {sum(f3)})")
    for name, tot in (("2D", sum(f2)), ("3D", sum(f3))):
        assert abs(tot - P) / P < 1.0e-6, (
            f"{name} summed interface normal force {tot:.8e} != applied "
            f"{P:.8e} -- equilibrium, not a tolerance question")
    assert abs(sum(f3) - sum(f2)) / sum(f2) <= tol, (
        f"summed force disagrees across the twin: 2D {sum(f2):.8e} vs 3D "
        f"{sum(f3):.8e}, budget {tol:.3e}")

    for ix in range(len(XS)):
        assert f2[ix] > 0.0, (
            f"station ix={ix}: the 2D pair carries no force -- it separated, so "
            "the deck is no longer measuring a two-point interface")
        assert abs(f3[ix] - f2[ix]) / abs(f2[ix]) <= tol, (
            f"station ix={ix}: per-node force 2D {f2[ix]:.8e} vs 3D "
            f"{f3[ix]:.8e} -- the load SPLIT across the interface differs, "
            f"budget {tol:.3e}")

    gscale = max(abs(v) for v in g2)
    assert gscale > 0.0, "the 2D interface has no gap at all"
    for ix in range(len(XS)):
        assert g2[ix] < 0.0, f"station ix={ix}: the 2D pair is not in contact"
        for iz in (0, 1):
            rel = abs(g3[ix][iz] - g2[ix]) / gscale
            assert rel <= tol, (
                f"station ix={ix}, iz={iz}: gap 2D {g2[ix]:.8e} vs 3D "
                f"{g3[ix][iz]:.8e} (rel {rel:.3e} > budget {tol:.3e}) -- the "
                "interface compliance differs between the twins")


# --------------------------------------------------------------------------
# PLANE STRESS -- G-T1b(d), the ADR's headline reason for the lane
# --------------------------------------------------------------------------
def test_2d_plane_stress_contact():
    """A plane-STRESS deck under contact, against the analytic thin-plate
    response.  2D versus ANALYTIC by necessity: there is no slice twin for this,
    which is the whole point (module docstring).

    Deck: two identical PlaneStress blocks of thickness TH, contact between
    them, uniform pressure on top, and the lateral boundary FREE -- one x-fix
    per block, enough to remove the rigid-body translation and no more, so the
    Poisson expansion is unconstrained.  That freedom is what makes the test
    discriminating: it is the lateral strain that separates the two states.

    The exact solution is a uniform uniaxial stress field, whose displacements
    are linear and therefore lie exactly in the bilinear quad's space; the
    penalty interface adds a uniform gap, i.e. a rigid translation of the top
    block, and changes no strain.  So the strain assertions are at 1e-6, not at
    a discretization tolerance.

        sigma = P_ps / (W * TH)              (compressive)
        plane STRESS   eps_yy = -sigma/E          eps_xx = +nu*sigma/E
        plane STRAIN   eps_yy = -(1-nu^2)*sigma/E eps_xx = +nu*(1+nu)*sigma/E

    Each measurement is asserted twice: it MATCHES the plane-stress value, and
    it is SEPARATED from the plane-strain one.  With nu = 0.3 the separations
    are EXACT numbers, not estimates: |1 - 1/(1-nu^2)| = 9.89% on eps_yy and
    |1 - 1/(1+nu)| = 23.08% on eps_xx, and the bounds below (5% and 15%) sit
    under both.  The second half is what says the slice could not have produced
    this measurement.

    The force checks pin ADR-85 How/7's thickness convention from both sides:
    the summed contact force equals the applied load (the 2D element carries TH
    inside its stiffness, so the load is an actual force), and the per-node
    penetration is exactly (P_ps/2)/KN -- kn is force/length on a nodal gap and
    is NOT scaled by thickness, so a kn silently multiplied by TH would show up
    here as a factor of 4 while leaving the force balance intact.
    """
    # TH is deliberately NOT 1.0: a kn silently multiplied by the element
    # thickness would be invisible at TH = 1.  The load is scaled down with it
    # so the uniaxial strain stays ~1% and the linear-elastic comparison is not
    # being made at a deformation nobody would run.
    TH, nu, P_ps = 0.25, 0.3, 25.0
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    for t, (x, y) in zip((1, 2, 3, 4),
                         ((0.0, 0.0), (W, 0.0), (W, H), (0.0, H))):
        ops.node(t, x, y)
    ops.fix(1, 1, 1)                    # base: y on both, x on ONE only
    ops.fix(2, 0, 1)
    ops.element("quad", 1, 1, 2, 3, 4, TH, "PlaneStress", 1)
    for t, (x, y) in zip((5, 6, 7, 8),
                         ((0.0, H - SEED), (W, H - SEED), (W, 2.0 * H), (0.0, 2.0 * H))):
        ops.node(t, x, y)
    ops.fix(5, 1, 0)                    # kills the top block's x rigid mode; the
                                        # two contact points kill y and rotation
    ops.element("quad", 2, 5, 6, 7, 8, TH, "PlaneStress", 1)
    ops.contactSurface(10, "-master", 2, 4, 3)
    ops.contactSurface(20, "-slave", 5, 6)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (7, 8):
        ops.load(t, 0.0, -P_ps / 2.0)
    _static(what="2D plane-stress contact")

    f = [ops.ladrunoContactForce(t) for t in (5, 6)]
    assert sum(f) > 0.0, f"the plane-stress interface transmitted nothing ({f})"
    assert abs(sum(f) - P_ps) / P_ps < 1.0e-6, (
        f"summed contact force {sum(f):.8e} != applied {P_ps:.8e}")
    for t, fi in zip((5, 6), f):
        pen = SEED - (ops.nodeDisp(t, 2) - ops.nodeDisp(4 if t == 5 else 3, 2))
        assert abs(fi - P_ps / 2.0) / (P_ps / 2.0) < 1.0e-6, (
            f"node {t}: contact force {fi:.8e} != P_ps/2 {P_ps / 2.0:.8e} (a uniform "
            "pressure on identical blocks splits evenly)")
        assert abs(pen - (P_ps / 2.0) / KN) / ((P_ps / 2.0) / KN) < 1.0e-6, (
            f"node {t}: penetration {pen:.8e} != (P_ps/2)/kn "
            f"{(P_ps / 2.0) / KN:.8e} -- kn was scaled by the element thickness, "
            "which ADR-85 How/7 says it must NOT be (kn is force/length on a "
            "nodal gap; the thickness lives inside the element stiffness)")

    sigma = P_ps / (W * TH)
    eps_yy = ((ops.nodeDisp(7, 2) + ops.nodeDisp(8, 2)) / 2.0
              - (ops.nodeDisp(5, 2) + ops.nodeDisp(6, 2)) / 2.0) / H
    eps_xx = ((ops.nodeDisp(6, 1) - ops.nodeDisp(5, 1))
              + (ops.nodeDisp(7, 1) - ops.nodeDisp(8, 1))) / (2.0 * W)
    ps_yy, ps_xx = -sigma / E, nu * sigma / E
    pn_yy, pn_xx = -(1.0 - nu * nu) * sigma / E, nu * (1.0 + nu) * sigma / E
    print(f"\n[plane stress] eps_yy {eps_yy:.8e} (stress {ps_yy:.8e}, "
          f"strain {pn_yy:.8e});  eps_xx {eps_xx:.8e} (stress {ps_xx:.8e}, "
          f"strain {pn_xx:.8e})")

    assert abs(eps_yy - ps_yy) / abs(ps_yy) < 1.0e-6, (
        f"eps_yy {eps_yy:.8e} != plane-stress -sigma/E {ps_yy:.8e}")
    assert abs(eps_xx - ps_xx) / abs(ps_xx) < 1.0e-6, (
        f"eps_xx {eps_xx:.8e} != plane-stress nu*sigma/E {ps_xx:.8e}")
    assert abs(eps_yy - pn_yy) / abs(pn_yy) > 0.05, (
        f"eps_yy {eps_yy:.8e} is indistinguishable from the plane-STRAIN value "
        f"{pn_yy:.8e}: the element is not condensing sigma_zz, so the deck is "
        "the slice it was supposed to replace")
    assert abs(eps_xx - pn_xx) / abs(pn_xx) > 0.15, (
        f"eps_xx {eps_xx:.8e} is indistinguishable from the plane-STRAIN value "
        f"{pn_xx:.8e}: the lateral Poisson response is the one-layer slice's, "
        "which is the capability ADR-85 Why says 2D exists to provide")
