"""ADR-85 T0 -- 2D rigid-plane contact acceptance (gate G-T0(b), plus the
back-compat half of G-T0(c)).

WHAT IS BEING GATED
    The `RIGID_PLANE` LadrunoContactFE mode is already ndm-PARAMETERIZED: its
    gap/residual/tangent loops run `d < ndm` (LadrunoContactFE.cpp:63-68,
    :655-671) and the handler derives ndm PER NODE from that node's coordinates
    (LadrunoContactHandler.cpp:1446-1489).  Unlike the NTS and mortar lanes it
    never calls the 3D-only `ladrunoSurfaceNodesOk` pre-flight
    (LadrunoContactHandler.cpp:100, "contact needs -ndm 3"), so a `-ndm 2` deck
    can reach it TODAY through the zero-padded 9-arg `contactPlane`.

    "Appears wired" is not a capability.  This file turns it into a tested,
    declared one -- and the ADR asks the question explicitly (Risks: "Does the
    T0 rigid-plane acceptance find real 2D bugs?").  If it does, T0 grows; the
    gate stays.

DECK
    One `quad` (upstream FourNodeQuad, PlaneStrain, nu = 0 so the uniaxial
    modulus is exactly E) of side L, x-DOFs fixed on every node so the column is
    uniaxial, riding on the analytical plane y = 0 with outward normal +y.  The
    two bottom nodes are the slave node-set.  Contact is the ONLY thing holding
    the block up in y, which is the point: the shipped 3D idiom seeds the slave
    slightly penetrating so the pair is active from step 1 (a slave held only by
    the penalty is singular while separated).

TOLERANCES ARE RELATIVE, per the ADR's tolerance policy (How/1): the gap is a
difference of coordinates and carries an eps*|x| noise floor, so every bound
below is expressed against a deck quantity (the penalty depth P/kn, the block
height L, the free-fall impact speed) and never as a bare absolute.

COVERAGE MAP -> ADR-85 G-T0
    (b) statics converge / g bounded / reactions balance
            test_2d_rigidplane_static
        explicit drop, bounded, no energy injection
            test_2d_rigidplane_explicit_drop
        SOFT variant   test_2d_rigidplane_soft
        visc variant   test_2d_rigidplane_visc
    (c) 9-arg contactPlane on a 2D surface still works (back-compat)
            test_2d_rigidplane_static (it IS the 9-arg form)
        7-arg 2D ergonomic form
            test_2d_rigidplane_7arg_form  -- EXPECTED TO FAIL until T0-C

The refusal matrix (the rest of G-T0(c)) lives in
tests/test_adr85_contact2d_t0_refusals.py.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# ---- deck constants (every tolerance below is expressed against one of these)
E = 1.0e3          # Young's modulus
NU = 0.0           # nu = 0 => PlaneStrain uniaxial modulus == E exactly
RHO = 1.0          # element-density mass (the quad's lumped mass, 4 equal nodes)
THICK = 1.0
L = 1.0            # block edge == the deck's length scale
KN = 1.0e6         # rigid-plane penalty
SEED = 1.0e-8      # initial penetration (contact active from step 1)
GRAV = 10.0

M_NODE = RHO * L * L * THICK / 4.0     # lumped nodal mass of one square quad
                                       # (FourNodeQuad lumps; rho comes from the
                                       # material when the element's own rho is 0)


# --------------------------------------------------------------------------
# builders
# --------------------------------------------------------------------------
def _block(y_bottom):
    """A single PlaneStrain quad whose bottom edge sits at y = y_bottom.  Nodes
    1,2 = bottom (the slave pair), 3,4 = top.  x fixed on all four (uniaxial
    column); y free everywhere -- only the contact holds the block up."""
    ops.node(1, 0.0, y_bottom)
    ops.node(2, L, y_bottom)
    ops.node(3, L, y_bottom + L)
    ops.node(4, 0.0, y_bottom + L)
    for t in (1, 2, 3, 4):
        ops.fix(t, 1, 0)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU, RHO)
    ops.element("quad", 1, 1, 2, 3, 4, THICK, "PlaneStrain", 1)


def _declare_plane(nine_arg=True, kn=KN, extra=()):
    """Declare the slave node-set + the analytical plane y = 0, normal +y.

    nine_arg=True  -> the zero-padded 3D form that works TODAY:
                          contactPlane tag ssTag  nx ny 0  px py 0  kn
    nine_arg=False -> the 2D ergonomic form that lands in T0-C:
                          contactPlane tag ssTag  nx ny  px py  kn
    Both must mean exactly the same thing; the ADR fixes the 9-arg form as
    PERMANENTLY valid (How/8), so this helper is also the back-compat gate.
    """
    ops.contactSurface(1, "-slave", 1, 2)
    if nine_arg:
        ops.contactPlane(1, 1, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, kn, *extra)
    else:
        ops.contactPlane(1, 1, 0.0, 1.0, 0.0, 0.0, kn, *extra)


def _explicit_setup():
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1.0e-12, 10, 0)
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")


def _gap(tag, y0):
    """Signed gap of a bottom node against the plane y = 0 (negative == in)."""
    return y0 + ops.nodeDisp(tag, 2)


def _cm_vy():
    return sum(ops.nodeVel(t, 2) for t in (1, 2, 3, 4)) / 4.0


# --------------------------------------------------------------------------
# 1. statics
# --------------------------------------------------------------------------
def _static_run(nine_arg, P=1.0e1, nsteps=5):
    """Press the block onto the plane with total load P, ramped over nsteps
    LoadControl steps.  Returns the observables the gate reads."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    _block(-SEED)
    _declare_plane(nine_arg=nine_arg)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (3, 4):
        ops.load(t, 0.0, -P / 2.0)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 40, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    ops.integrator("LoadControl", 1.0 / nsteps)
    for k in range(nsteps):
        assert ops.analyze(1) == 0, f"2D rigid-plane static step {k} did not converge"
    g1, g2 = _gap(1, -SEED), _gap(2, -SEED)
    return dict(
        P=P,
        g1=g1, g2=g2,
        pen1=-g1, pen2=-g2,
        compression=(ops.nodeDisp(1, 2) - ops.nodeDisp(4, 2)),
        uy=[ops.nodeDisp(t, 2) for t in (1, 2, 3, 4)],
        ux=[ops.nodeDisp(t, 1) for t in (1, 2, 3, 4)],
    )


def test_2d_rigidplane_static():
    """G-T0(b) statics + G-T0(c) 9-arg back-compat.

    Three independent claims on one converged solve:
      1. the penalty depth is the analytic P/(2*kn) at each of the two slave
         nodes (the load splits evenly by symmetry);
      2. REACTIONS BALANCE -- the summed contact force kn*(pen1+pen2) equals the
         applied P.  The rigid-plane adapter's force never reaches
         ops.reactions() (a contact FE's force variants are size 0), and
         ladrunoContactForce is fed only by the NTS branch, so the penalty depth
         IS the readout for this lane;
      3. the block still behaves like a block: with nu = 0 the PlaneStrain
         uniaxial modulus is E, so the column compresses P*L/(E*L*THICK).
    """
    r = _static_run(nine_arg=True)
    d_ref = r["P"] / (2.0 * KN)                    # the deck's penalty scale

    # 1. penetration, RELATIVE to the penalty scale
    for pen in (r["pen1"], r["pen2"]):
        assert pen > 0.0, f"contact never engaged (pen={pen:.3e}); gaps={r['g1']:.3e},{r['g2']:.3e}"
        assert abs(pen - d_ref) / d_ref < 1.0e-6, f"pen={pen:.6e} expect={d_ref:.6e}"
    assert abs(r["pen1"] - r["pen2"]) / d_ref < 1.0e-6, "the two slave nodes are not symmetric"

    # gap bound: the block sits ON the plane to within the penalty depth, and
    # nowhere near sinking through it (RELATIVE to the block height).
    for g in (r["g1"], r["g2"]):
        assert g >= -5.0 * d_ref, f"gap {g:.3e} deeper than 5x the penalty depth {d_ref:.3e}"
        assert abs(g) < 1.0e-3 * L, f"gap {g:.3e} is not small against the block height {L}"

    # 2. reaction balance
    react = KN * (r["pen1"] + r["pen2"])
    assert abs(react - r["P"]) / r["P"] < 1.0e-6, f"contact reaction {react} != applied {r['P']}"

    # 3. the block is really there
    comp_ref = r["P"] * L / (E * L * THICK)
    assert abs(r["compression"] - comp_ref) / comp_ref < 1.0e-6, (
        f"uniaxial compression {r['compression']:.6e} != {comp_ref:.6e}")

    # x is fixed everywhere: nothing may drift laterally
    assert max(abs(v) for v in r["ux"]) < 1.0e-14, f"lateral drift on a fixed-x deck: {r['ux']}"


# --------------------------------------------------------------------------
# 2. explicit drop
# --------------------------------------------------------------------------
def _drop(kn=KN, extra=(), dt=5.0e-5, nsteps=3000, gap0=1.0e-2):
    """Drop the block from gap0 under gravity onto the plane, central
    difference.  Returns the history summary the gates read."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    _block(gap0)
    _declare_plane(nine_arg=True, kn=kn, extra=extra)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t in (1, 2, 3, 4):                          # gravity, distributed by nodal mass
        ops.load(t, 0.0, -M_NODE * GRAV)
    _explicit_setup()
    min_gap, max_absu, max_absv, v_rebound, blew = gap0, 0.0, 0.0, 0.0, False
    touched = False
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            blew = True
            break
        g = min(_gap(1, gap0), _gap(2, gap0))
        min_gap = min(min_gap, g)
        max_absu = max(max_absu, max(abs(ops.nodeDisp(t, 2)) for t in (1, 2, 3, 4)))
        max_absv = max(max_absv, max(abs(ops.nodeVel(t, 2)) for t in (1, 2, 3, 4)))
        if g < 0.0:
            touched = True
        # rebound speed = the largest UPWARD centre-of-mass velocity reached
        # once the block has touched and is no longer penetrating.  (Not "when
        # it climbs back above half the drop height": a strongly damped run
        # never gets that high, which would silently degrade the -visc gate to
        # comparing 0.0 against something.)
        vcm = _cm_vy()
        if touched and g >= 0.0 and vcm > v_rebound:
            v_rebound = vcm
        if max_absu > 100.0 * gap0:                 # runaway == diverged
            blew = True
            break
    return dict(blew=blew, touched=touched, min_gap=min_gap, max_absu=max_absu,
                max_absv=max_absv, v_rebound=v_rebound,
                final_gap=min(_gap(1, gap0), _gap(2, gap0)), gap0=gap0, kn=kn)


def test_2d_rigidplane_explicit_drop():
    """G-T0(b) explicit.

    HONEST CLAIM.  A frictionless elastic penalty with no dissipation does not
    "come to rest" -- it bounces forever, and a test that asserted otherwise
    would be asserting a bug.  What IS gated:
      * the run never aborts and never runs away (the divergence signature),
      * the block engages and REBOUNDS (the active set releases -- a stuck
        active set would hold it down),
      * it never penetrates beyond the impact penalty depth v*sqrt(m/kn) by
        more than a factor, and ends above the plane to within that depth,
      * no energy injection: nodal speeds stay bounded by a small multiple of
        the free-fall impact speed.  (A few times, not 1x: on impact the bottom
        nodes decelerate while the top keeps falling, so the internal
        compression wave legitimately overshoots the rigid-body speed.  A
        DIVERGING contact does not overshoot by a factor, it grows without
        bound -- which is what this catches.)
    """
    r = _drop()
    v_impact = math.sqrt(2.0 * GRAV * r["gap0"])
    pen_ref = v_impact * math.sqrt(M_NODE / KN)     # the deck's impact penalty depth

    assert not r["blew"], "the 2D explicit drop diverged or aborted"
    assert r["touched"], f"the block never reached the plane (min gap {r['min_gap']:.3e})"
    assert r["v_rebound"] > 0.1 * v_impact, (
        f"no rebound: v_up={r['v_rebound']:.3e} vs impact {v_impact:.3e} "
        "(active set stuck ON?)")
    assert r["min_gap"] >= -10.0 * pen_ref, (
        f"penetrated {-r['min_gap']:.3e}, more than 10x the impact depth {pen_ref:.3e}")
    assert r["final_gap"] >= -10.0 * pen_ref, (
        f"ended below the plane: gap {r['final_gap']:.3e}")
    assert r["max_absv"] <= 4.0 * v_impact, (
        f"energy injected: max nodal |v| {r['max_absv']:.3e} vs impact {v_impact:.3e}")


# --------------------------------------------------------------------------
# 3. -soft (explicit SOFT=1 Courant-stable penalty) in 2D
# --------------------------------------------------------------------------
KN_STIFF = 1.0e9      # omega*dt = sqrt(kn/m)*dt = 63 >> 2 at DT_STRUCT
DT_STRUCT = 1.0e-3    # the element's own CFL is L/sqrt(E/rho) = 0.032, so this
                      # dt is stable for the BLOCK and unstable only for contact


def _held(kn, soft, dt=DT_STRUCT, nsteps=2000, v0=-0.5):
    """The ADR-39 P4 `_plane_loaded` rig in 2D: the block starts ON the plane
    under a SUSTAINED weight (so contact stays engaged and an unstable contact
    mode compounds) and is kicked downward.  Returns (diverged, max|uy|)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    _block(0.0)
    extra = ("-soft", soft) if soft > 0.0 else ()
    _declare_plane(nine_arg=True, kn=kn, extra=extra)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t in (1, 2, 3, 4):
        ops.load(t, 0.0, -M_NODE * GRAV)
    for t in (1, 2, 3, 4):
        ops.setNodeVel(t, 2, v0, "-commit")
    _explicit_setup()
    maxabs = 0.0
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return True, maxabs                     # NaN / Inf abort
        maxabs = max(maxabs, max(abs(ops.nodeDisp(t, 2)) for t in (1, 2, 3, 4)))
        if maxabs > 1.0 * L:                        # runaway, relative to the block
            return True, maxabs
    return False, maxabs


def test_2d_rigidplane_soft():
    """G-T0(b) SOFT variant.  The ADR-39 P4 headline, in 2D.

    k_soft = SOFSCL*4*m_eff/dt^2 is a DIMENSION-INDEPENDENT Courant argument
    (ADR-85 How/6), and the rigid-plane sizing path it rides -- gapModeInvMass /
    ladrunoNodeMass / ladrunoInvMassProj -- is already size-guarded for 2-DOF
    nodes.  So this test is a VERIFICATION, not a port: the same stiff kn that
    diverges at the structural dt must stay bounded once `-soft` sizes the
    penalty from the (element-density) nodal mass.  dt is identical across both
    legs; only SOFT differs.
    """
    div_stiff, _ = _held(KN_STIFF, 0.0)
    div_soft, maxabs_soft = _held(KN_STIFF, 0.1)
    assert div_stiff, ("a stiff penalty must DIVERGE at the structural dt "
                       "(omega*dt >> 2) -- if it does not, this deck no longer "
                       "gates anything")
    assert not div_soft, f"-soft must stabilize the 2D explicit contact (max|uy|={maxabs_soft:.3e})"
    # and bounded WELL below the runaway cutoff: the response is the soft
    # equilibrium W/(2*k_soft) plus the kick and the block's own elastic
    # ringing v0/omega_block, all of which are ~1e-2 L or smaller here.
    assert maxabs_soft < 0.1 * L, (
        f"soft response {maxabs_soft:.3e} not bounded near equilibrium "
        f"(< {0.1 * L:.3e}); k_soft = SOFSCL*4*m_eff/dt^2 = "
        f"{0.1 * 4.0 * M_NODE / (DT_STRUCT * DT_STRUCT):.3e}")


# --------------------------------------------------------------------------
# 4. -visc (viscous normal stabilization) in 2D
# --------------------------------------------------------------------------
def test_2d_rigidplane_visc():
    """G-T0(b) visc variant.

    The rigid-plane dashpot loops `d < ndm` (LadrunoContactFE.cpp:660-670), so
    2D inherits it for free -- ADR-85 How/6 lists it under "free (verify, don't
    port)".  This is that verification: the same drop with mu_c > 0 must rebound
    SLOWER than the frictionless-elastic one and stay bounded.  mu_c is sized
    from the deck (zeta ~ 0.25 per slave node) so the signal is unambiguous
    rather than a coin flip on numerical noise.
    """
    muc = 0.5 * 2.0 * math.sqrt(KN * M_NODE)        # zeta = 0.25 per slave node
    free = _drop()
    damped = _drop(extra=("-visc", muc))
    assert not free["blew"] and not damped["blew"], "a 2D -visc drop diverged"
    assert free["touched"] and damped["touched"], "a 2D -visc drop never engaged"
    assert free["v_rebound"] > 0.0, "the undamped baseline did not rebound"
    assert damped["v_rebound"] < 0.8 * free["v_rebound"], (
        f"-visc did not dissipate: rebound {damped['v_rebound']:.4e} vs "
        f"undamped {free['v_rebound']:.4e}")
    # a dashpot removes energy; it must never add any
    assert damped["max_absv"] <= free["max_absv"] * 1.05, (
        f"-visc AMPLIFIED the response: {damped['max_absv']:.3e} vs {free['max_absv']:.3e}")


# --------------------------------------------------------------------------
# 5. the 7-arg 2D form -- T0-C
# --------------------------------------------------------------------------
def test_2d_rigidplane_7arg_form():
    """G-T0(c) ergonomics -- EXPECTED TO FAIL UNTIL ADR-85 T0-C LANDS.

    The 7-arg 2D form `contactPlane tag ssTag nx ny px py kn` does not exist
    yet: today's parser demands 9 positional arguments
    (OpenSeesOutputCommands.cpp:1011) and refuses anything shorter, so this test
    currently dies inside the declaration.  That failure IS the pre-T0-C state
    and it is named so the mapping is unmistakable.

    After T0-C the two forms must be the SAME declaration -- the dimension comes
    from the referenced surface's node coordinates, never from interpreter ndm
    (ADR-85 How/8), so a 2D surface accepts both and they cannot disagree.  The
    gate is therefore equality of the whole solved state, not a re-derivation of
    the physics: the 9-arg leg is already pinned by
    test_2d_rigidplane_static above.

    NOT an xfail on purpose.  A non-strict xfail would report XPASS after T0-C
    (green either way -- it would stop gating), and a strict xfail would report
    a FAILURE the moment T0-C succeeds, forcing the T0-C author to edit the
    marker to record success.  A plain failing test says exactly one thing, in
    both directions.
    """
    ref = _static_run(nine_arg=True)
    got = _static_run(nine_arg=False)
    for k in ("pen1", "pen2", "compression"):
        scale = abs(ref[k])
        assert scale > 0.0, f"degenerate reference for {k}"
        assert abs(got[k] - ref[k]) / scale < 1.0e-12, (
            f"7-arg 2D contactPlane disagrees with the 9-arg zero-padded form on "
            f"{k}: {got[k]!r} vs {ref[k]!r}")
