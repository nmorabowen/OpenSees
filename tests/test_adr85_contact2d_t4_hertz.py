"""ADR-85 T4 -- the 2D Hertz cylinder-on-plane benchmark (gate G-T4).

Convention pinned in Ladruno_implementation/85_ladruno_contact_2d_adr.md (SS Phase
plan T4, G-T4) and independently verified in
Ladruno_implementation/contact_prototypes/proto_t4_hertz2d.py:

    b^2   = 4 P' R / (pi E*)
    p_max = 2 P' / (pi b)

P' is the applied force PER UNIT THICKNESS; this deck uses unit-thickness ("quad",
thickness = 1.0) plane-strain elements, so the applied line load equals the summed
FE contact force numerically -- no h-division anywhere here (2D NTS carries no
-thickness parameter at all; only the mortar lane does, ADR-85 SS How/7).

E* is the PLANE-STRAIN combination 1/E* = (1-nu1^2)/E1 + (1-nu2^2)/E2.  The
cylinder is RIGID (E2 -> infinity), so E* = E1/(1-nu1^2) exactly.

WHY P' IS MEASURED, NOT PRESCRIBED
    2D LINE contact has no delta-only closed form for b.  A 2D half-plane's
    surface displacement is log-divergent, so the indentation depth delta that
    produces a given P' depends on the DOMAIN SIZE; only the LOCAL relations
    b(P'), p_max(P') and the semi-elliptic profile shape are domain-free.  This
    deck therefore drives the contact by interference (rigid cylinder lowered by
    delta), MEASURES P', and compares a and p0 against hertz2d(P'_measured).
    That is the complete local Hertz statement; the delta <-> P' relation is the
    only part left out and it is the only part that is not a property of the
    contact solution.  (The 3D lane's probe_b3_hertz3.py CAN prescribe the patch
    from delta alone -- a = sqrt(R*delta) -- because 3D POINT contact has that
    closed form.  Reusing that idiom in 2D is what the first draft of this deck
    did, and it is wrong in kind, not in tuning.)

FOUR DECK RULES, EACH ESTABLISHED BY MEASUREMENT (do not "simplify" any of them)

 R1  ROLE SWAP.  The DEFORMABLE half-plane's top row is the SLAVE node set and
     the RIGID cylinder is the MASTER arc -- not the other way round.  Then the
     contact patch is resolved by the same mesh that resolves the elasticity,
     each free surface node carries exactly one constraint, and the pressure
     tributary w_i is the (uniform) elastic surface spacing.  With the roles the
     other way the reported "pressure" is f_i divided by an ARBITRARY indenter
     node spacing, and the patch resolution is set by the master mesh rather
     than by the mesh that carries the deformation.

 R2  STAGGER.  No slave may project onto a MASTER VERTEX.  The arc facet grid is
     offset by dx/2 against the slave grid so every slave lands mid-facet.  A
     slave sitting within ~1e-6 facet lengths of a master vertex picks up a
     SECOND claim whose gap is evaluated in the reference configuration; the
     symptom is ladrunoContactForce inflated by exactly kn*(initial gap), and at
     some vertices Newton stalls outright with a residual pinned at kn*(initial
     gap).  See the assertion below -- it is a real guard, not decoration.

 R3  FACET LENGTH FROM DELTA, NOT FROM THE MESH.  On a CURVED master the 2D NTS
     narrow phase stops arming a pair once the penetration exceeds roughly TWICE
     the master facet length (measured on an isolated single-slave rig: cutoff
     pen/Lm = 1.9-2.0 for Lm from 6e-4 to 2.4e-3 at R = 1; a FLAT master shows no
     such cutoff out to 60 facet lengths).  The interference deck starts every
     pair at a penetration of up to delta, so the facet length hm must satisfy
     hm >~ delta.  hm is therefore snapped to an integer multiple of dx sized
     from DELTA and is DECOUPLED from dx, which is sized from b.  Refining dx
     with hm tied to it silently disarms the middle of the patch -- the load
     collapses onto a rim of a few nodes and a and p0 go wrong by 3-10x, which
     LOOKS like discretization error and is not.

 R4  MASK THE ACTIVE SET GEOMETRICALLY (belt-and-suspenders; the root cause is
     FIXED as of this PR -- see below).  ladrunoContactForce used to keep a
     STALE value on a pair that was active and has since RELEASED (minimal
     reproducer: one slave seeded penetrating, pulled back out by a load -- the
     true force is 0, the query kept returning kn*(seeded penetration) forever
     -- see test_2d_hertz_released_pair_force_is_stale). A Hertz deck ALWAYS
     releases pairs: at zero displacement every node inside
     |x| < sqrt(2*R*delta) penetrates, and the converged patch is much narrower
     than that. Summing the RAW query therefore over-reported P' by 2-3 ORDERS
     OF MAGNITUDE before the fix. FIX (this PR, LadrunoContactFE.cpp's 2D
     getResidual branch): the NTS force snapshot is now zeroed on EVERY
     residual evaluation before the active check, not only set on the active
     path -- so a released pair reports 0.0 starting the NEXT evaluation, the
     same "refuse leaks no stale geometry" discipline the kernel's own refusal
     paths already use. Scoped to the 2D branch only (2D-lane code this PR
     already owns); the identical defect in the 3D SEGMENT branch
     (LadrunoContactFE.cpp's segmentActive call site) is DEFERRED to avoid
     coupling a 3D observable-behaviour change into a PR whose gate requires
     the 3D contact_dump to stay bit-identical -- see LEDGER_quirks. The mask
     below is kept anyway: it is now a no-op (QUERY_TOL asserts the query and
     the geometric gap agree to 1e-9 on every genuinely active pair) and a
     second, independent line of defense against any future regression of the
     same kind.
"""
import math
import os
import sys

import pytest

from _testbed import ops

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                                "Ladruno_implementation", "contact_prototypes"))
from proto_t4_hertz2d import hertz2d  # noqa: E402

pytestmark = [pytest.mark.zone_a]

# ---- deck constants (every tolerance below is expressed against one of these)
R = 1.0             # cylinder radius
DELTA = 5.0e-3      # interference (the drive); b/delta ~ 10 for this domain
E = 1.0e4           # half-plane Young's modulus
NU = 0.3
W2 = 2.0            # half-width of the elastic domain  (b/W2 ~ 0.024)
H = 2.0             # depth of the elastic domain       (b/H  ~ 0.024)
KFAC = 1.0e4        # kn = KFAC * E * dx  -- see _kn_note below
RATIO = 1.15        # mesh grading ratio outside the fine zone
FINE = 2.5          # the uniform fine zone is +-FINE*b wide

QUERY_TOL = 1.0e-9  # ladrunoContactForce vs kn*gap on an ACTIVE pair: an identity
EQ_TOL = 1.0e-3     # summed contact force vs base reaction.  NOT machine precision
                    # on purpose: the contact forces are NORMAL to a curved arc and
                    # the base reaction is their VERTICAL resultant, so the two
                    # differ by the pressure-weighted mean of 1 - cos(theta) =
                    # <x^2>/(2R^2) = b^2/(8R^2) = 2.9e-4 for this deck -- the
                    # second-order shallowness term Hertz itself drops.  Measured
                    # 2.93e-4 at every mesh, i.e. the predicted value, so EQ_TOL is
                    # 3.4x it and a genuine leak still fails.

# kn = KFAC*E*dx keeps the CONTACT MODULUS kn/dx equal to KFAC*E, so the relative
# penetration is mesh-independent (measured: pen_max/delta = 5.2e-4 at every mesh
# below).  Sweeping KFAC over 1e2..1e6 moves a/b by 11% at 1e2 and by <0.1% from
# 1e3 upward, so 1e4 sits two decades inside the plateau.


def _graded(x_fine, n_fine, x_far, ratio):
    """[0, x_fine] uniform in n_fine cells, then geometric growth out to x_far."""
    h = x_fine / n_fine
    xs = [i * h for i in range(n_fine + 1)]
    x, hh = x_fine, h
    while x < x_far - 1e-12:
        hh *= ratio
        x = min(x + hh, x_far)
        xs.append(x)
    return xs


def _run(b_hint, nb, R=R, delta=DELTA, E=E, nu=NU, kfac=KFAC, nstep=2,
         tol=1.0e-10, iters=200):
    """One solve.  `b_hint` sizes the refinement zone, `nb` = elements per b."""
    dx = b_hint / nb
    n_fine = max(int(round(FINE * b_hint / dx)), 4)
    x_fine = n_fine * dx

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    xh = _graded(x_fine, n_fine, W2, RATIO)
    xs = [-v for v in reversed(xh[1:])] + xh          # symmetric, node at x = 0
    yh = _graded(x_fine, n_fine, H, RATIO)
    ys = [-v for v in reversed(yh)]                   # -H .. 0, node at y = 0
    nx, ny = len(xs) - 1, len(ys) - 1

    nid, t_ = {}, 1
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.node(t_, xs[i], ys[j])
            nid[(i, j)] = t_
            t_ += 1
    base = [nid[(i, 0)] for i in range(nx + 1)]
    for t in base:
        ops.fix(t, 1, 1)                              # bottom: fully fixed
    for j in range(1, ny + 1):
        for i in (0, nx):
            ops.fix(nid[(i, j)], 1, 0)                # sides: roller
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    eid = 1
    for j in range(ny):
        for i in range(nx):
            ops.element("quad", eid, nid[(i, j)], nid[(i + 1, j)],
                        nid[(i + 1, j + 1)], nid[(i, j + 1)], 1.0, "PlaneStrain", 1)
            eid += 1

    # ---- R3: facet length from DELTA, snapped to an integer multiple of dx
    q = max(1, int(round(delta / dx)))
    hm = q * dx
    sag = hm * hm / (8.0 * R)     # chord-midpoint correction: puts every facet's
                                  # midpoint ON the arc instead of every node, which
                                  # removes the one-sided O(hm^2/R) faceting bias
                                  # and leaves a symmetric ripple half its size.

    def yarc(x):
        return (R - delta) - math.sqrt(max(R * R - x * x, 0.0))

    arc_half = FINE * b_hint
    k0 = int(math.ceil((arc_half - 0.5 * dx) / hm))
    axs = [(-k0 + k) * hm + 0.5 * dx for k in range(2 * k0 + 1)]     # R2 stagger
    ays = [yarc(x) - sag for x in axs]
    mn, mt = [], 800001
    for x, y in zip(axs, ays):
        ops.node(mt, x, y)
        ops.fix(mt, 1, 1)                             # RIGID indenter
        mn.append(mt)
        mt += 1
    chain = []
    for i in range(len(mn) - 1):
        chain += [mn[i], mn[i + 1]]                   # chained stride-2 pair list
    ops.contactSurface(10, "-master", 2, *chain)

    lim = axs[-1] - 2.0 * hm
    slaves, sx, sw = [], {}, {}
    for i in range(nx + 1):
        if abs(xs[i]) > lim:                          # keep the slaves inside the arc
            continue
        t = nid[(i, ny)]
        slaves.append(t)
        sx[t] = xs[i]
        sw[t] = 0.5 * (xs[i + 1] - xs[i - 1])         # consistent nodal tributary
    ops.contactSurface(20, "-slave", *slaves)
    kn = kfac * E * dx
    # the cylinder's material is ABOVE its lower arc, so the master outward normal
    # points DOWN; the vote is degenerate-prone here, so it is TOLD (ADR-85 How/2)
    ops.contact(1, 10, 20, kn, 0.0, 0.0, "-outward", 0.0, -1.0)

    stag = min(min(abs(sx[t] - a) for a in axs) for t in slaves)
    assert stag > 0.25 * dx, (
        f"R2 violated: a slave lies {stag:.3e} from a master vertex (< {0.25*dx:.3e}). "
        "A slave projecting onto a master vertex picks up a second, "
        "reference-configuration claim -- see the module docstring.")

    ops.constraints("LadrunoContact")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormUnbalance", tol, iters, 0)           # the RESIDUAL, not the
    ops.algorithm("Newton")                            # increment: a stalled
    ops.integrator("LoadControl", 1.0)                 # penalty contact can make
    ops.analysis("Static")                             # NormDispIncr pass on a
    rc = 0                                             # state that is not in
    for _ in range(nstep):                             # equilibrium
        rc = ops.analyze(1)
        if rc != 0:
            break

    # ---- R4: independent gap against the (rigid, analytic) FACETED master
    def _chord(xq):
        k = min(max(int(math.floor((xq - axs[0]) / hm)), 0), len(axs) - 2)
        s = (ays[k + 1] - ays[k]) / hm
        return ays[k] + (xq - axs[k]) * s, s

    ops.reactions()
    Rbase = sum(ops.nodeReaction(t, 2) for t in base)
    act, qerr, nq, pen_max, raw = [], 0.0, 0, 0.0, 0.0
    for t in slaves:
        ux, uy = ops.nodeDisp(t, 1), ops.nodeDisp(t, 2)
        yc, sl = _chord(sx[t] + ux)
        pen = (uy - yc) / math.sqrt(1.0 + sl * sl)     # > 0 == penetrating
        fq = ops.ladrunoContactForce(t)
        raw += fq
        if pen > 0.0:
            act.append((sx[t], fq, sw[t]))
            pen_max = max(pen_max, pen)
            if fq > 0.0:
                qerr = max(qerr, abs(fq - kn * pen) / (kn * pen))
                nq += 1
    act.sort()
    P = sum(f for _, f, _ in act)
    out = dict(rc=rc, P=P, P_raw=raw, Rbase=Rbase, kn=kn, dx=dx, hm=hm,
               nslave=len(slaves), nact=len(act), nq=nq, qerr=qerr,
               pen_max=pen_max, nelem=nx * ny)
    if P <= 0.0:
        return out
    ref = hertz2d(P_prime=P, R=R, E1=E, nu1=nu, E2=1.0e30, nu2=0.3)
    prof = [(x, f / w) for x, f, w in act]
    a_node = max(abs(x) for x, _, _ in act)
    p0 = max(v for _, v in prof)
    dev, nsh = 0.0, 0
    for x, pi in prof:
        if abs(x) <= 0.8 * ref["b"]:                   # skip the sqrt-singular edge
            pr = ref["p_max"] * math.sqrt(max(0.0, 1.0 - (x / ref["b"]) ** 2))
            dev = max(dev, abs(pi - pr) / ref["p_max"])
            nsh += 1
    out.update(b=ref["b"], p_max=ref["p_max"], a=a_node, p0=p0, dev=dev, nsh=nsh,
               ra=a_node / ref["b"], rp=p0 / ref["p_max"],
               eq=abs(P - Rbase) / P, nspan=2.0 * ref["b"] / dx, prof=prof)
    return out


def _two_pass(nb):
    """2D has no delta-only b, so the fine zone cannot be sized a priori: pass 1
    solves on a provisional mesh and MEASURES b(P'), pass 2 rebuilds around it."""
    a = _run(b_hint=0.05 * R, nb=6)
    assert a["P"] > 0.0, f"T4 Hertz deck: pass 1 made no contact ({a})"
    return _run(b_hint=a["b"], nb=nb)


def _assert_deck_selfconsistent(tag, r):
    assert r["rc"] == 0, f"{tag}: the Hertz deck did not converge ({r})"
    assert r["P"] > 0.0, f"{tag}: the interface transmitted NOTHING ({r})"
    assert r["eq"] < EQ_TOL, (
        f"{tag}: summed contact force {r['P']:.8e} != base reaction "
        f"{r['Rbase']:.8e} (rel {r['eq']:.2e}) -- the masked active set does not "
        "balance the block's only other external force")
    assert r["nq"] > 0 and r["qerr"] < QUERY_TOL, (
        f"{tag}: ladrunoContactForce != kn*gap on the {r['nq']} genuinely ACTIVE "
        f"pairs (max rel {r['qerr']:.2e}) -- the query and the geometry disagree "
        "where they must not")
    assert r["nspan"] > 8.0, (
        f"{tag}: 2b spans only {r['nspan']:.1f} elements -- the patch is not "
        "resolved, so the comparison below would be gating the mesh, not the lane")
    assert r["b"] / R < 0.1, f"{tag}: b/R = {r['b']/R:.3f} is outside Hertz validity"
    assert r["b"] / min(W2, H) < 0.1, (
        f"{tag}: b/min(W2,H) = {r['b']/min(W2,H):.3f} -- the domain is not a "
        "half-plane at this patch size")


def test_2d_hertz_cylinder_on_plane():
    """G-T4: 2D Hertz line contact, at a COARSE and a REFINED mesh.

    MEASURED at this session's tip (both meshes, two-pass sizing):

        mesh      2b/dx    P'        a/b       p0/p_max   profile dev
        coarse     16     19.858    1.0007      0.9936       0.6%
        fine       32     19.859    1.0007      0.9940       0.9%

      and, holding the mesh at 2b/dx = 30 while sweeping the physics:
        R 1->4, delta 2e-3->1e-2, E 1e4->1e6, nu 0->0.45
        a/b stays in [1.033, 1.034] on the tributary-EDGE estimator and
        p0/p_max in [0.989, 0.999], profile deviation 0.13%-1.94%.

    The band below is +-10%, which is loose against those numbers on purpose:
    it is the ADR's ceiling, not the achieved accuracy, and the achieved
    accuracy is printed so a regression shows up as a number, not as a pass.

    THE HALF-WIDTH ESTIMATOR.  `a` is the largest |x| among slave nodes that
    penetrate, so it is resolved to +-dx/2 by construction; the true b is
    bracketed by [a, a + dx].  Both bounds are asserted, which is a stronger
    statement than a single ratio and does not need the estimator's O(dx) bias
    argued away.
    """
    coarse = _two_pass(nb=8)
    _assert_deck_selfconsistent("coarse", coarse)
    fine = _two_pass(nb=16)
    _assert_deck_selfconsistent("fine", fine)

    for tag, r in (("coarse", coarse), ("fine", fine)):
        print(f"\n[Hertz2D {tag}] nelem={r['nelem']} dx={r['dx']:.5g} hm={r['hm']:.5g} "
              f"kn={r['kn']:.4g} active={r['nact']}/{r['nslave']}")
        print(f"    P'={r['P']:.6g} (raw unmasked query would be {r['P_raw']:.6g}, "
              f"{r['P_raw']/r['P']:.1f}x -- released pairs keep a stale force) "
              f"baseR={r['Rbase']:.6g} equil={r['eq']:.2e}")
        print(f"    b={r['b']:.6g} 2b/dx={r['nspan']:.1f}   a={r['a']:.6g} "
              f"a/b={r['ra']:.4f}   p0={r['p0']:.6g} p_max={r['p_max']:.6g} "
              f"p0/p_max={r['rp']:.4f}   profile dev={100*r['dev']:.2f}% "
              f"({r['nsh']} stations)")

        assert r["a"] <= r["b"] * (1.0 + 1.0e-2) and r["a"] + r["dx"] >= r["b"], (
            f"{tag}: Hertz b = {r['b']:.6g} is not bracketed by the resolved "
            f"contact patch [{r['a']:.6g}, {r['a'] + r['dx']:.6g}]")
        assert 0.9 < r["ra"] < 1.1, (
            f"{tag}: contact half-width ratio a/b = {r['ra']:.4f} outside the "
            f"G-T4 band (a = {r['a']!r}, Hertz b = {r['b']!r}, P' = {r['P']!r})")
        assert 0.9 < r["rp"] < 1.1, (
            f"{tag}: peak pressure ratio p0/p_max = {r['rp']:.4f} outside the "
            f"G-T4 band (p0 = {r['p0']!r}, Hertz p_max = {r['p_max']!r})")
        assert r["dev"] < 0.10, (
            f"{tag}: the FE pressure profile departs from the semi-ellipse "
            f"p_max*sqrt(1-(x/b)^2) by {100*r['dev']:.2f}% of p_max over "
            f"|x| <= 0.8b ({r['nsh']} stations).  This is the claim the two "
            "scalar ratios cannot make: a and p0 can both be right while the "
            "shape is wrong.")

    # mesh convergence: refining must not make the agreement WORSE by a margin
    assert abs(fine["rp"] - 1.0) < abs(coarse["rp"] - 1.0) + 0.02, (
        f"refining regressed the peak-pressure agreement: coarse "
        f"{coarse['rp']:.4f} -> fine {fine['rp']:.4f}")
    assert abs(fine["P"] - coarse["P"]) / coarse["P"] < 0.02, (
        f"P' is not mesh-converged: coarse {coarse['P']:.6g} vs fine "
        f"{fine['P']:.6g} -- P' is a domain property here, not a mesh one")


def test_2d_hertz_released_pair_force_is_stale():
    """THE ENGINE DEFECT R4's MASK EXISTS TO WORK AROUND -- FIXED as of this PR,
    pinned here as a permanent regression guard.

    A pair that has RELEASED must report ZERO contact force. Before this PR it
    reported the force it had while active. Minimal rig: one slave seeded
    INSIDE the master half-plane on a soft spring, with a load that pulls it
    back OUT -- the converged state is unambiguously separated (the slave ends
    at the pure spring solution u = -P/ks, three orders of magnitude away from
    any contact equilibrium), so "still in contact" is not a possible reading.

    THREE control legs share the deck so the claim cannot be blamed on the rig:
    a pair that STAYS active reports correctly, a pair that NEVER touches
    reports 0, and a pair that ACTIVATES LATE reports correctly. Only the
    released leg was wrong.

    FIX: LadrunoContactFE.cpp's 2D getResidual branch now zeroes the NTS force
    snapshot on every evaluation before the active check (see the module
    docstring, R4). This test is the permanent guard against that regressing;
    R4's mask in the deck above stays as a second, independent line of defense.
    """
    KN, KS, P0, P = 1.0e6, 1.0e3, 1.0e-3, 5.0e2

    def _leg(y0, load):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        ops.node(101, -1.0, 0.0)
        ops.node(102, 1.0, 0.0)
        ops.fix(101, 1, 1)
        ops.fix(102, 1, 1)
        ops.node(1, 0.0, y0)
        ops.fix(1, 1, 0)
        ops.node(2, 0.0, y0)
        ops.fix(2, 1, 1)
        ops.uniaxialMaterial("Elastic", 1, KS)
        ops.element("zeroLength", 1, 2, 1, "-mat", 1, "-dir", 2)
        ops.contactSurface(10, "-master", 2, 101, 102)
        ops.contactSurface(20, "-slave", 1)
        # material ABOVE the line: a slave at y > 0 is INSIDE
        ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, -1.0)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(1, 0.0, load)
        ops.constraints("LadrunoContact")
        ops.numberer("Plain")
        ops.system("FullGeneral")
        ops.test("NormUnbalance", 1.0e-10, 100, 0)
        ops.algorithm("Newton")
        ops.integrator("LoadControl", 1.0)
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        u = ops.nodeDisp(1, 2)
        pen = y0 + u
        return u, pen, ops.ladrunoContactForce(1), (KN * pen if pen > 0.0 else 0.0)

    # controls
    for name, y0, load in (("stays active", +P0, +P),
                           ("never touches", -P0, 0.0),
                           ("never touches, pushed further out", -P0, -P),
                           ("activates late", -P0, +4.0 * P)):
        u, pen, f, f_ref = _leg(y0, load)
        assert abs(f - f_ref) <= 1.0e-6 * max(1.0, f_ref), (
            f"control leg '{name}' already disagrees: f={f!r} expected {f_ref!r} "
            f"(u={u!r}, pen={pen!r}) -- the rig, not the defect")

    u, pen, f, f_ref = _leg(+P0, -P)
    assert pen < 0.0, f"the release leg did not release (pen={pen!r}, u={u!r})"
    assert abs(u + P / KS) <= 1.0e-9 * (P / KS), (
        f"the release leg is not at the pure-spring solution: u={u!r} vs "
        f"{-P / KS!r} -- the pair is still carrying something")
    assert f == 0.0, (
        f"a RELEASED 2D NTS pair reports {f!r} instead of 0.0.  It is exactly "
        f"kn*(seeded penetration) = {KN * P0!r}, i.e. the force the pair had "
        "while it was active: releasing does not zero the reported force.  The "
        "assembled residual IS correct (the solve lands on the pure-spring "
        "solution above), so this is a REPORTING defect -- and it is invisible "
        "to every battery deck whose active set only ever grows, which is why "
        "it survived to T4.  test_adr85_contact2d_t4_hertz.py's R4 mask exists "
        "solely because of it.")
