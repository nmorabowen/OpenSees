"""ADR-85 T3 -- 2D mortar `-tie` (mesh-tying) exactness (gate G-T3(c)).

WHAT THIS FILE IS
    Written against `_adr85_t3_design.md` oracle finding A5 and ADR-57/ADR-41's
    shared-node precedent, for the ADR-85 T3 C++.  Two deck defects were found
    and fixed by the gate run (12 failed / 67 passed against the built
    engine, both root-caused by hand-solve / instrumented diagnosis, both
    test-side) -- see the two notes immediately below before reading further.

    A5 (the design note): "the tie patch is exact for fields with a TANGENTIAL
    gradient; a normal-gradient component differs by exactly g.n". The deck
    below therefore drives u = c0 + c1*x (tangential -- x varies ALONG a
    horizontal master) on the master side and checks the FREE slave side
    reproduces it exactly across a non-matching interface -- the 2D analogue
    of the 3D `-tie` patch test (test_adr41_mortar_c4_1.py::
    test_c4_1_patch_test_split_matches_monolithic), but read at the tie's own
    machine-precision limit (via `analyze_augmented`) rather than at the
    O(1/epsTie) "penalty tolerance" that 3D file settles for.

DEFECT 1 (gate rerun) -- `-slave-segments 2` / `-master 2` take a FLAT
STRIDE-2 PAIR LIST, not a bare node-tag list.  Each consecutive PAIR of tags
is one segment, chained head-to-tail: 3 segments over 4 nodes [a,b,c,d] must
be declared [a,b, b,c, c,d] (6 tags), not the bare 4-tag list -- a bare list
silently declares len(tags)//2 DISJOINT segments with a hole in the middle
(even count, no accidental shared node, so the parser never objects).
Confirmed by the gate's instrumented diagnosis: a hand-solve of the HOLED
surface a bare list actually declares reproduced the "wrong" observed forces
exactly.  `_seg_pairs()` below is the fix, applied to every multi-segment
list in this file (the T-junction test already built its slave list from
explicit segment tuples and was unaffected).

DEFECT 2 (gate rerun) -- non-homogeneous `ops.sp` (an imposed nonzero
displacement) is NOT enforced under `constraints("LadrunoContact")`; the
engine warns "non-homogeneous SP ... NOT enforced ... Use the
DisplacementControl integrator" and silently leaves the DOF free.  This file
used to drive the master's tangential-gradient field via `sp` directly.  The
tie couples RELATIVE DISPLACEMENT (r_I = D*u_s - M*u_m, ADR-41 C4.1) rather
than absolute position, so simply pre-displacing the master's GEOMETRY (u=0
by construction) does not work either -- it would drive the tie toward
u_slave = 0, not toward the target field. The fix is `_pin_to_field()`: each
master node gets a STIFF x-direction spring to a fixed ground node PLUS a
balancing nodal force sized so the spring's equilibrium extension is exactly
the target displacement, absent the tie's own (much smaller) pull -- a
load-driven substitute for `sp` that stays fully within ordinary
force/stiffness equilibrium, which `LadrunoContact` handles correctly. See
its docstring for the anchor-stiffness sizing argument.

COVERAGE MAP -> ADR-85 G-T3(c)
    test_2d_tie_patch_tangential_gradient_exact   the headline: EXACT transfer
        of a tangential-gradient field, held-load Uzawa-augmented to
        machine precision (folds in the "ladrunoMortarTieResidual driven below
        tol by analyze_augmented" clause -- see its docstring for why using
        the augmented driver for BOTH claims at once is the honest way to
        state "exact").
    test_2d_tie_shared_node_t_junction            the ADR-41 C4 crux: a slave
        surface whose two segments share a node exactly at a master segment
        boundary; also the facet-declaration-order independence check
        (test_c4_1_shared_node_order_independent's 2D analogue).
    test_2d_tie_residual_epsTie_independent_augmented   the augmented-tie
        convergence gate in isolation, on a deck with REAL competing elastic
        stiffness on the slave side (so a pure-penalty tie leaves a genuine
        O(1/epsTie) residual and only augmentation removes it) -- swept over
        two decades of epsTie, plus the pure-penalty discriminator at the
        softer one.

WHY THE FIRST TWO TESTS USE FREE-FLOATING SLAVE NODES (NO ELEMENTS)
    The tie's own tangent is `epsTie * B~^T B~ (x) I` (ADR-41 C4.1, no clamp,
    no KKT). With NO other stiffness competing on the slave DOFs, the
    assembled equilibrium `epsTie*(B~^T B~)*u_s = epsTie*(B~^T M~)*u_m` is
    epsTie-INDEPENDENT algebraically (epsTie divides out on both sides) --
    which makes exact transfer the ONLY possible converged answer regardless
    of epsTie, and isolates the claim under test (the D/M weighted transfer of
    a tangential-gradient field across a non-matching interface is exact) from
    any O(1/epsTie) penalty residual. `analyze_augmented` is still run (rather
    than skipped as "unnecessary") to stay uniform with the rest of this
    battery and to exercise `ladrunoMortarTieResidual` directly, per the ADR's
    own three-part gate for G-T3(c).

    The third test deliberately does the OPPOSITE -- gives the slave side its
    own competing elastic stiffness -- specifically so a pure-penalty run
    LEAVES a residual and only the augmented loop removes it; that is what a
    discriminating "driven below tol" gate needs, and the first two tests
    cannot supply it (see their own reasoning above).
"""
import os
import sys

import pytest

from _testbed import ops

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "Ladruno_scripts"))
from analyze_augmented import analyze_augmented  # noqa: E402

pytestmark = [pytest.mark.zone_a]

EPSTIE = 1.0e6
AUGTOL = 1.0e-11
W = 2.0
C0 = 1.0e-3
C1 = 5.0e-4
K_ANCHOR = 1.0e14      # see _pin_to_field's docstring for the sizing argument


def _seg_pairs(tags):
    """`-master 2` / `-slave-segments 2` read a FLAT STRIDE-2 PAIR LIST (see
    the module docstring, DEFECT 1). Turns a chain of node tags into the
    correct flat pair list."""
    out = []
    for i in range(len(tags) - 1):
        out += [tags[i], tags[i + 1]]
    return out


def _pin_to_field(nodes, xs, c0, c1, k=K_ANCHOR, dof=1, tsTag=1, patTag=1):
    """Drive each of `nodes`' DOF `dof` to (c0 + c1*x) WITHOUT `ops.sp` (see
    the module docstring, DEFECT 2): a STIFF spring (uniaxial `Elastic`,
    stiffness `k`) from each node to a coincident FIXED ground node, plus a
    balancing nodal force `k*(c0+c1*x)`. In isolation this spring's own
    equilibrium is exactly `u = force/k = c0+c1*x`; the tie then pulls the
    node away from that target by `tie_reaction/k`, so the approximation
    error shrinks as 1/k.

    SIZING.  The tie's own reaction on a driven node is O(epsTie * target)
    (D/M are O(1) projection weights): with epsTie <= 1e6 in this file and
    target ~ 1e-3, that reaction is O(1e3), so `k=1e14` (the default) bounds
    the induced error at O(1e3/1e14) = O(1e-11) absolute, i.e. O(1e-8)
    relative to a 1e-3-scale target -- three-plus orders inside the 1e-6
    relative bound the exactness tests below actually assert, which is the
    margin that makes this a legitimate stand-in for an exact kinematic
    constraint rather than a tuned pass.
    """
    ops.timeSeries("Linear", tsTag)
    ops.pattern("Plain", patTag, tsTag)
    for idx, (tag, x) in enumerate(zip(nodes, xs)):
        gnd = 90000 + idx
        mat = 90000 + idx
        ops.node(gnd, x, 0.0)
        ops.fix(gnd, 1, 1)
        ops.uniaxialMaterial("Elastic", mat, k)
        ops.element("zeroLength", gnd, gnd, tag, "-mat", mat, "-dir", dof)
        target = c0 + c1 * x
        if dof == 1:
            ops.load(tag, k * target, 0.0)
        else:
            ops.load(tag, 0.0, k * target)


def _tie_gradient_model(c0=C0, c1=C1, ntop=3, nbot=2, epsTie=EPSTIE, w=W):
    """Master: nbot segments over [0,w], driven to c0+c1*x via `_pin_to_field`
    (tangential-gradient field, y held at 0 -- oracle A5's horizontal-master
    convention). Slave: ntop segments over the SAME span, non-matching,
    entirely free -- its equilibrium is governed by the tie alone (see the
    module docstring)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    master, mx = [], []
    for i in range(nbot + 1):
        x = w * i / nbot
        ops.node(100 + i, x, 0.0)
        ops.fix(100 + i, 0, 1)          # y pinned at 0; x free (anchor-driven)
        master.append(100 + i)
        mx.append(x)
    slave, sx = [], []
    for j in range(ntop + 1):
        x = w * j / ntop
        ops.node(200 + j, x, 0.0)
        slave.append(200 + j)
        sx.append(x)

    ops.contactSurface(10, "-master", 2, *_seg_pairs(master))
    ops.contactSurface(20, "-slave-segments", 2, *_seg_pairs(slave))
    # -outward required: the interface is FLUSH (coincident y=0 on both
    # surfaces), which degenerates the centroid vote (ADR-85 How/2).
    ops.contact(1, 10, 20, "-mortar", "-tie", "-epsTie", epsTie,
                "-outward", 0.0, 1.0)

    _pin_to_field(master, mx, c0, c1)
    return dict(master=master, slave=slave, ntop=ntop, nbot=nbot, w=w)


def _static_then_augment(augTol=AUGTOL, maxAug=100, tol=1.0e-12, iters=100):
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", tol, iters, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "2D tie: the anchor-driven step did not converge"
    status, n_aug, res = analyze_augmented(
        ops, maxAug=maxAug, augTol=augTol, query=ops.ladrunoMortarTieResidual)
    assert status == 0, f"2D tie augmentation failed to converge (residual {res:.3e})"
    assert res < augTol, f"2D tie residual {res:.3e} !< augTol {augTol:.0e}"
    return res


def test_2d_tie_patch_tangential_gradient_exact():
    """G-T3(c) headline: u = c0 + c1*x (tangential gradient, horizontal master
    -- oracle A5) transfers to the free slave across a non-matching interface,
    once driven to the tie's own machine-precision limit. Tolerance 1e-6
    relative -- see `_pin_to_field`'s docstring for the anchor-approximation
    error budget this leaves three-plus orders of headroom over.
    """
    m = _tie_gradient_model()
    _static_then_augment()

    for j, tag in enumerate(m["slave"]):
        x = W * j / m["ntop"]
        ref = C0 + C1 * x
        ux = ops.nodeDisp(tag, 1)
        uy = ops.nodeDisp(tag, 2)
        assert abs(ux - ref) <= 1.0e-6 * max(abs(ref), 1.0), (
            f"slave node {tag} (x={x:.4f}): ux={ux!r} != c0+c1*x={ref!r}")
        assert abs(uy) <= 1.0e-9, (
            f"slave node {tag} (x={x:.4f}): uy={uy!r} != 0 -- the field has no "
            "normal-direction component, so this should be exact too "
            "(oracle A5: a normal-gradient component would differ by g.n, but "
            "there is none here)")


def test_2d_tie_shared_node_t_junction():
    """The ADR-41 C4 crux, in 2D: a slave surface whose two SEGMENTS share a
    node exactly at a master segment boundary (a "T-junction" between the two
    discretizations). Master breaks at x=1 (nodes at 0,1,2); slave breaks at
    0, 0.5, 1 (SHARED with the master's own boundary), 2 -- so the slave node
    at x=1 is visited by TWO slave facets ([0.5,1] and [1,2]) while sitting
    exactly on the master's own internal node. This is the global-accumulator
    shared-node case (LEDGER_quirks MAJOR-1, ADR-41 C4.1) transplanted onto a
    junction that also touches the MASTER's discretization boundary, which the
    plain non-matching deck above never exercises (no slave node there ever
    coincides with a master node).

    Two claims: (1) the shared node's transferred value still matches the
    target field (order-independent double-visitation does not perturb it),
    and (2) declaring the two slave facets in FORWARD vs REVERSED order gives
    a BIT-IDENTICAL converged state (test_c4_1_shared_node_order_independent's
    own bound, machine precision -- not a tie/anchor tolerance: both legs
    solve the IDENTICAL system, just with the facets listed in a different
    order, so nothing here depends on the anchor's approximation quality).
    """
    def _build(swap):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        mx = [0.0, 1.0, 2.0]
        master = []
        for i, x in enumerate(mx):
            ops.node(100 + i, x, 0.0)
            ops.fix(100 + i, 0, 1)
            master.append(100 + i)
        sx = [0.0, 0.5, 1.0, 2.0]
        slave = []
        for j, x in enumerate(sx):
            ops.node(200 + j, x, 0.0)
            slave.append(200 + j)

        mseg = [master[0], master[1], master[1], master[2]]
        ops.contactSurface(10, "-master", 2, *mseg)
        segs = [(slave[0], slave[1]), (slave[1], slave[2]), (slave[2], slave[3])]
        if swap:
            segs = list(reversed(segs))
        sflat = [n for seg in segs for n in seg]
        ops.contactSurface(20, "-slave-segments", 2, *sflat)
        ops.contact(1, 10, 20, "-mortar", "-tie", "-epsTie", EPSTIE,
                    "-outward", 0.0, 1.0)

        _pin_to_field(master, mx, C0, C1)
        _static_then_augment()
        return {tag: (ops.nodeDisp(tag, 1), ops.nodeDisp(tag, 2))
                for tag in slave}, sx, slave

    fwd, sx, slave = _build(swap=False)
    rev, _sx2, _slave2 = _build(swap=True)

    for j, tag in enumerate(slave):
        x = sx[j]
        ref = C0 + C1 * x
        ux_f, uy_f = fwd[tag]
        assert abs(ux_f - ref) <= 1.0e-6 * max(abs(ref), 1.0), (
            f"T-junction slave node {tag} (x={x:.4f}): ux={ux_f!r} != {ref!r}")
        assert abs(uy_f) <= 1.0e-9

        ux_r, uy_r = rev[tag]
        assert abs(ux_f - ux_r) <= 1.0e-9 * max(abs(ux_f), 1.0), (
            f"facet-declaration order changed node {tag}'s ux: "
            f"forward {ux_f!r} vs reversed {ux_r!r} -- the T-junction / shared "
            "node must be order-independent (the global accumulator lambda_N "
            "already relies on, ADR-41 C4.1)")
        assert abs(uy_f - uy_r) <= 1.0e-9


# --------------------------------------------------------------------------
# the augmented-tie CONVERGENCE gate in isolation (real competing stiffness)
# --------------------------------------------------------------------------
E, NU = 2.0e4, 0.0
THICK = 1.0
PTIP = 20.0


def _split_strip(epsTie, nbot=2, ntop=3, w=2.0, h=1.0):
    """A 2D analogue of test_adr41_mortar_c4_1.py::_split_column: a fixed-base
    PlaneStrain strip (master) with a NON-matching, tension-loaded free strip
    (slave) tied on top. Unlike the gradient decks above, the slave here has
    its OWN elastic stiffness competing with the tie -- so a pure-penalty tie
    leaves a genuine O(1/epsTie) residual, and only augmentation removes it
    (see the module docstring). Purely LOAD-driven (a real applied tension on
    the top edge, no `sp` anywhere) -- unaffected by DEFECT 2."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)

    base, mid = [], []
    for i in range(nbot + 1):
        x = w * i / nbot
        ops.node(100 + i, x, 0.0)
        ops.fix(100 + i, 1, 1)
        ops.node(110 + i, x, h)
        ops.fix(110 + i, 1, 0)
        base.append(100 + i)
        mid.append(110 + i)
    for i in range(nbot):
        ops.element("quad", 300 + i, base[i], base[i + 1], mid[i + 1], mid[i],
                    THICK, "PlaneStrain", 1)

    slave, top = [], []
    for j in range(ntop + 1):
        x = w * j / ntop
        ops.node(200 + j, x, h)
        ops.fix(200 + j, 1, 0)
        ops.node(210 + j, x, 2.0 * h)
        ops.fix(210 + j, 1, 0)
        slave.append(200 + j)
        top.append(210 + j)
    for j in range(ntop):
        ops.element("quad", 400 + j, slave[j], slave[j + 1], top[j + 1], top[j],
                    THICK, "PlaneStrain", 1)

    ops.contactSurface(10, "-master", 2, *_seg_pairs(mid))
    ops.contactSurface(20, "-slave-segments", 2, *_seg_pairs(slave))
    ops.contact(1, 10, 20, "-mortar", "-tie", "-epsTie", epsTie,
                "-outward", 0.0, 1.0)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    le = w / ntop
    trib = [le / 2.0 if j in (0, ntop) else le for j in range(ntop + 1)]
    for t, wi in zip(top, trib):
        ops.load(t, 0.0, +PTIP * wi / w)     # TENSION away from the master

    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 60, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, f"epsTie={epsTie:.0e}: penalty tie step did not converge"


@pytest.mark.parametrize("epsTie", [1.0e6, 1.0e8])
def test_2d_tie_residual_epsTie_independent_augmented(epsTie):
    """G-T3(c): ladrunoMortarTieResidual driven below augTol by
    analyze_augmented, independent of epsTie (two decades), on a deck where
    the slave's OWN elastic stiffness genuinely competes with the tie penalty.
    """
    augTol = 1.0e-11
    _split_strip(epsTie)
    status, n_aug, res = analyze_augmented(
        ops, maxAug=80, augTol=augTol, query=ops.ladrunoMortarTieResidual)
    assert status == 0, f"epsTie={epsTie:.0e}: tie augmentation failed ({res:.3e})"
    assert res < augTol, f"epsTie={epsTie:.0e}: residual {res:.3e} !< augTol"


def test_2d_tie_residual_pure_penalty_does_not_meet_augtol():
    """THE DISCRIMINATOR: at the softer epsTie (1e6) a PURE-PENALTY tie (no
    augmentation) must NOT already satisfy augTol -- otherwise the gate above
    is vacuous. Mirrors test_2d_mortar_alm_pure_penalty_does_not_meet_augtol.
    """
    augTol = 1.0e-11
    _split_strip(1.0e6)
    res = ops.ladrunoMortarTieResidual()
    assert res >= augTol, (
        f"pure-penalty tie residual {res:.3e} already meets augTol {augTol:.0e} "
        "at the softer epsTie -- this deck cannot discriminate the gate above")
