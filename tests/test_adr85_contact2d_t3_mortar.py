"""ADR-85 T3 -- 2D mortar headline: machine-precision transfer + ALM (gate G-T3 a/b + sign).

WHAT THIS FILE IS
    T1b shipped the 2D NTS lane with a REAL, documented discretization error
    across a non-matching interface (test_adr85_contact2d_t1b_nts.py::
    test_2d_patch_compression measures 4.73% pressure deviation on a 2:5 mesh
    pair) -- that error is exactly what the ADR promises the T3 mortar lane
    removes: the 2D interval integrator on the slave-trace measure, 2-point
    Gauss-exact for straight segments (ADR-85 How/3, oracle CHECK 2/3), makes
    the SAME kind of non-matching interface transfer EXACT rather than merely
    "within a percent-scale band". Written against `_adr85_t3_design.md`
    (the T3 design note) and `85_ladruno_contact_2d_adr.md` (the ADR, esp.
    How/3 and G-T3(a)/(b)). Re-run and iterated against the live built
    engine (dist/bin) after the first gate pass -- see the two GATE RERUN
    notes below (SEED and tolerance-scaling test-side fixes) and the sign
    test's own docstring (a non-coincidental-numbers fix).

COVERAGE MAP -> ADR-85 G-T3
    (a) machine-precision uniform pressure transfer, non-matching 2D interface
            test_2d_mortar_patch_uniform_pressure_machine_precision
            test_2d_mortar_patch_uniform_pressure_far_field       (|x|~5e4 twin)
    (b) analyze_augmented drives penetration below an epsN-INDEPENDENT augTol
            test_2d_mortar_alm_epsN_independent            (3 decades, parametrized)
            test_2d_mortar_alm_pure_penalty_does_not_meet_augtol   (the discriminator)
    sign  penetration (g<0) <=> compression; a separated pair (g>=0) is inert
            test_2d_mortar_penetration_sign_compression_convention

WHY THE MASTER-ROW SETTLEMENT UNIFORMITY IS THE LOAD-BEARING CLAIM, NOT REACTIONS
    Summed reactions equalling the applied load is true of ANY converged model
    with the right boundary conditions -- it says nothing about whether the
    INTERFACE itself transferred the load correctly (a contact lane that
    "assembles but transmits nothing" plus a parallel path carrying everything
    would still balance). The patch deck here is the SAME laterally-confined
    (u_x == 0 everywhere, nu = 0) column construction test_adr85_contact2d_t1b_
    nts.py::_patch_model uses, for the SAME reason: with x fixed on every node,
    the exact elasticity solution is UNIFORM sigma_yy, which a bilinear quad
    reproduces EXACTLY -- so any measured non-uniformity is attributable to the
    INTERFACE alone. Where NTS leaves single-digit-percent non-uniformity (T1b),
    mortar's whole reason to exist is that the SAME measurement is exact to
    numerical noise: the mortar transfer is the CONSISTENT Galerkin nodal load
    for whatever traction the top block delivers (uniform, by the same laterally-
    confined argument, regardless of the top block's OWN mesh), so the master row
    must settle uniformly to the solver's own convergence floor, not to a percent
    band.

TOLERANCE POLICY
    NormDispIncr is driven to 1e-10 (GATE RERUN: originally 1e-13, but the
    far-field twin -- coordinates ~5e4 -- failed to converge there: floating-
    point noise on an absolute residual test scales with the coordinate
    magnitude, ~5e-12 at that scale, so 1e-13 was tighter than representable
    at |x|~5e4. 1e-10 clears that floor with two-plus orders of margin while
    staying far tighter than any assertion below) so the assertions below can
    be gated an order tighter than a basic statics identity would need, and
    still leave headroom under the "machine precision" claim: 1e-9 relative on
    the total-force equilibrium identity (three orders tighter than the 1e-6
    used for the equivalent NTS claim in T1b) and 1e-9 relative on master-row
    settlement uniformity (replacing T1b's 5% BAND with an EXACTNESS bound,
    which is the entire point of this gate). A literal 1e-12 is not used for
    the settlement claim because Newton residual noise compounds across the
    block's own DOFs, not just the interface; 1e-9 is still four-plus orders
    inside "a discretization error", the thing being ruled out.
"""
import os
import sys

import pytest

from _testbed import ops

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "Ladruno_scripts"))
from analyze_augmented import analyze_augmented  # noqa: E402

pytestmark = [pytest.mark.zone_a]

# ---- deck constants (2D quads -- the fork's upstream FourNodeQuad idiom, same
# element the T1b patch test uses, "element quad ... PlaneStrain")
E = 1.0e4
NU = 0.0            # nu = 0 -> the laterally-confined column is uniform sigma_yy,
                    # EXACT for a bilinear quad -- see the module docstring
THICK = 1.0
EPSN = 1.0e6
SEED = 1.0e-8       # seeded penetration: the pair is active from step 1.
                    # KEEP THIS SMALL (gate rerun finding, test-side): the top
                    # block's REFERENCE length is (2H) - (H-SEED) = H+SEED,
                    # not H, so comp_ref = eps*H is off by eps*SEED -- at the
                    # original SEED=1e-4 that is a 1e-4 RELATIVE error, well
                    # above the 1e-6 bound test_2d_mortar_patch_uniform_
                    # pressure_machine_precision asserts on comp_top. 1e-8
                    # (matching test_adr85_contact2d_t1b_nts.py's own SEED)
                    # pushes it to 1e-8 relative, three orders inside that
                    # bound. Confirmed empirically: the deviation scales
                    # EXACTLY linearly with SEED (a probe swept SEED over
                    # four decades and the relative comp_top error tracked it
                    # 1:1), so this is a closed-form-reference bookkeeping
                    # issue, not a physics or engine effect.
W = 2.0
H = 1.0


def _seg_pairs(tags):
    """`-master 2` / `-slave-segments 2` read a FLAT STRIDE-2 PAIR LIST: each
    consecutive PAIR of tags is one segment, chained head-to-tail -- e.g. 3
    segments over 4 nodes [a,b,c,d] must be declared as [a,b, b,c, c,d] (6
    tags), NOT the bare 4-tag node list. A bare node list silently declares
    len(tags)//2 DISJOINT segments with holes in between (the parser only
    requires an even count and no accidental shared node, so it never
    objects) -- this is a DEFECT-1 fix (gate rerun, root-caused by hand-solve
    against the exact discrete solution of the holed surface a bare list
    produces). This helper turns a chain of node tags into the correct flat
    pair list."""
    out = []
    for i in range(len(tags) - 1):
        out += [tags[i], tags[i + 1]]
    return out


def _tributary(n, width):
    """Consistent nodal tributary lengths of an edge meshed into n equal
    elements over `width` (half an element at each end, a whole one interior);
    copied from test_adr85_contact2d_t1b_nts.py::_tributary."""
    le = width / n
    return [le / 2.0 if j in (0, n) else le for j in range(n + 1)]


def _patch_model_mortar(ntop=3, nbot=2, epsN=EPSN, offset=(0.0, 0.0)):
    """Two stacked PlaneStrain blocks, NON-matching interface (3 slave segments
    vs 2 master segments over the same span W), joined by a `-mortar` contact
    instead of T1b's plain NTS penalty. `offset` translates the whole deck (the
    far-field twin) -- every geometric input below is built relative to it, so
    the deck is bit-for-bit the same shape, just relocated.
    """
    ox, oy = offset
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)

    base, mid = [], []
    for i in range(nbot + 1):
        x = ox + W * i / nbot
        ops.node(100 + i, x, oy + 0.0)
        ops.fix(100 + i, 1, 1)
        ops.node(110 + i, x, oy + H)
        ops.fix(110 + i, 1, 0)             # x fixed, y free
        base.append(100 + i)
        mid.append(110 + i)
    for i in range(nbot):
        ops.element("quad", 300 + i, base[i], base[i + 1], mid[i + 1], mid[i],
                    THICK, "PlaneStrain", 1)

    slave, top = [], []
    for j in range(ntop + 1):
        x = ox + W * j / ntop
        ops.node(200 + j, x, oy + H - SEED)
        ops.fix(200 + j, 1, 0)
        ops.node(210 + j, x, oy + 2.0 * H)
        ops.fix(210 + j, 1, 0)
        slave.append(200 + j)
        top.append(210 + j)
    for j in range(ntop):
        ops.element("quad", 400 + j, slave[j], slave[j + 1], top[j + 1], top[j],
                    THICK, "PlaneStrain", 1)

    ops.contactSurface(10, "-master", 2, *_seg_pairs(mid))
    ops.contactSurface(20, "-slave-segments", 2, *_seg_pairs(slave))
    ops.contact(1, 10, 20, "-mortar", "-epsN", epsN, "-outward", 0.0, 1.0)
    return dict(base=base, mid=mid, slave=slave, top=top, ntop=ntop, nbot=nbot)


def _load_uniform(m, P):
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t, w in zip(m["top"], _tributary(m["ntop"], W)):
        ops.load(t, 0.0, -P * w / W)


def _static(nsteps=2, tol=1.0e-10, iters=80):
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", tol, iters, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    for k in range(nsteps):
        assert ops.analyze(1) == 0, f"2D mortar patch: step {k} did not converge"


def _wmean(tags, n):
    w = _tributary(n, W)
    return sum(wi * ops.nodeDisp(t, 2) for t, wi in zip(tags, w)) / W


def _run_patch_and_assert(offset, P=2.0e2):
    m = _patch_model_mortar(ntop=3, nbot=2, offset=offset)
    _load_uniform(m, P)
    _static()

    ops.reactions()
    tot = sum(ops.nodeReaction(t, 2) for t in m["base"])
    assert abs(abs(tot) - P) / P < 1.0e-9, (
        f"summed base reaction {tot:.10e} != applied load {P:.10e} at offset "
        f"{offset}: the 2D mortar interface's equilibrium is off")

    # the load-bearing claim: master-row settlement uniform to near-machine
    # precision (see the module docstring -- this is what "the interface is
    # exact" MEANS on this laterally-confined, nu=0 deck).
    u_mid = [ops.nodeDisp(t, 2) for t in m["mid"]]
    u_mean = _wmean(m["mid"], m["nbot"])
    dev = max(abs(ui - u_mean) / abs(u_mean) for ui in u_mid)
    assert dev < 1.0e-9, (
        f"master-row settlement non-uniformity {dev:.3e} at offset {offset} "
        f"(u_mid={u_mid}, mean={u_mean:.10e}) -- mortar's whole reason to exist "
        "over NTS (T1b's 5% band) is that this is machine-precision uniform")

    # the top block's OWN compression must also match the closed-form uniaxial
    # answer (comp = sigma*H/E) -- otherwise "uniform" could pass on a
    # uniformly-WRONG interface.
    comp_ref = (P / (W * THICK)) * H / E
    comp_top = _wmean(m["slave"], m["ntop"]) - _wmean(m["top"], m["ntop"])
    assert abs(comp_top - comp_ref) / comp_ref < 1.0e-6, (
        f"top block compression {comp_top:.10e} != uniaxial sigma*H/E "
        f"{comp_ref:.10e} at offset {offset}")
    return dev


def test_2d_mortar_patch_uniform_pressure_machine_precision():
    """G-T3(a), the headline: a 3-segment slave over a 2-segment master (the SAME
    non-matching pairing T1b's NTS patch test uses) transmits a uniform pressure
    to near-machine precision -- not the 4.73%-band NTS gets on the analogous
    2:5 mesh in test_adr85_contact2d_t1b_nts.py::test_2d_patch_compression.
    """
    dev = _run_patch_and_assert((0.0, 0.0))
    print(f"\n[2D mortar patch] master-row settlement deviation {dev:.3e} "
          "(NTS's analogous band was 4.73%)")


def test_2d_mortar_patch_uniform_pressure_far_field():
    """The far-from-origin twin (~5e4 offset, ADR-85 How/1's own tolerance
    policy): the SAME deck, relocated, must hold the SAME *relative* bounds --
    the gap `g` is a difference of large coordinates and keeps an absolute noise
    floor, so this is what proves the assertions above are genuinely relative,
    not accidentally tight because the deck sits at the origin.
    """
    _run_patch_and_assert((5.0e4, 5.0e4))


# --------------------------------------------------------------------------
# G-T3(b): analyze_augmented drives penetration below an epsN-INDEPENDENT augTol
# --------------------------------------------------------------------------
def _alm_block(epsN, seed=1.0e-4, P=200.0):
    """A single matched 2D mortar segment pair -- the ADR-41 c2_1 idiom
    transcribed to 2D: one fixed master line, one free slave line, coincident,
    seeded penetrating, uniform load."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(101, 0.0, 0.0)
    ops.node(102, 1.0, 0.0)
    ops.fix(101, 1, 1)
    ops.fix(102, 1, 1)
    ops.node(1, 0.0, -seed)
    ops.node(2, 1.0, -seed)
    ops.fix(1, 1, 0)
    ops.fix(2, 1, 0)
    ops.contactSurface(10, "-master", 2, 101, 102)
    ops.contactSurface(20, "-slave-segments", 2, 1, 2)
    ops.contact(1, 10, 20, "-mortar", "-epsN", epsN, "-outward", 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, -P / 2.0)
    ops.load(2, 0.0, -P / 2.0)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 40, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, f"epsN={epsN:.0e}: penalty step did not converge"


@pytest.mark.parametrize("epsN", [1.0e4, 1.0e6, 1.0e8])
def test_2d_mortar_alm_epsN_independent(epsN):
    """G-T3(b): at a HELD load, held-load Uzawa augmentation (the ADR-41 D1
    `analyze_augmented` recipe, already the shipped 3D idiom) drives the 2D
    mortar penetration below augTol -- independent of which of THREE decades
    of epsN was declared, the epsN-INDEPENDENT clause the ADR promises.
    """
    augTol = 1.0e-9
    _alm_block(epsN)
    status, n_aug, pen = analyze_augmented(ops, maxAug=80, augTol=augTol)
    assert status == 0, (
        f"epsN={epsN:.0e}: within-step 2D mortar augmentation failed to "
        f"converge (measure {pen:.3e})")
    assert pen < augTol, (
        f"epsN={epsN:.0e}: penetration {pen:.3e} !< augTol {augTol:.0e}")


def test_2d_mortar_alm_pure_penalty_does_not_meet_augtol():
    """THE DISCRIMINATOR: at the SOFTEST epsN (1e4), a PURE-PENALTY solve (no
    augmentation at all) must NOT already satisfy augTol -- otherwise the gate
    above is vacuous (a broken/no-op augmentation loop would still pass if the
    penalty alone already met the bar). This proves the ALM loop is doing real
    work, mirroring test_d1_within_step_contact_convergence's
    `penalty_pen > alm_pen` check but stated as its own discriminating claim.
    """
    augTol = 1.0e-9
    _alm_block(1.0e4)
    pen = ops.ladrunoMortarPenetration()
    assert pen >= augTol, (
        f"pure-penalty penetration {pen:.3e} already meets augTol {augTol:.0e} "
        "at the softest epsN -- this deck cannot discriminate the ALM gate "
        "above (raise epsN's floor or lower augTol)")


# --------------------------------------------------------------------------
# sign convention: penetration (g<0) <=> compression
# --------------------------------------------------------------------------
def test_2d_mortar_penetration_sign_compression_convention():
    """A penetrating pair (g<0) engages the contact's own (much stiffer)
    penalty, collapsing the slave's displacement by orders of magnitude
    relative to a separated pair (g>=0, no penetration), which is inert and
    settles at the pure-spring answer. This ties the fork's sign convention
    (penetration <=> g<0, ADR-85 How/1) to an observable rather than trusting
    the internal sign alone.

    DISPLACEMENT, not a MASTER-side reaction (gate rerun finding, test-side):
    an earlier draft compared `ops.nodeReaction()` summed over the bare
    contact-only master nodes 101/102. Probed directly against the live
    engine: with NOTHING else attached to those nodes, `ops.reactions()`
    reads back EXACTLY 0.0 regardless of whether the pair is engaged --
    matching test_adr41_mortar_c2_1.py's own documented finding ("the contact
    FE's force variants are size-0, so contact forces don't flow into
    ops.reactions()"). The SLAVE's own displacement (used throughout the rest
    of this battery, e.g. `pair_2d_now_live` / `mortar_2d_now_live`) is
    exactly the observable the ADR's own "penalty-depth identity" guidance
    (How/7) points at instead.

    NON-COINCIDENTAL NUMBERS (gate rerun, test-side): an earlier draft used
    a SHARED `seed_sign * 1e-3` gap for both legs, with KN/KS/P chosen so
    that the pure-penalty equilibrium gap for the "penetrating" leg happened
    to land almost EXACTLY on that seed (both reduce to P_per_node /
    (epsN/2)) -- making the leg converge with ~0 NET displacement (correct,
    but a false negative: it looks inert when it is actually sitting AT its
    equilibrium already), and the "separated" leg's spring alone (KS too soft
    relative to P) pulled the node all the way THROUGH the master, silently
    re-engaging contact from the other side. Probed directly against the live
    engine to confirm the diagnosis (both effects reproduced and explained),
    then fixed by choosing values with a wide, explicit margin: a tiny
    (1e-8) penetrating seed far below any plausible equilibrium gap, a
    separated gap (2.0) far above what the control spring alone could ever
    close (KS=1e3 under P/2=500 moves at most 0.5), and asserting the
    separated leg against its OWN closed form (confirmed exact, not merely
    "small") rather than assuming zero.
    """
    KN, KS, P = 1.0e6, 1.0e3, 1.0e3

    def _leg(y0):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        ops.node(101, 0.0, 0.0)
        ops.node(102, 1.0, 0.0)
        ops.fix(101, 1, 1)
        ops.fix(102, 1, 1)
        ops.node(1, 0.0, y0)
        ops.node(2, 1.0, y0)
        ops.fix(1, 1, 0)
        ops.fix(2, 1, 0)
        for sn, gnd, x in ((1, 901, 0.0), (2, 902, 1.0)):
            ops.node(gnd, x, y0)
            ops.fix(gnd, 1, 1)
            ops.uniaxialMaterial("Elastic", gnd, KS)
            ops.element("zeroLength", 900 + sn, gnd, sn, "-mat", gnd, "-dir", 2)
        ops.contactSurface(10, "-master", 2, 101, 102)
        ops.contactSurface(20, "-slave-segments", 2, 1, 2)
        ops.contact(1, 10, 20, "-mortar", "-epsN", KN, "-outward", 0.0, 1.0)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(1, 0.0, -P / 2.0)
        ops.load(2, 0.0, -P / 2.0)
        ops.constraints("LadrunoContact")
        ops.numberer("Plain")
        ops.system("FullGeneral")
        ops.test("NormDispIncr", 1.0e-12, 40, 0)
        ops.algorithm("Newton")
        ops.integrator("LoadControl", 1.0)
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        return ops.nodeDisp(1, 2), ops.nodeDisp(2, 2)

    u1_pen, u2_pen = _leg(-1.0e-8)   # slave BELOW master: penetrating, g < 0
    u1_sep, u2_sep = _leg(+2.0)      # slave far ABOVE master: separated, g >= 0

    u_sep_ref = -(P / 2.0) / KS      # separated: the spring alone carries the load
    for name, u in (("u1_sep", u1_sep), ("u2_sep", u2_sep)):
        assert abs(u - u_sep_ref) / abs(u_sep_ref) < 1.0e-9, (
            f"the separated control ({name}) is not the pure spring: {u!r} "
            f"vs {u_sep_ref!r}")
    for name, u in (("u1_pen", u1_pen), ("u2_pen", u2_pen)):
        assert abs(u) < 1.0e-2 * abs(u_sep_ref), (
            f"penetrating pair ({name}) displaced almost as much as the "
            f"separated control ({u!r} vs {u_sep_ref!r}) -- compression did "
            "not engage the (much stiffer) penalty")
