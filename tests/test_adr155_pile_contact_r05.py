"""ADR-155 (pile-contact R0.5) -- three opt-in mortar controls, gated on closed-form oracles.

The pile-soil probe (apeGmsh piles-validation R1) found three mortar-lane gaps:

  N-2  the commit-cycle Uzawa update runs on EVERY commit, so even a LINEAR tied problem depends
       on the number of load steps  ->  `-augment commit|request|never`
  N-1  the brute-force mortar pairing ties a facet to the ANTIPODAL facets of a closed surface
       (one tie around a full cylinder is several times too stiff)  ->  `-maxGap <d>`
  G-9  no initial-gap control: a faceted skin in a faceted hole starts with spurious gaps and
       penetrations, and an interference fit must be built into the geometry
       ->  `-adjust [tol]` and `-gapOffset <g0>`

Gates (the closed forms are pinned build-free in
Ladruno_implementation/contact_prototypes/proto_adr155_r05.py):
  * N-2: a linear penalty tie with `-augment never` gives the SAME tip displacement for 1, 2 and 5
    load steps (to 1e-12 relative) and converges to the monolithic bar as epsTie -> inf (error
    ~ 1/epsTie). The default (`commit`) still drifts -- the finding, pinned.
  * N-2: `request` = pure penalty on physical steps, ALM inside the analyze_augmented bracket;
    `never` = pure penalty even inside the bracket.
  * N-1: one closed-cylinder tie with `-maxGap` == the 4-sector split (the R1 workaround) within
    1e-6 relative; without the guard it is far stiffer.
  * G-9: `-adjust` on a faceted tube inside a faceted ring (non-matching, spurious gaps AND
    penetrations as meshed) starts with EXACTLY zero displacement and zero penetration under zero
    load; the unadjusted interface does not.
  * G-9: `-gapOffset -delta` on two stacked blocks gives the series-spring interference pressure
    sigma = delta / (2L/E + 1/epsN) (penalty) and delta*E/(2L) after held-load augmentation.
  * Command surface: the flags are -mortar-only, gap shifts are refused on -tie, bad values refused.
"""
import math
import os
import sys

import pytest

from _testbed import ops

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "Ladruno_scripts"))
from analyze_augmented import analyze_augmented  # noqa: E402

pytestmark = [pytest.mark.zone_a]


def _static(tol=1.0e-12, maxit=60, system="FullGeneral", test="NormDispIncr"):
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system(system)
    ops.test(test, tol, maxit, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


# ------------------------------------------------------------------------------------------ N-2
E_COL, P_COL = 2.0e4, 20.0          # sigma = 20, u_tip(monolithic) = 2*P/E = 2e-3


def _split_column(epsTie, nsteps, extra=()):
    """The ADR-41 C4 split column: bottom brick [0,1]^2x[0,1] (MASTER face z=1) tied to a top block of
    TWO bricks split at x=0.5 (SLAVE face, non-matching). Uniaxial (x,y fixed, nu=0). End load P in
    `nsteps` equal LoadControl increments. Returns the tip z-displacement."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, E_COL, 0.0)
    bot = [(1, 0, 0, 0), (2, 1, 0, 0), (3, 1, 1, 0), (4, 0, 1, 0),
           (5, 0, 0, 1), (6, 1, 0, 1), (7, 1, 1, 1), (8, 0, 1, 1)]
    top = [(9, 0, 0, 1), (10, .5, 0, 1), (11, 1, 0, 1), (12, 0, 1, 1), (13, .5, 1, 1), (14, 1, 1, 1),
           (15, 0, 0, 2), (16, .5, 0, 2), (17, 1, 0, 2), (18, 0, 1, 2), (19, .5, 1, 2), (20, 1, 1, 2)]
    for t, x, y, z in bot + top:
        ops.node(t, float(x), float(y), float(z))
    for t in (1, 2, 3, 4):
        ops.fix(t, 1, 1, 1)
    for t in range(5, 21):
        ops.fix(t, 1, 1, 0)
    ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.element("LadrunoBrick", 2, 9, 10, 13, 12, 15, 16, 19, 18, 1)
    ops.element("LadrunoBrick", 3, 10, 11, 14, 13, 16, 17, 20, 19, 1)
    ops.contactSurface(1, "-master", 4, 5, 6, 7, 8)
    ops.contactSurface(2, "-slave-segments", 4, 9, 10, 13, 12, 10, 11, 14, 13)
    ops.contact(1, 1, 2, "-mortar", "-tie", "-epsTie", epsTie, "-outward", 0.0, 0.0, 1.0, *extra)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (15, 17, 18, 20):
        ops.load(t, 0.0, 0.0, P_COL / 8.0)
    for t in (16, 19):
        ops.load(t, 0.0, 0.0, P_COL / 4.0)
    _static()
    ops.integrator("LoadControl", 1.0 / nsteps)
    for _ in range(nsteps):
        assert ops.analyze(1) == 0
    return ops.nodeDisp(15, 3)


def test_n2_never_linear_tie_is_step_count_independent():
    """`-augment never`: a linear penalty tie is a LINEAR problem, so 1, 2 and 5 load steps give the
    same answer (to Newton round-off). The shipped default augments at every commit, so the same
    model drifts with the step count (R1 N-2: 0.7395 vs 0.7131 mm on the pile)."""
    u = [_split_column(1.0e5, n, ("-augment", "never")) for n in (1, 2, 5)]
    for ui in u[1:]:
        assert abs(ui - u[0]) <= 1.0e-12 * abs(u[0]), f"never: step-count dependent {u}"
    d = [_split_column(1.0e5, n) for n in (1, 2, 5)]
    assert d[0] == pytest.approx(u[0], rel=1e-14)       # 1 step: nothing to augment before it
    # the column is exactly 1-D (nu = 0), so the default's drift is the series-spring Uzawa recursion
    # of proto_adr155_r05.py T2 to round-off: 2.2e-3 / 2.1e-3 / 2.04e-3 for 1 / 2 / 5 steps.
    assert d == pytest.approx([2.2e-3, 2.1e-3, 2.04e-3], rel=1e-10), d
    assert u == pytest.approx([2.2e-3] * 3, rel=1e-12), u


def test_n2_never_converges_to_the_tie_limit():
    """Pure penalty: the tip error against the monolithic bar (2P/E) is O(1/epsTie) -- one decade
    of epsTie removes one decade of error -- and vanishes in the limit."""
    exact = 2.0 * P_COL / E_COL
    err = [abs(_split_column(eps, 1, ("-augment", "never")) - exact) for eps in (1e5, 1e6, 1e7, 1e8)]
    for a, b in zip(err, err[1:]):
        assert 8.0 < a / b < 12.0, f"penalty error not O(1/epsTie): {err}"
    assert err[-1] < 2.0e-4 * exact          # measured 1.0e-4 relative at epsTie = 1e8


# ------------------------------------------------------------------------------------ blocks
E_BLK, DELTA = 2.0e4, 1.0e-3        # two unit blocks, prescribed interference delta


def _blocks(extra, epsN=1.0e6, top_free_x=False, shear=0.0):
    """Two stacked unit bricks, bottom z=0 clamped, top z=2 held in z. Uniaxial in z (x,y fixed,
    nu=0) unless top_free_x (then the top block's x DOFs are free and `shear` is applied in x at the
    top face -- the interface carries it by cohesion). Frictionless mortar contact across the
    MATCHING, coincident z=1 interface, gap shifted by the options in `extra`."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, E_BLK, 0.0)
    sq = [(0, 0), (1, 0), (1, 1), (0, 1)]
    for k, z in enumerate((0.0, 1.0)):
        for i, (x, y) in enumerate(sq):
            ops.node(1 + 4 * k + i, float(x), float(y), z)
    for k, z in enumerate((1.0, 2.0)):
        for i, (x, y) in enumerate(sq):
            ops.node(9 + 4 * k + i, float(x), float(y), z)
    for t in (1, 2, 3, 4):
        ops.fix(t, 1, 1, 1)
    for t in (5, 6, 7, 8):
        ops.fix(t, 1, 1, 0)
    for t in (9, 10, 11, 12):
        ops.fix(t, 0 if top_free_x else 1, 1, 0)
    for t in (13, 14, 15, 16):
        ops.fix(t, 0 if top_free_x else 1, 1, 1)
    ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.element("LadrunoBrick", 2, 9, 10, 11, 12, 13, 14, 15, 16, 1)
    ops.contactSurface(1, "-master", 4, 5, 6, 7, 8)
    ops.contactSurface(2, "-slave-segments", 4, 9, 10, 11, 12)
    ops.contact(1, 1, 2, "-mortar", "-epsN", epsN, "-outward", 0.0, 0.0, 1.0, *extra)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (13, 14, 15, 16):
        ops.load(t, shear / 4.0, 0.0, 0.0)
    _static()


def _base_pressure():
    ops.reactions()
    return sum(ops.nodeReaction(t, 3) for t in (1, 2, 3, 4))   # = sigma*A, A = 1


def test_g9_gap_offset_interference_pressure():
    """-gapOffset -delta prescribes an interference fit on straight geometry. Two blocks in series
    with the penalty spring: sigma = delta/(2L/E + 1/epsN) exactly (constant stress, flat matched
    interface); held-load augmentation removes the penalty compliance -> sigma = delta*E/(2L)."""
    eps = 1.0e6
    pen = DELTA / (2.0 / E_BLK + 1.0 / eps)
    exact = DELTA * E_BLK / 2.0
    _blocks(("-gapOffset", -DELTA, "-augment", "request"), epsN=eps)
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    assert _base_pressure() == pytest.approx(pen, rel=1e-10)
    # a second PHYSICAL step under `request` does not augment: still the penalty value.
    assert ops.analyze(1) == 0
    assert _base_pressure() == pytest.approx(pen, rel=1e-10)
    # inside the analyze_augmented bracket `request` augments -> the exact interference pressure.
    st, _, _ = analyze_augmented(ops, maxAug=40, augTol=1.0e-12)
    assert st == 0
    assert _base_pressure() == pytest.approx(exact, rel=1e-8)


def test_n2_modes_on_the_interference_blocks():
    """The three -augment modes on one closed form. `commit` (the default) augments on every
    physical commit (drifts toward delta*E/2L step by step); `never` stays on the penalty value
    even inside the analyze_augmented bracket."""
    eps = 1.0e6
    pen = DELTA / (2.0 / E_BLK + 1.0 / eps)
    _blocks(("-gapOffset", -DELTA), epsN=eps)                      # default = commit
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    p1 = _base_pressure()
    assert ops.analyze(1) == 0
    assert p1 == pytest.approx(pen, rel=1e-10)
    assert _base_pressure() > p1 * (1.0 + 1e-4)                    # augmented by the commit
    _blocks(("-gapOffset", -DELTA, "-augment", "never"), epsN=eps)
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    st, _, _ = analyze_augmented(ops, maxAug=5, augTol=1.0e-12)
    assert st == 1                                                  # never converges: inert bracket
    assert _base_pressure() == pytest.approx(pen, rel=1e-10)


def test_g9_no_offset_is_inert():
    """A coincident interface with no shift carries no load (the baseline the offset acts on)."""
    _blocks(("-augment", "never"))
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    assert abs(_base_pressure()) < 1.0e-12


def test_n2_cohesive_bond_is_step_count_independent_under_never():
    """R1 N-2: a cohesion-only mortar bond (mu = 0) failed at step 2 under the default. With
    `-augment never` and an interference prestress, a shear ramp in 1 or 5 steps converges and gives
    the same answer (elastic stick is linear in the penalty)."""
    out = []
    for n in (1, 5):
        _blocks(("-gapOffset", -DELTA, "-augment", "never", "-cohesion", 100.0, "-epsT", 1.0e6),
                epsN=1.0e6, top_free_x=True, shear=5.0)
        ops.integrator("LoadControl", 1.0 / n)
        for _ in range(n):
            assert ops.analyze(1) == 0
        out.append(ops.nodeDisp(13, 1))
    assert out[0] != 0.0
    assert abs(out[1] - out[0]) <= 1.0e-10 * abs(out[0]), out


# --------------------------------------------------------------------------- tube in a ring
R_IN, R_IF, R_OUT = 0.3, 0.5, 1.2     # inner tube [R_IN, R_IF], ring [R_IF, R_OUT]
E_TUBE, NU_TUBE = 1.0e3, 0.3
EPS_C = 1.0e5                         # contact penalty ~ 20 E/h (the R1 window was ~25 E/h;
                                      # 200 E/h chatters on this faceted interface)


def _polar(r, th, z):
    return (r * math.cos(th), r * math.sin(th), z)


def _tube_in_ring(nti=16, nto=24):
    """A tube (one radial layer, `nti` sectors) inside a ring (two radial layers, `nto` sectors),
    one element through z = [0, 0.5]; uz fixed everywhere (a plane-strain slice), the ring's outer
    face clamped. The interface r = R_IF is NON-matching (nti != nto): both are chord facets of the
    same circle, so as meshed there are spurious gaps AND penetrations. Returns the facet lists
    (master = tube outer face, slave = ring inner face) with each facet's centroid angle."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, E_TUBE, NU_TUBE)
    H = 0.5

    def tag_t(i, j, k):          # tube: radial i in {0,1}, sector j, level k
        return 1 + i + 2 * (j % nti) + 2 * nti * k

    def tag_r(i, j, k):          # ring: radial i in {0,1,2}
        return 10001 + i + 3 * (j % nto) + 3 * nto * k

    for k in (0, 1):
        for j in range(nti):
            for i, r in enumerate((R_IN, R_IF)):
                ops.node(tag_t(i, j, k), *_polar(r, 2 * math.pi * j / nti, k * H))
                # uz fixed (plane-strain slice); uy fixed on the x axis (theta = 0, pi): the load is
                # along x, so this is the symmetry plane, and it removes the tube's free rigid
                # rotation about z (frictionless contact on a near-circle has no torsional stiffness).
                on_x_axis = (2 * j) % nti == 0
                ops.fix(tag_t(i, j, k), 0, 1 if on_x_axis else 0, 1)
        for j in range(nto):
            for i, r in enumerate((R_IF, 0.5 * (R_IF + R_OUT), R_OUT)):
                t = tag_r(i, j, k)
                ops.node(t, *_polar(r, 2 * math.pi * j / nto, k * H))
                ops.fix(t, 1, 1, 1) if i == 2 else ops.fix(t, 0, 0, 1)
    e = 1
    for j in range(nti):
        ops.element("LadrunoBrick", e, tag_t(0, j, 0), tag_t(1, j, 0), tag_t(1, j + 1, 0), tag_t(0, j + 1, 0),
                    tag_t(0, j, 1), tag_t(1, j, 1), tag_t(1, j + 1, 1), tag_t(0, j + 1, 1), 1)
        e += 1
    for j in range(nto):
        for i in (0, 1):
            ops.element("LadrunoBrick", e, tag_r(i, j, 0), tag_r(i + 1, j, 0), tag_r(i + 1, j + 1, 0),
                        tag_r(i, j + 1, 0), tag_r(i, j, 1), tag_r(i + 1, j, 1), tag_r(i + 1, j + 1, 1),
                        tag_r(i, j + 1, 1), 1)
            e += 1
    master = [((j + 0.5) * 2 * math.pi / nti,
               [tag_t(1, j, 0), tag_t(1, j + 1, 0), tag_t(1, j + 1, 1), tag_t(1, j, 1)])
              for j in range(nti)]
    slave = [((j + 0.5) * 2 * math.pi / nto,
              [tag_r(0, j, 0), tag_r(0, j + 1, 0), tag_r(0, j + 1, 1), tag_r(0, j, 1)])
             for j in range(nto)]
    loaded = [tag_t(0, j, k) for j in range(nti) for k in (0, 1)]
    tube_nodes = [tag_t(i, j, k) for i in (0, 1) for j in range(nti) for k in (0, 1)]
    return master, slave, loaded, tube_nodes


def _declare(master, slave, sectors, args, outward=False):
    """sectors=1: one contact over the whole closed interface; sectors=4: the R1 quadrant split
    (facets binned by centroid angle), optionally with the quadrant bisector as -outward."""
    stag = 1
    for q in range(sectors):
        lo, hi = q * 2 * math.pi / sectors, (q + 1) * 2 * math.pi / sectors
        m = [n for th, f in master if lo <= th < hi for n in f]
        s = [n for th, f in slave if lo <= th < hi for n in f]
        ops.contactSurface(stag, "-master", 4, *m)
        ops.contactSurface(stag + 1, "-slave-segments", 4, *s)
        ow = ()
        if outward:
            b = 0.5 * (lo + hi)
            ow = ("-outward", math.cos(b), math.sin(b), 0.0)
        ops.contact(q + 1, stag, stag + 1, "-mortar", *args, *ow)
        stag += 2


def _lateral(loaded, F):
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in loaded:
        ops.load(t, F / len(loaded), 0.0, 0.0)


def _mean_ux(nodes):
    return sum(ops.nodeDisp(t, 1) for t in nodes) / len(nodes)


def _tube_tie(sectors, extra):
    master, slave, loaded, tube = _tube_in_ring()
    _declare(master, slave, sectors, ("-tie", "-epsTie", 1.0e6) + tuple(extra))
    _lateral(loaded, 1.0)
    _static(tol=1e-12, system="UmfPack")
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    return _mean_ux(tube)


def test_n1_maxgap_closed_cylinder_tie_matches_sector_split():
    """One tie around the whole cylinder with -maxGap == the 4-sector split (the R1 workaround).
    Without the guard every tube facet is also tied to the far side of the ring (the clip accepts
    the anti-parallel facet), which stiffens the tube laterally."""
    split = _tube_tie(4, ("-augment", "never"))
    guarded = _tube_tie(1, ("-augment", "never", "-maxGap", 0.1))
    bare = _tube_tie(1, ("-augment", "never"))
    assert guarded == pytest.approx(split, rel=1e-6)
    # measured: bare 4.308e-4 vs split 4.456e-4 (-3.3 %) on this slice; on the R1 pile, where the
    # antipodal ties also lock the pile's rotation, the same defect was 3.65x.
    assert bare < 0.99 * split, f"unguarded closed tie no longer over-stiff? bare {bare}, split {split}"


def test_g9_adjust_faceted_tube_starts_stress_free():
    """A faceted tube in a faceted ring (16 vs 24 chords of one circle): as meshed the interface has
    penetrations (tube nodes poke past the ring's chords) and gaps. Zero load: the unadjusted contact
    pushes the bodies apart; with -adjust nothing moves and nothing penetrates -- exactly."""
    res = {}
    for key, extra in (("raw", ()), ("adjust", ("-adjust",))):
        master, slave, loaded, tube = _tube_in_ring()
        _declare(master, slave, 4, ("-epsN", EPS_C, "-augment", "never") + extra, outward=True)
        _static(tol=1e-10, system="UmfPack", test="NormUnbalance")
        ops.integrator("LoadControl", 1.0)
        assert ops.analyze(1) == 0
        umax = max(abs(ops.nodeDisp(t, d)) for t in tube for d in (1, 2))
        res[key] = (umax, ops.ladrunoMortarPenetration())
    assert res["raw"][0] > 1.0e-6, res
    assert res["adjust"] == (0.0, 0.0), res


def test_g9_adjust_then_load_and_interference():
    """With -adjust the faceted interface is usable: a lateral load converges (the tube bears on the
    front and opens at the back), and -adjust + -gapOffset gives a uniform radial shrink fit whose
    net lateral displacement is zero by symmetry."""
    master, slave, loaded, tube = _tube_in_ring()
    _declare(master, slave, 4, ("-epsN", EPS_C, "-augment", "never", "-adjust"), outward=True)
    _lateral(loaded, 1.0)
    _static(tol=1e-10, system="UmfPack", test="NormUnbalance")
    ops.integrator("LoadControl", 0.5)
    for _ in range(2):
        assert ops.analyze(1) == 0
    assert _mean_ux(tube) > 0.0
    # shrink fit: the tube moves radially INWARD everywhere, no net lateral motion.
    master, slave, loaded, tube = _tube_in_ring()
    _declare(master, slave, 4, ("-epsN", EPS_C, "-augment", "never", "-adjust",
                                "-gapOffset", -1.0e-4), outward=True)
    _static(tol=1e-10, system="UmfPack", test="NormUnbalance")
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    assert abs(_mean_ux(tube)) < 1.0e-12
    for t in tube:
        x, y = ops.nodeCoord(t, 1), ops.nodeCoord(t, 2)
        ur = (x * ops.nodeDisp(t, 1) + y * ops.nodeDisp(t, 2)) / math.hypot(x, y)
        assert ur < 0.0


# ------------------------------------------------------------------------------ command surface
def _two_facets():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(1 + i, float(x), float(y), 0.0)
        ops.node(5 + i, float(x), float(y), 0.0)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)


@pytest.mark.parametrize("bad", [
    ("-augment", "sometimes"),
    ("-maxGap", 0.0),
    ("-maxGap", -1.0),
    ("-adjust", -1.0),
])
def test_refuses_bad_values(bad):
    _two_facets()
    with pytest.raises(Exception):
        ops.contact(1, 1, 2, "-mortar", "-epsN", 1.0e5, *bad)


@pytest.mark.parametrize("flag", [("-augment", "never"), ("-maxGap", 0.1),
                                  ("-gapOffset", -1e-3), ("-adjust",)])
def test_flags_are_mortar_only(flag):
    _two_facets()
    ops.contactSurface(3, "-slave", 5, 6, 7, 8)
    with pytest.raises(Exception):
        ops.contact(1, 1, 3, 1.0e5, 0.0, 0.0, *flag)


@pytest.mark.parametrize("flag", [("-gapOffset", -1e-3), ("-adjust",), ("-adjust", 1e-3)])
def test_gap_shift_refused_on_tie(flag):
    _two_facets()
    with pytest.raises(Exception):
        ops.contact(1, 1, 2, "-mortar", "-tie", "-epsTie", 1.0e5, *flag)


def test_accepted_forms():
    """Every documented spelling parses (tie + augment/maxGap; contact + all four)."""
    _two_facets()
    ops.contact(1, 1, 2, "-mortar", "-tie", "-epsTie", 1.0e5, "-augment", "request", "-maxGap", 0.1)
    _two_facets()
    ops.contact(1, 1, 2, "-mortar", "-epsN", 1.0e5, "-adjust", 1.0e-3, "-gapOffset", -1e-4,
                "-maxGap", 0.1, "-augment", "commit", "-outward", 0.0, 0.0, 1.0)
    _two_facets()
    ops.contact(1, 1, 2, "-mortar", "-epsN", 1.0e5, "-adjust", "-augment", "never")


def test_oracle_proto_passes():
    """The build-free oracle (numpy) behind the closed forms above."""
    import runpy
    path = os.path.join(os.path.dirname(__file__), "..", "Ladruno_implementation",
                        "contact_prototypes", "proto_adr155_r05.py")
    mod = runpy.run_path(path, run_name="proto")
    for fn in ("t1_interference", "t2_step_count", "t3_adjust_exact_zero", "t4_maxgap_window"):
        mod[fn]()
    assert mod["FAILS"] == []


@pytest.mark.parametrize("flag", [("-maxGap", 0.1), ("-gapOffset", -1e-4), ("-adjust",)])
def test_2d_pair_refuses_3d_only_options(flag):
    """-maxGap/-gapOffset/-adjust are wired to the 3D mortar lane; a 2D pair draws a named FATAL at
    handle() (analyze returns < 0) instead of silently ignoring them. -augment works in 2D."""
    def block(extra):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        ops.node(101, 0.0, 0.0)
        ops.node(102, 1.0, 0.0)
        ops.fix(101, 1, 1)
        ops.fix(102, 1, 1)
        ops.node(1, 0.0, -1e-4)
        ops.node(2, 1.0, -1e-4)
        ops.fix(1, 1, 0)
        ops.fix(2, 1, 0)
        ops.contactSurface(10, "-master", 2, 101, 102)
        ops.contactSurface(20, "-slave-segments", 2, 1, 2)
        ops.contact(1, 10, 20, "-mortar", "-epsN", 1e6, "-outward", 0.0, 1.0, *extra)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(1, 0.0, -1.0)
        ops.load(2, 0.0, -1.0)
        _static()
        ops.integrator("LoadControl", 1.0)
        try:
            return ops.analyze(1)
        except Exception:
            return -1
    assert block(("-augment", "never")) == 0
    assert block(flag) < 0
