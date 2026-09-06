"""LadrunoShellModifier (WP-91) -- G9/G10 frame-vs-shell equivalence (Zone-A).

The physical question a user actually has: if I crack a flexural member 0.35 in a
FRAME model (ETABS frame modifiers scale A and I33 directly) and I crack "the same"
member 0.35 in a SHELL model, do I land in the same place?

Model: a flexure-controlled cantilever, L/d = 10, built twice --

    (a) elasticBeamColumn, with A and I scaled directly (the frame-modifier route)
    (b) a ShellMITC4 mesh of the same member in ITS OWN PLANE, with f11/f22/f12

In-plane bending of a shell is carried by MEMBRANE action -- sigma_xx varying over
the depth -- so the modifier that scales it is f11. G10 pins the corollary that
catches people: the `m` modifiers are OUT-OF-PLANE plate bending and do exactly
NOTHING to in-plane flexure, so cracking a wall with m11 is a silent no-op.

G9 deliberately asserts RATIOS, not absolute deflections. A bilinear membrane
locks in bending, so the shell is a few percent stiff at any practical mesh and
that gap is discretisation, not the modifier. The strong, mesh-independent claim
is that the modifier scales straight THROUGH that gap: the shell/frame ratio is
identical gross and cracked. Measured across 20x2 .. 120x24 the softening ratio
is 2.8571 at every density while the absolute error moves -12.6% -> -0.5%.

See test_ladrunoShellModifier_section.py for G1-G5, G7, G8.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


@pytest.fixture(scope="module", autouse=True)
def _build_stamp():
    if not hasattr(ops, "ladrunoBuild"):
        pytest.skip("openseespy wheel fallback in use -- LadrunoShellModifier needs the fork build")
    print(f"\n[ladrunoBuild] {ops.ladrunoBuild()}")


# Flexure-controlled by construction: L/d = 10, so the Timoshenko shear term is
# 0.72% of the Euler-Bernoulli deflection.
L, D, T = 6.0, 0.6, 0.3
E, NU = 25.0e6, 0.2
P = 100.0
A = T * D
I = T * D ** 3 / 12.0
CRACK = 0.35

_NX, _NY = 40, 8


def _solve():
    ops.system("BandGeneral")
    ops.numberer("RCM")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "linear static analysis failed"


def _frame(mod_A=1.0, mod_I=1.0, nele=40):
    """ETABS frame modifiers scale A and I33 directly; there is no fork frame
    decorator, so that is modelled by scaling the section properties."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    for i in range(nele + 1):
        ops.node(i + 1, i * L / nele, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    for i in range(nele):
        ops.element("elasticBeamColumn", i + 1, i + 1, i + 2,
                    A * mod_A, E, I * mod_I, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(nele + 1, 0.0, -P, 0.0)
    _solve()
    return -ops.nodeDisp(nele + 1, 2)


def _shell(mods=None, nx=_NX, ny=_NY):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)

    def nid(i, j):
        return j * (nx + 1) + i + 1

    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.node(nid(i, j), i * L / nx, j * D / ny - D / 2.0, 0.0)

    # One fix() per node -- a second SP on an already-constrained DOF is refused.
    # Built-in end fully fixed; everything else only held planar. Drilling (6) is
    # left free: ShellMITC4 carries its own fictitious drilling stiffness.
    for j in range(ny + 1):
        for i in range(nx + 1):
            if i == 0:
                ops.fix(nid(i, j), 1, 1, 1, 1, 1, 1)
            else:
                ops.fix(nid(i, j), 0, 0, 1, 1, 1, 0)

    ops.section("ElasticMembranePlateSection", 1, E, NU, T, 0.0)
    sec = 1
    if mods:
        args = []
        for k, v in mods.items():
            args += [f"-{k}", float(v)]
        ops.section("LadrunoShellModifier", 10, 1, *args)
        sec = 10

    e = 0
    for j in range(ny):
        for i in range(nx):
            e += 1
            ops.element("ShellMITC4", e,
                        nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1),
                        sec)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(ny + 1):
        ops.load(nid(nx, j), 0.0, -P / (ny + 1), 0.0, 0.0, 0.0, 0.0)
    _solve()
    return -ops.nodeDisp(nid(nx, ny // 2), 2)     # mid-depth tip node


# ---------------------------------------------------------------------------
# G9 -- the frame and the shell soften identically
# ---------------------------------------------------------------------------
def test_g9_frame_and_shell_soften_by_the_same_exact_ratio():
    """Cracking A,I in a frame and f11/f22/f12 in a shell must produce the SAME
    softening ratio, and it must be exactly 1/CRACK in both."""
    f_gross = _frame()
    f_crack = _frame(mod_A=CRACK, mod_I=CRACK)
    s_gross = _shell()
    s_crack = _shell({"f11": CRACK, "f22": CRACK, "f12": CRACK})

    # Tolerance is direct-solve round-off, not physics. MEASURED deviations from
    # 1/CRACK: frame 3.19e-10, shell 4.20e-11 -- the FRAME is the noisier of the
    # two, since scaling A and I changes the assembled system's conditioning.
    # 1e-8 keeps ~30x headroom while still asserting 8 significant figures; a
    # regression that actually broke the modifier would miss by orders of
    # magnitude, not by round-off.
    expected = 1.0 / CRACK
    assert f_crack / f_gross == pytest.approx(expected, rel=1e-8), \
        f"frame softening {f_crack/f_gross} != {expected}"
    assert s_crack / s_gross == pytest.approx(expected, rel=1e-8), \
        f"shell softening {s_crack/s_gross} != {expected}"


def test_g9_modifier_scales_through_the_frame_shell_discretisation_gap():
    """The shell/frame deflection ratio must be IDENTICAL gross and cracked.

    This is the mesh-independent statement. The absolute shell-vs-frame gap is
    membrane shear locking (a few percent, and it converges away under
    refinement); what must hold exactly is that the modifier does not distort
    the relationship -- it scales straight through it.
    """
    f_gross = _frame()
    f_crack = _frame(mod_A=CRACK, mod_I=CRACK)
    s_gross = _shell()
    s_crack = _shell({"f11": CRACK, "f22": CRACK, "f12": CRACK})

    # MEASURED: gap 0.973336798430770 gross vs 0.973336798699989 cracked, a
    # relative deviation of 2.77e-10 -- direct-solve round-off. See the note in
    # the ratio test above for why the tolerance is 1e-8.
    assert s_crack / f_crack == pytest.approx(s_gross / f_gross, rel=1e-8), \
        ("the modifier changed the shell-vs-frame relationship: "
         f"gross {s_gross/f_gross} vs cracked {s_crack/f_crack}")


def test_g9_shell_is_within_a_few_percent_of_beam_theory():
    """Sanity floor on the model itself, so G9's ratios are not being taken from
    a mesh that is physically meaningless. Bilinear membrane locking makes the
    shell STIFFER than Euler-Bernoulli at this density; refinement removes it
    (measured -3.4% at 40x8 -> -0.5% at 120x24)."""
    eb = P * L ** 3 / (3.0 * E * I)
    s_gross = _shell()
    assert abs(s_gross - eb) / eb < 0.06, \
        f"shell {s_gross} vs Euler-Bernoulli {eb} -- mesh or model is off"


# ---------------------------------------------------------------------------
# G10 -- the trap: `m` modifiers do NOTHING to in-plane flexure
# ---------------------------------------------------------------------------
def test_g10_plate_bending_modifiers_do_not_touch_in_plane_flexure():
    """A wall/deep beam cracked with m11/m22/m12 is a SILENT NO-OP in plane.

    In-plane bending is membrane action (sigma_xx over the depth), so it is
    governed by f11. The `m` modifiers are out-of-plane plate bending. An
    engineer who reaches for "the bending ones" gets no cracking at all and no
    warning, which is why this is pinned as a test rather than left to the guide.
    """
    s_gross = _shell()
    s_m = _shell({"m11": CRACK, "m22": CRACK, "m12": CRACK})
    assert s_m == pytest.approx(s_gross, rel=1e-12), \
        f"m-modifiers changed in-plane flexure: {s_gross} -> {s_m}"


def test_g10_transverse_shear_modifiers_also_do_not_touch_in_plane_flexure():
    """Same for v13/v23 -- transverse (out-of-plane) shear is not in the in-plane
    load path either."""
    s_gross = _shell()
    s_v = _shell({"v13": CRACK, "v23": CRACK})
    assert s_v == pytest.approx(s_gross, rel=1e-12), \
        f"v-modifiers changed in-plane flexure: {s_gross} -> {s_v}"


def test_g10_all_eight_equals_membrane_only_for_in_plane_bending():
    """Corollary: for a member loaded in its own plane, cracking all eight
    modifiers is indistinguishable from cracking only the three membrane ones."""
    s_f = _shell({"f11": CRACK, "f22": CRACK, "f12": CRACK})
    s_all = _shell({"f11": CRACK, "f22": CRACK, "f12": CRACK,
                    "m11": CRACK, "m22": CRACK, "m12": CRACK,
                    "v13": CRACK, "v23": CRACK})
    assert s_all == pytest.approx(s_f, rel=1e-12), \
        f"all-eight {s_all} != membrane-only {s_f}"
