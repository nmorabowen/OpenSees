"""ADR-85 T3 -- `-thickness` scaling (gate G-T3(f)).

WHAT THIS FILE IS
    Written against `_adr85_t3_design.md` oracle finding A6 ("the correct
    -thickness 2 twin's committed multiplier is h*lambda_1 ... the bitwise 2x
    traction claim holds; double-scaling is the gated h^2 error") and the ADR
    How/7 thickness convention (h multiplies epsN/epsT/muc/the friction
    clamps/the tie stiffness ONCE, at the 2D injection site; lambda_N is NEVER
    separately scaled -- it inherits h through epsN*gbar), for the ADR-85 T3
    C++. Iterated against the live built engine (dist/bin) after the first
    gate pass -- see the two GATE RERUN notes below.

DEFECT 2 (gate rerun finding) -- REDESIGNED, not patched. This file
originally imposed the penetration via `ops.sp`, which `constraints(
"LadrunoContact")` does not enforce for non-homogeneous values (the engine
warns and silently leaves the DOF free -- confirmed by the gate run). The
fix is the gate's own suggested option (c): build the slave nodes ALREADY at
their penetrated positions (an interference fit baked directly into the
reference geometry) and declare them fully `fix`ed (homogeneous, so always
enforced regardless of the constraint handler). With every DOF on both
surfaces fixed, the contact force depends on GEOMETRY x STIFFNESS alone -- no
load pattern, no equilibrium search, nothing for `-thickness` to interact
with except the declared penalty/tie densities themselves.

GATE RERUN #2 (found DURING the live-engine rerun, test-side) -- reactions
read back as EXACTLY 0.0 on a fully-`fix`ed, contact-only node. Probed
directly: `ops.reactions()` / `ops.nodeReaction()` does NOT flow the mortar
contact's own force at a node whose ONLY connection is the contact itself
(matching test_adr41_mortar_c2_1.py's documented finding, "the contact FE's
force variants are size-0, so contact forces don't flow into
ops.reactions()" -- confirmed here to apply on the MASTER side too, not just
free slave DOFs). The interference-fit deck therefore serves the GEOMETRY
claim (below) perfectly (confirmed byte-identical, exactly as designed), but
cannot supply a nonzero, readable FORCE observable by itself.

The fix, ALSO probed and confirmed: insert a genuinely ORDINARY element
(a very stiff `zeroLength` spring) between each slave node and a truly
`fix`ed ground node, and read the reaction at the GROUND node instead --
an ordinary element's force DOES flow into `ops.reactions()` correctly (only
the bare Contact FE's own contribution is the gap). With the spring several
orders stiffer than the contact's own effective penalty, the slave stays
"almost exactly" at its interference-fit position (probed ratio to the ideal
2x at K=1e16: 1.9999999998999993, i.e. ~5e-11 relative) while now exposing a
correctly-flowing reaction. This is why the geometry-byte-identity claim and
the reaction-doubling claim below live on TWO DIFFERENT deck constructions
rather than one: the fully-`fix`ed deck is exact for geometry but mute for
force; the spring-anchored deck is (very nearly, not exactly) exact for
force and reads a force. Each is used for the claim it can actually support.

WHY `system UmfPack` FOR THE GEOMETRY DECK, NOT `FullGeneral`
    Every DOF on both surfaces of `_thickness_deck_exact` is `fix`ed, so the
    model has ZERO free equations. tests/test_soe_zero_free_equations.py is
    this repo's own record that the classic dense SOEs (FullGeneral among
    them) NULL-POINTER-CRASH the whole process (`exit(-1)`, no Python
    traceback) on the FIRST zero-equation `setSize()` -- UmfPack is that
    file's own "always worked" control. The spring-anchored deck
    (`_thickness_deck_readable`) leaves the slave DOFs nominally free (only
    pinned via the stiff spring), so it has real free equations and uses
    `FullGeneral` like the rest of this battery.

WHY BITWISE 2x IS THE RIGHT BAR FOR GEOMETRY, AND A TIGHT-BUT-NOT-BITWISE BAR
FOR THE REACTION
    Multiplying an IEEE-754 double by exactly 2.0 only adjusts the exponent
    field -- no mantissa rounding occurs, and the interference-fit deck's
    KINEMATICS are identical between the h=1/h=2 twins by fiat (same fixed
    geometry, no solve-dependence at all), so `float.hex` equality on the
    penetration query is the right bar and is what test_2d_thickness_
    geometry_byte_identical uses. The REACTION comes through the finite-K
    anchor instead, which is (by construction, see above) an APPROXIMATION
    of an exact kinematic pin, not an exact one -- so its own doubling claim
    is gated at a tight RELATIVE tolerance (1e-6, four-plus orders inside the
    ~5e-11 relative deviation actually measured) rather than literal
    `float.hex` equality, with the reasoning stated here rather than silently
    loosened.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

PEN = 1.0e-3     # interference-fit penetration depth, IDENTICAL between twins
EPSN = 1.0e6
K_ANCHOR = 1.0e16   # see the module docstring's sizing/precision argument


def _thickness_deck_exact(h, epsN=EPSN, pen=PEN):
    """The fully-`fix`ed interference-fit pair: EXACT geometry, but reactions
    on the bare contact-only nodes read back 0.0 (see the module docstring)
    -- used ONLY for the geometry (ladrunoMortarPenetration) claim.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(101, -1.0, 0.0)
    ops.node(102, 1.0, 0.0)
    ops.fix(101, 1, 1)
    ops.fix(102, 1, 1)
    ops.contactSurface(10, "-master", 2, 101, 102)

    ops.node(1, -0.5, -pen)
    ops.node(2, 0.5, -pen)
    ops.fix(1, 1, 1)
    ops.fix(2, 1, 1)
    ops.contactSurface(20, "-slave-segments", 2, 1, 2)

    args = [1, 10, 20, "-mortar", "-epsN", epsN, "-outward", 0.0, 1.0]
    if h != 1.0:
        args += ["-thickness", h]
    ops.contact(*args)

    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("UmfPack")            # zero free DOFs -- see the module docstring
    ops.test("NormDispIncr", 1.0e-13, 40, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, f"h={h}: the interference-fit step did not converge"

    return ops.ladrunoMortarPenetration()


def _thickness_deck_readable(h, epsN=EPSN, pen=PEN, k=K_ANCHOR, force_flag=False):
    """The spring-anchored twin: slave DOFs are nominally FREE but pinned to
    the SAME interference-fit position via a very stiff `zeroLength` spring
    to a genuinely `fix`ed ground node -- an ORDINARY element, whose force
    DOES flow into `ops.reactions()` (see the module docstring). Returns the
    ground-node reactions (the readable force observable) and the
    penetration (a secondary sanity check -- not the byte-identity claim,
    which lives on `_thickness_deck_exact`)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(101, -1.0, 0.0)
    ops.node(102, 1.0, 0.0)
    ops.fix(101, 1, 1)
    ops.fix(102, 1, 1)
    ops.contactSurface(10, "-master", 2, 101, 102)

    ops.node(1, -0.5, -pen)
    ops.node(2, 0.5, -pen)
    ops.fix(1, 1, 0)                 # x fixed, y free (pinned by the spring)
    ops.fix(2, 1, 0)
    for tag, x in ((1, -0.5), (2, 0.5)):
        gnd = 950 + tag
        ops.node(gnd, x, -pen)
        ops.fix(gnd, 1, 1)
        ops.uniaxialMaterial("Elastic", gnd, k)
        ops.element("zeroLength", gnd, gnd, tag, "-mat", gnd, "-dir", 2)
    ops.contactSurface(20, "-slave-segments", 2, 1, 2)

    args = [1, 10, 20, "-mortar", "-epsN", epsN, "-outward", 0.0, 1.0]
    if h != 1.0 or force_flag:
        args += ["-thickness", h]
    ops.contact(*args)

    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-13, 60, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, f"h={h}: the spring-anchored step did not converge"

    ops.reactions()
    r1 = ops.nodeReaction(951, 2)
    r2 = ops.nodeReaction(952, 2)
    gap = ops.ladrunoMortarPenetration()
    return r1, r2, gap


def test_2d_thickness_doubles_reactions():
    """G-T3(f): `-thickness 2` scales the contact-induced reaction (read at
    the spring anchor's ground node, see the module docstring) by
    approximately EXACTLY 2 -- gated at 1e-6 relative, not literal bitwise
    equality, because the observable itself comes through a finite (if very
    stiff) approximation of an exact pin; see the module docstring for the
    measured ~5e-11 relative floor that bound leaves four-plus orders of
    headroom over.
    """
    r1_a, r2_a, _gap_a = _thickness_deck_readable(1.0)
    r1_b, r2_b, _gap_b = _thickness_deck_readable(2.0)

    assert r1_a != 0.0 and r2_a != 0.0, (
        "the default-thickness twin produced a zero reaction -- the deck is "
        "not exercising the contact force at all")

    for name, ra, rb in (("node 1", r1_a, r1_b), ("node 2", r2_a, r2_b)):
        assert abs(rb - 2.0 * ra) / abs(2.0 * ra) < 1.0e-6, (
            f"{name}: -thickness 2 did not scale the reaction ~2x: "
            f"h=1 reaction {ra!r}, h=2 reaction {rb!r} (ratio "
            f"{rb / ra!r}, expected ~2.0)")


def test_2d_thickness_geometry_byte_identical():
    """G-T3(f): the GEOMETRY observable (ladrunoMortarPenetration -- a length,
    dimension-blind per ADR-85 How/7) is byte-identical between the h=1 and
    h=2 twins, on the EXACT (fully-fixed) interference-fit deck: the
    kinematics are identical by construction, so -thickness must scale FORCE
    densities only, never the kinematics.
    """
    gap_a = _thickness_deck_exact(1.0)
    gap_b = _thickness_deck_exact(2.0)
    assert float.hex(gap_a) == float.hex(gap_b), (
        f"-thickness perturbed the geometry: gap(h=1)={gap_a!r} "
        f"({float.hex(gap_a)}) vs gap(h=2)={gap_b!r} ({float.hex(gap_b)}) -- "
        "the thickness multiplier must never reach lambda_N / the kinematics")


def test_2d_thickness_default_is_one_no_op():
    """A CONTROL on the two tests above: declaring `-thickness 1.0` EXPLICITLY
    must reproduce the same reactions/geometry as omitting the flag entirely
    (the declared default, ADR-85 design note: "default 1.0, loudly
    documented") -- otherwise "doubles" could pass on a build where the
    flag's mere PRESENCE (not its value) changes something. Geometry is
    checked on the EXACT deck (bitwise); the reaction is checked on the
    READABLE deck (tight relative, same reasoning as the doubling test).
    """
    gap_default = _thickness_deck_exact(1.0)
    r1_default, r2_default, _g = _thickness_deck_readable(1.0)

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(101, -1.0, 0.0)
    ops.node(102, 1.0, 0.0)
    ops.fix(101, 1, 1)
    ops.fix(102, 1, 1)
    ops.contactSurface(10, "-master", 2, 101, 102)
    ops.node(1, -0.5, -PEN)
    ops.node(2, 0.5, -PEN)
    ops.fix(1, 1, 1)
    ops.fix(2, 1, 1)
    ops.contactSurface(20, "-slave-segments", 2, 1, 2)
    ops.contact(1, 10, 20, "-mortar", "-epsN", EPSN, "-thickness", 1.0,
                "-outward", 0.0, 1.0)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-13, 40, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    gap_explicit = ops.ladrunoMortarPenetration()
    assert float.hex(gap_default) == float.hex(gap_explicit), (
        f"explicit -thickness 1.0 differs from the omitted-flag default "
        f"(geometry): {gap_default!r} vs {gap_explicit!r}")

    r1_explicit, r2_explicit, _g2 = _thickness_deck_readable(1.0, force_flag=True)
    assert r1_default == r1_explicit and r2_default == r2_explicit, (
        f"explicit -thickness 1.0 differs from the omitted-flag default "
        f"(reaction): ({r1_default!r},{r2_default!r}) vs "
        f"({r1_explicit!r},{r2_explicit!r})")
