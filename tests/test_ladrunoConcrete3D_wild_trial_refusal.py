"""WP concrete3d-hang-diagnosis review, defect 1 (CRITICAL) regression gate.

Before this fix, LadrunoConcrete3D::setTrialStrain returned 0 unconditionally regardless of
whether the return map actually converged, and commitState() copied that elastic-trial fallback
into the committed state unconditionally too -- so a wild, physically nonsensical trial strain was
silently ACCEPTED and COMMITTED as if it were a valid converged answer (review probe_absorb.cpp:
one such committed step left a nominal stress of hundreds of MPa at fc=30 and an effective stress
over a GPa, poisoning every ordinary small step at that Gauss point afterward).

Two element paths, because the fix has two independent seams (ADR-86b / WP-99,
SRC/material/LadrunoMaterialStatus.h):
  * TRIAL-time: setTrialStrain returns LADRUNO_MATERIAL_REFUSED; an element that FORWARDS this code
    (LadrunoBrick, via ladrunoBrickMustCut) sees it immediately during the Newton iteration and the
    analysis fails the step right there.
  * COMMIT-time: commitState() declares the refusal via ladrunoNoteCommitRefusal() (out of band --
    Domain::commit() drops elePtr->commitState()'s own return value for EVERY element, fork or
    vanilla). An element that DISCARDS the trial-time code entirely -- stock stdBrick, which never
    even looks at setTrialStrain's return value -- still gets the step failed at commit, because
    Domain::commit() itself aborts once it sees a pending refusal.

Both gates assert ops.analyze(1) != 0 (the step is cut) within a short wall-clock bound -- the
defect this guards against was NOT a hang, it was silent acceptance, so completion is expected;
the wall-clock bound is a cheap extra check that the (unrelated, separately fixed) sub-increment
attempt budget is doing its job too.
"""
import time

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_NODES = {
    1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0), 3: (1.0, 1.0, 0.0), 4: (0.0, 1.0, 0.0),
    5: (0.0, 0.0, 1.0), 6: (1.0, 0.0, 1.0), 7: (1.0, 1.0, 1.0), 8: (0.0, 1.0, 1.0),
}
_CONN = [1, 2, 3, 4, 5, 6, 7, 8]

# UNIFORM-STRAIN kinematics (u_i = E_ij x_j on every node, prescribed via sp -- no ops.fix() needed,
# the affine field is zero at node 1 (0,0,0) by construction), so both loading directions below are
# actual material strain states, not a hand-guessed nodal pattern. Every number is the review's own
# probe_absorb.cpp, unchanged EXCEPT the shear components are DOUBLED: probe_absorb calls the kernel
# directly, whose deps[3..5] are the TENSOR shear (LadrunoConcrete3D::setTrialStrain halves the
# element's ENGINEERING shear before handing it to the kernel -- see the "engineering -> tensor
# shear" comment there), so the ENGINEERING shear this affine field must impose is 2x probe_absorb's
# raw numbers for indices 3-5 to reproduce the SAME kernel-level deps. CRACK_DIR/CRACK_STEP/STEPS is
# probe_absorb's own 30 x 2e-5 realistic cracking direction/step -- reproduced here almost exactly
# (measured: this FE model's committed kappa_p=41.1086 after the crack phase vs probe_absorb's raw-
# kernel kp=41.11, sigEff matching to 4 significant figures) so the wild jump below starts from
# essentially the SAME state probe_absorb characterized. WILD_RAW is its own direction from a 4000-
# try random search (magnitude 0.0292, tries=420, the smallest jump it found that fails the KERNEL
# directly). At the FE level (full Newton + bbar + damage, not a bare kernel call) that same
# magnitude only pushes kappa_p further (measured: 41.1 -> 64.5, rc=0, still converges) -- so
# WILD_MAG here is a magnitude SWEEP finding on the FE model itself: scale=0.1 still converges
# (kp->83.1), scale=0.5 refuses (rc!=0, kp UNCHANGED at 83.1 -- the refused trial is not committed,
# which is exactly the defect-1 fix this test guards). 1.0 keeps margin past that measured boundary.
# e6 = (exx, eyy, ezz, gxy, gyz, gxz), ENGINEERING convention (as OpenSees sp/strain always is).
_CRACK_DIR = (1.0, -0.2, -0.2, 0.6, 0.0, 0.0)
_CRACK_STEP = 2.0e-5
_CRACK_STEPS = 30
_WILD_RAW = (0.174, -0.502, 0.956, 0.964, -0.732, -1.448)
_WILD_MAG = 1.0


def _affine_disp(x, y, z, e6):
    exx, eyy, ezz, gxy, gyz, gxz = e6
    ux = exx * x + 0.5 * gxy * y + 0.5 * gxz * z
    uy = 0.5 * gxy * x + eyy * y + 0.5 * gyz * z
    uz = 0.5 * gxz * x + 0.5 * gyz * y + ezz * z
    return ux, uy, uz


def _make_concrete_material(tag=1):
    # Every parameter matches the review's probe_absorb.cpp mkp() exactly (Hp=0.01 in particular --
    # the wrapper default Hp=0.5 gives a much slower kappa_p growth, so the same 30-step crack +
    # wild jump that fails the kernel directly does NOT reproduce at the FE level under the default
    # hardening modulus).
    ops.nDMaterial("LadrunoConcrete3D", tag, 30000.0, 0.2, 30.0, 3.0, 0.1, 30.0,
                   "-Df", 0.85, "-hardening", 0.3, 0.01,
                   "-tensionLaw", "bilinear", "-epsFc", 1.0e-3,
                   "-flowPotential", "cdpm2", "-compressionDrive", "cdpm2", "-tcTemper", "proj")


def _build(ele_name, extra_args, mat_tag=1):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _NODES.items():
        ops.node(tag, x, y, z)
    _make_concrete_material(mat_tag)
    ops.fix(1, 1, 1, 1)       # pins the (already-zero-displacement) origin against rigid body drift
    ops.element(ele_name, 1, *_CONN, mat_tag, *extra_args)

    # pattern 1: the CRACK_DIR affine field at its FULL 30-step magnitude, ramped in over 30 equal
    # LoadControl steps of dlam=1/30 -- each step is exactly one of probe_absorb.cpp's 2e-5 increments.
    e6_crack_full = tuple(c * _CRACK_STEP * _CRACK_STEPS for c in _CRACK_DIR)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n, (x, y, z) in _NODES.items():
        if n == 1:
            continue
        ux, uy, uz = _affine_disp(x, y, z, e6_crack_full)
        ops.sp(n, 1, ux); ops.sp(n, 2, uy); ops.sp(n, 3, uz)

    ops.system("UmfPack")     # LadrunoConcrete3D's tangent is non-symmetric (CDPM2 non-associated flow)
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-6, 20)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / _CRACK_STEPS)
    ops.analysis("Static")


_WALL_CLOCK_BOUND_S = 30.0   # generous; a correctly-refusing step should take well under a second


def _assert_stress_unchanged(sig_before, sig_after):
    # exact equality is too strict: re-querying "stress" after a REFUSED (reverted) trial still
    # exercises the same floating-point evaluation path a second time, and the deep noise floor
    # (~1e-16 to 1e-18, many orders below any physical stress here) is not bit-identical across two
    # separate evaluations. Assert the values are unchanged to a tolerance well below anything
    # physically meaningful, not bit-for-bit.
    assert len(sig_before) == len(sig_after)
    for i, (b, a) in enumerate(zip(sig_before, sig_after)):
        tol = 1.0e-9 * max(abs(b), abs(a), 1.0)
        assert abs(a - b) <= tol, (
            f"stress[{i}] changed after a refused step: {b!r} -> {a!r} (|d|={abs(a-b):.3e} > {tol:.3e})"
        )


def _crack_then_wild_jump():
    """Phase 1: probe_absorb.cpp's own realistic 30 x 2e-5 cracking sequence (must converge
    normally). Phase 2: its own found wild increment, added as a SEPARATE Constant-scaled pattern so
    it applies at full magnitude in exactly one step regardless of pattern-1's accumulated lambda.
    Returns (rc, sig_before_wild_jump, sig_after_wild_jump) -- the "before" snapshot is taken AFTER
    the (legitimate) crack phase, since that is the committed state the wild jump must not disturb."""
    for _ in range(_CRACK_STEPS):
        rc = ops.analyze(1)
        assert rc == 0, "the realistic 2e-5 cracking steps must converge normally (not the wild trial)"
    sig_before = list(ops.eleResponse(1, "stress"))

    e6_wild = tuple(c * _WILD_MAG for c in _WILD_RAW)
    ops.timeSeries("Constant", 2)
    ops.pattern("Plain", 2, 2)
    for n, (x, y, z) in _NODES.items():
        if n == 1:
            continue
        ux, uy, uz = _affine_disp(x, y, z, e6_wild)
        ops.sp(n, 1, ux); ops.sp(n, 2, uy); ops.sp(n, 3, uz)

    ops.integrator("LoadControl", 1.0)
    t0 = time.perf_counter()
    rc = ops.analyze(1)
    wall = time.perf_counter() - t0
    assert wall < _WALL_CLOCK_BOUND_S, f"the wild step took {wall:.1f}s (bound {_WALL_CLOCK_BOUND_S}s)"
    sig_after = list(ops.eleResponse(1, "stress"))
    return rc, sig_before, sig_after


def test_ladrunoBrick_wild_trial_refuses_and_cuts_the_step():
    """LadrunoBrick forwards setTrialStrain's return code (ladrunoBrickMustCut) -- the step must
    fail at the TRIAL, not silently commit the elastic-trial fallback."""
    _build("LadrunoBrick", ["-formulation", "bbar"])
    rc, sig_before, sig_after = _crack_then_wild_jump()
    assert rc != 0, "a wild trial strain must fail the step (LadrunoBrick propagates the refusal)"
    # the committed state must be UNTOUCHED by the refused step -- the defect this guards against
    # committed the garbage trial anyway.
    _assert_stress_unchanged(sig_before, sig_after)


def test_stdBrick_wild_trial_still_refuses_via_commit_latch():
    """stdBrick does not look at setTrialStrain's return code at all (a 'discarding' element per
    ADR-86b) -- the WP-99 commit-time latch (ladrunoNoteCommitRefusal -> Domain::commit() aborts)
    must still fail the step even though the element itself never sees LADRUNO_MATERIAL_REFUSED."""
    _build("stdBrick", [])
    rc, sig_before, sig_after = _crack_then_wild_jump()
    assert rc != 0, (
        "stdBrick discards setTrialStrain's return code, but the commit-time latch "
        "(ladrunoNoteCommitRefusal) must still cut the step"
    )
    _assert_stress_unchanged(sig_before, sig_after)
