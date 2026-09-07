"""ADR-95 P4 -- the UW `DruckerPrager` two-surface return map, repaired.

WHAT THIS PINS
--------------
P0/P1 established that the tension cutoff `f2 = I1 - T(alpha2)` of the UW
Drucker-Prager was **dead code upstream**: `Jact` carries flags 0/1 per surface,
but the residual/Jacobian assembly switched on the flag VALUE

    for (i = 0; i < 2; i++)
        if      (Jact(i) == 1) { R(0) = f1-residual; g(0,0) = ...; }
        else if (Jact(i) == 2) { R(1) = f2-residual; g(1,1) = ...; }

so the `== 2` arm was unreachable.  At the corner both loop passes wrote ROW 0,
leaving `R(1) = 0` and `g(1,1)` at the dummy `1` of the det(g) = 1
initialisation, hence `gamma(1) == 0` identically: the cutoff SELECTED a
different consistent tangent but applied NO stress return, and (with
`rho_bar = 0`, where a cone return cannot move `I1` at all) the Gauss point
committed with `I1` arbitrarily far above `T`.  On the f2-only set the same loop
assembled the f1 residual into row 0, i.e. a cone return for a pure cutoff step.
The wrong `g` then poisoned the consistent tangent through `g_contra`, which is
the 5-8 decade `detAmin` outlier P1 measured at the four corner GPs that open the
quadratic wall.

P4 makes the assembly INDEX-driven (the row index is the SURFACE index) and
fixes the one genuinely wrong factor in the tangent block: the radial-return
softening term must divide by the norm of the **trial** eta, not by the norm of
the returned eta that upstream recomputes just above it for `mState(1)`.

This module tests the repaired behaviour from the outside:

  1. a corner return actually returns `I1` to `T` (hydrostatic and non-hydrostatic
     trials, and both the degenerate and the hardened corner);
  2. the f2-only branch is reachable and does a pure volumetric return;
  3. the consistent tangent at a NON-degenerate corner matches a central finite
     difference of the return map at the same trial state -- the property that a
     mis-assembled `g` destroys, and the one the global Newton actually consumes.

THE DEGENERATE CORNER IS NOT A BUG
----------------------------------
`mTo` is not a free parameter: the constructor sets `T = sqrt(2/3)*sigma_y/rho`,
which is exactly the apex of the cone.  So for a NON-hardening deck (the ADR-95
campaign deck) `f1 = f2 = 0` forces `||eta|| = sqrt(2/3)*K(alpha1) - rho*T = 0`:
every corner return lands on the apex, the returned stress is pinned at
`sigma = (T/3)*1`, and the exact consistent tangent is therefore **zero**.  The
material's pre-existing `NormCep < 1e-10` floor catches that and hands the solver
`1e-3*Ce`.  That is the correct answer, not a defect, so the FD check is run on a
HARDENING deck, where `K(alpha1) > sigma_y` opens the corner out into a genuine
two-surface vertex with a non-zero tangent.

HOW THE STATES ARE REACHED
--------------------------
Same driver as `test_adr95_dp_branch_response.py`: one `LadrunoBrick` with a
KNOWN homogeneous strain prescribed on all 24 DOF via `sp()` under the
**Penalty** handler (Transformation would hand the solver a zero-size system).
Two-stage legs remove the first load pattern, reset pseudo-time, and re-prescribe
-- so the committed state entering stage 2 is exactly the stage-1 return.
"""
import math
import os
import sys

import pytest

# ADR-95 rule: every probe runs the binary built in THIS worktree.
# LADRUNO_DIST_BIN overrides the directory -- a long campaign run holds
# dist/bin/opensees.pyd open on Windows and silently blocks the next build's copy
# step, so a rebuild is verified from a staged copy instead.
_DIST_BIN = os.environ.get("LADRUNO_DIST_BIN") or os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dist", "bin")
if os.path.isdir(_DIST_BIN) and _DIST_BIN not in sys.path:
    sys.path.insert(0, _DIST_BIN)

from _testbed import ops  # noqa: E402

pytestmark = [pytest.mark.zone_a]

# --- the ADR-95 P1 campaign deck (h20_prandtl): nu = 0.45, SY = 0.2 kPa, --------
# --- phi_txc = 20 deg, rho_bar = 0 (fully non-associated), no hardening ---------
K_EL = 1.0e4
G_EL = 4.0e3
SY = 0.2
PHI_TXC = 20.0


def _rho_from_phi_txc(phi_deg):
    s = math.sin(math.radians(phi_deg))
    return 2.0 * s / (math.sqrt(3.0) * (3.0 - s))


RHO = _rho_from_phi_txc(PHI_TXC)
RHO_BAR = 0.0
T_CUT = math.sqrt(2.0 / 3.0) * SY / RHO      # tension cutoff, I1 units (1.09904)

# Hardening deck for the non-degenerate corner.  theta = 1 puts ALL of H into
# isotropic hardening K(alpha1); H > 3G is what makes the f2-only set reachable
# (see the module docstring's algebra).
H_HARD = 6.0e4
THETA = 1.0

_HEX = [
    (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
]
PENALTY = 1.0e14


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def _model(hard=0.0, theta=0.0):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y, z) in enumerate(_HEX, start=1):
        ops.node(i, x, y, z)
    # K G sigma_y rho rho_bar Kinf Ko delta1 delta2 H theta
    ops.nDMaterial("DruckerPrager", 1, K_EL, G_EL, SY, RHO, RHO_BAR,
                   0.0, 0.0, 0.0, 0.0, hard, theta, 0.0)
    ops.element("LadrunoBrick", 1, *range(1, 9), 1,
                "-geom", "linear", "-formulation", "bbar")


def _apply(eps, ptag):
    """Prescribe u_i = eps_ij x_j on all 24 DOF and take ONE step to it."""
    exx, eyy, ezz, gxy, gyz, gzx = eps
    ops.timeSeries("Linear", ptag)
    ops.pattern("Plain", ptag, ptag)
    for i, (x, y, z) in enumerate(_HEX, start=1):
        ops.sp(i, 1, exx * x + 0.5 * gxy * y + 0.5 * gzx * z)
        ops.sp(i, 2, 0.5 * gxy * x + eyy * y + 0.5 * gyz * z)
        ops.sp(i, 3, 0.5 * gzx * x + 0.5 * gyz * y + ezz * z)
    ops.constraints("Penalty", PENALTY, PENALTY)
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-13, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static", "-noWarnings")
    assert ops.analyze(1) == 0, f"uniform-strain driver failed on stage {ptag}"


def _drive(stages, hard=0.0, theta=0.0):
    """Run a list of (exx, eyy, ezz, gxy, gyz, gzx) total-strain stages."""
    _model(hard, theta)
    for k, eps in enumerate(stages, start=1):
        if k > 1:
            # the previous pattern owns the previous sp() set: drop it, rewind
            # pseudo-time so the new Linear series ramps 0 -> 1 again, and
            # rebuild the analysis around the new constraint set.  The MATERIAL
            # state is untouched, so stage k starts from stage k-1's return.
            ops.remove("loadPattern", k - 1)
            ops.setTime(0.0)
            ops.wipeAnalysis()
        _apply(eps, k)
    return 1


def _branch(tag=1, gp=1):
    r = ops.eleResponse(tag, "material", gp, "ladrunoBranch")
    r = list(r) if r else []
    assert len(r) == 8, f"ladrunoBranch width {len(r)} != 8"
    return r


def _tangent(tag=1, gp=1):
    r = ops.eleResponse(tag, "material", gp, "ladrunoTangent")
    r = list(r) if r else []
    assert len(r) == 36, f"ladrunoTangent width {len(r)} != 36"
    return [r[6 * i:6 * i + 6] for i in range(6)]


def _stress(tag=1, gp=1):
    r = ops.eleResponse(tag, "material", gp, "stress")
    r = list(r) if r else []
    assert len(r) == 6, f"stress width {len(r)} != 6"
    return r


# ---------------------------------------------------------------------------
def test_build_stamp_is_printed():
    assert hasattr(ops, "ladrunoBuild"), "ladrunoBuild() missing -- wrong build"
    build = ops.ladrunoBuild()
    print(f"ADR-95 P4 engine build (ladrunoBuild): {build}")
    assert isinstance(build, str) and len(build) >= 7


# ---------------------------------------------------------------------------
# 1. the corner now RETURNS STRESS
# ---------------------------------------------------------------------------
def test_hydrostatic_corner_returns_I1_to_the_cutoff():
    """Pure hydrostatic tension, campaign deck: the defining P4 assertion.

    `I1_trial = 3K*tr(eps) = 9.0`, more than eight times `T = 1.09904`.  Before
    P4 the returned `I1` was the trial value; it must now be `T`.

    A purely volumetric trial has `||eta_trial|| = 0`, so `n = 0` and `gamma(0)`
    comes out exactly zero: the whole return is the cutoff's.  That also
    exercises the tangent's divide-by-zero guard -- upstream's
    `4*G*G/norm_eta*gamma(0)` is `inf*0 = NaN` here.
    """
    _drive([(1.0e-4, 1.0e-4, 1.0e-4, 0.0, 0.0, 0.0)])
    i1_trial = 3.0 * K_EL * 3.0e-4
    for gp in range(1, 9):
        b = _branch(gp=gp)
        assert b[0] == 3.0, f"GP {gp}: expected the corner (3), got branch {b[0]}"
        assert b[4] > 0.0, f"GP {gp}: trial f2 = {b[4]} not positive"
        assert b[2] > 0.0, (
            f"GP {gp}: gamma1 = {b[2]} -- the tension-cutoff multiplier is still "
            "structurally zero, i.e. the P4 fix is not in this build"
        )
        assert b[2] == pytest.approx((i1_trial - T_CUT) / (9.0 * K_EL), rel=1e-8)
        assert abs(b[1]) < 1e-12, f"GP {gp}: gamma0 = {b[1]} on a volumetric trial"
        assert b[6] == pytest.approx(T_CUT, rel=1e-8), (
            f"GP {gp}: I1 = {b[6]} != T = {T_CUT}; the cutoff did not return"
        )
        assert math.isfinite(b[7]), f"GP {gp}: detAmin = {b[7]} not finite (NaN guard)"
        # n = 0 and gamma0 = 0 leave Cep = Ce - K*1(x)1 = 2G*IIdev, whose acoustic
        # tensor is G*I + (G/3)*n(x)n  =>  det/(2G)^3 = (4/3)/8 = 1/6, exactly.
        assert b[7] == pytest.approx(1.0 / 6.0, rel=1e-9), (
            f"GP {gp}: detAmin = {b[7]} != 1/6 -- the volumetric-return tangent "
            "should be exactly the deviatoric elastic operator"
        )
        assert 1e-3 <= abs(b[7]) <= 1e2, f"GP {gp}: |detAmin| = {abs(b[7])} not O(1)"

        C = _tangent(gp=gp)
        cmax = max(abs(v) for row in C for v in row)
        assert all(math.isfinite(v) for row in C for v in row), "tangent not finite"
        assert cmax / (2.0 * G_EL) < 100.0, f"GP {gp}: max|C|/(2G) = {cmax/(2*G_EL)}"


def test_deviatoric_corner_lands_on_the_apex_with_a_zero_tangent():
    """The P0 corner leg (3,1,1)e-4 -- now with a real return.

    Both multipliers are positive here.  Because `T` IS the cone apex on a
    non-hardening deck, the return lands exactly on `sigma = (T/3)*1` and the
    exact consistent tangent is ZERO -- caught by the material's own
    `NormCep < 1e-10` floor, which substitutes `1e-3*Ce`.  So `detAmin` is
    `1e-9 * 0.479...` here BY CONSTRUCTION and is asserted as such rather than
    as O(1): a non-degenerate corner needs hardening (next test).
    """
    _drive([(3.0e-4, 1.0e-4, 1.0e-4, 0.0, 0.0, 0.0)])
    i1_trial = 3.0 * K_EL * 5.0e-4
    for gp in range(1, 9):
        b = _branch(gp=gp)
        assert b[0] == 3.0, f"GP {gp}: expected the corner (3), got {b[0]}"
        assert b[1] > 0.0, f"GP {gp}: gamma0 = {b[1]} not positive"
        assert b[2] > 0.0, f"GP {gp}: gamma1 = {b[2]} not positive"
        assert b[2] == pytest.approx((i1_trial - T_CUT) / (9.0 * K_EL), rel=1e-8)
        assert b[6] == pytest.approx(T_CUT, rel=1e-8), (
            f"GP {gp}: I1 = {b[6]} != T = {T_CUT}"
        )
        C = _tangent(gp=gp)
        assert all(math.isfinite(v) for row in C for v in row), "tangent not finite"
        cmax = max(abs(v) for row in C for v in row)
        assert cmax / (2.0 * G_EL) < 100.0, f"GP {gp}: max|C|/(2G) = {cmax/(2*G_EL)}"
        # the floor, and NOT the 1e2..1e3 x 2G operator P1 measured pre-fix
        expect = 1e-9 * (K_EL + 4.0 * G_EL / 3.0) / (8.0 * G_EL)
        assert b[7] == pytest.approx(expect, rel=1e-6), (
            f"GP {gp}: detAmin = {b[7]}; the apex tangent should be the "
            f"1e-3*Ce floor ({expect}), not a pathological operator"
        )
        sig = _stress(gp=gp)
        for a in range(3):
            assert sig[a] == pytest.approx(T_CUT / 3.0, rel=1e-8)
        for a in range(3, 6):
            assert abs(sig[a]) < 1e-9


def test_f2_only_branch_is_reachable_and_returns_volumetrically():
    """Branch 2 -- unreachable on the campaign deck, reachable with hardening.

    On a NON-hardening deck `f1 <= 0 < f2` is algebraically impossible (it would
    need `||eta|| < 0`), which is why P1 never saw branch 2.  With `theta = 1`
    and `H > 3G`, `K(alpha1)` outruns the elastic unloading of `eta`, so a
    hardened point can sit inside the cone yet above the cutoff.  Pre-fix this
    branch assembled the **f1** residual into row 0, i.e. it did a cone return.
    """
    stage1 = (0.0, 0.0, 0.0, 8.0e-4, 0.0, 0.0)      # shear to yield + harden
    stage2 = (5.0e-5, 5.0e-5, 5.0e-5, 0.0, 0.0, 0.0)  # shear off, tension on
    _drive([stage1, stage2], hard=H_HARD, theta=THETA)
    i1_trial = 3.0 * K_EL * 1.5e-4
    for gp in range(1, 9):
        b = _branch(gp=gp)
        assert b[3] < 0.0, f"GP {gp}: trial f1 = {b[3]} -- expected inside the cone"
        assert b[4] > 0.0, f"GP {gp}: trial f2 = {b[4]} -- expected above the cutoff"
        assert b[0] == 2.0, f"GP {gp}: expected the f2-only branch (2), got {b[0]}"
        assert abs(b[1]) < 1e-14, f"GP {gp}: gamma0 = {b[1]} on a cutoff-only step"
        assert b[2] > 0.0, f"GP {gp}: gamma1 = {b[2]} not positive"
        assert b[2] == pytest.approx((i1_trial - T_CUT) / (9.0 * K_EL), rel=1e-8)
        assert b[6] == pytest.approx(T_CUT, rel=1e-8), f"GP {gp}: I1 = {b[6]}"
        # a pure volumetric return leaves the deviatoric elastic operator
        assert b[7] == pytest.approx(1.0 / 6.0, rel=1e-9), f"GP {gp}: detAmin {b[7]}"
        assert 1e-3 <= abs(b[7]) <= 1e2


# ---------------------------------------------------------------------------
# 2. the CONSISTENT TANGENT at a non-degenerate corner
# ---------------------------------------------------------------------------
_FD_STAGE1 = (0.0, 0.0, 0.0, 8.0e-4, 0.0, 0.0)
_FD_STAGE2 = (1.0e-4, 1.0e-4, 1.0e-4, 8.0e-4, 0.0, 0.0)
_FD_H = 1.0e-7


def _corner_state(perturb=None):
    """Run the hardened corner path, optionally with stage 2 perturbed."""
    s2 = list(_FD_STAGE2)
    if perturb is not None:
        k, dh = perturb
        s2[k] += dh
    _drive([_FD_STAGE1, tuple(s2)], hard=H_HARD, theta=THETA)
    return _branch(gp=1), _stress(gp=1), _tangent(gp=1)


def test_hardened_corner_is_non_degenerate():
    """Guard for the FD test: both surfaces active, both multipliers positive,
    and ||eta|| strictly positive at the return (i.e. NOT the apex)."""
    b, sig, C = _corner_state()
    assert b[0] == 3.0, f"expected the corner (3), got branch {b[0]}"
    assert b[1] > 1e-9, f"gamma0 = {b[1]} -- cone leg inactive"
    assert b[2] > 1e-9, f"gamma1 = {b[2]} -- cutoff leg inactive"
    assert b[6] == pytest.approx(T_CUT, rel=1e-8), f"I1 = {b[6]} != T"
    dev = [sig[i] - T_CUT / 3.0 for i in range(3)] + list(sig[3:])
    norm_eta = math.sqrt(sum(dev[i] ** 2 for i in range(3))
                         + 2.0 * sum(dev[i] ** 2 for i in range(3, 6)))
    assert norm_eta > 1e-2, (
        f"||eta|| = {norm_eta} at the return -- this is the degenerate apex, "
        "so the FD check below would be comparing against a zero tangent"
    )
    assert 1e-3 <= abs(b[7]) <= 1e2, f"|detAmin| = {abs(b[7])} not O(1)"
    cmax = max(abs(v) for row in C for v in row)
    assert cmax / (2.0 * G_EL) < 100.0, f"max|C|/(2G) = {cmax / (2 * G_EL)}"


@pytest.mark.parametrize("k", [0, 1, 2, 3, 4, 5])
def test_corner_tangent_matches_finite_difference(k):
    """Central FD of the return map, one Voigt direction per case.

    The consistent tangent of a return map is by definition the derivative of the
    returned stress with respect to the strain AT THE SAME committed state, so a
    central difference of the whole two-stage path (which re-runs stage 1
    identically, hence commits an identical state) must reproduce column k of C.
    A mis-assembled Jacobian -- the P4 defect -- breaks exactly this, because
    `g_contra` enters every rank-one term of the tangent.

    Voigt column 3..5 differentiates w.r.t. ENGINEERING shear, which is what the
    driver's `gxy/gyz/gzx` arguments and `mCep`'s columns both use.
    """
    _, sig0, C = _corner_state()
    _, sig_p, _ = _corner_state(perturb=(k, +_FD_H))
    _, sig_m, _ = _corner_state(perturb=(k, -_FD_H))

    fd = [(sig_p[i] - sig_m[i]) / (2.0 * _FD_H) for i in range(6)]
    an = [C[i][k] for i in range(6)]
    scale = max(max(abs(v) for v in an), max(abs(v) for v in fd), 1e-30)
    for i in range(6):
        assert abs(fd[i] - an[i]) <= 1e-3 * scale, (
            f"C[{i}][{k}] = {an[i]} but the finite difference of the return map "
            f"gives {fd[i]} (column scale {scale}); the consistent tangent does "
            "not correspond to the return actually performed"
        )
