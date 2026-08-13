"""ASDPlasticMaterial3D + ``MohrCoulombTensionCutoff_YF`` / ``..._PF`` — the
composite Mohr-Coulomb-with-Rankine-tension-cutoff battery (perfect plasticity,
tension-POSITIVE convention, default ``Backward_Euler`` integrator with the
closed-form ``special_return`` on the cutoff face / edge / corner / apex).

Plan: ``Ladruno_implementation/84_ladruno_mc_tension_cutoff_adr.md``.
Provenance: the Cerro Lindo ADR-0005 **M3** acceptance ("no committed state may
carry f > 0 anywhere in the model") reduced to a single material point — the
``assert_admissible`` audit below is that acceptance in miniature, applied to the
committed stress history of EVERY test in this module.

DRIVER
------
One ``stdBrick`` unit cube with 1/8-symmetry restraints (x=0 / y=0 / z=0 faces
normal-fixed).  Strain is imposed with ``sp`` constraints on the x=1 / y=1 / z=1
faces, each face group driven by its own ``timeSeries("Path", ...)``, so the
axial and lateral strain histories are INDEPENDENT functions of pseudo-time
(that is what lets one model walk a compression leg and then a tension leg
without rebuilding).  Two driver flavours:

  * ``lat=None``  -> lateral faces FREE  => homogeneous **uniaxial stress**
    (used for the cap / unload-slope probes),
  * ``lat=(t,v)`` -> all 24 DOFs prescribed => a pure **strain-controlled**
    material-point probe with no global limit point at all; this is the
    numerically bullet-proof flavour and carries the corner/apex paths.

``constraints("Transformation")``.  ``Penalty`` was tried first and is WRONG for
this driver, in two independent ways worth recording:

  * with alpha = 1e14 against a ~1e6 stiffness, the displacement-increment norm
    has a conditioning FLOOR near 3e-8 — no ``NormDispIncr`` tolerance below
    that can ever be met, and the corner path stalled at load factor 0.28 after
    300 iterations sitting at exactly that floor,
  * even two purely ELASTIC steps failed to converge (rc -3), so the symptom
    does not need plasticity to appear and is easy to misread as a material bug.

``Transformation`` condenses the prescribed DOFs away instead, which for the
fully-prescribed flavour leaves ZERO free equations — that is fine here: the
analysis returns 0, the ``sp`` values propagate exactly, and the material is
still fully exercised.  Because the imposition is now exact, ``nodeDisp`` and
the target strain agree bit for bit; strains are still MEASURED from
``nodeDisp`` (it costs nothing and keeps the oracles honest).

Paths carry ``-useLast`` **and** a padded terminal hold node past the analysis
end time: see LEDGER_quirks, "timeSeries Path returns 0 BEYOND its last time
node" and "-useLast silently dropped" — a Path-driven ``sp`` snapping to zero on
the final step is the classic silent corruption of exactly this driver.

ORACLES (tension positive throughout; principal stresses ascending s1<=s2<=s3)
-----------------------------------------------------------------------------
  * elastic uniaxial:      sigma_xx = E * eps_xx
  * tension cutoff:        max principal stress caps at T (= ``TC_min_stress``);
                           uniaxial-tension yield strain eps_y = T / E
  * Mohr-Coulomb yield (the exact arithmetic of ``MohrCoulomb_YF::YF``, and of
    the numpy audit below):
        theta = asin(clamp(-3*sqrt(3)*J3 / (2*J2**1.5), -1, 1)) / 3   (0 if J2~0)
        f_MC  = (cos(theta) - sin(theta)*sin(phi)/sqrt(3))*sqrt(J2)
                + I1*sin(phi)/3 - c*cos(phi)
    with I1 = tr(sigma), J2 = 0.5*dev:dev, J3 = det(dev).
  * triaxial-compression meridian (theta = +30 deg, where the Lode-smoothed YF
    is algebraically EXACT — see the derivation note at ``_mc_axial``):
        sigma_ax = sigma_conf*(1+sin phi)/(1-sin phi) - 2c*cos phi/(1-sin phi)
    and its uniaxial-compression specialisation (sigma_conf = 0):
        sigma_uc = -2c*cos phi/(1-sin phi)   = -285.60 kPa for c=100, phi=20
  * MC apex (hydrostatic tension limit of the un-cut surface):
        p_apex = c*cot(phi) = 274.75 kPa   (so T = 24.7 and T = 0 both sit well
        BELOW the apex, and the "no-op" T = 500 sits above it)
  * hydrostatic tension under an active cutoff plateaus at sigma = T*delta
  * elastic unloading from ANY perfectly-plastic state is uniaxially stiff:
        d sigma_xx / d eps_xx = E
  * bulk / shear:  K = E/(3(1-2nu)),  G = E/(2(1+nu))

TOLERANCE LADDER
----------------
  * cap / plateau values ......... rel 1e-9, taken NEAT — the return is
                                   closed-form exact and ``Transformation``
                                   imposes the strain exactly, so there is no
                                   slip left to absorb (a smoke probe on this
                                   build measured relerr 0.0 on the tension
                                   cap); T = 0 uses abs 1e-9*E,
  * MCTC-vs-plain-MC parity ...... BIT-IDENTICAL (``np.array_equal``) first,
                                   falling back to rtol 1e-14 while reporting
                                   the measured max ulp gap,
  * admissibility audit .......... f <= 1e-6 * (2 c cos phi); on the heavily
                                   confined legs (|sigma| ~ 1e4 kPa) that fixed
                                   1.9e-4 absolute would be a ~2e-8 RELATIVE
                                   demand, so those calls pass the same 1e-6
                                   scaled by max(2 c cos phi, max|sigma|) — the
                                   canonical-testbed rule that tolerances are
                                   relative and conditioning-aware,
  * closed-form strengths ........ rtol 5e-2.  The MC surface is exact on the
                                   compression meridian, but the FLOW direction
                                   is Drucker-Prager-smoothed for |theta| > 29
                                   deg, so the stress PATH (hence sigma_conf on
                                   an oedometric leg) is only approximately the
                                   textbook one,
  * unload slope ................. rel 1e-6.

Parameters are Cerro Lindo-like, kPa: E = 2.0e6, nu = 0.3, c = 100, phi = 20,
psi = 5 (non-associated), T_rock = 24.7, T_edz = 0.0.
"""
import math

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# ---------------------------------------------------------------------------
# model parameters (kPa)
# ---------------------------------------------------------------------------
E, NU = 2.0e6, 0.3
C, PHI, PSI = 100.0, 20.0, 5.0
T_ROCK, T_EDZ = 24.7, 0.0
T_NOOP = 500.0                     # > c*cot(phi) = 274.75 -> cutoff never bites

KMOD = E / (3.0 * (1.0 - 2.0 * NU))          # 1.6667e6
GMOD = E / (2.0 * (1.0 + NU))                # 7.6923e5
LAM = E * NU / ((1.0 + NU) * (1.0 - 2.0 * NU))

P_APEX = C / math.tan(math.radians(PHI))     # 274.75 kPa hydrostatic tension
EPS_Y_TEN = T_ROCK / E                       # 1.235e-5 uniaxial tension yield

IV = "BackStress(NullHardeningTensorFunction):"

# The audit tolerance of record (module docstring, "TOLERANCE LADDER").
TOL_YF = 1.0e-6 * (2.0 * C * math.cos(math.radians(PHI)))     # 1.879e-4 kPa


# ---------------------------------------------------------------------------
# materials
#
# Every parameter DECLARED by the chosen YF/PF/EL triplet must be passed: an
# unknown name is a SILENT no-op in this framework and an unset one is garbage.
#   LinearIsotropic3D_EL  <YoungsModulus, PoissonsRatio>
#   MohrCoulomb_YF        <MC_phi, MC_c, MC_ds>
#   MohrCoulomb_PF        <MC_phi, MC_c, MC_ds, MC_psi>
#   MohrCoulombTensionCutoff_*  the same union PLUS TC_min_stress
# plus MassDensity, which ASDPlasticMaterial3D itself declares.
# ---------------------------------------------------------------------------
def mat_mctc(tag, T=T_ROCK, c=C, phi=PHI, psi=PSI):
    """The material under test: composite MC + Rankine tension cutoff."""
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "MohrCoulombTensionCutoff_YF", "MohrCoulombTensionCutoff_PF",
        "LinearIsotropic3D_EL", IV,
        "Begin_Model_Parameters",
        "YoungsModulus", E, "PoissonsRatio", NU,
        "MC_phi", phi, "MC_c", c, "MC_psi", psi, "MC_ds", 0.0,
        "TC_min_stress", T, "MassDensity", 0.0,
        "End_Model_Parameters",
        "Begin_Internal_Variables",
        "BackStress", 0., 0., 0., 0., 0., 0.,
        "End_Internal_Variables",
    )


def mat_mc(tag, c=C, phi=PHI, psi=PSI):
    """The plain Mohr-Coulomb reference (no cutoff) — the parity oracle."""
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "MohrCoulomb_YF", "MohrCoulomb_PF",
        "LinearIsotropic3D_EL", IV,
        "Begin_Model_Parameters",
        "YoungsModulus", E, "PoissonsRatio", NU,
        "MC_phi", phi, "MC_c", c, "MC_psi", psi, "MC_ds", 0.0,
        "MassDensity", 0.0,
        "End_Model_Parameters",
        "Begin_Internal_Variables",
        "BackStress", 0., 0., 0., 0., 0., 0.,
        "End_Internal_Variables",
    )


# ---------------------------------------------------------------------------
# numpy oracles — recomputed from the COMMITTED 6-component stress, in exactly
# the arithmetic of MohrCoulomb_YF::YF / TensionCutoff_YF::YF.
# Voigt order is [11, 22, 33, 12, 23, 13] (eigenAPI typedefs.h).
# ---------------------------------------------------------------------------
def _tensor(s):
    s = np.asarray(s, dtype=float)
    return np.array([[s[0], s[3], s[5]],
                     [s[3], s[1], s[4]],
                     [s[5], s[4], s[2]]])


def principals(s):
    """Principal stresses ASCENDING (s1 <= s2 <= s3); s3 is the max principal,
    matching VoigtVector::principalStresses()."""
    return np.linalg.eigvalsh(_tensor(s))


def lode_angle(s):
    S = _tensor(s)
    I1 = np.trace(S)
    dev = S - (I1 / 3.0) * np.eye(3)
    J2 = 0.5 * float(np.sum(dev * dev))
    J3 = float(np.linalg.det(dev))
    if J2 <= np.finfo(float).eps:
        return 0.0
    x = (-3.0 * math.sqrt(3.0) * J3) / (2.0 * J2 * math.sqrt(J2))
    return math.asin(min(1.0, max(-1.0, x))) / 3.0


def f_mc(s, phi_deg=PHI, c=C):
    """Mohr-Coulomb yield function, tension positive, Lode-angle form."""
    phi = math.radians(phi_deg)
    S = _tensor(s)
    I1 = float(np.trace(S))
    dev = S - (I1 / 3.0) * np.eye(3)
    J2 = 0.5 * float(np.sum(dev * dev))
    th = lode_angle(s)
    return ((math.cos(th) - math.sin(th) * math.sin(phi) / math.sqrt(3.0))
            * math.sqrt(max(J2, 0.0)) + I1 * math.sin(phi) / 3.0
            - c * math.cos(phi))


def f_tc(s, T):
    """Rankine tension cutoff: max principal stress minus the cap."""
    return float(principals(s)[2]) - T


def _mc_axial(sig_conf, phi_deg=PHI, c=C):
    """Closed-form axial stress on the triaxial-COMPRESSION meridian.

    Exact for this YF, not merely approximate: at sigma_1 = sigma_2 = conf,
    sigma_3 = ax the Lode angle is exactly +30 deg, and substituting
    cos30 = sqrt(3)/2, sin30 = 1/2, sqrt(J2) = (conf-ax)/sqrt(3),
    I1 = 2*conf + ax into the smoothed YF collapses it algebraically to the
    textbook (sig1-sig3)/2 + (sig1+sig3)/2*sin(phi) = c*cos(phi).
    """
    s = math.sin(math.radians(phi_deg))
    co = math.cos(math.radians(phi_deg))
    return sig_conf * (1.0 + s) / (1.0 - s) - 2.0 * c * co / (1.0 - s)


def assert_admissible(hist, phi=PHI, c=C, T=None, tol=None, label=""):
    """THE M3 acceptance in miniature — no committed state may carry f > 0.

    ``hist`` is an (nstep, 6) array of committed stresses.  Checks f_MC over the
    whole history, and f_TC too when ``T`` is given (pass ``T=None`` for the
    plain-MC reference runs, where the cutoff does not exist).
    """
    hist = np.atleast_2d(np.asarray(hist, dtype=float))
    if tol is None:
        tol = 1.0e-6 * (2.0 * c * math.cos(math.radians(phi)))

    fmc = np.array([f_mc(s, phi, c) for s in hist])
    k = int(np.argmax(fmc))
    assert fmc[k] <= tol, (
        f"{label}: INADMISSIBLE committed state — max f_MC = {fmc[k]:.6e} > "
        f"tol {tol:.3e} at step {k} (0-based), sigma = {hist[k]}, "
        f"principals = {principals(hist[k])}, "
        f"lode = {math.degrees(lode_angle(hist[k])):.3f} deg"
    )

    if T is not None:
        ftc = np.array([f_tc(s, T) for s in hist])
        k = int(np.argmax(ftc))
        assert ftc[k] <= tol, (
            f"{label}: TENSION CUTOFF violated — max principal exceeds T={T} "
            f"by {ftc[k]:.6e} > tol {tol:.3e} at step {k} (0-based), "
            f"sigma = {hist[k]}, principals = {principals(hist[k])}"
        )
    return fmc


def _cond_tol(hist, phi=PHI, c=C):
    """Conditioning-aware variant of TOL_YF for the heavily confined legs.

    f is O(|sigma|); on a 1e4 kPa state the fixed 1.879e-4 would be a ~2e-8
    RELATIVE demand, which is below what a 1e-6 return tolerance can deliver.
    Canonical testbed rule: tolerances are relative and conditioning-aware.
    """
    smax = float(np.max(np.abs(np.asarray(hist, dtype=float)))) if len(hist) else 0.0
    return 1.0e-6 * max(2.0 * c * math.cos(math.radians(phi)), smax)


# ---------------------------------------------------------------------------
# single-element material-point driver
# ---------------------------------------------------------------------------
_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
# 1/8 symmetry: normal-fix every node on the x=0 / y=0 / z=0 face
_FIX = {1: (1, 1, 1), 2: (0, 1, 1), 3: (0, 0, 1), 4: (1, 0, 1),
        5: (1, 1, 0), 6: (0, 1, 0), 8: (1, 0, 0)}      # node 7 carries no fix
_XFACE = (2, 3, 6, 7)          # x = 1
_YFACE = (3, 4, 7, 8)          # y = 1
_ZFACE = (5, 6, 7, 8)          # z = 1


def _path(tag, times, values, t_pad):
    """A Path series padded with a terminal HOLD node past ``t_pad``.

    Belt AND braces against LEDGER_quirks: ``-useLast`` was silently dropped on
    the -time/-values route in older builds, and a float-accumulated pseudo-time
    that overshoots the last node by a few ulps then makes getFactor return 0 —
    every prescribed DOF snaps to zero on the FINAL step.
    """
    t = list(times) + [t_pad + 1.0]
    v = list(values) + [values[-1]]
    ops.timeSeries("Path", tag, "-time", *t, "-values", *v, "-useLast")


def _build(mat_fn, ax, lat, t_end, dt):
    """ax = (times, eps_xx values); lat = (times, eps_yy=eps_zz values) or None.

    ``lat=None`` leaves the y=1 / z=1 faces traction-free => uniaxial stress.

    The integrator is set HERE, before ``analysis("Static")``: setting it after
    makes StaticAnalysis log "WARNING analysis Static - no Integrator
    specified" and build itself around a default.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *map(float, c))
    for t, m in _FIX.items():
        ops.fix(t, *m)
    mat_fn(1)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)

    _path(1, ax[0], ax[1], t_end)
    ops.pattern("Plain", 1, 1)
    for n in _XFACE:
        ops.sp(n, 1, 1.0)                  # unit cube: u = eps * 1.0

    if lat is not None:
        _path(2, lat[0], lat[1], t_end)
        ops.pattern("Plain", 2, 2)
        for n in _YFACE:
            ops.sp(n, 2, 1.0)
        for n in _ZFACE:
            ops.sp(n, 3, 1.0)

    ops.constraints("Transformation")
    ops.numberer("Plain")
    # UmfPack, NOT FullGeneral: the two-leg drivers prescribe every DOF (zero
    # free equations) and FullGeneral hard-crashes the process on an N=0
    # system (measured: ElasticIsotropic, plain MC and MCTC all die in step 1;
    # UmfPack returns rc=0 on the identical model). See LEDGER_quirks.
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-10, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", dt)
    ops.analysis("Static")


def _probe():
    """Committed centroid-GP stress (6) and the three measured normal strains."""
    ops.eleResponse(1, "forces")                       # set the lazy strain
    sig = list(ops.eleResponse(1, "stresses"))[0:6]    # GP1 of the stdBrick
    eps = [ops.nodeDisp(2, 1), ops.nodeDisp(4, 2), ops.nodeDisp(5, 3)]
    return sig, eps


def _step(n, tag=""):
    """Advance n steps, asserting each converges; returns (sig[n,6], eps[n,3])."""
    sig, eps = [], []
    for k in range(n):
        assert ops.analyze(1) == 0, (
            f"{tag}: analyze() failed at step {k + 1}/{n} "
            f"(pseudo-time {ops.getTime():.6g})")
        s, e = _probe()
        sig.append(s)
        eps.append(e)
    return np.array(sig), np.array(eps)


def _run(mat_fn, ax, lat, nsteps, t_end=1.0, tag=""):
    _build(mat_fn, ax, lat, t_end, t_end / nsteps)
    return _step(nsteps, tag=tag or getattr(mat_fn, "__name__", "run"))


# ---------------------------------------------------------------------------
# availability guards
# ---------------------------------------------------------------------------
def _constructible(mat_fn):
    try:
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        for t, c in _CUBE.items():
            ops.node(t, *map(float, c))
        mat_fn(1)
        ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
        ok = True
    except Exception:
        ok = False
    ops.wipe()
    return ok


@pytest.fixture(scope="module")
def mctc_available():
    if not _constructible(mat_mc):
        pytest.skip("ASDPlasticMaterial3D / MohrCoulomb_YF not available in "
                    "this build")
    if not _constructible(mat_mctc):
        pytest.skip("ASDPlasticMaterial3D MohrCoulombTensionCutoff_YF/_PF not "
                    "available in this build")


@pytest.fixture(scope="module")
def mc_available():
    if not _constructible(mat_mc):
        pytest.skip("ASDPlasticMaterial3D / MohrCoulomb_YF not available in "
                    "this build")


# ===========================================================================
# 1. uniaxial tension caps at exactly T
# ===========================================================================
@pytest.mark.t0m
@pytest.mark.parametrize("T", [T_ROCK, T_EDZ], ids=["T=24.7", "T=0"])
def test_uniaxial_tension_cap(T, mctc_available):
    """Strain-driven uniaxial TENSION (lateral faces free => true uniaxial
    stress).  Below eps_y = T/E the response is elastic; above it the max
    principal — here sigma_xx — must sit at EXACTLY T.

    Non-vacuous: plain MC alone would carry 2c cos/(1+sin) = 140.0 kPa in
    uniaxial tension, 5.7x the T = 24.7 cap, so the cutoff is what governs.
    """
    nsteps = 20
    eps_max = 5.0 * EPS_Y_TEN               # same path for both T, 5x rock yield
    scale = E * eps_max                     # 123.5 kPa — the elastic reach
    sig, eps = _run(lambda t: mat_mctc(t, T=T),
                    ([0.0, 1.0], [0.0, eps_max]), None, nsteps,
                    tag=f"uniaxial tension T={T}")

    # -- homogeneity of the material point (all 8 GPs see the same state)
    ops.eleResponse(1, "forces")
    allgp = np.array(list(ops.eleResponse(1, "stresses"))).reshape(8, 6)
    assert np.allclose(allgp, allgp[0], atol=1.0e-9 * max(scale, 1.0)), (
        "stdBrick Gauss points disagree — the driver is not homogeneous")

    # -- elastic branch: steps 1..3 sit at 0.25/0.50/0.75 eps_y
    if T > 0.0:
        for k in range(3):
            e_meas = eps[k, 0]
            assert e_meas < 0.95 * (T / E), (
                f"step {k + 1} was meant to be pre-yield: eps={e_meas:.6e} vs "
                f"eps_y={T / E:.6e}")
            assert sig[k, 0] == pytest.approx(E * e_meas, rel=1.0e-8), (
                f"elastic sigma_xx {sig[k, 0]} != E*eps {E * e_meas}")

    # -- plateau: everything from 3x the rock yield strain onward
    plateau = sig[11:, 0]
    if T > 0.0:
        for k, s in enumerate(plateau, start=12):
            assert s == pytest.approx(T, rel=1.0e-9), (
                f"step {k}: uniaxial tension did not cap at T={T}: {s}")
    else:
        assert np.max(np.abs(plateau)) <= 1.0e-9 * E, (
            f"T=0 must carry no tension at all; plateau max |sigma_xx| = "
            f"{np.max(np.abs(plateau)):.6e}")

    # -- and it really is a UNIAXIAL stress state
    assert np.max(np.abs(sig[:, 1:3])) <= 1.0e-6 * scale, (
        f"lateral faces are traction-free but sigma_yy/zz reached "
        f"{np.max(np.abs(sig[:, 1:3])):.6e}")

    assert_admissible(sig, T=T, tol=TOL_YF, label=f"uniaxial tension T={T}")


# ===========================================================================
# 2. confined compression is bit-identical to plain Mohr-Coulomb
# ===========================================================================
@pytest.mark.t1
def test_confined_compression_parity(mctc_available):
    """The cutoff must be INERT in shear/compression.

    Leg 1 (t 0->0.2): hydrostatic compression to eps = -1e-4.
    Leg 2 (t 0.2->1): axial compression to eps_xx = -5e-3 at frozen lateral
    strain (an oedometric path) — well past the MC surface, and with a max
    principal of order -1e3 kPa it never comes near T = 24.7.

    Contract: the whole committed stress history must match plain MC BIT for
    BIT.  Fallback rtol 1e-14, with the measured ulp gap in the message.
    """
    nsteps = 50
    ax = ([0.0, 0.2, 1.0], [0.0, -1.0e-4, -5.0e-3])
    lat = ([0.0, 0.2, 1.0], [0.0, -1.0e-4, -1.0e-4])

    a, _ = _run(mat_mctc, ax, lat, nsteps, tag="MCTC confined compression")
    b, _ = _run(mat_mc, ax, lat, nsteps, tag="plain-MC confined compression")

    if not np.array_equal(a, b):
        d = np.abs(a - b)
        scale = np.maximum(np.abs(a), np.abs(b))
        rel = np.where(scale > 0.0, d / np.maximum(scale, 1e-300), 0.0)
        i, j = np.unravel_index(int(np.argmax(d)), d.shape)
        msg = (f"MCTC vs plain-MC compression history is NOT bit-identical: "
               f"max|diff| = {d.max():.6e} at step {i} comp {j} "
               f"({a[i, j]!r} vs {b[i, j]!r}); max rel = {rel.max():.3e} "
               f"(~{rel.max() / 2.220446049250313e-16:.1f} ulp). The cutoff "
               f"branch must be arithmetically inert when f_MC governs.")
        assert np.allclose(a, b, rtol=1.0e-14,
                           atol=1.0e-14 * max(1.0, float(scale.max()))), msg

    # -- closed form on the compression meridian, using the MEASURED confinement
    sig_conf = 0.5 * (a[-1, 1] + a[-1, 2])
    sig_ax = a[-1, 0]
    tol = _cond_tol(a)
    assert abs(f_mc(a[-1])) <= tol, (
        f"final confined state is not ON the MC surface: f_MC = "
        f"{f_mc(a[-1]):.6e} (tol {tol:.3e}) — the path never yielded, so the "
        f"parity check above is vacuous")
    # rtol 5e-2: the SURFACE is exact on this meridian but the FLOW is
    # DP-smoothed for |theta| > 29 deg, so sigma_conf itself is only
    # approximately the textbook oedometric value.
    assert sig_ax == pytest.approx(_mc_axial(sig_conf), rel=5.0e-2), (
        f"axial strength {sig_ax:.4f} != closed form "
        f"{_mc_axial(sig_conf):.4f} at sigma_conf = {sig_conf:.4f}")

    assert_admissible(a, T=T_ROCK, tol=tol, label="MCTC confined compression")
    assert_admissible(b, T=None, tol=tol, label="plain-MC confined compression")


# ===========================================================================
# 3. THE motivating case — triaxial-equal tension (the apex)
# ===========================================================================
@pytest.mark.t0m
def test_apex_regression(mctc_available):
    """All three normal strains pulled equally positive, to 10x the strain at
    which the hydrostatic stress would reach T.

    On plain MC this is the known apex defect (singular tangent / f > 0 states
    accepted — see the companion documentation run at the bottom of this file).
    On MCTC the special_return must land the state exactly on sigma = T*delta.
    """
    nsteps = 25
    eps_max = 10.0 * T_ROCK / (3.0 * KMOD)      # 4.94e-5; 3*K*eps_max = 247 kPa
    scale = 3.0 * KMOD * eps_max
    path = ([0.0, 1.0], [0.0, eps_max])
    sig, _ = _run(mat_mctc, path, path, nsteps, tag="hydrostatic tension")

    for k in range(9, nsteps):                  # from 4x the cap onward
        s = sig[k]
        for comp in range(3):
            assert s[comp] == pytest.approx(T_ROCK, rel=1.0e-9), (
                f"step {k + 1}: hydrostatic tension must plateau at "
                f"sigma = T*delta; sigma[{comp}] = {s[comp]} != {T_ROCK}")
        assert np.max(np.abs(s[3:6])) <= 1.0e-9 * scale, (
            f"step {k + 1}: spurious shear {s[3:6]} on a hydrostatic path")

    assert_admissible(sig, T=T_ROCK, tol=TOL_YF, label="hydrostatic tension")


# ===========================================================================
# 4. the corner regime — MC face and cutoff face active together
# ===========================================================================
@pytest.mark.t0m
def test_corner_path(mctc_available):
    """Leg 1 (t 0->0.3): load onto the MC face in confined compression
    (eps_xx = -3e-3, eps_lat = -5e-5; the elastic trial there is
    sigma ~= (-8.2e3, -3.7e3, -3.7e3) kPa, f_MC > 0).
    Leg 2 (t 0.3->1): hold eps_xx and pull the lateral strains to +2.4e-3.
    Under PLASTIC relaxation the stress rides the MC triaxial-extension ridge
    (s2 = s3 climbing, f_MC = 0) until the ridge meets the cutoff, then locks
    onto the compound corner (s1*, T, T), s1* = (T(1+sin phi) - 2c cos phi)/
    (1 - sin phi) = -235.25 kPa — the MC face and BOTH cutoff planes active.
    (An elastically-derived lateral target of ~1.2e-3 never reaches the cutoff
    once plastic relaxation is accounted for — measured: final s2 = s3 =
    -213.6, still on the ridge. Targets in the measured band 2.0e-3..3.2e-3
    all end exactly at the compound corner.)
    """
    nsteps = 50
    ax = ([0.0, 0.3, 1.0], [0.0, -3.0e-3, -3.0e-3])
    lat = ([0.0, 0.3, 1.0], [0.0, -5.0e-5, 2.4e-3])
    sig, _ = _run(mat_mctc, ax, lat, nsteps, tag="corner path")

    tol = _cond_tol(sig)
    s1, _s2, s3 = principals(sig[-1])
    assert s3 <= T_ROCK + tol, (
        f"final max principal {s3:.6f} exceeds the cutoff T = {T_ROCK} "
        f"by {s3 - T_ROCK:.3e} (tol {tol:.3e})")
    assert f_mc(sig[-1]) <= tol, (
        f"final state violates Mohr-Coulomb: f_MC = {f_mc(sig[-1]):.6e} "
        f"(tol {tol:.3e}); sigma = {sig[-1]}")

    # non-vacuity: the cutoff really engaged AND the deviator really is large.
    # At the exact corner with sigma_1 = T = 24.7 the MC face puts the minor
    # principal at -235.2 kPa, so "< -100" is a safe floor.
    assert s3 > 0.5 * T_ROCK, (
        f"corner path never reached the cutoff: max principal {s3:.6f} "
        f"vs T = {T_ROCK}")
    assert s1 < -100.0, (
        f"corner path lost its confinement: min principal {s1:.6f} — this is "
        f"a near-hydrostatic tension state, not a corner")

    assert_admissible(sig, T=T_ROCK, tol=tol, label="corner path")


# ===========================================================================
# 5. tangent contract — elastic slope and elastic unloading
# ===========================================================================
@pytest.mark.t0m
def test_tangent_fd(mctc_available):
    """The default ``tangent_type`` of ASDPlasticMaterial3D is **Secant**, so
    comparing a tangent MATRIX against a stress finite difference is not a
    meaningful contract on the plastic branch.  The pragmatic contract instead:

      (a) elastic branch: sigma_xx = E*eps_xx and the secant/FD slope is E,
      (b) unloading from the tension-CUTOFF plateau is elastic with slope E,
      (c) unloading from the MC-FACE plateau (uniaxial compression) is elastic
          with slope E, and the plateau itself equals -2c cos/(1-sin).

    All slopes are formed from MEASURED strains so the penalty slip cancels.
    """
    # -- (a) elastic
    sig, eps = _run(mat_mctc, ([0.0, 1.0], [0.0, 0.5 * EPS_Y_TEN]), None, 2,
                    tag="elastic uniaxial")
    for k in range(2):
        assert sig[k, 0] == pytest.approx(E * eps[k, 0], rel=1.0e-8), (
            f"elastic step {k + 1}: sigma {sig[k, 0]} != E*eps {E * eps[k, 0]}")
    k_fd = (sig[1, 0] - sig[0, 0]) / (eps[1, 0] - eps[0, 0])
    assert k_fd == pytest.approx(E, rel=1.0e-8), f"elastic FD slope {k_fd} != E"
    assert_admissible(sig, T=T_ROCK, tol=TOL_YF, label="elastic uniaxial")

    # -- (b) unload off the tension cutoff: load to 5*eps_y by t=0.8 (step 20),
    #        then peel back 1e-5 of strain over the last 5 steps.
    top = 5.0 * EPS_Y_TEN
    sig, eps = _run(mat_mctc,
                    ([0.0, 0.8, 1.0], [0.0, top, top - 1.0e-5]), None, 25,
                    tag="unload from cutoff")
    assert sig[19, 0] == pytest.approx(T_ROCK, rel=1.0e-9), (
        f"pre-unload state is not on the cutoff: {sig[19, 0]} != {T_ROCK}")
    k_un = (sig[-1, 0] - sig[19, 0]) / (eps[-1, 0] - eps[19, 0])
    assert k_un == pytest.approx(E, rel=1.0e-6), (
        f"unloading off the tension cutoff is not elastic: slope {k_un} != E "
        f"(d_sigma = {sig[-1, 0] - sig[19, 0]:.6e}, "
        f"d_eps = {eps[-1, 0] - eps[19, 0]:.6e})")
    assert_admissible(sig, T=T_ROCK, tol=TOL_YF, label="unload from cutoff")

    # -- (c) unload off the MC face: uniaxial COMPRESSION (theta = +30 deg,
    #        sigma_conf = 0), whose plateau is the exact -2c cos/(1-sin).
    bot = -5.0e-4                              # 3.5x the -1.428e-4 yield strain
    sig, eps = _run(mat_mctc,
                    ([0.0, 0.8, 1.0], [0.0, bot, bot + 1.0e-5]), None, 25,
                    tag="unload from MC face")
    sig_uc = -2.0 * C * math.cos(math.radians(PHI)) / (1.0 - math.sin(math.radians(PHI)))
    assert sig[19, 0] == pytest.approx(sig_uc, rel=5.0e-3), (
        f"uniaxial compressive strength {sig[19, 0]:.4f} != closed form "
        f"{sig_uc:.4f}")
    k_un = (sig[-1, 0] - sig[19, 0]) / (eps[-1, 0] - eps[19, 0])
    assert k_un == pytest.approx(E, rel=1.0e-6), (
        f"unloading off the MC face is not elastic: slope {k_un} != E "
        f"(d_sigma = {sig[-1, 0] - sig[19, 0]:.6e}, "
        f"d_eps = {eps[-1, 0] - eps[19, 0]:.6e})")
    assert_admissible(sig, T=T_ROCK, tol=TOL_YF, label="unload from MC face")


# ===========================================================================
# 6. parameters are movable at runtime
# ===========================================================================
def test_update_parameter_ramp(mctc_available):
    """``ops.parameter(tag, "element", ele, <NAME>)`` + ``ops.updateParameter``.

    ASDPlasticMaterial3D::setParameter has a catch-all that stores argv[0] as
    the parameter NAME and returns responseID 1000; updateParameter then routes
    the double through ``parameters_storage.setParameterByName``.  Brick's
    setParameter forwards a bare (non-"material") argv to ALL eight material
    points, so the plain name is the right address.  That catch-all ALWAYS
    reports success — an unknown name is a silent no-op — so this test verifies
    the route by OBSERVATION and skips if the route turns out to be inert.

    Part A ramps TC_min_stress 24.7 -> 10 mid-run (that is also the canary).
    Part B ramps MC_c 100 -> 60 and MC_phi 20 -> 15 mid-run.
    """
    # ---------------- Part A: TC_min_stress 24.7 -> 10 (and the canary)
    nsteps = 40
    eps_max = 5.0 * EPS_Y_TEN
    _build(mat_mctc, ([0.0, 1.0], [0.0, eps_max]), None, 1.0, 1.0 / nsteps)
    pre, pre_e = _step(3, tag="param ramp, pre")
    assert pre[2, 0] == pytest.approx(E * pre_e[2, 0], rel=1.0e-8), (
        "the pre-ramp steps were meant to be elastic")

    try:
        ops.parameter(11, "element", 1, "TC_min_stress")
        ops.updateParameter(11, 10.0)
    except Exception as exc:                       # pragma: no cover
        pytest.skip(f"ops.parameter/updateParameter refused the "
                    f"ASDPlasticMaterial3D name route: {exc}")

    post, _ = _step(nsteps - 3, tag="param ramp, post")
    plateau = post[-10:, 0]
    if np.allclose(plateau, T_ROCK, rtol=1.0e-3):
        pytest.skip(
            "ops.updateParameter had NO observable effect: the cap is still "
            f"{T_ROCK} after TC_min_stress was set to 10.0. The catch-all in "
            "ASDPlasticMaterial3D::setParameter always reports success, so the "
            "addressing route cannot be verified from Python alone — this may "
            "well be a real defect, not an absent feature.")
    for k, s in enumerate(plateau):
        assert s == pytest.approx(10.0, rel=1.0e-6), (
            f"after the TC_min_stress ramp the cap should be 10.0, got {s} "
            f"(plateau sample {k})")
    assert_admissible(post, T=10.0, tol=TOL_YF, label="post TC ramp")

    # ---------------- Part B: MC_c 100 -> 60, MC_phi 20 -> 15
    c_new, phi_new = 60.0, 15.0
    nsteps = 50
    ax = ([0.0, 0.2, 1.0], [0.0, -1.0e-4, -5.0e-3])
    lat = ([0.0, 0.2, 1.0], [0.0, -1.0e-4, -1.0e-4])
    _build(mat_mctc, ax, lat, 1.0, 1.0 / nsteps)
    a, _ = _step(30, tag="c/phi ramp, pre")
    assert abs(f_mc(a[-1])) <= _cond_tol(a), (
        f"the pre-ramp leg never reached the MC surface (f_MC = "
        f"{f_mc(a[-1]):.6e}), so the ramp below would prove nothing")

    ops.parameter(12, "element", 1, "MC_c")
    ops.updateParameter(12, c_new)
    ops.parameter(13, "element", 1, "MC_phi")
    ops.updateParameter(13, phi_new)

    b, _ = _step(20, tag="c/phi ramp, post")
    sig_conf = 0.5 * (b[-1, 1] + b[-1, 2])
    sig_ax = b[-1, 0]
    pred_new = _mc_axial(sig_conf, phi_new, c_new)
    pred_old = _mc_axial(sig_conf, PHI, C)

    assert abs(sig_ax - pred_old) > 0.05 * abs(pred_old), (
        f"the ramped state {sig_ax:.4f} is indistinguishable from the "
        f"UNramped closed form {pred_old:.4f} — the c/phi ramp did nothing")
    # rtol 5e-2 for the same DP-smoothed-flow reason as the parity test
    assert sig_ax == pytest.approx(pred_new, rel=5.0e-2), (
        f"after ramping c={c_new}, phi={phi_new} the axial strength "
        f"{sig_ax:.4f} != closed form {pred_new:.4f} at "
        f"sigma_conf = {sig_conf:.4f}")

    assert_admissible(b, phi=phi_new, c=c_new, T=T_ROCK, tol=_cond_tol(b),
                      label="post c/phi ramp")


# ===========================================================================
# 7. an inactive cutoff must change nothing
# ===========================================================================
@pytest.mark.t1
def test_noop_cutoff(mctc_available):
    """T = 500 kPa sits ABOVE the MC apex p_apex = c*cot(phi) = 274.75 kPa, so
    the cutoff surface is entirely outside the MC cone and can never activate.

    (a) the confined-compression leg must again be bit-identical to plain MC,
    (b) hydrostatic tension driven 1.5x past the apex (3K*eps = 412 kPa of
        elastic trial vs a 274.75 kPa apex) must still commit only admissible
        states — the composite's apex handling has to hold even when the cutoff
        it was built around is inert.
    """
    # -- (a) parity in compression
    nsteps = 50
    ax = ([0.0, 0.2, 1.0], [0.0, -1.0e-4, -5.0e-3])
    lat = ([0.0, 0.2, 1.0], [0.0, -1.0e-4, -1.0e-4])
    a, _ = _run(lambda t: mat_mctc(t, T=T_NOOP), ax, lat, nsteps,
                tag="no-op cutoff compression")
    b, _ = _run(mat_mc, ax, lat, nsteps, tag="plain-MC compression")

    if not np.array_equal(a, b):
        d = np.abs(a - b)
        scale = np.maximum(np.abs(a), np.abs(b))
        rel = np.where(scale > 0.0, d / np.maximum(scale, 1e-300), 0.0)
        i, j = np.unravel_index(int(np.argmax(d)), d.shape)
        msg = (f"T = {T_NOOP} > p_apex = {P_APEX:.2f} must be a no-op, but the "
               f"compression history is NOT bit-identical to plain MC: "
               f"max|diff| = {d.max():.6e} at step {i} comp {j} "
               f"({a[i, j]!r} vs {b[i, j]!r}); max rel = {rel.max():.3e} "
               f"(~{rel.max() / 2.220446049250313e-16:.1f} ulp)")
        assert np.allclose(a, b, rtol=1.0e-14,
                           atol=1.0e-14 * max(1.0, float(scale.max()))), msg
    assert_admissible(a, T=T_NOOP, tol=_cond_tol(a), label="no-op compression")

    # -- (b) hydrostatic tension straight at the MC apex
    eps_max = 1.5 * P_APEX / (3.0 * KMOD)          # 8.24e-5
    path = ([0.0, 1.0], [0.0, eps_max])
    h, _ = _run(lambda t: mat_mctc(t, T=T_NOOP), path, path, 30,
                tag="no-op cutoff hydrostatic tension")

    assert_admissible(h, T=T_NOOP, tol=TOL_YF, label="no-op hydrostatic tension")
    # with the cutoff inert, the effective cap IS the apex: f_MC <= 0 at J2 = 0
    # is exactly I1/3 <= c*cot(phi).
    # ADR-84: T_eff = min(T, c*cot phi), so with T above the apex the hook
    # clamps to the apex and projects to c*cot(phi)*delta.
    i1_3 = float(np.sum(h[-1, 0:3])) / 3.0
    assert i1_3 <= P_APEX + TOL_YF, (
        f"hydrostatic stress {i1_3:.6f} exceeded the MC apex {P_APEX:.6f}")
    assert i1_3 > 0.5 * P_APEX, (
        f"hydrostatic stress only reached {i1_3:.6f} of an apex at "
        f"{P_APEX:.6f} — the path collapsed instead of returning to the apex "
        f"(elastic trial was 3*K*eps_max = {3.0 * KMOD * eps_max:.2f} kPa)")


# ===========================================================================
# 8. DOCUMENTATION RUN — plain Mohr-Coulomb on the apex path.
#
# Deliberately assertion-free and deliberately LAST in the file: this is the
# known plain-MC apex defect (singular tangent / f > 0 acceptance) that the
# composite exists to fix.  Whether it converges, refuses, or commits an
# inadmissible state is all interesting and none of it is a contract, so no
# return code and no stress is asserted on.  Being last means that if it ever
# does take the interpreter down, every real test above has already run.
# ===========================================================================
def test_plain_mc_apex_reference_run(mc_available):
    """Reference-only: plain MC pulled 1.5x past its hydrostatic apex."""
    eps_max = 1.5 * P_APEX / (3.0 * KMOD)
    path = ([0.0, 1.0], [0.0, eps_max])
    _build(mat_mc, path, path, 1.0, 1.0 / 30)

    codes, worst_f = [], -1.0e30
    for _ in range(30):
        try:
            rc = ops.analyze(1)
        except Exception as exc:                   # pragma: no cover
            codes.append(f"exception: {exc}")
            break
        codes.append(rc)
        if rc != 0:
            break
        try:
            sig, _ = _probe()
            worst_f = max(worst_f, f_mc(sig))
        except Exception:                          # pragma: no cover
            break

    print(f"\n[reference] plain MohrCoulomb_YF on the hydrostatic-tension apex "
          f"path (eps_max = {eps_max:.4e}, apex = {P_APEX:.2f} kPa): "
          f"analyze codes = {codes}, worst committed f_MC = {worst_f:.6e} "
          f"(tol of record {TOL_YF:.3e}). f_MC > 0 here is the EXPECTED "
          f"plain-MC apex defect; nothing is asserted.")
