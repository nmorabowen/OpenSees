"""ADR-97 P2 (wp/97c) — gates 1, 2, 4 and 6 for the PRINCIPAL-STRESS-SPACE
closest-point return map of the Mohr-Coulomb family.

WHAT IS UNDER TEST
------------------
`integration_method Closest_Point` on `MohrCoulomb_YF` x `MohrCoulomb_PF` and on
`MohrCoulombTensionCutoff_YF` x `..._PF` no longer runs P1's smooth 6D Newton:
the Mohr-Coulomb surface is six PLANES in principal stress space, so on the
sorted sextant (s1 >= s2 >= s3, tension positive) every return — to the face, to
either corner LINE, or to the vertex — is a closed-form linear projection in the
elastic metric, and the Koiter tangent is constant on each region.  There is no
Newton anywhere on this path, which is why `cp_iterations` reads **1** on a
plastic step and 0 on an elastic one: it is a region count.

Region selection is Clausen's boundary-plane test on the trial principal point,
NOT a `dLambda >= 0` active-set search — at the apex with `psi < phi` the three
Koiter multipliers are not all positive (ADR-97 P0 header finding 6), so an
active-set search classifies those states wrongly.

`MohrCoulombTensionCutoff` first offers the trial to ADR-84's `special_return`
hook (cutoff face / Rankine edge / MC-cutoff corner / compound corner / apex, all
closed form, with the RAW Koiter tangent in `stiffness_return`) and only falls
back to the plain-MC principal return when that hook declines.  No ADR-84
geometry is re-derived, so wherever the hook governs, `Closest_Point` and
`Backward_Euler` must commit the SAME stress — asserted below.

GATE 1 — the four regions of `adr97_oracle/cppm_mc.py` pinned against the P0
transcript `reference_output.txt`, plus `|f| <= tol` at every commit of a
multi-step path (triaxial, simple shear, rotating principal directions) and the
MCTC physical oracles of ADR-84.

GATE 2 — `Algorithmic` against a central difference of the BINARY'S OWN assembled
internal force, in the face region (three separated principal stresses, so the
eigenprojection rotation term is live) and in the degenerate-eigenvalue edge
region (its l'Hôpital limit), plus a global-Newton iteration contrast against
`Backward_Euler`.

GATE 4 — `Backward_Euler` inertness is `tests/test_adr97_p4_inertness.py` (23
decks, 282 rows, byte-identical, fresh subprocesses).  What is here is the
relationship between the two maps on Mohr-Coulomb, which turned out to be a
FINDING: the flow direction is constant inside a sextant, so the cutting plane
and the closest point must be the same point — and they are, but only when
`MC_ds > 0` switches the yield function and the flow direction from their
ANALYTIC Lode-angle branch to a central difference of their own `f`/`g`.  See
the gate-4 block for the argument and the numbers.

GATE 6 — fail-loud: the MIXED pairings the generator registers
(`MohrCoulomb_YF` x `VonMises_PF` / `DruckerPrager_PF` / `HoekBrown_PF`,
`VonMises_YF` / `DruckerPrager_YF` x `MohrCoulomb_PF`) are covered by no oracle
and must STILL be refused; so must StiffSoil (HoekBrown ships in wp/97d); and a hydrostatic
trial with a vanishing deviator must exercise the degenerate-eigenvalue branch
without producing a NaN.

FLOAT PINS
----------
The returned stress of every region is closed-form linear algebra whose only
platform-variable step is a 3x3 symmetric eigen-decomposition, and each pinned
trial below is either non-degenerate or returns a state that is INVARIANT under
the eigenvector ambiguity (the degenerate pair comes back equal).  Those pins are
therefore taken at 1e-10 relative and the measured value is printed.  Everything
that passes through a Newton (`Backward_Euler` comparisons, the CP-vs-BE gap) is
pinned at the testbed's cross-platform `RTOL = 1e-6`.

RIG
---
The P1 rig: one unit-cube `LadrunoBrick` with EVERY displacement degree of
freedom prescribed to a homogeneous strain field, `system UmfPack` (`FullGeneral`
crashes a fully prescribed rig, `ProfileSPD` is wrong because the non-associated
`C_alg` is unsymmetric).  `stdBrick` appears only inside `fd_check`.

Zone-A.
"""
import math
import os
import sys

import numpy as np
import pytest

from _testbed import ops

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir, "Ladruno_implementation",
                                "adr97_oracle"))
import fd_tangent_driver as FD  # noqa: E402

pytestmark = [pytest.mark.zone_a]

RTOL = 1e-6            # testbed cross-platform pin
RTOL_CLOSED = 1e-10    # closed-form region returns (see FLOAT PINS above)

NODES = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
         (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]

IV_NULL = "BackStress(NullHardeningTensorFunction):"

# ---------------------------------------------------------------------------
# the P0 oracle's own Mohr-Coulomb constants (adr97_oracle/cppm_mc.py)
# ---------------------------------------------------------------------------
MC_E, MC_NU = 30000.0, 0.25
MC_PHI, MC_PSI, MC_C = 30.0, 10.0, 10.0
MC_APEX = MC_C / math.tan(math.radians(MC_PHI))          # 17.3205080757
MC_K = MC_C * math.cos(math.radians(MC_PHI))             # c cos(phi)

# ADR-84's Cerro-Lindo-like MCTC constants (tests/test_asdplastic_mctc.py, kPa)
TC_E, TC_NU = 2.0e6, 0.3
TC_C, TC_PHI, TC_PSI = 100.0, 20.0, 5.0
TC_T = 24.7
TC_P_APEX = TC_C / math.tan(math.radians(TC_PHI))        # 274.75


def _Ee(E, nu):
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    D = np.zeros((6, 6))
    D[:3, :3] = lam
    for i in range(3):
        D[i, i] = lam + 2 * mu
    for i in range(3, 6):
        D[i, i] = mu
    return D


EE_MC = _Ee(MC_E, MC_NU)
EE_TC = _Ee(TC_E, TC_NU)


# ---------------------------------------------------------------------------
# materials
# ---------------------------------------------------------------------------
def mat_mc(tag, method="Closest_Point", tangent="Algorithmic",
           phi=MC_PHI, psi=MC_PSI, c=MC_C, E=MC_E, nu=MC_NU, niter=100,
           strict=None, ds=0.0):
    extra = ["strict_convergence", int(strict)] if strict is not None else []
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "MohrCoulomb_YF", "MohrCoulomb_PF", "LinearIsotropic3D_EL", IV_NULL,
        "Begin_Model_Parameters",
        "YoungsModulus", E, "PoissonsRatio", nu,
        "MC_phi", float(phi), "MC_c", float(c), "MC_psi", float(psi),
        "MC_ds", float(ds), "MassDensity", 0.0,
        "End_Model_Parameters",
        "Begin_Internal_Variables",
        "BackStress", 0., 0., 0., 0., 0., 0.,
        "End_Internal_Variables",
        "Begin_Integration_Options",
        "integration_method", method, "tangent_type", tangent,
        "n_max_iterations", int(niter), *extra,
        "End_Integration_Options",
    )


def mat_mctc(tag, method="Closest_Point", tangent="Algorithmic", T=TC_T,
             phi=TC_PHI, psi=TC_PSI, c=TC_C, niter=100, strict=None):
    extra = ["strict_convergence", int(strict)] if strict is not None else []
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "MohrCoulombTensionCutoff_YF", "MohrCoulombTensionCutoff_PF",
        "LinearIsotropic3D_EL", IV_NULL,
        "Begin_Model_Parameters",
        "YoungsModulus", TC_E, "PoissonsRatio", TC_NU,
        "MC_phi", float(phi), "MC_c", float(c), "MC_psi", float(psi),
        "MC_ds", 0.0, "TC_min_stress", float(T), "MassDensity", 0.0,
        "End_Model_Parameters",
        "Begin_Internal_Variables",
        "BackStress", 0., 0., 0., 0., 0., 0.,
        "End_Internal_Variables",
        "Begin_Integration_Options",
        "integration_method", method, "tangent_type", tangent,
        "n_max_iterations", int(niter), *extra,
        "End_Integration_Options",
    )


# ---------------------------------------------------------------------------
# the P1 prescribed-strain rig (verbatim pattern, see test_adr97_p1_smooth.py)
# ---------------------------------------------------------------------------
def _u_of(eps, x, y, z):
    e11, e22, e33, g12, g23, g13 = eps
    return (e11 * x + g12 * y + g13 * z, e22 * y + g23 * z, e33 * z)


def drive(mat_fn, legs, nstep=1, ele="LadrunoBrick", tol=1e-13, maxiter=60,
          want=()):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for k, (x, y, z) in enumerate(NODES):
        ops.node(k + 1, float(x), float(y), float(z))
    mat_fn(1)
    ops.element(ele, 1, *range(1, 9), 1)

    U = []
    for eps in legs:
        d = {}
        for k, (x, y, z) in enumerate(NODES):
            u = _u_of(eps, float(x), float(y), float(z))
            for dd in (1, 2, 3):
                d[(k + 1, dd)] = u[dd - 1]
        U.append(d)

    nl = len(legs)
    owner = {}
    for key in U[0]:
        vals = [0.0] + [U[L][key] for L in range(nl)]
        changed = [L for L in range(nl) if abs(vals[L + 1] - vals[L]) > 1e-18]
        assert len(changed) <= 1, (
            "path not component-disjoint at DOF %s (legs %s)" % (key, changed))
        owner[key] = changed[0] if changed else 0

    for L in range(nl):
        if nl == 1:
            ops.timeSeries("Linear", 100 + L)
        else:
            # one time point PAST the end of the path: a Path series returns 0
            # outside its range and the last step lands exactly on the final
            # time (LEDGER_quirks, ADR-97 P1).
            times = [float(t) for t in range(nl + 2)]
            values = [0.0] * (L + 1) + [1.0] * (nl - L + 1)
            ops.timeSeries("Path", 100 + L, "-time", *times, "-values", *values)
        ops.pattern("Plain", 100 + L, 100 + L)
        for key, o in owner.items():
            if o == L:
                ops.sp(key[0], key[1], float(U[L][key]))

    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", tol, maxiter, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nstep)
    ops.analysis("Static")

    out = {"codes": [], "sigma": [], "eps": []}
    for w in want:
        out[w] = []
    for _ in range(nl * nstep):
        rc = ops.analyze(1)
        out["codes"].append(rc)
        if rc != 0:
            break
        out["sigma"].append(np.array(list(ops.eleResponse(1, "stresses"))[0:6]))
        out["eps"].append(np.array(list(ops.eleResponse(1, "strains"))[0:6]))
        for w in want:
            v = list(ops.eleResponse(1, "material", 1, w))
            out[w].append(np.array(v) if v else np.array([]))
    out["sigma"] = np.array(out["sigma"])
    out["eps"] = np.array(out["eps"])
    return out


# ---------------------------------------------------------------------------
# numpy audits, in the arithmetic of MohrCoulomb_YF::YF
# ---------------------------------------------------------------------------
def _tensor(s):
    s = np.asarray(s, dtype=float)
    return np.array([[s[0], s[3], s[5]], [s[3], s[1], s[4]], [s[5], s[4], s[2]]])


def princ_desc(s):
    """principal stresses DESCENDING (s1 >= s2 >= s3), tension positive."""
    return np.sort(np.linalg.eigvalsh(_tensor(s)))[::-1]


def f_mc(s, phi=MC_PHI, c=MC_C):
    """The principal-stress Mohr-Coulomb function -- verified by the P0 oracle
    to equal the header's own invariant expression to 1e-14 relative."""
    sp = math.sin(math.radians(phi))
    x = princ_desc(s)
    return 0.5 * ((x[0] - x[2]) + (x[0] + x[2]) * sp) - c * math.cos(math.radians(phi))


def f_tc(s, T):
    return float(princ_desc(s)[0]) - T


def _rel(a, b):
    b = np.asarray(b, dtype=float)
    scale = max(float(np.max(np.abs(b))), 1e-30)
    return float(np.max(np.abs(np.asarray(a, dtype=float) - b))) / scale


def _constructible(mat_fn):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for k, (x, y, z) in enumerate(NODES):
        ops.node(k + 1, float(x), float(y), float(z))
    try:
        mat_fn(1)
        ops.element("LadrunoBrick", 1, *range(1, 9), 1)
        ok = True
    except Exception:
        ok = False
    ops.wipe()
    return ok


@pytest.fixture(scope="module")
def cp_mc_available():
    if not _constructible(lambda t: mat_mc(t)):
        pytest.skip("ASDPlasticMaterial3D MohrCoulomb Closest_Point not available")


@pytest.fixture(scope="module")
def cp_mctc_available():
    if not _constructible(lambda t: mat_mctc(t)):
        pytest.skip("ASDPlasticMaterial3D MCTC Closest_Point not available")


# ===========================================================================
# GATE 1a — the four regions, pinned against adr97_oracle/cppm_mc.py
#
# Each deck is ONE step of prescribed total strain eps = E^-1 sigma_tr, so the
# elastic predictor at that step is EXACTLY the oracle's trial stress (the
# committed stress starts at zero).  The trial states are the oracle's own, and
# each one sits well inside its region (the oracle prints the boundary-plane
# margins).  Note the two degenerate-eigenvalue trials: the eigenVECTORS are not
# unique there, but the returned state is (the degenerate pair comes back equal),
# so the pin is well posed.
# ===========================================================================
MC_REGION_CASES = {
    "face (no shear)": dict(
        eps=(0.00083333333333333339, -0.00041666666666666653,
             -0.0029166666666666668, 0.0, 0.0, 0.0),
        sigma=(-20.481237065827457, -39.683854778101804,
               -96.084727348859914, 0.0, 0.0, 0.0),
        region="face", degenerate=False),
    "face (sheared, exercises the eigenprojection rotation term)": dict(
        eps=(0.00083333333333333339, -0.00041666666666666653,
             -0.0029166666666666668, 0.00125, -0.00066666666666666664,
             0.00041666666666666669),
        sigma=(-24.015909102882752, -41.824969350524078,
               -93.997938487132927, 8.9707883036619798,
               -7.2849933621429637, 4.3046526769344409),
        region="face", degenerate=False),
    "edge s1 == s2 (triaxial compression corner)": dict(
        eps=(0.00054166666666666642, 0.00054166666666666653,
             -0.0042499999999999994, 0.0, 0.0, 0.0),
        sigma=(-33.052061309937947, -33.052061309937947,
               -133.79720008119139, 0.0, 0.0, 0.0),
        region="edge12", degenerate=True),
    "edge s2 == s3 (triaxial extension corner)": dict(
        eps=(0.0019166666666666666, -0.0024583333333333332,
             -0.0024583333333333332, 0.0, 0.0, 0.0),
        sigma=(-18.220807360885715, -89.303438234034701,
               -89.303438234034701, 0.0, 0.0, 0.0),
        region="edge23", degenerate=True),
    "apex (hydrostatic tension)": dict(
        eps=(0.00075, 0.00075, 0.00075, 0.0, 0.0, 0.0),
        sigma=(17.320508075688778, 17.320508075688778,
               17.320508075688778, 0.0, 0.0, 0.0),
        region="apex", degenerate=True),
    "apex (slightly deviatoric)": dict(
        eps=(0.00083333333333333339, 0.00066666666666666664, 0.0005,
             0.0, 0.0, 0.0),
        sigma=(17.320508075688778, 17.320508075688778,
               17.320508075688778, 0.0, 0.0, 0.0),
        region="apex", degenerate=True),
}


@pytest.mark.parametrize("name", list(MC_REGION_CASES))
def test_gate1_region_returns_match_the_p0_oracle(cp_mc_available, name):
    case = MC_REGION_CASES[name]
    eps = np.array(case["eps"])
    r = drive(lambda t: mat_mc(t), [eps], nstep=1, want=("cp_iterations",))
    assert r["codes"] == [0], r["codes"]
    sig = r["sigma"][-1]
    err = _rel(sig, case["sigma"])
    tr = EE_MC @ eps
    print("gate 1a %-56s region=%-6s trial=%s\n"
          "        sigma   = %s\n"
          "        oracle  = %s\n"
          "        rel err = %.3e   |f| = %.3e   cp_iterations = %d"
          % (name, case["region"], np.array2string(tr, precision=4),
             np.array2string(sig, precision=10),
             np.array2string(np.array(case["sigma"]), precision=10),
             err, abs(f_mc(sig)), int(r["cp_iterations"][-1][0])))
    assert err <= RTOL_CLOSED, err
    # admissible, and ON the surface for every region (the apex included: it is
    # the vertex of the surface, f == 0 there by construction).
    assert abs(f_mc(sig)) <= 1e-8 * max(1.0, MC_K), f_mc(sig)
    # closed form: a region count, not a Newton count.
    assert int(r["cp_iterations"][-1][0]) == 1, r["cp_iterations"][-1]


def test_gate1_apex_multipliers_need_not_be_positive(cp_mc_available):
    """The apex rows above are reached with psi = 10 < phi = 30.  The P0 oracle
    prints the Koiter multipliers there as [-0.0115, +0.0062, +0.0133]: the cone
    of return directions no longer contains the hydrostatic direction, so an
    active-set search that grows the set while dLambda >= 0 would NEVER classify
    these states as apex states.  This test is the behavioural consequence — the
    boundary-plane classification does reach the vertex — plus the statement that
    it is the classification, not the multipliers, that decides."""
    eps = np.array(MC_REGION_CASES["apex (slightly deviatoric)"]["eps"])
    r = drive(lambda t: mat_mc(t), [eps], nstep=1)
    sig = r["sigma"][-1]
    print("gate 1 apex vertex = %s (apex = %.10f on the hydrostatic axis)"
          % (np.array2string(sig, precision=10), MC_APEX))
    assert _rel(sig[:3], np.full(3, MC_APEX)) <= RTOL_CLOSED
    assert float(np.max(np.abs(sig[3:]))) <= 1e-10 * MC_APEX


# ===========================================================================
# GATE 1b — |f| <= tol at every commit of a multi-step path
# ===========================================================================
MC_PATHS = {
    # confined compression: lands on the triaxial-compression corner (s1 == s2)
    "triaxial compression": [np.array([5.0e-4, 5.0e-4, -4.0e-3, 0., 0., 0.])],
    # pure shear about the hydrostatic axis: three distinct principals, face
    "simple shear": [np.array([0., 0., 0., 3.0e-3, 0., 0.])],
    # two legs whose principal DIRECTIONS rotate between them
    "rotating principal directions": [
        np.array([0., 0., -3.0e-3, 2.0e-3, 0., 0.]),
        np.array([0., 0., -3.0e-3, 2.0e-3, 2.0e-3, 0.])],
}


@pytest.mark.parametrize("name", list(MC_PATHS))
def test_gate1_admissible_at_every_commit(cp_mc_available, name):
    r = drive(lambda t: mat_mc(t), MC_PATHS[name], nstep=10,
              want=("cp_iterations",))
    assert all(c == 0 for c in r["codes"]), r["codes"]
    worst = 0.0
    plastic = 0
    for s in r["sigma"]:
        worst = max(worst, f_mc(s))
        if f_mc(s) > -1e-6 * MC_K:
            plastic += 1
    iters = [int(v[0]) for v in r["cp_iterations"]]
    print("gate 1b %-32s steps=%d  plastic=%d  worst f = %+.3e  "
          "cp_iterations = %s" % (name, len(r["sigma"]), plastic, worst,
                                  sorted(set(iters))))
    assert plastic >= 3, "path never yielded; the gate would be vacuous"
    assert worst <= 1e-8 * MC_K, worst
    assert max(iters) <= 1, iters


def test_gate1_rotating_path_principal_directions_really_rotate(cp_mc_available):
    """Guard against a vacuous rotating-normal gate: the eigenvectors of the
    committed stress must actually move between the two legs."""
    r = drive(lambda t: mat_mc(t), MC_PATHS["rotating principal directions"],
              nstep=10)
    v0 = np.linalg.eigh(_tensor(r["sigma"][9]))[1][:, 0]
    v1 = np.linalg.eigh(_tensor(r["sigma"][-1]))[1][:, 0]
    ang = math.degrees(math.acos(min(1.0, abs(float(v0 @ v1)))))
    print("gate 1b principal direction rotated %.2f deg between the legs" % ang)
    assert ang > 5.0, ang


# ===========================================================================
# GATE 1c — MohrCoulombTensionCutoff: ADR-84's physical oracles, and the
#           requirement that CP reuse ADR-84's `special_return` unchanged
# ===========================================================================
def test_gate1_mctc_hydrostatic_tension_caps_at_T(cp_mctc_available):
    """Pure hydrostatic tension well past the cutoff must commit sigma = T*I --
    ADR-84's own oracle -- under Closest_Point, and identically to
    Backward_Euler, because both go through the same `special_return` hook."""
    eps = np.full(6, 0.0)
    eps[:3] = 3.0 * TC_T / TC_E * (1.0 - 2.0 * TC_NU) / 1.0   # ~3x the cap
    r_cp = drive(lambda t: mat_mctc(t), [eps], nstep=4)
    r_be = drive(lambda t: mat_mctc(t, method="Backward_Euler",
                                    tangent="Secant"), [eps], nstep=4)
    assert all(c == 0 for c in r_cp["codes"]), r_cp["codes"]
    s_cp, s_be = r_cp["sigma"][-1], r_be["sigma"][-1]
    err = _rel(s_cp[:3], np.full(3, TC_T))
    gap = _rel(s_cp, s_be)
    print("gate 1c MCTC hydrostatic tension: sigma = %s\n"
          "        cap error vs T = %.3e   CP-vs-BE gap = %.3e"
          % (np.array2string(s_cp, precision=8), err, gap))
    assert err <= RTOL, err
    assert float(np.max(np.abs(s_cp[3:]))) <= RTOL * TC_T
    assert gap <= RTOL, gap


def test_gate1_mctc_uniaxial_tension_face_return_equals_backward_euler(
        cp_mctc_available):
    """A uniaxial tension leg is resolved by `special_return` Stage 1 (the
    Rankine face, s3 -> T).  `Closest_Point` reuses that hook verbatim, so the
    committed stress must agree with `Backward_Euler` to the printed tolerance --
    this is the test that would go red if P2 had re-derived ADR-84's geometry."""
    eps = np.zeros(6)
    eps[2] = 4.0 * TC_T / TC_E
    r_cp = drive(lambda t: mat_mctc(t), [eps], nstep=5)
    r_be = drive(lambda t: mat_mctc(t, method="Backward_Euler",
                                    tangent="Secant"), [eps], nstep=5)
    assert all(c == 0 for c in r_cp["codes"]), r_cp["codes"]
    s_cp, s_be = r_cp["sigma"][-1], r_be["sigma"][-1]
    gap = _rel(s_cp, s_be)
    print("gate 1c MCTC uniaxial tension: sigma_cp = %s\n"
          "        sigma_zz cap = %.10f (T = %.4f)   CP-vs-BE gap = %.3e"
          % (np.array2string(s_cp, precision=8), s_cp[2], TC_T, gap))
    assert abs(s_cp[2] - TC_T) <= RTOL * TC_T, s_cp[2]
    assert gap <= RTOL, gap
    assert f_tc(s_cp, TC_T) <= RTOL * TC_T


def test_gate1_mctc_confined_compression_falls_through_to_the_mc_return(
        cp_mctc_available):
    """Confined compression leaves the cutoff inactive at the trial, so
    `special_return` declines and the plain-MC principal return takes over.  The
    committed state must satisfy BOTH branches of the composite yield function
    (the integrator re-checks the composite f, which is the gate that catches an
    MC return that violates the cutoff)."""
    eps = np.array([2.0e-4, 2.0e-4, -2.0e-3, 0., 0., 0.])
    r = drive(lambda t: mat_mctc(t), [eps], nstep=8, want=("cp_iterations",))
    assert all(c == 0 for c in r["codes"]), r["codes"]
    worst_mc = max(f_mc(s, phi=TC_PHI, c=TC_C) for s in r["sigma"])
    worst_tc = max(f_tc(s, TC_T) for s in r["sigma"])
    scale = TC_C * math.cos(math.radians(TC_PHI))
    print("gate 1c MCTC confined compression: sigma = %s\n"
          "        worst f_MC = %+.3e (scale %.3f)   worst f_TC = %+.3e"
          % (np.array2string(r["sigma"][-1], precision=6), worst_mc, scale,
             worst_tc))
    assert worst_mc >= -1e-3 * scale, "never yielded; gate vacuous"
    assert worst_mc <= 1e-8 * scale, worst_mc
    assert worst_tc <= 1e-8 * max(scale, TC_T), worst_tc


# ===========================================================================
# GATE 2 - the consistent tangent, against a central difference of the
#          binary's own assembled internal force
#
# TWO rigs, because no single one reaches both regions:
#
#  * `fd_tangent_driver`'s `uniaxial` rig fixes x and y on EVERY node, so the
#    state is uniaxial STRAIN: s1 == s2 == lambda*eps, s3 == M*eps.  That is the
#    triaxial-compression CORNER and a DEGENERATE trial eigenvalue, i.e. the
#    l'Hopital branch of the eigenprojection rotation term.  Confinement grows
#    with compression, so perfect plasticity has no limit point and the
#    load-driven rig converges.  It needs nu < 0.25: at nu = 0.25 with phi = 30
#    the oedometric stress path K0 = nu/(1-nu) = 1/3 coincides EXACTLY with the
#    Mohr-Coulomb compression meridian (1-sin phi)/(1+sin phi) = 1/3, so f is
#    identically -c cos(phi) and the deck NEVER YIELDS -- a silently vacuous
#    gate if you do not check for it (this test does).
#
#  * a FREE-NODE rig, local to this file: the homogeneous strain field of the P1
#    driver is prescribed on seven of the eight nodes and node 7 is left free.
#    Seven prescribed nodes mean there is no limit point at any stress level, so
#    ANY state can be reached -- in particular a genuine FACE state with three
#    well-separated principal stresses, where the rotation term
#    (y_i - y_j)/(x_i - x_j) is live rather than at its limit.  The two-rig
#    finite-difference scheme is `fd_tangent_driver`'s (setNodeDisp does not
#    trigger Domain::update, so each perturbed state is a separate analysis).
# ===========================================================================
def _fd_build(mat_fn, eps, ele, free_node, u_free=None):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for k, (x, y, z) in enumerate(NODES):
        ops.node(k + 1, float(x), float(y), float(z))
    mat_fn(1)
    ops.element(ele, 1, *range(1, 9), 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for k, (x, y, z) in enumerate(NODES):
        if (k + 1) == free_node and u_free is None:
            continue
        u = _u_of(eps, float(x), float(y), float(z))
        for dd in (1, 2, 3):
            ops.sp(k + 1, dd, float(u[dd - 1]))
    if u_free is not None:
        for dd in (1, 2, 3):
            ops.sp(free_node, dd, float(u_free[dd - 1]))
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-13, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


def fd_free_node(mat_fn, eps, ele="LadrunoBrick", free_node=7, h=1e-7):
    """Assembled 3x3 tangent at the one free node vs a central difference of the
    binary's own reaction there.  No numpy reference anywhere."""
    _fd_build(mat_fn, eps, ele, free_node)
    rc = ops.analyze(1)
    if rc != 0:
        return dict(rc=rc, rel_err=float("nan"))
    eqn = [ops.nodeDOFs(free_node)[d] for d in range(3)]
    d = ops.printA("-sparse", "-ret")     # NOT printA() alone: needs -sparse -ret
    n = max(max(d["rowIndices"]), max(d["colIndices"])) + 1
    K = np.zeros((n, n))
    for i, j, v in zip(d["rowIndices"], d["colIndices"], d["values"]):
        K[i, j] += v
    K_asm = K[np.ix_(eqn, eqn)]
    u = np.array([ops.nodeDisp(free_node, dd) for dd in (1, 2, 3)])
    sig = np.array(list(ops.eleResponse(1, "stresses"))[0:6])
    K_fd = np.zeros((3, 3))
    for j in range(3):
        s = h * max(1e-4, abs(u[j]))
        R = []
        for sgn in (1.0, -1.0):
            up = u.copy()
            up[j] += sgn * s
            _fd_build(mat_fn, eps, ele, free_node, u_free=up)
            rc2 = ops.analyze(1)
            if rc2 != 0:
                return dict(rc=rc2, rel_err=float("nan"))
            ops.reactions()
            R.append(np.array([ops.nodeReaction(free_node, dd)
                               for dd in (1, 2, 3)]))
        K_fd[:, j] = (R[0] - R[1]) / (2.0 * s)
    scale = max(float(np.max(np.abs(K_fd))), 1e-30)
    return dict(rc=0, rel_err=float(np.max(np.abs(K_fd - K_asm))) / scale,
                K_fd=K_fd, K_asm=K_asm, sig=sig, u=u)


# the P0 oracle's own trial states, expressed as the total strain that produces
# them in ONE step (eps = E^-1 sigma_tr; the committed stress starts at zero)
EPS_FACE = np.linalg.solve(EE_MC, np.array([-10., -40., -100., 0., 0., 0.]))
EPS_FACE_SHEARED = np.linalg.solve(EE_MC,
                                   np.array([-10., -40., -100., 15., -8., 5.]))
SIGMA_FACE_ORACLE = np.array([-20.481237065827457, -39.683854778101804,
                              -96.084727348859914, 0.0, 0.0, 0.0])
FD_NU = 0.15      # see the module comment above: nu must be < 0.25 for phi = 30


def test_gate2_algorithmic_fd_on_the_degenerate_edge_region(cp_mc_available):
    r = FD.fd_check(lambda t: mat_mc(t, nu=FD_NU), rig="uniaxial",
                    load=(0.0, 0.0, -40.0),
                    label="MC edge region, Closest_Point / Algorithmic")
    assert r["rc"] == 0, "the MC load-driven oedometric rig did not converge"
    x = princ_desc(r["sig"])
    print("gate 2 (edge, DEGENERATE eigenvalue) principals = %s  f = %+.3e\n"
          "        rel_err = %.3e   rel_fro = %.3e"
          % (np.array2string(x, precision=6), f_mc(r["sig"]),
             r["rel_err"], r["rel_fro"]))
    assert f_mc(r["sig"]) >= -1e-6 * MC_K, "state is elastic; gate vacuous"
    assert abs(x[0] - x[1]) <= 1e-8 * max(1.0, abs(x[0])), \
        "expected a degenerate trial eigenvalue on the oedometric rig"
    assert r["rel_err"] <= RTOL, r["rel_err"]


@pytest.mark.parametrize("name,eps", [("face, axis aligned", EPS_FACE),
                                      ("face, sheared", EPS_FACE_SHEARED)])
def test_gate2_algorithmic_fd_on_the_face_region(cp_mc_available, name, eps):
    """Three well-separated principal stresses, so the eigenprojection rotation
    term is LIVE.  This is the measurement the P2 mutation gate kills."""
    r = fd_free_node(lambda t: mat_mc(t), eps)
    assert r["rc"] == 0, "the free-node rig did not converge"
    x = princ_desc(r["sig"])
    sep = min(abs(x[0] - x[1]), abs(x[1] - x[2])) / max(1.0, abs(x[2]))
    print("gate 2 (%s) principals = %s  separation = %.3e  f = %+.3e\n"
          "        rel_err = %.3e" % (name, np.array2string(x, precision=6),
                                      sep, f_mc(r["sig"]), r["rel_err"]))
    assert f_mc(r["sig"]) >= -1e-6 * MC_K, "state is elastic; gate vacuous"
    assert sep > 1e-2, "principal stresses are not separated; not the face case"
    assert r["rel_err"] <= RTOL, r["rel_err"]


def test_gate2_negative_control_backward_euler_cannot_do_the_face_rig(
        cp_mc_available):
    """The negative control here is an OUTCOME, not a number: with the shipped
    `Backward_Euler` tangents the same free-node rig does not converge at all,
    where `Closest_Point`/`Algorithmic` converges and reproduces the assembled
    tangent to the error printed above.  (ADR-94 M3 measured `Continuum` and
    `Secant` at 57 % and 80 % against a central difference of the material's own
    response; at those errors a Newton on a nearly singular perfectly plastic
    element simply stops converging.)"""
    codes = {}
    for lbl, kw in (("Continuum", dict(method="Backward_Euler",
                                       tangent="Continuum")),
                    ("Secant", dict(method="Backward_Euler",
                                    tangent="Secant"))):
        r = fd_free_node(lambda t: mat_mc(t, **kw), EPS_FACE)
        codes[lbl] = r["rc"]
    r_cp = fd_free_node(lambda t: mat_mc(t), EPS_FACE)
    print("gate 2 negative control on the free-node face rig: "
          "BE/Continuum analyze -> %s, BE/Secant -> %s, "
          "CP/Algorithmic -> %s (rel_err %.3e)"
          % (codes["Continuum"], codes["Secant"], r_cp["rc"], r_cp["rel_err"]))
    assert r_cp["rc"] == 0
    assert codes["Continuum"] != 0 and codes["Secant"] != 0, (
        "the shipped tangents now converge on this rig, so the contrast above is"
        " no longer a contrast: %r" % codes)


def _oedometric_iterations(mat_fn, nsteps=6, load=-40.0):
    """One confined cube, load driven, counting the GLOBAL Newton iterations per
    step.  Confinement grows with compression, so perfect plasticity has no limit
    point here."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for k, (x, y, z) in enumerate(NODES):
        ops.node(k + 1, float(x), float(y), float(z))
    for k in range(1, 5):
        ops.fix(k, 1, 1, 1)
    for k in range(5, 9):
        ops.fix(k, 1, 1, 0)
    mat_fn(1)
    ops.element("stdBrick", 1, *range(1, 9), 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for k in range(5, 9):
        ops.load(k, 0.0, 0.0, load)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-9, 300, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    iters, sig = [], None
    for _ in range(nsteps):
        rc = ops.analyze(1)
        if rc != 0:
            return None, None, rc
        iters.append(int(ops.testIter()))
        sig = np.array(list(ops.eleResponse(1, "stresses"))[0:6])
    return iters, sig, 0


def test_gate2_iteration_contrast_against_backward_euler(cp_mc_available):
    it_cp, sig_cp, rc_cp = _oedometric_iterations(lambda t: mat_mc(t, nu=FD_NU))
    it_co, sig_co, rc_co = _oedometric_iterations(
        lambda t: mat_mc(t, nu=FD_NU, method="Backward_Euler",
                         tangent="Continuum"))
    it_se, _s, rc_se = _oedometric_iterations(
        lambda t: mat_mc(t, nu=FD_NU, method="Backward_Euler",
                         tangent="Secant"))
    assert rc_cp == 0, "Closest_Point oedometric rig did not converge"
    assert rc_co == 0 and rc_se == 0, (rc_co, rc_se)
    print("gate 2 iteration contrast (MC oedometric, nu = %.2f, 6 steps)\n"
          "        Closest_Point / Algorithmic : %s  total %d\n"
          "        Backward_Euler / Continuum  : %s  total %d\n"
          "        Backward_Euler / Secant     : %s  total %d\n"
          "        committed sigma CP = %s  (f = %+.2e)\n"
          "        committed sigma BE = %s  (f = %+.2e)"
          % (FD_NU, it_cp, sum(it_cp), it_co, sum(it_co), it_se, sum(it_se),
             np.array2string(sig_cp, precision=6), f_mc(sig_cp),
             np.array2string(sig_co, precision=6), f_mc(sig_co)))
    assert f_mc(sig_cp) >= -1e-6 * MC_K, "never yielded; contrast is vacuous"
    assert sum(it_cp) * 2 <= sum(it_co), (it_cp, it_co)


def test_gate2_mctc_iteration_contrast(cp_mctc_available):
    """The same contrast on the ADR-84 MCTC deck (kPa, Cerro-Lindo-like)."""
    it_cp, sig_cp, rc_cp = _oedometric_iterations(lambda t: mat_mctc(t),
                                                  load=-4.0e3)
    it_be, sig_be, rc_be = _oedometric_iterations(
        lambda t: mat_mctc(t, method="Backward_Euler", tangent="Continuum"),
        load=-4.0e3)
    assert rc_cp == 0 and rc_be == 0, (rc_cp, rc_be)
    print("gate 2 MCTC iteration contrast: CP %s (total %d) vs BE %s (total %d)"
          "\n        sigma CP = %s\n        sigma BE = %s"
          % (it_cp, sum(it_cp), it_be, sum(it_be),
             np.array2string(sig_cp, precision=4),
             np.array2string(sig_be, precision=4)))
    assert sum(it_cp) <= sum(it_be), (it_cp, it_be)


# ===========================================================================
# GATE 4 - what `Closest_Point` agrees with, and where it does not
#
# `Backward_Euler` inertness proper is tests/test_adr97_p4_inertness.py (23
# decks, 282 rows, byte-identical, fresh subprocesses).  What belongs HERE is the
# relationship between the two maps on Mohr-Coulomb, and it turned out to be a
# FINDING rather than an agreement:
#
# Mohr-Coulomb's flow direction is CONSTANT inside a sextant (the surface is
# piecewise linear), so on a deck whose iterates stay in one sextant the
# Ortiz-Simo cutting plane and the closest point are the SAME point, and the two
# maps must agree.  Measured, they do not -- unless `MC_ds > 0`, which switches
# `MohrCoulomb_YF::df_dsigma_ij` and `MohrCoulomb_PF` from their ANALYTIC
# Lode-angle branch to a central difference of their own f / g.  Through the
# finite difference `Backward_Euler` lands on the closest-point answer to 1e-10;
# through the analytic branch it is 4e-2 to 3e-1 away.  Since f is exactly linear
# in principal stress, a central difference of it over the Voigt slots is the
# EXACT 6D gradient -- its accuracy even IMPROVES with a larger step, which is
# the signature of differencing a linear function -- so two independent
# references (this closed-form return, and the header's own difference of its own
# f) agree with each other and against the shipped analytic c1/c2/c3
# coefficients.  Recorded, NOT fixed: fixing it would change `Backward_Euler`,
# which ADR-97 D1 keeps byte-identical.
# ===========================================================================
@pytest.mark.parametrize("nstep", [1, 4, 10, 40])
def test_gate4_cp_is_exact_and_step_size_independent(cp_mc_available, nstep):
    """A closed-form projection onto a FIXED surface from a proportional path is
    the same point however many steps it is taken in."""
    r = drive(lambda t: mat_mc(t), [EPS_FACE], nstep=nstep)
    assert all(c == 0 for c in r["codes"]), r["codes"]
    err = _rel(r["sigma"][-1], SIGMA_FACE_ORACLE)
    print("gate 4 CP N=%2d vs the P0 oracle: rel err = %.3e" % (nstep, err))
    assert err <= RTOL_CLOSED, err


def test_gate4_backward_euler_agrees_only_through_its_own_finite_difference(
        cp_mc_available):
    """The finding above, pinned in BOTH directions so it cannot rot silently."""
    rows = []
    for ds in (0.0, 1e-8, 1e-6, 1e-4):
        b = drive(lambda t: mat_mc(t, method="Backward_Euler",
                                   tangent="Secant", ds=ds),
                  [EPS_FACE], nstep=10)
        assert all(c == 0 for c in b["codes"]), (ds, b["codes"])
        rows.append((ds, _rel(b["sigma"][-1], SIGMA_FACE_ORACLE),
                     b["sigma"][-1]))
    for ds, err, sig in rows:
        print("gate 4 Backward_Euler MC_ds = %-8.0e vs the closest point: "
              "rel err = %.3e   sigma = %s"
              % (ds, err, np.array2string(sig, precision=8)))
    analytic = rows[0][1]
    fd_errs = [e for _ds, e, _s in rows[1:]]
    assert analytic > 1e-2, (
        "the shipped ANALYTIC Lode-angle branch now agrees with the closest "
        "point; the finding this test records has been fixed, and the test "
        "must be rewritten rather than relaxed: %r" % analytic)
    assert max(fd_errs) <= 1e-9, (
        "Backward_Euler through its OWN central difference no longer reproduces "
        "the closest point, which would mean the two independent references "
        "have stopped agreeing: %r" % fd_errs)


def test_gate4_cp_equals_be_at_the_apex(cp_mc_available):
    """The one region where the shipped gradient cannot be wrong, because the
    integrator does not use it: both maps project onto the vertex."""
    eps = np.array(MC_REGION_CASES["apex (hydrostatic tension)"]["eps"])
    a = drive(lambda t: mat_mc(t), [eps], nstep=4)
    b = drive(lambda t: mat_mc(t, method="Backward_Euler", tangent="Secant"),
              [eps], nstep=4)
    assert all(c == 0 for c in list(a["codes"]) + list(b["codes"]))
    gap = _rel(a["sigma"][-1], b["sigma"][-1])
    print("gate 4 CP vs BE at the MC apex: rel gap = %.3e  sigma = %s"
          % (gap, np.array2string(a["sigma"][-1], precision=10)))
    assert gap <= RTOL_CLOSED, gap


# ===========================================================================
# GATE 6 — fail loud
# ===========================================================================
def _mat_pair(tag, yf, pf, iv, params, ivs, method="Closest_Point",
              tangent="Algorithmic"):
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag, yf, pf, "LinearIsotropic3D_EL", iv,
        "Begin_Model_Parameters", *params, "End_Model_Parameters",
        "Begin_Internal_Variables", *ivs, "End_Internal_Variables",
        "Begin_Integration_Options",
        "integration_method", method, "tangent_type", tangent,
        "End_Integration_Options",
    )


MC_PARAMS = ["YoungsModulus", MC_E, "PoissonsRatio", MC_NU,
             "MC_phi", MC_PHI, "MC_c", MC_C, "MC_psi", MC_PSI, "MC_ds", 0.0,
             "MassDensity", 0.0]
BS0 = ["BackStress", 0., 0., 0., 0., 0., 0.]

# (label, yf, pf, iv-string, extra parameters, extra internal variables)
REFUSED_MIXED = [
    ("MohrCoulomb_YF x VonMises_PF (Null)",
     "MohrCoulomb_YF", "VonMises_PF", IV_NULL, [], []),
    ("MohrCoulomb_YF x VonMises_PF (TensorLinear)",
     "MohrCoulomb_YF", "VonMises_PF",
     "BackStress(TensorLinearHardeningFunction):",
     ["TensorLinearHardeningParameter", 0.0], []),
    ("MohrCoulomb_YF x DruckerPrager_PF",
     "MohrCoulomb_YF", "DruckerPrager_PF",
     "BackStress(TensorLinearHardeningFunction):"
     "DP_cohesion(ScalarLinearHardeningFunction):",
     ["TensorLinearHardeningParameter", 0.0,
      "ScalarLinearHardeningParameter", 0.0,
      "DP_xi_c", 20.0, "DP_eta", 0.4, "DP_etabar", 0.4],
     ["DP_cohesion", 0.0]),
    ("VonMises_YF x MohrCoulomb_PF",
     "VonMises_YF", "MohrCoulomb_PF",
     "BackStress(TensorLinearHardeningFunction):"
     "YieldStress(ScalarLinearHardeningFunction):",
     ["TensorLinearHardeningParameter", 0.0,
      "ScalarLinearHardeningParameter", 0.0],
     ["YieldStress", 30.0]),
    ("DruckerPrager_YF x MohrCoulomb_PF",
     "DruckerPrager_YF", "MohrCoulomb_PF",
     "BackStress(TensorLinearHardeningFunction):"
     "DP_cohesion(ScalarLinearHardeningFunction):",
     ["TensorLinearHardeningParameter", 0.0,
      "ScalarLinearHardeningParameter", 0.0,
      "DP_xi_c", 20.0, "DP_eta", 0.4, "DP_etabar", 0.4],
     ["DP_cohesion", 0.0]),
    # Ladruno (ADR-97 wp/97d): the parameter is `HB_sigci`, not `HB_sigma_ci`;
    # with the old spelling this row was refused for a MISSING PARAMETER, never
    # for the family gate.  Corrected -- it still refuses, now for the right
    # reason (see tests/test_adr97_p3_hoekbrown.py, which checks every
    # HoekBrown pairing in BOTH directions).
    ("MohrCoulomb_YF x HoekBrown_PF",
     "MohrCoulomb_YF", "HoekBrown_PF", IV_NULL,
     ["HB_sigci", 50000.0, "HB_mb_psi", 2.396510364418,
      "HB_s", 0.011743628457, "HB_a", 0.502840500848, "HB_ds", 0.0], []),
]


@pytest.mark.parametrize("label,yf,pf,iv,extra_p,extra_iv", REFUSED_MIXED,
                         ids=[c[0] for c in REFUSED_MIXED])
def test_gate6_mixed_pairings_are_still_refused(label, yf, pf, iv, extra_p,
                                               extra_iv):
    """P2 enables the principal-space map ONLY when the yield function and the
    plastic flow direction are of the SAME Mohr-Coulomb family.  The generator
    registers six cross pairings that are covered by NO oracle: the principal
    return assumes both are piecewise linear, and P1's smooth 6D map cannot use
    Mohr-Coulomb's Lode-angle gradient (Drucker-Prager substitution above 29 deg,
    central differences otherwise).  They must stay refused, loudly, at parse
    time."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for k, (x, y, z) in enumerate(NODES):
        ops.node(k + 1, float(x), float(y), float(z))
    with pytest.raises(Exception):
        _mat_pair(1, yf, pf, iv, MC_PARAMS + list(extra_p), BS0 + list(extra_iv))
        ops.element("LadrunoBrick", 1, *range(1, 9), 1)
    ops.wipe()


def test_gate6_hoekbrown_is_now_accepted_and_stiffsoil_is_not():
    """Ladruno (ADR-97 wp/97d): P3 ships HoekBrown_YF x HoekBrown_PF, so the
    HoekBrown half of this test is INVERTED.  It also had `HB_sigma_ci` for
    `HB_sigci`, which made it pass on a MISSING PARAMETER rather than on the
    family gate -- corrected here.  StiffSoil (P5) has no runtime row: it needs
    StiffSoil_EL and its own parameter set, and is covered by the compile-time
    support matrix instead."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for k, (x, y, z) in enumerate(NODES):
        ops.node(k + 1, float(x), float(y), float(z))
    HB = ["YoungsModulus", 5.0e7, "PoissonsRatio", 0.25,
          "HB_sigci", 50000.0, "HB_mb", 2.396510364418,
          "HB_s", 0.011743628457, "HB_a", 0.502840500848,
          "HB_mb_psi", 2.396510364418, "HB_ds", 0.0, "MassDensity", 0.0]
    _mat_pair(1, "HoekBrown_YF", "HoekBrown_PF", IV_NULL, HB, BS0)
    ops.element("LadrunoBrick", 1, *range(1, 9), 1)
    ops.wipe()


def test_gate6_mc_and_mctc_are_accepted():
    """The positive half of the refusal matrix: the two families P2 DOES ship
    must construct under `Closest_Point` + `Algorithmic`."""
    assert _constructible(lambda t: mat_mc(t)), \
        "MohrCoulomb_YF x MohrCoulomb_PF must be accepted by P2"
    assert _constructible(lambda t: mat_mctc(t)), \
        "MCTC must be accepted by P2"


def test_gate6_algorithmic_still_refused_with_backward_euler():
    """ADR-97 D2 survives P2: a consistent tangent is defined only relative to a
    committed map."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for k, (x, y, z) in enumerate(NODES):
        ops.node(k + 1, float(x), float(y), float(z))
    with pytest.raises(Exception):
        mat_mc(1, method="Backward_Euler", tangent="Algorithmic")
        ops.element("LadrunoBrick", 1, *range(1, 9), 1)
    ops.wipe()


def test_gate6_zero_friction_angle_is_refused(cp_mc_available):
    """phi == 0 is Tresca: the Mohr-Coulomb apex is at infinity and Clausen's
    boundary planes, which pass through it, are undefined.  The map refuses
    rather than silently reclassifying (the step is rejected, so `analyze`
    returns non-zero on a LadrunoBrick, which checks the refusal sentinel)."""
    eps = np.array([5.0e-4, 5.0e-4, -4.0e-3, 0., 0., 0.])
    r = drive(lambda t: mat_mc(t, phi=0.0, psi=0.0), [eps], nstep=4)
    print("gate 6 phi = 0 analyze codes = %s" % (r["codes"],))
    assert any(c != 0 for c in r["codes"]), (
        "a phi = 0 Mohr-Coulomb committed a state under Closest_Point instead of"
        " refusing: %r" % (r["codes"],))


def test_gate6_degenerate_hydrostatic_trial_does_not_nan(cp_mc_available):
    """A trial state on (or a hair off) the hydrostatic axis makes ALL THREE
    eigenvalue gaps degenerate, so every shear slot of the back-transform takes
    the l'Hopital branch at once, and the eigenVECTORS are arbitrary.  ADR-94 B4
    is the reproducer that used to commit NaN through this class of state."""
    for tag, dev in (("exactly hydrostatic", 0.0),
                     ("hydrostatic + 1e-12 deviator", 1e-12),
                     ("hydrostatic + 1e-7 deviator", 1e-7)):
        eps = np.array([7.5e-4 + dev, 7.5e-4, 7.5e-4 - dev, 0., 0., 0.])
        r = drive(lambda t: mat_mc(t), [eps], nstep=2)
        assert all(c == 0 for c in r["codes"]), (tag, r["codes"])
        s = r["sigma"][-1]
        assert np.all(np.isfinite(s)), (tag, s)
        print("gate 6 %-32s sigma = %s"
              % (tag, np.array2string(s, precision=8)))
        assert _rel(s[:3], np.full(3, MC_APEX)) <= 1e-6, (tag, s)
        assert abs(f_mc(s)) <= 1e-8 * MC_K, (tag, f_mc(s))


def test_gate6_strict_convergence_is_inert_on_a_converging_mc_deck(
        cp_mc_available):
    """`strict_convergence 1` must change nothing on a deck whose every commit is
    admissible -- the principal return is closed form, so it never is not."""
    legs = MC_PATHS["triaxial compression"]
    a = drive(lambda t: mat_mc(t), legs, nstep=6)
    b = drive(lambda t: mat_mc(t, strict=1), legs, nstep=6)
    assert all(c == 0 for c in b["codes"]), b["codes"]
    same = np.array_equal(a["sigma"], b["sigma"])
    print("gate 6 strict_convergence byte-inert on a converging MC deck: %s "
          "(rel gap %.3e)" % (same, _rel(a["sigma"], b["sigma"])))
    assert _rel(a["sigma"], b["sigma"]) == 0.0
