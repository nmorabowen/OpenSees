"""ADR-97 P3 (`wp/97d-cp-hoekbrown`) -- `integration_method Closest_Point` +
`tangent_type Algorithmic` for the HOEK-BROWN family.

The map is a principal-stress-space return to a CURVED surface (Clausen &
Damkilde 2008): a 4x4 Newton in the surface's own variable on the face or on
either curved edge, and a closed-form vertex at the apex.  It reuses P2's
principal machinery verbatim (spectral decomposition, eigenprojection
back-transform with the rotation term and its l'Hopital limit, analytic
`Rs^-1`, tangent policy, strict-convergence contract) and replaces only the
projection.

WHAT THE ORACLE IS.  `Ladruno_implementation/adr97_oracle/cppm_hb.py`, a numpy
closest-point solver for the same surface whose region classification is
CLOSED by a 400-trial KKT scan (max return-direction residual 1.672e-14 with
every Koiter multiplier >= 0).  The elastic domain is convex, so a vanishing
return-direction residual with non-negative multipliers is a PROOF of global
optimality, not a spot check.

FLOAT PINS.  The oracle-region returns are pinned at 1e-10 RELATIVE, tighter
than the fork's usual 1e-6 cross-platform rule, for the same reason P2's are:
the only platform-variable step is a 3x3 symmetric eigen-decomposition, and
every pinned trial is either non-degenerate or returns a state INVARIANT under
the eigenvector ambiguity (the degenerate pair comes back equal).  The
transcription check that preceded the build reproduced all seven regions at
<= 1.9e-13 relative with the same Newton counts.  Everything that goes through
the global Newton or a finite difference is pinned at 1e-6.

THE POTENTIAL.  `HoekBrown_PF::g` is evaluated in the un-negated frame and so
collapses to a TRESCA potential (`HB_mb_psi` inert, zero dilatancy, 32.86 deg
off the normal even at mb_psi == mb).  `Closest_Point` uses a frame-consistent
Hoek-Brown potential in its own code path; the shipped `g` is untouched so
`Backward_Euler` stays byte-identical (ADR-97 D1).  Gate 4 PINS that difference
rather than hiding it -- CP and BE are NOT expected to agree on a Hoek-Brown
deck, and the test says by how much.
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
RTOL_REGION = 1e-10    # oracle-region returns (see FLOAT PINS above)

NODES = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
         (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]

IV_NULL = "BackStress(NullHardeningTensorFunction):"

# ---------------------------------------------------------------------------
# the P0 oracle's Hoek-Brown constants -- the SAME material as
# tests/test_adr94_hlist_hb.py (Hoek & Brown 2018 GSI formulas, kPa)
# ---------------------------------------------------------------------------
HB_SIGCI, HB_MI, HB_GSI, HB_D = 50000.0, 10.0, 60.0, 0.0
HB_MB = HB_MI * math.exp((HB_GSI - 100.0) / (28.0 - 14.0 * HB_D))
HB_S = math.exp((HB_GSI - 100.0) / (9.0 - 3.0 * HB_D))
HB_A = 0.5 + (1.0 / 6.0) * (math.exp(-HB_GSI / 15.0) - math.exp(-20.0 / 3.0))
HB_T = HB_S * HB_SIGCI / HB_MB              # 245.0151818950, == APEX_STRESS
HB_SCALE = HB_SIGCI * HB_S ** HB_A          # 5350.4268083793, YF_STRENGTH_SCALE
HB_E, HB_NU = 5.0e7, 0.25


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


EE_HB = _Ee(HB_E, HB_NU)


def mat_hb(tag, method="Closest_Point", tangent="Algorithmic",
           mb_psi=None, niter=100, strict=None, ds=0.0):
    extra = ["strict_convergence", int(strict)] if strict is not None else []
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "HoekBrown_YF", "HoekBrown_PF", "LinearIsotropic3D_EL", IV_NULL,
        "Begin_Model_Parameters",
        "YoungsModulus", HB_E, "PoissonsRatio", HB_NU,
        "HB_sigci", HB_SIGCI, "HB_mb", HB_MB, "HB_s", HB_S, "HB_a", HB_A,
        "HB_mb_psi", float(HB_MB if mb_psi is None else mb_psi),
        "HB_ds", float(ds), "MassDensity", 0.0,
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
# the P1/P2 prescribed-strain rig (verbatim pattern)
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
    ops.system("UmfPack")          # NOT FullGeneral: it crashes prescribed rigs
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
# numpy audits, in the arithmetic of HoekBrown_YF::YF (tension positive)
# ---------------------------------------------------------------------------
def _tensor(s):
    s = np.asarray(s, dtype=float)
    return np.array([[s[0], s[3], s[5]], [s[3], s[1], s[4]], [s[5], s[4], s[2]]])


def princ_desc(s):
    """principal stresses DESCENDING (y1 >= y2 >= y3), tension positive."""
    return np.sort(np.linalg.eigvalsh(_tensor(s)))[::-1]


def f_hb(s):
    """the header's own composite max(f_shear, f_tension), in principal space."""
    y = princ_desc(s)
    arg = max(HB_S - HB_MB * y[0] / HB_SIGCI, 0.0)
    return max(y[0] - y[2] - HB_SIGCI * arg ** HB_A, y[0] - HB_T)


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
def cp_hb_available():
    if not _constructible(lambda t: mat_hb(t)):
        pytest.skip("ASDPlasticMaterial3D HoekBrown Closest_Point not available")


# ===========================================================================
# GATE 1a -- the four regions, pinned against adr97_oracle/cppm_hb.py
#
# Each deck is ONE step of prescribed total strain eps = E^-1 sigma_tr, so the
# elastic predictor at that step is EXACTLY the oracle's trial stress (the
# committed stress starts at zero).  `sigma` below is Q diag(y_oracle) Q^T with
# Q the trial eigenvectors -- for the sheared trial all three eigenvalues are
# distinct, so that reconstruction is unique.
# ===========================================================================
HB_REGION_CASES = {
    "face (no shear)": dict(
        eps=(0.000115, 1.5000000000000016e-05, -0.00046, 0.0, 0.0, 0.0),
        sigma=(-3642.7736831410998, -6441.6222402965004, -25123.715278045001,
               0.0, 0.0, 0.0),
        region="face", its_max=5, degenerate=False),
    "face (sheared, exercises the eigenprojection rotation term)": dict(
        eps=(0.000115, 1.5000000000000016e-05, -0.00046,
             7.4999999999999993e-05, -4.0000000000000003e-05,
             2.5000000000000001e-05),
        sigma=(-4010.3137595367621, -6738.4186262415187, -25163.13544737453,
               1026.2415898229497, -786.30508932484497, 475.42935678468655),
        region="face", its_max=5, degenerate=False),
    "edge y1 == y2 (triaxial compression corner)": dict(
        eps=(8.9999999999999965e-05, 8.9999999999999979e-05,
             -0.00055999999999999995, 0.0, 0.0, 0.0),
        sigma=(-5028.7623527101996, -5028.7623527101996, -30069.032969829899,
               0.0, 0.0, 0.0),
        region="edge12", its_max=5, degenerate=True),
    "edge y2 == y3 (triaxial extension corner)": dict(
        eps=(0.00016000000000000001, -0.00028999999999999995, -0.00029,
             0.0, 0.0, 0.0),
        sigma=(-2395.1014970300998, -20077.293773552599, -20077.293773552599,
               0.0, 0.0, 0.0),
        region="edge23", its_max=5, degenerate=True),
    "apex (hydrostatic tension)": dict(
        eps=(3.9999999999999998e-06, 3.9999999999999998e-06,
             3.9999999999999998e-06, 0.0, 0.0, 0.0),
        sigma=(245.01518189499933, 245.01518189499933, 245.01518189499933,
               0.0, 0.0, 0.0),
        region="apex", its_max=1, degenerate=True),
    "apex (deviatoric, still in the elastic-metric cone)": dict(
        eps=(5.75e-06, 4.5000000000000001e-06, 3.2499999999999998e-06,
             0.0, 0.0, 0.0),
        sigma=(245.01518189499933, 245.01518189499933, 245.01518189499933,
               0.0, 0.0, 0.0),
        region="apex", its_max=1, degenerate=True),
    "near-apex FACE (the header's CHECK_APEX_REGION says apex)": dict(
        eps=(1.03e-05, 7.9999999999999996e-07, 5.5000000000000003e-07,
             0.0, 0.0, 0.0),
        sigma=(244.88118682640001, 131.8202302997, 122.39973437259999,
               0.0, 0.0, 0.0),
        region="face", its_max=5, degenerate=False),
}


@pytest.mark.parametrize("name", list(HB_REGION_CASES))
def test_gate1_region_returns_match_the_p0_oracle(cp_hb_available, name):
    case = HB_REGION_CASES[name]
    r = drive(lambda t: mat_hb(t), [case["eps"]], want=("cp_iterations",))
    assert r["codes"] == [0], "the one-step region deck did not converge: %s" % (
        r["codes"],)
    sig = r["sigma"][0]
    rel = _rel(sig, case["sigma"])
    its = int(r["cp_iterations"][0][0])
    print("gate 1a [%s]\n"
          "        region       = %s\n"
          "        principals   = %s\n"
          "        rel. error   = %.3e   (pin %.0e)\n"
          "        |f| committed= %.3e   cp_iterations = %d"
          % (name, case["region"], np.array2string(princ_desc(sig), precision=8),
             rel, RTOL_REGION, abs(f_hb(sig)), its))
    assert rel <= RTOL_REGION, rel
    # the composite must be satisfied at the committed state; the gradient
    # scaled floor is the HB instance of ADR-94's f_relative_tol lesson
    assert f_hb(sig) <= 1e-8 * HB_SCALE, f_hb(sig)
    assert 1 <= its <= case["its_max"], its


def test_gate1_curved_face_newton_stays_within_five_iterations(cp_hb_available):
    """The <= 5 gate is the reason the map is written in the surface's OWN
    variable with a NORMALIZED flow direction: the P0 oracle measured 6 without
    the normalization and 6 without the cutting-plane seed."""
    counts = {}
    for name, case in HB_REGION_CASES.items():
        if case["region"] != "face":
            continue
        r = drive(lambda t: mat_hb(t), [case["eps"]], want=("cp_iterations",))
        assert r["codes"] == [0]
        counts[name] = int(r["cp_iterations"][0][0])
    print("gate 1a curved-face Newton counts: %s" % counts)
    assert counts, "no face case ran"
    assert max(counts.values()) <= 5, counts


# ---------------------------------------------------------------------------
# GATE 1b -- admissible at EVERY commit, on three multi-step paths
# ---------------------------------------------------------------------------
EPS_TRIAX = (0.00023, 3.0e-05, -0.00092, 0.0, 0.0, 0.0)
EPS_SHEAR = (0.0, 0.0, 0.0, 0.0011, 0.0, 0.0)
# two legs whose principal DIRECTIONS rotate between them, in the P2 shape: the
# second leg may only switch on a component the first left at zero, or `drive`'s
# component-disjointness check (one Path series per leg) cannot be satisfied.
EPS_ROT_A = (0.0, 0.0, -0.0015, 0.001, 0.0, 0.0)
EPS_ROT_B = (0.0, 0.0, -0.0015, 0.001, 0.001, 0.0)

PATHS = {
    "triaxial compression": ([EPS_TRIAX], 10),
    "simple shear": ([EPS_SHEAR], 10),
    "rotating principal directions": ([EPS_ROT_A, EPS_ROT_B], 10),
}


@pytest.mark.parametrize("name", list(PATHS))
def test_gate1_admissible_at_every_commit(cp_hb_available, name):
    legs, nstep = PATHS[name]
    r = drive(lambda t: mat_hb(t), legs, nstep=nstep, want=("cp_iterations",))
    assert all(c == 0 for c in r["codes"]), (name, r["codes"])
    fs = [f_hb(s) for s in r["sigma"]]
    its = sorted({int(v[0]) for v in r["cp_iterations"]})
    plastic = sum(1 for v in r["cp_iterations"] if int(v[0]) > 0)
    print("gate 1b [%s] steps=%d plastic=%d  worst f = %+.3e  "
          "(1e-8*scale = %.2e)  cp_iterations seen = %s"
          % (name, len(fs), plastic, max(fs), 1e-8 * HB_SCALE, its))
    assert max(fs) <= 1e-8 * HB_SCALE, max(fs)
    assert plastic > 0, "path never yielded; the gate is vacuous"
    assert max(its) <= 5, its


def test_gate1_rotating_path_principal_directions_really_rotate(cp_hb_available):
    legs, nstep = PATHS["rotating principal directions"]
    r = drive(lambda t: mat_hb(t), legs, nstep=nstep)
    assert all(c == 0 for c in r["codes"])
    v0 = np.linalg.eigh(_tensor(r["sigma"][nstep - 1]))[1][:, -1]
    v1 = np.linalg.eigh(_tensor(r["sigma"][-1]))[1][:, -1]
    ang = math.degrees(math.acos(min(1.0, abs(float(v0 @ v1)))))
    print("gate 1b rotating path: first eigenvector moved %.2f deg" % ang)
    assert ang > 5.0, ang


# ---------------------------------------------------------------------------
# GATE 1c -- the ADR-94 tension plateau, under Closest_Point
#
# Uniaxial STRESS tension (lateral faces free), so the stress path is exactly
# sigma = (sxx, 0, 0) and the limit is the root of
#     sxx = sigma_ci * (s - mb*sxx/sigma_ci)^a
# -- 244.4854419069 kPa, 0.216 % BELOW the apex T = 245.0151818950 (at y1 = T
# the clamp leaves f_shear = y1 > 0, so the surface is reached strictly before
# the vertex).  ADR-94's own test pins BE's plateau at |sigma_t| = T with
# rel = 5e-2, which does not separate the two; this one pins the exact root.
# ---------------------------------------------------------------------------
_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
_FIX = {1: (1, 1, 1), 2: (0, 1, 1), 3: (0, 0, 1), 4: (1, 0, 1),
        5: (1, 1, 0), 6: (0, 1, 0), 8: (1, 0, 0)}      # node 7 carries no fix
_XFACE = (2, 3, 6, 7)


def _uniaxial_stress_root():
    lo, hi = 0.0, HB_T
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        arg = max(HB_S - HB_MB * mid / HB_SIGCI, 0.0)
        if mid - HB_SIGCI * arg ** HB_A < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


HB_UNIAXIAL_TENSION_LIMIT = _uniaxial_stress_root()   # 244.4854419069


def _uniaxial_tension(mat_fn, nsteps=60, factor=3.0):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *map(float, c))
    for t, m in _FIX.items():                 # ADR-94's own restraint pattern
        ops.fix(t, *m)
    mat_fn(1)
    ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in _XFACE:
        ops.sp(n, 1, factor * HB_T / HB_E)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-10, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    codes, hist = [], []
    for _ in range(nsteps):
        rc = ops.analyze(1)
        codes.append(rc)
        if rc != 0:
            break
        ops.eleResponse(1, "forces")
        hist.append(list(ops.eleResponse(1, "stresses"))[0:6])
    return codes, (np.array(hist) if hist else np.zeros((0, 6)))


def test_gate1_tension_plateau_under_closest_point(cp_hb_available):
    codes, hist = _uniaxial_tension(lambda t: mat_hb(t))
    assert len(hist) > 0, "no step committed at all"
    sxx = float(hist[-1, 0])
    ok = sum(1 for c in codes if c == 0)
    print("gate 1c uniaxial tension under Closest_Point:\n"
          "        steps committed  = %d / %d   (codes tail %s)\n"
          "        last sigma_xx    = %.7f kPa\n"
          "        exact limit      = %.7f kPa   (rel %.3e)\n"
          "        apex T           = %.7f kPa   (ADR-94's own pin)"
          % (ok, 60, codes[-3:], sxx, HB_UNIAXIAL_TENSION_LIMIT,
             abs(sxx - HB_UNIAXIAL_TENSION_LIMIT) / HB_UNIAXIAL_TENSION_LIMIT,
             HB_T))
    # ADR-94's own tolerance on the plateau, kept so this test and that one
    # agree on the physics; the tighter statement is the print above.
    assert sxx == pytest.approx(HB_T, rel=5.0e-2), sxx
    assert sxx <= HB_T * (1.0 + 1e-9), (
        "the plateau is ABOVE the apex, i.e. outside the surface: %r" % sxx)
    assert sxx == pytest.approx(HB_UNIAXIAL_TENSION_LIMIT, rel=1.0e-3), (
        "the Closest_Point plateau is not the exact uniaxial-stress limit")


def test_gate1_tension_past_the_corner_is_finite(cp_hb_available):
    """ADR-94's residual: the corner is where Backward_Euler stalls.  The
    frame-consistent potential HAS a return to the apex (the shipped Tresca `g`
    does not -- all six of its flow directions have negative trace there), so a
    trial pushed WELL past the corner must commit a finite, admissible state."""
    for label, eps in (("2x past the corner",
                        tuple(2.0 * HB_T / HB_E for _ in range(3)) + (0., 0., 0.)),
                       ("20x past the corner",
                        tuple(20.0 * HB_T / HB_E for _ in range(3)) + (0., 0., 0.)),
                       ("hydrostatic + tiny deviator",
                        (10.0 * HB_T / HB_E, 10.0 * HB_T / HB_E + 1e-15,
                         10.0 * HB_T / HB_E - 1e-15, 0., 0., 0.))):
        r = drive(lambda t: mat_hb(t), [eps])
        assert r["codes"] == [0], (label, r["codes"])
        sig = r["sigma"][0]
        assert np.all(np.isfinite(sig)), (label, sig)
        print("gate 1c [%s] committed principals = %s  f = %+.3e"
              % (label, np.array2string(princ_desc(sig), precision=8), f_hb(sig)))
        assert _rel(sig, np.array([HB_T, HB_T, HB_T, 0., 0., 0.])) <= 1e-9, sig


# ===========================================================================
# GATE 2 -- the consistent tangent against a central difference of the
# BINARY's own assembled internal force.  No numpy reference anywhere.
#
# Two rigs, as in P2: `fd_tangent_driver`'s load-driven oedometric rig reaches
# a DEGENERATE-eigenvalue state (y1 == y2, the l'Hopital branch of the rotation
# term), and P2's free-node rig -- seven of eight nodes prescribed, node 7 free
# -- is the only one that reaches a genuine FACE state with three separated
# principal stresses.  The APEX tangent is rank 0 by construction (the vertex
# does not move), so there is NO finite-difference gate on it and none is
# missing; that statement is the same one P1 makes for the Drucker-Prager apex
# and P2 for the Mohr-Coulomb vertex.
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
    _fd_build(mat_fn, eps, ele, free_node)
    rc = ops.analyze(1)
    if rc != 0:
        return dict(rc=rc, rel_err=float("nan"))
    eqn = [ops.nodeDOFs(free_node)[d] for d in range(3)]
    d = ops.printA("-sparse", "-ret")     # NOT printA() alone
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


FD_FACE = np.array(HB_REGION_CASES["face (no shear)"]["eps"])
FD_FACE_SHEARED = np.array(
    HB_REGION_CASES["face (sheared, exercises the eigenprojection rotation "
                    "term)"]["eps"])


@pytest.mark.parametrize("name,eps", [("face, axis aligned", FD_FACE),
                                      ("face, sheared", FD_FACE_SHEARED)])
def test_gate2_algorithmic_fd_on_the_face_region(cp_hb_available, name, eps):
    """Three well-separated principal stresses, so BOTH the curvature term
    dl*dm/dy inside the Jacobian AND the eigenprojection rotation term are
    live.  This is the measurement the P3 mutation gate kills."""
    r = fd_free_node(lambda t: mat_hb(t), eps)
    assert r["rc"] == 0, "the free-node rig did not converge"
    x = princ_desc(r["sig"])
    sep = min(abs(x[0] - x[1]), abs(x[1] - x[2])) / max(1.0, abs(x[2]))
    print("gate 2 (%s) principals = %s\n"
          "        separation = %.3e   f = %+.3e   rel_err = %.3e"
          % (name, np.array2string(x, precision=6), sep, f_hb(r["sig"]),
             r["rel_err"]))
    assert f_hb(r["sig"]) >= -1e-6 * HB_SCALE, "state is elastic; gate vacuous"
    assert sep > 1e-3, "principal stresses are not separated; not the face case"
    assert r["rel_err"] <= RTOL, r["rel_err"]


def test_gate2_algorithmic_fd_on_the_load_driven_rig(cp_hb_available):
    """`fd_tangent_driver`'s own oedometric rig.  It lands on the DEGENERATE
    y1 == y2 edge, i.e. the l'Hopital branch of the rotation term -- the region
    the free-node rig above cannot produce."""
    r = FD.fd_check(lambda t: mat_hb(t), rig="uniaxial",
                    load=(0.0, 0.0, -3.0e4),
                    label="HB edge region, Closest_Point / Algorithmic")
    assert r["rc"] == 0, "the HB load-driven oedometric rig did not converge"
    x = princ_desc(r["sig"])
    print("gate 2 (load-driven rig) principals = %s  f = %+.3e\n"
          "        rel_err = %.3e   rel_fro = %.3e   degenerate = %s"
          % (np.array2string(x, precision=6), f_hb(r["sig"]), r["rel_err"],
             r["rel_fro"], abs(x[0] - x[1]) <= 1e-8 * max(1.0, abs(x[0]))))
    assert f_hb(r["sig"]) >= -1e-6 * HB_SCALE, "state is elastic; gate vacuous"
    assert r["rel_err"] <= RTOL, r["rel_err"]


def _oedometric_iterations(mat_fn, nsteps=6, load=-3.0e4):
    """One confined cube, load driven, counting the GLOBAL Newton iterations per
    step.  Confinement grows with compression, so perfect plasticity has no
    limit point here."""
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


def test_gate2_iteration_contrast_against_backward_euler(cp_hb_available):
    """On this deck the contrast is not a ratio, it is an OUTCOME: with the
    shipped tangents `Backward_Euler` does not converge at all where
    `Closest_Point`/`Algorithmic` does.  (ADR-94 M3 measured `Continuum` and
    `Secant` at 57 % and 80 % against a central difference of the material's own
    response; on a perfectly plastic Hoek-Brown edge state, where the consistent
    tangent is rank deficient anyway, a Newton at that error simply stops.)  The
    assertion is written so that it still holds -- as a >= 2x iteration
    reduction -- if a future change makes Backward_Euler converge here."""
    it_cp, sig_cp, rc_cp = _oedometric_iterations(lambda t: mat_hb(t))
    it_co, sig_co, rc_co = _oedometric_iterations(
        lambda t: mat_hb(t, method="Backward_Euler", tangent="Continuum"))
    it_se, sig_se, rc_se = _oedometric_iterations(
        lambda t: mat_hb(t, method="Backward_Euler", tangent="Secant"))
    assert rc_cp == 0, "Closest_Point HB oedometric rig did not converge (%s)" % rc_cp
    print("gate 2 iteration contrast (HB oedometric, 6 steps, %.1e/node)\n"
          "        Closest_Point / Algorithmic : %s  total %s   (analyze -> 0)\n"
          "        Backward_Euler / Continuum  : %s   (analyze -> %s)\n"
          "        Backward_Euler / Secant     : %s   (analyze -> %s)\n"
          "        committed sigma CP = %s  (f = %+.3e)"
          % (-3.0e4, it_cp, sum(it_cp),
             (it_co if rc_co == 0 else "did not converge"), rc_co,
             (it_se if rc_se == 0 else "did not converge"), rc_se,
             np.array2string(sig_cp, precision=4), f_hb(sig_cp)))
    assert f_hb(sig_cp) >= -1e-6 * HB_SCALE, "never yielded; contrast is vacuous"
    assert f_hb(sig_cp) <= 1e-8 * HB_SCALE
    for lbl, rc, it in (("Continuum", rc_co, it_co), ("Secant", rc_se, it_se)):
        if rc == 0:
            assert 2 * sum(it_cp) <= sum(it), (lbl, it_cp, it)


# ===========================================================================
# GATE 3 -- hardening.  There is none to gate: `HoekBrown_YF` and
# `HoekBrown_PF` are registered with `BackStress<NullHardeningTensorFunction>`
# and nothing else, and `HoekBrown_YF::YIELD_FUNCTION_HARDENING` returns 0.0.
# The principal-space map has no q-row, and `asdp_all_ivs_are_inert` is folded
# into `ladruno_cp_principal_family` at COMPILE time, so a hypothetical
# HoekBrown_YF<ArmstrongFrederick..> would be refused at parse time rather than
# silently integrated with its hardening ignored.  The test below is that
# statement, made executable.
# ===========================================================================
def test_gate3_the_hoek_brown_family_is_perfectly_plastic(cp_hb_available):
    r = drive(lambda t: mat_hb(t), [EPS_TRIAX], nstep=10)
    assert all(c == 0 for c in r["codes"])
    bs = np.array(list(ops.eleResponse(1, "material", 1, "BackStress")))
    print("gate 3: committed BackStress after 10 plastic steps = %s"
          % np.array2string(bs, precision=3))
    assert np.max(np.abs(bs)) == 0.0, (
        "the Hoek-Brown internal variable moved; it is registered Null-hardening"
        " and the principal-space map has no q-row")


# ===========================================================================
# GATE 4 -- Backward_Euler is untouched, and the potential gap is PINNED
# ===========================================================================
def test_gate4_cp_and_be_disagree_by_the_measured_potential_gap(cp_hb_available):
    """CP and BE are NOT expected to agree on a Hoek-Brown deck, and this test
    says by how much and WHY: `HoekBrown_PF::g` is a Tresca potential (zero
    dilatancy, HB_mb_psi inert) while `Closest_Point` uses the frame-consistent
    Hoek-Brown potential.  The P0 oracle measured the same trial at plastic
    volumetric strain +2.208e-05 (frame consistent) against -6.8e-21 (header).

    This test turns RED the day the shipped `g` is fixed -- which is the point:
    that fix is the owner's separate PR, and this is its warrant."""
    eps = HB_REGION_CASES["face (no shear)"]["eps"]
    r_cp = drive(lambda t: mat_hb(t), [eps], want=("pstrain",))
    r_be = drive(lambda t: mat_hb(t, method="Backward_Euler",
                                  tangent="Continuum"), [eps], want=("pstrain",))
    assert r_cp["codes"] == [0] and r_be["codes"] == [0], (r_cp["codes"],
                                                           r_be["codes"])
    s_cp, s_be = r_cp["sigma"][0], r_be["sigma"][0]
    ev_cp = float(np.sum(r_cp["pstrain"][0][:3]))
    ev_be = float(np.sum(r_be["pstrain"][0][:3]))
    gap = _rel(s_cp, s_be)
    print("gate 4 potential gap on the oracle's face trial:\n"
          "        CP principals = %s   plastic eps_vol = %+.6e\n"
          "        BE principals = %s   plastic eps_vol = %+.6e\n"
          "        |CP - BE| / |BE| = %.3e   (%.2f %% of the strength scale)"
          % (np.array2string(princ_desc(s_cp), precision=6), ev_cp,
             np.array2string(princ_desc(s_be), precision=6), ev_be,
             gap, 100.0 * np.max(np.abs(s_cp - s_be)) / HB_SCALE))
    # the shipped Tresca `g` is exactly non-dilatant ...
    assert abs(ev_be) <= 1e-12, (
        "Backward_Euler's flow is no longer non-dilatant -- HoekBrown_PF::g may"
        " have been fixed; see ADR-97 P3 and re-derive this pin (%r)" % ev_be)
    # ... and the frame-consistent one is not
    assert ev_cp >= 1e-5, ev_cp
    assert gap >= 1e-2, (
        "CP and BE now agree on a Hoek-Brown deck; either the shipped potential"
        " was fixed or Closest_Point silently fell back to it (%r)" % gap)
    # both must still be admissible
    assert f_hb(s_cp) <= 1e-8 * HB_SCALE and f_hb(s_be) <= 1e-6 * HB_SCALE


def test_gate4_step_refinement_on_a_curved_surface(cp_hb_available):
    """A closest point onto a PLANE reached along a proportional path is step
    independent (P2 measured 1.5e-16 to 5.2e-16 over N = 1..40).  A CURVED
    surface cannot be: after step k the state sits ON the surface, so step
    k+1's trial starts from a curved point and its closest point is a different
    point.  That is a property of the map, not a defect -- and it is worth
    pinning, because a silent fall-back to an incremental cutting plane would
    look exactly like a LARGE version of it.

    N = 1 is the oracle's own one-step return and must be exact; the refinement
    sequence must then CONVERGE (successive differences shrinking) rather than
    wander."""
    ref = HB_REGION_CASES["face (no shear)"]["sigma"]
    eps = HB_REGION_CASES["face (no shear)"]["eps"]
    sig = {}
    for nstep in (1, 4, 10, 40, 160):
        r = drive(lambda t: mat_hb(t), [eps], nstep=nstep)
        assert all(c == 0 for c in r["codes"]), (nstep, r["codes"])
        sig[nstep] = r["sigma"][-1]
    print("gate 4 step refinement on the curved face (oracle = the N = 1 map):")
    for nstep in (1, 4, 10, 40, 160):
        print("        N = %3d   rel vs the oracle = %.3e   rel vs N = 160 = %.3e"
              % (nstep, _rel(sig[nstep], ref), _rel(sig[nstep], sig[160])))
    assert _rel(sig[1], ref) <= 1e-9, _rel(sig[1], ref)
    d = [_rel(sig[n], sig[160]) for n in (1, 4, 10, 40)]
    assert d[0] > d[1] > d[2] > d[3], (
        "the step-refinement sequence is not converging: %s" % d)
    assert d[0] <= 1e-2, ("the one-step and incremental limits differ by more"
                          " than 1 %% -- that is not curvature, that is a"
                          " different map: %r" % d[0])


# ===========================================================================
# GATE 6 -- fail loud
# ===========================================================================
HB_YF_P = ["HB_sigci", HB_SIGCI, "HB_mb", HB_MB, "HB_s", HB_S, "HB_a", HB_A,
           "HB_ds", 0.0]
HB_PF_P = ["HB_sigci", HB_SIGCI, "HB_mb_psi", HB_MB, "HB_s", HB_S, "HB_a", HB_A,
           "HB_ds", 0.0]
EL_P = ["YoungsModulus", HB_E, "PoissonsRatio", HB_NU, "MassDensity", 0.0]
BS0 = ["BackStress", 0., 0., 0., 0., 0., 0.]
# `utuple_concat_unique_type` de-duplicates internal variables by TYPE, so a
# pairing whose YF and PF carry BackStress with DIFFERENT hardening laws gets
# TWO back stresses (ADR-97 P1 finding).  The IV strings below spell that out;
# a wrong one makes the material unfindable and the refusal test vacuous, which
# is why every row is checked under Backward_Euler as well.
IV_LIN = "BackStress(TensorLinearHardeningFunction):"
IV_DPC = "DP_cohesion(ScalarLinearHardeningFunction):"
IV_YS = "YieldStress(ScalarLinearHardeningFunction):"
LIN_P = ["TensorLinearHardeningParameter", 0.0,
         "ScalarLinearHardeningParameter", 0.0]


def _mat_pair(tag, yf, pf, iv, params, ivs, method="Closest_Point",
              tangent="Algorithmic"):
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag, yf, pf, "LinearIsotropic3D_EL", iv,
        "Begin_Model_Parameters", *params, "End_Model_Parameters",
        "Begin_Internal_Variables", *ivs, "End_Internal_Variables",
        "Begin_Integration_Options",
        "integration_method", method, "tangent_type", tangent,
        "n_max_iterations", 100,
        "End_Integration_Options",
    )


# every MIXED pairing the generator registers that involves Hoek-Brown.  Each
# row carries its OWN parameter and internal-variable list, because ADR-94 made
# a missing model parameter an error: a refusal test built on a wrong parameter
# set is VACUOUS (it refuses for the wrong reason), which is why every row below
# is checked in BOTH directions -- it must CONSTRUCT under Backward_Euler and be
# REFUSED under Closest_Point.
REFUSED_HB_MIXED = [
    ("MohrCoulomb_YF x HoekBrown_PF", "MohrCoulomb_YF", "HoekBrown_PF",
     IV_NULL, ["MC_phi", 30.0, "MC_c", 100.0, "MC_ds", 0.0] + HB_PF_P, []),
    ("HoekBrown_YF x VonMises_PF (Null)", "HoekBrown_YF", "VonMises_PF",
     IV_NULL, HB_YF_P, []),
    ("HoekBrown_YF x VonMises_PF (TensorLinear)", "HoekBrown_YF", "VonMises_PF",
     IV_NULL + IV_LIN, HB_YF_P + ["TensorLinearHardeningParameter", 0.0], BS0),
    ("HoekBrown_YF x MohrCoulomb_PF", "HoekBrown_YF", "MohrCoulomb_PF",
     IV_NULL, HB_YF_P + ["MC_phi", 30.0, "MC_c", 100.0, "MC_ds", 0.0,
                         "MC_psi", 10.0], []),
    ("HoekBrown_YF x DruckerPrager_PF", "HoekBrown_YF", "DruckerPrager_PF",
     IV_NULL + IV_LIN + IV_DPC, HB_YF_P + LIN_P + ["DP_etabar", 0.4],
     BS0 + ["DP_cohesion", 0.0]),
    ("VonMises_YF x HoekBrown_PF", "VonMises_YF", "HoekBrown_PF",
     IV_LIN + IV_YS + IV_NULL, HB_PF_P + LIN_P,
     ["YieldStress", 3000.0] + BS0),
    ("DruckerPrager_YF x HoekBrown_PF", "DruckerPrager_YF", "HoekBrown_PF",
     IV_LIN + IV_DPC + IV_NULL,
     HB_PF_P + LIN_P + ["DP_xi_c", 2000.0, "DP_eta", 0.4],
     ["DP_cohesion", 0.0] + BS0),
]


def _try_build(yf, pf, iv, params, ivs, method):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for k, (x, y, z) in enumerate(NODES):
        ops.node(k + 1, float(x), float(y), float(z))
    try:
        _mat_pair(1, yf, pf, iv, EL_P + list(params), BS0 + list(ivs),
                  method=method,
                  tangent=("Algorithmic" if method == "Closest_Point"
                           else "Continuum"))
        ops.element("LadrunoBrick", 1, *range(1, 9), 1)
        ok = True
    except Exception:
        ok = False
    ops.wipe()
    return ok


@pytest.mark.parametrize("label,yf,pf,iv,params,ivs", REFUSED_HB_MIXED,
                         ids=[c[0] for c in REFUSED_HB_MIXED])
def test_gate6_mixed_hoek_brown_pairings_stay_refused(cp_hb_available, label,
                                                      yf, pf, iv, params, ivs):
    """The family-marker EQUALITY rule P2 introduced, applied to Hoek-Brown: a
    mixed pairing has no verified map at all -- the principal return assumes
    BOTH the surface and the potential are of that family -- so it is refused at
    PARSE time rather than silently approximated.  Checked in both directions so
    the refusal is attributable to `Closest_Point` and not to a typo in the
    parameter list."""
    be = _try_build(yf, pf, iv, params, ivs, "Backward_Euler")
    cp = _try_build(yf, pf, iv, params, ivs, "Closest_Point")
    print("gate 6: %-42s  Backward_Euler -> %s   Closest_Point -> %s"
          % (label, "built" if be else "FAILED",
             "ACCEPTED" if cp else "refused"))
    assert be, ("%s does not even construct under Backward_Euler, so the"
                " refusal below would be vacuous -- fix the parameter list"
                % label)
    assert not cp, "%s must stay refused under Closest_Point (ADR-97 P3)" % label


def test_gate6_matched_hoek_brown_pairing_is_accepted(cp_hb_available):
    """The positive half of the same matrix -- without it the tests above pass
    for a build in which Closest_Point is refused for everything."""
    assert _try_build("HoekBrown_YF", "HoekBrown_PF", IV_NULL,
                      HB_YF_P + ["HB_mb_psi", HB_MB], [], "Closest_Point"),         "HoekBrown_YF x HoekBrown_PF must be ACCEPTED by ADR-97 P3"
    assert _constructible(lambda t: mat_hb(t))


def test_gate6_algorithmic_still_refused_with_backward_euler(cp_hb_available):
    """ADR-97 D2 survives P3: a consistent tangent is defined only relative to
    a specific committed map."""
    assert not _constructible(
        lambda t: mat_hb(t, method="Backward_Euler", tangent="Algorithmic"))


def test_gate6_strict_convergence_is_byte_inert_on_a_converging_deck(
        cp_hb_available):
    eps = HB_REGION_CASES["face (no shear)"]["eps"]
    a = drive(lambda t: mat_hb(t, strict=0), [eps], nstep=4)
    b = drive(lambda t: mat_hb(t, strict=1), [eps], nstep=4)
    assert a["codes"] == b["codes"] == [0, 0, 0, 0]
    gap = float(np.max(np.abs(a["sigma"] - b["sigma"])))
    print("gate 6 strict_convergence byte-inertness: max |delta sigma| = %r" % gap)
    assert gap == 0.0, gap


def test_gate6_degenerate_and_hydrostatic_states_commit_no_nan(cp_hb_available):
    """The class of state ADR-94 B4 used to commit NaN through: exactly
    hydrostatic, and hydrostatic plus a deviator at the round-off scale."""
    base = 5.0 * HB_T / HB_E
    for label, eps in (("exactly hydrostatic tension", (base, base, base, 0., 0., 0.)),
                       ("hydrostatic + 1e-15", (base, base + 1e-15, base - 1e-15,
                                                0., 0., 0.)),
                       ("hydrostatic + 1e-9", (base, base + 1e-9, base - 1e-9,
                                               0., 0., 0.)),
                       ("hydrostatic COMPRESSION", (-base, -base, -base, 0., 0., 0.))):
        r = drive(lambda t: mat_hb(t), [eps])
        assert r["codes"] == [0], (label, r["codes"])
        sig = r["sigma"][0]
        assert np.all(np.isfinite(sig)), (label, sig)
        print("gate 6 [%s] -> principals %s   f = %+.3e"
              % (label, np.array2string(princ_desc(sig), precision=8), f_hb(sig)))
        assert f_hb(sig) <= 1e-8 * HB_SCALE, (label, f_hb(sig))
