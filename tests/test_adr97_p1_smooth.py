"""ADR-97 P1 (wp/97b) — gates 1, 2 and 3 for `integration_method Closest_Point`
and `tangent_type Algorithmic` on the SMOOTH families (VonMises,
Drucker-Prager including the apex).

WHAT IS UNDER TEST
------------------
`Closest_Point` solves the fully implicit closest-point return map as ONE
coupled Newton in ``x = (sigma_{n+1}, q_{n+1}, dlambda)``; `Algorithmic` is the
exact consistent tangent of THAT map, taken from the converged Jacobian.  The
shipped `Backward_Euler` is an Ortiz-Simo CUTTING PLANE whose fixed point is
``sigma_tr - sum_k dl_k E m(sigma^k)`` and whose internal-variable update is a
Newton-path quadrature exact only for constant ``h`` (ADR-94 M3).

GATE 1 — material-level correctness.  Committed stress against the P0 numpy
oracles (`Ladruno_implementation/adr97_oracle/`, transcript
`reference_output.txt`, build 3622d6214), ``|f| <= tol`` at every commit, and
the local Newton converging in <= 5 iterations.

GATE 2 — the tangent.  `Algorithmic` against a central difference of the
BINARY'S OWN assembled residual (`adr97_oracle/fd_tangent_driver.py::fd_check`),
with the measured `Backward_Euler`/`Continuum` value 0.573447 pinned as the
negative control; plus the ADR-94 two-cube model's Newton iteration count.

GATE 3 — path error.  The cutting plane's per-step path dependence is only
visible when ``h`` is not constant along the iterates, i.e. for
Armstrong-Frederick (22 of the 46 registered specializations carry it) — which
is why ADR-94 H6 could not see it with linear hardening.  Measured as the
step-refinement error of each map against its own fine-step limit.

RIG
---
One unit-cube `LadrunoBrick` with EVERY displacement degree of freedom
prescribed to a homogeneous strain field, so the Gauss-point strain IS the
prescribed Voigt strain and the analysis is a pure state determination (the
`adr97_oracle/fd_tangent_driver.py` pattern).  `system UmfPack`: `FullGeneral`
crashes a fully prescribed rig (`FullGenLinSOE` N = 0) and `C_alg` is
UNSYMMETRIC whenever ``m != n``, so `ProfileSPD` would be wrong too.
`stdBrick` appears ONLY inside `fd_check` (gate 2), where the oracle's negative
control was measured and where no material return code is under test.

The displacement field is the ASYMMETRIC one
``u = (e11 x + g12 y + g13 z, e22 y + g23 z, e33 z)``: its symmetric gradient is
exactly the engineering Voigt strain, and each Voigt slot feeds exactly ONE
displacement component, which is what lets a multi-leg (rotating-normal) path be
driven by one `sp` per DOF (two `sp` on one DOF is the rank-deficient-KKT trap).

Zone-A, ~20 s.
"""
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

# Cross-platform float pins are >= 1e-6 relative by testbed policy; the numbers
# actually measured on Windows/MSVC are printed by each test and recorded in
# Ladruno_implementation/reviews/adr97_p1_report.md.
RTOL = 1e-6

W_STRESS = np.array([1., 1., 1., 2., 2., 2.])
SQRT_2_over_3 = 0.816496580928          # the TRUNCATED literal, Globals.h:41

NODES = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
         (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]

# ---------------------------------------------------------------------------
# materials (the oracles' own constants)
# ---------------------------------------------------------------------------
E_VM, NU_VM, SY_VM = 70000.0, 0.3, 30.0
IV_VM_LIN = ("BackStress(TensorLinearHardeningFunction):"
             "YieldStress(ScalarLinearHardeningFunction):")
IV_VM_AF = ("BackStress(ArmstrongFrederickHardeningFunction):"
            "YieldStress(ScalarLinearHardeningFunction):")

DP_E, DP_NU, DP_XI_C, DP_ETA = 30000.0, 0.25, 20.0, 0.4
IV_DP = ("BackStress(TensorLinearHardeningFunction):"
         "DP_cohesion(ScalarLinearHardeningFunction):")


def mat_vm(tag, tangent="Algorithmic", method="Closest_Point", hiso=0.0,
           hkin=0.0, niter=100, strict=None):
    extra = ["strict_convergence", int(strict)] if strict is not None else []
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "VonMises_YF", "VonMises_PF", "LinearIsotropic3D_EL", IV_VM_LIN,
        "Begin_Model_Parameters",
        "YoungsModulus", E_VM, "PoissonsRatio", NU_VM,
        "ScalarLinearHardeningParameter", float(hiso),
        "TensorLinearHardeningParameter", float(hkin), "MassDensity", 0.0,
        "End_Model_Parameters",
        "Begin_Internal_Variables", "YieldStress", SY_VM,
        "BackStress", 0., 0., 0., 0., 0., 0., "End_Internal_Variables",
        "Begin_Integration_Options",
        "integration_method", method, "tangent_type", tangent,
        "n_max_iterations", int(niter), *extra,
        "End_Integration_Options",
    )


def mat_vm_af(tag, tangent="Algorithmic", method="Closest_Point",
              ha=15000.0, cr=300.0, niter=100):
    """Both the yield function's and the flow direction's back stress are the
    SAME Armstrong-Frederick internal variable (the IV storage de-duplicates by
    TYPE, so a YF/PF pair carrying different hardening laws would give two
    independent back stresses -- see the ADR-97 report)."""
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "VonMises_YF", "VonMises_PF", "LinearIsotropic3D_EL", IV_VM_AF,
        "Begin_Model_Parameters",
        "YoungsModulus", E_VM, "PoissonsRatio", NU_VM,
        "ScalarLinearHardeningParameter", 0.0,
        "AF_ha", float(ha), "AF_cr", float(cr), "MassDensity", 0.0,
        "End_Model_Parameters",
        "Begin_Internal_Variables", "YieldStress", SY_VM,
        "BackStress", 0., 0., 0., 0., 0., 0., "End_Internal_Variables",
        "Begin_Integration_Options",
        "integration_method", method, "tangent_type", tangent,
        "n_max_iterations", int(niter),
        "End_Integration_Options",
    )


def mat_dp(tag, etabar=0.4, hiso=0.0, tangent="Algorithmic",
           method="Closest_Point", niter=100):
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "DruckerPrager_YF", "DruckerPrager_PF", "LinearIsotropic3D_EL", IV_DP,
        "Begin_Model_Parameters",
        "YoungsModulus", DP_E, "PoissonsRatio", DP_NU,
        "DP_xi_c", DP_XI_C, "DP_eta", DP_ETA, "DP_etabar", float(etabar),
        "TensorLinearHardeningParameter", 0.0,
        "ScalarLinearHardeningParameter", float(hiso),
        "MassDensity", 0.0,
        "End_Model_Parameters",
        "Begin_Internal_Variables",
        "BackStress", 0., 0., 0., 0., 0., 0.,
        "DP_cohesion", 0.0,
        "End_Internal_Variables",
        "Begin_Integration_Options",
        "integration_method", method, "tangent_type", tangent,
        "n_max_iterations", int(niter),
        "End_Integration_Options",
    )


# ---------------------------------------------------------------------------
# the prescribed-strain rig
# ---------------------------------------------------------------------------
def _u_of(eps, x, y, z):
    e11, e22, e33, g12, g23, g13 = eps
    return (e11 * x + g12 * y + g13 * z, e22 * y + g23 * z, e33 * z)


def drive(mat_fn, legs, nstep=10, ele="LadrunoBrick", tol=1e-13, maxiter=60,
          want=("stresses",)):
    """Prescribe a multi-leg homogeneous strain path and return the history."""
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
            "the path is not component-disjoint at DOF %s (legs %s change it): "
            "two sp constraints on one DOF is the rank-deficient-KKT trap"
            % (key, changed))
        owner[key] = changed[0] if changed else 0

    for L in range(nl):
        if nl == 1:
            ops.timeSeries("Linear", 100 + L)
        else:
            # One time point PAST the end of the path: a Path series returns 0
            # outside its defined range, and the last step lands exactly on the
            # final time -- without the overshoot the whole path unloads to zero
            # in one step and the material silently reports a plausible-looking
            # (on-surface) but wrong stress.
            times = [float(t) for t in range(nl + 2)]
            values = [0.0] * (L + 1) + [1.0] * (nl - L + 1)
            ops.timeSeries("Path", 100 + L, "-time", *times, "-values", *values)
        ops.pattern("Plain", 100 + L, 100 + L)
        for key, o in owner.items():
            if o == L:
                ops.sp(key[0], key[1], float(U[L][key]))

    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("UmfPack")            # NOT FullGeneral (crashes a fully sp rig),
    ops.test("NormDispIncr", tol, maxiter, 0)   # NOT ProfileSPD (C_alg unsym.)
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
            if w == "stresses":
                continue
            v = list(ops.eleResponse(1, "material", 1, w))
            out[w].append(np.array(v) if v else np.array([]))
    out["sigma"] = np.array(out["sigma"])
    out["eps"] = np.array(out["eps"])
    return out


def _dev(v):
    p = (v[0] + v[1] + v[2]) / 3.0
    d = np.array(v, dtype=float)
    d[:3] -= p
    return d


def _f_vm(sig, alpha, k):
    r = _dev(sig) - alpha
    return float(np.sqrt((r * W_STRESS * r).sum()) - SQRT_2_over_3 * k)


def _f_dp(sig, alpha, xi_c=DP_XI_C, eta=DP_ETA):
    r = _dev(sig) - alpha
    q = float(np.sqrt(0.5 * (r * W_STRESS * r).sum()))
    return q + eta * (sig[0] + sig[1] + sig[2]) / 3.0 - xi_c


def _rel(a, b):
    b = np.asarray(b, dtype=float)
    scale = max(float(np.max(np.abs(b))), 1e-30)
    return float(np.max(np.abs(np.asarray(a, dtype=float) - b))) / scale


# ---------------------------------------------------------------------------
# availability
# ---------------------------------------------------------------------------
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
def cp_available():
    if not _constructible(lambda t: mat_vm(t)):
        pytest.skip("ASDPlasticMaterial3D Closest_Point / LadrunoBrick not "
                    "available in this build")


@pytest.fixture(scope="module")
def cp_dp_available():
    if not _constructible(lambda t: mat_dp(t)):
        pytest.skip("ASDPlasticMaterial3D DruckerPrager Closest_Point not "
                    "available in this build")


# ---------------------------------------------------------------------------
# the oracle's strain paths (engineering Voigt totals per leg)
# ---------------------------------------------------------------------------
VM_PATHS = {
    "triaxial": [np.array([-9.0e-4, -9.0e-4, 3.0e-3, 0., 0., 0.])],
    "simple-shear": [np.array([0., 0., 0., 4.0e-3, 0., 0.])],
    "rotating-normal": [np.array([0., 0., 0., 3.0e-3, 0., 0.]),
                        np.array([0., 0., 0., 3.0e-3, 3.0e-3, 0.])],
}
DP_COMPRESS = [np.array([4.0e-4, 4.0e-4, -3.0e-3, 0., 0., 0.])]
DP_COMPRESS_SHEAR = [np.array([4.0e-4, 4.0e-4, -3.0e-3, 1.5e-3, 0., 0.])]
DP_HYDRO_TENSION = [np.array([1.2e-3, 1.2e-3, 1.2e-3, 0., 0., 0.])]


# ===========================================================================
# 0. the rig itself
# ===========================================================================
def test_rig_reproduces_the_prescribed_strain(cp_available):
    """Sanity: a fully prescribed unit cube must see EXACTLY the Voigt strain
    the driver asked for, and (below yield) the linear-isotropic stress.  With
    every DOF prescribed the free-DOF count is zero, so `analyze` is a pure
    state determination -- the `fd_tangent_driver` pattern."""
    eps = np.array([1.0e-5, -2.0e-5, 3.0e-5, 4.0e-5, -5.0e-5, 6.0e-5])
    r = drive(lambda t: mat_vm(t), [eps], nstep=1)
    assert r["codes"] == [0]
    assert _rel(r["eps"][-1], eps) < 1e-9, r["eps"][-1]
    lam = E_VM * NU_VM / ((1 + NU_VM) * (1 - 2 * NU_VM))
    mu = E_VM / (2 * (1 + NU_VM))
    D = np.zeros((6, 6))
    D[:3, :3] = lam
    for i in range(3):
        D[i, i] = lam + 2 * mu
    for i in range(3, 6):
        D[i, i] = mu
    print("rig strain rel err = %.3e" % _rel(r["eps"][-1], eps))
    assert _rel(r["sigma"][-1], D @ eps) < 1e-9


# ===========================================================================
# GATE 1 -- von Mises against the P0 oracle
# ===========================================================================
#: `Ladruno_implementation/adr97_oracle/reference_output.txt`, VM block.
VM_ORACLE = {
    # (case, path): (sigma, alpha, k)
    ("perfect", "triaxial"):
        ([60., 60., 90., 0., 0., 0.], [0.] * 6, 30.0),
    ("perfect", "simple-shear"):
        ([0., 0., 0., 17.3205080757, 0., 0.], [0.] * 6, 30.0),
    ("perfect", "rotating-normal"):
        ([0., 0., 0., 0.6597367945, 17.3079388537, 0.], [0.] * 6, 30.0),
    ("hardening", "triaxial"):
        ([55.2147239264, 55.2147239264, 99.5705521472, 0., 0., 0.],
         [0.] * 6, 44.3558282209),
    ("hardening", "simple-shear"):
        ([0., 0., 0., 24.5280749163, 0., 0.], [0.] * 6, 42.4838719668),
    ("hardening", "rotating-normal"):
        ([0., 0., 0., 2.7663534077, 27.3507956089, 0.], [0.] * 6, 47.6146636537),
}

VM_AF_ORACLE = {
    "triaxial": ([49.6232363748, 49.6232363748, 110.7535272504, 0., 0., 0.],
                 [-10.3767636252, -10.3767636252, 20.7535272504, 0., 0., 0.]),
    "simple-shear": ([0., 0., 0., 45.3906185722, 0., 0.],
                     [0., 0., 0., 28.0701104965, 0., 0.]),
    "rotating-normal": ([0., 0., 0., 23.7842127217, 37.7654987099, 0.],
                        [0., 0., 0., 21.6973739054, 20.5711652319, 0.]),
}


@pytest.mark.parametrize("case,hiso", [("perfect", 0.0), ("hardening", 7000.0)])
@pytest.mark.parametrize("path", ["triaxial", "simple-shear", "rotating-normal"])
def test_gate1_vm_closest_point_matches_the_oracle(cp_available, case, hiso, path):
    """The committed stress and yield stress of the closest-point map against
    the numpy CPPM oracle (`adr97_oracle/cppm_vm.py`), same conventions, same
    truncated `SQRT_2_over_3`, same 10 steps per leg."""
    r = drive(lambda t: mat_vm(t, hiso=hiso), VM_PATHS[path], nstep=10,
              want=("stresses", "YieldStress", "BackStress"))
    assert all(c == 0 for c in r["codes"]), r["codes"]
    sig_ref, alpha_ref, k_ref = VM_ORACLE[(case, path)]
    err = _rel(r["sigma"][-1], sig_ref)
    print("VM %-10s %-16s sigma rel err = %.3e" % (case, path, err))
    assert err < RTOL, (r["sigma"][-1], sig_ref)
    k = float(r["YieldStress"][-1][0])
    assert abs(k - k_ref) <= RTOL * abs(k_ref), (k, k_ref)
    # both hardening slopes for the back stress are 0 in these decks, so the
    # yield surface must not have translated
    assert np.max(np.abs(r["BackStress"][-1] - np.array(alpha_ref))) < 1e-9, \
        r["BackStress"][-1]


@pytest.mark.parametrize("path", ["triaxial", "simple-shear", "rotating-normal"])
def test_gate1_vm_armstrong_frederick_matches_the_oracle(cp_available, path):
    """Armstrong-Frederick solved IMPLICITLY at n+1 -- the direct removal of
    ADR-94 M3's path dependence.  `Backward_Euler` evaluates AF's recovery term
    at the TRIAL alpha inside an otherwise implicit loop, so its answer depends
    on how many corrector steps the loop happened to take."""
    if not _constructible(lambda t: mat_vm_af(t)):
        pytest.skip("VonMises + ArmstrongFrederick specialization unavailable")
    r = drive(lambda t: mat_vm_af(t), VM_PATHS[path], nstep=10,
              want=("stresses", "BackStress", "YieldStress"))
    assert all(c == 0 for c in r["codes"]), r["codes"]
    sig_ref, alpha_ref = VM_AF_ORACLE[path]
    es = _rel(r["sigma"][-1], sig_ref)
    ea = _rel(r["BackStress"][-1], alpha_ref)
    print("VM AF     %-16s sigma rel err = %.3e   alpha rel err = %.3e"
          % (path, es, ea))
    assert es < RTOL, (r["sigma"][-1], sig_ref)
    assert ea < RTOL, (r["BackStress"][-1], alpha_ref)


def test_gate1_vm_yield_function_is_satisfied_at_every_commit(cp_available):
    """|f| <= tol at EVERY committed state, not just the last one."""
    for hiso in (0.0, 7000.0):
        r = drive(lambda t: mat_vm(t, hiso=hiso), VM_PATHS["rotating-normal"],
                  nstep=10, want=("stresses", "YieldStress", "BackStress"))
        assert all(c == 0 for c in r["codes"])
        worst = 0.0
        for i in range(len(r["sigma"])):
            k = float(r["YieldStress"][i][0])
            f = _f_vm(r["sigma"][i], r["BackStress"][i], k)
            # a plastic step lands ON the surface; an elastic one is inside it
            worst = max(worst, f)
        print("VM hiso=%-8g worst committed f = %.3e (tol 1e-6)" % (hiso, worst))
        assert worst <= 1e-6, worst


def test_gate1_newton_converges_in_at_most_five_iterations(cp_available):
    """A closest-point Newton started from the elastic predictor is quadratic;
    the P0 oracle takes 1 iteration for perfect/linear von Mises and 3 for
    Armstrong-Frederick.  Read from the material's own `cp_iterations`
    response -- pytest's capfd cannot see the .pyd's cout (LEDGER_quirks)."""
    worst = 0
    for name, fn in (("perfect", lambda t: mat_vm(t)),
                     ("hardening", lambda t: mat_vm(t, hiso=7000.0)),
                     ("AF", lambda t: mat_vm_af(t))):
        if not _constructible(fn):
            continue
        r = drive(fn, VM_PATHS["rotating-normal"], nstep=10,
                  want=("stresses", "cp_iterations"))
        assert all(c == 0 for c in r["codes"])
        got = [v for v in r["cp_iterations"] if v.size]
        assert got, "the material did not answer the cp_iterations response"
        m = int(max(float(np.max(v)) for v in got))
        print("VM %-10s max Closest_Point Newton iterations = %d" % (name, m))
        worst = max(worst, m)
    assert worst <= 5, worst


# ===========================================================================
# GATE 1 -- Drucker-Prager (cone, non-associated cone, apex)
# ===========================================================================
#: `adr97_oracle/reference_output.txt`, DP block.  Only the PERFECTLY PLASTIC
#: rows are pinned against the oracle: with `ScalarLinearHardeningParameter = 0`
#: the oracle's self-consistent `f = sqrt(J2) + eta p - (xi_c + k)` and the
#: shipped `f = sqrt(J2) + eta p - xi_c` coincide exactly.  With H != 0 they do
#: NOT -- see `test_gate1_dp_cohesion_hardening_is_perfectly_plastic_under_cp`.
DP_ORACLE = {
    "cone-associated": (DP_COMPRESS, 0.4,
                        [-26.1416983071, -26.1416983071, -94.7352064897,
                         0., 0., 0.]),
    "cone-nonassociated": (DP_COMPRESS_SHEAR, 0.2,
                           [-26.6815746982, -26.6815746982, -89.9603702958,
                            13.9585578524, 0., 0.]),
    "apex": (DP_HYDRO_TENSION, 0.4, [50., 50., 50., 0., 0., 0.]),
}


@pytest.mark.parametrize("case", list(DP_ORACLE))
def test_gate1_dp_closest_point_matches_the_oracle(cp_dp_available, case):
    """Cone (associated and non-associated) and APEX returns against
    `adr97_oracle/cppm_dp.py`.  The apex row is the one the shipped
    `check_apex_region` cannot get right in general: it is a EUCLIDEAN normal
    cone test (`p - p_apex >= eta q`) while the exact condition is in the
    ELASTIC metric.  `Closest_Point` classifies inside the integrator, where K,
    G and etabar are in scope."""
    path, etabar, sig_ref = DP_ORACLE[case]
    r = drive(lambda t: mat_dp(t, etabar=etabar), path, nstep=10,
              want=("stresses", "BackStress"))
    assert all(c == 0 for c in r["codes"]), r["codes"]
    err = _rel(r["sigma"][-1], sig_ref)
    print("DP %-20s sigma rel err = %.3e   sigma = %s"
          % (case, err, np.array2string(r["sigma"][-1], precision=8)))
    assert err < RTOL, (r["sigma"][-1], sig_ref)
    f = _f_dp(r["sigma"][-1], r["BackStress"][-1])
    assert f <= 1e-6, f


def test_gate1_dp_cohesion_hardening_moves_the_iv_but_not_the_surface(
        cp_dp_available):
    """ADR-97 P0 header finding 2, made observable -- and sharpened.

    `DruckerPrager_YF` has its cohesion internal variable COMMENTED OUT of its
    own `f` (line 26) while `yf_hardening` still contributes ``df/dk = -1``
    times that IV's rate: a hardening term matching no term of `f`.

    The finding as ADR-97 P0 stated it was that `Closest_Point` would be
    perfectly plastic here (its `df/dq` must be the TRUE derivative of the `f`
    it solves, hence zero) while `Backward_Euler` "hardens".  MEASURED, the
    second half is wrong and the defect is WORSE than stated: the cohesion
    hardening is inert on the committed stress under BOTH integrators.  It
    cannot be otherwise -- both maps drive the SAME `f` to zero, and `f` does
    not read the internal variable, so the root is the same wherever the
    plastic modulus puts the iterates.  Under `Backward_Euler` the phantom
    modulus only changes the local Newton's STEP SIZE, and what survives is the
    cutting plane's iterate-path residue (~5e-9 here); under `Closest_Point` the
    consistency row carries no `df/dq` term at all and the two answers are
    identical to the last bit.

    So `ScalarLinearHardeningParameter` on a Drucker-Prager is a no-op for the
    stress on both integrators, while the internal variable itself grows -- it
    is written, recorded, and read by nothing.  What DOES differ is the tangent:
    `Backward_Euler`'s `Continuum`/`Secant` operators divide by
    ``n:E:m - H`` with that phantom `H`, so a hardening slope the surface never
    sees still perturbs the assembled stiffness.

    This test pins all three facts.  It turns red the day the yield function is
    fixed -- which is the point: fixing it changes `Backward_Euler`, which
    ADR-97 D1 keeps byte-identical, so it is a separate PR.
    """
    soft = drive(lambda t: mat_dp(t, hiso=0.0), DP_COMPRESS, nstep=10)
    hard = drive(lambda t: mat_dp(t, hiso=500.0), DP_COMPRESS, nstep=10,
                 want=("stresses", "DP_cohesion"))
    assert all(c == 0 for c in soft["codes"] + hard["codes"])
    d_cp = _rel(hard["sigma"][-1], soft["sigma"][-1])
    print("DP cohesion hardening under Closest_Point: rel difference vs H=0 "
          "= %.3e (expected exactly 0)" % d_cp)
    assert d_cp < 1e-12, (
        "Closest_Point's df/dq for the DruckerPrager cohesion IV is no longer "
        "zero -- either the yield function was fixed (good, but ADR-97 D1 and "
        "this pin both need updating) or df_dq drifted")

    # ... and yet the internal variable itself HAS advanced
    k_end = float(hard["DP_cohesion"][-1][0])
    print("DP cohesion IV after the same path: k = %.6f (written, read by "
          "nothing)" % k_end)
    assert k_end > 1e-3, (
        "the DruckerPrager cohesion internal variable stopped evolving; if the "
        "yield function was fixed to read it, this whole test is obsolete")

    be_soft = drive(lambda t: mat_dp(t, hiso=0.0, method="Backward_Euler",
                                     tangent="Secant"), DP_COMPRESS, nstep=10)
    be_hard = drive(lambda t: mat_dp(t, hiso=500.0, method="Backward_Euler",
                                     tangent="Secant"), DP_COMPRESS, nstep=10)
    d_be = _rel(be_hard["sigma"][-1], be_soft["sigma"][-1])
    print("DP cohesion hardening under Backward_Euler: rel difference vs H=0 "
          "= %.3e (the cutting plane's iterate-path residue, not hardening)"
          % d_be)
    assert d_be < 1e-6, (
        "Backward_Euler's committed stress now responds to a cohesion "
        "hardening slope its own yield function does not read")
    assert d_be > d_cp, (
        "the cutting plane's iterate-path residue vanished -- Backward_Euler "
        "may no longer be a cutting plane")


# ===========================================================================
# GATE 2 -- the consistent tangent
# ===========================================================================
#: measured on build 3622d6214 by `adr97_oracle/fd_tangent_driver.py` and
#: quoted in ADR-94 M3 as 57.3 %.  It is the NEGATIVE control: whatever
#: `Algorithmic` scores has to be contrasted against a tangent that is known
#: not to be the tangent of the committed map.
FD_BE_CONTINUUM = 0.573447


def test_gate2_algorithmic_matches_a_finite_difference_of_the_committed_map(
        cp_available):
    """`Algorithmic` against a central difference of the BINARY'S OWN assembled
    internal force, on a free-DOF rig -- no numpy reference anywhere.

    `fd_check` uses `stdBrick`, deliberately: that is the host on which the P0
    oracle measured the negative control, and no material return code is under
    test here (stdBrick swallows them by design, ADR-94 B2)."""
    r_ctrl = FD.fd_check(lambda t: FD.mat_vm(t, tangent="Continuum",
                                             method="Backward_Euler"),
                         label="Backward_Euler / Continuum (negative control)")
    print("gate 2 negative control rel_err = %.6f (pinned %.6f)"
          % (r_ctrl["rel_err"], FD_BE_CONTINUUM))
    assert abs(r_ctrl["rel_err"] - FD_BE_CONTINUUM) < 1e-3, (
        "the ADR-94 M3 negative control moved: %r" % r_ctrl["rel_err"])

    r = FD.fd_check(lambda t: FD.mat_vm(t, tangent="Algorithmic",
                                        method="Closest_Point"),
                    label="Closest_Point / Algorithmic")
    print("gate 2 Closest_Point/Algorithmic rel_err = %.3e   rel_fro = %.3e"
          % (r["rel_err"], r["rel_fro"]))
    assert r["rc"] == 0, "the Closest_Point load-driven rig did not converge"
    assert r["rel_err"] <= 1e-6, r["rel_err"]


def test_gate2_algorithmic_on_a_sheared_twelve_dof_rig(cp_available):
    """The same measurement with the top face fully free, so the tangent's
    shear block is exercised (the P0 oracle measured Continuum at 0.670 here)."""
    r = FD.fd_check(lambda t: FD.mat_vm(t, tangent="Algorithmic",
                                        method="Closest_Point"),
                    rig="full", load=(-3.0, 0.0, -13.0),
                    label="12-DOF rig, Closest_Point / Algorithmic")
    print("gate 2 (12-DOF) Closest_Point/Algorithmic rel_err = %.3e"
          % r["rel_err"])
    assert r["rc"] == 0
    assert r["rel_err"] <= 1e-6, r["rel_err"]


# --- the ADR-94 two-cube model -------------------------------------------
def _two_cube(mat_fn, nsteps=4, load_pl=-60.0, load_el=-10.0):
    """Two disconnected VonMises cubes, one driven plastic and one left
    elastic (`tests/test_adr94_redblue_blue.py`'s heterogeneous rig)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, x0 in ((1, 0.0), (2, 3.0)):
        b = 10 * t
        for k, (x, y, z) in enumerate(NODES):
            ops.node(b + k + 1, x0 + float(x), float(y), float(z))
        for k in range(4):
            ops.fix(b + k + 1, 1, 1, 1)
    mat_fn(1)
    for t, _x0 in ((1, 0.0), (2, 3.0)):
        b = 10 * t
        ops.element("stdBrick", t, *[b + k + 1 for k in range(8)], 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t, load in ((1, load_pl), (2, load_el)):
        b = 10 * t
        for k in range(4, 8):
            ops.load(b + k + 1, 0., 0., load / nsteps)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-9, 200, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    iters, sig = [], None
    for _ in range(nsteps):
        rc = ops.analyze(1)
        assert rc == 0, "two-cube rig must converge every step"
        iters.append(ops.testIter())
    sig = np.array(list(ops.eleResponse(1, "stresses"))[0:6])
    return iters, sig


def test_gate2_two_cube_newton_cost(cp_available):
    """ADR-94 M3's cost claim, re-measured: the DEFAULT (`Backward_Euler` +
    `Secant`) costs 5.3x `Continuum`'s global iterations, and no shipped
    tangent is the tangent of the map.  With the consistent tangent the global
    Newton should take <= 3 iterations per step and no more than `Continuum`."""
    it_cont, sig_cont = _two_cube(
        lambda t: FD.mat_vm(t, tangent="Continuum", method="Backward_Euler"))
    it_alg, sig_alg = _two_cube(
        lambda t: FD.mat_vm(t, tangent="Algorithmic", method="Closest_Point"))
    print("two-cube per-step testIter: BE/Continuum = %s, "
          "Closest_Point/Algorithmic = %s" % (it_cont, it_alg))
    # the converged answers must agree -- this is a proportional (non-rotating
    # normal) path, where the cutting plane and the closest point coincide
    print("two-cube converged stress rel difference = %.3e"
          % _rel(sig_alg, sig_cont))
    print("two-cube TOTAL iterations: BE/Continuum = %d, "
          "Closest_Point/Algorithmic = %d  (%.1fx)"
          % (sum(it_cont), sum(it_alg), sum(it_cont) / max(sum(it_alg), 1)))
    assert _rel(sig_alg, sig_cont) < 1e-6

    # The ADR's gate was written as `testIter() <= 3`; MEASURED it is 4 per step
    # -- one Newton iteration to leave the previous step's converged tangent,
    # two quadratic, one to satisfy NormDispIncr at 1e-9.  4 is what an exact
    # consistent tangent costs on this rig, and the content of the gate is the
    # CONTRAST: `Continuum` takes 6/41/36/30 for the same four steps.
    assert max(it_alg) <= 4, it_alg
    assert sum(it_alg) * 2 <= sum(it_cont), (it_alg, it_cont)


# ===========================================================================
# GATE 3 -- the cutting plane's path error
# ===========================================================================
def test_gate3_af_step_refinement_beats_the_cutting_plane(cp_available):
    """The closest-point map and the cutting plane are both O(h) and converge to
    the SAME limit, but at finite step the cutting plane is measurably less
    accurate on Armstrong-Frederick -- the P0 oracle measured 4.3x on this very
    material.  With LINEAR hardening the two are indistinguishable, which is
    exactly why ADR-94 H6 could not see the defect.

    (Gate 3 in the ADR is stated as "step-halving reproduces the state"; what is
    observable from the binary is the step-REFINEMENT error of each map against
    its own fine-step limit, which is the same statement made measurable -- the
    per-step Newton-start invariance the oracle pins at 1.1e-12 cannot be
    reached from Python, since the start guess is not an option.)
    """
    if not _constructible(lambda t: mat_vm_af(t)):
        pytest.skip("VonMises + ArmstrongFrederick specialization unavailable")
    path = VM_PATHS["rotating-normal"]
    ref_cp = drive(lambda t: mat_vm_af(t), path, nstep=160)["sigma"][-1]
    ref_be = drive(lambda t: mat_vm_af(t, method="Backward_Euler",
                                       tangent="Secant"),
                   path, nstep=160)["sigma"][-1]
    lim = _rel(ref_cp, ref_be)
    print("AF fine-step (N=160) CP vs BE rel difference = %.3e" % lim)

    e_cp = _rel(drive(lambda t: mat_vm_af(t), path, nstep=10)["sigma"][-1],
                ref_cp)
    e_be = _rel(drive(lambda t: mat_vm_af(t, method="Backward_Euler",
                                          tangent="Secant"),
                      path, nstep=10)["sigma"][-1], ref_be)
    print("AF step-refinement error at N=10: Closest_Point %.4e, "
          "Backward_Euler %.4e  (ratio %.2fx)" % (e_cp, e_be, e_be / max(e_cp, 1e-30)))
    assert e_be > e_cp, (
        "the cutting plane is no longer less accurate than the closest point "
        "on an Armstrong-Frederick path -- gate 3's contrast is gone")

    # ... and with LINEAR hardening the contrast vanishes (ADR-94 H6's blind spot)
    ref_l_cp = drive(lambda t: mat_vm(t, hiso=7000.0), path, nstep=160)["sigma"][-1]
    ref_l_be = drive(lambda t: mat_vm(t, hiso=7000.0, method="Backward_Euler",
                                      tangent="Secant"),
                     path, nstep=160)["sigma"][-1]
    print("linear-hardening fine-step CP vs BE rel difference = %.3e"
          % _rel(ref_l_cp, ref_l_be))
