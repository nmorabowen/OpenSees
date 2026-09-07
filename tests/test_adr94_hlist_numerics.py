"""ADR-94 R1-A (numerics lane) — sentinel tests for H1, H6, H7, H8 of the
``ASDPlasticMaterial3D`` review.

Every assertion here pins a DEFECT **as currently observed**, except where the
docstring says ``REFUTED``.  A future fix therefore turns the test RED and
announces itself; it is not silently absorbed.  Verdicts and the deciding
numbers live in ``Ladruno_implementation/_adr94_hlist_R1A.md``; the numpy
references live in ``Ladruno_implementation/adr94_oracle/hex8_tangent.py``.

WHAT EACH H IS
--------------
* **H1 — class-static ``Stiffness``.**  ``Stiffness``, ``dsigma``,
  ``depsilon_elpl``, ``intersection_stress/strain`` are ``static`` members of
  the ``ASDPlasticMaterial3D<E,Y,P,tag>`` specialization (declared 4116-4120,
  defined 4150+), so ONE 6x6 is shared by every Gauss point of every element of
  every material tag of that YF/PF/EL combination.  ``getTangent()`` (698) just
  copies it.  Host elements ``setTrialStrain`` all Gauss points in ``update()``
  and call ``getTangent()`` in a SEPARATE later loop (``Brick.cpp`` 1069 vs
  1201), so the whole model is assembled with the tangent of whichever Gauss
  point was integrated last.  CONFIRMED, blocker.
* **H6 — ``Backward_Euler`` is a cutting plane, not a closest-point map.**
  ``n``, ``m``, ``H`` are re-evaluated at the current iterate and the stress is
  corrected incrementally (2244-2320).  The accuracy half of the claim is
  REFUTED on VonMises (where the two maps coincide exactly); the tangent half
  is CONFIRMED — no ``tangent_type`` reproduces the consistent tangent.
* **H7 — the ``dLambda + deltaLambda < 0`` fallback (2298-2303)** prints
  ``PLASTIC INCONSISTENCY - ELASTIC STEP!``, sets ``Stiffness = Eelastic`` and
  ``return 0``.  CONFIRMED and worse than hypothesised: it fires at iteration
  0, so the commit is the ELASTIC PREDICTOR (not a partial iterate), it is
  committed with ``f >> tol``, and ``strict_convergence`` does not gate it.
* **H8 — ``Backward_Euler_LineSearch``** hardcodes ``max_iter = 30`` and
  ``tol_rel = 1e-8``, tests a LINEAR prediction of Phi (so a Newton direction
  always accepts alpha = 1), never sees ``strict_convergence``, and its
  "substepping" loop halves ``dEps`` and solves ONE reduced increment rather
  than a chain, overwriting ``TrialStrain`` with ``CommitStrain + dEps/2^k``.
  CONFIRMED.

TRAPS OBEYED (ADR-94 §8)
------------------------
* ``system("UmfPack")`` everywhere; the tangent is read with
  ``printA -sparse -ret`` (the dense ``-ret`` path needs ``FullGeneral``,
  whose ``getA()`` is the only non-null one — see ``OpenSeesCommands.cpp``
  2718 — and ``FullGeneral`` is banned here).
* Every rig leaves free DOFs, so a wrong tangent is observable.
* Refusal gates use ``TenNodeTetrahedron`` (``stdBrick`` swallows return
  codes).  ``stdBrick`` is used only where no refusal is under test — it is
  the host the H1 row cites.
"""
import math
import os
import sys

import numpy as np
import pytest

from _testbed import ops

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir, "Ladruno_implementation",
                                "adr94_oracle"))
import hex8_tangent as O  # noqa: E402

import test_asdplastic_mctc as M  # noqa: E402
import test_adr84_p2a_strict_convergence as P  # noqa: E402

pytestmark = [pytest.mark.zone_a]

# ---------------------------------------------------------------------------
# VonMises control model, kPa-like magnitudes
# ---------------------------------------------------------------------------
E_VM, NU_VM, SY_VM, H_VM = 70000.0, 0.3, 30.0, 7000.0
# Softening steeper than n:E:m = 2G = E/(1+nu) = 53846 -- the H7 reproducer.
H_SOFT = -120000.0
IV_VM = ("BackStress(TensorLinearHardeningFunction):"
         "YieldStress(ScalarLinearHardeningFunction):")

# Free DOFs of the uniaxial-strain hex rig, in the 24-DOF oracle numbering:
# nodes 5..8, z component.  Verified against OpenSees by
# ``test_H6_oracle_matches_opensees_elastic_stiffness``.
FREE = [14, 17, 20, 23]


def mat_vm(tag, tangent="Continuum", method="Backward_Euler",
           hiso=H_VM, niter=100, strict=None):
    extra = ["strict_convergence", int(strict)] if strict is not None else []
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "VonMises_YF", "VonMises_PF", "LinearIsotropic3D_EL", IV_VM,
        "Begin_Model_Parameters",
        "YoungsModulus", E_VM, "PoissonsRatio", NU_VM,
        "ScalarLinearHardeningParameter", hiso,
        "TensorLinearHardeningParameter", 0.0, "MassDensity", 0.0,
        "End_Model_Parameters",
        "Begin_Internal_Variables", "YieldStress", SY_VM,
        "BackStress", 0., 0., 0., 0., 0., 0., "End_Internal_Variables",
        "Begin_Integration_Options",
        "integration_method", method, "tangent_type", tangent,
        "n_max_iterations", int(niter), *extra,
        "End_Integration_Options",
    )


def _sparse_K(n):
    """Assembled tangent as a dense ``n x n`` (``printA`` re-forms it first)."""
    d = ops.printA("-sparse", "-ret")
    K = np.zeros((n, n))
    for i, j, v in zip(d["rowIndices"], d["colIndices"], d["values"]):
        K[i, j] += v
    return K


def _rel(a, b):
    return float(np.max(np.abs(a - b)) / max(float(np.max(np.abs(b))), 1e-30))


# ---------------------------------------------------------------------------
# rig A — uniaxial-strain unit cube.  x and y fixed everywhere, z fixed on the
# base, z LOADED (not prescribed) on the top: 4 free DOFs, and the strain field
# is homogeneous, so H1 cannot contaminate an H6 measurement.
# ---------------------------------------------------------------------------
def _cube_build(mat_fn, load_z, nsteps, tol=1e-12):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for k, (x, y, z) in enumerate(O.NODES):
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
        ops.load(k, 0., 0., load_z)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", tol, 60, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")


# Load per top node: 4 P over unit area.  Yield at |sigma_zz| = 52.5, so -60
# per node (sigma_zz = -240) is deep in the plastic range.
P_PLASTIC = -60.0


def _cube_run(mat_fn, nsteps, load_z=P_PLASTIC):
    _cube_build(mat_fn, load_z, nsteps)
    iters = 0
    for _ in range(nsteps):
        rc = ops.analyze(1)
        iters += ops.testIter()
        if rc != 0:
            return rc, iters, None, None
    eps = np.array(list(ops.eleResponse(1, "strains"))[0:6])
    sig = np.array(list(ops.eleResponse(1, "stresses"))[0:6])
    return 0, iters, eps, sig


# ---------------------------------------------------------------------------
# rig B — two DISCONNECTED unit cubes in one model, one plastic one elastic.
# ---------------------------------------------------------------------------
def _cubes_build(tags_x0, loads, tangent="Continuum"):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, x0 in tags_x0:
        b = 10 * t
        for k, (x, y, z) in enumerate(O.NODES):
            ops.node(b + k + 1, x0 + float(x), float(y), float(z))
        for k in range(4):
            ops.fix(b + k + 1, 1, 1, 1)
    mat_vm(1, tangent)
    for t, x0 in tags_x0:
        b = 10 * t
        ops.element("stdBrick", t, *[b + k + 1 for k in range(8)], 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t, load in loads:
        b = 10 * t
        for k in range(4, 8):
            ops.load(b + k + 1, 0., 0., load)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-11, 80, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


# 4 x (-10) = -40 kPa: elastic (yield needs 52.5).  4 x (-60): plastic.
LOAD_EL, LOAD_PL = -10.0, -60.0
NDOF_CUBE = 12          # a free-standing cube: 4 top nodes x 3 DOF


def _cubes_K(tags_x0, loads):
    _cubes_build(tags_x0, loads)
    rc = ops.analyze(1)
    return rc, _sparse_K(NDOF_CUBE * len(tags_x0))


# ---------------------------------------------------------------------------
# availability
# ---------------------------------------------------------------------------
def _vm_constructible():
    try:
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        for k, (x, y, z) in enumerate(O.NODES):
            ops.node(k + 1, float(x), float(y), float(z))
        mat_vm(1)
        ops.element("stdBrick", 1, *range(1, 9), 1)
        ok = True
    except Exception:
        ok = False
    ops.wipe()
    return ok


@pytest.fixture(scope="module")
def vm_available():
    if not _vm_constructible():
        pytest.skip("ASDPlasticMaterial3D / VonMises_YF not available")


@pytest.fixture(scope="module")
def mc_available():
    if not P._tet_constructible(lambda t: P.mat_mc(t)):
        pytest.skip("ASDPlasticMaterial3D / MohrCoulomb_YF / "
                    "TenNodeTetrahedron not available")


# ===========================================================================
# H1 — the class-static Stiffness
# ===========================================================================
@pytest.mark.t0m
def test_H1_oracle_matches_opensees_elastic_stiffness():
    """Validate the numpy hex-8 oracle before anything is measured with it.

    ``hex8_K(C_elastic)`` must equal the assembled OpenSees ``stdBrick``
    tangent to machine precision, which also pins the ``FREE`` DOF map used by
    the H6 tests.  Measured: 4.4e-16.
    """
    _cube_build(lambda t: ops.nDMaterial("ElasticIsotropic", t, E_VM, NU_VM),
                -1.0e-4, 1)
    ops.analyze(1)
    K_ops = _sparse_K(4)
    K_ref = O.hex8_K(O.C_elastic(E_VM, NU_VM))[np.ix_(FREE, FREE)]
    assert _rel(K_ops, K_ref) < 1e-12


@pytest.mark.t0m
def test_H1_one_static_tangent_is_shared_by_every_element(vm_available):
    """CONFIRMED (blocker).  Two disconnected cubes of the SAME ASDP
    specialization, one plastic one elastic, are assembled with ONE tangent.

    Measured on ``52314165a``: the two elements' stand-alone tangents differ by
    13.6%, yet inside the combined model the two diagonal blocks are equal to
    round-off, and the plastic element's block is 13.2% away from its own
    correct tangent.  Assemble the block sum of the two singles and you get a
    DIFFERENT matrix from the two-element model -- which is the whole claim.

    A fix (per-instance ``Stiffness``) makes ``blocks_identical`` false and
    turns this test red.
    """
    _, K_pl = _cubes_K([(1, 0.0)], [(1, LOAD_PL)])
    _, K_el = _cubes_K([(2, 3.0)], [(2, LOAD_EL)])
    assert _rel(K_pl, K_el) > 0.1, "the two states must be genuinely different"

    _, K_both = _cubes_K([(1, 0.0), (2, 3.0)], [(1, LOAD_PL), (2, LOAD_EL)])
    n = NDOF_CUBE
    blk_pl, blk_el = K_both[:n, :n], K_both[n:, n:]

    # DEFECT: both elements got the same 6x6 tangent.
    assert _rel(blk_pl, blk_el) < 1e-9, "H1 fixed? the blocks now differ"
    # DEFECT: the plastic element's block is NOT its own tangent...
    assert _rel(blk_pl, K_pl) > 0.1
    # ...it is the ELASTIC element's, i.e. the last one integrated.
    assert _rel(blk_el, K_el) < 1e-9


@pytest.mark.t0m
def test_H1_shared_tangent_follows_the_last_element_integrated(vm_available):
    """CONFIRMED (blocker).  WHICH tangent everybody gets is decided by domain
    iteration order, not by the element: swap which cube is plastic and the
    wrong block swaps with it.

    Together with the previous test this is the ADR-75b blocker: the result of
    an assembly depends on the order Gauss points were integrated in, so the
    material can never be safe under threaded assembly.
    """
    n = NDOF_CUBE
    # element 1 plastic, element 2 elastic -> everyone gets ELEMENT 2's tangent
    _, K_a = _cubes_K([(1, 0.0), (2, 3.0)], [(1, LOAD_PL), (2, LOAD_EL)])
    # element 1 elastic, element 2 plastic -> everyone gets ELEMENT 2's again
    _, K_b = _cubes_K([(1, 0.0), (2, 3.0)], [(1, LOAD_EL), (2, LOAD_PL)])

    _, K_el_alone = _cubes_K([(2, 3.0)], [(2, LOAD_EL)])
    _, K_pl_alone = _cubes_K([(1, 0.0)], [(1, LOAD_PL)])

    # case a: the shared tangent is the ELASTIC one (element 2 ran last)
    assert _rel(K_a[n:, n:], K_el_alone) < 1e-9
    # case b: the shared tangent is now the PLASTIC one, and it is imposed on
    # the ELASTIC element (block 1) as well
    assert _rel(K_b[:n, :n], K_b[n:, n:]) < 1e-9
    assert _rel(K_b[:n, :n], K_el_alone) > 0.1
    assert _rel(K_b[:n, :n], K_pl_alone) < 0.05


# ===========================================================================
# H6 — cutting plane vs closest point, and the tangent options
# ===========================================================================
@pytest.mark.t0m
def test_H6_backward_euler_is_exact_for_von_mises(vm_available):
    """REFUTED (the accuracy half, on VonMises).

    For ``f = ||dev s|| - sqrt(2/3) sy`` with associated flow and LINEAR
    isotropic hardening the return direction ``n = s/||s||`` does not rotate
    during the correction, so ``Phi(lambda)`` is exactly linear and the cutting
    plane lands on the closest-point answer in ONE Newton step.  Measured
    against the numpy closed form (``vm_radial_return``) at ``dlam = 2.3e-3``:
    relative stress error 3.7e-14.

    Consequence for the ADR: the H6 row's proposed "convergence-order study on
    VM+linear hardening" is DEGENERATE -- it can only ever measure backward
    Euler's own first order, never cutting-plane-vs-closest-point.  A YF whose
    normal rotates (Lode-dependent MC, kinematic hardening, stress-dependent
    elasticity) is required for that.
    """
    rc, _, eps, sig = _cube_run(lambda t: mat_vm(t, "Continuum"), 1)
    assert rc == 0
    sig_ref, dlam, _ = O.vm_radial_return(np.zeros(6), eps, E_VM, NU_VM,
                                          SY_VM, H_VM)
    assert dlam > 1e-4, "the step must be genuinely plastic"
    assert _rel(sig, sig_ref) < 1e-10


@pytest.mark.t0m
def test_H6_no_tangent_option_reproduces_the_consistent_tangent(vm_available):
    """CONFIRMED (major).  On a homogeneous plastic state (so H1 is invisible)
    NONE of the five ``tangent_type`` options is the consistent tangent of the
    return map.  Measured vs ``hex8_K(vm_consistent_tangent)`` at
    ``dEps_zz = -3.6e-3``:

        Continuum                          57.3%   3 global Newton iters
        Secant  (the DEFAULT)              79.9%  16
        Elastic                           102.5%  23
        Numerical_Algorithmic_FirstOrder   31.0%   4
        Numerical_Algorithmic_SecondOrder  31.0%   4

    The numerical pair are closest yet still 31% out because they differentiate
    ``compute_local_stress()`` -- a THIRD map that is not the one committed.
    The operational consequence is the iteration count: the shipped default
    costs 5.3x the iterations of ``Continuum`` on this step.
    """
    errs, iters = {}, {}
    for tg in ("Continuum", "Secant", "Elastic",
               "Numerical_Algorithmic_FirstOrder",
               "Numerical_Algorithmic_SecondOrder"):
        rc, it, eps, _ = _cube_run(lambda t, g=tg: mat_vm(t, g), 1)
        assert rc == 0
        C_alg = O.vm_consistent_tangent(np.zeros(6), eps, E_VM, NU_VM,
                                        SY_VM, H_VM)
        K_ref = O.hex8_K(C_alg)[np.ix_(FREE, FREE)]
        errs[tg] = _rel(_sparse_K(4), K_ref)
        iters[tg] = it

    assert min(errs.values()) > 0.25, "a consistent tangent appeared: %r" % errs
    assert errs["Continuum"] > 0.4
    assert errs["Secant"] > 0.6
    assert errs["Elastic"] > 0.9
    assert errs["Numerical_Algorithmic_FirstOrder"] > 0.25
    # the DEFAULT (Secant) is the expensive one
    assert iters["Secant"] >= 4 * iters["Continuum"]


@pytest.mark.t0m
def test_H6_continuum_tangent_error_is_first_order_in_the_step(vm_available):
    """CONFIRMED (major).  The ``Continuum`` operator is the ``dlam -> 0``
    limit of the consistent tangent, so its error is O(dEps) and only vanishes
    as the step does -- it is never the tangent of the map actually used.
    Measured (nsteps, dEps_zz, K error): 1/-3.6e-3/57.3%, 2/-1.8e-3/46.2%,
    4/-9.1e-4/31.3%, 8/-4.5e-4/19.1%, 16/-2.3e-4/10.7% -- a halving per
    halving of the step.
    """
    errs = {}
    for ns in (2, 8):
        rc, _, eps, _ = _cube_run(lambda t: mat_vm(t, "Continuum"), ns)
        assert rc == 0
        eprev = eps * (ns - 1) / ns
        s_prev, _, sy_prev = O.vm_radial_return(np.zeros(6), eprev, E_VM,
                                                NU_VM, SY_VM, H_VM)
        C_alg = O.vm_consistent_tangent(s_prev, eps - eprev, E_VM, NU_VM,
                                        sy_prev, H_VM)
        errs[ns] = _rel(_sparse_K(4), O.hex8_K(C_alg)[np.ix_(FREE, FREE)])
    assert errs[2] > 0.35
    assert errs[8] < 0.5 * errs[2] + 0.05     # first-order decay, not O(1)
    assert errs[8] > 0.10                     # but still far from consistent


# ===========================================================================
# H7 — the "PLASTIC INCONSISTENCY - ELASTIC STEP!" fallback
# ===========================================================================
def _f_vm(s, sy=SY_VM):
    S = np.array([[s[0], s[3], s[5]], [s[3], s[1], s[4]], [s[5], s[4], s[2]]])
    d = S - np.trace(S) / 3.0 * np.eye(3)
    return math.sqrt(float(np.sum(d * d))) - math.sqrt(2.0 / 3.0) * sy


def _soft_history(strict=None, nsteps=4):
    """Softening VonMises with ``H_iso`` steeper than ``n:E:m``: at the first
    local iteration ``dPhi/dlambda = H - n:E:m > 0`` with ``Phi > 0``, so
    ``deltaLambda < 0`` and the 2298 branch fires immediately."""
    _cube_build(lambda t: mat_vm(t, "Continuum", hiso=H_SOFT, strict=strict),
                P_PLASTIC, nsteps, tol=1e-10)
    rows = []
    for _ in range(nsteps):
        rc = ops.analyze(1)
        eps = np.array(list(ops.eleResponse(1, "strains"))[0:6])
        sig = np.array(list(ops.eleResponse(1, "stresses"))[0:6])
        rows.append((rc, eps, sig))
        if rc != 0:
            break
    return rows


@pytest.mark.t0m
def test_H7_inconsistency_branch_commits_the_elastic_predictor(vm_available):
    """CONFIRMED, and WORSE than the H7 row states (major).

    The row predicts the commit is a "partially corrected iterate".  Measured:
    the branch fires at iteration 0, where ``dLambda == 0`` and nothing has
    been corrected, so the commit is the ELASTIC PREDICTOR **exactly**
    (difference 0.000e+00) -- and the step is accepted with ``return 0``.  The
    material then behaves purely elastically for the whole run while ``f_VM``
    grows monotonically: measured ``f_VM`` at the four commits = +3.50, +31.49,
    +59.49, +87.48 kPa against ``f_absolute_tol = 1e-6``.

    Blast radius: any specialization whose hardening modulus softens faster
    than ``n:E:m`` (here ``H_iso = -120000`` vs ``2G = 53846``) -- VonMises,
    DruckerPrager and MohrCoulomb all accept a negative
    ``ScalarLinearHardeningParameter``.

    (The branch also prints ``PLASTIC INCONSISTENCY - ELASTIC STEP!`` on
    ``cout``.  That line is NOT asserted on: it is one of the 161 ``cout``
    sites of H12 and does not survive pytest's fd capture reliably.  The
    signature asserted instead is unique to this branch -- no other path
    returns 0 with a commit that is the elastic predictor to machine precision
    while ``f > 0``.)
    """
    rows = _soft_history()

    Ce = O.C_elastic(E_VM, NU_VM)
    s_prev, e_prev, fs = np.zeros(6), np.zeros(6), []
    for rc, eps, sig in rows:
        assert rc == 0, "the branch reports SUCCESS -- that is the defect"
        predictor = s_prev + Ce @ (eps - e_prev)
        assert float(np.max(np.abs(sig - predictor))) < 1e-9 * max(
            1.0, float(np.max(np.abs(predictor))))
        fs.append(_f_vm(sig))
        s_prev, e_prev = sig, eps

    assert fs[-1] > 50.0, "committed states are far outside the surface"
    assert all(b > a for a, b in zip(fs, fs[1:])), "and drifting further out"


@pytest.mark.t0m
def test_H7_strict_convergence_does_not_gate_the_inconsistency_branch(
        vm_available):
    """CONFIRMED (major).  ADR-84 P2a gated the exhaustion-accept and the
    f-decreasing exit; this is an EIGHTH silent-accept site in the same
    (default) integrator and ``strict_convergence`` never sees it -- the
    branch returns 0 before the ``be_strict`` check at 2338.  Measured: with
    ``strict_convergence 1`` all four steps still return 0 with
    ``f_VM = +87.5``, bit-identical to flag-off.
    """
    off = _soft_history(strict=None)
    on = _soft_history(strict=1)
    assert [r[0] for r in on] == [0, 0, 0, 0]
    assert _rel(on[-1][2], off[-1][2]) < 1e-12
    assert _f_vm(on[-1][2]) > 50.0


# ===========================================================================
# H8 — Backward_Euler_LineSearch
# ===========================================================================
def _mc_tet(method, strict=None, niter=None):
    """The ADR-84 tet rig (a genuinely plastic MC leg) on a chosen
    integrator.  ``TenNodeTetrahedron`` propagates the material return code."""
    def mat(tag):
        extra = []
        if strict is not None:
            extra += ["strict_convergence", int(strict)]
        if niter is not None:
            extra += ["n_max_iterations", int(niter)]
        ops.nDMaterial(
            "ASDPlasticMaterial3D", tag,
            "MohrCoulomb_YF", "MohrCoulomb_PF", "LinearIsotropic3D_EL", M.IV,
            "Begin_Model_Parameters",
            "YoungsModulus", M.E, "PoissonsRatio", M.NU,
            "MC_phi", M.PHI, "MC_c", M.C, "MC_psi", M.PSI, "MC_ds", 0.0,
            "MassDensity", 0.0, "End_Model_Parameters",
            "Begin_Internal_Variables",
            "BackStress", 0., 0., 0., 0., 0., 0., "End_Internal_Variables",
            "Begin_Integration_Options",
            "integration_method", method, *extra, "End_Integration_Options",
        )
    codes, hist = P._drive(lambda: P._tet_build(mat), P.TET_NSTEPS)
    return codes, hist


@pytest.fixture(scope="module")
def ls_runs(mc_available):
    return {
        ("BE", 100, None): _mc_tet("Backward_Euler", niter=100),
        ("LS", 2, None): _mc_tet("Backward_Euler_LineSearch", niter=2),
        ("LS", 100, None): _mc_tet("Backward_Euler_LineSearch", niter=100),
        ("LS", 100, 1): _mc_tet("Backward_Euler_LineSearch", niter=100,
                                strict=1),
    }


@pytest.mark.t0m
def test_H8_line_search_ignores_n_max_iterations(ls_runs):
    """CONFIRMED (major).  ``Backward_Euler_LineSearch`` hardcodes
    ``max_iter = 30`` (2385) instead of reading ``INT_OPT_n_max_iterations``,
    so the option is silently inert: ``n_max_iterations`` 2 and 100 give
    BIT-IDENTICAL histories.  On the same rig ``Backward_Euler`` with
    ``n_max_iterations 2`` is the ADR-84 exhaustion reproducer (worst
    ``f_MC = 77.6`` vs 6.3e-4 at 100), so the option is not inert in general --
    only in this integrator.
    """
    c2, h2 = ls_runs[("LS", 2, None)]
    c100, h100 = ls_runs[("LS", 100, None)]
    assert c2 == c100
    assert h2.shape == h100.shape
    assert np.array_equal(h2, h100)


@pytest.mark.t0m
def test_H8_line_search_ignores_strict_convergence(ls_runs):
    """CONFIRMED (major).  ``strict_convergence`` is read only inside
    ``Backward_Euler`` (2081, 2086, 2184, 2338); ``Backward_Euler_LineSearch``
    has no ``be_strict`` at all, so the fork's one loud-failure switch is a
    no-op the moment a user selects this integrator.  Measured: flag on and
    flag off are bit-identical.
    """
    c_off, h_off = ls_runs[("LS", 100, None)]
    c_on, h_on = ls_runs[("LS", 100, 1)]
    assert c_off == c_on
    assert np.array_equal(h_off, h_on)


@pytest.mark.t0m
def test_H8_line_search_is_less_robust_than_plain_backward_euler(ls_runs):
    """CONFIRMED (major).  The name promises robustness; measured on the ADR-84
    MC tet leg it completes 2 of 20 steps where plain ``Backward_Euler``
    completes 20 of 20.

    Reading the code, the "line search" cannot help: the acceptance test is on
    a LINEAR prediction ``Phi + dPhi/dlambda * dl`` (2477-2486), so for an
    unclipped Newton direction ``dl = -alpha*Phi/dPhi`` it reduces to
    ``|1-alpha| <= 1 - 1e-4*alpha``, true for every ``alpha`` in (0,1] --
    ``alpha = 1`` is always accepted on the first try and no backtracking ever
    happens.  Nor does the "substepping": on failure it halves ``dEps`` and
    solves ONE reduced increment (2578-2596) instead of chaining substeps,
    overwriting ``TrialStrain`` with ``CommitStrain + dEps/2^k`` and returning
    SUCCESS for a strain the element never asked for.
    """
    c_be, h_be = ls_runs[("BE", 100, None)]
    c_ls, h_ls = ls_runs[("LS", 100, None)]
    ok_be = sum(1 for c in c_be if c == 0)
    ok_ls = sum(1 for c in c_ls if c == 0)
    assert ok_be == P.TET_NSTEPS
    assert ok_ls < 5, "BE_LS became robust -- re-read the H8 verdict"
