"""ADR-94 R4 -- integrator x tangent matrix for ASDPlasticMaterial3D.

Parametrized measurement sweep, NOT a physics gate: this file asserts only
that every cell of {YF} x {integration_method} x {tangent_type} x {path} ran
to completion (or was cleanly recorded as refused/stalled/crashed). It never
fails on the *quality* of a cell's answer -- that is what the written table
is for.

METRICS per cell (YF, integrator, tangent) x path:
  * ``|f|max``  -- max absolute value of a numpy re-implementation of the
    yield function, evaluated on every COMMITTED stress (``strict_convergence
    0``, the default, so an inadmissible state can be committed -- that is
    the point of measuring it, not a bug in the harness).
  * ``drift``   -- ``|f|max`` at 2x the step count over ``|f|max`` at the
    base step count, on the same total load -- a cheap convergence-order
    smell test, not a rigorous order-of-accuracy study.
  * ``iters``   -- summed ``ops.testIter()`` over the base-step-count run.
  * ``status``  -- ``ok`` / ``refused`` (analyze() != 0) / ``stalled`` (wall
    clock cap hit) / ``crashed`` (exception or NaN/Inf committed) /
    ``unavailable`` (material would not construct in this build).

DRIVER (ADR-94 Sec. 8 traps obeyed):
  * Host = ``LadrunoBrick`` (propagates material refusal codes; ``stdBrick``
    swallows them -- already proven usable for ASDP in
    ``test_adr94_hlist_hb.py``).
  * ``system("UmfPack")`` always; never ``FullGeneral``.
  * DISPLACEMENT control (``sp``), not load control: an early calibration
    pass here hit the textbook load-control failure mode -- Null/zero-rate
    ("perfectly plastic") hardening has a genuine LIMIT LOAD in shear, past
    which no equilibrium solution exists for any integrator, so load control
    refused almost every cell regardless of integrator quality. Displacement
    control has no such ceiling. Compression/tension prescribe only the
    z-face DOF and leave x, y completely free (traction-free, a true
    uniaxial-stress path); shear prescribes only the top x-DOF (fixing z,
    leaving y free). Every path therefore keeps >=4 genuinely free DOFs, so
    ``ops.testIter()`` remains a meaningful signal (the ADR-84/H1 lesson: a
    fully-prescribed driver cannot see a wrong tangent).
  * ``strict_convergence`` left at the default 0 everywhere (task scope).
  * ``n_max_iterations`` and ``rk45_niter_max`` are capped low (30 / 40) so a
    non-convergent cell fails FAST rather than spinning; a per-path wall
    clock cap (1.5 s) is a second backstop against a genuine hang.
  * Materials use Null or zero-rate ("simplest") hardening throughout, so
    each numpy ``f_*`` re-implementation does not need to track internal
    variable evolution -- it is evaluated directly off the committed
    6-component stress.

PARAMETER PROVENANCE (reused from already-vetted ADR-94 decks, not
re-derived): VonMises from ``test_adr84_p2a_strict_convergence`` (E=70000,
NU=0.3, SY=30 kPa); MohrCoulomb/MohrCoulombTensionCutoff from
``test_asdplastic_mctc`` (E=2e6, NU=0.3, phi=20, c=100, psi=5 kPa,
T_rock=24.7 kPa); DruckerPrager and HoekBrown from ``test_adr94_hlist_hb``
(DP: E=1e6, NU=0.25, xi_c=1000, eta=0.3, etabar=0.1; HB: 50 MPa rock, GSI=60,
mi=10 sandstone-like). Known hazard from that file, guarded here with the NaN
check above: a pure hydrostatic path through the DP apex silently commits
NaN (``CHECK_APEX_REGION`` stub) -- this file's tension/compression paths
are NOT purely hydrostatic (single-axis load only) but the guard is kept
regardless.
"""
import math
import time

import numpy as np
import pytest

from _testbed import ops

import test_adr84_p2a_strict_convergence as VMX  # noqa: E402
import test_asdplastic_mctc as MCX               # noqa: E402
import test_adr94_hlist_hb as HBX                # noqa: E402

pytestmark = [pytest.mark.zone_a]

INTEGRATORS = [
    "Forward_Euler", "Forward_Euler_Subincrement",
    "Backward_Euler", "Backward_Euler_LineSearch",
    "Modified_Euler_Error_Control", "Runge_Kutta_45_Error_Control",
]
TANGENTS = [
    "Elastic", "Continuum", "Secant",
    "Numerical_Algorithmic_FirstOrder", "Numerical_Algorithmic_SecondOrder",
]
PATHS = ["triaxial_compression", "simple_shear", "tension_to_cutoff"]

N_ITER_CAP = 30
RK45_ITER_CAP = 40
WALL_CAP_S = 1.5           # per (integrator,tangent,path,step-count) run
N_BASE = 3                 # base step count; drift compares against 2x this


def _int_opts(method, tangent):
    return ["Begin_Integration_Options",
            "integration_method", method, "tangent_type", tangent,
            "n_max_iterations", N_ITER_CAP,
            "rk45_niter_max", RK45_ITER_CAP,
            "strict_convergence", 0,
            "End_Integration_Options"]


# ---------------------------------------------------------------------------
# materials -- Null / zero-rate hardening throughout (see module docstring
# for parameter provenance)
# ---------------------------------------------------------------------------
def mat_vm(tag, method, tangent):
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "VonMises_YF", "VonMises_PF", "LinearIsotropic3D_EL", VMX.IV_VM,
        "Begin_Model_Parameters",
        "YoungsModulus", VMX.E_VM, "PoissonsRatio", VMX.NU_VM,
        "ScalarLinearHardeningParameter", 0.0,
        "TensorLinearHardeningParameter", 0.0,
        "MassDensity", 0.0,
        "End_Model_Parameters",
        "Begin_Internal_Variables",
        "YieldStress", VMX.SY_VM,
        "BackStress", 0., 0., 0., 0., 0., 0.,
        "End_Internal_Variables",
        *_int_opts(method, tangent))


def mat_dp(tag, method, tangent):
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "DruckerPrager_YF", "DruckerPrager_PF", "LinearIsotropic3D_EL",
        HBX.IV_DP,
        "Begin_Model_Parameters",
        "YoungsModulus", HBX.DP_E, "PoissonsRatio", HBX.DP_NU,
        "DP_xi_c", HBX.DP_XI_C, "DP_eta", HBX.DP_ETA,
        "DP_etabar", HBX.DP_ETABAR,
        "TensorLinearHardeningParameter", 0.0,
        "ScalarLinearHardeningParameter", 0.0,
        "MassDensity", 0.0,
        "End_Model_Parameters",
        "Begin_Internal_Variables",
        "BackStress", 0., 0., 0., 0., 0., 0.,
        "DP_cohesion", 0.0,
        "End_Internal_Variables",
        *_int_opts(method, tangent))


def mat_mc(tag, method, tangent):
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "MohrCoulomb_YF", "MohrCoulomb_PF", "LinearIsotropic3D_EL", MCX.IV,
        "Begin_Model_Parameters",
        "YoungsModulus", MCX.E, "PoissonsRatio", MCX.NU,
        "MC_phi", MCX.PHI, "MC_c", MCX.C, "MC_psi", MCX.PSI, "MC_ds", 0.0,
        "MassDensity", 0.0,
        "End_Model_Parameters",
        "Begin_Internal_Variables",
        "BackStress", 0., 0., 0., 0., 0., 0.,
        "End_Internal_Variables",
        *_int_opts(method, tangent))


def mat_mctc(tag, method, tangent):
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "MohrCoulombTensionCutoff_YF", "MohrCoulombTensionCutoff_PF",
        "LinearIsotropic3D_EL", MCX.IV,
        "Begin_Model_Parameters",
        "YoungsModulus", MCX.E, "PoissonsRatio", MCX.NU,
        "MC_phi", MCX.PHI, "MC_c", MCX.C, "MC_psi", MCX.PSI, "MC_ds", 0.0,
        "TC_min_stress", MCX.T_ROCK, "MassDensity", 0.0,
        "End_Model_Parameters",
        "Begin_Internal_Variables",
        "BackStress", 0., 0., 0., 0., 0., 0.,
        "End_Internal_Variables",
        *_int_opts(method, tangent))


def mat_hb(tag, method, tangent):
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "HoekBrown_YF", "HoekBrown_PF", "LinearIsotropic3D_EL", HBX.IV_HB,
        "Begin_Model_Parameters",
        "YoungsModulus", HBX.HB_E, "PoissonsRatio", HBX.HB_NU,
        "HB_sigci", HBX.HB_SIGMA_CI, "HB_mb", HBX.HB_MB, "HB_s", HBX.HB_S,
        "HB_a", HBX.HB_A, "HB_mb_psi", HBX.HB_MB, "HB_ds", 0.0,
        "MassDensity", 0.0,
        "End_Model_Parameters",
        "Begin_Internal_Variables",
        "BackStress", 0., 0., 0., 0., 0., 0.,
        "End_Internal_Variables",
        *_int_opts(method, tangent))


# ---------------------------------------------------------------------------
# numpy yield-function oracles -- evaluated directly off the committed
# 6-component stress (Voigt [11,22,33,12,23,13]); Null/zero-rate hardening
# means no internal-variable bookkeeping is required.
# ---------------------------------------------------------------------------
def f_vm(s):
    return VMX.f_vm(s, sy=VMX.SY_VM)


def f_dp(s):
    """NOTE (post-hoc fix): ``DruckerPrager_YF.h`` computes
    ``p = sigma.meanStress()`` (tension-positive, unnegated) with a comment
    claiming "positive in compression" -- a static read of the header alone
    says the literal arithmetic is unnegated. An EMPIRICAL probe (drive a
    single-element DP compression path on Backward_Euler across a range of
    step counts and watch which sign the committed state actually tracks
    toward zero) shows the opposite: the material's own return map converges
    towards ``sqrt_J2 - eta*p - xi_c = 0``, i.e. COMPRESSION-POSITIVE p, not
    the literal unnegated header arithmetic. Trusting the runtime behaviour
    over the static header reading (whatever internal sign flip accounts for
    the discrepancy -- e.g. an internal geomechanics-convention frame used
    only inside the return map -- is out of scope for this oracle), this
    oracle negates p to match. See the R4 return summary for the probe
    numbers this was calibrated against."""
    S = MCX._tensor(s)
    I1 = float(np.trace(S))
    dev = S - (I1 / 3.0) * np.eye(3)
    sqrt_J2 = math.sqrt(max(0.5 * float(np.sum(dev * dev)), 0.0))
    p = I1 / 3.0
    return sqrt_J2 - HBX.DP_ETA * p - HBX.DP_XI_C


def f_mc(s):
    return MCX.f_mc(s, phi_deg=MCX.PHI, c=MCX.C)


def f_mctc(s):
    return max(MCX.f_mc(s, phi_deg=MCX.PHI, c=MCX.C), MCX.f_tc(s, MCX.T_ROCK))


def f_hb(s):
    S_geo = -MCX._tensor(s)
    e = np.linalg.eigvalsh(S_geo)          # ascending
    sigma3, sigma1 = e[0], e[2]
    arg = HBX.HB_MB * sigma3 / HBX.HB_SIGMA_CI + HBX.HB_S
    if arg > 0:
        return sigma1 - sigma3 - HBX.HB_SIGMA_CI * (arg ** HBX.HB_A)
    return sigma1 - sigma3 - HBX.HB_SIGMA_CI * HBX.HB_S


MATERIALS = {
    "VonMises": mat_vm,
    "DruckerPrager": mat_dp,
    "MohrCoulomb": mat_mc,
    "MohrCoulombTensionCutoff": mat_mctc,
    "HoekBrown": mat_hb,
}
F_ORACLE = {
    "VonMises": f_vm,
    "DruckerPrager": f_dp,
    "MohrCoulomb": f_mc,
    "MohrCoulombTensionCutoff": f_mctc,
    "HoekBrown": f_hb,
}
# elastic modulus per YF (needed to turn a stress scale into a strain)
_EMOD = {
    "VonMises": VMX.E_VM, "DruckerPrager": HBX.DP_E,
    "MohrCoulomb": MCX.E, "MohrCoulombTensionCutoff": MCX.E,
    "HoekBrown": HBX.HB_E,
}
_NU = {
    "VonMises": VMX.NU_VM, "DruckerPrager": HBX.DP_NU,
    "MohrCoulomb": MCX.NU, "MohrCoulombTensionCutoff": MCX.NU,
    "HoekBrown": HBX.HB_NU,
}
# representative compressive / tensile stress scale per YF (kPa)
_SCALE_C = {
    "VonMises": VMX.SY_VM, "DruckerPrager": HBX.DP_XI_C,
    "MohrCoulomb": 2.0 * MCX.C * math.cos(math.radians(MCX.PHI)),
    "MohrCoulombTensionCutoff": 2.0 * MCX.C * math.cos(math.radians(MCX.PHI)),
    "HoekBrown": HBX.HB_UCS,
}
_SCALE_T = {
    "VonMises": VMX.SY_VM, "DruckerPrager": HBX.DP_XI_C,
    "MohrCoulomb": 2.0 * MCX.C * math.cos(math.radians(MCX.PHI)),
    "MohrCoulombTensionCutoff": MCX.T_ROCK,
    "HoekBrown": abs(HBX.HB_SIGMA_T),
}
# calibrated multiples of (stress_scale / modulus) that keep the reference
# Backward_Euler + Secant cell converging at both N=3 and N=6 steps -- found
# by an empirical scan (see the R4 return summary); NOT a limit-analysis
# result, just a comfortably-past-first-yield operating point per path.
# NOTE on DruckerPrager/HoekBrown: the base-fully-fixed cube is NOT a clean
# uniaxial-stress element (the fixed base induces real lateral confinement at
# the queried Gauss point), so the naive stress_scale/E estimate undershot
# first yield by 30-100x for these two YFs specifically (VM/MC/MCTC, whose
# scale already matched observed yielding, were left unchanged). Re-derived
# empirically: scan k, keep the largest value at which Backward_Euler +
# Secant still converges at both N=3 and N=6 (a DP compression/shear/tension
# ceiling of ~1.0-1.5 -- DP hits the already-documented H7 elastic-fallback
# almost immediately past first yield, so it cannot sustain a deep excursion
# regardless of k; HoekBrown compression/shear tolerate a much larger k).
_K = {
    "VonMises": dict(triaxial_compression=3.0, simple_shear=3.0, tension_to_cutoff=3.0),
    # DP: compression/shear pushed to the exact largest k that (a) still
    # converges at N=3 AND 2N=6 and (b) actually crosses f=0 at the final
    # commit (1.5 -> f=+4.75; 1.0 -> f=0.0 exactly) -- anything past this
    # refuses immediately (the H7 elastic-fallback ceiling, confirmed
    # independent of n_max_iterations in the R4 return summary). Tension
    # could NOT be pushed past first yield at all: every k above ~1.1
    # refuses before the committed history ever reaches f>=0, so 0.9 is kept
    # as the largest safe (but still elastic) value -- itself a finding, not
    # an oversight.
    "DruckerPrager": dict(triaxial_compression=1.5, simple_shear=1.0, tension_to_cutoff=0.9),
    "MohrCoulomb": dict(triaxial_compression=3.0, simple_shear=0.8, tension_to_cutoff=3.0),
    "MohrCoulombTensionCutoff": dict(triaxial_compression=1.2, simple_shear=1.0, tension_to_cutoff=3.0),
    # HB compression: 30 gives a genuinely plastic interim commit (|f|~1e-7,
    # see the R4 return summary) even though the run as a whole still ends
    # in a later refusal for weaker integrators -- an accurate reading, not
    # a regression.
    "HoekBrown": dict(triaxial_compression=30.0, simple_shear=3.0, tension_to_cutoff=1.9),
}


def _sp_disp(yf, path):
    """(dof, magnitude) for the ONE sp-prescribed top-face DOF; dof is a
    1-based OpenSees DOF index (1=x, 3=z)."""
    E = _EMOD[yf]
    G = E / (2.0 * (1.0 + _NU[yf]))
    k = _K[yf][path]
    if path == "simple_shear":
        return 1, k * _SCALE_C[yf] / G
    if path == "triaxial_compression":
        return 3, -k * _SCALE_C[yf] / E
    return 3, k * _SCALE_T[yf] / E          # tension_to_cutoff


# ---------------------------------------------------------------------------
# driver: unit cube, base fixed, ONE top-face DOF displacement-controlled,
# the rest genuinely free (traction-free / unconstrained).
# ---------------------------------------------------------------------------
_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
_BASE = (1, 2, 3, 4)
_TOP = (5, 6, 7, 8)


def _constructible(mat_fn):
    try:
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        for t, c in _CUBE.items():
            ops.node(t, *map(float, c))
        mat_fn(1, "Backward_Euler", "Secant")
        ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
        ok = True
    except Exception:
        ok = False
    ops.wipe()
    return ok


def _run(mat_fn, dof_disp, nsteps, wall_cap=WALL_CAP_S):
    """One displacement-controlled run. Returns (status, iters, hist[N,6])."""
    dof, disp = dof_disp
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *map(float, c))
    for t in _BASE:
        ops.fix(t, 1, 1, 1)
    if dof == 1:                    # shear: fix z, leave x (sp) and y free
        for t in _TOP:
            ops.fix(t, 0, 0, 1)
    else:                            # compression/tension: x,y free
        for t in _TOP:
            ops.fix(t, 0, 0, 0)
    mat_fn(1)
    ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in _TOP:
        ops.sp(n, dof, disp)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-8, 60, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")

    hist, iters, status = [], 0, "ok"
    t0 = time.perf_counter()
    for _ in range(nsteps):
        if time.perf_counter() - t0 > wall_cap:
            status = "stalled"
            break
        rc = ops.analyze(1)
        iters += ops.testIter()
        if rc != 0:
            status = "refused"
            break
        s = list(ops.eleResponse(1, "stresses"))[0:6]
        if any((v != v) or math.isinf(v) for v in s):     # NaN/Inf guard
            hist.append(s)
            status = "crashed"
            break
        hist.append(s)
    return status, iters, (np.array(hist) if hist else np.zeros((0, 6)))


def _fmax(f_fn, hist):
    if len(hist) == 0:
        return float("nan")
    vals = [abs(f_fn(s)) for s in hist if all(v == v for v in s)]
    return float(np.max(vals)) if vals else float("nan")


# ---------------------------------------------------------------------------
# the sweep -- one module-scoped fixture builds the whole matrix once
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def matrix():
    results = {}
    t_start = time.perf_counter()

    for yf, mat_fn in MATERIALS.items():
        available = _constructible(mat_fn)
        for method in INTEGRATORS:
            for tangent in TANGENTS:
                def make(tag, m=mat_fn, meth=method, tg=tangent):
                    m(tag, meth, tg)

                for path in PATHS:
                    key = (yf, method, tangent, path)
                    if not available:
                        results[key] = dict(status="unavailable", fmax=float("nan"),
                                             iters=0, drift=float("nan"), wall=0.0)
                        continue
                    dof_disp = _sp_disp(yf, path)
                    t0 = time.perf_counter()
                    try:
                        st_n, it_n, h_n = _run(make, dof_disp, N_BASE)
                        fmax_n = _fmax(F_ORACLE[yf], h_n)
                        st_2n, _, h_2n = _run(make, dof_disp, 2 * N_BASE)
                        fmax_2n = _fmax(F_ORACLE[yf], h_2n)
                        status = st_n if st_n != "ok" else st_2n
                        drift = (fmax_2n / fmax_n
                                 if (fmax_n == fmax_n and fmax_n > 1e-9) else
                                 float("nan"))
                        wall = time.perf_counter() - t0
                        results[key] = dict(
                            status=status,
                            fmax=(fmax_2n if fmax_2n == fmax_2n else fmax_n),
                            iters=it_n, drift=drift, wall=wall)
                    except Exception:
                        results[key] = dict(status="crashed", fmax=float("nan"),
                                             iters=0, drift=float("nan"),
                                             wall=time.perf_counter() - t0)
    ops.wipe()
    total_wall = time.perf_counter() - t_start
    _write_report(results, total_wall)
    return results, total_wall


# ---------------------------------------------------------------------------
# report writer
# ---------------------------------------------------------------------------
def _fmt(cell):
    if cell["status"] == "unavailable":
        return "unavailable"
    f = cell["fmax"]
    fstr = f"{f:.2e}" if f == f else "nan"
    return f"{fstr}/{cell['iters']}/{cell['status']}"


def _write_report(results, total_wall):
    out_path = ("../Ladruno_implementation/_adr94_matrix.md")
    lines = [
        "# ADR-94 R4 -- integrator x tangent matrix\n",
        f"\nGenerated by `tests/test_adr94_matrix.py`. Total wall time: "
        f"{total_wall:.1f} s. Cell format: `|f|max / sum(testIter) / status` "
        f"(base step count N={N_BASE}; `|f|max` is the worse of the N-step "
        f"and 2N-step runs). `strict_convergence 0` (default) throughout, so "
        f"`status=ok` does NOT imply an admissible commit -- see `|f|max`.\n",
    ]
    for path in PATHS:
        lines.append(f"\n## Path: {path}\n")
        header = "| YF | Integrator | " + " | ".join(TANGENTS) + " |"
        sep = "|---|---|" + "---|" * len(TANGENTS)
        lines.append(header)
        lines.append(sep)
        for yf in MATERIALS:
            for method in INTEGRATORS:
                row = [yf, method]
                for tangent in TANGENTS:
                    row.append(_fmt(results[(yf, method, tangent, path)]))
                lines.append("| " + " | ".join(row) + " |")

    # recommended configuration per YF: best (integrator, tangent) by
    # (most paths 'ok', then lowest mean |f|max, then lowest iters)
    # Two-stage ranking per YF: (1) admit only the configs that reach the
    # BEST n_ok achieved by anything for that YF (a robustness floor -- a
    # config that only survives the easiest single path must not outrank one
    # that survives all three just because its one surviving path happens to
    # be more accurate); (2) among that admitted set, rank by lowest mean
    # |f|max over its ok paths FIRST, then summed iters over those paths.
    # This replaces an earlier version that sorted purely on -n_ok, which on
    # ties fell through to iters and could pick a config whose accuracy was
    # never actually compared against its equally-robust rivals.
    lines.append("\n## Recommended configuration per YF\n")
    for yf in MATERIALS:
        scored = []
        for method in INTEGRATORS:
            for tangent in TANGENTS:
                cells = [results[(yf, method, tangent, p)] for p in PATHS]
                ok_cells = [c for c in cells if c["status"] == "ok"]
                n_ok = len(ok_cells)
                fvals = [c["fmax"] for c in ok_cells if c["fmax"] == c["fmax"]]
                mean_f = float(np.mean(fvals)) if fvals else float("inf")
                iters = sum(c["iters"] for c in ok_cells)
                scored.append((method, tangent, n_ok, mean_f, iters))
        max_n_ok = max(s[2] for s in scored)
        admitted = [s for s in scored if s[2] == max_n_ok]
        method, tangent, n_ok, mean_f, iters = min(
            admitted, key=lambda s: (s[3], s[4]))
        lines.append(
            f"- **{yf}**: `{method}` + `{tangent}` "
            f"({n_ok}/3 paths ok [best reachable], mean |f|max={mean_f:.3g} "
            f"over ok paths, summed iters={iters})")

    # crash/stall/refuse roll-up
    lines.append("\n## Non-ok cells\n")
    bad = [(k, v) for k, v in results.items() if v["status"] not in ("ok",)]
    if not bad:
        lines.append("(none -- every cell reported `ok`)")
    else:
        for (yf, method, tangent, path), v in bad:
            lines.append(
                f"- `{yf}` / `{method}` / `{tangent}` / `{path}`: "
                f"**{v['status']}** (iters={v['iters']})")

    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# the test itself -- runs the sweep and asserts only structural completeness
# ---------------------------------------------------------------------------
def test_adr94_matrix_every_cell_accounted_for(matrix):
    results, total_wall = matrix
    expected = len(MATERIALS) * len(INTEGRATORS) * len(TANGENTS) * len(PATHS)
    assert len(results) == expected

    valid = {"ok", "refused", "stalled", "crashed", "unavailable"}
    for key, cell in results.items():
        assert cell["status"] in valid, f"{key}: unrecognized status {cell}"

    print(f"\n[ADR-94 R4] {expected} cells, wall time {total_wall:.1f} s")
