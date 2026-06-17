"""
Zone-A regression: black-box finite-difference DAMPING sensitivity.

Pins the capability scoped in Ladruno_implementation/damping_sensitivity/:
the response gradient w.r.t. a damping parameter via finite differences, with
NO DDM derivative and NO vanilla edits. Uses vanilla OpenSees only (rayleigh,
Viscous, modalDamping, UniformDamping), so this is upstreamable (zone_a).

Two layers:
  * FDSensitivity unit tests (pure math, no engine): central/forward exactness on
    a quadratic, caching + honest n_solves, step-plateau monotonicity.
  * Damping-channel tests: SDOF at resonance, gradient of steady amplitude vs the
    parameter-free oracle  dU/dp = -U/p  (exact because U(p) = A/p when the damping
    coefficient scales linearly with the parameter). Covers element-Rayleigh,
    Viscous material, modal damping, nodal alphaM, and a UniformDamping object.

Run:  <py3.12> -m pytest tests/test_fd_damping_sensitivity.py -v
(the dist/bin opensees.pyd is built for CPython 3.12)
"""
import os
import sys
import math

import pytest

# --- bootstrap the locally-built engine -----------------------------------
DIST = r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\dist\bin"
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")
if os.path.isdir(DIST):
    try:
        os.add_dll_directory(DIST)
    except (FileNotFoundError, OSError):
        pass
    if DIST not in sys.path:
        sys.path.insert(0, DIST)

# --- import the helper from the damping_sensitivity bundle -----------------
_BUNDLE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "Ladruno_implementation", "damping_sensitivity",
)
if _BUNDLE not in sys.path:
    sys.path.insert(0, _BUNDLE)
from fd_sensitivity import FDSensitivity, fd_gradient  # noqa: E402

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


# ==========================================================================
# Layer 1 — FDSensitivity unit tests (pure math, no engine)
# ==========================================================================
zone_a = pytest.mark.zone_a


@zone_a
def test_central_difference_exact_on_quadratic():
    # f(x) = x^2 -> f'(3) = 6 ; central difference is exact for quadratics.
    fd = FDSensitivity(lambda p: p[0] ** 2)
    assert math.isclose(fd.gradient([3.0])[0], 6.0, rel_tol=1e-7)


@zone_a
def test_gradient_vector_multiparam():
    # f(x,y) = x^2 + 3y -> grad = [2x, 3]
    fd = FDSensitivity(lambda p: p[0] ** 2 + 3.0 * p[1])
    g = fd.gradient([4.0, 10.0])
    assert math.isclose(g[0], 8.0, rel_tol=1e-6)
    assert math.isclose(g[1], 3.0, rel_tol=1e-6)


@zone_a
def test_forward_scheme_reuses_base_solve():
    # forward scheme: 1 base + 1 per component = 3 solves for 2 params.
    fd = FDSensitivity(lambda p: p[0] ** 2 + p[1] ** 2, scheme="forward")
    fd.gradient([1.0, 1.0])
    assert fd.n_solves == 3


@zone_a
def test_cache_avoids_resolves():
    # Re-running gradient at the same point must hit the cache (no new solves).
    fd = FDSensitivity(lambda p: p[0] ** 2)
    fd.gradient([2.0])
    n_after_first = fd.n_solves
    fd.gradient([2.0])
    assert fd.n_solves == n_after_first  # all cache hits the second time


@zone_a
def test_step_study_plateau_monotone_improvement():
    # On a cubic, central-diff error ~ h^2: shrinking the step toward the plateau
    # must not get worse. f=x^3 -> f'(1)=3.
    fd = FDSensitivity(lambda p: p[0] ** 3)
    rows = fd.step_study([1.0], rel_steps=(1e-1, 1e-2, 1e-3))
    errs = [abs(g - 3.0) for _, g in rows]
    assert errs[0] > errs[1] > errs[2]


# ==========================================================================
# Layer 2 — damping-channel sensitivity (SDOF at resonance)
# ==========================================================================
M = 1.0
F_N = 1.0
WN = 2.0 * math.pi * F_N
K = WN * WN * M
F0 = 1.0
DT = 0.002
N_CYCLES = 120
T = N_CYCLES / F_N
NSTEPS = int(T / DT)


def _sdof(extra_doRayleigh=True):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0); ops.fix(1, 1)
    ops.node(2, 0.0); ops.mass(2, M)
    ops.uniaxialMaterial("Elastic", 1, K)
    dr = ("-doRayleigh", 1) if extra_doRayleigh else ()
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1, *dr)
    ops.timeSeries("Trig", 1, 0.0, T, 1.0 / F_N, "-factor", F0)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0)


def _amplitude():
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("ProfileSPD")
    ops.test("NormDispIncr", 1e-10, 25)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    tail = T - 10.0 / F_N
    umax = 0.0
    for _ in range(NSTEPS):
        ops.analyze(1, DT)
        if ops.getTime() >= tail:
            umax = max(umax, abs(ops.nodeDisp(2, 1)))
    return umax


def _fwd_rayleigh_alphaM(p):
    (aM,) = p
    _sdof()
    ops.rayleigh(aM, 0.0, 0.0, 0.0)
    return _amplitude()


def _fwd_viscous(p):
    (C,) = p
    # built directly: spring (Elastic) + dashpot (Viscous) in parallel in one
    # zeroLength (two materials on dir 1). Damping is the material's own force.
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0); ops.fix(1, 1)
    ops.node(2, 0.0); ops.mass(2, M)
    ops.uniaxialMaterial("Elastic", 1, K)
    ops.uniaxialMaterial("Viscous", 2, C, 1.0)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, 2, "-dir", 1, 1)
    ops.timeSeries("Trig", 1, 0.0, T, 1.0 / F_N, "-factor", F0)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0)
    return _amplitude()


def _fwd_modal(p):
    (xi,) = p
    _sdof()
    ops.eigen("-fullGenLapack", 1)
    ops.modalDamping(xi)
    return _amplitude()


def _fwd_uniform(p):
    (xi,) = p
    _sdof()
    ops.damping("Uniform", 1, xi, 0.5 * WN, 2.0 * WN)
    ops.region(1, "-ele", 1, "-damp", 1)
    return _amplitude()


# parameter-free oracle: U(p)=A/p at resonance  =>  dU/dp = -U/p (exact)
def _oracle(U, p):
    return -U / p


@zone_a
@pytest.mark.parametrize("name,fwd,p0", [
    ("rayleigh_alphaM", _fwd_rayleigh_alphaM, 2.0 * 0.05 * WN),
    ("viscous_C",       _fwd_viscous,        2.0 * 0.08 * M * WN),
    ("modal_xi",        _fwd_modal,          0.05),
    ("nodal_alphaM",    _fwd_rayleigh_alphaM, 2.0 * 0.07 * WN),
    ("uniform_object",  _fwd_uniform,        0.05),
])
def test_damping_channel_gradient_matches_oracle(name, fwd, p0):
    fd = FDSensitivity(fwd, rel_step=3e-3, scheme="central")
    U = fwd([p0])
    g = fd.gradient([p0])[0]
    g_oracle = _oracle(U, p0)
    assert g < 0.0, f"{name}: more damping must reduce resonant amplitude"
    rel = abs(g - g_oracle) / abs(g_oracle)
    assert rel < 5e-3, f"{name}: FD grad {g:.6e} vs oracle {g_oracle:.6e} (rel {rel:.2%})"


@zone_a
def test_uniform_object_input_zeta_is_not_realised_xi():
    # Honest pin of the Damping-object quirk: input zeta != realised SDOF xi.
    # (Documents WHY the simple xi-oracle does not apply to UniformDamping.)
    U = _fwd_uniform([0.05])
    xi_eff = F0 / (2.0 * K * U)
    assert xi_eff < 0.05, f"realised xi_eff={xi_eff:.4f} should be < input 0.05"
