"""
Sensitivity of peak inter-story drift to DAMPING, on a 2-DOF shear building.
Black-box finite differences via fd_sensitivity.fd_gradient. Zero vanilla edits.

Model: 2-story shear frame (ndm=1, ndf=1). Floor masses m1,m2; story springs
k1,k2 (zeroLength Elastic). Classical Rayleigh damping matched to a target
ratio xi at modes 1 and 2. Driven by harmonic ground acceleration near mode 1
(UniformExcitation) so damping strongly shapes the response.

Two studies:
  (A) d(peak story-1 drift)/d(xi)        -- scalar damping-ratio parameter
  (B) gradient w.r.t. [alphaM, betaK]    -- the raw Rayleigh coefficients (vector)
"""
import os, sys, math
from fd_sensitivity import fd_gradient, step_study

DIST = r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\dist\bin"
os.environ["LADRUNO_OPENSEES_QUIET"] = "1"
os.add_dll_directory(DIST)
sys.path.insert(0, DIST)
import opensees as ops  # noqa: E402

# ---- structure ---------------------------------------------------------
m1, m2 = 2.0, 1.5
k1, k2 = 1200.0, 800.0          # story stiffnesses
dt = 0.005
n_cycles = 60


def _build_frame():
    """Build the bare frame (no damping, no analysis). Returns (w1, w2)."""
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(0, 0.0); ops.fix(0, 1)     # ground
    ops.node(1, 0.0); ops.mass(1, m1)
    ops.node(2, 0.0); ops.mass(2, m2)
    ops.uniaxialMaterial("Elastic", 1, k1)
    ops.uniaxialMaterial("Elastic", 2, k2)
    # NOTE: zeroLength defaults to -doRayleigh 0, which SILENTLY drops
    # stiffness-proportional (betaK) Rayleigh damping. Must set -doRayleigh 1
    # or betaK has zero effect (and d(response)/d(betaK) is honestly 0).
    ops.element("zeroLength", 1, 0, 1, "-mat", 1, "-dir", 1, "-doRayleigh", 1)  # story 1
    ops.element("zeroLength", 2, 1, 2, "-mat", 2, "-dir", 1, "-doRayleigh", 1)  # story 2
    eig = ops.eigen("-fullGenLapack", 2)
    w1, w2 = math.sqrt(eig[0]), math.sqrt(eig[1])
    return w1, w2


def _rayleigh_from_xi(xi, w1, w2):
    """Classical Rayleigh coeffs giving damping ratio xi at modes 1 and 2."""
    alphaM = xi * 2.0 * w1 * w2 / (w1 + w2)
    betaK = xi * 2.0 / (w1 + w2)
    return alphaM, betaK


def _run(alphaM, betaK, w_drive):
    """Forward solve with explicit Rayleigh coeffs. Returns peak story-1 drift."""
    ops.rayleigh(alphaM, betaK, 0.0, 0.0)
    # harmonic ground acceleration near mode 1
    T = n_cycles * (2.0 * math.pi / w_drive)
    period = 2.0 * math.pi / w_drive
    ops.timeSeries("Trig", 1, 0.0, T, period, "-factor", 1.0)
    ops.pattern("UniformExcitation", 1, 1, "-accel", 1)

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1e-10, 25)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")

    nsteps = int(T / dt)
    peak = 0.0
    for _ in range(nsteps):
        ops.analyze(1, dt)
        d1 = abs(ops.nodeDisp(1, 1))          # story-1 drift (floor1 - ground)
        if d1 > peak:
            peak = d1
    return peak


def forward_xi(params):
    """forward(params=[xi]) -> peak story-1 drift. Rebuilds modes each call."""
    (xi,) = params
    w1, w2 = _build_frame()
    aM, bK = _rayleigh_from_xi(xi, w1, w2)
    return _run(aM, bK, w_drive=w1)


def forward_coeffs(params):
    """forward(params=[alphaM, betaK]) -> peak story-1 drift."""
    aM, bK = params
    w1, w2 = _build_frame()
    return _run(aM, bK, w_drive=w1)


if __name__ == "__main__":
    # baseline modal info
    w1, w2 = _build_frame()
    print("=" * 68)
    print(f" 2-DOF shear frame: w1={w1:.4f} rad/s ({w1/2/math.pi:.3f} Hz), "
          f"w2={w2:.4f} rad/s")
    print(f" driven by harmonic ground accel at mode 1; response = peak story-1 drift")
    print("=" * 68)

    # ---- Study A: d(peak drift)/d(xi) ---------------------------------
    xi0 = 0.05
    base = forward_xi([xi0])
    g_xi = fd_gradient(forward_xi, [xi0], rel_step=1e-2, scheme="central")[0]
    print(f"\n STUDY A  scalar damping-ratio parameter, xi0={xi0}")
    print(f"   peak story-1 drift at xi0 = {base:.6e}")
    print(f"   d(peak drift)/d(xi)       = {g_xi:.6e}   (negative: more damping -> less drift)")
    print(f"   step-size plateau (central FD):")
    print(f"     {'rel.step':>10} {'d(drift)/d(xi)':>18}")
    for r, g in step_study(forward_xi, [xi0], comp=0):
        print(f"     {r:>10.0e} {g:>18.6e}")

    # ---- Study B: gradient w.r.t. [alphaM, betaK] ---------------------
    aM0, bK0 = _rayleigh_from_xi(xi0, w1, w2)
    grad = fd_gradient(forward_coeffs, [aM0, bK0], rel_step=1e-2, scheme="central")
    print(f"\n STUDY B  raw Rayleigh coefficients, [alphaM, betaK]=[{aM0:.5f}, {bK0:.6f}]")
    print(f"   d(peak drift)/d(alphaM) = {grad[0]:.6e}")
    print(f"   d(peak drift)/d(betaK)  = {grad[1]:.6e}")
    print(f"   -> gradient VECTOR returned in one call (2 params -> 4 forward solves)")
    print("=" * 68)
