"""
Finite-difference damping sensitivity WITHOUT touching vanilla OpenSees.

Proves: dResponse/d(damping parameter) via black-box finite differences
(perturb the parameter, re-run the forward solve, difference the response).
No C++ changes, no DDM derivative methods, no reliability module needed.

Self-validation: a resonant SDOF has a closed-form damping sensitivity.
  steady-state amplitude at resonance   U(xi) = F0 / (2 k xi)
  exact sensitivity                     dU/dxi = -F0 / (2 k xi^2)
If the FD gradient matches the analytic one, the black-box route is sound.
"""
import os, sys, math

DIST = r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\dist\bin"
os.environ["LADRUNO_OPENSEES_QUIET"] = "1"
os.add_dll_directory(DIST)
sys.path.insert(0, DIST)
import opensees as ops  # noqa: E402

# ---- SDOF at resonance -------------------------------------------------
m = 1.0
f_n = 1.0                      # Hz
wn = 2.0 * math.pi * f_n       # rad/s
k = wn * wn * m                # so that sqrt(k/m) = wn
F0 = 1.0                       # harmonic force amplitude, driven AT resonance

dt = 0.002
n_cycles = 120                 # long enough for the transient to die out
T = n_cycles / f_n
nsteps = int(T / dt)


def run(xi):
    """Forward solve. Returns steady-state amplitude (max|u| over last 10 cycles).

    Damping ratio xi is injected as MASS-proportional Rayleigh damping:
        C = alphaM * M,  alphaM = 2*xi*wn  ->  c = 2*xi*wn*m = 2*xi*m*wn  (exact SDOF c)
    xi is the 'parameter' we differentiate the response with respect to.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 0.0)
    ops.fix(1, 1)
    ops.uniaxialMaterial("Elastic", 1, k)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1)
    ops.mass(2, m)

    # harmonic load at the resonant frequency
    ops.timeSeries("Trig", 1, 0.0, T, 1.0 / f_n, "-factor", F0)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0)

    ops.rayleigh(2.0 * xi * wn, 0.0, 0.0, 0.0)   # alphaM only

    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("ProfileSPD")
    ops.test("NormDispIncr", 1e-10, 20)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")

    tail_start = T - 10.0 / f_n
    umax = 0.0
    for _ in range(nsteps):
        ops.analyze(1, dt)
        t = ops.getTime()
        if t >= tail_start:
            u = abs(ops.nodeDisp(2, 1))
            if u > umax:
                umax = u
    return umax


def analytic_U(xi):
    return F0 / (2.0 * k * xi)


def analytic_dUdxi(xi):
    return -F0 / (2.0 * k * xi * xi)


# ---- Finite-difference gradient at xi0 ---------------------------------
xi0 = 0.10

print(f"{'='*64}")
print(f" Resonant SDOF: f_n={f_n} Hz, k={k:.4f}, F0={F0}, xi0={xi0}")
print(f" steady-state amplitude (FD vs analytic):")
U_fd = run(xi0)
U_an = analytic_U(xi0)
print(f"   U_fd      = {U_fd:.6e}")
print(f"   U_analytic= {U_an:.6e}   (rel.err {abs(U_fd-U_an)/U_an:.2%})")
print(f"{'='*64}")
print(f" dU/dxi : central finite difference vs analytic -F0/(2 k xi^2)")
print(f" analytic dU/dxi = {analytic_dUdxi(xi0):.6e}\n")
print(f"   {'rel.step':>10} {'FD dU/dxi':>16} {'rel.err':>12}")

for rel in (1e-1, 3e-2, 1e-2, 3e-3, 1e-3):
    h = rel * xi0
    g_fd = (run(xi0 + h) - run(xi0 - h)) / (2.0 * h)
    g_an = analytic_dUdxi(xi0)
    print(f"   {rel:>10.0e} {g_fd:>16.6e} {abs(g_fd-g_an)/abs(g_an):>11.2%}")

print(f"\n Negative sign = more damping reduces resonant amplitude. As expected.")
print(f"{'='*64}")
