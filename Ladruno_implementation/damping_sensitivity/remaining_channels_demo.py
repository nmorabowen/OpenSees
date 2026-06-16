"""
Close the "all sources of damping" claim: FD sensitivity for the three damping
channels not yet demonstrated -- MODAL damping, a DAMPING OBJECT (UniformDamping),
and NODAL mass-proportional alphaM. Uses the FDSensitivity class. Zero vanilla edits.

All three are SDOF at resonance, so each has the SAME closed-form oracle as a
generic damping ratio xi (amplitude U = F0/(2*k*xi), dU/dxi = -F0/(2*k*xi^2)),
except nodal-alphaM which is parameterised directly by alphaM
(U = F0/(alphaM*m*wn), dU/dalphaM = -F0/(alphaM*m*wn^2)).
"""
import os, sys, math
from fd_sensitivity import FDSensitivity

DIST = r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\dist\bin"
os.environ["LADRUNO_OPENSEES_QUIET"] = "1"
os.add_dll_directory(DIST)
sys.path.insert(0, DIST)
import opensees as ops  # noqa: E402

m = 1.0
f_n = 1.0
wn = 2.0 * math.pi * f_n
k = wn * wn * m
F0 = 1.0
dt = 0.002
n_cycles = 120
T = n_cycles / f_n
nsteps = int(T / dt)


def _common_model():
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0); ops.fix(1, 1)
    ops.node(2, 0.0); ops.mass(2, m)
    ops.uniaxialMaterial("Elastic", 1, k)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1, "-doRayleigh", 1)
    ops.timeSeries("Trig", 1, 0.0, T, 1.0 / f_n, "-factor", F0)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0)


def _run_to_amplitude():
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("ProfileSPD")
    ops.test("NormDispIncr", 1e-10, 25)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    tail = T - 10.0 / f_n
    umax = 0.0
    for _ in range(nsteps):
        ops.analyze(1, dt)
        if ops.getTime() >= tail:
            umax = max(umax, abs(ops.nodeDisp(2, 1)))
    return umax


# ---- channel: MODAL damping -------------------------------------------
def forward_modal(params):
    (xi,) = params
    _common_model()
    ops.eigen("-fullGenLapack", 1)
    ops.modalDamping(xi)
    return _run_to_amplitude()


# ---- channel: DAMPING OBJECT (UniformDamping over a band around wn) ----
def forward_uniform(params):
    (xi,) = params
    _common_model()
    # damping('Uniform', tag, zeta, fLower, fUpper) ; band brackets wn (rad/s)
    ops.damping("Uniform", 1, xi, 0.5 * wn, 2.0 * wn)
    ops.region(1, "-ele", 1, "-damp", 1)
    return _run_to_amplitude()


# ---- channel: NODAL mass-proportional alphaM --------------------------
def forward_nodal_alphaM(params):
    (alphaM,) = params
    _common_model()
    ops.rayleigh(alphaM, 0.0, 0.0, 0.0)   # alphaM hits the nodal mass -> C = alphaM*m
    return _run_to_amplitude()


# Parameter-free oracle: at resonance the steady amplitude is U(p) = A/p for a
# constant A, because the damping coefficient scales linearly with the parameter
# p (modal xi, generic xi, nodal alphaM via c=alphaM*m, and a Damping object's
# input zeta -> realized xi_eff proportional to it). Hence EXACTLY
#     dU/dp = -A/p^2 = -U/p .
# This sidesteps per-channel closed forms (and the typo class of bug they invite)
# AND the fact that a Damping object's input zeta != the realized SDOF xi.
def proportional_oracle(U, p):
    return -U / p


def report(name, fwd, x0, note=""):
    fd = FDSensitivity(fwd, rel_step=3e-3, scheme="central")
    U_fd = fwd([x0])              # one direct solve for the amplitude value
    g = fd.gradient([x0])[0]
    g_oracle = proportional_oracle(U_fd, x0)
    xi_eff = F0 / (2.0 * k * U_fd)   # damping ratio actually realised at resonance
    print(f"\n--- {name} ---")
    if note:
        print(f"  {note}")
    print(f"  amplitude U(p)            = {U_fd:.6e}   (realised xi_eff = {xi_eff:.4f})")
    print(f"  gradient  FD              = {g:.6e}")
    print(f"  gradient  oracle (-U/p)   = {g_oracle:.6e}   (rel.err {abs(g-g_oracle)/abs(g_oracle):.2%})")
    print(f"  plateau:")
    for r, gg in fd.step_study([x0]):
        print(f"     {r:>9.0e}  {gg:>15.6e}")
    print(f"  forward solves used (cache on): {fd.n_solves}")


if __name__ == "__main__":
    print("=" * 68)
    print(f" Remaining damping channels, SDOF @ resonance (f_n={f_n} Hz, k={k:.4f})")
    print(" FD sensitivity via FDSensitivity class; validated vs oracle dU/dp=-U/p")
    print("=" * 68)

    xi0 = 0.05
    report("MODAL damping  d(U)/d(xi)", forward_modal, xi0)
    report("UniformDamping object  d(U)/d(xi)", forward_uniform, xi0,
           note="NB: object input zeta != realised xi_eff (band-fit operator); "
                "FD still exact via -U/p.")
    aM0 = 2.0 * xi0 * wn          # alphaM giving ~5% damping
    report("NODAL alphaM  d(U)/d(alphaM)", forward_nodal_alphaM, aM0)
    print("\n" + "=" * 68)
