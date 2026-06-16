"""
FD sensitivity w.r.t. a VISCOUS DAMPER coefficient -- a damping channel that
has NO DDM derivative anywhere in OpenSees (ViscousMaterial exposes no
*Sensitivity methods). Black-box FD differentiates it anyway, with zero
vanilla edits.

Self-validation: SDOF spring + linear dashpot, driven at resonance.
    m u'' + C u' + k u = F0 sin(wn t)
    steady amplitude     U(C)  = F0 / (C * wn)
    exact sensitivity    dU/dC = -F0 / (C^2 * wn)
"""
import os, sys, math
from fd_sensitivity import fd_gradient, step_study

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


def forward(params):
    """forward(params=[C]) -> steady resonant amplitude.

    C is the linear viscous dashpot coefficient (Viscous material, alpha=1),
    in PARALLEL with the elastic spring via two zeroLength materials on dir 1.
    No Rayleigh damping at all -- damping is the material's own rate-dependent
    force, which has no DDM derivative. FD differentiates it regardless.
    """
    (C,) = params
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0); ops.fix(1, 1)
    ops.node(2, 0.0); ops.mass(2, m)
    ops.uniaxialMaterial("Elastic", 1, k)
    ops.uniaxialMaterial("Viscous", 2, C, 1.0)          # linear dashpot: f = C*v
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, 2, "-dir", 1, 1)  # spring + dashpot

    ops.timeSeries("Trig", 1, 0.0, T, 1.0 / f_n, "-factor", F0)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0)

    ops.constraints("Plain")
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


def analytic_U(C):    return F0 / (C * wn)
def analytic_dUdC(C): return -F0 / (C * C * wn)


if __name__ == "__main__":
    xi0 = 0.08
    C0 = 2.0 * xi0 * m * wn          # dashpot coeff giving ~8% damping
    print("=" * 64)
    print(f" SDOF spring+dashpot at resonance. C0={C0:.5f} (xi~{xi0})")
    print(f" Viscous material -> NO DDM derivative exists; FD handles it.")
    print("=" * 64)
    U_fd, U_an = forward([C0]), analytic_U(C0)
    print(f" amplitude:  FD={U_fd:.6e}  analytic={U_an:.6e}  "
          f"(rel.err {abs(U_fd-U_an)/U_an:.2%})")
    print(f"\n dU/dC : central FD vs analytic -F0/(C^2 wn) = {analytic_dUdC(C0):.6e}")
    print(f"   {'rel.step':>10} {'FD dU/dC':>16} {'rel.err':>10}")
    g_an = analytic_dUdC(C0)
    for r, g in step_study(forward, [C0]):
        print(f"   {r:>10.0e} {g:>16.6e} {abs(g-g_an)/abs(g_an):>9.2%}")
    print("=" * 64)
