"""WP-199 regression: late-time stability of a free-surface half-space over LadrunoCPML.

A 2-D plane-strain half-space (E = 2e8, nu = 0.25, rho = 2000; the apeGmsh
example-04 soil) on a 1.5 m Q4 grid: FourNodeQuad interior, a 9 m LadrunoCPML
layer on the base and on one side, free surface at y = 0 (unstretched, box top
at +1e6), half-symmetry at x = 0. A 10 Hz Ricker vertical surface load, Newmark
average acceleration, dt = 1.5 ms, 2000 steps (3 s). The element flags match the
example (`-m 2 -R0 1e-4 -f0 10 -lumped`).

Before the fix, with `-lumped` the layer's acceleration term used the HRZ mass
while its velocity / displacement-inertial memory terms used the consistent
N_a N_b coupling. The semi-discrete layer was then not a stretched copy of the
lumped lattice, and that lattice's grid modes that are backward in the layer
normal (u_x, odd-even in x near the tangential Nyquist) were amplified by the
PML. Interior kinetic energy fell to 3e-4 of its peak by 0.6 s, then grew about
x4 per 0.075 s at a rate independent of dt (semi-discrete, not a Newmark issue);
nodal displacements reached 6e13 m at 3 s. With one lumping rule for the whole
inertial operator the energy decays monotonically to ~1e-10 of the peak.

The existing CPML anchors did not catch it: the Lamb / Green anchors run
165-190 steps, and the 2500-step Lyapunov case is a 2x2-element toy with no
free surface. Wall time: ~5 s per case on the dev box.
"""
import math

import numpy as np
import pytest

import opensees as ops

E, NU, RHO = 2.0e8, 0.25, 2000.0
DT, N_STEPS, F0, T0 = 1.5e-3, 2000, 10.0, 0.12


def _run_halfspace(h: float = 1.5, lumped: bool = True, extra=()) -> tuple[np.ndarray, np.ndarray, float]:
    """Return (t, interior kinetic energy) sampled every 25 steps, plus max |u| at the end."""
    Lx, Ly, Lp = 30.0, 21.0, 9.0
    nx, ny = int(round((Lx + Lp) / h)), int(round((Ly + Lp) / h))
    nxi, nyi = int(round(Lx / h)), int(round(Ly / h))

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU, RHO)

    def nid(i: int, j: int) -> int:
        return 1 + j * (nx + 1) + i

    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.node(nid(i, j), i * h, -j * h)
            fx = 1 if (i == 0 or i == nx or j == ny) else 0
            fy = 1 if (i == nx or j == ny) else 0
            if fx or fy:
                ops.fix(nid(i, j), fx, fy)

    cpml_flags = ["-box", -1.0e6, Lx, -Ly, 1.0e6, "-L", Lp, "-m", 2.0,
                  "-R0", 1.0e-4, "-f0", F0, *(["-lumped"] if lumped else []), *extra]
    etag = 1
    for j in range(ny):
        for i in range(nx):
            n = (nid(i, j + 1), nid(i + 1, j + 1), nid(i + 1, j), nid(i, j))
            if i >= nxi or j >= nyi:
                ops.element("LadrunoCPML", etag, *n, "-E", E, "-nu", NU, "-rho", RHO, *cpml_flags)
            else:
                ops.element("quad", etag, *n, 1.0, "PlaneStrain", 1, 0.0, RHO, 0.0, 0.0)
            etag += 1

    vals = []
    for k in range(N_STEPS + 1):
        a = (math.pi * F0 * (k * DT - T0)) ** 2
        vals.append(-1.0e5 * (1.0 - 2.0 * a) * math.exp(-a))
    ops.timeSeries("Path", 1, "-dt", DT, "-values", *vals)
    ops.pattern("Plain", 1, 1)
    ops.load(nid(0, 0), 0.0, 1.0)

    ops.system("UmfPack")
    ops.numberer("RCM")
    ops.constraints("Plain")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.algorithm("Linear", "-FactorOnce")
    ops.analysis("Transient")

    interior = [nid(i, j) for j in range(nyi + 1) for i in range(nxi + 1)]
    m_node = RHO * h * h
    t_hist, ke_hist = [], []
    for k in range(1, N_STEPS + 1):
        assert ops.analyze(1, DT) == 0, f"analysis failed at step {k}"
        if k % 25 == 0:
            ke = 0.5 * m_node * sum(ops.nodeVel(n, 1) ** 2 + ops.nodeVel(n, 2) ** 2 for n in interior)
            t_hist.append(k * DT)
            ke_hist.append(ke)
    u_max = max(abs(ops.nodeDisp(n, d)) for n in ops.getNodeTags() for d in (1, 2))
    ops.wipe()
    return np.array(t_hist), np.array(ke_hist), u_max


@pytest.mark.parametrize(
    "label, lumped, extra",
    [
        ("lumped_example04", True, ()),         # the apeGmsh example-04 configuration
        ("lumped_kappa1", True, ("-kappaMax", 1.0)),  # fastest pre-fix growth (~x10 / 0.04 s)
        ("consistent", False, ()),               # was already stable; keep it so
    ],
)
def test_cpml_halfspace_energy_never_grows_late(label, lumped, extra):
    t, ke, u_max = _run_halfspace(lumped=lumped, extra=tuple(extra))
    peak = ke[t <= 0.5].max()  # source-driven peak (a growing mode cannot hide it)
    rel = ke / peak

    early = rel[(t >= 0.75) & (t < 1.5)].max()
    late = rel[t >= 1.5].max()
    final = rel[-1]
    print(f"[{label}] KE/peak: max(0.75-1.5 s)={early:.3e} max(1.5-3 s)={late:.3e} "
          f"final={final:.3e} u_max={u_max:.3e}")

    # Pre-fix (lumped): KE/peak = 1 (the growing mode dominates) from ~1.1 s on,
    # u_max ~ 6e13 m. Post-fix: monotone decay to ~1e-10.
    assert np.all(np.isfinite(ke))
    assert early < 1.0e-4, f"{label}: energy not absorbed by 0.75-1.5 s ({early:.3e})"
    assert late < 1.0e-6, f"{label}: late-time energy {late:.3e} of peak (growing mode?)"
    assert late <= early, f"{label}: energy grew between windows ({early:.3e} -> {late:.3e})"
    assert final < 1.0e-7, f"{label}: final energy {final:.3e} of peak"
    assert u_max < 1.0e-6, f"{label}: residual displacement {u_max:.3e} m"
