"""LadrunoBrick — EXPLICIT-DYNAMICS validation under BOTH fork integrators.

`uri -lumped` is the cheap explicit hex LadrunoBrick adds. This file validates it
dynamically under the two explicit time steppers the fork ships:

    * CentralDifferenceLadruno  (leap-frog central difference, non-dissipative)
    * ExplicitBathe             (Noh-Bathe composite, sub-step p=0.54)

The load-bearing idea: if the element's mass matrix, internal force, and
critical-time-step contribution are all consistent, **two independent explicit
integrators must agree on the element's dynamics**. That cross-check is stronger
than either integrator alone, and it is exactly the gap the static tests cannot
see (mass + dt_cr never enter a static solve).

Gates:
  1. free axial vibration period agrees across CDL and ExplicitBathe, and both
     land near the continuum fixed-free value T1 = 4L/c, c = sqrt(E/rho);
  2. ops.criticalTimeStep() is positive, ~ le/c, and scales correctly with rho;
  3. the response stays bounded with no zero-energy growth over a long run;
  4. a transverse (hourglass-exciting) transient is controlled dynamically by the
     FB stiffness hourglass under both integrators.

numpy-free (math only), like the other LadrunoBrick Zone-A files.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

CDL = ("CentralDifferenceLadruno",)
EB = ("ExplicitBathe", 0.54)          # p = 0.54, the source-battery sub-step value

_E, _NU, _RHO = 1000.0, 0.0, 1.0      # nu=0 -> clean 1-D axial speed c = sqrt(E/rho)
_C = math.sqrt(_E / _RHO)             # = sqrt(M/rho) with M=E at nu=0
_LE = 1.0                             # unit hexes


def _nid(i, j, k, nx):
    return 1 + i + (nx + 1) * j + (nx + 1) * 2 * k


def _build_column(nx, formulation="uri", hg="stiffness", lateral_fixed=True,
                  lumped=True):
    """A column of nx unit hexes along x (1x1 cross-section), x=0 face clamped.
    lateral_fixed: pin all y,z DOFs (pure 1-D axial). Reuses the proven valid-hex
    connectivity from the bending file."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for k in range(2):
        for j in range(2):
            for i in range(nx + 1):
                ops.node(_nid(i, j, k, nx), i * _LE, j * _LE, k * _LE)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU, _RHO)

    args = ["-formulation", formulation]
    if formulation == "uri":
        args += ["-hourglass", hg]
    if lumped:
        args += ["-lumped"]

    e = 1
    for i in range(nx):
        ops.element("LadrunoBrick", e,
                    _nid(i, 0, 0, nx), _nid(i + 1, 0, 0, nx), _nid(i + 1, 1, 0, nx), _nid(i, 1, 0, nx),
                    _nid(i, 0, 1, nx), _nid(i + 1, 0, 1, nx), _nid(i + 1, 1, 1, nx), _nid(i, 1, 1, nx),
                    1, *args)
        e += 1

    for k in range(2):
        for j in range(2):
            ops.fix(_nid(0, j, k, nx), 1, 1, 1)            # clamp x=0 face
    if lateral_fixed:
        for k in range(2):
            for j in range(2):
                for i in range(1, nx + 1):
                    ops.fix(_nid(i, j, k, nx), 0, 1, 1)    # x free, y,z pinned


def _free_axial(integrator_args, nx=8, A=0.05, periods=2.2, safety=0.1):
    """Pull the bar into a linear axial ramp (excites the fundamental), release,
    free-vibrate. Returns (t[], tip_x[], T1_analytic)."""
    _build_column(nx, formulation="uri", hg="stiffness", lateral_fixed=True)
    L = nx * _LE
    # initial displacement u_x(x) = A x / L on the free (x) DOFs
    for k in range(2):
        for j in range(2):
            for i in range(1, nx + 1):
                ops.setNodeDisp(_nid(i, j, k, nx), 1, A * (i * _LE) / L, "-commit")

    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(*integrator_args)
    ops.analysis("Transient")

    dt = safety * _LE / _C
    T1 = 4.0 * L / _C
    nsteps = int(periods * T1 / dt)
    tip = _nid(nx, 0, 0, nx)
    t, u = [ops.getTime()], [ops.nodeDisp(tip, 1)]
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            break
        t.append(ops.getTime())
        u.append(ops.nodeDisp(tip, 1))
    return t, u, T1


def _period(t, u):
    """Period from successive upward zero-crossings of u about its mean."""
    um = sum(u) / len(u)
    c = [x - um for x in u]
    cross = []
    for i in range(1, len(c)):
        if c[i - 1] <= 0.0 < c[i]:
            frac = -c[i - 1] / (c[i] - c[i - 1])
            cross.append(t[i - 1] + frac * (t[i] - t[i - 1]))
    if len(cross) < 2:
        return None
    diffs = [cross[k + 1] - cross[k] for k in range(len(cross) - 1)]
    return sum(diffs) / len(diffs)


def _finite(xs):
    return all(math.isfinite(v) for v in xs)


def _maxabs(xs):
    return max(abs(v) for v in xs)


# --------------------------------------------------------------------------
# 1. cross-integrator agreement (the headline dynamic gate)
# --------------------------------------------------------------------------
def test_axial_period_agrees_across_integrators():
    """Free axial vibration period must agree between CentralDifferenceLadruno and
    ExplicitBathe (same discrete M and K → same dynamics, up to time-integration
    dispersion), and both must land near the continuum T1 = 4L/c."""
    t_cdl, u_cdl, T1 = _free_axial(CDL)
    t_eb, u_eb, _ = _free_axial(EB)
    assert _finite(u_cdl) and _finite(u_eb), "explicit axial run went non-finite"

    P_cdl = _period(t_cdl, u_cdl)
    P_eb = _period(t_eb, u_eb)
    assert P_cdl and P_eb, "could not extract a vibration period (too few crossings)"

    # the two integrators must agree closely on the period
    assert abs(P_cdl - P_eb) / P_cdl < 0.08, (
        f"period mismatch across integrators: CDL={P_cdl:.4f}, EB={P_eb:.4f}"
    )
    # ... and both near the continuum fundamental (coarse-mesh discretisation error)
    for nm, P in (("CDL", P_cdl), ("EB", P_eb)):
        assert abs(P - T1) / T1 < 0.20, f"{nm} period {P:.4f} vs analytic T1={T1:.4f}"


# --------------------------------------------------------------------------
# 2. critical time step queryable + physically scaled
# --------------------------------------------------------------------------
def _dtcr(rho):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    cube = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
            5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
    for tg, cc in cube.items():
        ops.node(tg, *cc)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU, rho)
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1)
    ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1,
                "-formulation", "uri", "-hourglass", "stiffness", "-lumped")
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(7, 0.0, 0.0, 1.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator("ExplicitBathe", 0.54, "-cfl")
    ops.analysis("Transient")
    ops.analyze(1, 1.0e-6)                  # priming step triggers the dt_cr compute
    return ops.criticalTimeStep()


def test_critical_timestep_sane_and_scales():
    """The element must feed the CFL machinery a sane critical step: positive, of
    order le/c, and scaling like sqrt(rho) (4x density -> ~2x dt_cr). This is the
    only gate that exercises the element MASS matrix + dt_cr contribution."""
    dt1 = _dtcr(_RHO)
    dt4 = _dtcr(4.0 * _RHO)
    assert dt1 > 0.0 and math.isfinite(dt1), f"dt_cr not positive/finite ({dt1})"
    # order of magnitude: a unit hex with c=sqrt(E/rho) has dt_cr ~ le/c
    le_over_c = _LE / _C
    assert 0.05 * le_over_c < dt1 < 3.0 * le_over_c, f"dt_cr={dt1:.3e} not ~le/c={le_over_c:.3e}"
    # density scaling: dt_cr ∝ sqrt(rho), so 4x rho -> ~2x dt_cr
    assert 1.8 < dt4 / dt1 < 2.2, f"dt_cr should ~double for 4x rho (got {dt4/dt1:.3f})"


# --------------------------------------------------------------------------
# 3. long-run stability / no zero-energy growth
# --------------------------------------------------------------------------
@pytest.mark.parametrize("integ", [CDL, EB], ids=["CDL", "ExplicitBathe"])
def test_axial_no_growth(integ):
    """Undamped free vibration must not grow: a stable, energy-respecting explicit
    run keeps the late-window amplitude near the early-window amplitude. Growth
    would flag an instability or a broken hourglass restoring force."""
    t, u, _ = _free_axial(integ, periods=4.0)
    assert _finite(u), "run went non-finite"
    n = len(u)
    early = _maxabs(u[: n // 5])
    late = _maxabs(u[-n // 5:])
    assert late <= 1.5 * early, f"amplitude grew (early={early:.3e}, late={late:.3e})"


# --------------------------------------------------------------------------
# 4. transverse hourglass controlled under explicit dynamics
# --------------------------------------------------------------------------
@pytest.mark.parametrize("integ", [CDL, EB], ids=["CDL", "ExplicitBathe"])
def test_transverse_hourglass_controlled(integ):
    """A transverse tip kick on a one-element-deep uri bar lives in the hourglass
    (bending) space. With FB `stiffness` hourglass control it must stay bounded
    under both explicit integrators — i.e. the zero-energy modes are stabilised
    dynamically, not just statically."""
    nx = 6
    _build_column(nx, formulation="uri", hg="stiffness", lateral_fixed=False)
    # initial transverse (z) velocity on the free-end face
    tip_face = [_nid(nx, j, k, nx) for k in range(2) for j in range(2)]
    for n in tip_face:
        ops.setNodeVel(n, 3, 1.0, "-commit")

    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(*integ)
    ops.analysis("Transient")

    dt = 0.1 * _LE / _C
    maxd = 0.0
    for _ in range(500):
        assert ops.analyze(1, dt) == 0
        d = ops.nodeDisp(tip_face[0], 3)
        assert math.isfinite(d), "transverse explicit run went non-finite"
        maxd = max(maxd, abs(d))
    assert maxd < 1.0e2, f"transverse hourglass response unbounded (max={maxd:.3e})"
