"""LadrunoDynamicRelaxation (Ladruno quasi-static DR integrator, classTag 33005)
— Zone-A battery.

Quasi-static dynamic relaxation: explicit matrix-free leap-frog with a Gershgorin
fictitious mass + Cundall kinetic damping, driven to STATIC equilibrium under the
stock transient driver (`algorithm Linear` + `system Diagonal`). The von Mises
snap-through truss (PR-1 / tests/_proto_arch_snapthrough.py; limit load ≈3.80) is
the shared oracle.

  DR-1  sub-critical relax → same equilibrium `uy` stock LoadControl reaches.
  DR-6  snap-through — super-critical constant load relaxes to the FAR stable
        branch (`uy ≈ −0.217`, the AL-6 oracle), where LoadControl fails.
  DR-8  M* independence — the (massless) truss still relaxes via the Gershgorin
        stiffness probe (NOT getMass, which is ~0 here — the ADR-20 BLOCKER).
  DR-S  parser smoke for the documented argument forms.

Theory / design: Ladruno_implementation/21_ladruno_dynamic_relaxation_adr.md.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

DR = "LadrunoDynamicRelaxation"

_H = 0.10
_SPAN = 1.0
_E = 1.0e4
_A = 1.0
_APEX = 2
_LIMIT = 3.80          # von Mises snap-through limit load (PR-1)


def _build_geom():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, -_SPAN, 0.0)
    ops.node(_APEX, 0.0, _H)
    ops.node(3, _SPAN, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(3, 1, 1)
    ops.fix(_APEX, 1, 0)
    ops.uniaxialMaterial("Elastic", 1, _E)
    ops.element("corotTruss", 1, 1, _APEX, _A, 1)
    ops.element("corotTruss", 2, _APEX, 3, _A, 1)


def _loadcontrol_uy(target_lambda, dlam=0.05):
    """Stable rising-branch equilibrium uy at a sub-critical load, via LoadControl."""
    _build_geom()
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(_APEX, 0.0, -1.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1e-10, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", dlam)
    ops.analysis("Static")
    while ops.getTime() < target_lambda - 1e-9:
        if ops.analyze(1) != 0:
            break
    return ops.nodeDisp(_APEX, 2)


def _relax(P, integrator_args, chunk=50, dt=1.0, max_iter=400000,
           tol_vel=1e-6, tol_du=1e-8):
    """Build the truss under a CONSTANT downward load P and relax to rest with the
    DR integrator. Returns (converged, uy, steps)."""
    _build_geom()
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(_APEX, 0.0, -P)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormUnbalance", 1e-6, 1, 0)
    ops.algorithm("Linear")
    ops.integrator(*integrator_args)
    ops.analysis("Transient")

    prev = ops.nodeDisp(_APEX, 2)
    done = 0
    while done < max_iter:
        if ops.analyze(chunk, dt) != 0:
            return False, ops.nodeDisp(_APEX, 2), done
        done += chunk
        vmax = abs(ops.nodeVel(_APEX, 2))
        u = ops.nodeDisp(_APEX, 2)
        dumax = abs(u - prev)
        prev = u
        if vmax < tol_vel and dumax < tol_du:
            return True, u, done
    return False, prev, done


# --------------------------------------------------------------------------
# DR-1 : sub-critical relax reaches the LoadControl equilibrium
# --------------------------------------------------------------------------
def test_dr1_subcritical_matches_loadcontrol():
    P = 3.0                      # < limit 3.80 -> rising-branch stable equilibrium
    ref_uy = _loadcontrol_uy(P)
    conv, uy, steps = _relax(P, (DR,))
    assert conv, f"DR did not relax to rest (steps={steps}, uy={uy})"
    assert abs(uy - ref_uy) < 1e-3, (
        f"DR rest uy={uy:.6f} != LoadControl uy={ref_uy:.6f}"
    )


# --------------------------------------------------------------------------
# DR-6 : snap-through to the far stable branch (LoadControl fails here)
# --------------------------------------------------------------------------
def test_dr6_snap_through_far_branch():
    P = 1.15 * _LIMIT            # super-critical -> must snap through
    # the far stable branch is MUCH stiffer than the near-flat start config, so a
    # M* frozen at the start would violate omega*dt<2 there -> -recompute refreshes
    # the Gershgorin mass as the truss stiffens (the snap-through analogue of the
    # softening-refresh the ADR notes).
    conv, uy, steps = _relax(P, (DR, "-recompute", 50), max_iter=800000)
    assert conv, f"DR did not settle after snap-through (steps={steps}, uy={uy})"
    # far stable branch (the AL-6 prototype oracle): uy well below the far flat -2H
    assert uy < -2.0 * _H, f"DR did not snap through: uy={uy:.4f} (expected < {-2*_H})"
    assert abs(uy - (-0.217)) < 0.03, f"far-branch uy={uy:.4f} off the AL-6 oracle (-0.217)"


# --------------------------------------------------------------------------
# DR-8 : M* independence — massless truss still relaxes (Gershgorin, not getMass)
# --------------------------------------------------------------------------
def test_dr8_massless_relaxes_via_gershgorin():
    # the truss carries NO mass (Elastic, no -rho) so getMass ~ 0; gershgorin builds
    # M* from the stiffness probe, so DR must still converge.
    P = 2.5
    conv, uy, steps = _relax(P, (DR, "-mass", "gershgorin"))
    assert conv, f"massless DR (gershgorin) failed to relax (steps={steps})"
    assert math.isfinite(uy) and uy < 0.0


# --------------------------------------------------------------------------
# DR-S : parser smoke
# --------------------------------------------------------------------------
@pytest.mark.parametrize("args", [
    (DR,),
    (DR, "-mass", "gershgorin"),
    (DR, "-mass", "unity"),
    (DR, "-mass", "lumped", 1.0),
    (DR, "-dt", 1.0, "-recompute", 100),
])
def test_drs_parser_smoke(args):
    _build_geom()
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(_APEX, 0.0, -1.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormUnbalance", 1e-6, 1, 0)
    ops.algorithm("Linear")
    ops.integrator(*args)
    ops.analysis("Transient")
    assert ops.analyze(5, 1.0) == 0, f"first steps failed for {args}"
    assert math.isfinite(ops.nodeDisp(_APEX, 2))
