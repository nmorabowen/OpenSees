"""ADR-52 W2-E1 — explicit bulk viscosity on LadrunoBrick (runtime validation).

A single 8-node LadrunoBrick (unit cube), bottom face fixed, top face free in z,
x/y fixed everywhere (uniaxial strain -> nonzero volumetric strain rate, so bulk
viscosity engages). Explicit central difference (CentralDifferenceLadruno), kicked
with an initial top-face velocity so it oscillates axially.

Two decisive checks (the adversarial-review "decisive validation"):
  1. OFF-PATH IDENTITY: `-bulkViscosity 0 0` is bit-identical to no option (the
     bvActive guard skips the block) -- the regression gate for the default path.
  2. DISSIPATION + SIGN: with b1>0 the oscillation amplitude DECAYS (energy removed)
     and never AMPLIFIES vs the undamped run (confirms the dissipative sign).
"""
import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

L, E, NU, RHO = 1.0, 1000.0, 0.0, 1.0   # nu=0 -> D00=E, c_d=sqrt(E/rho)=31.6
DT, NSTEPS, V0 = 1.0e-3, 1600, -1.0     # dt << dt_cr (~0.045); ~11 axial periods


def _run(bulk):
    """Build + integrate; return the top-node (5) z-displacement history.

    bulk = None -> no -bulkViscosity option; (b1,b2) -> with the option.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    pts = [(0, 0, 0), (L, 0, 0), (L, L, 0), (0, L, 0),
           (0, 0, L), (L, 0, L), (L, L, L), (0, L, L)]
    tags = list(range(1, 9))
    for t, (x, y, z) in zip(tags, pts):
        ops.node(t, float(x), float(y), float(z))
    for t in tags[:4]:
        ops.fix(t, 1, 1, 1)               # bottom face: fully fixed
    for t in tags[4:]:
        ops.fix(t, 1, 1, 0)               # top face: z free, x/y fixed (uniaxial)

    ops.nDMaterial("ElasticIsotropic", 1, E, NU, RHO)
    if bulk is None:
        ops.element("LadrunoBrick", 1, *tags, 1, "-lumped")
    else:
        ops.element("LadrunoBrick", 1, *tags, 1, "-lumped",
                    "-bulkViscosity", float(bulk[0]), float(bulk[1]))

    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10)
    ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno")
    ops.analysis("Transient")

    for t in tags[4:]:
        ops.setNodeVel(t, 3, V0, "-commit")   # kick the top face axially

    hist = np.empty(NSTEPS)
    for k in range(NSTEPS):
        assert ops.analyze(1, DT) == 0, f"explicit step {k} failed"
        hist[k] = ops.nodeDisp(5, 3)
    return hist


def test_bulk_viscosity_off_path_bit_identical():
    """`-bulkViscosity 0 0` must reproduce the no-option run exactly (guard off)."""
    base = _run(None)
    zero = _run((0.0, 0.0))
    assert np.array_equal(base, zero), \
        "b1=b2=0 changed the trajectory -- the off-default path is not bit-identical"


def test_bulk_viscosity_dissipates_and_sign_is_correct():
    """b1>0 decays the oscillation (energy removed) and never amplifies it."""
    undamped = _run(None)
    damped = _run((1.0, 1.2))   # strong b1 for an unambiguous decay signal

    q = NSTEPS // 4
    amp_early = np.max(np.abs(damped[:q]))
    amp_late = np.max(np.abs(damped[-q:]))
    amp_undamped_late = np.max(np.abs(undamped[-q:]))

    # the undamped central-difference run is ~energy-conserving (no decay)
    assert np.max(np.abs(undamped[-q:])) > 0.5 * np.max(np.abs(undamped[:q])), \
        "undamped baseline unexpectedly decayed -- check stability/setup"
    # self-decay: late amplitude well below early (energy is being removed)
    assert amp_late < 0.9 * amp_early, \
        f"no decay with bulk viscosity: early={amp_early:.4e} late={amp_late:.4e}"
    # dissipation, not injection: damped late amplitude below the undamped late
    assert amp_late < amp_undamped_late, \
        f"damped late {amp_late:.4e} not below undamped late {amp_undamped_late:.4e}"
    # dissipative SIGN: the damped run never amplifies beyond the undamped envelope
    assert np.max(np.abs(damped)) <= 1.05 * np.max(np.abs(undamped)), \
        "bulk viscosity AMPLIFIED the response -- wrong sign (anti-damping)"
