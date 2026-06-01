"""LadrunoBrick — cumulative VISCOUS-hourglass dissipation reporting.

`uri -hourglass viscous` controls the zero-energy (hourglass) modes with a
Flanagan-Belytschko rate-form damping force (Belytschko 8.7.10). That force
stores NO energy — it dissipates — so the stateless `hourglassEnergy()` cannot
report it instantaneously. Instead the element integrates the work done against
the damping force over committed steps into a committed accumulator, and reports
that running total through the same `"hourglassEnergy"` response (aliases
`"hgDissipation"` / `"hourglassDissipation"`):

    ΔE = c_visc · Σ_{a,i} q̇_aι·Δq_aι ,   q̇ from getTrialVel, Δq from Δu.

This is the GLSTAT-style "hourglass energy" LS-DYNA reports for viscous HG
control — a diagnostic that the spurious modes are being damped, and how much
energy that cost.

Test design — exciting a CLEAN hourglass mode. uri-viscous has no hourglass
stiffness, so its STATIC tangent is rank-deficient (explicit-only). To probe it
we pin only the x,y DOFs (leaving z free, so the z-hourglass mode is untouched)
and kick with the exact ξη hourglass pattern  v_z = (2x−1)(2y−1) = ±1 — which is
zero-strain at the 1-pt centroid (a genuine zero-energy mode) and carries zero
net z-momentum (no rigid drift). Under CentralDifferenceLadruno this mode has no
restoring force, only the viscous damper, so its kinetic energy is dissipated
cleanly and entirely.

Gates:
  1. the accumulator is monotonically non-decreasing and ends strictly positive;
  2. it is BOUNDED by the hourglass kinetic energy imparted (the only sink is the
     viscous damper, so it can never report more energy than was there) AND with
     light damping it recovers a substantial fraction of it;
  3. that energy-conservation bound holds across a RANGE of damping levels — the
     accumulator never fabricates energy regardless of eps;
  4. a rigid / constant-velocity field has zero hourglass content (γ ⟂ the linear
     field) -> exactly zero dissipation;
  5. the `"hgDissipation"` alias returns the same value as `"hourglassEnergy"`;
  6. only uri-viscous reports a growing dissipation — stiffness reports its
     (instantaneous, oscillating) stored energy, physical reports 0.

Note on the work integral: hgDissipated is the discrete work ∫f·du done against
the FB rate damper — a monotone, energy-bounded DIAGNOSTIC, not a machine-precise
energy balance. Under a leapfrog stepper it tracks the dissipated energy well for
light damping; for very strong damping the per-step velocity collapse makes f·Δu
under-count (and the explicit step itself needs a smaller dt to stay stable), so
the gates below assert the robust properties (monotone / positive / bounded),
not exact convergence.

numpy-free (math only), like the other LadrunoBrick Zone-A files.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_E, _NU, _RHO = 1000.0, 0.0, 1.0       # nu=0: clean shear modulus, no vol coupling
_C = math.sqrt(_E / _RHO)
_LE = 1.0

_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
_CONN = [1, 2, 3, 4, 5, 6, 7, 8]

# the ξη hourglass pattern in z:  h = (2x−1)(2y−1)  evaluated at each node.
# Zero-strain at the centroid (a true zero-energy mode), Σh = 0 (no z-momentum).
_HG_Z = {t: (2 * c[0] - 1) * (2 * c[1] - 1) for t, c in _CUBE.items()}


def _zfree_cube(eps=0.0, hg="viscous"):
    """A single unit hex with x,y pinned and z free (so the ξη z-hourglass mode is
    unconstrained), uri + the given hourglass flavour, lumped mass. eps is the
    hourglass coefficient (0 -> element default)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *c)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU, _RHO)
    for t in _CUBE:
        ops.fix(t, 1, 1, 0)                       # x,y pinned; z free
    args = ["-formulation", "uri", "-hourglass", hg]
    if eps > 0.0:
        args += [eps]
    args += ["-lumped"]
    ops.element("LadrunoBrick", 1, *_CONN, 1, *args)


def _hg(tag=1, name="hourglassEnergy"):
    ops.eleResponse(tag, "forces")               # ensure trial state is current
    return list(ops.eleResponse(tag, name))[0]


def _setup_transient():
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno")
    ops.analysis("Transient")


def _run_hg_kick(eps=0.0, vz=1.0, nsteps=800, safety=0.03, record=False):
    """Kick the z DOFs with the pure ξη hourglass pattern and integrate with
    CentralDifferenceLadruno. safety<<1 keeps the viscous explicit step stable at
    the eps levels used here. Returns the final accumulator, or its history."""
    _zfree_cube(eps)
    for t in _CUBE:
        ops.setNodeVel(t, 3, _HG_Z[t] * vz, "-commit")
    _setup_transient()
    dt = safety * _LE / _C
    hist = [_hg()]
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            break
        if record:
            hist.append(_hg())
    return hist if record else _hg()


def _hourglass_kinetic_energy(vz=1.0):
    """Lumped KE of the ξη kick: each node carries m = rho*V/8, v = ±vz."""
    m = _RHO * (_LE ** 3) / 8.0
    return sum(0.5 * m * (_HG_Z[t] * vz) ** 2 for t in _CUBE)


# --------------------------------------------------------------------------
# 1. monotonic accumulation + positivity (headline)
# --------------------------------------------------------------------------
def test_dissipation_monotonic_and_positive():
    """Under a hourglass-exciting explicit run, the viscous-hourglass dissipation
    accumulator only ever grows (work against a damper is non-negative) and ends
    strictly positive — the spurious mode is being bled off, and reported."""
    hist = _run_hg_kick(eps=0.05, record=True, nsteps=400)
    assert all(math.isfinite(v) for v in hist), "dissipation history went non-finite"
    for a, b in zip(hist, hist[1:]):
        assert b >= a - 1e-12, f"dissipation decreased: {a:.6e} -> {b:.6e}"
    assert hist[0] == 0.0, "accumulator should start at exactly 0"
    assert hist[-1] > 0.0, f"no viscous dissipation accumulated (end={hist[-1]:.3e})"


# --------------------------------------------------------------------------
# 2. bounded by the hourglass kinetic energy, and a substantial fraction of it
# --------------------------------------------------------------------------
def test_dissipation_bounded_and_substantial():
    """The only sink is the viscous damper and the kick is a pure zero-energy mode
    (no restoring force), so the work done against the damper can never exceed the
    hourglass kinetic energy imparted. With light damping (where the discrete
    work integral tracks the energy closely) it should also recover a substantial
    fraction of it over a long run."""
    E = _run_hg_kick(eps=0.05, nsteps=1500, safety=0.02)
    KE0 = _hourglass_kinetic_energy()
    assert E > 0.0, "expected positive dissipation"
    assert E <= KE0 * 1.01, f"dissipation {E:.4e} exceeds imparted HG energy {KE0:.4e}"
    assert E > 0.5 * KE0, (
        f"light-damping dissipation {E:.4e} should recover a substantial fraction "
        f"of the imparted HG energy {KE0:.4e} (got {E / KE0:.3f})"
    )


# --------------------------------------------------------------------------
# 3. energy is never fabricated — the bound holds across damping levels
# --------------------------------------------------------------------------
@pytest.mark.parametrize("eps", [0.02, 0.05, 0.1])
def test_dissipation_never_exceeds_energy(eps):
    """Whatever the damping level, the reported dissipation is a positive quantity
    bounded above by the kinetic energy imparted — the work accumulator never
    fabricates energy (it integrates work against a real, energy-removing force)."""
    E = _run_hg_kick(eps=eps, nsteps=600, safety=0.02)
    KE0 = _hourglass_kinetic_energy()
    assert 0.0 < E <= KE0 * 1.01, (
        f"eps={eps}: dissipation {E:.4e} must be positive and <= imparted "
        f"HG energy {KE0:.4e}"
    )


# --------------------------------------------------------------------------
# 4. rigid / constant velocity has no hourglass content -> zero dissipation
# --------------------------------------------------------------------------
def test_rigid_velocity_no_dissipation():
    """A uniform z-velocity is a rigid translation: γ ⟂ the constant field, so the
    generalized hourglass velocity q̇ = Σ γ·v is exactly zero and nothing is
    dissipated, even though the body is moving."""
    _zfree_cube()
    for t in _CUBE:
        ops.setNodeVel(t, 3, 1.0, "-commit")     # uniform -> pure rigid translation
    _setup_transient()
    dt = 0.1 * _LE / _C
    for _ in range(100):
        assert ops.analyze(1, dt) == 0
    assert _hg() == 0.0, "rigid translation must dissipate exactly 0 (γ orthogonality)"


# --------------------------------------------------------------------------
# 5. the dissipation alias mirrors hourglassEnergy
# --------------------------------------------------------------------------
def test_dissipation_alias_matches():
    """`hgDissipation` / `hourglassDissipation` are explicit aliases for the same
    id-8 response; for the viscous flavour all return the cumulative dissipation."""
    _zfree_cube()
    for t in _CUBE:
        ops.setNodeVel(t, 3, _HG_Z[t] * 1.0, "-commit")
    _setup_transient()
    dt = 0.1 * _LE / _C
    for _ in range(50):
        assert ops.analyze(1, dt) == 0
    e1 = _hg(1, "hourglassEnergy")
    e2 = _hg(1, "hgDissipation")
    e3 = _hg(1, "hourglassDissipation")
    assert e1 > 0.0 and e1 == e2 == e3, f"alias mismatch: {e1}, {e2}, {e3}"


# --------------------------------------------------------------------------
# 6. only viscous reports a growing accumulator (others store / report 0)
# --------------------------------------------------------------------------
@pytest.mark.parametrize("hg,expect", [
    ("viscous", "grows"),       # dissipation accumulator > 0 and non-decreasing
    ("stiffness", "stored"),    # instantaneous stored energy, oscillating, > 0 somewhere
    ("physical", "zero"),       # assumed strain -> no separable hourglass energy
])
def test_only_viscous_accumulates(hg, expect):
    """Same ξη hourglass kick under explicit dynamics: viscous accumulates a
    positive dissipation total; stiffness reports its (instantaneous, nonzero)
    stored hourglass energy; physical reports 0. Separates the dissipation report
    from the stored-energy report."""
    _zfree_cube(hg=hg)
    for t in _CUBE:
        ops.setNodeVel(t, 3, _HG_Z[t] * 1.0, "-commit")
    _setup_transient()
    dt = 0.1 * _LE / _C
    early = None
    emax = 0.0
    for s in range(120):
        assert ops.analyze(1, dt) == 0
        e = _hg()
        emax = max(emax, e)
        if s == 5:
            early = e
    end = _hg()

    if expect == "grows":
        assert end > 0.0 and end >= early, f"viscous should accumulate (got {early}->{end})"
    elif expect == "stored":
        assert emax > 0.0, f"stiffness should report stored hourglass energy (max={emax})"
    else:  # zero
        assert emax == 0.0, f"physical should report 0 hourglass energy (max={emax})"
