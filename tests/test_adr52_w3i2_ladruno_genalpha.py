"""ADR-52 W3-I2 — LadrunoGeneralizedAlpha: sensitivity/DDM subclass of GeneralizedAlpha.

Strict superset of LadrunoHHT: the Chung-Hulbert generalized-alpha method carries TWO
spectral parameters — alphaF weights stiffness/damping (Ualpha / Ualphadot), alphaM
weights inertia (Ualphadotdot) — where HHT has a single alpha on K/C and full-step M.
Because U/Udot/Udotdot still follow the Newmark recurrence, the sensitivity propagation
(saveSensitivity / commitSensitivity / formSensitivityRHS / computeSensitivities) is
identical to Newmark; only formEleResidual/formNodUnbalance change:
  - the mass multiplicator picks up the alphaM-weighting of d(Ualphadotdot)/dh,
  - the damping multiplicator picks up the alphaF-weighting of d(Ualphadot)/dh,
  - dM/dh acts on Ualphadotdot, dC/dh on Ualphadot,
  - the extra element-only term is -K*(1-alphaF)*dU_n.

Decisive checks (alphaM != alphaF, both != 1 so every weighting is exercised):
  1. EXISTS + SENSITIVITY-OFF BYTE-IDENTITY vs vanilla GeneralizedAlpha.
  2. DDM vs CENTRAL FD of d(U)/dE, undamped: exercises the alphaF stiffness weighting,
     the new addK_Force term, and the alphaM-weighted mass multiplicator.
  3. DDM vs CENTRAL FD, mass-proportional Rayleigh (C = alphaM_ray*M != 0): exercises the
     alphaF-weighted damping multiplicator (addD_Force) too.
  4. REDUCES TO NEWMARK DDM AT alphaM == alphaF == 1: both spectral params at 1 collapse
     every weighting (beta=0.25, gamma=0.5) -> must equal Newmark's DDM.

Scope note: as for LadrunoHHT, the dM/dh-on-Ualphadotdot and dC/dh-on-Ualphadot terms can
only be FD-exercised with a nonzero dM/dh / dC/dh, which needs an element implementing
getMassSensitivity / getTangentStiffSensitivity. Truss implements neither (Element's base
returns zero), so those two terms are validated by derivation + the alphaM=alphaF=1 ->
Newmark reduction; the mass-proportional case still validates the alphaF-weighted damping
multiplicator. param=E (dM/dh=dC/dh=0) keeps the FD oracle clean.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# --- SDOF truss: fixed node 1 -- truss -- mass node 2 (x only) ---------------
E0, A, L, M = 2.0e3, 1.5, 1.0, 1.0
DT, NSTEPS = 5.0e-4, 80
ALPHA_M, ALPHA_F = 0.9, 0.8         # both != 1 and != each other -> two-alpha terms live


def _build(E, int_args, with_sens, alphaM_ray=0.0):
    """Build the SDOF truss, run a transient sine pulse; return (U2x history, sens)."""
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)                 # free in x only
    ops.mass(2, M, M)

    ops.uniaxialMaterial('Elastic', 1, E)
    ops.element('Truss', 1, 1, 2, A, 1)

    ops.timeSeries('Trig', 1, 0.0, DT * NSTEPS, DT * NSTEPS, '-factor', 100.0)
    ops.pattern('Plain', 1, 1)
    ops.load(2, 1.0, 0.0)

    if with_sens:
        ops.parameter(1, 'element', 1, 'E')

    if alphaM_ray > 0.0:
        ops.rayleigh(alphaM_ray, 0.0, 0.0, 0.0)   # C = alphaM_ray * M (E-independent)

    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', 1.0e-12, 50)
    ops.algorithm('Newton')
    ops.integrator(*int_args)
    ops.analysis('Transient')
    if with_sens:
        ops.sensitivityAlgorithm('-computeAtEachStep')

    hist = []
    for _ in range(NSTEPS):
        assert ops.analyze(1, DT) == 0
        hist.append(ops.nodeDisp(2, 1))

    sens = ops.sensNodeDisp(2, 1, 1) if with_sens else None
    return hist, sens


def test_ladruno_genalpha_sensitivity_off_byte_identical():
    """Plain LadrunoGeneralizedAlpha == vanilla GeneralizedAlpha, bit for bit."""
    ref, _ = _build(E0, ('GeneralizedAlpha', ALPHA_M, ALPHA_F), with_sens=False)
    lad, _ = _build(E0, ('LadrunoGeneralizedAlpha', ALPHA_M, ALPHA_F), with_sens=False)
    assert lad == ref, "LadrunoGeneralizedAlpha sensitivity-off path must be byte-identical"


def _fd_gradient(alphaM_ray):
    dE = E0 * 1.0e-6
    hp, _ = _build(E0 + dE, ('LadrunoGeneralizedAlpha', ALPHA_M, ALPHA_F), with_sens=False, alphaM_ray=alphaM_ray)
    hm, _ = _build(E0 - dE, ('LadrunoGeneralizedAlpha', ALPHA_M, ALPHA_F), with_sens=False, alphaM_ray=alphaM_ray)
    return (hp[-1] - hm[-1]) / (2.0 * dE)


@pytest.mark.parametrize("alphaM_ray", [0.0, 4.0], ids=["undamped", "massPropDamped"])
def test_ladruno_genalpha_ddm_matches_finite_difference(alphaM_ray):
    """DDM d(U2x)/dE from LadrunoGeneralizedAlpha matches a central FD of the response."""
    _, ddm = _build(E0, ('LadrunoGeneralizedAlpha', ALPHA_M, ALPHA_F), with_sens=True, alphaM_ray=alphaM_ray)
    fd = _fd_gradient(alphaM_ray)
    assert abs(fd) > 0.0
    rel_err = abs(ddm - fd) / abs(fd)
    assert rel_err < 1.0e-5, (
        f"DDM vs FD mismatch (alphaM_ray={alphaM_ray}): ddm={ddm:.6e} fd={fd:.6e} "
        f"rel_err={rel_err:.3e}"
    )


def test_ladruno_genalpha_reduces_to_newmark_ddm():
    """At alphaM == alphaF == 1 every weighting collapses (beta=0.25, gamma=0.5):
    LadrunoGeneralizedAlpha's DDM must equal Newmark's DDM on the same model."""
    _, ddm_lad = _build(E0, ('LadrunoGeneralizedAlpha', 1.0, 1.0), with_sens=True)
    _, ddm_nm = _build(E0, ('Newmark', 0.5, 0.25), with_sens=True)
    assert abs(ddm_nm) > 0.0
    rel_err = abs(ddm_lad - ddm_nm) / abs(ddm_nm)
    assert rel_err < 1.0e-8, (
        f"LadrunoGeneralizedAlpha(1,1) DDM must match Newmark DDM: lad={ddm_lad:.8e} "
        f"nm={ddm_nm:.8e} rel_err={rel_err:.3e}"
    )
