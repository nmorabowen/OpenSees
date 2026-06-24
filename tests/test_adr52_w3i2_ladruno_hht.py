"""ADR-52 W3-I2 — LadrunoHHT: sensitivity-carrying (DDM) subclass of HHT.

LadrunoHHT reuses the entire HHT algorithm and adds the Direct Differentiation
Method (DDM) seam that vanilla OpenSees ships only on Newmark / DisplacementControl
/ ArcLength / MinUnbalDispNorm — never on the numerically-damped HHT. Because HHT's
U/Udot/Udotdot evolve by the SAME Newmark relations, the sensitivity propagation is
identical to Newmark; only the dynamic *residual* is re-derived at the HHT
alpha-intermediate state (Ualpha / Ualphadot), which:
  - alpha-weights the damping multiplicator d(Ualphadot)/dh,
  - evaluates dC/dh on Ualphadot (not Udot),
  - adds an EXTRA element-only stiffness term  -K*(1-alpha)*dU_n  that vanishes for
    Newmark (alpha == 1).

Decisive checks (alpha != 1 so every HHT-specific term is exercised):
  1. EXISTS + SENSITIVITY-OFF BYTE-IDENTITY: a plain transient run with LadrunoHHT
     is bit-for-bit identical to vanilla HHT (the no-sensitivity path delegates to
     the base TransientIntegrator).
  2. DDM vs CENTRAL FINITE DIFFERENCE, undamped: d(U)/dE from the integrator matches
     a central FD of the response. Exercises the stiffness alpha-weighting, the new
     addK_Force term, and the (Newmark-identical) mass terms.
  3. DDM vs CENTRAL FD, mass-proportional Rayleigh damping (C = alphaM*M != 0):
     additionally exercises the alpha-weighted damping MULTIPLICATOR (the addD_Force
     term, tmp2) — the velocity-sensitivity flow through a nonzero C. alphaM is
     E-independent so dC/dh == 0, hence the *dC/dh-on-Ualphadot* term
     (addD_ForceSensitivity) is exercised but contributes zero here.
  4. REDUCES TO NEWMARK DDM AT alpha == 1: with alpha == 1 the HHT-specific (1-alpha)
     terms drop out, so LadrunoHHT's DDM must equal Newmark's DDM on the same model —
     a cross-check that the alpha-weighting is the *only* difference from Newmark.

Scope note (the one term NOT independently FD-validated): the dC/dh-on-Ualphadot
contribution (addD_ForceSensitivity, where Ualphadot — not Udot — is the velocity the
HHT residual evaluates C at) can only be exercised numerically with dC/dh != 0, which
requires an element implementing getTangentStiffSensitivity / getMassSensitivity. Truss
implements neither (Element's base returns zero and warns for betaK*K DDM), so a
stiffness-proportional Rayleigh FD case would false-fail for an element limitation, not
an integrator bug. That term is therefore validated by derivation + the alpha==1 ->
Newmark reduction (which pins the surrounding algebra); revisit with a DDM-mass/stiffness
-capable element if one is added to the fork.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# --- SDOF truss: fixed node 1 -- truss -- mass node 2 (x only) ---------------
E0, A, L, M = 2.0e3, 1.5, 1.0, 1.0   # k = EA/L = 3000, w = sqrt(k/m) ~ 54.8, T ~ 0.115
DT, NSTEPS = 5.0e-4, 80
ALPHA = 0.8                          # HHT numerical damping (alpha != 1 -> HHT terms live)


def _build(E, int_args, with_sens, alphaM=0.0):
    """Build the SDOF truss, run a transient sine pulse; return (U2x history, sens).

    int_args : full argument tuple passed to ops.integrator (name + numeric args).
    with_sens: if True, set a parameter on E + sensitivityAlgorithm and return the
               final DDM sensitivity d(U2x)/dE as the second item (else None).
    alphaM   : mass-proportional Rayleigh coefficient (C = alphaM*M, E-independent).
    """
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)                 # free in x only
    ops.mass(2, M, M)

    ops.uniaxialMaterial('Elastic', 1, E)
    ops.element('Truss', 1, 1, 2, A, 1)

    # one sine pulse in x at node 2 over the whole run
    ops.timeSeries('Trig', 1, 0.0, DT * NSTEPS, DT * NSTEPS, '-factor', 100.0)
    ops.pattern('Plain', 1, 1)
    ops.load(2, 1.0, 0.0)

    if with_sens:
        ops.parameter(1, 'element', 1, 'E')

    if alphaM > 0.0:
        ops.rayleigh(alphaM, 0.0, 0.0, 0.0)   # C = alphaM * M (E-independent)

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


def test_ladrunohht_sensitivity_off_byte_identical():
    """Plain LadrunoHHT == vanilla HHT, bit for bit (no-sensitivity path)."""
    ref, _ = _build(E0, ('HHT', ALPHA), with_sens=False)
    lad, _ = _build(E0, ('LadrunoHHT', ALPHA), with_sens=False)
    assert lad == ref, "LadrunoHHT sensitivity-off path must be byte-identical to HHT"


def _fd_gradient(alphaM):
    """Central finite difference of final U2x w.r.t. E (DDM-independent)."""
    dE = E0 * 1.0e-6
    hp, _ = _build(E0 + dE, ('LadrunoHHT', ALPHA), with_sens=False, alphaM=alphaM)
    hm, _ = _build(E0 - dE, ('LadrunoHHT', ALPHA), with_sens=False, alphaM=alphaM)
    return (hp[-1] - hm[-1]) / (2.0 * dE)


@pytest.mark.parametrize("alphaM", [0.0, 4.0], ids=["undamped", "massPropDamped"])
def test_ladrunohht_ddm_matches_finite_difference(alphaM):
    """DDM d(U2x)/dE from LadrunoHHT matches a central FD of the response."""
    _, ddm = _build(E0, ('LadrunoHHT', ALPHA), with_sens=True, alphaM=alphaM)
    fd = _fd_gradient(alphaM)
    assert abs(fd) > 0.0
    rel_err = abs(ddm - fd) / abs(fd)
    assert rel_err < 1.0e-5, (
        f"DDM vs FD mismatch (alphaM={alphaM}): ddm={ddm:.6e} fd={fd:.6e} "
        f"rel_err={rel_err:.3e}"
    )


def test_ladrunohht_reduces_to_newmark_ddm_at_alpha1():
    """At alpha == 1 the HHT-specific (1-alpha) terms drop out: LadrunoHHT's DDM
    must equal Newmark's DDM (avg-accel: gamma=0.5, beta=0.25) on the same model."""
    _, ddm_lad = _build(E0, ('LadrunoHHT', 1.0), with_sens=True)
    _, ddm_nm = _build(E0, ('Newmark', 0.5, 0.25), with_sens=True)
    assert abs(ddm_nm) > 0.0
    rel_err = abs(ddm_lad - ddm_nm) / abs(ddm_nm)
    assert rel_err < 1.0e-8, (
        f"LadrunoHHT(alpha=1) DDM must match Newmark DDM: lad={ddm_lad:.8e} "
        f"nm={ddm_nm:.8e} rel_err={rel_err:.3e}"
    )
