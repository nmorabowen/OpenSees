"""ADR-41 C3.3 — finish C3: tangential λ_T Uzawa + the non-symmetric mortar Csl.

C3.1/C3.2 shipped PENALTY friction (force + symmetric tangent). C3.3 adds:
  - the tangential augmented-Lagrange multiplier λ_T (the friction analogue of λ_N): the residual
    injects it via the OFFSET TRICK gTeff_eff = gTeff + λ_T/epsT (reusing the penalty return map),
    and one Uzawa step per Domain::commit (λ_T ← −tFric) drives the spurious STICK elastic creep →
    0 at FINITE epsT (epsT-INDEPENDENT tangential position) — the held-load analyze_augmented proc
    augments it across commits, exactly as it augments λ_N;
  - `-consistanttan` is now FUNCTIONAL for mortar: the non-symmetric Coulomb Csl (the pressure
    coupling −μ·epsN·t̂⊗n) is assembled (needs a non-symmetric solver) for quadratic Coulomb Newton.

Gates:
  - held-load augmentation drives the stick tangential creep → ~0, epsT-INDEPENDENT (vs the penalty
    creep ∝ 1/epsT) — the λ_T win;
  - `-consistanttan` on a slipping Coulomb config converges on FullGeneral and reaches the SAME
    converged state as the symmetric tangent (the tangent only changes the Newton path, not the root).

Mechanics validated oracle-first in proto_c3_mortar_friction.py (T7 consistent Csl FD ~1e-11; T8
λ_T Uzawa stick-creep → 0 epsT-independent).
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

DELTA = 1.0e-4
COH = 50.0        # Tresca cohesion ⇒ cap = COH·A, pressure-independent


def _master():
    for t, x, y in [(1, 0, 0), (2, 1, 0), (3, 1, 1), (4, 0, 1)]:
        ops.node(t, float(x), float(y), 0.0)
        ops.fix(t, 1, 1, 1)


def _slave(z, fix_z):
    for t, x, y in [(5, 0, 0), (6, 1, 0), (7, 1, 1), (8, 0, 1)]:
        ops.node(t, float(x), float(y), z)
        ops.fix(t, 0, 1, 1 if fix_z else 0)            # x free; y fixed; z per arg


def _lambdaT_creep(epsN, augment):
    """A slave facet in normal contact (z fixed at −DELTA ⇒ Tresca cap engaged) under a tangential
    stick load Q < cap. Returns the slave x-disp after analyze(1) [penalty creep] then, if augment,
    after held-load re-commits [λ_T Uzawa drives the creep → 0]. Tresca ⇒ cap is N-independent."""
    Q = 20.0                                           # < cap = COH·A = 50 ⇒ stick
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _master()
    _slave(-DELTA, fix_z=True)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    ops.contact(1, 1, 2, "-mortar", "-epsN", epsN, "-epsT", epsN, "-cohesion", COH,
                "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, Q / 4.0, 0.0, 0.0)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-11, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0, "penalty stick solve did not converge"
    u_penalty = ops.nodeDisp(5, 1)
    if not augment:
        return u_penalty, None
    # held-load augmentation: each zero-increment commit performs one λ_T Uzawa step.
    ops.integrator("LoadControl", 0.0)
    for _ in range(15):
        assert ops.analyze(1) == 0, "held-load λ_T augmentation step did not converge"
    return u_penalty, ops.nodeDisp(5, 1)


def test_c3_3_lambdaT_uzawa_reduces_stick_creep_epsN_independent():
    """THE λ_T WIN: penalty friction leaves a stick elastic creep ∝ 1/epsT; the tangential Uzawa
    drives it → ~0 at FINITE epsT, and the augmented tangential position is epsT-INDEPENDENT
    (penalty alone cannot reach it). Mirrors the C2.2 normal λ_N headline, on the tangential side."""
    u_pen_lo, u_aug_lo = _lambdaT_creep(epsN=1.0e4, augment=True)
    u_pen_hi, u_aug_hi = _lambdaT_creep(epsN=1.0e6, augment=True)

    # penalty creep is epsT-dependent (∝ 1/epsT): 100× stiffer epsT ⇒ ~100× less creep.
    assert u_pen_lo > 20.0 * u_pen_hi > 0.0, (
        f"penalty creep not ∝1/epsT: {u_pen_lo:.2e} vs {u_pen_hi:.2e}")
    # λ_T augmentation drives the creep → ~0 at BOTH epsT, and the result is epsT-INDEPENDENT.
    assert abs(u_aug_lo) < 1.0e-6, f"λ_T did not remove the creep at epsT=1e4 (u={u_aug_lo:.2e})"
    assert abs(u_aug_hi) < 1.0e-6, f"λ_T did not remove the creep at epsT=1e6 (u={u_aug_hi:.2e})"
    assert abs(u_aug_lo - u_aug_hi) < 1.0e-7, (
        f"augmented position epsT-dependent: {u_aug_lo:.2e} vs {u_aug_hi:.2e}")
    # and the augmentation actually SHRANK the creep (not a no-op).
    assert abs(u_aug_lo) < 1.0e-2 * abs(u_pen_lo), "augmentation did not reduce the creep"


def _coulomb_slip(consistent, system):
    """A spring-restrained slave in Coulomb SLIP (Q > μN). Returns (rc, x-disp). With -consistanttan
    + a non-symmetric solver the full Coulomb tangent is used; the converged state is solver/tangent-
    independent (the tangent only changes the Newton path)."""
    Q, Pz, kx, mu = 120.0, 200.0, 1.0e3, 0.5
    epsN = 1.0e4
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _master()
    _slave(-DELTA, fix_z=True)
    ops.uniaxialMaterial("Elastic", 1, kx)
    for i, t in enumerate((5, 6, 7, 8)):
        gt = 200 + t
        ops.node(gt, *[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0)][i])
        ops.fix(gt, 1, 1, 1)
        ops.element("zeroLength", 300 + t, gt, t, "-mat", 1, "-dir", 1)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    opts = ["-mortar", "-epsN", epsN, "-mu", mu, "-epsT", epsN, "-outward", 0.0, 0.0, 1.0]
    if consistent:
        opts += ["-consistanttan"]
    ops.contact(1, 1, 2, *opts)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, Q / 4.0, 0.0, -Pz / 4.0)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system(system)
    ops.test("NormDispIncr", 1.0e-10, 60, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    ops.integrator("LoadControl", 1.0)
    rc = ops.analyze(1)
    return rc, ops.nodeDisp(5, 1), ops.testIter()


def test_c3_3_consistanttan_coulomb_converges_same_root():
    """`-consistanttan` (non-symmetric Coulomb Csl) is now FUNCTIONAL for mortar: a slipping Coulomb
    config converges on a non-symmetric solver (FullGeneral) to the SAME converged state as the
    default symmetric tangent (the tangent changes only the Newton path, not the root)."""
    rc_sym, ux_sym, it_sym = _coulomb_slip(consistent=False, system="FullGeneral")
    rc_cons, ux_cons, it_cons = _coulomb_slip(consistent=True, system="FullGeneral")
    assert rc_sym == 0, "symmetric-tangent Coulomb slip did not converge"
    assert rc_cons == 0, "-consistanttan Coulomb slip did not converge on FullGeneral"
    assert abs(ux_sym - ux_cons) < 1.0e-7, (
        f"consistent vs symmetric reached different roots: {ux_cons} vs {ux_sym}")
    # the consistent (non-symmetric Coulomb) tangent is the EXACT derivative ⇒ it must not need
    # MORE Newton iterations than the symmetric approximation (quadratic vs superlinear) — the
    # justification for offering -consistanttan at all.
    assert it_cons <= it_sym, f"consistent tangent took MORE iters ({it_cons}) than symmetric ({it_sym})"
