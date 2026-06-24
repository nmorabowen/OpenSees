"""B2 follow-up — frictional SOFT=2 (Courant-stable Coulomb/Tresca on the segment-based explicit penalty).

B2 (#406) shipped the SOFT=2 segment-based explicit NORMAL penalty as a frictionless MVP. This adds
friction on the SAME clipped segment overlap: per IN-CONTACT slave node (the soft-penalty mask p_I<0)
the SHIPPED LadrunoFrictionKernel return map runs with the SOFT normal pressure N_I=|p_I| and a
Courant-stable STICK penalty k_soft_t,I = SOFSCL·4·m_eff_t,I/dt² sized from the SEGMENT-BASED
tangential gap-mode mass (the B2 normal m_eff with n→t, worst-cased over the two tangents) — so a
stiff friction kt never throttles the explicit dt_cr via the stick mode (the B1-kt rule, generalized
to the segment operator). Penalty friction (λ_T≡0, single-pass — the explicit SOFT=2 philosophy);
the slip is the DISPLACEMENT-based weighted tangential gap (the C3.1 closest-point lesson).

Gates:
  1. FRICTION FORCE — a SOFT=2 contact dragged tangentially is arrested by Coulomb friction (the
     tangential travel collapses vs the frictionless free slide) and the travel decreases monotonically
     with μ (the force opposes motion and scales with μN).
  2. softKt STICK STABILITY — a STIFF `-epsT` (which would diverge via ω_t·dt≫2) runs STABLE under
     explicit because softKt sizes the stick penalty Courant-stably (the headline — the stick mode no
     longer throttles dt_cr).
  3. μ=0 BYTE-IDENTITY — no friction params ⇒ the friction block short-circuits ⇒ bit-identical to the
     frictionless SOFT=2 trajectory.
  4. (command accept is covered by test_adr39_contact_p5_soft2: `-mortar -soft -mu` now succeeds.)

Mechanics (the tangential m_eff_t + the cone over the soft pressure) validated oracle-first in
contact_prototypes/proto_soft2_friction.py (12/12).
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

M = 1.0
DT = 1.0e-3
KN = 1.0e6


def _fixed_master():
    for t, (x, y) in zip((1, 2, 3, 4), [(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(t, float(x), float(y), 0.0)
        ops.fix(t, 1, 1, 1)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)


def _slave_xz_free(z):
    """A matched slave facet with x,z FREE (y fixed) so it can be pressed down AND dragged tangentially."""
    for t, (x, y) in zip((5, 6, 7, 8), [(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(t, float(x), float(y), z)
        ops.mass(t, M, M, M)
        ops.fix(t, 0, 1, 0)


def _drag(mu, epsT=None, soft=0.1, P=200.0, vx=0.5, nsteps=400, force_mu=False,
         cohesion=None, tauMax=None):
    """Press the matched slave facet onto the fixed master (downward load P, normal contact), give it a
    tangential x-velocity, run SOFT=2 explicit. Coulomb friction opposes the slide ⇒ the x-travel of
    node 5 collapses vs the frictionless free slide vx·t. Returns (stable, x_travel). force_mu adds an
    explicit `-mu 0.0` even when μ=0 (to test the friction-param-present-but-zero opt-in inertness)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _fixed_master()
    _slave_xz_free(0.0)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    args = ["-mortar", "-epsN", KN, "-soft", soft]
    if mu > 0.0 or force_mu:
        args += ["-mu", mu]
    if cohesion is not None:
        args += ["-cohesion", cohesion]
    if tauMax is not None:
        args += ["-tauMax", tauMax]
    if epsT is not None:
        args += ["-epsT", epsT]
    args += ["-outward", 0.0, 0.0, 1.0]
    ops.contact(1, 1, 2, *args)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, 0.0, 0.0, -P)
        ops.setNodeVel(t, 1, vx, "-commit")     # tangential (x) velocity
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    stable = True
    for _ in range(nsteps):
        if ops.analyze(1, DT) != 0:
            stable = False
            break
        if abs(ops.nodeDisp(5, 1)) > 5.0:
            stable = False
            break
    return stable, ops.nodeDisp(5, 1)


def test_soft2_friction_arrests_tangential_drag():
    """GATE 1 — the headline. A SOFT=2 contact dragged tangentially at vx slides FREELY without
    friction (x ≈ vx·t) but is ARRESTED by Coulomb friction (the travel collapses), and the travel
    decreases MONOTONICALLY with μ — proving the friction force is active on the SOFT=2 segment overlap
    and opposes the slip ∝ μN. Same geometry + dt — only μ differs."""
    st0, x0 = _drag(0.0)            # frictionless: free slide
    st3, x3 = _drag(0.3)
    st6, x6 = _drag(0.6)
    assert st0 and st3 and st6, "SOFT=2 friction drag should stay stable at the structural dt"
    # frictionless free slide ≈ vx·t = 0.5·(400·1e-3) = 0.20
    assert x0 > 0.15, f"frictionless SOFT=2 should slide freely (x={x0:.4f})"
    # friction collapses the travel (an order of magnitude+) and it shrinks with μ
    assert x3 < 0.25 * x0, f"-mu must arrest the drag (x(0.3)={x3:.4f} vs free {x0:.4f})"
    assert x6 < x3, f"tangential travel must decrease with μ (x(0.6)={x6:.4f} !< x(0.3)={x3:.4f})"


def _stick_stability(epsT, mu=2.0, soft=0.1, P=200.0, nsteps=2000):
    """A SOFT=2 contact in STICK (tiny tangential kick) at a high μ. Under explicit the softKt sizing
    REPLACES the configured (possibly stiff) -epsT with the Courant-stable k_soft_t, so a stiff epsT
    that would diverge (ω_t·dt≫2) stays BOUNDED. Returns (stable, max|x|)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _fixed_master()
    _slave_xz_free(0.0)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    ops.contact(1, 1, 2, "-mortar", "-epsN", KN, "-soft", soft, "-mu", mu, "-epsT", epsT,
                "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, 0.0, 0.0, -P)
        ops.setNodeVel(t, 1, 0.01, "-commit")    # small tangential kick to excite the stick mode
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    mx = 0.0
    for _ in range(nsteps):
        if ops.analyze(1, DT) != 0:
            return False, mx
        mx = max(mx, abs(ops.nodeDisp(5, 1)))
        if mx > 1.0:
            return False, mx
    return True, mx


def test_soft2_friction_softkt_stick_stable():
    """GATE 2 — the headline stability win. A STIFF `-epsT 1e9` stick penalty would give ω_t·dt ≈ 30 ≫ 2
    (divergent under central difference at the structural dt), but softKt sizes the stick penalty
    Courant-stably (k_soft_t = SOFSCL·4·m_eff_t/dt² ⇒ ω_t·dt = 2√SOFSCL < 2), so the SOFT=2 friction
    stays BOUNDED — the friction kt no longer throttles dt_cr (the B1-kt rule, segment lane)."""
    stable, mx = _stick_stability(epsT=1.0e9)
    assert stable, f"softKt must keep a stiff friction kt stable under explicit (max|x|={mx:.3e})"
    assert mx < 1.0e-3, f"the stick response should stay bounded near the origin (max|x|={mx:.3e})"


def test_soft2_friction_cohesion_and_tresca_cone():
    """GATE — the unified cone `cap = min(μ|N|+c, τmax)` reaches the FE through `-cohesion`/`-tauMax`:
    (a) COHESION (μ=0, c>0) gives a shear strength INDEPENDENT of the normal pressure ⇒ it arrests the
        drag even with no friction coefficient, and more c ⇒ more arrest;
    (b) a TRESCA cap (μ>0, small τmax) LIMITS the friction below pure Coulomb ⇒ MORE slide than the
        uncapped Coulomb cone (the cap is an upper bound on the shear)."""
    _, x_free = _drag(0.0)
    _, x_c50 = _drag(0.0, cohesion=50.0)
    _, x_c200 = _drag(0.0, cohesion=200.0)
    assert x_c50 < 0.25 * x_free, f"cohesion must arrest the drag at μ=0 (x={x_c50:.4f} vs free {x_free:.4f})"
    assert x_c200 < x_c50, f"more cohesion ⇒ more arrest (x(200)={x_c200:.4f} !< x(50)={x_c50:.4f})"
    # Tresca cap < pure Coulomb at the same μ ⇒ less friction ⇒ more travel
    _, x_coulomb = _drag(0.6)
    _, x_tresca = _drag(0.6, tauMax=20.0)
    assert x_tresca > x_coulomb, (
        f"a Tresca cap must limit friction below Coulomb (x_capped={x_tresca:.4f} !> x_coulomb={x_coulomb:.4f})")


def test_soft2_friction_mu0_byte_identical():
    """GATE 3 — μ=0 short-circuits the friction block (`mu>0 ∧ c>0 ∧ τmax>0` all false) ⇒ a contact
    with an explicit `-mu 0.0` is BIT-IDENTICAL to one with no friction params at all (the friction
    addition is strictly opt-in / byte-identical when off)."""
    _, x_nofric = _drag(0.0)                    # no -mu at all
    _, x_mu0 = _drag(0.0, force_mu=True)        # explicit -mu 0.0 ⇒ short-circuits ⇒ must match exactly
    assert x_nofric == x_mu0, f"explicit -mu 0.0 perturbed the frictionless SOFT=2: {x_nofric!r} vs {x_mu0!r}"
    assert x_nofric > 0.15, f"the frictionless reference should slide freely (x={x_nofric:.4f})"
