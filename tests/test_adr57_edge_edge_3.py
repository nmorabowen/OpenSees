"""ADR-57 E5 — edge-edge explicit Courant-stable SOFT penalty.

The EDGE_EDGE mode gains the SOFT lane (the B1/B2 SOFT analogue for the edge operator): under the
explicit CentralDifferenceLadruno, `-edgeSoft <SOFSCL>` REPLACES the edge penalty εN with the
Courant-stable k_soft = SOFSCL·4·m_eff/dt² (m_eff = the 4-node edge gap-mode generalized mass,
1/((1−s)²invMproj_a0 + s²invMproj_a1 + (1−t)²invMproj_b0 + t²invMproj_b1)), so the edge contact's
own central-difference mode ω·dt = 2√SOFSCL ≤ 2 NEVER throttles dt_cr — edge-on impact runs at the
STRUCTURAL dt. Inert under implicit / -edgeSoft-absent ⇒ byte-identical. Mechanics validated
oracle-first in contact_prototypes/proto_e5_soft.py.

Two crossed bars: a slave bar's bottom edge (∥x at z=+0.02) drops with a downward velocity onto a
fixed master bar's top edge (∥y at z=0), n=ẑ. m_eff = 2 (two unit slave masses, fixed master, the
symmetric crossing). At the STRUCTURAL dt=1e-3 the configured stiff εN=1e8 has ω·dt = √(εN/2)·dt ≈
7 ≫ 2 (UNSTABLE) while SOFT has ω·dt = 2√0.1 ≈ 0.63 ≤ 2 (stable).

Gates:
  - EXPLICIT SOFT impact at the structural dt: bounded penetration, the bar REBOUNDS (no blow-up);
  - the SAME dt with the configured stiff εN DIVERGES (ω·dt ≫ 2) — SOFT lets you run where stiff can't;
  - IMPLICIT byte-identity: under a static (non-CDL) solve, -edgeSoft is BITWISE identical to its absence.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

EPS = 1.0e8       # the configured (stiff) edge-edge penalty εN — throttles dt_cr under explicit
ZC = 0.02         # slave bottom edge seeded just ABOVE the master (separated ⇒ a real impact)
DT = 1.0e-3       # the STRUCTURAL dt; SOFT must run stably at it, the stiff εN must not
V0 = -2.0         # downward approach velocity of the slave bottom edge


def _build_bars():
    """Master (fixed) bar in the y-z plane, top edge m1-m2 ∥ y at z=0. Slave bar in the x-z plane,
    bottom edge s1-s2 ∥ x at z=ZC (above). n = ê_s×ê_m = ẑ. Slave top edge (7,8) anchored."""
    for t, x, y, z in [(1, 0.0, -1.0, 0.0), (2, 0.0, 1.0, 0.0),
                       (3, 0.0, 1.0, -0.3), (4, 0.0, -1.0, -0.3)]:
        ops.node(t, x, y, z)
        ops.fix(t, 1, 1, 1)
    for t, x, y, z in [(5, -1.0, 0.0, ZC), (6, 1.0, 0.0, ZC),
                       (7, 1.0, 0.0, ZC + 0.3), (8, -1.0, 0.0, ZC + 0.3)]:
        ops.node(t, x, y, z)
    for t in (7, 8):
        ops.fix(t, 1, 1, 1)


def _edge_contact(soft):
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    args = [1, 1, 2, "-mortar", "-epsN", EPS, "-edgeedge", "-edgeBand", 0.15,
            "-outward", 0.0, 0.0, 1.0]
    if soft:
        args += ["-edgeSoft", 0.1]
    ops.contact(*args)


def _impact_history(soft, nsteps=45):
    """EXPLICIT central-difference impact at the structural dt. Returns the absolute-z history of
    node 5 (z = ZC + disp). The bottom edge is free in z with mass + a downward velocity; x,y fixed."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _build_bars()
    for t in (5, 6):
        ops.fix(t, 1, 1, 0)               # x,y fixed; z free (the impact DOF)
        ops.mass(t, 1.0, 1.0, 1.0)
        ops.setNodeVel(t, 3, V0, "-commit")
    _edge_contact(soft)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    z = [ZC + ops.nodeDisp(5, 3)]
    for _ in range(nsteps):
        ops.analyze(1, DT)                # explicit: no convergence check (diverging runs still "pass")
        zi = ZC + ops.nodeDisp(5, 3)
        z.append(zi)
        if not math.isfinite(zi) or abs(zi) > 1.0e6:
            break                         # diverged ⇒ stop feeding NaN/huge back into the integrator
    return z


def test_e5_soft_runs_at_structural_dt():
    """SOFT impact at the structural dt: ω·dt = 2√0.1 ≤ 2 ⇒ bounded penetration + a clean REBOUND
    (the bar comes back up), no blow-up. The analytic penetration is d = |v0|·√(m_eff/k_soft)."""
    z = _impact_history(soft=True)
    assert all(math.isfinite(zi) for zi in z), "SOFT impact produced non-finite z (diverged?)"
    deepest = min(z)
    # analytic SOFT penetration: m_eff=2, k_soft=0.1·4·2/dt²=8e5 ⇒ d=|v0|√(m_eff/k_soft)
    k_soft = 0.1 * 4.0 * 2.0 / (DT * DT)
    d_analytic = abs(V0) * math.sqrt(2.0 / k_soft)
    assert -deepest < 5.0 * d_analytic, f"SOFT penetration {-deepest} unbounded vs analytic {d_analytic}"
    assert deepest < 0.0, "the bar never penetrated (impact did not occur)"
    # rebound: the bar recovers ABOVE the deepest point and ends moving up (energy-restituting, bounded)
    assert z[-1] > deepest + 0.5 * (-deepest), f"SOFT bar did not rebound (z_end={z[-1]}, deepest={deepest})"
    assert ops.nodeVel(5, 3) > 0.0, "SOFT bar not moving up at the end (no rebound)"


def test_e5_stiff_throttles_diverges():
    """THE CONTRAST: the SAME structural dt with the configured stiff εN=1e8 has ω·dt = √(εN/2)·dt ≈
    7 ≫ 2 ⇒ the central-difference contact mode is UNSTABLE and the integration blows up. SOFT runs
    where the stiff penalty cannot (it would force dt_cr ~ 2/ω ≈ 2.8e-4, ~4× smaller)."""
    z = _impact_history(soft=False)
    omega_dt = math.sqrt(EPS / 2.0) * DT
    assert omega_dt > 2.0, f"sanity: stiff εN should be unstable at this dt (ω·dt={omega_dt})"
    # divergence ⇒ a non-finite z OR an absurd excursion (|z| ≫ the 0.02 drop) the SOFT run never shows
    diverged = any((not math.isfinite(zi)) or abs(zi) > 1.0 for zi in z)
    assert diverged, f"stiff εN did NOT diverge at the structural dt (max|z|={max(abs(zi) for zi in z)})"


def test_e5_implicit_byte_identical():
    """IMPLICIT byte-identity: under a static (non-CDL) solve the SOFT path falls through to the
    configured εN (the dynamic_cast to CentralDifferenceLadruno fails), so -edgeSoft is BITWISE
    identical to its absence. A pressed, contact-ACTIVE crossing (pre-penetrating, modest load —
    the proven E2 static setup) compared node-for-node, with -edgeSoft toggled."""
    ZC_PEN = -1.0e-4          # pre-penetrating ⇒ the contact is ACTIVE from step 1 (the soft path runs)
    EPS_IMPL = 1.0e6
    KSPRING = 500.0

    def _static_push(soft):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        for t, x, y, z in [(1, 0.0, -1.0, 0.0), (2, 0.0, 1.0, 0.0),
                           (3, 0.0, 1.0, -0.3), (4, 0.0, -1.0, -0.3)]:
            ops.node(t, x, y, z)
            ops.fix(t, 1, 1, 1)
        for t, x, y, z in [(5, -1.0, 0.0, ZC_PEN), (6, 1.0, 0.0, ZC_PEN),
                           (7, 1.0, 0.0, ZC_PEN + 0.3), (8, -1.0, 0.0, ZC_PEN + 0.3)]:
            ops.node(t, x, y, z)
        for t in (7, 8):
            ops.fix(t, 1, 1, 1)
        for t in (5, 6):
            ops.fix(t, 1, 1, 0)           # x,y fixed; z free
        # soft regularizing z-springs (the point-pair differential mode, as in E2)
        ops.uniaxialMaterial("Elastic", 1, KSPRING)
        for gnd, sn, x in [(305, 5, -1.0), (306, 6, 1.0)]:
            ops.node(gnd, x, 0.0, ZC_PEN)
            ops.fix(gnd, 1, 1, 1)
            ops.element("zeroLength", 900 + sn, gnd, sn, "-mat", 1, "-dir", 3)
        ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
        ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
        args = [1, 1, 2, "-mortar", "-epsN", EPS_IMPL, "-edgeedge", "-edgeBand", 0.15,
                "-outward", 0.0, 0.0, 1.0]
        if soft:
            args += ["-edgeSoft", 0.1]
        ops.contact(*args)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        for t in (5, 6):
            ops.load(t, 0.0, 0.0, -200.0)             # press the bottom edge DOWN into the master
        ops.constraints("LadrunoContact")
        ops.numberer("Plain")
        ops.system("FullGeneral")
        ops.integrator("LoadControl", 1.0)
        ops.test("NormDispIncr", 1.0e-11, 40, 0)
        ops.algorithm("Newton")
        ops.analysis("Static")
        out = []
        for _ in range(3):
            assert ops.analyze(1) == 0, "implicit edge-edge push did not converge"
            out.append(tuple(ops.nodeDisp(t, 3) for t in (5, 6)))
        return out
    soft, plain = _static_push(True), _static_push(False)
    assert soft == plain, "implicit -edgeSoft not bitwise identical to its absence"
    # RIGOR: prove the EDGE-EDGE contact actually GOVERNS (so the soft-fallback path is genuinely
    # exercised, not trivially identical via spring-dominated displacements). εN=1e6 ≫ 2k=1e3 ⇒ the
    # contact carries ~all the load: the bottom edge settles at the εN-dominated depth (~P/εN ≈ 4e-4
    # from the ZC_PEN seed), ~10³× shallower than the pure-spring depth P/2k = 0.4. If the springs
    # governed, |z| would be ~0.4 and the contact (hence the soft fallback) would be inert.
    z5_abs = ZC_PEN + soft[-1][0]
    # 3 LoadControl steps × (200 N × 2 nodes) = 1200 N total; the two springs (Σk = 2·KSPRING = 1000)
    # alone would deflect 1200/1000 = 1.2. εN=1e6 ≫ 2k ⇒ the contact carries ~all of it: |z5| ≈ 1.2e-3,
    # ~10³× shallower. If the springs governed (soft path inert) |z5| would be ~1.2.
    spring_only_depth = 3.0 * (2.0 * 200.0) / (2.0 * KSPRING)
    assert abs(z5_abs) < 0.01 * spring_only_depth, \
        f"contact did not govern (|z5|={z5_abs} vs spring-only {spring_only_depth}); soft path not exercised"
