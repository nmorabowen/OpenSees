"""LadrunoUP (ADR-71 P4, WP4.C) — PDMY liquefaction-column battery.

Owns the P4 free-field liquefaction gates (plan §WP4.A, 4.C pins): upstream
quadUP reference vs LadrunoUP Q4, PressureDependMultiYield (PDMY) both sides,
the gravity(elastic stage 0) -> updateMaterialStage 1 -> base-shaking idiom,
and THE novel bit — the §3.2 obligation-3 gate that the stage flip reaches the
materials THROUGH the element AND dirties the auto-alpha stabilization cache.

PDMY PARAMETER SET (module constant PDMY below): a medium-loose sand in the
UCSD modern 16-argument signature
    nDMaterial PressureDependMultiYield tag nd rho refShear refBulk phi
        gammaPeak refPress d phaseTransAng contract dil1 dil2 liq1 liq2 liq3 NYS
following the staging idiom of the upstream McGann/Arduino workshop example
(Workshops/OpenSeesDays/site/freeFieldDepend.tcl: build stage-0 elastic ->
gravity -> updateMaterialStage -stage 1 -> setTime 0/wipeAnalysis -> shake) with
liquefaction switched ON (liq1=10 kPa, liq2=0.02, liq3=1 — the workshop file
zeroes them; contract=0.17 accelerates contraction so the short record
liquefies). rho=2 Mg/m3, G=5.5e4, B=1.5e5 kPa, phi=33, refPress=80 kPa, d=0.5.

ARG-MAPPING ADDENDUM to the P1 equiv battery (verified here, one-element
bit-exact cross-check during scaffolding): upstream quadUP's <b1 b2> drive BOTH
the solid-momentum rows AND the fluid seepage source, so LadrunoUP parity needs
`-body b1 b2 -fluidBody b1 b2` (setting -fluidBody 0 zeroes the hydrostatic
steady p and shifts the effective-stress baseline by exactly p_hydro — the
first scaffolding run had sigma'_v off by 2x through this).

STAGING QUIRK (documented, drives the leg-1 design): the always-on ctor body
force is a STEP load at t=0; under the undamped gamma=1/2 tight leg it rings
and the p-as-disp vs p-as-vel parameterizations disagree O(1). The tight leg
therefore uses the P1 methodology's RAMPED top surcharge (body force zero);
the body-force column runs under the production gamma=0.6 set only.

CACHE-DIRTY GATE (§3.2 obligation 3, the 4.C headline) — instruments:
  * p-p-block alpha meter: in the assembled transient tangent the p-row/p-col
    entries are c1*H + c2*(S + alpha*Htilde) — ALL state-independent except
    through alpha. Two manual-alpha runs calibrate the per-entry scale
    (apre*c2*Htilde); any run's alpha is then read EXACTLY from its p-p
    entries, immune to trajectory divergence (K-block contamination sank the
    naive whole-matrix extraction during scaffolding).
  * twin identity (the P2 trick, pin-literal): -stab auto 0.25 vs
    -stab <manual alpha_post> where alpha_post = alpha0*h^2/M_post and M_post
    = the material tangent probed through the element AFTER the flip
    (eleResponse 'material' tangent — equals PDMY getInitialTangent at the
    settled post-flip state, loadStage-1 branch: factor(currentStress)*M0).
    Identical shaking trajectories prove auto refreshed to the POST-stage
    moduli; the PRE-stage control (alpha from M0 = refBulk + 4/3 refShear,
    stage-0 initial tangent) must disagree loudly or the cache never dirtied.
  First-passing-run numbers (single element): f = alpha_auto/alpha_pre =
  1.405093 constant over 120 shaking steps (pp-meter, all entries, all probe
  steps), == M0/M_probe to 5 digits; twin rel 2.4e-8; pre-stage control rel
  1.1e-1 (4.6e6x separation).

BUG (reported to MAIN, xfail'd here with repro + workaround):
  On MULTI-ELEMENT models updateMaterialStage dirties the stabilization/K0
  caches of ONLY ONE element. Chain: MaterialStageParameter::setDomain
  (SRC/domain/component/MaterialStageParameter.cpp:76) stops its domain scan at
  the FIRST element that accepts the parameter ("only need to find one" — true
  for the UCSD materials' static shared loadStage, NOT for per-element caches),
  so exactly one LadrunoUP co-registers (LadrunoUP.cpp:1270). Measured with
  the pp-meter on a 2-element shared-tag stack: f1 = 1.4051 (scan-found
  element, correct post-stage refresh — identical to the single-element gate)
  but f2 = 1.0000 (STALE stage-0 alpha). The materials themselves DO all flip
  (copy ctor shares the static matN slot, PressureDependMultiYield.cpp:426 —
  verified below), so the physics runs plastic; only the per-element caches go
  silently stale, i.e. the stabilization is off by the stage-modulus ratio
  (1.4x here) on every element except one.
  WORKAROUND (gated green below): after updateMaterialStage, poke every
  element with `setParameter -val <same k> -ele <e> xPerm` — updateParameter
  case 3 re-dirties both caches (LadrunoUP.cpp:1300) without changing physics
  (measured: f2 1.0000 -> 1.4051 with a same-value poke).

METER QUIRK (cost half a day of contradictory readings during scaffolding):
  nodeDOFs() on a Transformation-handler SLAVE node (equalDOF u-ties here)
  returns an equation map whose p-slot lands on a K-scale row of printA —
  garbage for the pp-meter. Use MASTER/untied nodes' p-DOFs only.

Run mechanics (worktree): hermetic bootstrap of THIS worktree's dist/bin —
invoke with `py -3.12 -S` and PYTHONPATH carrying pytest's site-packages
(project memory "OpenSees Python test env"). Battery runtime ~= 45 s.
"""
import math
import os
import sys

# --------------------------------------------------------------------------
# hermetic import: this worktree's build, not a stale preloaded pyd
# --------------------------------------------------------------------------
_WT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_BIN = os.path.join(_WT, "dist", "bin")
if os.path.isfile(os.path.join(_BIN, "opensees.pyd")) or os.path.isfile(
        os.path.join(_BIN, "opensees.so")):
    os.environ["PATH"] = _BIN + os.pathsep + os.environ.get("PATH", "")
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(_BIN)
    sys.path.insert(0, _BIN)
    import opensees as ops
    assert os.path.normcase(os.path.dirname(ops.__file__)) == \
        os.path.normcase(_BIN), (
        f"stale opensees import: {ops.__file__} (want {_BIN}); "
        "run with py -3.12 -S so the boot .pth cannot preload another build")
else:  # CI layout: opensees.so copied next to tests/
    from _testbed import ops  # noqa: F401

import numpy as np
import pytest

pytestmark = [pytest.mark.zone_a]

GRAV = 9.81
ALPHA0 = 0.25

# medium-loose sand, UCSD modern 16-arg signature (module docstring cites the
# workshop precedent; liquefaction params switched on, contract raised)
PDMY = dict(rho=2.0, G=5.5e4, B=1.5e5, phi=33.0, gammaPeak=0.1, refPress=80.0,
            d=0.5, PT=27.0, contract=0.17, dil1=0.0, dil2=0.0,
            liq1=10.0, liq2=0.02, liq3=1.0, NYS=20)
# stage-0 (elastic) constrained modulus M0 = B + 4G/3 — PDMY getInitialTangent
# at loadStage 0 is exactly the (G, B) elastic tangent (factor = 1).
M0 = PDMY["B"] + 4.0 * PDMY["G"] / 3.0

COL = dict(nE=10, b=1.0, dy=1.0, perm=1.0e-5, Kf=2.2e5, poro=0.4)


def _mk_pdmy(tag):
    p = PDMY
    ops.nDMaterial("PressureDependMultiYield", tag, 2, p["rho"], p["G"],
                   p["B"], p["phi"], p["gammaPeak"], p["refPress"], p["d"],
                   p["PT"], p["contract"], p["dil1"], p["dil2"], p["liq1"],
                   p["liq2"], p["liq3"], p["NYS"])


def _rel_inf(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    assert np.all(np.isfinite(a)) and np.all(np.isfinite(b))
    return float(np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-30))


# ==========================================================================
# shared free-field column (buoyant self-weight variant, quadUP-parity BCs)
# ==========================================================================
def _column_build(which, stab=("off",)):
    """1 x nE stack, simple-shear ties, base fixed+sealed, top drained.

    which in {'ladruno','quadup'}. bY = -g*(rho-1) — with rho=2 numerically
    equal to full g on the fluid; quadUP's b2 drives solid AND seepage rows,
    LadrunoUP parity passes it to BOTH -body and -fluidBody (module header).
    """
    c = COL
    nE = c["nE"]
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    nid = {}
    k = 1
    for j in range(nE + 1):
        for i in range(2):
            ops.node(k, i * c["b"], j * c["dy"])
            nid[(i, j)] = k
            k += 1
    _mk_pdmy(1)
    bY = -GRAV * (PDMY["rho"] - 1.0)
    for j in range(nE):
        n1, n2 = nid[(0, j)], nid[(1, j)]
        n3, n4 = nid[(1, j + 1)], nid[(0, j + 1)]
        if which == "ladruno":
            args = ["LadrunoUP", j + 1, n1, n2, n3, n4, 1,
                    "-Kf", c["Kf"], "-poro", c["poro"], "-rhoF", 1.0,
                    "-perm", c["perm"], c["perm"], "-dynSeepage", "off",
                    "-body", 0.0, bY, "-fluidBody", 0.0, bY]
            if stab[0] == "auto":
                args += ["-stab", "auto", stab[1]]
            elif stab[0] == "off":
                args += ["-stab", "off"]
            else:
                args += ["-stab", float(stab[1])]
            ops.element(*args)
        else:
            # quadUP: bulk = pre-combined Qbar ~= Kf/n (P1 equiv arg map)
            ops.element("quadUP", j + 1, n1, n2, n3, n4, 1.0, 1,
                        c["Kf"] / c["poro"], 1.0, c["perm"], c["perm"],
                        0.0, bY, 0.0)
    ops.fix(nid[(0, 0)], 1, 1, 0)
    ops.fix(nid[(1, 0)], 1, 1, 0)
    ops.fix(nid[(0, nE)], 0, 0, 1)
    ops.fix(nid[(1, nE)], 0, 0, 1)
    for j in range(1, nE + 1):
        ops.equalDOF(nid[(0, j)], nid[(1, j)], 1, 2)
    return nid


def _getp(which, n):
    """Contract conversion: honest-p disp slot vs upstream rate-DOF vel."""
    return ops.nodeDisp(n, 3) if which == "ladruno" else ops.nodeVel(n, 3)


def _sigv_eff(ele):
    """-sigma'_yy at GP1 (compression positive), material committed stress."""
    return -ops.eleResponse(ele, "material", 1, "stress")[1]


def _column_gravity(flip=True):
    """Upstream staging idiom: elastic gravity (big-dt Newmark 0.6 transient),
    updateMaterialStage 1, plastic re-equilibration."""
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-7, 50)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")
    for _ in range(15):
        assert ops.analyze(1, 200.0) == 0
    if flip:
        ops.updateMaterialStage("-material", 1, "-stage", 1)
        for _ in range(15):
            assert ops.analyze(1, 200.0) == 0


def _column_shake(which, nid, amp, tsh, dt=0.01, rayleigh=True):
    """setTime 0 / wipeAnalysis / UniformExcitation sine — returns per-step
    (surface ux, mid-depth r_u) plus sigma'_v0; identical Newton + substep
    fallback policy on both legs."""
    nE = COL["nE"]
    mid = nE // 2
    sv0 = _sigv_eff(mid)
    p0 = _getp(which, nid[(0, mid)])
    ops.setTime(0.0)
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.integrator("Newmark", 0.6, 0.3025)
    if rayleigh:
        w1, w2, z = 2 * math.pi * 1.0, 2 * math.pi * 3.0, 0.05
        ops.rayleigh(2 * z * w1 * w2 / (w1 + w2), 2 * z / (w1 + w2), 0.0, 0.0)
    ops.analysis("Transient")
    ops.timeSeries("Trig", 2, 0.0, tsh, 1.0 / 1.5, "-factor", amp)
    ops.pattern("UniformExcitation", 3, 1, "-accel", 2)
    surf, ru = [], []
    ns = int(round(tsh / dt))
    for s in range(ns):
        ops.test("NormDispIncr", 1e-5, 50)
        ops.algorithm("Newton")
        if ops.analyze(1, dt) != 0:
            ok = True
            ops.test("NormDispIncr", 1e-4, 100)
            ops.algorithm("KrylovNewton")
            for _ in range(10):
                if ops.analyze(1, dt / 10) != 0:
                    ok = False
                    break
            assert ok, f"{which} shaking hard-failed at t={s * dt:.2f}"
        surf.append(ops.nodeDisp(nid[(0, nE)], 1))
        ru.append((_getp(which, nid[(0, mid)]) - p0) / sv0)
    return np.array(surf), np.array(ru), sv0


# ==========================================================================
# 0. the transport itself: the stage flip reaches the materials THROUGH the
#    element (setParameter catch-all forwarding, §3.2 obligation 3 first half)
# ==========================================================================
def test_stage_flip_reaches_materials_through_element():
    """All 10 elements' GP materials report the stage-0 ELASTIC tangent
    (M0 = B + 4G/3 exactly) before the flip and a pressure-dependent stage-1
    tangent (!= M0, depth-graded) after — the flip travelled
    MaterialStageParameter -> element setParameter catch-all -> per-GP copies
    (all copies share the UCSD static loadStage, so every element flips even
    though the domain scan registers only one — module header)."""
    nid = _column_build("ladruno", stab=("auto", ALPHA0))
    nE = COL["nE"]
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-7, 50)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")
    for _ in range(10):
        assert ops.analyze(1, 200.0) == 0
    m_pre = [ops.eleResponse(e, "material", 1, "tangent")[0]
             for e in range(1, nE + 1)]
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    for _ in range(5):
        assert ops.analyze(1, 200.0) == 0
    m_post = [ops.eleResponse(e, "material", 1, "tangent")[0]
              for e in range(1, nE + 1)]
    for e, (a, b) in enumerate(zip(m_pre, m_post), 1):
        assert abs(a - M0) <= 1e-6 * M0, f"ele {e} pre-flip tangent {a} != M0"
        assert b < 0.95 * M0, (
            f"ele {e} post-flip tangent {b} still elastic — stage flip did "
            "not reach this element's materials")
    # depth grading: deeper (higher confinement) => stiffer post-flip tangent
    assert m_post[0] > m_post[-1]


# ==========================================================================
# 1a. gate (a) leg 1 (tight): gamma=1/2, beta=1/4 gravity-stage equivalence
# ==========================================================================
def test_quadup_equivalence_leg1_gravity_tight():
    """PDMY stage-0 consolidation column (RAMPED top surcharge — the step
    body force rings undamped, module header) vs FourNodeQuadUP under Newmark
    gamma=1/2 beta=1/4: full u AND p trajectories <= 1e-6 rel.
    First-passing-run: u 2.9e-15, p 4.1e-15 — machine-identical one-step rule,
    pinning the PDMY-lane arg mapping + contract conversion in one shot."""
    nE, b, dy = 4, COL["b"], COL["dy"]
    perm, Kf, poro = COL["perm"], COL["Kf"], COL["poro"]

    def run(which):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        nid = {}
        k = 1
        for j in range(nE + 1):
            for i in range(2):
                ops.node(k, i * b, j * dy)
                nid[(i, j)] = k
                k += 1
        _mk_pdmy(1)
        for j in range(nE):
            n1, n2 = nid[(0, j)], nid[(1, j)]
            n3, n4 = nid[(1, j + 1)], nid[(0, j + 1)]
            if which == "ladruno":
                ops.element("LadrunoUP", j + 1, n1, n2, n3, n4, 1,
                            "-Kf", Kf, "-poro", poro, "-rhoF", 1.0,
                            "-perm", perm, perm, "-stab", "off",
                            "-dynSeepage", "off")
            else:
                ops.element("quadUP", j + 1, n1, n2, n3, n4, 1.0, 1,
                            Kf / poro, 1.0, perm, perm, 0.0, 0.0, 0.0)
        for j in range(nE + 1):
            fy = 1 if j == 0 else 0
            fp = 1 if j == nE else 0
            ops.fix(nid[(0, j)], 1, fy, fp)
            ops.fix(nid[(1, j)], 1, fy, fp)
        ops.timeSeries("Path", 1, "-time", 0.0, 0.05, 1e6,
                       "-values", 0.0, 1.0, 1.0)
        ops.pattern("Plain", 1, 1)
        ops.load(nid[(0, nE)], 0.0, -5.0, 0.0)
        ops.load(nid[(1, nE)], 0.0, -5.0, 0.0)
        ops.constraints("Transformation")
        ops.numberer("Plain")
        ops.system("UmfPack")
        ops.test("NormUnbalance", 1e-11, 50)
        ops.algorithm("Newton")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient")
        U, P = [], []
        for _ in range(40):
            assert ops.analyze(1, 0.05) == 0
            U.append([ops.nodeDisp(nid[(i, j)], d)
                      for j in range(nE + 1) for i in range(2) for d in (1, 2)])
            P.append([_getp(which, nid[(i, j)])
                      for j in range(nE + 1) for i in range(2)])
        return np.array(U), np.array(P)

    uL, pL = run("ladruno")
    uU, pU = run("quadup")
    assert np.max(np.abs(uU)) > 1e-6 and np.max(np.abs(pU)) > 1e-2
    rel_u = _rel_inf(uL, uU)
    rel_p = _rel_inf(pL, pU)
    assert rel_u <= 1e-6, f"leg-1 u trajectories diverge: rel {rel_u:.3e}"
    assert rel_p <= 1e-6, f"leg-1 p trajectories diverge: rel {rel_p:.3e}"


# ==========================================================================
# 1b. gate (a) leg 2 (production): full staging + shaking, response band
# ==========================================================================
def test_quadup_equivalence_leg2_shaking_band():
    """Full idiom both sides (elastic gravity -> flip -> 4 s / 0.4 m/s^2 sine,
    gamma=0.6/beta=0.3025, 5% Rayleigh): pre-shake effective stress and the
    shaking response must agree at response level. Two different nonlinear
    element formulations => generous 10-15%% pinned bands, but with the
    -fluidBody parity mapping the measured gaps are tiny:
    sigma'_v0 rel 6.2e-6, surface-ux peak rel 9.7e-5, trajectory max 5.4e-4,
    r_u peak abs diff 6.1e-5 (first passing run)."""
    nid = _column_build("ladruno", stab=("off",))
    _column_gravity()
    sL, rL, svL = _column_shake("ladruno", nid, 0.4, 4.0)
    nid = _column_build("quadup")
    _column_gravity()
    sU, rU, svU = _column_shake("quadup", nid, 0.4, 4.0)
    assert abs(svL - svU) / svU <= 1e-3, f"sigma'_v0 mismatch {svL} vs {svU}"
    n = min(len(sL), len(sU))
    pkL, pkU = np.max(np.abs(sL[:n])), np.max(np.abs(sU[:n]))
    assert pkU > 1e-4  # non-trivial response guard
    assert abs(pkL - pkU) / pkU <= 0.10, (
        f"surface peak band: L {pkL:.4e} U {pkU:.4e}")
    assert _rel_inf(sL[:n], sU[:n]) <= 0.15, "surface trajectory band"
    ruL, ruU = np.max(rL[:n]), np.max(rU[:n])
    assert ruU > 0.05  # the shaking must actually build excess pressure
    assert abs(ruL - ruU) <= 0.05, f"r_u peak band: L {ruL:.3f} U {ruU:.3f}"


# ==========================================================================
# 2. gate (b): liquefaction onset — r_u rises toward 1 at mid-depth
# ==========================================================================
def test_liquefaction_onset_ru():
    """Production config (-stab auto 0.25), strong shaking (2.0 m/s^2, 6 s):
    mid-depth excess-pore-pressure ratio r_u = (p - p0)/sigma'_v0 rises toward
    1. Gates: r_u max >= 0.75 (measured 0.914), monotone build-up (quantile
    staircase), r_u holds near its peak at the end (undrained time scale >>
    record: perm 1e-5, 6 s). sigma'_v0 at mid-depth ~= gamma'*z = 54 kPa
    (buoyant, hydrostatic baseline) — sanity-pinned +-10%."""
    nid = _column_build("ladruno", stab=("auto", ALPHA0))
    _column_gravity()
    _, ru, sv0 = _column_shake("ladruno", nid, 2.0, 6.0)
    assert abs(sv0 - 54.0) <= 5.4, f"sigma'_v0 {sv0:.2f} != ~54 kPa"
    ru_max = float(np.max(ru))
    assert ru_max >= 0.75, f"no liquefaction onset: r_u max {ru_max:.3f}"
    q10, q50 = ru[int(len(ru) * 0.10)], ru[int(len(ru) * 0.50)]
    assert q10 < q50 <= ru_max + 1e-12, (
        f"r_u not building up: q10 {q10:.3f} q50 {q50:.3f} max {ru_max:.3f}")
    assert ru[-1] >= 0.8 * ru_max, "r_u dissipated implausibly fast"


# ==========================================================================
# 3. gate (c): cache-dirty — single-element twin + pp-block alpha meter
# ==========================================================================
# Single-element rig with a p-inert gravity stage: ALL p-DOFs fixed during
# gravity+flip (ptilde rows condensed out; Htilde cannot act), released for
# the shake — so the shake response depends ONLY on the post-flip alpha and a
# manual-alpha twin can close to machine level (measured 9.3e-16 in the
# elastic no-flip control during scaffolding).
_C1 = dict(perm=1.0e-5, Kf=2.2e5, poro=0.4, sur=-60.0)


def _one_build(stab):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    for i, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], 1):
        ops.node(i, float(x), float(y))
    _mk_pdmy(1)
    args = ["LadrunoUP", 1, 1, 2, 3, 4, 1, "-Kf", _C1["Kf"],
            "-poro", _C1["poro"], "-rhoF", 1.0,
            "-perm", _C1["perm"], _C1["perm"], "-dynSeepage", "off"]
    if stab[0] == "auto":
        args += ["-stab", "auto", stab[1]]
    else:
        args += ["-stab", float(stab[1])]
    ops.element(*args)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 1)
    ops.fix(3, 0, 0, 1)
    ops.fix(4, 0, 0, 1)
    ops.equalDOF(3, 4, 1, 2)


def _one_gravity_flip(flip=True):
    """Ramped surcharge, all-p-fixed consolidation, optional stage flip +
    settle. Returns the post-phase material tangent M (GP1) probed THROUGH
    the element — the exact modulus class the auto-alpha re-probe consumes
    (PDMY getInitialTangent stage-1 == elastic-branch getTangent here)."""
    ops.timeSeries("Path", 1, "-time", 0.0, 10.0, 1e9, "-values", 0.0, 1.0, 1.0)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 0.0, _C1["sur"] / 2.0, 0.0)
    ops.load(4, 0.0, _C1["sur"] / 2.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-10, 50)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")
    for _ in range(10):
        assert ops.analyze(1, 5.0) == 0
    if flip:
        ops.updateMaterialStage("-material", 1, "-stage", 1)
        for _ in range(3):
            assert ops.analyze(1, 5.0) == 0
    return ops.eleResponse(1, "material", 1, "tangent")[0]


def _one_release_and_shake(collect_pp_at=(), amp=0.15, nst=120, dt=0.01):
    """Release base p, shake; optionally collect the 2x2 p-p block of the
    assembled tangent (printA COLUMN-major — P1 quirk) at given steps."""
    ops.remove("sp", 1, 3)
    ops.remove("sp", 2, 3)
    ops.setTime(0.0)
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-10, 150)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")
    ops.timeSeries("Trig", 2, 0.0, nst * dt, 0.5, "-factor", amp)
    ops.pattern("UniformExcitation", 3, 1, "-accel", 2)
    traj, pp, pdofs = [], {}, None
    for s in range(1, nst + 1):
        assert ops.analyze(1, dt) == 0
        if pdofs is None:
            pdofs = [ops.nodeDOFs(1)[2], ops.nodeDOFs(2)[2]]
        traj.append([ops.nodeDisp(n, 3) for n in (1, 2)]
                    + [ops.nodeDisp(3, 1)])
        if s in collect_pp_at:
            A = np.array(ops.printA("-ret"))
            N = int(round(math.sqrt(len(A))))
            A = A.reshape(N, N, order="F")
            pp[s] = np.array([[A[r, c] for c in pdofs] for r in pdofs])
    return np.array(traj), pp


def test_cache_dirty_stage_flip_single_element():
    """THE 4.C headline (§3.2 obligation 3, second half): the stage flip must
    dirty the auto-alpha cache so it re-probes the POST-stage initial tangent.

    Instruments (module header):
      (1) pp-block alpha meter: manual-apre and manual-2apre runs calibrate
          the per-entry Htilde scale; the auto run's alpha is then read
          exactly from its p-p entries at several shaking steps.
      (2) twin identity + control (pin-literal).

    Gates (first-passing-run values in brackets):
      * meter self-check: manual leg reads f = 1 at every probe step [1.000000]
      * f_auto == M0 / M_probe(post-flip) +-2% at EVERY probe step [1.405093
        vs 1.405093, constant across steps 1..120]
      * f_auto >= 1.25 — decisively NOT the stale stage-0 value 1.0 (this is
        the control's teeth: a broken transport reproduces the PRE-stage
        twin instead, which the multi-element xfail below demonstrates live)
      * twin: auto vs manual(alpha_post) trajectories <= 1e-6 [2.4e-8]
      * control: auto vs manual(alpha_pre) >= 1e-2 [1.1e-1]
    """
    apre = ALPHA0 / M0  # h = largest edge = 1 m exactly on this element
    probe_steps = (1, 40, 120)

    _one_build(("auto", ALPHA0))
    m_probe = _one_gravity_flip()
    tA, ppA = _one_release_and_shake(collect_pp_at=probe_steps)
    f_expect = M0 / m_probe
    assert f_expect > 1.25, f"probe M {m_probe} did not drop vs M0 {M0}"

    apost = apre * f_expect
    _one_build(("manual", apost))
    _one_gravity_flip()
    tM, _ = _one_release_and_shake()

    _one_build(("manual", apre))
    _one_gravity_flip()
    tP, ppP = _one_release_and_shake(collect_pp_at=probe_steps)

    _one_build(("manual", 2 * apre))
    _one_gravity_flip()
    _, pp2 = _one_release_and_shake(collect_pp_at=probe_steps)

    # (1) pp-meter: scale from the two manual runs at step 1 (state-blind)
    hscale = pp2[1] - ppP[1]
    assert np.min(np.abs(hscale)) > 0.0
    for s in probe_steps:
        f_man = 1.0 + np.median((ppP[s] - ppP[1]) / hscale)
        assert abs(f_man - 1.0) <= 1e-6, (
            f"meter self-check failed at step {s}: manual f {f_man}")
        f_auto = 1.0 + float(np.median((ppA[s] - ppP[1]) / hscale))
        assert abs(f_auto - f_expect) <= 0.02 * f_expect, (
            f"step {s}: auto alpha f={f_auto:.6f} != post-stage {f_expect:.6f}"
            " — cache refreshed to the wrong moduli")
        assert f_auto >= 1.25, (
            f"step {s}: auto alpha still at the stage-0 value (f={f_auto:.4f})"
            " — the flip never dirtied the cache (ELEMENT BUG)")

    # (2) twin identity + control
    rel_twin = _rel_inf(tA, tM)
    rel_ctrl = _rel_inf(tA, tP)
    assert np.max(np.abs(tA[:, :2])) > 1.0  # non-trivial p response
    assert rel_twin <= 1e-6, (
        f"post-stage twin broke: rel {rel_twin:.3e} (auto != manual "
        "alpha_post => the auto cache holds something else)")
    assert rel_ctrl >= 1e-2, (
        f"PRE-stage control matched the auto run (rel {rel_ctrl:.3e}) — the "
        "cache did NOT refresh on the stage flip: element bug")


# ==========================================================================
# 4. the multi-element cache-dirty BUG (xfail, loud) + verified workaround
# ==========================================================================
def _two_build(stab):
    """2-element shared-material-tag stack, all p fixed during gravity."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    nid = {}
    k = 1
    for j in range(3):
        for i in range(2):
            ops.node(k, i * 1.0, j * 1.0)
            nid[(i, j)] = k
            k += 1
    _mk_pdmy(1)
    for j in range(2):
        n1, n2 = nid[(0, j)], nid[(1, j)]
        n3, n4 = nid[(1, j + 1)], nid[(0, j + 1)]
        args = ["LadrunoUP", j + 1, n1, n2, n3, n4, 1, "-Kf", _C1["Kf"],
                "-poro", _C1["poro"], "-rhoF", 1.0,
                "-perm", _C1["perm"], _C1["perm"], "-dynSeepage", "off"]
        if stab[0] == "auto":
            args += ["-stab", "auto", stab[1]]
        else:
            args += ["-stab", float(stab[1])]
        ops.element(*args)
    ops.fix(nid[(0, 0)], 1, 1, 1)
    ops.fix(nid[(1, 0)], 1, 1, 1)
    for j in (1, 2):
        ops.fix(nid[(0, j)], 0, 0, 1)
        ops.fix(nid[(1, j)], 0, 0, 1)
        ops.equalDOF(nid[(0, j)], nid[(1, j)], 1, 2)
    return nid


def _two_run(stab, poke=False):
    """Gravity + flip (+ optional same-value xPerm poke = the workaround),
    release non-top p, one transient step, return the 4x4 p-p block
    (p-DOFs of nodes j=0 and j=1, left+right)."""
    nid = _two_build(stab)
    ops.timeSeries("Path", 1, "-time", 0.0, 10.0, 1e9, "-values", 0.0, 1.0, 1.0)
    ops.pattern("Plain", 1, 1)
    ops.load(nid[(0, 2)], 0.0, _C1["sur"] / 2.0, 0.0)
    ops.load(nid[(1, 2)], 0.0, _C1["sur"] / 2.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-10, 50)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")
    for _ in range(10):
        assert ops.analyze(1, 5.0) == 0
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    if poke:
        for e in (1, 2):
            ops.setParameter("-val", _C1["perm"], "-ele", e, "xPerm")
    for _ in range(3):
        assert ops.analyze(1, 5.0) == 0
    m_probe = ops.eleResponse(1, "material", 1, "tangent")[0]
    for j in (0, 1):
        for i in range(2):
            ops.remove("sp", nid[(i, j)], 3)
    ops.setTime(0.0)
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-10, 50)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")
    assert ops.analyze(1, 0.01) == 0
    # MASTER-node p-DOFs only — the slave (right) nodes' nodeDOFs p-slot is
    # garbage under Transformation ties (module-header meter quirk); the
    # left j=0 node is element-1-pure, left j=1 node mixes both elements.
    pdofs = [ops.nodeDOFs(nid[(0, 0)])[2], ops.nodeDOFs(nid[(1, 0)])[2],
             ops.nodeDOFs(nid[(0, 1)])[2]]
    A = np.array(ops.printA("-ret"))
    N = int(round(math.sqrt(len(A))))
    A = A.reshape(N, N, order="F")
    return np.array([[A[r, c] for c in pdofs] for r in pdofs]), m_probe


def _two_alpha_f(poke):
    """Per-element alpha factors (f1, f2) from the pp meter: rows/cols
    {p(0,0), p(1,0)} carry element 1 only; the p(0,1) diagonal carries both
    elements (equal geometry => equal Htilde diagonals:
    d22 = ((f1-1)+(f2-1))/2)."""
    ppA, m_probe = _two_run(("auto", ALPHA0), poke=poke)
    apre = ALPHA0 / M0
    ppP, _ = _two_run(("manual", apre), poke=poke)
    pp2, _ = _two_run(("manual", 2 * apre), poke=poke)
    hs = pp2 - ppP
    d = (ppA - ppP) / hs
    f1 = 1.0 + float(np.median([d[0, 0], d[0, 1], d[1, 1], d[0, 2], d[1, 2]]))
    f2 = 1.0 + 2.0 * float(d[2, 2]) - (f1 - 1.0)
    return f1, f2, M0 / m_probe


def test_cache_dirty_multi_element_workaround_xperm_poke():
    """WORKAROUND gate (passes): after updateMaterialStage, a same-value
    `setParameter -val k -ele e xPerm` poke on every element re-dirties the
    stabilization cache (LadrunoUP.cpp updateParameter case 3) — both elements
    of the shared-tag stack then read the POST-stage alpha
    (measured f1 = f2 = 1.4051 = the probe-expected post-stage value; without
    the poke f2 stays 1.0000 — the xfail below)."""
    f1, f2, f_expect = _two_alpha_f(poke=True)
    assert f_expect > 1.25
    for name, f in (("ele1", f1), ("ele2", f2)):
        assert f >= 1.25, (
            f"{name} alpha stale (f={f:.4f}) even WITH the xPerm poke — "
            "workaround broken")
        assert abs(f - f_expect) <= 0.10 * f_expect, (
            f"{name} refreshed alpha f={f:.4f} not in the post-stage band "
            f"(expected ~{f_expect:.4f})")


def test_cache_dirty_multi_element_shared_tag():
    """CONTRACT gate (§3.2 obligation 3): every element's auto-alpha re-probes
    the post-stage moduli after updateMaterialStage — on the 2-element
    shared-tag stack both f1 and f2 sit in the post-stage band like the
    single-element gate.

    HISTORY (P4 4.C finding, FIXED at MAIN integration): originally a strict
    xfail — MaterialStageParameter::setDomain stops its scan at the FIRST
    accepting element (MaterialStageParameter.cpp:76 'only need to find one',
    a materials-side assumption that holds for the UCSD static loadStage but
    not for per-element caches), so exactly one LadrunoUP co-registered and
    the rest kept stale auto-alpha (measured f1 = 1.4051 / f2 = 1.0000, a
    ~1.4x stabilization error). FIX: LadrunoUP::updateParameter(1) now
    BROADCASTS dirtyStageCaches() to every LadrunoUP in the domain (flags
    only, lazy recompute — over-dirtying is harmless). This test is the
    regression gate; the xPerm-poke workaround test above remains as a
    plumbing check of updateParameter case 3."""
    f1, f2, f_expect = _two_alpha_f(poke=False)
    assert f_expect > 1.25
    assert f1 >= 1.25, f"ele1 alpha stale after stage flip (f={f1:.4f})"
    assert f2 >= 1.25, f"ele2 alpha stale after stage flip (f={f2:.4f})"
