"""LadrunoConcrete3D (CDPM2-grade solid concrete, classTag 33017) — Zone-A openseespy battery.

This is the OpenSees material-in-an-element gate (the numpy oracle + the g++ kernel byte-check
live in tests/test_ladrunoConcrete3D_material.py). A single ``stdBrick`` unit cube under
statically-determinate 1/8-symmetry restraints develops a homogeneous uniaxial-STRESS state with
free lateral contraction; node 2 dof 1 is driven by displacement control and the centroid Gauss
point is read as a single-point constitutive probe.

The consistent tangent is NON-SYMMETRIC (non-associated flow + the semi-implicit theta-freeze), so
the analysis uses ``system FullGeneral``.

Gates (the element-level shadow of the oracle D0/D1/C1 + the wiring):
  * elastic         : below onset (eps0 = ft/E) sigma_xx = E eps_xx
  * tension peak    : nominal sigma_xx peaks at ~ft then SOFTENS (the damage peak; the effective
                      stress is monotonic) with omega_t -> high
  * compression peak: nominal |sigma_xx| peaks at ~fc (the MW failure surface at kappa_p=1)
  * damage routing  : in tension softening the effective stress exceeds the (damaged) nominal
  * response wiring : stress/strain/tangent/effectiveStress/damage/kappaP/plasticStrain reachable
  * state cycle     : FE_Datastore round-trip recovers the committed state

Plan/ADR: Ladruno_implementation/31_ladruno_concrete3d_adr.md.
"""
import math
import os
import sys

import pytest

from _testbed import ops
from _testbed.roundtrip import setresponse_smoke, database_roundtrip

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "_testbed"))
import concrete3d_ref as ref  # noqa: E402  the numpy oracle (the verified spec)

pytestmark = [pytest.mark.zone_a]

_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}

# oracle-gate concrete (E,nu,fc,ft,Gf,Gc); compression is NEGATIVE, fc/ft POSITIVE magnitudes
_E, _NU, _FC, _FT, _GF, _GC = 30000.0, 0.2, 30.0, 3.0, 0.1, 5.0


def _mat(tag, **kw):
    """nDMaterial LadrunoConcrete3D with the oracle-gate params (overridable via kwargs)."""
    E = kw.get("E", _E); nu = kw.get("nu", _NU)
    fc = kw.get("fc", _FC); ft = kw.get("ft", _FT)
    Gf = kw.get("Gf", _GF); Gc = kw.get("Gc", _GC)
    args = ["LadrunoConcrete3D", tag, E, nu, fc, ft, Gf, Gc]
    if "rho" in kw:
        args += ["-rho", kw["rho"]]
    if "lch" in kw:
        args += ["-lch", kw["lch"]]
    if kw.get("implex"):
        args += ["-implex"]
    if "eta" in kw:
        args += ["-eta", kw["eta"]]
    if "ctTemper" in kw:
        args += ["-ctTemper", kw["ctTemper"]]
    if "hoop" in kw:
        args += ["-hoop", kw["hoop"]]
    if "hoopFy" in kw:
        args += ["-hoopFy", kw["hoopFy"]]
    ops.nDMaterial(*args)


def _build(mat_fn, solver="FullGeneral"):
    """Single stdBrick unit cube, 1/8-symmetry determinate restraints. UNSYMMETRIC solver by default
    (the CDPM2 tangent is non-symmetric); pass solver="ProfileSPD" for the IMPL-EX symmetric-solver test."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *c)
    mat_fn(1)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 0, 1, 1)
    ops.fix(3, 0, 0, 1)
    ops.fix(4, 1, 0, 1)
    ops.fix(5, 1, 1, 0)
    ops.fix(6, 0, 1, 0)
    ops.fix(8, 1, 0, 0)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (2, 3, 6, 7):
        ops.load(n, 0.25, 0.0, 0.0)        # unit x-traction reference (DisplacementControl drives)
    ops.system(solver)                     # the CDPM2 tangent is NON-SYMMETRIC (Tier-1)
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-8, 200, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


def _sig_xx(gp=1):
    ops.eleResponse(1, "forces")                          # trigger the lazy strain push
    return list(ops.eleResponse(1, "stresses"))[(gp - 1) * 6]


def _wt(gp=1):
    return list(ops.eleResponse(1, "material", gp, "damage"))[0]


def _run(mat_fn, eps_target, nsteps):
    """Drive node-2 dof-1 displacement (=eps_xx on the unit cube) and return [(eps, sig_xx)]."""
    _build(mat_fn)
    out = []
    d = eps_target / nsteps
    ops.integrator("DisplacementControl", 2, 1, d)
    for _ in range(nsteps):
        if ops.analyze(1) != 0:
            break                                         # softening may stall; keep what converged
        out.append((ops.nodeDisp(2, 1), _sig_xx()))
    return out


def _drive_adaptive(mat_fn, eps_target, base_steps, max_cuts=7, solver="FullGeneral"):
    """Displacement-control driver with step-CUTTING through the softening limit point — the only
    way a single implicit element gets past an unconfined tension/compression peak (the snap-back
    regime). Returns [(eps_xx, sig_xx, omega_t)] for every converged increment."""
    _build(mat_fn, solver=solver)
    out = []
    step = eps_target / base_steps
    cuts = 0
    guard = 0
    while abs(ops.nodeDisp(2, 1)) < abs(eps_target) and guard < base_steps * 40:
        guard += 1
        ops.integrator("DisplacementControl", 2, 1, step)
        if ops.analyze(1) == 0:
            out.append((ops.nodeDisp(2, 1), _sig_xx(), _wt()))
            if cuts > 0 and guard % 4 == 0:               # tentatively grow the step back up
                step *= 2.0; cuts -= 1
        else:
            step *= 0.5; cuts += 1
            if cuts > max_cuts:
                break                                     # genuinely stuck — keep what converged
    return out


# ---------------------------------------------------------------------------
# elastic — below the damage onset eps0 = ft/E, the response is linear elastic
# ---------------------------------------------------------------------------
@pytest.mark.t0m
def test_elastic_uniaxial():
    eps = 0.2 * (_FT / _E)            # 20% of the tensile onset strain — safely elastic
    res = _run(lambda t: _mat(t), eps, 1)
    e, s = res[-1]
    assert s == pytest.approx(_E * e, rel=1e-5), f"elastic sigma {s} != E eps {_E * e}"


# ---------------------------------------------------------------------------
# tension — the nominal stress PEAKS at ~ft then softens (the damage peak)
# ---------------------------------------------------------------------------
@pytest.mark.t1
def test_uniaxial_tension_peak_and_softening():
    # tension damage onset is immediate (eps0 = ft/E ~ 1e-4) so the analysis hits the softening
    # limit point at once — drive with step-CUTTING through it (the snap-back regime; cf. ADR's
    # "unconfined softening is snap-backy -> Tier-3 explicit" note).
    res = _drive_adaptive(lambda t: _mat(t), 0.03, 300)
    assert len(res) > 12, f"only {len(res)} steps converged"
    sig = [s for _, s, _ in res]
    peak = max(sig)
    ipk = sig.index(peak)
    assert peak == pytest.approx(_FT, rel=0.05), f"tensile peak {peak} != ft {_FT}"
    assert ipk > 0, "peak at the very first step (no pre-peak rise captured)"
    # post-peak: the stress degrades and omega_t develops past the limit point
    assert sig[-1] < peak, f"did not degrade past peak (end {sig[-1]} vs peak {peak})"
    wt_max = max(w for _, _, w in res)
    assert wt_max > 0.1, f"omega_t {wt_max} did not develop in tension softening"


# ---------------------------------------------------------------------------
# compression — the nominal stress peaks at ~fc (the MW failure surface at kappa_p=1)
# ---------------------------------------------------------------------------
@pytest.mark.t1
def test_uniaxial_compression_peak():
    res = _run(lambda t: _mat(t), -0.006, 300)
    assert len(res) > 50
    mag = [abs(s) for _, s in res]
    peak = max(mag)
    assert peak == pytest.approx(_FC, rel=0.08), f"compressive peak {peak} != fc {_FC}"
    # effective-plasticity is monotonic (no plastic peak); the nominal peak is the DAMAGE peak
    assert mag.index(peak) > 0


# ---------------------------------------------------------------------------
# damage routing — in tension softening the (undamaged) effective stress exceeds the nominal
# ---------------------------------------------------------------------------
@pytest.mark.t1
def test_effective_exceeds_nominal_in_tension_softening():
    res = _drive_adaptive(lambda t: _mat(t), 0.03, 300)
    assert len(res) > 12, f"only {len(res)} steps converged"
    wt_max = max(w for _, _, w in res)
    assert wt_max > 0.05, f"expected tensile damage, omega_t_max={wt_max}"
    # at the last converged (damaged) state the undamaged effective stress exceeds the nominal
    nominal = _sig_xx()
    seff = list(ops.eleResponse(1, "material", 1, "effectiveStress"))[0]
    assert seff > nominal + 1e-6, f"effective {seff} should exceed damaged nominal {nominal}"


# ---------------------------------------------------------------------------
# response wiring + state round-trip
# ---------------------------------------------------------------------------
@pytest.mark.t1
def test_response_wiring():
    _run(lambda t: _mat(t), 0.01, 100)
    setresponse_smoke(1, [
        ("material", 1, "stress"),
        ("material", 1, "strain"),
        ("material", 1, "tangent"),
        ("material", 1, "effectiveStress"),
        ("material", 1, "damage"),
        ("material", 1, "kappaP"),
        ("material", 1, "plasticStrain"),
    ])


@pytest.mark.t1
def test_database_roundtrip():
    def build():
        _build(lambda t: _mat(t))
        ops.integrator("DisplacementControl", 2, 1, 0.5 * (_FT / _E))
        assert ops.analyze(1) == 0
    database_roundtrip(build, probe_nodes=[2], ndf=3,
                       probe_fn=lambda: list(ops.eleResponse(1, "material", 1, "stress")))


# ---------------------------------------------------------------------------
# P2g — CYCLIC no-heal: damage is MONOTONE on elastic unload (omega does NOT heal) and the nominal
# stress unloads along the degraded damage secant. This is the end-to-end check of the wrapper commit
# cycle carrying the new monotone drive history (sigtmax/sigcmax) under a real DisplacementControl
# analysis (the material/g++ levels are gated by test_p2g_monotone_damage_gate + the dmg_cyclic_unload
# byte-check). Tension: damage onset is immediate, so load with step-cutting through the limit point,
# then elastically unload (stable, plain DisplacementControl).
# ---------------------------------------------------------------------------
@pytest.mark.t1
def test_cyclic_no_heal_unload():
    _build(lambda t: _mat(t))
    out = []
    eps_peak = 0.0015
    step = eps_peak / 200; cuts = 0; guard = 0
    while ops.nodeDisp(2, 1) < eps_peak and guard < 2000:
        guard += 1
        ops.integrator("DisplacementControl", 2, 1, step)
        if ops.analyze(1) == 0:
            out.append((ops.nodeDisp(2, 1), _sig_xx(), _wt()))
            if cuts > 0 and guard % 4 == 0:
                step *= 2.0; cuts -= 1
        else:
            step *= 0.5; cuts += 1
            if cuts > 10:
                break
    assert len(out) > 10, f"tension load did not converge ({len(out)} steps)"
    eps_rev, _, wt_rev = out[-1]
    assert wt_rev > 0.2, f"omega_t at reversal {wt_rev} too small to test no-heal"

    # --- ELASTIC UNLOAD back toward half the reversal strain (plain DisplacementControl) ---
    nun = 60
    ops.integrator("DisplacementControl", 2, 1, (0.5 * eps_rev - eps_rev) / nun)
    unl = []
    for _ in range(nun):
        if ops.analyze(1) != 0:
            break
        unl.append((ops.nodeDisp(2, 1), _sig_xx(), _wt()))
    assert len(unl) > 10, f"unload did not converge ({len(unl)} steps)"

    # (1) NO HEALING — omega_t is FROZEN across the whole elastic unload (== the reversal value).
    wt_unl = [w for _, _, w in unl]
    assert min(wt_unl) >= wt_rev - 1.0e-9, f"omega_t healed on unload: {min(wt_unl)} < {wt_rev}"
    assert max(wt_unl) <= wt_rev + 1.0e-9, f"omega_t grew on elastic unload: {max(wt_unl)} > {wt_rev}"

    # (2) SECANT UNLOAD — the nominal stress drops toward 0 along the degraded secant; the unload
    #     stiffness d(sig)/d(eps) is positive, below the elastic E, and ~ (1-omega_t)*E (uniaxial stress).
    e0, s0, _ = unl[0]; e1, s1, _ = unl[-1]
    assert s1 < s0, f"stress did not reduce on unload ({s1} !< {s0})"
    sec_stiff = (s0 - s1) / (e0 - e1)
    assert 0.0 < sec_stiff < _E, f"unload secant stiffness {sec_stiff} not in (0, E={_E})"
    assert sec_stiff == pytest.approx((1.0 - wt_rev) * _E, rel=0.3), \
        f"secant stiffness {sec_stiff} != (1-wt)E {(1.0 - wt_rev) * _E}"


# ---------------------------------------------------------------------------
# P2h — the -ctTemper {none|alphat|proj} compression->tension damage-coupling modes parse, run, and
# round-trip. In PURE tension none and alphat are byte-identical (alpha_c=0 => w_t=1) and proj is lightly
# softer; all peak at ~ft. The discriminating compression->tension RESTORATION is gated at the material
# level (test_p2h_cttemper_gate + the g++ dmg_cttemper_* byte-checks); here we confirm the wrapper
# plumbing (parser + commit cycle + serialization) for every mode end-to-end.
# ---------------------------------------------------------------------------
@pytest.mark.t1
@pytest.mark.parametrize("mode", ["none", "alphat", "proj"])
def test_cttemper_parses_and_runs(mode):
    res = _drive_adaptive(lambda t: _mat(t, ctTemper=mode), 0.02, 300)
    assert len(res) > 12, f"{mode}: only {len(res)} steps converged"
    peak = max(s for _, s, _ in res)
    assert peak == pytest.approx(_FT, rel=0.05), f"{mode} tension peak {peak} != ft {_FT}"


@pytest.mark.t1
def test_cttemper_database_roundtrip():
    """A -ctTemper mode survives a FE_Datastore round-trip (the ctTemper int is serialized)."""
    def build():
        _build(lambda t: _mat(t, ctTemper="alphat"))
        ops.integrator("DisplacementControl", 2, 1, 0.5 * (_FT / _E))
        assert ops.analyze(1) == 0
    database_roundtrip(build, probe_nodes=[2], ndf=3,
                       probe_fn=lambda: list(ops.eleResponse(1, "material", 1, "stress")))


# ---------------------------------------------------------------------------
# P3 Tier-3 — EXPLICIT dynamics. Tier-1 implicit STALLS at the unconfined tension-softening limit point
# (the snap-back; test_uniaxial_tension_peak_and_softening needs step-cutting to crawl past it) because the
# consistent tangent goes INDEFINITE there (gate TD2). An EXPLICIT integrator advances on the MASS matrix
# and NEVER factorizes that tangent, so the same material runs straight through tensile cracking. The
# committed stress is byte-identical to Tier-1 (g++ B7: t3=0 — do_tangent=false == do_tangent=true). Here we
# DEMONSTRATE the explicit PATH end-to-end: a free-dynamics tension yank (initial face velocity — the proven
# explicit idiom, NOT a prescribed-SP ramp which is an implicit/static construct) drives the cube into the
# damaging regime under CentralDifferenceLadruno; the run COMPLETES (no stall — the Tier-3 value vs the
# implicit limit-point stall), the stress stays BOUNDED (no instability), and omega_t develops (the inelastic
# regime was integrated under explicit). rho is picked for a convenient explicit dt. (NB the elastic
# strain-to-peak eps0=ft/E~1e-4 is a blink under free dynamics, so this is a "runs through cracking" demo, not
# a clean quasi-static loading curve — that lives in the numpy oracle + the g++ B7 correctness gate.)
# ---------------------------------------------------------------------------
@pytest.mark.t1
def test_tier3_explicit_path_runs():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *c)
    _mat(1, rho=2.4e-3)                           # -rho => element mass for the explicit (Diagonal) solve
    ops.fix(1, 1, 1, 1); ops.fix(2, 0, 1, 1); ops.fix(3, 0, 0, 1); ops.fix(4, 1, 0, 1)
    ops.fix(5, 1, 1, 0); ops.fix(6, 0, 1, 0); ops.fix(8, 1, 0, 0)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    for n in (2, 3, 6, 7):                        # initial +x velocity on the x=1 face => free tension stretch
        ops.setNodeVel(n, 1, 0.4, "-commit")
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")                        # lumped-mass explicit solve (never factorizes the tangent)
    ops.test("NormDispIncr", 1.0e-12, 1)
    ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno")    # auto dt_cr; prime once then a sub-critical fixed dt
    ops.analysis("Transient")
    assert ops.analyze(1, 1.0e-9) == 0, "explicit prime step failed"
    dt = 0.3 * ops.criticalTimeStep()
    sig, done = [], 0
    for _ in range(400):
        if ops.analyze(1, dt) != 0:
            break
        done += 1
        sig.append(_sig_xx())
    # (1) COMPLETES — the material integrates under explicit (mass via getRho + transient + getStress),
    #     the path that makes softening safe (the global solve never factorizes the indefinite tangent).
    assert done == 400, f"explicit run stalled at step {done}/400"
    # (2) BOUNDED + LOADED — the free tension stretch develops a sensible stress (no instability/blow-up);
    #     the inelastic-softening correctness in explicit mode (do_tangent-independence) is the g++ B7 gate.
    assert all(abs(s) < 5.0 * _FT for s in sig), "stress blew up (explicit instability)"
    assert max(abs(s) for s in sig) > 0.2 * _FT, "material did not load under explicit"


# ---------------------------------------------------------------------------
# MULTIAXIAL / SHEAR coverage — pins the wrapper's engineering<->tensor shear conversions
# (strain x0.5 in, x2 out; tangent shear-COLUMNS x0.5; stress unscaled) against the numpy oracle.
# The element gates above drive only node-2 dof-1, so every engineering shear component is identically
# zero and a wrong shear factor would stay green (the adversarial-review false-green gap). These probe
# the material OBJECT directly with arbitrary 6-strain vectors (incl. shear) via the `NDTest` facility
# (SetStrain -> setTrialStrain; GetStress -> getStress; GetTangentStiffness -> getTangent), so the
# conversion is exercised at non-zero shear with no element/solve in the way.
# ---------------------------------------------------------------------------
def _matobj(**kw):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _mat(1, **kw)
    return 1


def _nd_stress(tag, eps_eng):
    ops.NDTest("SetStrain", tag, *[float(x) for x in eps_eng])
    return list(ops.NDTest("GetStress", tag))


def _nd_tangent(tag, eps_eng):
    ops.NDTest("SetStrain", tag, *[float(x) for x in eps_eng])
    return list(ops.NDTest("GetTangentStiffness", tag))   # 36, row-major


def _eps_tensor(eps_eng):
    """engineering {.., gxy, gyz, gzx} -> tensor Voigt {.., exy, eyz, ezx} (the oracle convention)."""
    return ref.np.array([eps_eng[0], eps_eng[1], eps_eng[2],
                         0.5 * eps_eng[3], 0.5 * eps_eng[4], 0.5 * eps_eng[5]])


@pytest.mark.t1
def test_elastic_simple_shear():
    """Pure elastic simple shear gamma_xy: sig_xy = G*gamma_xy and the normal stresses ~ 0 — the
    decisive check of the shear factor. strain x0.5 IN (eps_xy_tensor = gamma/2) and stress shear
    UNscaled (sig_xy = 2G*eps_xy_tensor = G*gamma)."""
    G = _E / (2.0 * (1.0 + _NU))
    gamma = 0.3 * (_FT / _E)                            # safely elastic
    tag = _matobj()
    sig = _nd_stress(tag, [0, 0, 0, gamma, 0, 0])
    assert sig[3] == pytest.approx(G * gamma, rel=1e-5), f"sig_xy {sig[3]} != G*gamma {G * gamma}"
    for i in (0, 1, 2, 4, 5):
        assert abs(sig[i]) < 1e-6 * G * gamma + 1e-9, f"spurious sig[{i}] = {sig[i]}"
    # getStrain round-trips engineering shear (x2 out): SetStrain gamma -> GetStrain gamma
    eps_back = list(ops.NDTest("GetStrain", tag))
    assert eps_back[3] == pytest.approx(gamma, rel=1e-6), f"getStrain shear {eps_back[3]} != {gamma}"


@pytest.mark.t1
def test_multiaxial_elastic_stress_vs_oracle():
    """A small multiaxial strain WITH shear, elastic: all 6 stress comps == the oracle C0:eps (true
    tensor stress). Pins strain-in shear x0.5 + Poisson coupling + stress-out unscaled."""
    mp = ref.make_material(_E, _NU, _FC, _FT)
    s = _FT / _E
    eps_eng = [0.2 * s, -0.1 * s, 0.05 * s, 0.25 * s, -0.15 * s, 0.1 * s]
    sig = _nd_stress(_matobj(), eps_eng)
    sig_ref = ref.elastic_C(mp) @ _eps_tensor(eps_eng)   # tensor-convention C @ tensor strain = true stress
    for i in range(6):
        assert sig[i] == pytest.approx(float(sig_ref[i]), rel=1e-5, abs=1e-7), \
            f"stress[{i}] {sig[i]} != oracle {sig_ref[i]}"


@pytest.mark.t1
def test_elastic_tangent_vs_oracle():
    """The material's 6x6 tangent (wrapper engineering convention) == the oracle elastic C0 with shear
    COLUMNS halved: T[i][j] = C0[i][j] (j<3) or 0.5*C0[i][j] (j>=3). Pins getTangent's column-only
    x0.5 + the lambda (Poisson) coupling + the 2G->G shear stiffness — untouched by the uniaxial axial
    asserts. (Elastic state, so getTangent == the elastic operator.)"""
    mp = ref.make_material(_E, _NU, _FC, _FT)
    s = _FT / _E
    T = _nd_tangent(_matobj(), [0.15 * s, -0.05 * s, 0.0, 0.2 * s, 0.0, 0.1 * s])
    C0 = ref.elastic_C(mp)
    cmax = max(abs(float(C0[i][j])) for i in range(6) for j in range(6))
    for i in range(6):
        for j in range(6):
            expect = float(C0[i][j]) * (0.5 if j >= 3 else 1.0)
            assert T[i * 6 + j] == pytest.approx(expect, rel=1e-5, abs=1e-6 * cmax), \
                f"tangent[{i}][{j}] {T[i*6+j]} != expected {expect} (col-halved oracle C0)"


# ===========================================================================
# PHASE-2 REDUCED VIEWS — PlaneStrain / PlaneStress / AxiSymmetric / PlateFiber
#
# One class serves every dimensional view via a `dim` mode (the LadrunoJ2 doctrine,
# PR #90): the kernel return map always runs on the full 6-comp tensor; the element
# vector is mapped to the reduced ordering, and PlaneStress/PlateFiber enforce
# sigma_zz = 0 by a nested Newton on eps_zz + static condensation of the 33 dof.
#
# CDPM2 has NO upstream OpenSees peer, so the gates are SELF-REFERENTIAL: a reduced
# view created by an element's getCopy(type) is checked against the shipped 3D
# material under the matching out-of-plane constraint —
#   PlaneStrain : in-plane stress/tangent == the 3D slice with eps_zz = gxz = gyz = 0
#   PlaneStress : in-plane stress/tangent == the 3D state with eps_zz solved so sig_zz = 0
# The elastic tangents are additionally pinned to the CLOSED-FORM plane-strain /
# plane-stress isotropic moduli (independent of the wrapper). NDTest can only feed a
# 6-comp strain (3D), so reduced views are reached through the `quad` element
# (PlaneStrain/PlaneStress) and `bbarQuad` (AxiSymmetric2D), exactly as #90 did.
# AxiSym shares the PlaneStrain mapping machinery; PlateFiber shares the PlaneStress
# condensation machinery (so the two gated archetypes cover all four views' code).
# ===========================================================================
_QNODES = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0), 4: (0.0, 1.0)}


def _build_quad(qtype):
    """A single unit quad of the reduced-view material; the element calls getCopy(qtype)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for t, (x, y) in _QNODES.items():
        ops.node(t, x, y)
    _mat(1)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.fix(4, 1, 0)
    ops.element("quad", 1, 1, 2, 3, 4, 1.0, qtype, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.system("FullGeneral")                                  # CDPM2 tangent is NON-SYMMETRIC
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-10, 100, 0)
    ops.algorithm("Newton")


def _gp1_strain():
    return list(ops.eleResponse(1, "strain"))[0:3]             # (exx, eyy, gxy) engineering, GP1


def _gp1_stress():
    return list(ops.eleResponse(1, "stress"))[0:3]             # (sxx, syy, sxy), GP1


def _gp1_mat_tangent():
    return list(ops.eleResponse(1, "material", 1, "tangent"))  # reduced ncomp x ncomp, row-major


def _solve_elastic_quad(qtype, fx=0.4, fy=0.3):
    """Drive a quad into a small (sub-onset) multiaxial+shear elastic state and read GP1."""
    _build_quad(qtype)
    ops.load(2, fx, 0.0)
    ops.load(3, fx, fy)                                         # mixed -> in-plane shear
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, f"{qtype} elastic quad did not converge"
    return _gp1_strain(), _gp1_stress(), _gp1_mat_tangent()


def _planestrain_closed_form_D():
    E, nu = _E, _NU
    k = E / ((1.0 + nu) * (1.0 - 2.0 * nu))
    return [[k * (1.0 - nu), k * nu, 0.0],
            [k * nu, k * (1.0 - nu), 0.0],
            [0.0, 0.0, k * (1.0 - 2.0 * nu) / 2.0]]            # sig = D * eps_eng (gamma)


def _planestress_closed_form_D():
    E, nu = _E, _NU
    k = E / (1.0 - nu * nu)
    return [[k, k * nu, 0.0],
            [k * nu, k, 0.0],
            [0.0, 0.0, k * (1.0 - nu) / 2.0]]


# ---------------------------------------------------------------------------
# PlaneStrain — the MAPPING path (no condensation): the reduced view is exactly
# the 3D slice with eps_zz = gxz = gyz = 0.
# ---------------------------------------------------------------------------
@pytest.mark.t1
def test_planestrain_stress_matches_3d_slice():
    eps, sig, _ = _solve_elastic_quad("PlaneStrain")
    exx, eyy, gxy = eps
    assert max(abs(exx), abs(eyy), abs(gxy)) < _FT / _E, "drove past the damage onset (not elastic)"
    # the shipped 3D material fed the same in-plane strain with eps_zz = gxz = gyz = 0
    sig3d = _nd_stress(_matobj(), [exx, eyy, 0.0, gxy, 0.0, 0.0])
    for k, full in enumerate((0, 1, 3)):
        assert sig[k] == pytest.approx(sig3d[full], rel=1e-6, abs=1e-7), \
            f"PlaneStrain stress[{k}] {sig[k]} != 3D slice {sig3d[full]}"


@pytest.mark.t1
def test_planestrain_tangent_closed_form():
    _, _, T = _solve_elastic_quad("PlaneStrain")
    D = _planestrain_closed_form_D()
    dmax = max(abs(D[i][j]) for i in range(3) for j in range(3))
    for i in range(3):
        for j in range(3):
            assert T[i * 3 + j] == pytest.approx(D[i][j], rel=1e-5, abs=1e-6 * dmax), \
                f"PlaneStrain tangent[{i}][{j}] {T[i*3+j]} != closed form {D[i][j]}"


# ---------------------------------------------------------------------------
# PlaneStress — the CONDENSATION path: eps_zz solved so sig_zz = 0.
# ---------------------------------------------------------------------------
def _solve_szz_zero_3d(tag, exx, eyy, gxy, ezz0=0.0, commit=False):
    """Newton on eps_zz so the 3D material's sigma_zz = 0; returns (in-plane stress, eps_zz)."""
    ezz = ezz0
    s = None
    for _ in range(60):
        s = _nd_stress(tag, [exx, eyy, ezz, gxy, 0.0, 0.0])
        T = _nd_tangent(tag, [exx, eyy, ezz, gxy, 0.0, 0.0])    # current 6x6 (row-major)
        smag = sum(v * v for v in s) ** 0.5
        if abs(s[2]) <= 1e-11 * (smag if smag > 1.0 else 1.0):
            break
        d22 = T[2 * 6 + 2]                                       # d(sig_zz)/d(eps_zz), normal col unscaled
        ezz -= s[2] / d22
    if commit:
        ops.NDTest("CommitState", tag)
    return [s[0], s[1], s[3]], ezz


@pytest.mark.t1
def test_planestress_stress_enforces_szz_zero():
    eps, sig, _ = _solve_elastic_quad("PlaneStress")
    exx, eyy, gxy = eps
    assert max(abs(exx), abs(eyy), abs(gxy)) < _FT / _E, "drove past the damage onset (not elastic)"
    sig3d, _ = _solve_szz_zero_3d(_matobj(), exx, eyy, gxy)
    for k in range(3):
        assert sig[k] == pytest.approx(sig3d[k], rel=1e-6, abs=1e-7), \
            f"PlaneStress stress[{k}] {sig[k]} != sigma_zz=0 3D solve {sig3d[k]}"


@pytest.mark.t1
def test_planestress_tangent_closed_form():
    _, _, T = _solve_elastic_quad("PlaneStress")
    D = _planestress_closed_form_D()
    dmax = max(abs(D[i][j]) for i in range(3) for j in range(3))
    for i in range(3):
        for j in range(3):
            assert T[i * 3 + j] == pytest.approx(D[i][j], rel=1e-5, abs=1e-6 * dmax), \
                f"PlaneStress tangent[{i}][{j}] {T[i*3+j]} != closed form {D[i][j]}"


# ---------------------------------------------------------------------------
# PlaneStress NONLINEAR — the condensation must hold under plasticity + damage.
# Drive the quad in (pre-peak, hardening) compression so it converges robustly,
# record the GP1 strain path, then REPLAY it through the 3D material with the
# sigma_zz = 0 Newton + commit each step: the final in-plane stress must match.
# This validates the condensed return map (state + tangent) end-to-end.
# ---------------------------------------------------------------------------
@pytest.mark.t1
def test_planestress_nonlinear_replay_matches_3d():
    # uniaxial in-plane COMPRESSION (hardening pre-peak ⇒ robust convergence), driven by
    # DisplacementControl on the right edge (nodes 2,3 tied in x via equalDOF) so the GP
    # state is a clean, well-controlled uniaxial-stress compression that goes plastic+damaged.
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for t, (x, y) in _QNODES.items():
        ops.node(t, x, y)
    _mat(1)
    ops.fix(1, 1, 1)                                            # origin pinned
    ops.fix(2, 0, 1)                                            # bottom-right: y fixed
    ops.fix(4, 1, 0)                                            # top-left: x fixed
    ops.equalDOF(2, 3, 1)                                       # right edge moves uniformly in x
    ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStress", 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, -1.0, 0.0)                                      # reference; DisplacementControl drives
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-9, 100, 0)
    ops.algorithm("Newton")
    nsteps = 50
    eps_target = -1.2e-3                                        # past onset, pre-peak (converges)
    ops.integrator("DisplacementControl", 2, 1, eps_target / nsteps)
    ops.analysis("Static")
    path = []
    for _ in range(nsteps):
        assert ops.analyze(1) == 0, "PlaneStress compression quad did not converge"
        path.append(_gp1_strain())
    sig_quad = _gp1_stress()
    # must have actually gone nonlinear (guard against a vacuous elastic comparison)
    assert max(abs(c) for c in path[-1]) > _FT / _E, "stayed elastic — replay would be trivial"

    # replay the SAME GP1 strain path through a fresh 3D material, committing each step
    tag = _matobj()
    ezz = 0.0
    sig3d = None
    for (exx, eyy, gxy) in path:
        sig3d, ezz = _solve_szz_zero_3d(tag, exx, eyy, gxy, ezz0=ezz, commit=True)
    for k in range(3):
        assert sig_quad[k] == pytest.approx(sig3d[k], rel=2e-5, abs=1e-6), \
            f"nonlinear PlaneStress stress[{k}] {sig_quad[k]} != replayed 3D {sig3d[k]}"


# ---------------------------------------------------------------------------
# AxiSymmetric — getCopy("AxiSymmetric2D") via bbarQuad. The order-4 mapping shares
# the (passing) PlaneStrain non-condense machinery; the only AxiSym-specific datum is
# vmap={00,11,22(hoop),01} — confirmed against the element's own documented strain
# order (ConstantPressureVolumeQuad.cpp:367-370). We verify the getCopy path + order-4
# + the tangent mapping SOLVE-FREE: building the element calls getCopy("AxiSymmetric2D")
# (it would error if it returned null/wrong order), and the material's elastic 4×4
# tangent (read via eleResponse, no analysis) must equal the closed-form axisymmetric
# isotropic operator. (A full bbarQuad SOLVE is NOT gated: the mixed u-p element is
# numerically fragile with the non-symmetric CDPM2 tangent — an element-interaction
# issue, not a reduced-view bug; cf. the snap-back caveat.)
# ---------------------------------------------------------------------------
def _axisym_closed_form_C():
    E, nu = _E, _NU
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    mu = E / (2.0 * (1.0 + nu))
    a = lam + 2.0 * mu
    return [[a, lam, lam, 0.0],
            [lam, a, lam, 0.0],
            [lam, lam, a, 0.0],
            [0.0, 0.0, 0.0, mu]]                               # rz shear column halved (2mu->mu)


@pytest.mark.t1
def test_axisymmetric_getcopy_and_tangent():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 1.0, 0.0); ops.node(2, 2.0, 0.0)
    ops.node(3, 2.0, 1.0); ops.node(4, 1.0, 1.0)               # off the symmetry axis (r > 0)
    _mat(1)
    ops.element("bbarQuad", 1, 1, 2, 3, 4, 1.0, 1)             # getCopy("AxiSymmetric2D") -> order 4
    T = list(ops.eleResponse(1, "material", 1, "tangent"))     # elastic 4x4, no solve
    assert len(T) == 16, f"AxiSym tangent has {len(T)} entries, expected 16 (order 4)"
    C = _axisym_closed_form_C()
    cmax = max(abs(C[i][j]) for i in range(4) for j in range(4))
    for i in range(4):
        for j in range(4):
            assert T[i * 4 + j] == pytest.approx(C[i][j], rel=1e-5, abs=1e-6 * cmax), \
                f"AxiSym tangent[{i}][{j}] {T[i*4+j]} != closed form {C[i][j]}"


# ===========================================================================
# P3 Tier-2 IMPL-EX (`-implex`) — element-level plumbing. The IMPL-EX algorithm + its SPD/consistency
# falsification battery live in the numpy oracle (test_p3_implex_gate) and the g++ kernel byte-check
# (B5, reported stress == oracle to ~1e-16). These gate the WRAPPER wiring end-to-end: the parser flag,
# integrate() passing dt=ops_Dt + reading the IMPL-EX state, commitState, and send/recvSelf of the new
# committed fields. Plus the user-facing payoff: a single-sign softening run on a SYMMETRIC solver.
# ===========================================================================
def _run_loadcontrol(mat_fn, fnode, nsteps):
    """Uniform LoadControl compression (the IMPL-EX-friendly stepping — a monotone pseudo-time, unlike
    DisplacementControl's lambda-based dt near a limit point). Returns (converged_steps, sig_xx, kappaP)."""
    _build(mat_fn)
    ops.remove("loadPattern", 1)
    ops.pattern("Plain", 1, 1)
    for n in (2, 3, 6, 7):
        ops.load(n, fnode, 0.0, 0.0)
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    ok = 0
    for _ in range(nsteps):
        if ops.analyze(1) != 0:
            break
        ok += 1
    s = _sig_xx()
    kp = list(ops.eleResponse(1, "material", 1, "kappaP"))[0]
    return ok, s, kp


@pytest.mark.t1
def test_implex_loadcontrol_matches_tier1():
    """-implex end-to-end under uniform LoadControl compression into the HARDENING regime: converges,
    develops plasticity (kappaP>0), and the committed stress matches the Tier-1 run (pre-peak omega_c=0,
    so the explicit/implicit gap is the plastic-flow extrapolation only — tight). Exercises the parser,
    integrate()+dt=ops_Dt, the IMPL-EX commit, and confirms IMPL-EX does not corrupt the implicit state.
    (NB IMPL-EX wants ~uniform pseudo-time; DisplacementControl past a limit point feeds a non-uniform
    lambda-dt that over-extrapolates — use Tier-1 / transient / uniform stepping there.)"""
    o1, s1, k1 = _run_loadcontrol(lambda t: _mat(t), -5.0, 300)
    o2, s2, k2 = _run_loadcontrol(lambda t: _mat(t, implex=True), -5.0, 300)
    assert o1 == 300 and o2 == 300, f"did not fully converge (tier1 {o1}/300, implex {o2}/300)"
    assert k2 > 0.05, f"-implex did not go plastic (kappaP={k2})"
    assert s2 == pytest.approx(s1, rel=1e-4), f"-implex sig {s2} != Tier-1 {s1}"


@pytest.mark.t1
def test_implex_tangent_better_conditioned_in_softening():
    """The IMPL-EX payoff at the material level, via NDTest (strain-controlled ⇒ no global limit point).
    Drive uniaxial TENSION past onset into softening (SetStrain + CommitState each step) for a Tier-1 and
    an -implex material, then read the axial tangent C[0][0] (index 0 is a NORMAL component — unaffected
    by the engineering shear-column halving, so its sign is convention-robust): the Tier-1 damaged tangent
    is INDEFINITE (C[0][0] < 0, the snap-back) while the -implex degraded-elastic secant D_dam(w~):C0 is
    POSITIVE (C[0][0] > 0). Both report a softened axial stress."""
    def drive(implex):
        ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
        _mat(1, implex=implex)
        for e in [i * 4.0e-4 / 240 for i in range(1, 241)]:   # 0 -> 4e-4 (well past onset eps0=1e-4)
            ops.NDTest("SetStrain", 1, e, 0, 0, 0, 0, 0)
            ops.NDTest("CommitState", 1)
        sig = list(ops.NDTest("GetStress", 1))[0]
        T = list(ops.NDTest("GetTangentStiffness", 1))        # 36, row-major
        return sig, T[0]                                       # sig_xx, C[0][0]
    s_t1, c_t1 = drive(False)
    s_ix, c_ix = drive(True)
    # both softened well below the tensile peak ft
    assert s_t1 < _FT and s_ix < _FT, f"did not soften (tier1 {s_t1}, implex {s_ix})"
    # the decisive contrast: Tier-1 axial tangent INDEFINITE, IMPL-EX secant POSITIVE
    assert c_t1 < 0.0, f"expected Tier-1 softening axial tangent < 0, got {c_t1}"
    assert c_ix > 0.0, f"expected -implex secant axial tangent > 0, got {c_ix}"


@pytest.mark.t1
def test_implex_database_roundtrip():
    """FE_Datastore round-trip of a partially-damaged -implex material exercises the new send/recvSelf
    IMPL-EX committed fields (implex flag + wt/wc/dwt/dwc/dtn + depl[6])."""
    def build():
        _build(lambda t: _mat(t, implex=True))
        ops.integrator("DisplacementControl", 2, 1, 0.5 * (_FT / _E))
        assert ops.analyze(1) == 0
    database_roundtrip(build, probe_nodes=[2], ndf=3,
                       probe_fn=lambda: list(ops.eleResponse(1, "material", 1, "stress")))


# ===========================================================================
# P3 Duvaut-Lions viscoplasticity (`-eta`) — element-level plumbing. The relaxation MATH (the
# Simo-Hughes blend) + the closed-form 1-D overstress oracle live in the numpy oracle (test_p3_eta_gate)
# and the g++ kernel byte-check (B6: NETA, reported stress == oracle to ~1e-14, eta=0 == inviscid). These
# gate the WRAPPER wiring end-to-end: the parser flag, integrate() passing dt=ops_Dt into the kernel
# relaxation, commit, and send/recvSelf of the new eta field.
# ===========================================================================
def _run_loadcontrol_disp(mat_fn, fnode, nsteps):
    """Uniform LoadControl compression (monotone pseudo-time ⇒ dt=ops_Dt=1/nsteps, uniform — the
    well-defined Duvaut-Lions dt). Returns (converged_steps, nodeDisp(2,1), sig_xx, kappaP)."""
    _build(mat_fn)
    ops.remove("loadPattern", 1)
    ops.pattern("Plain", 1, 1)
    for n in (2, 3, 6, 7):
        ops.load(n, fnode, 0.0, 0.0)
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    ok = 0
    for _ in range(nsteps):
        if ops.analyze(1) != 0:
            break
        ok += 1
    return ok, ops.nodeDisp(2, 1), _sig_xx(), list(ops.eleResponse(1, "material", 1, "kappaP"))[0]


@pytest.mark.t1
def test_eta_loadcontrol_inviscid_limit_and_effect():
    """-eta end-to-end under uniform LoadControl compression into the HARDENING regime. Gates the WRAPPER
    wiring: the parser flag, integrate() reading dt=ops_Dt and applying the Duvaut-Lions relaxation, and
    commit. Under FORCE control the converged gauss stress is fixed by equilibrium (≈ load), so the
    relaxation shows in the STRAIN: a vanishing eta (eta << dt) recovers Tier-1, while a finite eta
    (eta ~ dt ⇒ less plastic flow per step ⇒ stiffer over the step) changes the strain at a given load."""
    o0, d0, s0, k0 = _run_loadcontrol_disp(lambda t: _mat(t), -5.0, 300)                 # Tier-1 inviscid
    oe, de, se, ke = _run_loadcontrol_disp(lambda t: _mat(t, eta=1.0e-9), -5.0, 300)     # eta << dt => inviscid
    ob, db, sb, kb = _run_loadcontrol_disp(lambda t: _mat(t, eta=0.05), -5.0, 300)       # eta >> dt => relaxed
    assert o0 == 300 and oe == 300 and ob == 300, f"did not fully converge ({o0},{oe},{ob})"
    assert k0 > 0.05, f"Tier-1 did not go plastic (kappaP={k0})"
    # inviscid limit: a vanishing eta recovers Tier-1 (strain AND stress)
    assert de == pytest.approx(d0, rel=1.0e-4), f"eta->0 disp {de} != Tier-1 {d0}"
    assert se == pytest.approx(s0, rel=1.0e-4), f"eta->0 sig {se} != Tier-1 {s0}"
    # genuine effect: a real eta changes the response (less plastic strain per step => smaller |disp|)
    assert abs(db - d0) > 1.0e-3 * abs(d0), f"-eta had no end-to-end effect (visc disp {db} vs tier1 {d0})"


@pytest.mark.t1
def test_eta_database_roundtrip():
    """FE_Datastore round-trip of a viscous (-eta) material exercises the new send/recvSelf eta field."""
    def build():
        _build(lambda t: _mat(t, eta=0.5))
        ops.integrator("DisplacementControl", 2, 1, 0.5 * (_FT / _E))
        assert ops.analyze(1) == 0
    database_roundtrip(build, probe_nodes=[2], ndf=3,
                       probe_fn=lambda: list(ops.eleResponse(1, "material", 1, "stress")))


# ===========================================================================
# Tier-3 EXPLICIT-DYNAMICS demo — the ADR §4.4 "softening is a non-issue" payoff.
#
# CDPM2 unconfined softening is snap-backy: the Tier-1 implicit path needs an
# UNSYMMETRIC solver (FullGeneral) AND displacement-control step-CUTTING through
# the limit point (see _drive_adaptive / test_uniaxial_tension_peak_and_softening).
# An EXPLICIT integrator sidesteps both: it never forms or factorizes a global
# tangent (so the indefinite/non-symmetric damaged tangent, TD2/TD3, is irrelevant)
# and a displacement-driven boundary makes the response kinematically controlled —
# there is NO global limit point to track. The material kernel already runs with
# doTangent=false; this section validates that end-to-end in OpenSees under BOTH
# fork explicit steppers, cross-checked against the SAME (oracle-verified) static
# softening backbone.
#
# Specimen: the 1/8-symmetry unit cube (a uniaxial-stress constitutive probe), now
# given mass via -rho. The driven face x-DOF is prescribed by a Linear timeSeries
# (support motion) so the axial strain is imposed kinematically; the free lateral
# (Poisson) DOFs carry the inertia. Loading is slow (hundreds of wave transits over
# the run) so the gauss stress tracks the quasi-static backbone.
# ===========================================================================
CDL = ("CentralDifferenceLadruno",)
EB = ("ExplicitBathe", 0.54)            # p = 0.54, the source-battery sub-step value
_RHO = 2.4                              # density (consistent units) -> finite element mass


def _build_explicit(mat_fn, rate, alphaM=0.0):
    """1/8-symmetry unit cube (same determinate restraints as _build) with MASS, the driven face
    x-DOF prescribed as u_x(t) = rate * t via a Linear timeSeries (support motion). system Diagonal
    + algorithm Linear => no global tangent is ever formed/factorized (the Tier-3 point)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *c)
    mat_fn(1)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 0, 1, 1)
    ops.fix(3, 0, 0, 1)
    ops.fix(4, 1, 0, 1)
    ops.fix(5, 1, 1, 0)
    ops.fix(6, 0, 1, 0)
    ops.fix(8, 1, 0, 0)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Linear", 1)                       # factor(t) = t
    ops.pattern("Plain", 1, 1)
    for n in (2, 3, 6, 7):
        ops.sp(n, 1, rate)                            # prescribed u_x = rate * t (the imposed strain ramp)
    if alphaM > 0.0:
        ops.rayleigh(alphaM, 0.0, 0.0, 0.0)           # mass-proportional only => stays tangent-free
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1.0e-12, 1)
    ops.algorithm("Linear")


def _run_explicit(integrator_args, eps_target, nsteps, dt, alphaM=0.0, **mat_kw):
    """Drive the cube to eps_target axial strain over nsteps explicit steps of size dt. Returns
    [(eps_xx, sig_xx, omega_t)] per step (eps_xx = nodeDisp(2,1) on the unit cube). Extra mat_kw
    (e.g. lch=...) are forwarded to the material."""
    t_end = nsteps * dt
    rate = eps_target / t_end                          # u_x(t) = rate * t reaches eps_target at t_end
    _build_explicit(lambda tag: _mat(tag, rho=_RHO, **mat_kw), rate, alphaM=alphaM)
    ops.integrator(*integrator_args)
    ops.analysis("Transient")
    out = []
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            break
        out.append((ops.nodeDisp(2, 1), _sig_xx(), _wt()))
    return out


# the explicit step: ~0.35 x the unit-cube dilatational CFL limit le/c_d, c_d = sqrt((K+4G/3)/rho).
# K=E/(3(1-2nu)), G=E/(2(1+nu)) for E=30000,nu=0.2 => c_d~118, le/c_d~0.0085 -> dt=0.003 is safe.
_DT_EXPL = 0.003
_NSTEPS_EXPL = 2000                                    # T = 6.0 -> ~700 wave transits (quasi-static)


@pytest.mark.t1
@pytest.mark.parametrize("integ", [CDL, EB], ids=["CDL", "ExplicitBathe"])
def test_explicit_tension_softening_runs_through(integ):
    """The headline Tier-3 gate: an explicit run marches a CDPM2 specimen straight THROUGH tension
    softening — no unsymmetric solver, no step-cutting, no convergence test — and stays finite. The
    nominal gauss stress peaks at ~ft then degrades while omega_t -> ~1, exactly the static backbone
    (test_uniaxial_tension_peak_and_softening) but reached by forward time-marching."""
    res = _run_explicit(integ, 0.03, _NSTEPS_EXPL, _DT_EXPL, alphaM=20.0)
    assert len(res) == _NSTEPS_EXPL, f"explicit run did not complete ({len(res)}/{_NSTEPS_EXPL} steps)"
    sig = [s for _, s, _ in res]
    assert all(math.isfinite(s) for s in sig), "explicit softening run went non-finite"
    peak = max(sig)
    ipk = sig.index(peak)
    assert peak == pytest.approx(_FT, rel=0.10), f"explicit tensile peak {peak} != ft {_FT}"
    assert ipk > 0, "peak at the very first step (no pre-peak elastic rise captured)"
    # softens substantially past the peak (oracle: ~1.22 at eps=0.03, i.e. ~0.4 ft)
    assert sig[-1] < 0.6 * peak, f"did not soften past peak (end {sig[-1]:.3f} vs peak {peak:.3f})"
    wt_max = max(w for _, _, w in res)
    assert wt_max > 0.5, f"omega_t {wt_max} did not develop in explicit tension softening"


def test_explicit_backbone_matches_oracle():
    """Cross-validation against the VERIFIED SPEC: the cube (default lch=1.0, no -autoRegularization)
    is a uniaxial-stress probe, so the explicit nominal-stress backbone must track the numpy oracle's
    drive_uniaxial_tension_damaged (lch=1.0) — pinning that mass/rate do not corrupt the rate-INDEPENDENT
    (eta=0) CDPM2 response. Compared at the PEAK (~ft) and at the final softened strain; tolerances
    absorb the explicit dynamic ripple about the quasi-static backbone."""
    import numpy as np
    mp = ref.make_material(_E, _NU, _FC, _FT)
    bb = ref.drive_uniaxial_tension_damaged(mp, np.linspace(0.0, 0.03, 600), Gf=_GF, lch=1.0)
    sig_o = bb["sig11"]
    peak_o = float(sig_o.max())
    end_o = float(sig_o[-1])                                       # ~1.22 at eps=0.03 (softened, omega_t~1)

    expl = _run_explicit(CDL, 0.03, _NSTEPS_EXPL, _DT_EXPL, alphaM=20.0)
    assert len(expl) == _NSTEPS_EXPL
    sig_e = [s for _, s, _ in expl]
    assert max(sig_e) == pytest.approx(peak_o, rel=0.10), (
        f"explicit peak {max(sig_e):.3f} != oracle peak {peak_o:.3f}")
    # final softened plateau: average the last 5% of steps to suppress the explicit dynamic ripple
    tail = sig_e[-_NSTEPS_EXPL // 20:]
    end_e = sum(tail) / len(tail)
    assert end_e == pytest.approx(end_o, rel=0.25), (
        f"explicit softened backbone {end_e:.3f} != oracle {end_o:.3f} at eps=0.03")


def test_explicit_completes_where_fixedstep_implicit_stalls():
    """The Tier-3 payoff, made concrete: the SAME single cube with a STEEP (snap-backy) softening law
    (lch=50 ⇒ eps_f=Gf/(ft·lch)≈6.7e-4, a sharp post-peak drop) under fixed-step implicit Newton
    (DisplacementControl, no step-cutting) STALLS at the immediate softening limit point (eps0=ft/E),
    converging only a handful of steps; the explicit run marches the full ramp. (Tier-1 needs the
    _drive_adaptive step-cutting + FullGeneral to get through — see the static softening test.) NB the
    contrast needs a genuinely snap-backy law: with the unit cube's default autoReg lch≈1 the softening
    is gradual (eps_f≈0.03) and plain DisplacementControl can path-follow many steps — platform-dependent
    and NOT the Tier-3 point. The steep law makes the implicit stall unambiguous on every platform."""
    lch = 50.0
    implicit_fixed = _run(lambda t: _mat(t, lch=lch), 0.03, 300)  # plain DisplacementControl, breaks on stall
    explicit = _run_explicit(CDL, 0.03, _NSTEPS_EXPL, _DT_EXPL, alphaM=20.0, lch=lch)
    assert len(implicit_fixed) < 60, (
        f"fixed-step implicit unexpectedly survived steep softening ({len(implicit_fixed)} steps) — "
        "contrast no longer demonstrates the Tier-3 payoff")
    assert len(explicit) == _NSTEPS_EXPL, (
        f"explicit failed to complete ({len(explicit)}/{_NSTEPS_EXPL}) — the Tier-3 claim")


# ===========================================================================
#  P5 — CONFINED-FIBER VIEW (DIM_BEAMFIBER, §4.6 "Mander by mechanism").
#  A single-fiber NDFiber section (area 1, at the centroid) in a zeroLengthSection element: the
#  imposed axial nodal displacement IS the section axial strain (unit length), and the section axial
#  force P == the fiber's nominal axial stress (A=1). getCopy("BeamFiber") routes the confined view
#  (lateral 11,22,12 condensed vs the passive hoop). Drives match the oracle confined_step path, so
#  P(eps) is byte-faithful to the numpy oracle; hoop=0 reduces to the free uniaxial-stress backbone.
# ===========================================================================
_CFIB_LCH = 1.0   # the material's default lch (no -lch / -autoRegularization) => gradual, convergent softening


def _build_confined_zls(hoopK, hoopFy=None):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 0.0, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.fix(2, 0, 1, 1, 1, 1, 1)              # free axial (dof 1) only => pure axial fiber strain
    kw = {"hoop": hoopK}
    if hoopFy is not None:
        kw["hoopFy"] = hoopFy
    _mat(1, **kw)
    ops.section("NDFiber", 1)
    ops.fiber(0.0, 0.0, 1.0, 1)               # single fiber, unit area => P == axial stress
    ops.element("zeroLengthSection", 1, 1, 2, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0, 0, 0, 0, 0, 0)
    ops.system("FullGeneral")                 # the confined tangent is non-symmetric
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-8, 200, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


def _Paxial():
    ops.eleResponse(1, "forces")
    return float(ops.eleResponse(1, "section", "force")[0])   # section axial resultant P (= sig_axial, A=1)


def _run_confined(hoopK, eps_target, nsteps, hoopFy=None):
    """DisplacementControl the axial dof to eps_target over nsteps; collect (eps, P). Stops on the first
    non-converged step (records what it reached) so a snap-back never fails the harness."""
    _build_confined_zls(hoopK, hoopFy)
    deps = eps_target / nsteps
    out = []
    for _ in range(nsteps):
        ops.integrator("DisplacementControl", 2, 1, deps)
        if ops.analyze(1) != 0:
            break
        out.append((ops.nodeDisp(2, 1), _Paxial()))
    return out


def test_confined_fiber_reduce_to_oracle():
    """hoop=0 BeamFiber == the free uniaxial-stress backbone: the section axial force P(eps) matches the
    numpy oracle drive_confined_fiber(hoop_K=0) along a compression path (the element-level reduce-to-1D,
    the §4.6 analog of the cube uniaxial gate). Pure axial => the fiber sees no shear, the lateral block
    is condensed exactly as the oracle."""
    eps_t = -1.8e-3
    n = 180
    rec = _run_confined(0.0, eps_t, n)
    assert len(rec) == n, f"hoop=0 confined fiber stalled at step {len(rec)}/{n}"
    mp = ref.make_material(_E, _NU, _FC, _FT)
    path = ref.np.array([e for e, _ in rec])
    o = ref.drive_confined_fiber(mp, path, _GF, _GC, lch=_CFIB_LCH, As=2.0, hoop_K=0.0)
    P = ref.np.array([p for _, p in rec])
    err = float(ref.np.max(ref.np.abs(P - o["sig11"])) / _FC)
    assert err < 5.0e-3, f"hoop=0 P(eps) vs oracle free backbone rel-fc err {err:.2e}"


def test_confinement_raises_capacity():
    """Mander-by-mechanism at the element level: a passive hoop (K>0) raises the peak axial capacity
    above the unconfined fiber. The non-associated dilatancy mobilizes the hoop self-consistently => the
    confined peak |P| exceeds the unconfined peak |P| (and the confining pressure is genuinely active)."""
    eps_t = -3.0e-3
    n = 300
    plain = _run_confined(0.0, eps_t, n)
    conf = _run_confined(1200.0, eps_t, n)
    assert len(plain) > 50 and len(conf) > 50, "confined-fiber runs stalled too early to compare peaks"
    pk_plain = max(abs(p) for _, p in plain)
    pk_conf = max(abs(p) for _, p in conf)
    assert pk_conf > 1.02 * pk_plain, (
        f"confined peak |P|={pk_conf:.3f} should exceed the unconfined peak |P|={pk_plain:.3f}")


def test_confined_fiber_view_builds():
    """getCopy(\"BeamFiber\") is reachable + the confined view loads through an NDFiber section without
    exit(-1) (order 3, type BeamFiber); a couple of elastic steps stay finite. Cheap reachability gate
    independent of the softening-convergence path."""
    rec = _run_confined(1200.0, -2.0e-4, 8, hoopFy=3.0)   # tiny strain (elastic-ish), hoop yield set
    assert len(rec) == 8
    for _, p in rec:
        assert math.isfinite(p)
