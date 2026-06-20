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
    if "lch" in kw:
        args += ["-lch", kw["lch"]]
    if kw.get("implex"):
        args += ["-implex"]
    if "eta" in kw:
        args += ["-eta", kw["eta"]]
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
