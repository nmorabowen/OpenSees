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
    ops.nDMaterial(*args)


def _build(mat_fn):
    """Single stdBrick unit cube, 1/8-symmetry determinate restraints. UNSYMMETRIC solver."""
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
    ops.system("FullGeneral")              # the CDPM2 tangent is NON-SYMMETRIC
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


def _drive_adaptive(mat_fn, eps_target, base_steps, max_cuts=7):
    """Displacement-control driver with step-CUTTING through the softening limit point — the only
    way a single implicit element gets past an unconfined tension/compression peak (the snap-back
    regime). Returns [(eps_xx, sig_xx, omega_t)] for every converged increment."""
    _build(mat_fn)
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
