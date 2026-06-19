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
import pytest

from _testbed import ops
from _testbed.roundtrip import setresponse_smoke, database_roundtrip

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
    res = _run(lambda t: _mat(t), 0.05, 400)
    assert len(res) > 50
    sig = [s for _, s in res]
    peak = max(sig)
    assert peak == pytest.approx(_FT, rel=0.05), f"tensile peak {peak} != ft {_FT}"
    assert sig[-1] < 0.5 * peak, f"did not soften: end {sig[-1]} vs peak {peak}"
    # the peak is NOT at the very first step (it develops, then degrades)
    assert sig.index(peak) > 0
    # damage variable omega_t climbs toward 1 in the softening tail
    wt, wc = list(ops.eleResponse(1, "material", 1, "damage"))
    assert wt > 0.5, f"omega_t {wt} did not develop in tension softening"


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
    _build(lambda t: _mat(t))
    ops.integrator("DisplacementControl", 2, 1, 0.02 / 200)
    ok = 0
    for _ in range(200):
        if ops.analyze(1) != 0:
            break
        ok += 1
    assert ok > 50
    nominal = _sig_xx()
    seff = list(ops.eleResponse(1, "material", 1, "effectiveStress"))[0]
    wt, wc = list(ops.eleResponse(1, "material", 1, "damage"))
    assert wt > 0.3, f"expected tensile damage, omega_t={wt}"
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
