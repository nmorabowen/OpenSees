"""LadrunoBondSlip (1D bond-slip tau-s) — Zone-A material battery.

Verifies the CEB-FIP MC2010 backbone + the ADR-20 D4 must-fix regularizations
against a pure-Python oracle, via the direct uniaxial harness
(testUniaxialMaterial / setStrain / getStress / getTangent; setStrain COMMITS each
step). See SRC/material/uniaxial/LadrunoBondSlip.{h,cpp} and
Ladruno_implementation/20_ladruno_embedded_reinforcement_adr.md (D4).
"""
import math
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# a representative bond law (consistent units)
TAU_MAX, S1, S2, S3, TAU_F, ALPHA = 10.0, 1.0e-3, 3.0e-3, 10.0e-3, 2.0, 0.4


def _params(Gf=0.0, s0=0.0):
    return dict(tau_max=TAU_MAX, s1=S1, s2=S2, s3=S3, tau_f=TAU_F,
                alpha=ALPHA, Gf=Gf, s0=s0)


def _k0(P):
    s0 = P["s0"] if P["s0"] > 0 else 0.1 * P["s1"]
    return P["tau_max"] * (s0 / P["s1"]) ** P["alpha"] / s0


def _s3eff(P):
    if P["Gf"] > 0 and P["tau_max"] > P["tau_f"]:
        return P["s2"] + 2.0 * P["Gf"] / (P["tau_max"] - P["tau_f"])
    return P["s3"]


def _env(a, P):
    """oracle backbone for a = |slip| >= 0 (mirrors the C++ envelope())."""
    s0 = P["s0"] if P["s0"] > 0 else 0.1 * P["s1"]
    k0, s3 = _k0(P), _s3eff(P)
    tm, s1, s2, tf, al = P["tau_max"], P["s1"], P["s2"], P["tau_f"], P["alpha"]
    if a <= s0:
        return k0 * a
    if a <= s1:
        return tm * (a / s1) ** al
    if a <= s2:
        return tm
    if a <= s3:
        return tm - (tm - tf) * (a - s2) / (s3 - s2)
    return tf


def _build(tag, P):
    ops.uniaxialMaterial("LadrunoBondSlip", tag, P["tau_max"], P["s1"], P["s2"],
                         P["s3"], P["tau_f"], P["alpha"], "-Gf", P["Gf"],
                         "-s0", P["s0"])


def _drive(P, strains):
    ops.wipe()
    _build(1, P)
    ops.testUniaxialMaterial(1)
    out = []
    for e in strains:
        ops.setStrain(float(e))
        out.append((ops.getStress(), ops.getTangent()))
    return out


# ---------------------------------------------------------------- 1. backbone
def test_backbone_matches_oracle():
    """Monotonic loading walks the full backbone (linear / power / plateau /
    softening / residual); each point must match the oracle."""
    P = _params()
    s3 = _s3eff(P)
    strains = [0.5 * 0.1 * S1,            # linear initial segment
               math.sqrt(0.1 * S1 * S1),  # ascending power law
               0.5 * (S1 + S2),           # plateau
               0.5 * (S2 + s3),           # softening
               2.0 * s3]                  # residual
    got = _drive(P, strains)
    for e, (sig, _t) in zip(strains, got):
        assert sig == pytest.approx(_env(e, P), rel=1e-6, abs=1e-9)


def test_plateau_equals_tau_max():
    P = _params()
    (sig, _t), = _drive(P, [0.5 * (S1 + S2)])
    assert sig == pytest.approx(TAU_MAX, rel=1e-9)


def test_residual_equals_tau_f():
    P = _params()
    (sig, _t), = _drive(P, [3.0 * S3])
    assert sig == pytest.approx(TAU_F, rel=1e-9)


# -------------------------------------------------- 2. D4.1 initial-slip reg.
def test_initial_tangent_finite_equals_k0():
    """The power law has dtau/ds -> inf at s->0; the linear s<s0 segment fixes it.
    A first tiny step must give the finite k0, not a blow-up."""
    P = _params()
    (sig, tan), = _drive(P, [0.25 * 0.1 * S1])   # inside [0, s0]
    assert math.isfinite(tan)
    assert tan == pytest.approx(_k0(P), rel=1e-6)
    assert sig == pytest.approx(_k0(P) * 0.25 * 0.1 * S1, rel=1e-6)


# ----------------------------------------------------------- 3. sign symmetry
@pytest.mark.parametrize("a", [0.5e-3, 2.0e-3, 6.0e-3])
def test_sign_symmetry(a):
    P = _params()
    (sp, _tp), = _drive(P, [a])
    (sn, _tn), = _drive(P, [-a])
    assert sn == pytest.approx(-sp, rel=1e-9, abs=1e-12)


# ------------------------------------------------- 4. D4.3 Gf overrides s3
def test_Gf_overrides_s3():
    """With Gf>0 the softening branch dissipates Gf per unit area, so it reaches
    the residual exactly at s3 = s2 + 2 Gf/(tau_max - tau_f)."""
    Gf = 0.5 * (TAU_MAX - TAU_F) * (6.0e-3)      # -> s3eff = s2 + 6e-3
    P = _params(Gf=Gf)
    s3 = _s3eff(P)
    assert s3 == pytest.approx(S2 + 6.0e-3, rel=1e-9)
    (below, _), (above, _) = _drive(P, [0.99 * s3, 1.01 * s3])
    assert below > TAU_F + 1e-6          # still softening just below s3
    assert above == pytest.approx(TAU_F, rel=1e-9)   # residual just above


# ----------------------------------------------- 5. elastic unload at k0
def test_unload_is_elastic_k0():
    """Load into the plateau, then unload one step: stress drops by k0*ds and the
    tangent is k0 (finite, no negative-stiffness surprise on unload)."""
    P = _params()
    e_peak = 0.5 * (S1 + S2)             # plateau
    ds = 0.2e-3
    res = _drive(P, [e_peak, e_peak - ds])
    (sig0, _t0), (sig1, tan1) = res
    assert sig0 == pytest.approx(TAU_MAX, rel=1e-9)
    assert tan1 == pytest.approx(_k0(P), rel=1e-6)
    assert (sig0 - sig1) == pytest.approx(_k0(P) * ds, rel=1e-5)
