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


# ------------------------------------------- 6. reversed slip: friction slider
def _slider_oracle(P, strains):
    """Oracle for the whole history: backbone past the peak |slip|; inside it, k0
    from the committed point with |tau| capped at the envelope of the peak."""
    k0 = _k0(P)
    slip, tau, peak = 0.0, 0.0, 0.0
    out = []
    for e in strains:
        a = abs(e)
        if a >= peak:
            tau, peak = math.copysign(_env(a, P), e), a
        else:
            cap = _env(peak, P)
            tau = max(-cap, min(cap, tau + k0 * (e - slip)))
        slip = e
        out.append(tau)
    return out


def _walk(targets, ds):
    """Slip history through `targets` in steps no larger than ds."""
    path, here = [], 0.0
    for target in targets:
        n = max(1, math.ceil(abs(target - here) / ds))
        path += [here + (target - here) * (i + 1) / n for i in range(n)]
        here = target
    return path


# past s3, then fully reversed twice; inner reversals sit on the friction branch
CYCLE = [1.5 * S3, 0.0, -0.5 * S3, 0.25 * S3, -1.5 * S3, 0.0, 2.0 * S3]


def test_reversal_follows_slider_oracle():
    P = _params()
    path = _walk(CYCLE, 0.05e-3)
    got = _drive(P, path)
    for e, want, (sig, _t) in zip(path, _slider_oracle(P, path), got):
        assert sig == pytest.approx(want, rel=1e-9, abs=1e-9), f"slip {e}"


def test_reversal_is_bounded_and_has_no_jump():
    """The defect: an unbounded unload line (about 28*tau_max), then a sign jump
    the step after the slip crossed zero. Every step moves at most k0*ds."""
    P = _params()
    path = _walk(CYCLE, 0.05e-3)
    stress = [sig for sig, _t in _drive(P, path)]
    assert max(abs(t) for t in stress) <= TAU_MAX * (1.0 + 1e-12)
    prev_e, prev_t = 0.0, 0.0
    for e, t in zip(path, stress):
        assert abs(t - prev_t) <= _k0(P) * abs(e - prev_e) * (1.0 + 1e-9), f"slip {e}"
        prev_e, prev_t = e, t


@pytest.mark.parametrize("peak, cap", [(0.5 * (S1 + S2), TAU_MAX),
                                       (2.0 * S3, TAU_F)])
def test_reversed_friction_stress_is_peak_envelope(peak, cap):
    """Reversing from the peak slides at -envelope(peak), through zero slip and on
    to just inside -peak, with the finite k0 tangent."""
    P = _params()
    got = _drive(P, [peak, 0.0, -0.5 * peak, -0.99 * peak])
    for sig, tan in got[1:]:
        assert sig == pytest.approx(-cap, rel=1e-9)
        assert tan == pytest.approx(_k0(P), rel=1e-9)


def test_reversal_from_friction_unloads_at_k0():
    """Turning back while sliding leaves the cap at once, at k0."""
    P = _params()
    peak, ds = 2.0 * S3, 0.01e-3
    got = _drive(P, [peak, 0.0, ds])
    assert got[2][0] == pytest.approx(-TAU_F + _k0(P) * ds, rel=1e-9)


def test_reversal_rejoins_backbone_past_opposite_peak():
    P = _params()
    peak = 0.5 * (S2 + S3)                # on the softening branch
    beyond = 1.1 * peak
    got = _drive(P, [peak, 0.0, -peak, -beyond])
    assert got[2][0] == pytest.approx(-_env(peak, P), rel=1e-9)
    assert got[3][0] == pytest.approx(-_env(beyond, P), rel=1e-9)
    assert got[3][1] < 0.0                # softening tangent again


def test_closed_slip_loop_dissipates():
    """A full residual-friction loop dissipates about 4*tau_f*peak."""
    P = _params()
    peak = 2.0 * S3
    path = _walk([peak, -peak, peak], 0.05e-3)
    stress = [sig for sig, _t in _drive(P, path)]
    start = path.index(peak)
    work = sum(0.5 * (stress[i] + stress[i + 1]) * (path[i + 1] - path[i])
               for i in range(start, len(path) - 1))
    assert work == pytest.approx(4.0 * TAU_F * peak, rel=0.01)
