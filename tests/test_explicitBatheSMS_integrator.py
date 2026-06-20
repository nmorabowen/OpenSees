"""ExplicitBathe SMS family (Ladruno) — Zone-A validation.

Selective mass scaling on the Noh-Bathe ExplicitBathe integrator, the explicit-Bathe
siblings of CentralDifferenceSMS / CentralDifferenceSMSConsistent:
  * ExplicitBatheSMS            (classTag 33009) — LUMPED (additive nodal mass; ADR 36).
  * ExplicitBatheSMSConsistent  (classTag 33010) — CONSISTENT/Olovsson (matrix-free PCG
                                 at BOTH Noh-Bathe sub-steps; ADR 38).

Command takes the Noh-Bathe sub-step parameter p first:
  integrator ExplicitBatheSMS $p $dtTarget <opts>
  integrator ExplicitBatheSMSConsistent $p $dtTarget <opts>

Each test carries a control that FAILS on a broken / no-op / lumped-behaving integrator.
Theory: Ladruno_implementation/{36,38}_*.md; numpy oracle:
Ladruno_implementation/mass_scaling_consistent/oracle_olovsson_sms.py.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

EB = "ExplicitBathe"
EBSMS = "ExplicitBatheSMS"
EBCONS = "ExplicitBatheSMSConsistent"
P = 0.54                  # Noh-Bathe sub-step parameter (good default)

LENGTHS = [1.0] * 5 + [0.05] + [1.0] * 5   # oracle Case-A bar (one tiny interior element)
E, RHO, A = 1.0e4, 1.0, 1.0
DT_TARGET = 0.012
DT_REF = 0.002            # < global stable step -> stable plain-EB reference


def _build_bar():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    x = 0.0
    ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
    coords = {1: 0.0}
    for i, L in enumerate(LENGTHS):
        x += L
        nd = i + 2
        ops.node(nd, x, 0.0); ops.fix(nd, 0, 1)
        coords[nd] = x
    ops.uniaxialMaterial("Elastic", 1, E)
    for i in range(len(LENGTHS)):
        ops.element("Truss", i + 1, i + 1, i + 2, A, 1, "-rho", RHO)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    return coords


def _free_vibration(integrator_args, dt, t_end, ic_slope=1.0e-3):
    coords = _build_bar()
    tip = max(coords)
    for nd, x in coords.items():
        if x > 0.0:
            ops.setNodeDisp(nd, 1, ic_slope * x, "-commit")
    ops.integrator(*integrator_args)
    ops.analysis("Transient")
    n = int(round(t_end / dt))
    us = [ops.nodeDisp(tip, 1)]
    for _ in range(n):
        if ops.analyze(1, dt) != 0:
            return dt, None
        u = ops.nodeDisp(tip, 1)
        if not math.isfinite(u):
            return dt, None
        us.append(u)
    return dt, us


def _dominant_freq(dt, us):
    import numpy as np
    x = np.asarray(us, float)
    x = x - x.mean()
    n = x.size
    spec = np.abs(np.fft.rfft(x * np.hanning(n)))
    freqs = np.fft.rfftfreq(n, d=dt)
    spec[0] = 0.0
    return float(freqs[int(np.argmax(spec))])


# --------------------------------------------------------------------------
# EBSMS-1 (lumped): constructs + runs, and is supra-stable where plain EB diverges.
# --------------------------------------------------------------------------
def test_eb_sms_lumped_supra_stable():
    _, u_sms = _free_vibration((EBSMS, P, DT_TARGET), DT_TARGET, 4.0)
    _, u_eb = _free_vibration((EB, P), DT_TARGET, 4.0)
    assert u_sms is not None and max(abs(u) for u in u_sms) < 1.0, "lumped EB-SMS not stable"
    assert u_eb is None or max(abs(u) for u in u_eb) > 1.0e2, (
        "plain ExplicitBathe was expected to DIVERGE at dt=%g" % DT_TARGET
    )


# --------------------------------------------------------------------------
# EBSMS-2 (lumped): injects nodal mass (the lumped mechanism; RHS assembly sees it).
# --------------------------------------------------------------------------
def test_eb_sms_lumped_injects_nodal_mass():
    _build_bar()
    base = ops.nodeMass(2)
    assert abs(base[0]) < 1e-9
    ops.integrator(EBSMS, P, DT_TARGET)
    ops.analysis("Transient")
    ops.analyze(1, DT_TARGET)
    after = ops.nodeMass(2)
    assert after[0] > 0.05, "lumped EB-SMS should inject nodal mass; got %r" % after


# --------------------------------------------------------------------------
# EBCONS-1 (consistent): supra-stable + injects NO nodal mass (matrix-free).
# --------------------------------------------------------------------------
def test_eb_consistent_supra_stable_no_nodal_mass():
    _, u_cons = _free_vibration((EBCONS, P, DT_TARGET), DT_TARGET, 4.0)
    assert u_cons is not None and max(abs(u) for u in u_cons) < 1.0, "EB-consistent not stable"

    _build_bar()
    ops.integrator(EBCONS, P, DT_TARGET)
    ops.analysis("Transient")
    ops.analyze(1, DT_TARGET)
    after = ops.nodeMass(2)
    assert abs(after[0]) < 1e-9, (
        "consistent EB-SMS must NOT mutate nodal mass (matrix-free); got %r" % after
    )


# --------------------------------------------------------------------------
# EBCONS-2 (the selling point): consistent preserves f1; lumped shifts it.
# --------------------------------------------------------------------------
def test_eb_consistent_preserves_frequency_vs_lumped():
    pytest.importorskip("numpy")
    T_END = 6.0
    dt_r, u_ref = _free_vibration((EB, P), DT_REF, T_END)
    dt_c, u_cons = _free_vibration((EBCONS, P, DT_TARGET), DT_TARGET, T_END)
    dt_l, u_lump = _free_vibration((EBSMS, P, DT_TARGET), DT_TARGET, T_END)
    assert u_ref and u_cons and u_lump, "a leg failed to run"

    f_ref = _dominant_freq(dt_r, u_ref)
    f_cons = _dominant_freq(dt_c, u_cons)
    f_lump = _dominant_freq(dt_l, u_lump)
    assert f_ref > 0.0
    e_cons = abs(f_cons - f_ref) / f_ref
    e_lump = abs(f_lump - f_ref) / f_ref

    assert e_lump > 0.15, "control: lumped EB-SMS f1 error %.3f expected large" % e_lump
    assert e_cons < 0.05, "consistent EB-SMS f1 error %.3f too large" % e_cons
    assert e_cons < 0.3 * e_lump, "consistent (%.3f) should beat lumped (%.3f)" % (e_cons, e_lump)


# --------------------------------------------------------------------------
# EBCONS-3: reduce-to-base — dtTarget below min dt_e => bit-identical to ExplicitBathe.
# --------------------------------------------------------------------------
def test_eb_consistent_reduces_to_base():
    dt = 5.0e-5
    _, u_cons = _free_vibration((EBCONS, P, 1.0e-6), dt, 0.02)
    _, u_eb = _free_vibration((EB, P), dt, 0.02)
    assert u_cons and u_eb and len(u_cons) == len(u_eb)
    scale = max(abs(b) for b in u_eb) or 1.0
    worst = max(abs(a - b) for a, b in zip(u_cons, u_eb))
    assert worst / scale < 1.0e-10, "reduce-to-base not bit-identical: rel %.3e" % (worst / scale)


# --------------------------------------------------------------------------
# EBCONS-4: the matrix-free PCG converges fast at BOTH Noh-Bathe sub-steps.
# --------------------------------------------------------------------------
def test_eb_consistent_pcg_converges_fast(capfd):
    import re
    _free_vibration((EBCONS, P, DT_TARGET, "-verbose"), DT_TARGET, 0.2)
    out = capfd.readouterr()
    iters = [int(m) for m in re.findall(r"PCG\s+(\d+)\s+iters", out.err + out.out)]
    assert iters, "no PCG iteration report under -verbose"
    assert max(iters) <= 30, "PCG slow to converge: %r" % iters
