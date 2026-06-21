"""ExplicitBatheLNVD SMS family (Ladruno) — Zone-A validation.

Selective mass scaling on the Noh-Bathe + FLAC-local-non-viscous-damping ExplicitBatheLNVD
integrator, the LNVD siblings of the ExplicitBathe SMS family:
  * ExplicitBatheLNVDSMS            (classTag 33011) — LUMPED (ADR 36).
  * ExplicitBatheLNVDSMSConsistent  (classTag 33012) — CONSISTENT/Olovsson (ADR 38).

Command takes the Noh-Bathe p AND the FLAC alpha first:
  integrator ExplicitBatheLNVDSMS $p $alpha $dtTarget <opts>
  integrator ExplicitBatheLNVDSMSConsistent $p $alpha $dtTarget <opts>

Theory: Ladruno_implementation/{36,38}_*.md.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

EBL = "ExplicitBatheLNVD"
EBLSMS = "ExplicitBatheLNVDSMS"
EBLCONS = "ExplicitBatheLNVDSMSConsistent"
P = 0.54
ALPHA = 0.1               # FLAC local damping (exercises the LNVD path)

LENGTHS = [1.0] * 5 + [0.05] + [1.0] * 5
E, RHO, A = 1.0e4, 1.0, 1.0
DT_TARGET = 0.012
DT_REF = 0.002


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


# EBL-SMS-1 (lumped): supra-stable where plain LNVD diverges.
def test_ebl_sms_lumped_supra_stable():
    _, u_sms = _free_vibration((EBLSMS, P, ALPHA, DT_TARGET), DT_TARGET, 4.0)
    _, u_ebl = _free_vibration((EBL, P, ALPHA), DT_TARGET, 4.0)
    assert u_sms is not None and max(abs(u) for u in u_sms) < 1.0, "lumped LNVD-SMS not stable"
    assert u_ebl is None or max(abs(u) for u in u_ebl) > 1.0e2, (
        "plain ExplicitBatheLNVD was expected to DIVERGE at dt=%g" % DT_TARGET
    )


# EBL-SMS-2 (lumped): injects nodal mass.
def test_ebl_sms_lumped_injects_nodal_mass():
    _build_bar()
    assert abs(ops.nodeMass(2)[0]) < 1e-9
    ops.integrator(EBLSMS, P, ALPHA, DT_TARGET)
    ops.analysis("Transient")
    ops.analyze(1, DT_TARGET)
    assert ops.nodeMass(2)[0] > 0.05, "lumped LNVD-SMS should inject nodal mass"


# EBL-CONS-1 (consistent): supra-stable + NO nodal mass mutation.
def test_ebl_consistent_supra_stable_no_nodal_mass():
    _, u_cons = _free_vibration((EBLCONS, P, ALPHA, DT_TARGET), DT_TARGET, 4.0)
    assert u_cons is not None and max(abs(u) for u in u_cons) < 1.0, "LNVD-consistent not stable"
    _build_bar()
    ops.integrator(EBLCONS, P, ALPHA, DT_TARGET)
    ops.analysis("Transient")
    ops.analyze(1, DT_TARGET)
    assert abs(ops.nodeMass(2)[0]) < 1e-9, "consistent LNVD-SMS must NOT mutate nodal mass"


# EBL-CONS-2 (selling point): consistent preserves f1; lumped shifts it. alpha=0 for a
# clean FFT peak (FLAC damping would smear it; both legs would still share it, but 0 is
# cleanest for the frequency read).
def test_ebl_consistent_preserves_frequency_vs_lumped():
    pytest.importorskip("numpy")
    a0 = 0.0
    T_END = 6.0
    dt_r, u_ref = _free_vibration((EBL, P, a0), DT_REF, T_END)
    dt_c, u_cons = _free_vibration((EBLCONS, P, a0, DT_TARGET), DT_TARGET, T_END)
    dt_l, u_lump = _free_vibration((EBLSMS, P, a0, DT_TARGET), DT_TARGET, T_END)
    assert u_ref and u_cons and u_lump
    f_ref = _dominant_freq(dt_r, u_ref)
    e_cons = abs(_dominant_freq(dt_c, u_cons) - f_ref) / f_ref
    e_lump = abs(_dominant_freq(dt_l, u_lump) - f_ref) / f_ref
    assert e_lump > 0.15, "control: lumped LNVD-SMS f1 error %.3f expected large" % e_lump
    assert e_cons < 0.05, "consistent LNVD-SMS f1 error %.3f too large" % e_cons
    assert e_cons < 0.3 * e_lump, "consistent (%.3f) should beat lumped (%.3f)" % (e_cons, e_lump)


# EBL-CONS-3: reduce-to-base bit-identical to ExplicitBatheLNVD (same alpha).
def test_ebl_consistent_reduces_to_base():
    dt = 5.0e-5
    _, u_cons = _free_vibration((EBLCONS, P, ALPHA, 1.0e-6), dt, 0.02)
    _, u_ebl = _free_vibration((EBL, P, ALPHA), dt, 0.02)
    assert u_cons and u_ebl and len(u_cons) == len(u_ebl)
    scale = max(abs(b) for b in u_ebl) or 1.0
    worst = max(abs(a - b) for a, b in zip(u_cons, u_ebl))
    assert worst / scale < 1.0e-10, "reduce-to-base not bit-identical: rel %.3e" % (worst / scale)


# EBL-CONS-4: PCG converges fast at both sub-steps.
def test_ebl_consistent_pcg_converges_fast(capfd):
    import re
    _free_vibration((EBLCONS, P, ALPHA, DT_TARGET, "-verbose"), DT_TARGET, 0.2)
    out = capfd.readouterr()
    iters = [int(m) for m in re.findall(r"PCG\s+(\d+)\s+iters", out.err + out.out)]
    assert iters, "no PCG iteration report under -verbose"
    assert max(iters) <= 30, "PCG slow to converge: %r" % iters
