"""CentralDifferenceSMSConsistent (Ladruno CONSISTENT / Olovsson mass scaling,
classTag 33008) — Zone-A validation. ADR 38 (the v2 of CentralDifferenceSMS 33007).

Conventional (lumped) mass scaling adds an isotropic diagonal increment to the
throttling elements' nodes, which shifts ALL global frequencies — including the
fundamental f1 that seismic SSI cannot afford to move. The consistent path injects the
centroidal Olovsson scaling mass  M_bar_e = beta[diag(m_a) - m m^T/M_e]  whose row sums
are zero, so rigid-body translation gets no added inertia and f1 is PRESERVED while the
relative/deformation modes that govern the explicit step are loaded. M_tilde is
non-diagonal, solved matrix-free by lumped-preconditioned CG inside the step.

Each test has a control that FAILS on a broken / lumped-behaving / no-op integrator:
  TCONS-1  constructs + runs (parse + lifecycle).
  TCONS-2  supra-stable: plain CD diverges at dtTarget, consistent SMS holds (necessity).
  TCONS-3  injects NO nodal mass (distinguishes it from the lumped sibling, which does).
  TCONS-4  f1 PRESERVATION: consistent tracks the fine-dt reference frequency far better
           than lumped SMS at the SAME aggressive dtTarget (the v2 selling point); the
           lumped leg is the non-vacuity control (it MUST shift f1 a lot).
  TCONS-5  reduce-to-base: dtTarget below min dt_e => nothing scaled => trajectory
           bit-identical to CentralDifferenceLadruno.
  TCONS-6  PCG health: -verbose reports a small iteration count (capfd).

Oracle (numpy): Ladruno_implementation/mass_scaling_consistent/oracle_olovsson_sms.py.
Theory: Ladruno_implementation/38_ladruno_consistent_mass_scaling_adr.md.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

CONS = "CentralDifferenceSMSConsistent"
SMS = "CentralDifferenceSMS"
CDL = "CentralDifferenceLadruno"


# --------------------------------------------------------------------------
# A fixed–free axial bar: 5 bulk elements (L=1), ONE tiny interior element
# (L=0.05, ~20x smaller -> ~400x smaller dt_e), then 5 more bulk elements. This is
# the oracle's Case A. dtTarget=0.012 is aggressive (above even the bulk dt_e ~0.01),
# so EVERY element is scaled: lumped scaling then loads the whole bar and craters f1,
# while Olovsson loads only the relative modes and keeps f1.
# --------------------------------------------------------------------------
LENGTHS = [1.0] * 5 + [0.05] + [1.0] * 5
E, RHO, A = 1.0e4, 1.0, 1.0
DT_TARGET = 0.012
DT_REF = 0.002          # < global stable step (~0.0085) -> stable plain-CD reference


def _build_bar():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    x = 0.0
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 1, 1)
    coords = {1: 0.0}
    for i, L in enumerate(LENGTHS):
        x += L
        nd = i + 2
        ops.node(nd, x, 0.0)
        ops.fix(nd, 0, 1)          # axial only (y fixed)
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
    """Seed the bar with the static-like LINEAR axial field u(x)=ic_slope*x (mode-1
    dominant) and free-vibrate. Returns (dt, [u_tip(t)]) or (dt, None) if it diverged."""
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
    """Dominant (peak-power) frequency of the tip time history via real FFT."""
    import numpy as np
    x = np.asarray(us, float)
    x = x - x.mean()
    n = x.size
    spec = np.abs(np.fft.rfft(x * np.hanning(n)))
    freqs = np.fft.rfftfreq(n, d=dt)
    spec[0] = 0.0
    return float(freqs[int(np.argmax(spec))])


# --------------------------------------------------------------------------
# TCONS-1: parse + lifecycle smoke.
# --------------------------------------------------------------------------
def test_consistent_constructs_and_runs():
    dt, us = _free_vibration((CONS, DT_TARGET, "-verbose"), DT_TARGET, 0.2)
    assert us is not None and len(us) > 5, "CentralDifferenceSMSConsistent failed to run"


# --------------------------------------------------------------------------
# TCONS-2 (necessity): plain CD DIVERGES at dtTarget (the tiny element's dt_cr is
# ~400x smaller); consistent SMS stays bounded.
# --------------------------------------------------------------------------
def test_consistent_supra_stable():
    _, u_cons = _free_vibration((CONS, DT_TARGET), DT_TARGET, 4.0)
    _, u_cdl = _free_vibration((CDL,), DT_TARGET, 4.0)
    assert u_cons is not None and max(abs(u) for u in u_cons) < 1.0, (
        "consistent SMS not stable/bounded at supra-stable dt"
    )
    assert u_cdl is None or max(abs(u) for u in u_cdl) > 1.0e2, (
        "plain CentralDifferenceLadruno was expected to DIVERGE at dt=%g" % DT_TARGET
    )


# --------------------------------------------------------------------------
# TCONS-3 (distinguishing): the consistent path injects NO nodal mass — the scaling
# mass lives only in the integrator's blocks (matrix-free). The lumped sibling, by
# contrast, mutates nodeMass. This fails if the integrator accidentally took the
# lumped path.
# --------------------------------------------------------------------------
def test_consistent_injects_no_nodal_mass():
    _build_bar()
    base = ops.nodeMass(2)
    assert abs(base[0]) < 1e-9, "expected zero baseline nodeMass, got %r" % base

    ops.integrator(CONS, DT_TARGET)
    ops.analysis("Transient")
    ops.analyze(1, DT_TARGET)
    after = ops.nodeMass(2)
    assert abs(after[0]) < 1e-9, (
        "consistent SMS must NOT mutate nodal mass (it is matrix-free); got %r — "
        "this looks like the lumped path" % after
    )

    # control: the LUMPED sibling DOES inject nodal mass on the same model.
    _build_bar()
    ops.integrator(SMS, DT_TARGET)
    ops.analysis("Transient")
    ops.analyze(1, DT_TARGET)
    lumped = ops.nodeMass(2)
    assert lumped[0] > 0.05, (
        "control: lumped CentralDifferenceSMS should inject nodal mass; got %r" % lumped
    )


# --------------------------------------------------------------------------
# TCONS-4 (the v2 selling point): consistent SMS preserves f1; lumped SMS shifts it.
# --------------------------------------------------------------------------
def test_consistent_preserves_frequency_vs_lumped():
    pytest.importorskip("numpy")
    T_END = 6.0
    dt_r, u_ref = _free_vibration((CDL,), DT_REF, T_END)
    dt_c, u_cons = _free_vibration((CONS, DT_TARGET), DT_TARGET, T_END)
    dt_l, u_lump = _free_vibration((SMS, DT_TARGET), DT_TARGET, T_END)
    assert u_ref and u_cons and u_lump, "a leg failed to run"

    f_ref = _dominant_freq(dt_r, u_ref)
    f_cons = _dominant_freq(dt_c, u_cons)
    f_lump = _dominant_freq(dt_l, u_lump)
    assert f_ref > 0.0

    e_cons = abs(f_cons - f_ref) / f_ref
    e_lump = abs(f_lump - f_ref) / f_ref

    # non-vacuity: lumped scaling MUST visibly shift f1 (added mass everywhere).
    assert e_lump > 0.15, (
        "control failed: lumped SMS f1 error %.3f expected large (>0.15); "
        "f_ref=%.4f f_lump=%.4f" % (e_lump, f_ref, f_lump)
    )
    # the claim: consistent preserves f1 and is markedly better than lumped.
    assert e_cons < 0.05, (
        "consistent SMS f1 error %.3f too large (>0.05); f_ref=%.4f f_cons=%.4f"
        % (e_cons, f_ref, f_cons)
    )
    assert e_cons < 0.3 * e_lump, (
        "consistent (%.3f) should be >3x better than lumped (%.3f)" % (e_cons, e_lump)
    )


# --------------------------------------------------------------------------
# TCONS-5: dtTarget below min dt_e => nothing scaled => the refineAccel hook is a
# no-op => trajectory bit-identical to CentralDifferenceLadruno.
# --------------------------------------------------------------------------
def test_consistent_reduces_to_base():
    dt = 0.0009                      # < the tiny element's dt_e (~5e-4? no: choose tiny)
    dt = 5.0e-5                      # safely below every element's dt_e -> no scaling
    _, u_cons = _free_vibration((CONS, 1.0e-6), dt, 0.02)   # target unreachably small
    _, u_cdl = _free_vibration((CDL,), dt, 0.02)
    assert u_cons and u_cdl and len(u_cons) == len(u_cdl)
    worst = max(abs(a - b) for a, b in zip(u_cons, u_cdl))
    scale = max(abs(b) for b in u_cdl) or 1.0
    assert worst / scale < 1.0e-10, (
        "reduce-to-base not bit-identical to CentralDifferenceLadruno: worst rel %.3e"
        % (worst / scale)
    )


# --------------------------------------------------------------------------
# TCONS-6: the matrix-free PCG converges in a few iterations (lumped diagonal is an
# excellent preconditioner). -verbose prints "PCG <k> iters"; assert k is small.
# --------------------------------------------------------------------------
def test_consistent_pcg_converges_fast(capfd):
    import re
    _free_vibration((CONS, DT_TARGET, "-verbose"), DT_TARGET, 0.2)
    out = capfd.readouterr()
    text = out.err + out.out
    # the per-step report is exactly:  "... PCG <k> iters, rel.resid <r>"
    iters = [int(m) for m in re.findall(r"PCG\s+(\d+)\s+iters", text)]
    assert iters, "no PCG iteration report found under -verbose"
    assert max(iters) <= 30, "PCG slow to converge: iters=%r" % iters


# --------------------------------------------------------------------------
# TCONS-7 (adversarial-review regression): under the DEFAULT Transformation constraint
# handler, an element touching an MP-constrained node has a TRANSFORMED FE id whose size
# can EXCEED the element's own n (a retained node's TransformationDOF_Group absorbs the
# slave's DOFs). Pairing that oversized id with the n×n M_bar was an out-of-bounds read.
# The consistent path now excludes elements touching EITHER a slave or a master MP node
# (+ a defensive getID().Size()==n guard). This builds exactly that hazard — a scaled
# tiny element whose node is the equalDOF master — and asserts: no crash/garbage, the run
# stays bounded, and the excluded element is reported as still-governing (nConstrained>0).
# --------------------------------------------------------------------------
def test_consistent_excludes_constraint_touching_elements(capfd):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
    ops.node(2, 0.1, 0.0); ops.fix(2, 0, 1)   # tiny stiff element 1-2 (sub-target)
    ops.node(3, 1.1, 0.0); ops.fix(3, 0, 1)
    ops.node(4, 1.1, 0.5); ops.fix(4, 0, 1)   # slave node, tied to node 3 in x
    ops.uniaxialMaterial("Elastic", 1, 1.0e6)  # stiff/short -> tiny dt_e
    ops.uniaxialMaterial("Elastic", 2, 1.0e4)
    ops.element("Truss", 1, 1, 2, 1.0, 1, "-rho", 1.0)   # tiny, touches node 2
    ops.element("Truss", 2, 2, 3, 1.0, 2, "-rho", 1.0)   # bulk, touches master node 3
    ops.equalDOF(3, 4, 1)                       # node 3 = retained/master, node 4 = slave
    ops.constraints("Transformation")           # the DEFAULT handler that expands the FE id
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.setNodeDisp(2, 1, 1.0e-3, "-commit")
    ops.integrator(CONS, 5.0e-3, "-verbose")    # dtTarget scales the tiny element
    ops.analysis("Transient")

    umax = 0.0
    for _ in range(300):
        if ops.analyze(1, 5.0e-3) != 0:
            pytest.fail("consistent SMS step failed (possible OOB/garbage from the "
                        "transformed FE id) on a constrained model")
        u = abs(ops.nodeDisp(2, 1))
        if not math.isfinite(u):
            pytest.fail("non-finite response — out-of-bounds M_bar read corrupted the solve")
        umax = max(umax, u)
    assert umax < 1.0, "run not bounded: max|u|=%r" % umax

    out = capfd.readouterr()
    text = out.err + out.out
    assert "EXCLUDED" in text or "MP-constrained" in text, (
        "expected a constraint-exclusion report (the constraint-touching element should be "
        "excluded from scaling, not silently OOB-scaled)"
    )
