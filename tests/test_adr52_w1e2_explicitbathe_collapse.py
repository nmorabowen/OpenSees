"""ADR-52 W1-E2 — ExplicitBathe* 6 -> 1 collapse: Zone-A byte-identity + cross-validation.

The six historical explicit-Bathe class tags
  ExplicitBathe (33000) / ExplicitBatheLNVD (33002) / ExplicitBatheSMS (33009) /
  ExplicitBatheSMSConsistent (33010) / ExplicitBatheLNVDSMS (33011) /
  ExplicitBatheLNVDSMSConsistent (33012)
are now ONE class selected by orthogonal flags on the integrator command:
  integrator ExplicitBathe $p [-lnvd <alpha>] [-sms <dtTarget>] [-consistent] ...

The per-alias batteries (test_explicitBathe*.py) still run under the DEPRECATED command
names and pin the behaviour. This file pins the *collapse contract* the aliases cannot:

  1. each new FLAG form is BYTE-IDENTICAL to its legacy alias (same trace, rel < 1e-12)
     across all five non-trivial flag combos — the proof the merge changed no numerics;
  2. the LNVD x CONSISTENT cross-layer (the combo never exercised before W1-E2): FLAC
     local damping + the matrix-free Olovsson PCG compose — stable, no nodal-mass
     mutation, and the damping still relaxes a loaded bar to its static solution;
  3. flag GATING reduces to the base path: -lnvd 0.0 (alpha=0) == plain ExplicitBathe.

Model: the oracle Case-A bar from test_explicitBatheSMS_integrator.py (one tiny interior
element so dtTarget=0.012 throttles -> SMS has work to do); system Diagonal + algorithm
Linear (the explicit recipe). See Ladruno_implementation/52_*adr.md (W1-E2).
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

P = 0.54                  # Noh-Bathe sub-step parameter
ALPHA = 0.8               # FLAC local-damping coefficient (classic default)
DT_TARGET = 0.012         # SMS target step (throttled by the tiny interior element)

LENGTHS = [1.0] * 5 + [0.05] + [1.0] * 5   # oracle Case-A bar
E, RHO, A = 1.0e4, 1.0, 1.0


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
    """Seed a linear IC, run the explicit transient, return the tip-disp trace."""
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
            return None
        u = ops.nodeDisp(tip, 1)
        if not math.isfinite(u):
            return None
        us.append(u)
    return us


def _relax_to_static(integrator_args, P_load, dt, n):
    """Apply a tip load and damp the transient out; return the settled tip disp."""
    coords = _build_bar()
    tip = max(coords)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(tip, P_load, 0.0)
    ops.integrator(*integrator_args)
    ops.analysis("Transient")
    for _ in range(n):
        if ops.analyze(1, dt) != 0:
            return None
        if not math.isfinite(ops.nodeDisp(tip, 1)):
            return None
    return ops.nodeDisp(tip, 1)


def _max_rel_diff(a, b):
    assert a is not None and b is not None, "a leg failed to run"
    assert len(a) == len(b), "trace length mismatch (%d vs %d)" % (len(a), len(b))
    scale = max((abs(v) for v in b), default=0.0) or 1.0
    return max(abs(x - y) for x, y in zip(a, b)) / scale


# ---------------------------------------------------------------------------
# 1. Each new FLAG form == its DEPRECATED alias, byte-for-byte (rel < 1e-12).
#    (combo, flag-form integrator args, alias-form integrator args)
# ---------------------------------------------------------------------------
COMBOS = [
    ("lnvd",
     ("ExplicitBathe", P, "-lnvd", ALPHA),
     ("ExplicitBatheLNVD", P, ALPHA)),
    ("sms",
     ("ExplicitBathe", P, "-sms", DT_TARGET),
     ("ExplicitBatheSMS", P, DT_TARGET)),
    ("sms_consistent",
     ("ExplicitBathe", P, "-sms", DT_TARGET, "-consistent"),
     ("ExplicitBatheSMSConsistent", P, DT_TARGET)),
    ("lnvd_sms",
     ("ExplicitBathe", P, "-lnvd", ALPHA, "-sms", DT_TARGET),
     ("ExplicitBatheLNVDSMS", P, ALPHA, DT_TARGET)),
    ("lnvd_sms_consistent",
     ("ExplicitBathe", P, "-lnvd", ALPHA, "-sms", DT_TARGET, "-consistent"),
     ("ExplicitBatheLNVDSMSConsistent", P, ALPHA, DT_TARGET)),
]


@pytest.mark.parametrize("name,flag_args,alias_args", COMBOS, ids=[c[0] for c in COMBOS])
def test_flag_form_byte_identical_to_legacy_alias(name, flag_args, alias_args):
    dt = DT_TARGET if "sms" in name else 5.0e-4
    t_end = 4.0 if "sms" in name else 0.4
    u_flag = _free_vibration(flag_args, dt, t_end)
    u_alias = _free_vibration(alias_args, dt, t_end)
    rel = _max_rel_diff(u_flag, u_alias)
    assert rel < 1.0e-12, (
        "%s: -flag form must equal the legacy alias byte-for-byte; rel=%.3e" % (name, rel)
    )


# ---------------------------------------------------------------------------
# 2. LNVD x CONSISTENT cross-validation (the previously-untested combo, 33012):
#    FLAC local damping + matrix-free Olovsson PCG must COMPOSE.
# ---------------------------------------------------------------------------
def test_lnvd_consistent_supra_stable_no_nodal_mass():
    # (a) supra-stable at dtTarget (the un-scaled bar would diverge here).
    u = _free_vibration(("ExplicitBathe", P, "-lnvd", ALPHA, "-sms", DT_TARGET, "-consistent"),
                        DT_TARGET, 4.0)
    assert u is not None and max(abs(v) for v in u) < 1.0, "LNVD+consistent not stable"

    # (b) consistent path mutates NO nodal mass (matrix-free), even with LNVD active.
    _build_bar()
    ops.integrator("ExplicitBathe", P, "-lnvd", ALPHA, "-sms", DT_TARGET, "-consistent")
    ops.analysis("Transient")
    ops.analyze(1, DT_TARGET)
    after = ops.nodeMass(2)
    assert abs(after[0]) < 1e-9, (
        "consistent must NOT inject nodal mass even with -lnvd; got %r" % after
    )


def test_lnvd_consistent_relaxes_loaded_bar_to_static():
    """FLAC damping (through the consistent solve) drives a tip-loaded bar to PL/EA."""
    P_load = 2.0
    # static closed form: sum of element extensions = P*L_i/(E*A) along the chain.
    u_exact = sum(P_load * L / (E * A) for L in LENGTHS)
    u = _relax_to_static(
        ("ExplicitBathe", P, "-lnvd", ALPHA, "-sms", DT_TARGET, "-consistent"),
        P_load, DT_TARGET, 4000)
    assert u is not None, "LNVD+consistent relaxation did not run"
    assert abs(u - u_exact) / abs(u_exact) < 1.0e-3, (
        "LNVD+consistent settled at %.6g, expected static %.6g" % (u, u_exact)
    )


def test_lnvd_consistent_pcg_active(capfd):
    """The matrix-free PCG runs at the sub-steps even with LNVD layered on top."""
    import re
    _free_vibration(("ExplicitBathe", P, "-lnvd", ALPHA, "-sms", DT_TARGET,
                     "-consistent", "-verbose"), DT_TARGET, 0.2)
    out = capfd.readouterr()
    iters = [int(m) for m in re.findall(r"PCG\s+(\d+)\s+iters", out.err + out.out)]
    assert iters, "no PCG iteration report under -verbose (consistent PCG inactive?)"
    assert max(iters) <= 30, "PCG slow to converge with LNVD: %r" % iters


# ---------------------------------------------------------------------------
# 3. Flag gating reduces to the base path.
# ---------------------------------------------------------------------------
def test_lnvd_zero_alpha_reduces_to_base():
    """-lnvd 0.0 routes through the override but adds no damping => == plain ExplicitBathe."""
    dt, t_end = 5.0e-4, 0.4
    u_lnvd0 = _free_vibration(("ExplicitBathe", P, "-lnvd", 0.0), dt, t_end)
    u_base = _free_vibration(("ExplicitBathe", P), dt, t_end)
    rel = _max_rel_diff(u_lnvd0, u_base)
    assert rel < 1.0e-12, "-lnvd 0.0 must be byte-identical to plain ExplicitBathe; rel=%.3e" % rel


def test_consistent_without_sms_is_ignored():
    """-consistent alone (no -sms) is an invalid combo -> ignored, runs as plain EB."""
    dt, t_end = 5.0e-4, 0.4
    u_cons = _free_vibration(("ExplicitBathe", P, "-consistent"), dt, t_end)
    u_base = _free_vibration(("ExplicitBathe", P), dt, t_end)
    rel = _max_rel_diff(u_cons, u_base)
    assert rel < 1.0e-12, "-consistent without -sms should reduce to plain EB; rel=%.3e" % rel
