"""WP-267: DDM sensitivities of LadrunoBondSlip, end to end.

Before WP-267 LadrunoBondSlip had no setParameter and inherited
UniaxialMaterial's sensitivity stubs, so a bond law in a DDM run contributed
d(tau)/d(theta) = 0 (quirk WP-262-uniaxial-ddm-stubs-silent-zero). WP-267 adds
the parameters (tau_max s1 s2 s3 tau_f alpha Gf s0) and implements the
sensitivities with dual numbers on the templated kernel
(LadrunoBondSlipKernel.h): one Dual<1> pass seeded on the active parameter and
on the committed (Cslip, Cstress, CslipMaxAbs) sensitivities.

The gate: a 1-DOF zeroLength (the DOF is the slip, the force the bond stress)
with the bond law in parallel with an elastic zeroLength (10 > the softening
slope 4), under LOAD control: up into softening, an elastic unload, and a
reload past the old peak. sensNodeDisp(2, 1, 1) must equal a central FD of the
run at every step, for tau_max, s1, tau_f, alpha and s3. A Steel01 run of the
same harness is the control.

The protocol stays off the two branches where the material reports k0 instead
of the consistent tangent 0 (the plateau, here 1 < s <= 1.05, and sliding at
the cap): the DDM solve uses the reported tangent, so there it would not be
the derivative of the run. See tests/test_ladrunoBondSlip_dual_cpp.py (TAN).
"""
import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a, pytest.mark.t1]

BASE = dict(tau_max=10.0, s1=1.0, s2=1.05, s3=3.0, tau_f=2.0, alpha=0.4)
KE = 10.0


def _protocol():
    # F = tau(s) + KE s; F = 20 at the peak s1, 20.5 at the plateau end s2: no
    # load level lands in [19.9, 20.6] on a virgin branch, and the reload passes
    # the old peak (26.6) without revisiting it.
    up = 1.9 * np.arange(1, 15)                     # 1.9 .. 26.6
    down = 26.6 - 1.8 * np.arange(1, 8)             # 24.8 .. 14.0 (elastic unload)
    re = 14.0 + 1.7 * np.arange(1, 11)              # 15.7 .. 31.0 (past the old peak)
    return np.concatenate([up, down, re])


def _material(kind, vals):
    if kind == "Steel01":
        ops.uniaxialMaterial("Steel01", 1, vals["Fy"], vals["E"], 0.02)
        return
    ops.uniaxialMaterial("LadrunoBondSlip", 1, vals["tau_max"], vals["s1"], vals["s2"],
                         vals["s3"], vals["tau_f"], vals["alpha"])


def _run(kind, vals, sens_param=None):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 0.0)
    ops.fix(1, 1)
    _material(kind, vals)
    ops.uniaxialMaterial("Elastic", 2, KE)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1)
    ops.element("zeroLength", 2, 1, 2, "-mat", 2, "-dir", 1)
    prot = _protocol()
    ops.timeSeries("Path", 1, "-time", *range(len(prot) + 1),
                   "-values", 0.0, *[float(v) for v in prot])
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0)
    if sens_param is not None:
        ops.parameter(1, "element", 1, "material", 1, sens_param)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-14, 200)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    if sens_param is not None:
        ops.sensitivityAlgorithm("-computeAtEachStep")
    us, ss = [], []
    for k in range(len(prot)):
        assert ops.analyze(1) == 0, f"{kind}: step {k} failed"
        us.append(ops.nodeDisp(2, 1))
        ss.append(ops.sensNodeDisp(2, 1, 1) if sens_param is not None else 0.0)
    ops.wipe()
    return np.array(us), np.array(ss)


def _check(kind, base, par, rel_h=1.0e-6):
    _, ddm = _run(kind, base, par)
    h = rel_h * abs(base[par])
    vp = dict(base); vp[par] += h
    vm = dict(base); vm[par] -= h
    up, _ = _run(kind, vp)
    um, _ = _run(kind, vm)
    fd = (up - um) / (2.0 * h)
    scale = np.abs(fd).max()
    assert scale > 0.0, "the response does not depend on the parameter"
    return np.abs(ddm - fd).max() / scale, ddm, fd


BASE_S01 = dict(Fy=10.0, E=30.0)     # yields at F = 13.33: between two load levels


def test_protocol_softens_and_unloads_elastically():
    u, _ = _run("LadrunoBondSlip", BASE)
    prot = _protocol()
    tau = prot - KE * u                              # equilibrium: the bond stress
    i_pk = 13                                        # first peak, F = 26.6
    assert BASE["s2"] < u[i_pk] < BASE["s3"], "first peak not in softening"
    assert tau[i_pk] < 0.8 * BASE["tau_max"], "no real softening"
    assert np.all(np.diff(u[14:21]) < 0.0), "the unload does not unload"
    assert tau[20] > -tau[i_pk], "the unload reached the friction cap (inconsistent tangent)"
    assert u[-1] > u[i_pk] and u[-1] < BASE["s3"], "the reload does not pass the old peak"
    # never converged on the plateau (k0 reported there, not the consistent 0)
    assert not np.any((u > BASE["s1"]) & (u <= BASE["s2"]))


@pytest.mark.parametrize("par", ["Fy", "E"])
def test_control_steel01_harness(par):
    err, *_ = _check("Steel01", BASE_S01, par)
    assert err <= 1e-5, f"Steel01 control: DDM vs FD {err:.3e} -- the harness is wrong"


@pytest.mark.parametrize("par", ["tau_max", "s1", "tau_f", "alpha", "s3"])
def test_ladruno_bond_slip_ddm_matches_fd(par):
    err, ddm, fd = _check("LadrunoBondSlip", BASE, par)
    assert np.abs(ddm).max() > 0.0, "DDM sensitivity is identically zero (the pre-WP-267 stub)"
    assert err <= 1e-5, f"LadrunoBondSlip d(u)/d({par}): DDM vs FD {err:.3e}"
