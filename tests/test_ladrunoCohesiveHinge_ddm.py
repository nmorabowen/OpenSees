"""WP-267: DDM sensitivities of LadrunoCohesiveHinge, end to end.

Before WP-267 LadrunoCohesiveHinge inherited UniaxialMaterial's sensitivity
stubs (getStressSensitivity = 0, commitSensitivity = -1): a DDM run returned
d(u)/d(theta) = 0 for Mc, Gf or Kpen without a warning (quirk
WP-262-uniaxial-ddm-stubs-silent-zero). WP-267 implements them with dual
numbers on the templated kernel (LadrunoCohesiveHingeKernel.h): one Dual<1>
pass seeded on the active parameter and on the committed kappaMax /
kappaMax_old sensitivities.

The gate: a rotational spring (1-DOF zeroLength: the DOF is the jump, the force
the moment) carrying the hinge, in parallel with an elastic zeroLength stiff
enough (5e4 > max |dM/dk| = aExp Mc = 2e4) to keep the total tangent positive,
so LOAD control can drive it deep into softening and back: up to +300, then
down to -400 (new damage on the negative side). The DDM sensitivity
sensNodeDisp(2, 1, 1) must equal a central FD of the whole run at every step,
for Mc, Gf and Kpen, on the exponential, linear and IMPL-EX laws. A Steel01
run of the same harness is the control: if it fails, the harness is wrong.
"""
import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a, pytest.mark.t1]

MC, GF = 100.0, 0.5                # kappa0 = 1e-5, exp aExp ~ 200, linear kappaf ~ 0.01
KPEN0 = 1000.0 * MC * MC / (2.0 * GF)
KE = 5.0e4                         # parallel elastic spring
P_UP, P_DOWN = 300.0, -400.0
N_UP, N_DOWN = 20, 40


def _protocol():
    up = np.linspace(P_UP / N_UP, P_UP, N_UP)
    down = np.linspace(P_UP + (P_DOWN - P_UP) / N_DOWN, P_DOWN, N_DOWN)
    return np.concatenate([up, down])


def _material(kind, vals):
    if kind == "Steel01":
        ops.uniaxialMaterial("Steel01", 1, vals["Fy"], vals["E"], 0.01)
        return
    args = [vals["Mc"], vals["Gf"]]
    if "Kpen" in vals:
        args += ["-penalty", vals["Kpen"]]
    shape, implex = kind
    args += ["-linear" if shape == "linear" else "-exp"]
    if implex:
        args += ["-implex"]
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 1, *args)


def _run(kind, vals, sens_param=None):
    """Run the cyclic protocol. Returns the jump per step and, if sens_param,
    the DDM d(jump)/d(param)."""
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
    ops.integrator("LoadControl", 1.0)       # pseudo-time; the Path series is the load
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
    vp = dict(base)
    vp[par] += h
    vm = dict(base)
    vm[par] -= h
    up, _ = _run(kind, vp)
    um, _ = _run(kind, vm)
    fd = (up - um) / (2.0 * h)
    scale = np.abs(fd).max()
    assert scale > 0.0, "the response does not depend on the parameter"
    return np.abs(ddm - fd).max() / scale, ddm, fd


BASE = dict(Mc=MC, Gf=GF)
BASE_K = dict(Mc=MC, Gf=GF, Kpen=KPEN0)    # explicit penalty: FD of Kpen = updateParameter(Kpen)
BASE_S01 = dict(Fy=MC, E=2.0e5)     # yields at F = 125: between two load levels
KINDS = {"exp": ("exp", False), "linear": ("linear", False), "exp_implex": ("exp", True)}


@pytest.mark.parametrize("name", sorted(KINDS))
def test_protocol_softens_both_ways(name):
    u, _ = _run(KINDS[name], BASE)
    prot = _protocol()
    m_hinge = prot - KE * u                 # equilibrium: the hinge moment
    kappa0 = MC / KPEN0
    assert u.max() > 100 * kappa0 and u.min() < -100 * kappa0, "the hinge never opens"
    i_up = int(np.argmax(u))
    assert m_hinge[i_up] < 0.6 * MC, "no real softening at the positive peak"
    assert u.min() < -1.2 * u.max(), "no new damage on the negative side"


@pytest.mark.parametrize("par", ["Fy", "E"])
def test_control_steel01_harness(par):
    """Steel01 has hand-coded DDM: the harness must reproduce its FD."""
    err, *_ = _check("Steel01", BASE_S01, par)
    assert err <= 1e-5, f"Steel01 control: DDM vs FD {err:.3e} -- the harness is wrong"


@pytest.mark.parametrize("name,par", [("exp", "Mc"), ("exp", "Gf"), ("exp", "Kpen"),
                                      ("linear", "Mc"), ("linear", "Gf"),
                                      ("exp_implex", "Mc"), ("exp_implex", "Gf")])
def test_ladruno_cohesive_hinge_ddm_matches_fd(name, par):
    base = BASE_K if par == "Kpen" else BASE
    err, ddm, fd = _check(KINDS[name], base, par)
    assert np.abs(ddm).max() > 0.0, "DDM sensitivity is identically zero (the pre-WP-267 stub)"
    assert err <= 1e-5, f"LadrunoCohesiveHinge {name} d(u)/d({par}): DDM vs FD {err:.3e}"
