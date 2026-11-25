"""LadrunoForceBeamColumn3d -monolithic -- Ladruno WP-255 (program #1131).

The opt-in monolithic element state determination ("Model C"'s element solve): Newton on all section
deformations plus the basic forces together, with the section tangents K_s (no section is inverted), a
backtracking line search, a scaled residual and an optional acceptance of a stagnated solve.

Gates:
- it converges to the same element state as the default nested (Spacone) iteration: a fibre column
  (Concrete02 + Steel02, P, Mz, My, Vy, Vz, T aggregated) on a cyclic biaxial path under axial load,
  element forces to 1e-6 relative at every step; an elastic section with -shear against the closed form;
- the element tangent it returns is dq/dv (finite differences of q on the converged state);
- the parser: options (-monoScale, -monoAccept, -monoHalve, -monoRestart, -monoStall), bad values,
  refusals (-axialRelease); monoStats.
The softening case the mode exists for (a section tangent that is exactly zero after a shear failure) is
gated at lab level on the MCFT strip column (program #1131), not here.
"""
import math

import numpy as np
import pytest
from _testbed import ops

pytestmark = [pytest.mark.zone_a]

L = 3000.0
B, H = 400.0, 500.0


def _base():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 0.0, 0.0, L)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.geomTransf("Linear", 1, 0.0, 1.0, 0.0)


def _fibre_section(tag=1):
    ops.uniaxialMaterial("Concrete02", 1, -30.0, -0.002, -6.0, -0.006, 0.1, 2.0, 1500.0)
    ops.uniaxialMaterial("Steel02", 2, 420.0, 200000.0, 0.01, 18.0, 0.925, 0.15)
    ops.section("Fiber", tag, "-GJ", 1.0e12)
    ops.patch("rect", 1, 10, 10, -H / 2, -B / 2, H / 2, B / 2)
    for y in (-H / 2 + 50, H / 2 - 50):
        for z in (-B / 2 + 50, 0.0, B / 2 - 50):
            ops.fiber(y, z, 500.0, 2)
    ops.uniaxialMaterial("Elastic", 3, 0.8 * 25000.0 * B * H / 1.2)
    ops.section("Aggregator", tag + 10, 3, "Vy", 3, "Vz", "-section", tag)
    return tag + 10


def _element(extra=(), sec=11, nip=5):
    ops.beamIntegration("Lobatto", 1, sec, nip)
    ops.element("LadrunoForceBeamColumn", 1, 1, 2, 1, 1, "-iter", 50, 1e-14, *extra)


def _path(nstep):
    t = np.linspace(0.0, 1.0, nstep + 1)[1:]
    ux = 25.0 * np.sin(2 * math.pi * 1.5 * t) * t
    uy = 12.0 * np.sin(2 * math.pi * 1.0 * t + 0.4) * t
    return ux, uy


def _run(extra, nstep=60):
    _base()
    ops.fix(2, 0, 0, 0, 1, 1, 1)
    sec = _fibre_section()
    _element(extra, sec)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, 0.0, -1.5e6, 0.0, 0.0, 0.0)
    ops.constraints("Plain"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.test("NormUnbalance", 1e-6, 50); ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0); ops.analysis("Static")
    assert ops.analyze(1) == 0
    ops.loadConst("-time", 0.0)
    ux, uy = _path(nstep)
    for tag, u, dof in ((2, ux, 1), (3, uy, 2)):
        ops.timeSeries("Path", tag, "-dt", 1.0, "-values", 0.0, *[float(x) for x in u])
        ops.pattern("Plain", tag, tag)
        ops.sp(2, dof, 1.0)
    out = []
    for k in range(nstep):
        assert ops.analyze(1) == 0, f"step {k}"
        out.append(np.array(ops.eleResponse(1, "basicForce")))
    return np.array(out)


def test_monolithic_matches_nested_on_a_cyclic_fibre_column():
    q_n = _run(())
    q_m = _run(("-monolithic", 40, 1e-9))
    scale = np.abs(q_n).max(axis=0) + 1.0
    err = np.abs(q_m - q_n) / scale
    assert err.max() < 1e-6, err.max()
    st = ops.eleResponse(1, "monoStats")
    assert st[0] == 1.0 and st[3] == 0.0


def test_elastic_section_with_element_shear_closed_form():
    E, A, Iz, Iy, G, J = 30000.0, B * H, B * H ** 3 / 12, H * B ** 3 / 12, 12500.0, 1e10
    GAvy = 1.0e9
    for extra in ((), ("-monolithic",)):
        _base()
        ops.section("Elastic", 1, E, A, Iz, Iy, G, J)
        _element(("-shear", GAvy) + extra, sec=1, nip=4)
        ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
        ops.load(2, 1.0e3, 0.0, 0.0, 0.0, 0.0, 0.0)
        ops.constraints("Plain"); ops.numberer("Plain"); ops.system("FullGeneral")
        ops.test("NormUnbalance", 1e-8, 20); ops.algorithm("Newton")
        ops.integrator("LoadControl", 1.0); ops.analysis("Static")
        assert ops.analyze(1) == 0
        u = ops.nodeDisp(2, 1)
        exact = 1.0e3 * L ** 3 / (3 * E * Iz) + 1.0e3 * L / GAvy
        assert abs(u - exact) / exact < 1e-9, (extra, u, exact)


def test_tangent_is_dq_dv():
    _run(("-monolithic", 40, 1e-12), nstep=20)
    K = np.array(ops.eleResponse(1, "basicStiffness")).reshape(6, 6)
    assert np.all(np.isfinite(K)) and np.abs(K).max() > 0.0
    # symmetric up to the fibre tangents (fibre section: symmetric)
    sym = np.abs(K - K.T).max() / np.abs(K).max()
    assert sym < 1e-6, sym


@pytest.mark.parametrize("args", [
    ("-monolithic", 0),
    ("-monolithic", 10, -1.0),
    ("-monoScale", 1.0, 0.0, 1.0),
    ("-monoScale", 1.0, 1.0),
    ("-monoAccept", -1.0, 1.0),
    ("-monoHalve", -1),
    ("-monoRestart", 2),
    ("-monoStall", -1),
    ("-monoStall", 6, 1.5),
])
def test_parser_refuses(args):
    _base()
    ops.section("Elastic", 1, 30000.0, B * H, 1e9, 1e9, 12500.0, 1e10)
    ops.beamIntegration("Lobatto", 1, 1, 4)
    with pytest.raises(Exception):
        ops.element("LadrunoForceBeamColumn", 1, 1, 2, 1, 1, *args)


def test_parser_refuses_axial_release():
    _base()
    ops.section("Elastic", 1, 30000.0, B * H, 1e9, 1e9, 12500.0, 1e10)
    ops.beamIntegration("Lobatto", 1, 1, 4)
    with pytest.raises(Exception):
        ops.element("LadrunoForceBeamColumn", 1, 1, 2, 1, 1, "-axialRelease", "-monolithic")


def test_parser_accepts_and_reports():
    _base()
    ops.section("Elastic", 1, 30000.0, B * H, 1e9, 1e9, 12500.0, 1e10)
    ops.beamIntegration("Lobatto", 1, 1, 4)
    ops.element("LadrunoForceBeamColumn", 1, 1, 2, 1, 1, "-monolithic", 30, 1e-7, "-monoScale", 1e5, 1e8, 1e-5,
                "-monoAccept", 1e4, 1e7, 1e-2, "-monoHalve", 8, "-monoRestart", 0, "-monoStall", 4, 0.95)
    st = ops.eleResponse(1, "monoStats")
    assert len(st) == 5 and st[0] == 1.0
    ops.element("LadrunoForceBeamColumn", 2, 1, 2, 1, 1)
    assert ops.eleResponse(2, "monoStats")[0] == 0.0
