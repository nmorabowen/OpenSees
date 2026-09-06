"""LadrunoShellModifier (WP-91) -- G6 mass gate (Zone-A).

The `-mass` flag scales ONLY getRho() (mass per unit area): rho_outer = massMod *
rho_inner. It must leave the stiffness untouched, so for a cantilevered shell patch
with all other modifiers = 1:

    massMod = 0.5  =>  M_global -> 0.5*M_global,  K_global unchanged
                   =>  lambda = K/M doubles          (eigen() returns lambda = omega^2)
                   =>  omega (and hence frequency) scales by EXACTLY sqrt(2)

See test_ladrunoShellModifier_section.py for G1-G5, G7, G8 and the shared modifier/
congruence contract notes.
"""
import math

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


@pytest.fixture(scope="module", autouse=True)
def _build_stamp():
    if not hasattr(ops, "ladrunoBuild"):
        pytest.skip("openseespy wheel fallback in use -- LadrunoShellModifier needs the fork build")
    print(f"\n[ladrunoBuild] {ops.ladrunoBuild()}")


E0, NU0, H0, RHO0 = 200.0, 0.25, 0.2, 2.5
_NX = 3  # small cantilever strip: enough modes to be a meaningful eigen problem


def _mod_args(**mods):
    args = []
    for k, v in mods.items():
        args += [f"-{k}", float(v)]
    return args


def _bare_section(tag, E=E0, nu=NU0, h=H0, rho=RHO0):
    ops.section("ElasticMembranePlateSection", tag, E, nu, h, rho, 1.0)


def _wrap(tag, inner_tag, **mods):
    ops.section("LadrunoShellModifier", tag, inner_tag, *_mod_args(**mods))


def _strip_nodes(nx=_NX):
    return {(i, j): i * 2 + j + 1 for i in range(nx + 1) for j in (0, 1)}


def _build_cantilever(sec_builder, nx=_NX, ele="ShellMITC4"):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    nodes = _strip_nodes(nx)
    for (i, j), tag in nodes.items():
        ops.node(tag, float(i), float(j), 0.0)
    sec_builder(1)
    for i in range(nx):
        n1, n2, n3, n4 = nodes[(i, 0)], nodes[(i + 1, 0)], nodes[(i + 1, 1)], nodes[(i, 1)]
        ops.element(ele, i + 1, n1, n2, n3, n4, 1)
    for (i, j), tag in nodes.items():
        if i == 0:
            ops.fix(tag, 1, 1, 1, 1, 1, 1)
    return nodes


def _free_dof_K(nodes):
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    tip = [nodes[(_NX, 0)], nodes[(_NX, 1)]]
    for t in tip:
        ops.load(t, 0.0, 0.0, 0.001, 0.0, 0.0, 0.0)  # tiny load, only to trigger formTangent
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "G6 static tangent solve failed"
    K = np.array(ops.printA("-ret"))
    n = int(round(len(K) ** 0.5))
    return K.reshape(n, n)


def _first_frequency():
    """Matches the established repo idiom (tests/test_ladrunoplane_dynamics.py::_eigs):
    ops.eigen() runs against its own default analysis machinery, no explicit
    numberer/system/constraints call needed beforehand."""
    lam = ops.eigen("-fullGenLapack", 1)[0]
    return math.sqrt(lam)


def _run_mass_case(mass_mod, nx=_NX):
    def build(tag):
        _bare_section(tag + 100)
        _wrap(tag, tag + 100, f11=1.0, f22=1.0, f12=1.0, m11=1.0, m22=1.0, m12=1.0,
              v13=1.0, v23=1.0, mass=mass_mod)

    _build_cantilever(build, nx)
    omega1 = _first_frequency()

    # rebuild fresh (eigen() leaves the SOE in a state unsuited for printA) for the
    # stiffness-untouched check
    nodes = _build_cantilever(build, nx)
    K = _free_dof_K(nodes)
    return omega1, K


def test_g6_mass_modifier_scales_frequency_by_sqrt2():
    omega_full, K_full = _run_mass_case(1.0)
    omega_half, K_half = _run_mass_case(0.5)

    assert omega_full > 0.0
    ratio = omega_half / omega_full
    assert ratio == pytest.approx(math.sqrt(2.0), rel=1.0e-6), (
        f"massMod=0.5 frequency ratio {ratio} != sqrt(2) (omega_full={omega_full}, "
        f"omega_half={omega_half})"
    )
    # stiffness must be untouched by -mass
    np.testing.assert_allclose(K_half, K_full, rtol=1.0e-12, atol=1.0e-12,
                                err_msg="-mass changed the assembled stiffness")


def test_g6_mass_modifier_leaves_bare_section_frequency_unchanged_at_1():
    """massMod = 1.0 must reproduce the bare (unwrapped) inner section's frequency
    exactly -- a narrower restatement of the G1 identity property specialized to the
    dynamic (mass) path."""
    def wrapped(tag):
        _bare_section(tag + 100)
        _wrap(tag, tag + 100, mass=1.0)

    def bare(tag):
        _bare_section(tag)

    _build_cantilever(wrapped)
    omega_wrap = _first_frequency()
    _build_cantilever(bare)
    omega_bare = _first_frequency()

    assert omega_wrap == pytest.approx(omega_bare, rel=1.0e-11)
