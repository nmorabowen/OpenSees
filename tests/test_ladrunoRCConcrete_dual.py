"""WP-263: `-dualTangent` (tangentMode 3) end to end, through the binary.

The kernel-level oracle is tests/test_ladrunoRCConcrete_dual_cpp.py (g++, no binary). This file
checks that the option reaches the production path in every host:

  * both materials of the family (LadrunoRCConcrete, LadrunoRCFiniteStrain) accept it;
  * the ELEMENT stiffness of a damaged LadrunoBrick equals the central FD of its resisting force
    with -dualTangent (small strain, -geom linear; and finite strain, LadrunoRCFiniteStrain
    -geom finite, through the Hencky / LogStrain seam), while the default tangent misses it by
    >= 1e-2 (break-the-gate: the default omits the projector rotation and the damage rate);
  * under -implex the prescribed stress path is the -implex path (the tangent changes, not the
    stress), and the returned material tangent is not the frozen-damage secant;
  * the PlateFiber view (ASDShellQ4 + LayeredShell, the sigma_33 condensation Newton) runs the
    cyclic interlock membrane-shear path with -dualTangent to the same resultants as the default.

Orchestrator: run after building (WP-263 does not run build.bat).
"""
import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

E, NU, KC = 30000.0, 0.2, 2.0 / 3.0
CE = [0.0, 0.0007, 0.0020, 0.0100]
CS = [0.0, 24.0,   30.0,   5.0]
CD = [0.0, 0.0,    0.25,   1.0 - 5.0 / 45.0]
TE = [0.0, 0.0001, 0.0010, 0.0040]
TS = [0.0, 3.0,    0.5,    0.0]
TD = [0.0, 0.0,    1.0 - 0.5 / 5.0, 0.999]
BACKBONE = ["-Ce", *CE, "-Cs", *CS, "-Cd", *CD, "-Te", *TE, "-Ts", *TS, "-Td", *TD, "-Kc", KC]

# distorted brick (the RCFiniteStrain battery geometry)
_NODES = {
    1: (0.00, 0.00, 0.00), 2: (1.00, 0.10, 0.05), 3: (1.10, 1.00, 0.00),
    4: (0.05, 0.95, 0.10), 5: (0.00, 0.05, 1.00), 6: (1.00, 0.00, 1.05),
    7: (1.05, 1.00, 1.10), 8: (0.00, 1.00, 0.95),
}
_CONN = [1, 2, 3, 4, 5, 6, 7, 8]
# tension + shear past the crack (dt ~ 0.5-0.7): principal axes off the frame, damage loading
FBAR_SMALL = [[1.0 + 3.0e-4, 1.0e-4, 0.0], [1.0e-4, 1.0 - 5.0e-5, 0.0], [0.0, 0.0, 1.0 - 5.0e-5]]


def _mat(name, tag, dual, extra=()):
    a = [name, tag, E, NU, *BACKBONE, *extra]
    if dual:
        a.append("-dualTangent")
    ops.nDMaterial(*a)


def _affine_disp(Fbar):
    u = np.zeros(24)
    for tag, (x, y, z) in _NODES.items():
        u[(tag - 1) * 3:(tag - 1) * 3 + 3] = (np.asarray(Fbar) - np.eye(3)) @ np.array([x, y, z])
    return u


def _impose(u, name, geom, dual):
    """Virgin brick, u imposed in ONE step, so the committed state is the trial at u from zero
    history: the resisting force is a pure function of u and its FD is the map's derivative."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _NODES.items():
        ops.node(tag, x, y, z)
    _mat(name, 1, dual)
    ops.element("LadrunoBrick", 1, *_CONN, 1, "-formulation", "std", "-geom", geom)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for tag in _NODES:
        for d in range(3):
            ops.sp(tag, d + 1, float(u[(tag - 1) * 3 + d]))
    ops.constraints("Lagrange")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 60, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "imposed-displacement solve failed"


def _k_vs_fd(name, geom, dual, h=1.0e-8):
    u = _affine_disp(FBAR_SMALL)
    _impose(u, name, geom, dual)
    dmg = list(ops.eleResponse(1, "material", 1, "damage"))
    K = np.array(ops.eleResponse(1, "stiff"), dtype=float).reshape(24, 24)
    Kfd = np.zeros((24, 24))
    for d in range(24):
        up = u.copy(); up[d] += h
        um = u.copy(); um[d] -= h
        _impose(up, name, geom, dual); fp = np.array(ops.eleForce(1), dtype=float)
        _impose(um, name, geom, dual); fm = np.array(ops.eleForce(1), dtype=float)
        Kfd[:, d] = (fp - fm) / (2.0 * h)
    return dmg, np.abs(K - Kfd).max() / np.abs(Kfd).max()


@pytest.mark.parametrize("name", ["LadrunoRCConcrete", "LadrunoRCFiniteStrain"])
def test_dual_tangent_parses(name):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _mat(name, 7, True)
    _mat(name, 8, True, ("-implex",))


@pytest.mark.parametrize("name,geom", [("LadrunoRCConcrete", "linear"), ("LadrunoRCFiniteStrain", "finite")])
def test_damaged_element_stiffness_matches_fd(name, geom):
    dmg, err_dual = _k_vs_fd(name, geom, dual=True)
    assert dmg[0] > 0.3, f"{name}: the state is not damaged (dt = {dmg[0]:.3f})"
    _, err_def = _k_vs_fd(name, geom, dual=False)
    print(f"\n{name} -geom {geom}: dt={dmg[0]:.3f}  |K - K_fd|/|K_fd|: dual {err_dual:.2e}, default {err_def:.2e}")
    assert err_dual <= 1.0e-5, f"{name}: dual element stiffness != FD ({err_dual:.3e})"
    assert err_def >= 1.0e-2, f"{name}: the default tangent matches the FD ({err_def:.3e}): the gate cannot bite"


# ---- IMPL-EX: the prescribed path is tangent-independent -----------------------------
_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}


def _drive_implex(dual, nsteps=40, exx=1.3e-3, gxy=0.8e-3):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y, z) in _CUBE.items():
        ops.node(t, float(x), float(y), float(z))
    _mat("LadrunoRCConcrete", 1, dual, ("-implex",))
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _CUBE.items():
        ops.sp(t, 1, exx * float(x) + gxy * float(y))   # u_x = exx X + gxy Y -> eps_xx, gamma_xy
        ops.sp(t, 2, 0.0); ops.sp(t, 3, 0.0)
    ops.constraints("Penalty", 1e16, 1e16); ops.numberer("Plain"); ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-10, 60, 0); ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps); ops.analysis("Static")
    sig = []
    for _ in range(nsteps):
        assert ops.analyze(1) == 0
        sig.append(list(ops.eleResponse(1, "material", 1, "stress")))
    C = np.array(ops.eleResponse(1, "material", 1, "tangent"), dtype=float).reshape(6, 6)
    dmg = list(ops.eleResponse(1, "material", 1, "damage"))
    return np.array(sig), C, dmg


def test_implex_dual_keeps_the_implex_stress_path():
    s0, C0, d0 = _drive_implex(dual=False)
    s3, C3, d3 = _drive_implex(dual=True)
    assert d3[0] > 0.3, f"no damage on the IMPL-EX path ({d3})"
    scale = np.abs(s0).max()
    assert np.abs(s3 - s0).max() <= 1e-9 * scale, "the dual tangent changed the IMPL-EX stress path"
    rel = np.abs(C3 - C0).max() / np.abs(C3).max()
    print(f"\nIMPL-EX: |C_dual - C_secant|/|C_dual| = {rel:.3e}")
    assert rel >= 1e-2, "the IMPL-EX dual tangent is the frozen-damage secant: mode 3 is not wired"


# ---- PlateFiber view: ASDShellQ4 + LayeredShell, the sigma_33 condensation --------------
_QUAD = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0), 4: (0.0, 1.0)}


def _shell_cyclic_shear(dual, exx=5.0e-4, gamma=2.0e-3, n1=40, nleg=50):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for t, (x, y) in _QUAD.items():
        ops.node(t, x, y, 0.0)
    _mat("LadrunoRCConcrete", 1, dual, ("-interlock", "-cyclic", "-agg", 16.0, "-crackSpacing", 50.0, "-beta"))
    ops.section("LayeredShell", 1, 4, 1, 0.025, 1, 0.025, 1, 0.025, 1, 0.025)
    ops.element("ASDShellQ4", 1, 1, 2, 3, 4, 1)
    vals = (0, 0, 1, -1, 1)
    times = [float(i) for i in range(len(vals) + 1)]
    ops.timeSeries("Path", 1, "-time", *times, "-values", *([0.0] + [1.0] * (len(vals) - 1) + [1.0]))
    ops.timeSeries("Path", 2, "-time", *times, "-values", *(list(vals) + [vals[-1]]))
    ops.pattern("Plain", 1, 1)
    for t, (x, y) in _QUAD.items():
        ops.sp(t, 1, exx * x)
        ops.sp(t, 3, 0.0); ops.sp(t, 4, 0.0); ops.sp(t, 5, 0.0); ops.sp(t, 6, 0.0)
    ops.pattern("Plain", 2, 2)
    for t, (x, y) in _QUAD.items():
        ops.sp(t, 2, gamma * x)
    ops.system("FullGeneral"); ops.numberer("Plain"); ops.constraints("Penalty", 1.0e12, 1.0e12)
    ops.test("NormDispIncr", 1.0e-9, 100, 0); ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / n1); ops.analysis("Static")
    for _ in range(n1):
        assert ops.analyze(1) == 0, "shell tension stage failed"
    out = []
    ops.integrator("LoadControl", 1.0 / nleg)
    for stage in range(1, len(vals) - 1):
        for _ in range(nleg):
            assert ops.analyze(1) == 0, f"shell shear stage {stage} failed (dual={dual})"
            ops.eleResponse(1, "forces")
            st = list(ops.eleResponse(1, "stresses"))
            out.append(st[:3])
    return np.array(out)


def test_platefiber_view_runs_cyclic_interlock_with_dual():
    a = _shell_cyclic_shear(dual=False)
    b = _shell_cyclic_shear(dual=True)
    assert a.shape == b.shape and len(a) > 50
    scale = np.abs(a).max()
    assert np.abs(a - b).max() <= 1e-6 * scale, "the PlateFiber resultants moved with the tangent"
