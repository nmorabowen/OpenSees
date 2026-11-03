"""CentralDifferenceExplicit -commitSolveState: second state determination skipped.

Rate-independent models: bit-identical trajectories with and without the flag.
Materials with rate-dependent memory: a different, bounded trajectory, because the
committed state reflects v_{n+1/2} instead of v_{n+1}. Memoryless rate materials:
bit-identical, because the committed state is never read back. Velocity read at
force time from the nodes (Rayleigh force) and nodal velocity output stay identical,
because the nodes still receive the full-step velocity.
"""
import math

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

CDE = "CentralDifferenceExplicit"
DT = 0.004
N = 20
E, A, RHO = 1.0e4, 1.0, 1.0
TIP = N + 1
FAPP = 1.0


def _build_bar(int_args, viscous=None):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.uniaxialMaterial("Elastic", 1, E)
    ops.uniaxialMaterial("Viscous", 2, 5.0, 1.0)               # memoryless dashpot
    ops.uniaxialMaterial("ViscousDamper", 3, 1.0e4, 5.0, 1.0)  # Maxwell, with memory
    for n in range(1, N + 2):
        ops.node(n, float(n - 1))
    ops.fix(1, 1)
    n_ele = N - 1 if viscous == "truss-damper" else N
    for e in range(1, n_ele + 1):
        ops.element("Truss", e, e, e + 1, A, 1, "-rho", RHO)
    if viscous == "truss-damper":
        ops.element("Truss", N, N, N + 1, A, 3, "-rho", RHO)
    elif viscous in ("zl-damper", "zl-viscous"):
        ops.node(900, float(N))
        ops.fix(900, 1)
        mat = 3 if viscous == "zl-damper" else 2
        ops.element("zeroLength", 901, 900, TIP, "-mat", mat, "-dir", 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(TIP, FAPP)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(CDE, "-cfl", *int_args)
    ops.analysis("Transient")


def _run_bar(int_args, viscous=None, nsteps=300):
    _build_bar(int_args, viscous)
    traj = []
    for _ in range(nsteps):
        assert ops.analyze(1, DT) == 0
        traj.append(ops.nodeDisp(TIP, 1))
    return traj


def _assert_identical(a, b, label):
    bad = [i for i, (x, y) in enumerate(zip(a, b)) if x != y]
    assert not bad, "%s: diverges at steps %s" % (label, bad[:5])


def test_rate_independent_bit_identity():
    t_on = _run_bar(("-commitSolveState",))
    t_off = _run_bar(())
    _assert_identical(t_on, t_off, "elastic bar")
    assert max(abs(u) for u in t_on) > 0.0


@pytest.mark.parametrize("path,expect_delta", [
    ("zl-damper", True),
    ("truss-damper", True),
    ("zl-viscous", False),
])
def test_rate_dependent_characterization(path, expect_delta):
    t_on = _run_bar(("-commitSolveState",), viscous=path)
    t_off = _run_bar((), viscous=path)
    differs = any(x != y for x, y in zip(t_on, t_off))
    if expect_delta:
        assert differs, "%s: a delta is expected with rate-dependent memory" % path
        denom = max(abs(u) for u in t_off)
        delta = max(abs(x - y) for x, y in zip(t_on, t_off))
        assert delta < 0.05 * denom, "%s: delta %.3e, peak %.3e" % (path, delta, denom)
    else:
        assert not differs, "%s: a memoryless dashpot is trajectory-identical" % path


def _run_brick(int_args, nsteps=200):
    """Plastic stdBrick block with Rayleigh damping: the damping force reads nodal
    velocity at force time, the J2 state is rate independent, and the trajectory
    records nodal velocity as well, so the full-step output is covered."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    Eb, nu, rho = 2.0e5, 0.3, 8.0e-9
    K = Eb / (3 * (1 - 2 * nu))
    G = Eb / (2 * (1 + nu))
    ops.nDMaterial("J2Plasticity", 1, K, G, 250.0, 250.0, 0.0, 2000.0, rho)
    nx = ny = nz = 3
    L = 10.0

    def nid(i, j, k):
        return 1 + i + (nx + 1) * (j + (ny + 1) * k)

    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                ops.node(nid(i, j, k), i * L, j * L, k * L)
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.fix(nid(i, j, 0), 1, 1, 1)
    e = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                e += 1
                ops.element("stdBrick", e,
                            nid(i, j, k), nid(i+1, j, k), nid(i+1, j+1, k), nid(i, j+1, k),
                            nid(i, j, k+1), nid(i+1, j, k+1), nid(i+1, j+1, k+1),
                            nid(i, j+1, k+1), 1, "-lumped")
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    top = nid(nx // 2, ny // 2, nz)
    ops.load(top, 100.0, 0.0, -100.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.rayleigh(100.0, 0.0, 1.0e-8, 0.0)
    ops.integrator(CDE, "-cfl", *int_args)
    ops.analysis("Transient")
    c_dil = math.sqrt((K + 4.0 * G / 3.0) / rho)
    dt = 0.25 * L / c_dil
    traj = []
    for _ in range(nsteps):
        assert ops.analyze(1, dt) == 0
        traj.append(ops.nodeDisp(top, 3))
        traj.append(ops.nodeVel(top, 3))
    return traj


def test_velocity_reading_brick_bit_identity():
    t_on = _run_brick(("-commitSolveState",))
    t_off = _run_brick(())
    _assert_identical(t_on, t_off, "stdBrick")
    assert max(abs(u) for u in t_on) > 0.0
    assert all(math.isfinite(u) for u in t_on)


def test_topology_change_and_revert():
    """Element removal and a failed step keep the identity with and without the flag."""
    def run(int_args):
        _build_bar(("-divergence", 5.0) + tuple(int_args))
        ops.element("Truss", 1000, N // 2, N // 2 + 1, A, 1, "-rho", RHO)
        ops.domainChange()
        traj = []
        for _ in range(80):
            assert ops.analyze(1, DT) == 0
            traj.append(ops.nodeDisp(TIP, 1))
        ops.remove("element", 1000)
        for _ in range(80):
            assert ops.analyze(1, DT) == 0
            traj.append(ops.nodeDisp(TIP, 1))
        failed = 0
        for _ in range(25):
            if ops.analyze(1, 60.0 * DT) != 0:
                failed += 1
                break
        assert failed == 1
        for _ in range(80):
            assert ops.analyze(1, DT) == 0
            traj.append(ops.nodeDisp(TIP, 1))
        return traj

    t_on = run(("-commitSolveState",))
    t_off = run(())
    _assert_identical(t_on, t_off, "topology and revert")
