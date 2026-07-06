"""CentralDifferenceLadruno -commitSolveState (ADR-67 P-NEW-2) — Zone-A gate.

Opt-in skip of the SECOND element constitutive pass in update() (the pass runs
at the SAME displacement, only node v/a differ — measured ~22% of explicit
wall, 40b lane-D). Adversarial-gate contract:

  CS-1  rate-INdependent deck: flag ON vs OFF bit-identical (the G-BYTE claim);
  CS-2  rate-DEPENDENT deck (ZeroLength+Viscous AND Truss+Viscous — the two
        gate-flagged trigger paths): flag ON must produce a DIFFERENT but
        bounded trajectory (documents the G-EQUIV surface — committed damper
        state reflects v_{n+1/2} instead of v_{n+1});
  CS-3  velocity-READING but rate-independent fork element (LadrunoBrick with
        -bulkViscosity + uri -hourglass viscous): flag ON vs OFF bit-identical —
        proves the retained setResponse invariant (velocity is consumed at
        force-time / commit from NODES, which still get v_{n+1});
  CS-4  mid-run remove element + failed-step revert/retry under the flag:
        ON vs OFF bit-identity is preserved through domainChanged and
        revertToLastStep (composition with the P-NEW-1 mass cache).

RemoveRecorder/collapse criteria were gate-verified CLEAR by source analysis
(they read trial DISPLACEMENT only, identical in both passes) — covered by
CS-1/CS-4 rather than a dedicated Collapse deck.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

CDL = "CentralDifferenceLadruno"
DT = 0.004
N = 20
E, A, RHO = 1.0e4, 1.0, 1.0
TIP = N + 1
FAPP = 1.0


# ---------------------------------------------------------------- 1-D bars
def _build_bar(int_args, viscous=None):
    """Fixed-free truss bar. viscous=None -> elastic (rate-independent);
    'zl-damper'/'truss-damper' -> Maxwell ViscousDamper (rate material WITH
    internal memory -> the G-EQUIV delta reaches the trajectory);
    'zl-viscous' -> memoryless pure dashpot (Viscous): its stress depends only
    on the CURRENT trial rate, so the committed state never feeds back and the
    trajectory stays bit-identical even under the flag."""
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.uniaxialMaterial("Elastic", 1, E)
    ops.uniaxialMaterial("Viscous", 2, 5.0, 1.0)          # memoryless dashpot
    ops.uniaxialMaterial("ViscousDamper", 3, 1.0e4, 5.0, 1.0)  # Maxwell (memory)
    for n in range(1, N + 2):
        ops.node(n, float(n - 1))
    ops.fix(1, 1)
    n_ele = N - 1 if viscous == "truss-damper" else N
    for e in range(1, n_ele + 1):
        ops.element("Truss", e, e, e + 1, A, 1, "-rho", RHO)
    if viscous == "truss-damper":
        # gate finding: Truss passes a velocity-derived rate to its material
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
    ops.integrator(CDL, "-cfl", *int_args)
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
    assert not bad, f"{label}: diverges at steps {bad[:5]}"


def test_cs1_rate_independent_bit_identity():
    t_on = _run_bar(("-commitSolveState",))
    t_off = _run_bar(())
    _assert_identical(t_on, t_off, "CS-1 elastic bar")
    assert max(abs(u) for u in t_on) > 0.0


@pytest.mark.parametrize("path,expect_delta", [
    ("zl-damper", True),     # Maxwell damper: committed internal state feeds back
    ("truss-damper", True),  # same via the Truss rate path
    ("zl-viscous", False),   # memoryless dashpot: committed state never read back
])
def test_cs2_rate_dependent_characterization(path, expect_delta):
    t_on = _run_bar(("-commitSolveState",), viscous=path)
    t_off = _run_bar((), viscous=path)
    differs = any(x != y for x, y in zip(t_on, t_off))
    if expect_delta:
        # the G-EQUIV surface is REAL for rate materials WITH memory...
        assert differs, (
            f"CS-2 {path}: expected a (small) G-EQUIV delta for a rate-memory deck")
        # ...but bounded: same physics to engineering accuracy (the committed
        # state shifts by O(dt * da * dF/dv) per step)
        denom = max(abs(u) for u in t_off)
        delta = max(abs(x - y) for x, y in zip(t_on, t_off))
        assert delta < 0.05 * denom, (
            f"CS-2 {path}: delta {delta:.3e} exceeds 5% of peak {denom:.3e} — "
            "the lagged-rate commit should be a small perturbation")
    else:
        # memoryless rate material: bit-identical even under the flag
        assert not differs, (
            f"CS-2 {path}: memoryless dashpot should be trajectory-identical")


# ---------------------------------------------------------- 3-D brick block
def _run_brick(int_args, flavor, nsteps=200):
    """Small LadrunoBrick block with a velocity-READING option the gate flagged.
    flavor='std-bv' -> std formulation + -bulkViscosity (force-time term);
    flavor='uri-hgv' -> uri formulation + -hourglass viscous (rate-form damping).
    (-bulkViscosity is std/bbar-only and -hourglass is uri-only, so they are
    exercised as separate variants.)"""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    Eb, nu, rho = 2.0e5, 0.3, 8.0e-9
    K = Eb / (3 * (1 - 2 * nu))
    G = Eb / (2 * (1 + nu))
    ops.nDMaterial("LadrunoJ2", 1, K, G, "-iso", "voce", 250.0, 0.0, 0.0, 2000.0,
                   "-rho", rho)
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
    if flavor == "std-bv":
        opts = ("-formulation", "std", "-bulkViscosity", 0.06, 1.2)
    else:
        opts = ("-formulation", "uri", "-hourglass", "viscous", 0.05)
    e = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                e += 1
                ops.element("LadrunoBrick", e,
                            nid(i, j, k), nid(i+1, j, k), nid(i+1, j+1, k), nid(i, j+1, k),
                            nid(i, j, k+1), nid(i+1, j, k+1), nid(i+1, j+1, k+1),
                            nid(i, j+1, k+1), 1, "-lumped", *opts)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    top = nid(nx // 2, ny // 2, nz)
    ops.load(top, 100.0, 0.0, -100.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(CDL, "-cfl", *int_args)
    ops.analysis("Transient")
    # dt from the dilatational wave speed (criticalTimeStep() is only available
    # after the first newStep computes it)
    import math
    c_dil = math.sqrt((K + 4.0 * G / 3.0) / rho)
    dt = 0.5 * L / c_dil
    traj = []
    for _ in range(nsteps):
        assert ops.analyze(1, dt) == 0
        traj.append(ops.nodeDisp(top, 3))
    return traj


@pytest.mark.parametrize("flavor", ["std-bv", "uri-hgv"])
def test_cs3_velocity_reading_brick_bit_identity(flavor):
    t_on = _run_brick(("-commitSolveState",), flavor)
    t_off = _run_brick((), flavor)
    _assert_identical(t_on, t_off, f"CS-3 brick {flavor}")
    assert max(abs(u) for u in t_on) > 0.0


def test_cs4_topology_and_revert_composition():
    """Remove-element + failed-step revert under the flag (composes P-NEW-1+2)."""
    def run(int_args):
        _build_bar(("-divergence", 5.0) + tuple(int_args))
        # redundant twin so removal keeps integrity
        ops.element("Truss", 1000, N // 2, N // 2 + 1, A, 1, "-rho", RHO)
        # rebuild analysis after adding the element (keeps both runs identical)
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
    _assert_identical(t_on, t_off, "CS-4 topology+revert")
