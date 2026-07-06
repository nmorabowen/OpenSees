"""CentralDifferenceLadruno constant-mass tangent cache (ADR-67 P-NEW-1) —
Zone-A gate: the cache must be BIT-IDENTICAL to every-step reassembly
(-noMassCache) on the plain path AND across both invalidation channels the
adversarial gate flagged:

  MC-1  quiet-path bit-identity (the G-BYTE claim itself);
  MC-2  mid-run topology change (remove element -> domainChanged -> the cache
        must reassemble the new mass, not solve against the stale diagonal);
  MC-3  failed-step revert + retry at a smaller dt (revertToLastStep re-arms
        the first-step starter; the cache must be cleared WITH massBuilt so the
        starter really assembles — the gate's top-priority guard).

Each test runs the identical load/step sequence twice (cache default-ON vs
-noMassCache) and compares the full recorded displacement trajectory with
EXACT float equality (bit-identity, not a tolerance).

Theory/measurements: Ladruno_implementation/40b_phase0_dominance_report.md
(lane-D addendum: -factorOnce A/B, md5-identical, -17.9% wall) and the ADR-67
implementation log (P-NEW-1).
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

CDL = "CentralDifferenceLadruno"

# 1-D fixed-free truss bar: N elements, redundant DOUBLED element mid-bar so
# MC-2 can remove one without losing structural integrity.
N = 20
E, A, RHO = 1.0e4, 1.0, 1.0
H = 1.0
TIP = N + 1
REDUNDANT_TAG = 1000        # parallel twin of element N//2
DT = 0.004                  # well under dt_cr for this bar (~2h/c = 2/sqrt(E/rho)*h)
FAPP = 1.0


def _build(extra_int_args=()):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.uniaxialMaterial("Elastic", 1, E)
    for n in range(1, N + 2):
        ops.node(n, (n - 1) * H)
    ops.fix(1, 1)
    for e in range(1, N + 1):
        ops.element("Truss", e, e, e + 1, A, 1, "-rho", RHO)
    mid = N // 2
    ops.element("Truss", REDUNDANT_TAG, mid, mid + 1, A, 1, "-rho", RHO)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(TIP, FAPP)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(CDL, "-cfl", *extra_int_args)
    ops.analysis("Transient")


def _run_plain(int_args, nsteps=300):
    _build(int_args)
    traj = []
    for _ in range(nsteps):
        assert ops.analyze(1, DT) == 0
        traj.append(ops.nodeDisp(TIP, 1))
    return traj


def _run_topology(int_args, n1=100, n2=200):
    _build(int_args)
    traj = []
    for _ in range(n1):
        assert ops.analyze(1, DT) == 0
        traj.append(ops.nodeDisp(TIP, 1))
    # mid-run stiffness+mass change -> domainChange -> the cache MUST reassemble
    ops.remove("element", REDUNDANT_TAG)
    for _ in range(n2):
        assert ops.analyze(1, DT) == 0
        traj.append(ops.nodeDisp(TIP, 1))
    return traj


def _run_revert_retry(int_args, n1=100, n2=200, dt_big_factor=60.0):
    # arm the KE-growth circuit breaker so the oversized step FAILS within a
    # couple of attempts (pure NaN overflow needs ~100 steps to reach Inf);
    # a failed analyze drives Domain::revertToLastCommit + revertToLastStep —
    # the exact path the cache guard protects.
    _build(("-divergence", 5.0) + tuple(int_args))
    traj = []
    for _ in range(n1):
        assert ops.analyze(1, DT) == 0
        traj.append(ops.nodeDisp(TIP, 1))
    # drive the step far over dt_cr until the breaker trips the step
    # (deterministic: both runs fail at the same attempt), then retry small.
    failed = 0
    for _ in range(25):
        if ops.analyze(1, dt_big_factor * DT) != 0:
            failed += 1
            break
    assert failed == 1, "expected the oversized-dt step to fail (breaker/instability)"
    for _ in range(n2):
        assert ops.analyze(1, DT) == 0
        traj.append(ops.nodeDisp(TIP, 1))
    return traj


def _assert_bit_identical(t_cache, t_nocache, label):
    assert len(t_cache) == len(t_nocache)
    bad = [i for i, (a, b) in enumerate(zip(t_cache, t_nocache)) if a != b]
    assert not bad, (
        f"{label}: cache vs -noMassCache diverge at steps {bad[:5]} "
        f"(first: {t_cache[bad[0]]!r} vs {t_nocache[bad[0]]!r})")


def test_mc1_quiet_path_bit_identity():
    t_cache = _run_plain(())
    t_nocache = _run_plain(("-noMassCache",))
    _assert_bit_identical(t_cache, t_nocache, "MC-1 plain")
    # sanity: the bar actually moves
    assert max(abs(u) for u in t_cache) > 0.0


def test_mc2_topology_change_invalidates():
    t_cache = _run_topology(())
    t_nocache = _run_topology(("-noMassCache",))
    _assert_bit_identical(t_cache, t_nocache, "MC-2 remove-element")
    # sanity: the removal changed the response (pre vs post windows differ)
    pre = t_cache[50]
    assert t_cache[-1] != pre


def test_mc3_revert_retry_bit_identity():
    t_cache = _run_revert_retry(())
    t_nocache = _run_revert_retry(("-noMassCache",))
    _assert_bit_identical(t_cache, t_nocache, "MC-3 revert-retry")
