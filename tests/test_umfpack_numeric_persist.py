"""UMFPACK Numeric-persist factorization reuse (ADR-40 rank 8/10) — Zone-A gate:
persisting the UMFPACK Numeric LU across solves (default ON) and reusing it while
the assembled matrix is unchanged must be BIT-IDENTICAL to refactorizing every
solve (`-noNumericPersist`, the legacy behavior) on every algorithm path and
across every invalidation channel the adversarial gate flagged:

  NP-1  full Newton (A reassembled every iter -> refactor every iter): reuse
        never fires, must still be identical (no regression on the no-reuse path);
  NP-2  ModifiedNewton (tangent formed once/step, iters trisolve-only): the
        PRIMARY reuse lane — identical trajectory + the run really iterated;
  NP-3  DisplacementControl (the dUhat intra-iteration reuse lane, lane E):
        identical load-displacement path through snap-through, incl. the
        newStep/update reference-displacement solves;
  NP-4  mid-run `remove element` -> domainChanged -> setSize frees Numeric +
        rebuilds Symbolic: identical, and the response actually changes across
        the removal (invalidation truly happened, not a stale-factor solve);
  NP-5  a step driven to non-convergence + retry: identical, proving the error
        path leaves factored=false / Numeric freed so the retry refactorizes.

Each test runs the identical sequence twice (persist default-ON vs
`-noNumericPersist`) and compares the full recorded trajectory with EXACT float
equality (bit-identity, not a tolerance).

Design / measurements: Ladruno_implementation/40c_umfpack_numeric_persist_design.md
and 40b_phase0_dominance_report.md (rank-8/10 factor-vs-trisolve addenda).
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# Shallow von Mises truss (corotational -> geometric nonlinearity forces Newton
# to iterate, so ModifiedNewton has a factorization to reuse). Snap-through
# limit load ~3.80. A redundant twin of bar 1 (parallel, same nodes) lets NP-4
# remove an element without losing structural integrity.
_H = 0.10
_SPAN = 1.0
_E = 1.0e4
_A = 1.0
_APEX = 2
_TWIN = 1000            # redundant parallel twin of element 1


def _build(redundant=False):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, -_SPAN, 0.0)
    ops.node(_APEX, 0.0, _H)
    ops.node(3, _SPAN, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(3, 1, 1)
    ops.fix(_APEX, 1, 0)                 # apex: X fixed by symmetry, Y free
    ops.uniaxialMaterial("Elastic", 1, _E)
    ops.element("corotTruss", 1, 1, _APEX, _A, 1)
    ops.element("corotTruss", 2, _APEX, 3, _A, 1)
    if redundant:
        ops.element("corotTruss", _TWIN, 1, _APEX, _A, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(_APEX, 0.0, -1.0)


def _solver(persist):
    # persist default-ON vs the legacy free-every-solve escape.
    if persist:
        ops.system("UmfPack")
    else:
        ops.system("UmfPack", "-noNumericPersist")
    ops.numberer("Plain")
    ops.constraints("Plain")


def _run_newton(persist, algo, lam_step=0.1, nsteps=30, redundant=False):
    _build(redundant)
    _solver(persist)
    ops.test("NormDispIncr", 1.0e-12, 50, 0)
    ops.algorithm(algo)
    ops.integrator("LoadControl", lam_step)
    ops.analysis("Static")
    traj = []
    for _ in range(nsteps):
        assert ops.analyze(1) == 0
        traj.append(ops.nodeDisp(_APEX, 2))
    return traj


def _run_dispcontrol(persist, du=-1.0e-3, nsteps=250):
    _build()
    _solver(persist)
    ops.test("NormDispIncr", 1.0e-12, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("DisplacementControl", _APEX, 2, du)
    ops.analysis("Static")
    traj = []
    for _ in range(nsteps):
        assert ops.analyze(1) == 0
        traj.append((ops.nodeDisp(_APEX, 2), ops.getLoadFactor(1)))
    return traj


def _run_topology(persist, n1=15, n2=15):
    _build(redundant=True)
    _solver(persist)
    ops.test("NormDispIncr", 1.0e-12, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.1)
    ops.analysis("Static")
    traj = []
    for _ in range(n1):
        assert ops.analyze(1) == 0
        traj.append(ops.nodeDisp(_APEX, 2))
    # mid-run stiffness change -> domainChange -> setSize -> Numeric MUST be
    # freed + Symbolic rebuilt (the sparsity pattern shrinks).
    ops.remove("element", _TWIN)
    for _ in range(n2):
        assert ops.analyze(1) == 0
        traj.append(ops.nodeDisp(_APEX, 2))
    return traj


def _run_revert_retry(persist, n1=20):
    # climb near the ~3.80 limit load with load control, then attempt one
    # oversized step that overshoots the limit (no equilibrium exists -> the
    # algorithm fails deterministically), then retry with small steps. A failed
    # analyze leaves the domain at the last committed state; the next tangent
    # reassembly (zeroA) must re-arm the factorization.
    _build()
    _solver(persist)
    ops.test("NormDispIncr", 1.0e-12, 12, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.1)
    ops.analysis("Static")
    traj = []
    for _ in range(n1):
        assert ops.analyze(1) == 0
        traj.append(ops.nodeDisp(_APEX, 2))
    # oversized step: jump the load far past the limit load in one increment.
    ops.integrator("LoadControl", 5.0)
    failed = ops.analyze(1)
    assert failed != 0, "expected the oversized load step to fail (past limit load)"
    # retry small
    ops.integrator("LoadControl", 0.05)
    for _ in range(n1):
        assert ops.analyze(1) == 0
        traj.append(ops.nodeDisp(_APEX, 2))
    return traj


def _assert_bit_identical(t_on, t_off, label):
    assert len(t_on) == len(t_off)
    bad = [i for i, (a, b) in enumerate(zip(t_on, t_off)) if a != b]
    assert not bad, (
        f"{label}: persist vs -noNumericPersist diverge at steps {bad[:5]} "
        f"(first: {t_on[bad[0]]!r} vs {t_off[bad[0]]!r})")


def test_np1_newton_bit_identity():
    t_on = _run_newton(True, "Newton")
    t_off = _run_newton(False, "Newton")
    _assert_bit_identical(t_on, t_off, "NP-1 Newton")
    assert min(t_on) < -1e-4                      # the apex actually moves down


def test_np2_modified_newton_reuse_bit_identity():
    t_on = _run_newton(True, "ModifiedNewton")
    t_off = _run_newton(False, "ModifiedNewton")
    _assert_bit_identical(t_on, t_off, "NP-2 ModifiedNewton")
    assert min(t_on) < -1e-4


def test_np3_displacement_control_bit_identity():
    t_on = _run_dispcontrol(True)
    t_off = _run_dispcontrol(False)
    _assert_bit_identical(t_on, t_off, "NP-3 DisplacementControl")
    # sanity: DisplacementControl drove the apex the prescribed distance and a
    # substantial NONLINEAR load developed (decreasing load increments => Newton
    # iterated => the dUhat newStep/update reference-displacement solves ran).
    disps = [u for u, _ in t_on]
    lams = [lf for _, lf in t_on]
    assert abs(disps[-1]) > 0.2 and max(lams) > 3.0
    inc = [lams[i + 1] - lams[i] for i in range(len(lams) - 1)]
    assert inc[0] != inc[-1]          # load-displacement curve is nonlinear


def test_np4_topology_change_invalidates():
    t_on = _run_topology(True)
    t_off = _run_topology(False)
    _assert_bit_identical(t_on, t_off, "NP-4 remove-element")
    # sanity: removing the redundant bar softened the structure (post-removal
    # increment is larger in magnitude than the matched pre-removal increment)
    pre = t_on[14] - t_on[13]
    post = t_on[15] - t_on[14]
    assert abs(post) > abs(pre)


def test_np5_revert_retry_bit_identity():
    t_on = _run_revert_retry(True)
    t_off = _run_revert_retry(False)
    _assert_bit_identical(t_on, t_off, "NP-5 revert-retry")
