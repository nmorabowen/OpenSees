"""WP-260b — Explicit wave hygiene + LadrunoCPML threading (Zone-A gates).

Verifies every ADR-68 gate for WP-260b:
  1. ExplicitBathe P-NEW-1 mass cache vs -noMassCache (G-BYTE bit-identical on
     plain run, mid-run removeElement invalidation, failed-step revert+retry,
     and LadrunoProjection constraint-projector latch coupling).
  2. ExplicitBathe P-NEW-2 -commitSolveState (G-BYTE on rate-independent bar
     and velocity-reading LadrunoBrick -bulkViscosity / uri -hourglass viscous;
     bounded nonzero G-EQUIV delta on Maxwell ViscousDamper).
  3. ExplicitBathe -lnvd hoisted V_fake / lnvdDamped quasi-static convergence.
  4. ADR-67 P1.1 lifetime-invariant gate: mid-run removeElement + continue on
     both CentralDifferenceSMSConsistent and ExplicitBathe -sms -consistent
     (exercises ConsistentPCGScratch resize + rank-1 centroidal matvec).
  5. ADR-67 P1.2' massless-element short-circuit in CriticalTimeStep.cpp.
  6. LadrunoBrick formUri(0), bulk-viscosity c_d cache (incl. setParameter
     invalidation on both LadrunoBrick and LadrunoQuad), LadrunoBrick20 T7
     inertia skip, and lumped pure-alphaM O(n) fast paths.
  7. LadrunoConstraintProjector pre-factored Cholesky L (nRet == 1 and nRet > 1
     rigid-link momentum conservation and constraint manifold satisfaction).
  8. LadrunoCPML OpenMP threading (1T vs 4T bit-identity under both CDL and
     ExplicitBathe) + lazy K_eff_/C_eff_ materialization under implicit Newmark
     and revertToLastCommit.
"""
import math
import os
import tempfile
import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


def _bit_equal(a: float, b: float) -> bool:
    if a != b:
        return False
    return math.copysign(1.0, a) == math.copysign(1.0, b)


def _assert_bit_identical(a, b, label):
    assert len(a) == len(b), f"{label}: length mismatch {len(a)} vs {len(b)}"
    bad = [i for i, (x, y) in enumerate(zip(a, b)) if not _bit_equal(x, y)]
    assert not bad, (
        f"{label}: diverge at steps {bad[:5]} "
        f"(first: {a[bad[0]]!r} vs {b[bad[0]]!r})"
    )


# ----------------------------------------------------------------------------
# 1. ExplicitBathe P-NEW-1 mass cache (G-BYTE)
# ----------------------------------------------------------------------------
N_BAR = 20
E_BAR, A_BAR, RHO_BAR = 1.0e4, 1.0, 1.0
TIP_BAR = N_BAR + 1
REDUNDANT_ELE = 1000
DT_EB = 0.005


def _build_eb_bar(extra_int_args=(), viscous=None):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.uniaxialMaterial("Elastic", 1, E_BAR)
    ops.uniaxialMaterial("Viscous", 2, 5.0, 1.0)
    ops.uniaxialMaterial("ViscousDamper", 3, 1.0e4, 5.0, 1.0)
    for n in range(1, N_BAR + 2):
        ops.node(n, float(n - 1))
    ops.fix(1, 1)
    n_ele = N_BAR - 1 if viscous == "truss-damper" else N_BAR
    for e in range(1, n_ele + 1):
        ops.element("Truss", e, e, e + 1, A_BAR, 1, "-rho", RHO_BAR)
    if viscous is None:
        mid = N_BAR // 2
        ops.element("Truss", REDUNDANT_ELE, mid, mid + 1, A_BAR, 1, "-rho", RHO_BAR)
    elif viscous == "truss-damper":
        ops.element("Truss", N_BAR, N_BAR, N_BAR + 1, A_BAR, 3, "-rho", RHO_BAR)
    elif viscous == "zl-damper":
        ops.node(900, float(N_BAR))
        ops.fix(900, 1)
        ops.element("zeroLength", 901, 900, TIP_BAR, "-mat", 3, "-dir", 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(TIP_BAR, 1.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator("ExplicitBathe", 0.54, "-cfl", *extra_int_args)
    ops.analysis("Transient")


def test_eb_mc1_quiet_path_bit_identity():
    def run(args):
        _build_eb_bar(args)
        out = []
        for _ in range(200):
            assert ops.analyze(1, DT_EB) == 0
            out.append(ops.nodeDisp(TIP_BAR, 1))
        return out

    t_on = run(())
    t_off = run(("-noMassCache",))
    _assert_bit_identical(t_on, t_off, "EB-MC-1 plain")
    assert max(abs(u) for u in t_on) > 0.0


def test_eb_mc2_remove_element_invalidates_mass_cache():
    def run(args):
        _build_eb_bar(args)
        out = []
        for _ in range(80):
            assert ops.analyze(1, DT_EB) == 0
            out.append(ops.nodeDisp(TIP_BAR, 1))
        ops.remove("element", REDUNDANT_ELE)
        for _ in range(120):
            assert ops.analyze(1, DT_EB) == 0
            out.append(ops.nodeDisp(TIP_BAR, 1))
        return out

    t_on = run(())
    t_off = run(("-noMassCache",))
    _assert_bit_identical(t_on, t_off, "EB-MC-2 remove-element")
    assert t_on[-1] != t_on[40]


def test_eb_mc3_revert_retry_bit_identity():
    def run(args):
        _build_eb_bar(("-divergence", 50.0) + tuple(args))
        out = []
        for _ in range(80):
            assert ops.analyze(1, DT_EB) == 0
            out.append(ops.nodeDisp(TIP_BAR, 1))
        failed = 0
        for _ in range(30):
            if ops.analyze(1, 80.0 * DT_EB) != 0:
                failed += 1
                break
        assert failed == 1, "expected oversized ExplicitBathe step to trip breaker"
        for _ in range(120):
            assert ops.analyze(1, DT_EB) == 0
            out.append(ops.nodeDisp(TIP_BAR, 1))
        return out

    t_on = run(())
    t_off = run(("-noMassCache",))
    _assert_bit_identical(t_on, t_off, "EB-MC-3 revert-retry")


def test_eb_mc4_with_projection_handler_bit_identity():
    def run(args):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        ops.uniaxialMaterial("Elastic", 1, 1.0e4)
        ops.node(1, 0.0, 0.0)
        ops.node(2, 1.0, 0.0)
        ops.node(3, 1.0, 1.0)
        ops.node(4, 0.0, 1.0)
        ops.fix(1, 1, 1)
        ops.fix(4, 1, 1)
        ops.mass(2, 2.0, 2.0)
        ops.mass(3, 3.0, 3.0)
        ops.element("Truss", 1, 1, 2, 1.0, 1, "-rho", 1.0)
        ops.element("Truss", 2, 4, 3, 1.0, 1, "-rho", 1.0)
        ops.equalDOF(2, 3, 1)
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 10.0, 0.0)
        ops.constraints("LadrunoProjection")
        ops.numberer("Plain")
        ops.system("Diagonal")
        ops.test("NormDispIncr", 1e-12, 1)
        ops.algorithm("Linear")
        ops.integrator("ExplicitBathe", 0.54, *args)
        ops.analysis("Transient")
        out = []
        for _ in range(150):
            assert ops.analyze(1, 0.002) == 0
            u2 = ops.nodeDisp(2, 1)
            u3 = ops.nodeDisp(3, 1)
            assert abs(u2 - u3) < 1e-13
            out.append(u2)
        return out

    t_on = run(())
    t_off = run(("-noMassCache",))
    _assert_bit_identical(t_on, t_off, "EB-MC-4 projection")
    assert max(abs(u) for u in t_on) > 0.0


# ----------------------------------------------------------------------------
# 2. ExplicitBathe P-NEW-2 -commitSolveState
# ----------------------------------------------------------------------------
def test_eb_cs1_rate_independent_bit_identity():
    def run(args):
        _build_eb_bar(args)
        out = []
        for _ in range(150):
            assert ops.analyze(1, DT_EB) == 0
            out.append(ops.nodeDisp(TIP_BAR, 1))
        return out

    t_on = run(("-commitSolveState",))
    t_off = run(())
    _assert_bit_identical(t_on, t_off, "EB-CS-1 elastic bar")


def test_eb_cs2_rate_dependent_bounded_delta():
    def run(args):
        _build_eb_bar(args, viscous="zl-damper")
        out = []
        for _ in range(150):
            assert ops.analyze(1, DT_EB) == 0
            out.append(ops.nodeDisp(TIP_BAR, 1))
        return out

    t_on = run(("-commitSolveState",))
    t_off = run(())
    delta = max(abs(a - b) for a, b in zip(t_on, t_off))
    peak = max(abs(b) for b in t_off)
    assert delta > 0.0, "EB-CS-2 expected nonzero G-EQUIV delta on Maxwell damper"
    assert delta < 0.05 * peak, f"EB-CS-2 delta {delta:.3e} too large vs peak {peak:.3e}"


def test_eb_cs3_brick_velocity_reading_bit_identity():
    def run(int_args, flavor):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        Eb, nu, rho = 2.0e5, 0.3, 8.0e-9
        K = Eb / (3 * (1 - 2 * nu))
        G = Eb / (2 * (1 + nu))
        ops.nDMaterial("LadrunoJ2", 1, K, G, "-iso", "voce", 250.0, 0.0, 0.0, 2000.0, "-rho", rho)
        nx = ny = nz = 2
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
        etag = 1
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    ns = (
                        nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                        nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1),
                    )
                    if flavor == "std-bv":
                        ops.element("LadrunoBrick", etag, *ns, 1, "-formulation", "std",
                                    "-lumped", "-bulkViscosity", 0.06, 1.2)
                    else:
                        ops.element("LadrunoBrick", etag, *ns, 1, "-formulation", "uri",
                                    "-hourglass", "viscous", 0.1, "-lumped")
                    etag += 1
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        top = nid(1, 1, nz)
        ops.load(top, 5000.0, 0.0, 0.0)
        ops.constraints("Transformation")
        ops.numberer("Plain")
        ops.system("Diagonal")
        ops.test("NormDispIncr", 1e-12, 1)
        ops.algorithm("Linear")
        ops.integrator("ExplicitBathe", 0.54, *int_args)
        ops.analysis("Transient")
        out = []
        for _ in range(100):
            assert ops.analyze(1, 1.0e-6) == 0
            out.append(ops.nodeDisp(top, 1))
        return out

    for flavor in ("std-bv", "uri-hgv"):
        t_on = run(("-commitSolveState",), flavor)
        t_off = run((), flavor)
        _assert_bit_identical(t_on, t_off, f"EB-CS-3 {flavor}")
        assert max(abs(u) for u in t_on) > 0.0


# ----------------------------------------------------------------------------
# 3. ExplicitBathe -lnvd hoisted V_fake / lnvdDamped
# ----------------------------------------------------------------------------
def test_eb_lnvd_hoisted_converges_to_static():
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    k, m, f = 100.0, 1.0, 10.0
    ops.uniaxialMaterial("Elastic", 1, k)
    ops.node(1, 0.0)
    ops.node(2, 1.0)
    ops.fix(1, 1)
    ops.mass(2, m)
    ops.element("Truss", 1, 1, 2, 1.0, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, f)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator("ExplicitBathe", 0.54, "-lnvd", 0.8)
    ops.analysis("Transient")
    for _ in range(400):
        assert ops.analyze(1, 0.05) == 0
    u = ops.nodeDisp(2, 1)
    assert abs(u - f / k) < 1e-3, f"LNVD static convergence failed: u={u} vs {f/k}"


def test_eb_lnvd_transformation_equaldof():
    """ExplicitBathe -lnvd hoisted V_fake fast read under constraints Transformation
    with equalDOF (also verified 0-ULP vs base build by Ladruno_scripts/wp260b_ab_md5.py)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    k, m, f = 100.0, 1.0, 10.0
    ops.uniaxialMaterial("Elastic", 1, k)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 1.0, 0.0)
    ops.node(3, 1.0, 1.0)
    ops.node(4, 0.0, 1.0)
    ops.fix(1, 1, 1)
    ops.fix(4, 1, 1)
    ops.equalDOF(2, 3, 1, 2)
    ops.mass(2, m, m)
    ops.mass(3, m, m)
    ops.element("Truss", 1, 1, 2, 1.0, 1)
    ops.element("Truss", 2, 4, 3, 1.0, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, f, 0.0)
    ops.load(3, f, 0.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator("ExplicitBathe", 0.54, "-lnvd", 0.8)
    ops.analysis("Transient")
    for _ in range(400):
        assert ops.analyze(1, 0.05) == 0
    u2 = ops.nodeDisp(2)
    u3 = ops.nodeDisp(3)
    assert abs(u2[0] - u3[0]) < 1e-14 and abs(u2[1] - u3[1]) < 1e-14
    # Two parallel trusses of stiffness k under total load 2*f -> static u = (2*f)/(2*k) = f/k
    assert abs(u2[0] - f / k) < 1e-3, f"LNVD Transformation+equalDOF u={u2[0]} vs {f/k}"


# ----------------------------------------------------------------------------
# 4. ADR-67 P1.1 lifetime-invariant gate: consistent-SMS + mid-run removeElement
# ----------------------------------------------------------------------------
@pytest.mark.parametrize("integrator_kind", ["cdl_sms_cons", "eb_sms_cons"])
@pytest.mark.parametrize("change_neq", [False, True])
def test_sms_consistent_remove_element_lifetime_invariant(integrator_kind, change_neq):
    """Verify mid-run removal of a SCALED element (and optionally its free node so
    neq shrinks and ConsistentPCGScratch reallocates) produces a post-removal
    trajectory 0-ULP bit-identical to a fresh model built without it."""
    dt_target = 2.5e-4
    dt_run = 0.8 * dt_target
    n_warmup = 5
    n_post = 40
    # Delayed step load: zero during warmup (so blocks + ConsistentPCGScratch build
    # with element 99 active while u=v=a=0), then nonzero after removal.
    load_vals = [0.0] * (n_warmup + 1) + [1.0] * (n_post + 10)

    def run_case(mode):
        # mode: "remove_mid" (build with 99, warmup, remove 99, load),
        #       "fresh_without" (never add 99, same warmup + load),
        #       "keep_control" (keep 99 throughout — non-vacuity check).
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        ops.nDMaterial("ElasticIsotropic", 1, 2.0e7, 0.25, 2.0)

        xs = [0.0, 1.0, 1.2, 2.2]
        nid = {}
        t = 1
        for ix, x in enumerate(xs):
            for iy in range(2):
                for iz in range(2):
                    nid[(ix, iy, iz)] = t
                    ops.node(t, x, float(iy), float(iz))
                    if ix == 0:
                        ops.fix(t, 1, 1, 1)
                    t += 1
        for ix in range(3):
            ns = (
                nid[(ix, 0, 0)], nid[(ix + 1, 0, 0)], nid[(ix + 1, 1, 0)], nid[(ix, 1, 0)],
                nid[(ix, 0, 1)], nid[(ix + 1, 0, 1)], nid[(ix + 1, 1, 1)], nid[(ix, 1, 1)],
            )
            ops.element("LadrunoBrick", ix + 1, *ns, 1, "-formulation", "std", "-lumped")

        if mode != "fresh_without":
            # Stiff, short truss whose dt_cr = L / sqrt(E/rho) << dt_target so
            # buildMassScalingConsistent genuinely scales element 99!
            ops.uniaxialMaterial("Elastic", 2, 1.0e8)
            if change_neq:
                free_node = 99
                ops.node(free_node, 2.25, 0.0, 0.0)
                ops.element("Truss", 99, nid[(3, 0, 0)], free_node, 0.2, 2, "-rho", 2.0)
            else:
                ops.element("Truss", 99, nid[(1, 0, 0)], nid[(2, 0, 0)], 0.2, 2, "-rho", 2.0)

        ops.timeSeries("Path", 1, "-dt", dt_run, "-values", *load_vals)
        ops.pattern("Plain", 1, 1)
        ops.load(nid[(3, 1, 1)], 1000.0, 0.0, 0.0)
        ops.load(nid[(3, 0, 0)], 1000.0, 0.0, 0.0)

        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("Diagonal")
        ops.test("NormDispIncr", 1e-12, 1)
        ops.algorithm("Linear")
        if integrator_kind == "cdl_sms_cons":
            ops.integrator("CentralDifferenceSMSConsistent", dt_target,
                           "-maxAddedMass", 10.0, "-pcgTol", 1e-14)
        else:
            ops.integrator("ExplicitBathe", 0.54, "-sms", dt_target,
                           "-consistent", "-maxAddedMass", 10.0, "-pcgTol", 1e-14)
        ops.analysis("Transient")

        for _ in range(n_warmup):
            assert ops.analyze(1, dt_run) == 0

        if mode == "remove_mid":
            ops.remove("element", 99)
            if change_neq:
                ops.remove("node", 99)

        traj = []
        for _ in range(n_post):
            assert ops.analyze(1, dt_run) == 0
            traj.append(ops.nodeDisp(nid[(3, 0, 0)], 1))
            traj.append(ops.nodeDisp(nid[(3, 1, 1)], 1))
        return traj

    t_removed = run_case("remove_mid")
    t_fresh = run_case("fresh_without")
    t_kept = run_case("keep_control")

    _assert_bit_identical(
        t_removed, t_fresh,
        f"consistent-SMS removeElement ({integrator_kind}, change_neq={change_neq})",
    )
    peak = max(abs(x) for x in t_fresh)
    assert peak > 0.0
    max_ctrl_diff = max(abs(a - b) for a, b in zip(t_kept, t_fresh))
    assert max_ctrl_diff / peak > 1e-3, (
        f"non-vacuity: scaled element 99 must materially affect trajectory "
        f"(got rel diff {max_ctrl_diff / peak:.3e})"
    )


# ----------------------------------------------------------------------------
# 4b. Independent numpy Olovsson dense-matrix oracle for rank-1 consistent-SMS (G-EQUIV)
# ----------------------------------------------------------------------------
def test_sms_consistent_rank1_vs_independent_olovsson_oracle():
    """G-EQUIV oracle for default-ON rank-1 consistentBlockAccumulate:
    1) 1D 11-truss bar (oracle_olovsson_sms.py Case A) at -pcgTol 1e-14: compares
       step-by-step accelerations & displacements over 30 steps against dense
       np.linalg.solve(M_tilde, -K u).
    2) 3D 3-brick chain (ndmE=3, nn=8) at -pcgTol 1e-14: compares both
       CentralDifferenceSMSConsistent and ExplicitBathe -sms -consistent initial
       accelerations across all 36 free DOFs against an independent 8-pt Gauss
       numpy hex8 dense Olovsson solve."""
    # --- Part 1: 1D 11-truss bar multi-step oracle ---
    lengths = [1.0] * 5 + [0.05] + [1.0] * 5
    E, rho, A = 1.0e4, 1.0, 1.0
    dt_target = 0.012
    nfree = len(lengths)

    K = np.zeros((nfree, nfree))
    mlump = np.zeros(nfree)
    dM = np.zeros((nfree, nfree))
    for e, L in enumerate(lengths):
        ke = (E * A / L) * np.array([[1.0, -1.0], [-1.0, 1.0]])
        me = np.array([0.5 * rho * A * L, 0.5 * rho * A * L])
        dt_e = L / math.sqrt(E / rho)
        beta = (dt_target / dt_e) ** 2 - 1.0 if dt_e < dt_target else 0.0
        Mbar = beta * (np.diag(me) - np.outer(me, me) / me.sum())
        dofs = [e - 1, e]
        for a in range(2):
            if dofs[a] >= 0:
                mlump[dofs[a]] += me[a]
            for b in range(2):
                if dofs[a] >= 0 and dofs[b] >= 0:
                    K[dofs[a], dofs[b]] += ke[a, b]
                    dM[dofs[a], dofs[b]] += Mbar[a, b]
    M_tilde = np.diag(mlump) + dM

    xs_bar = np.cumsum([0.0] + lengths)
    u_ref = 1.0e-3 * xs_bar[1:]
    a_ref = np.linalg.solve(M_tilde, -(K @ u_ref))
    v_half = 0.5 * dt_target * a_ref

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 1, 1)
    for i, x in enumerate(xs_bar[1:], start=2):
        ops.node(i, float(x), 0.0)
        ops.fix(i, 0, 1)
    ops.uniaxialMaterial("Elastic", 1, E)
    for i in range(nfree):
        ops.element("Truss", i + 1, i + 1, i + 2, A, 1, "-rho", rho)
    # Apply initial displacement AFTER element creation so Truss::setDomain uses
    # the undeformed reference length L.
    for i, x in enumerate(xs_bar[1:], start=2):
        ops.setNodeDisp(i, 1, 1.0e-3 * float(x), "-commit")
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator("CentralDifferenceSMSConsistent", dt_target,
                   "-pcgTol", 1e-14, "-pcgMaxIt", 500)
    ops.analysis("Transient")

    max_rel_a_1d = 0.0
    max_rel_u_1d = 0.0
    for _ in range(30):
        u_ref = u_ref + dt_target * v_half
        a_ref = np.linalg.solve(M_tilde, -(K @ u_ref))
        v_half = v_half + dt_target * a_ref

        assert ops.analyze(1, dt_target) == 0
        u_ops = np.array([ops.nodeDisp(i + 2, 1) for i in range(nfree)])
        a_ops = np.array([ops.nodeAccel(i + 2, 1) for i in range(nfree)])
        max_rel_u_1d = max(max_rel_u_1d, np.linalg.norm(u_ops - u_ref) / np.linalg.norm(u_ref))
        max_rel_a_1d = max(max_rel_a_1d, np.linalg.norm(a_ops - a_ref) / np.linalg.norm(a_ref))

    assert max_rel_a_1d < 1e-12, f"1D Olovsson accel rel err {max_rel_a_1d:.3e} >= 1e-12"
    assert max_rel_u_1d < 1e-12, f"1D Olovsson disp rel err {max_rel_u_1d:.3e} >= 1e-12"
    # Non-vacuity: unscaled lumped solve differs hugely from Olovsson
    a_lump_only = -(K @ u_ref) / mlump
    assert np.linalg.norm(a_ref - a_lump_only) / np.linalg.norm(a_ref) > 0.1

    # --- Part 2: 3D Hex8 (ndmE=3, nn=8) dense Olovsson oracle ---
    E3, nu3, rho3 = 2.0e7, 0.25, 2.0
    lam3 = E3 * nu3 / ((1.0 + nu3) * (1.0 - 2.0 * nu3))
    mu3 = E3 / (2.0 * (1.0 + nu3))
    C3 = np.array([
        [lam3 + 2 * mu3, lam3, lam3, 0, 0, 0],
        [lam3, lam3 + 2 * mu3, lam3, 0, 0, 0],
        [lam3, lam3, lam3 + 2 * mu3, 0, 0, 0],
        [0, 0, 0, mu3, 0, 0],
        [0, 0, 0, 0, mu3, 0],
        [0, 0, 0, 0, 0, mu3],
    ], dtype=float)

    def hex8_ke_me(lx, ly=1.0, lz=1.0):
        # Reference node order matching LadrunoBrick:
        # (-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),(-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1)
        signs = np.array([
            [-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
            [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1],
        ], dtype=float)
        gp = 1.0 / math.sqrt(3.0)
        Ke = np.zeros((24, 24))
        detJ = 0.125 * lx * ly * lz
        for xi in (-gp, gp):
            for eta in (-gp, gp):
                for zeta in (-gp, gp):
                    B = np.zeros((6, 24))
                    for a in range(8):
                        sx, sy, sz = signs[a]
                        dNdx = 0.125 * sx * (1 + sy * eta) * (1 + sz * zeta) * (2.0 / lx)
                        dNdy = 0.125 * (1 + sx * xi) * sy * (1 + sz * zeta) * (2.0 / ly)
                        dNdz = 0.125 * (1 + sx * xi) * (1 + sy * eta) * sz * (2.0 / lz)
                        c = 3 * a
                        B[0, c] = dNdx
                        B[1, c + 1] = dNdy
                        B[2, c + 2] = dNdz
                        B[3, c] = dNdy
                        B[3, c + 1] = dNdx
                        B[4, c + 1] = dNdz
                        B[4, c + 2] = dNdy
                        B[5, c] = dNdz
                        B[5, c + 2] = dNdx
                    Ke += (B.T @ C3 @ B) * detJ
        m_node = rho3 * lx * ly * lz / 8.0
        me = np.full(24, m_node)
        s = 1.0 / np.sqrt(me)
        lam_max = np.linalg.eigvalsh((Ke * s[:, None]) * s[None, :]).max()
        dt_e = 2.0 / math.sqrt(lam_max)
        return me, dt_e

    xs3 = [0.0, 1.0, 0.2, 1.0]  # lengths: 1.0, 0.2 (throttling!), 1.0
    xcoords = [0.0, 1.0, 1.2, 2.2]
    dt_target3 = 2.5e-4
    # Map (ix, iy, iz) -> node tag 1..16; ix==0 fixed -> free nodes are tags 5..16 (36 DOFs)
    nid = {}
    t = 1
    for ix in range(4):
        for iy in range(2):
            for iz in range(2):
                nid[(ix, iy, iz)] = t
                t += 1

    nfree3 = 36
    M_tilde3 = np.zeros((nfree3, nfree3))
    for ix, lx in enumerate([1.0, 0.2, 1.0]):
        me24, dt_e = hex8_ke_me(lx)
        beta = (dt_target3 / dt_e) ** 2 - 1.0 if dt_e < dt_target3 else 0.0
        ns = (
            nid[(ix, 0, 0)], nid[(ix + 1, 0, 0)], nid[(ix + 1, 1, 0)], nid[(ix, 1, 0)],
            nid[(ix, 0, 1)], nid[(ix + 1, 0, 1)], nid[(ix + 1, 1, 1)], nid[(ix, 1, 1)],
        )
        m_node = me24[0]
        m_tot = 8.0 * m_node
        for d in range(3):
            for a, na in enumerate(ns):
                ea = 3 * (na - 5) + d if na >= 5 else -1
                if ea >= 0:
                    M_tilde3[ea, ea] += m_node
                if beta > 0.0 and ea >= 0:
                    for b, nb in enumerate(ns):
                        eb = 3 * (nb - 5) + d if nb >= 5 else -1
                        if eb >= 0:
                            M_tilde3[ea, eb] += beta * ((m_node if a == b else 0.0) - m_node * m_node / m_tot)

    F_ext = np.zeros(nfree3)
    loads_3d = {
        nid[(3, 0, 0)]: (1200.0, -300.0, 150.0),
        nid[(3, 1, 0)]: (800.0, 400.0, -200.0),
        nid[(3, 1, 1)]: (1500.0, -100.0, 500.0),
        nid[(3, 0, 1)]: (600.0, 250.0, -350.0),
    }
    for nt, fv in loads_3d.items():
        for d in range(3):
            F_ext[3 * (nt - 5) + d] = fv[d]

    a_oracle3 = np.linalg.solve(M_tilde3, F_ext)
    dt_step3 = 1.0e-4

    for kind in ("cdl_sms_cons", "eb_sms_cons"):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        ops.nDMaterial("ElasticIsotropic", 1, E3, nu3, rho3)
        for ix, x in enumerate(xcoords):
            for iy in range(2):
                for iz in range(2):
                    nt = nid[(ix, iy, iz)]
                    ops.node(nt, x, float(iy), float(iz))
                    if ix == 0:
                        ops.fix(nt, 1, 1, 1)
        for ix in range(3):
            ns = (
                nid[(ix, 0, 0)], nid[(ix + 1, 0, 0)], nid[(ix + 1, 1, 0)], nid[(ix, 1, 0)],
                nid[(ix, 0, 1)], nid[(ix + 1, 0, 1)], nid[(ix + 1, 1, 1)], nid[(ix, 1, 1)],
            )
            ops.element("LadrunoBrick", ix + 1, *ns, 1, "-formulation", "std", "-lumped")
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        for nt, fv in loads_3d.items():
            ops.load(nt, *fv)
        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("Diagonal")
        ops.test("NormDispIncr", 1e-12, 1)
        ops.algorithm("Linear")
        p_bathe = 0.54
        if kind == "cdl_sms_cons":
            ops.integrator("CentralDifferenceSMSConsistent", dt_target3,
                           "-maxAddedMass", 10.0, "-pcgTol", 1e-14, "-pcgMaxIt", 500)
        else:
            ops.integrator("ExplicitBathe", p_bathe, "-sms", dt_target3,
                           "-consistent", "-maxAddedMass", 10.0, "-pcgTol", 1e-14, "-pcgMaxIt", 500)
        ops.analysis("Transient")
        assert ops.analyze(1, dt_step3) == 0

        u1 = np.zeros(nfree3)
        for nt in range(5, 17):
            u1[3 * (nt - 5):3 * (nt - 5) + 3] = ops.nodeDisp(nt)
        if kind == "cdl_sms_cons":
            # From rest (u0=v0=0): u1 = 0.5 * dt^2 * a0
            a_ops3 = 2.0 * u1 / (dt_step3 ** 2)
        else:
            # ExplicitBathe from rest (u0=v0=a0=0): substep 1 solves M_tilde a_p = F_ext,
            # then u1 = 0.5 * (1 - p) * dt^2 * a_p
            a_ops3 = 2.0 * u1 / ((1.0 - p_bathe) * (dt_step3 ** 2))
        rel_err3 = np.linalg.norm(a_ops3 - a_oracle3) / np.linalg.norm(a_oracle3)
        assert rel_err3 < 1e-12, f"3D Hex8 {kind} Olovsson oracle rel err {rel_err3:.3e} >= 1e-12"


# ----------------------------------------------------------------------------
# 5. ADR-67 P1.2' massless / non-positive-mass short-circuit in CriticalTimeStep
# ----------------------------------------------------------------------------
def test_cts_massless_element_short_circuit():
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    E, A, rho, L = 1.0e4, 1.0, 1.0, 1.0
    ops.uniaxialMaterial("Elastic", 1, E)
    ops.node(1, 0.0)
    ops.node(2, L)
    ops.node(3, L)
    ops.fix(1, 1)
    ops.fix(3, 1)
    ops.element("Truss", 1, 1, 2, A, 1, "-rho", rho)
    # Massless zeroLength springs — hit the P1.2' short-circuit before DGGEV
    for k in range(10):
        ops.element("zeroLength", 100 + k, 2, 3, "-mat", 1, "-dir", 1)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno", "-cfl")
    ops.analysis("Transient")
    ops.initialize()
    dt_cr = ops.criticalTimeStep()
    expected = L / math.sqrt(E / rho)
    assert abs(dt_cr - expected) / expected < 1e-10, f"dt_cr={dt_cr} vs expected={expected}"

    # Prove the short-circuit is ACTUALLY taken (fails on base build's DGGEV path):
    # an element with all-negative mdiag (-rho < 0) and negative stiffness (-E < 0)
    # has a spurious positive generalized eigenvalue in DGGEV (dt = 0.001 < 0.01),
    # whereas the non-positive-mass short-circuit skips it and keeps element #1's dt=0.01.
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.uniaxialMaterial("Elastic", 1, E)
    ops.uniaxialMaterial("Elastic", 2, -1.0e6)
    ops.node(1, 0.0)
    ops.node(2, L)
    ops.node(3, 2.0 * L)
    ops.fix(1, 1)
    ops.fix(3, 1)
    ops.element("Truss", 1, 1, 2, A, 1, "-rho", rho)
    ops.element("Truss", 2, 2, 3, A, 2, "-rho", -rho)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno", "-cfl")
    ops.analysis("Transient")
    ops.initialize()
    dt_cr_neg = ops.criticalTimeStep()
    assert abs(dt_cr_neg - expected) / expected < 1e-10, (
        f"CTS short-circuit not taken: got dt_cr={dt_cr_neg} instead of {expected}"
    )


# ----------------------------------------------------------------------------
# 6. Bulk-viscosity c_d cache + setParameter invalidation + alphaM fast path
# ----------------------------------------------------------------------------
def test_quad_bulk_viscosity_and_alphaM_fast_path():
    def run_quad(with_param_update=False, alpha_m=0.0):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        ops.nDMaterial("ElasticIsotropic", 1, 1.0e7, 0.25, 2.0)
        ops.node(1, 0.0, 0.0)
        ops.node(2, 1.0, 0.0)
        ops.node(3, 1.0, 1.0)
        ops.node(4, 0.0, 1.0)
        ops.fix(1, 1, 1)
        ops.fix(2, 1, 1)
        ops.element("LadrunoQuad", 1, 1, 2, 3, 4, 1, "-type", "PlaneStrain", "-thick", 1.0,
                    "-formulation", "std", "-bulkViscosity", 0.06, 1.2)
        if with_param_update:
            ops.parameter(1, "element", 1, "E")
        if alpha_m > 0.0:
            ops.rayleigh(alpha_m, 0.0, 0.0, 0.0)
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(3, 100.0, 0.0)
        ops.load(4, 100.0, 0.0)
        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("Diagonal")
        ops.test("NormDispIncr", 1e-12, 1)
        ops.algorithm("Linear")
        ops.integrator("CentralDifferenceLadruno")
        ops.analysis("Transient")
        for _ in range(30):
            assert ops.analyze(1, 1.0e-4) == 0
        if with_param_update:
            ops.updateParameter(1, 4.0e7)
        for _ in range(30):
            assert ops.analyze(1, 1.0e-4) == 0
        return ops.nodeDisp(3, 1)

    u_plain = run_quad(False, 0.0)
    u_updated = run_quad(True, 0.0)
    u_damped = run_quad(False, 20.0)
    assert math.isfinite(u_plain) and abs(u_plain) > 0.0
    assert u_updated != u_plain, "updateParameter(E) must invalidate Quad c_d / stiffness"
    assert abs(u_damped) < abs(u_plain), "alphaM fast path must damp response"


def test_brick20_t7_and_alphaM_lumped():
    def run_b20(alpha_m):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        ops.nDMaterial("ElasticIsotropic", 1, 1.0e7, 0.25, 2.0)
        # Standard 20-node hex on [0,1]^3
        corners = [
            (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
            (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1),
        ]
        mids = [
            (0.5, 0, 0), (1, 0.5, 0), (0.5, 1, 0), (0, 0.5, 0),
            (0.5, 0, 1), (1, 0.5, 1), (0.5, 1, 1), (0, 0.5, 1),
            (0, 0, 0.5), (1, 0, 0.5), (1, 1, 0.5), (0, 1, 0.5),
        ]
        for idx, (x, y, z) in enumerate(corners + mids, start=1):
            ops.node(idx, float(x), float(y), float(z))
            if z == 0.0:
                ops.fix(idx, 1, 1, 1)
        ops.element("LadrunoBrick20", 1, *range(1, 21), 1, "-formulation", "std", "-lumped")
        if alpha_m > 0.0:
            ops.rayleigh(alpha_m, 0.0, 0.0, 0.0)
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(7, 100.0, 0.0, 0.0)
        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("Diagonal")
        ops.test("NormDispIncr", 1e-12, 1)
        ops.algorithm("Linear")
        ops.integrator("CentralDifferenceLadruno")
        ops.analysis("Transient")
        for _ in range(40):
            assert ops.analyze(1, 5.0e-5) == 0
        return ops.nodeDisp(7, 1)

    u0 = run_b20(0.0)
    ud = run_b20(50.0)
    assert math.isfinite(u0) and abs(u0) > 0.0
    assert abs(ud) < abs(u0), "Brick20 lumped alphaM fast path must damp response"


# ----------------------------------------------------------------------------
# 7. LadrunoConstraintProjector pre-factored Cholesky L (G-EQUIV for nRet > 1)
# ----------------------------------------------------------------------------
def test_projector_cholesky_rigid_link_manifold():
    """Verify nRet=6 Cholesky solve in LadrunoConstraintProjector against:
    (1) independent numpy np.linalg.solve(L^T M L, L^T F_ext) at step 1 (< 1e-13 rel);
    (2) rigid-body generalized momentum balance L^T (M a) = L^T (F_ext - F_beam)
        after 100 steps (< 1e-11 rel);
    (3) rigid-link kinematic manifold."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 2.0, 0.0, 0.0)
    ops.node(3, 2.0, 1.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    m2 = np.array([5.0, 5.0, 5.0, 0.5, 0.5, 0.5])
    m3 = np.array([3.0, 3.0, 3.0, 0.3, 0.3, 0.3])
    ops.mass(2, *m2.tolist())
    ops.mass(3, *m3.tolist())
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.element("elasticBeamColumn", 1, 1, 2, 1.0, 1.0e5, 4.0e4, 1.0, 0.1, 0.1, 1)
    # Rigid beam link couples 6 retained DOFs of node 2 -> exercises nRet=6 Cholesky L!
    ops.rigidLink("beam", 2, 3)
    f3 = np.array([10.0, 20.0, 5.0, 0.0, 0.0, 0.0])
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, *f3.tolist())
    ops.constraints("LadrunoProjection")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno")
    ops.analysis("Transient")

    # Independent numpy oracle for step 1 from rest (u0=0, v0=0 => F_int=0):
    # r = x3 - x2 = (0, 1, 0) => a3 = [a2x - r_y*rz2, a2y, a2z + r_y*rx2, rx2, ry2, rz2]
    L = np.zeros((12, 6))
    L[:6, :] = np.eye(6)
    L[6:, :] = np.eye(6)
    L[6, 5] = -1.0
    L[8, 3] = 1.0
    M = np.diag(np.concatenate([m2, m3]))
    F_ext = np.concatenate([np.zeros(6), f3])
    ar_ref = np.linalg.solve(L.T @ M @ L, L.T @ F_ext)
    a_ref = L @ ar_ref

    dt = 1.0e-4
    assert ops.analyze(1, dt) == 0
    u_step1 = np.concatenate([ops.nodeDisp(2), ops.nodeDisp(3)])
    a_step1_ops = 2.0 * u_step1 / (dt ** 2)
    rel_step1 = np.linalg.norm(a_step1_ops - a_ref) / np.linalg.norm(a_ref)
    assert rel_step1 < 1e-13, (
        f"Projector Cholesky nRet=6 step-1 accel mismatch vs numpy solve: {rel_step1:.3e}"
    )

    for _ in range(99):
        assert ops.analyze(1, dt) == 0
    u2 = ops.nodeDisp(2)
    u3 = ops.nodeDisp(3)
    assert abs(u3[0] - (u2[0] - 1.0 * u2[5])) < 1e-12
    assert abs(u3[1] - u2[1]) < 1e-12
    assert abs(u3[2] - (u2[2] + 1.0 * u2[3])) < 1e-12
    assert max(abs(x) for x in u2) > 0.0

    # Generalized momentum balance at step 100: L^T (M a) == L^T (F_ext - F_beam_node2)
    a_now = np.concatenate([ops.nodeAccel(2), ops.nodeAccel(3)])
    f_beam_n2 = np.array(ops.eleForce(1)[6:12])
    F_net = np.concatenate([-f_beam_n2, f3])
    lhs = L.T @ (M @ a_now)
    rhs = L.T @ F_net
    rel_mom = np.linalg.norm(lhs - rhs) / np.linalg.norm(rhs)
    assert rel_mom < 1e-11, (
        f"Projector Cholesky nRet=6 generalized momentum balance error: {rel_mom:.3e}"
    )


# ----------------------------------------------------------------------------
# 8. LadrunoCPML threading (1T vs 4T bit-identity) + lazy tangent materialization
# ----------------------------------------------------------------------------
def _run_cpml_wave(integrator_kind, nthreads, nsteps=30):
    ops.wipe()
    ops.ladrunoThreads(nthreads)
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    E, nu, rho = 1.0e8, 0.25, 2000.0
    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)

    nx, ny, nz = 4, 4, 4
    L = 1.0

    def nid(i, j, k):
        return 1 + i + (nx + 1) * (j + (ny + 1) * k)

    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                ops.node(nid(i, j, k), i * L, j * L, k * L)
                if i == 0 or i == nx or j == 0 or j == ny or k == 0:
                    ops.fix(nid(i, j, k), 1, 1, 1)

    etag = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                ns = (
                    nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                    nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1),
                )
                in_core = (1 <= i <= nx - 2) and (1 <= j <= ny - 2) and (k >= 1)
                if in_core:
                    ops.element("LadrunoBrick", etag, *ns, 1, "-formulation", "uri", "-lumped")
                else:
                    ops.element(
                        "LadrunoCPML", etag, *ns,
                        "-E", E, "-nu", nu, "-rho", rho,
                        "-box", 1.0, float(nx - 1), 1.0, float(ny - 1), 1.0, float(nz),
                        1.0, 1.0, 1.0,
                        "-m", 2.0, "-R0", 1e-4, "-kappaMax", 1.0, "-alphaMax", 0.0,
                        "-lumped",
                    )
                etag += 1

    ops.timeSeries("Trig", 1, 0.0, 0.02, 0.01)
    ops.pattern("Plain", 1, 1)
    center = nid(nx // 2, ny // 2, nz)
    ops.load(center, 1.0e4, 0.0, 0.0)

    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("Diagonal" if integrator_kind != "implicit" else "BandGeneral")
    ops.test("NormDispIncr", 1e-10, 10)
    if integrator_kind == "cdl":
        ops.algorithm("Linear")
        ops.integrator("CentralDifferenceLadruno")
        dt = 5.0e-4
    elif integrator_kind == "eb":
        ops.algorithm("Linear")
        ops.integrator("ExplicitBathe", 0.54)
        dt = 5.0e-4
    else:
        ops.algorithm("Newton")
        ops.integrator("Newmark", 0.5, 0.25)
        dt = 1.0e-3
    ops.analysis("Transient")

    traj = []
    for _ in range(nsteps):
        assert ops.analyze(1, dt) == 0
        traj.append(ops.nodeDisp(center, 1))
    ops.ladrunoThreads(1)
    return traj


@pytest.mark.parametrize("kind", ["cdl", "eb", "implicit"])
def test_cpml_threaded_bit_identity(kind):
    t1 = _run_cpml_wave(kind, nthreads=1)
    t4 = _run_cpml_wave(kind, nthreads=4)
    _assert_bit_identical(t1, t4, f"CPML 1T vs 4T ({kind})")
    assert max(abs(u) for u in t1) > 0.0


# ----------------------------------------------------------------------------
# 9. Stress-dependent NDMaterial exclusion from bvCdCache + Brick setParameter
# ----------------------------------------------------------------------------
def test_bv_cd_cache_excludes_stress_dependent_hs():
    """LadrunoHS getInitialTangent() varies with committed stress & stage; bvCdCache
    must NOT cache c_d for it (0-ULP bit-identical to setParameter-disabled oracle)."""
    def run_brick_hs(disable_cache, stage_switch=False):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        ops.nDMaterial(
            "LadrunoHS", 1, 35000.0, 35000.0, 105000.0, 0.5, 10.0, 38.0, 8.0,
            "-nu", 0.2, "-pref", 100.0, "-Rf", 0.9, "-rho", 2.0,
            "-stage", 0 if stage_switch else 1,
        )
        coords = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                  (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
        for i, (x, y, z) in enumerate(coords, 1):
            ops.node(i, float(x), float(y), float(z))
            if z == 0:
                ops.fix(i, 1, 1, 1)
        ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1,
                    "-formulation", "std", "-lumped", "-bulkViscosity", 0.5, 1.2)
        if disable_cache:
            ops.parameter(1, "element", 1, "material", 1, "E50_ref")
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        for n in (5, 6, 7, 8):
            ops.load(n, 0.0, 0.0, -500.0)
        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("Diagonal")
        ops.test("NormDispIncr", 1e-12, 1)
        ops.algorithm("Linear")
        ops.integrator("CentralDifferenceLadruno")
        ops.analysis("Transient")
        out = []
        n1 = 20 if stage_switch else 50
        for _ in range(n1):
            assert ops.analyze(1, 5e-4) == 0
            out.append(ops.nodeDisp(5, 3))
        if stage_switch:
            ops.updateMaterialStage("-material", 1, "-stage", 1)
            for _ in range(30):
                assert ops.analyze(1, 5e-4) == 0
                out.append(ops.nodeDisp(5, 3))
        return out

    _assert_bit_identical(
        run_brick_hs(False, False), run_brick_hs(True, False),
        "LadrunoBrick + LadrunoHS (stage 1) bvCdCache exclusion",
    )
    _assert_bit_identical(
        run_brick_hs(False, True), run_brick_hs(True, True),
        "LadrunoBrick + LadrunoHS (stage 0->1) bvCdCache exclusion",
    )

    def run_quad_hs(disable_cache):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        ops.nDMaterial(
            "LadrunoHS", 1, 35000.0, 35000.0, 105000.0, 0.5, 10.0, 38.0, 8.0,
            "-nu", 0.2, "-pref", 100.0, "-Rf", 0.9, "-rho", 2.0, "-stage", 1,
        )
        for i, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], 1):
            ops.node(i, float(x), float(y))
            if y == 0:
                ops.fix(i, 1, 1)
        ops.element("LadrunoQuad", 1, 1, 2, 3, 4, 1, "-type", "PlaneStrain", "-thick", 1.0,
                    "-formulation", "std", "-bulkViscosity", 0.5, 1.2)
        if disable_cache:
            ops.parameter(1, "element", 1, "material", 1, "E50_ref")
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(3, 0.0, -500.0)
        ops.load(4, 0.0, -500.0)
        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("Diagonal")
        ops.test("NormDispIncr", 1e-12, 1)
        ops.algorithm("Linear")
        ops.integrator("CentralDifferenceLadruno")
        ops.analysis("Transient")
        out = []
        for _ in range(50):
            assert ops.analyze(1, 5e-4) == 0
            out.append(ops.nodeDisp(3, 2))
        return out

    _assert_bit_identical(
        run_quad_hs(False), run_quad_hs(True),
        "LadrunoQuad + LadrunoHS bvCdCache exclusion",
    )


def test_brick_bv_cd_cache_set_parameter_invalidation():
    """Verify LadrunoBrick ElasticIsotropic bvCdCache is 0-ULP identical to uncached
    before and after an element-level updateParameter(E)."""
    def run(attach_before_steps):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        ops.nDMaterial("ElasticIsotropic", 1, 1.0e7, 0.25, 2.0)
        coords = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                  (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
        for i, (x, y, z) in enumerate(coords, 1):
            ops.node(i, float(x), float(y), float(z))
            if z == 0:
                ops.fix(i, 1, 1, 1)
        ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1,
                    "-formulation", "std", "-lumped", "-bulkViscosity", 0.5, 1.2)
        if attach_before_steps:
            ops.parameter(1, "element", 1, "E")
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        for n in (5, 6, 7, 8):
            ops.load(n, 0.0, 0.0, -500.0)
        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("Diagonal")
        ops.test("NormDispIncr", 1e-12, 1)
        ops.algorithm("Linear")
        ops.integrator("CentralDifferenceLadruno")
        ops.analysis("Transient")
        out = []
        for _ in range(20):
            assert ops.analyze(1, 5e-5) == 0
            out.append(ops.nodeDisp(5, 3))
        if not attach_before_steps:
            # Attach mid-run after bvCdCache has been populated for 20 steps!
            ops.parameter(1, "element", 1, "E")
        ops.updateParameter(1, 4.0e7)
        for _ in range(30):
            assert ops.analyze(1, 5e-5) == 0
            out.append(ops.nodeDisp(5, 3))
        return out

    u_pre = run(True)
    u_mid = run(False)
    _assert_bit_identical(u_pre, u_mid, "LadrunoBrick bvCdCache mid-run setParameter")
    assert max(abs(x) for x in u_pre) > 0.0


def test_quad_bv_cd_cache_set_parameter_invalidation():
    """Verify LadrunoQuad ElasticIsotropic bvCdCache is 0-ULP identical to uncached
    before and after a mid-run element-level updateParameter(E)."""
    def run(attach_before_steps):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        ops.nDMaterial("ElasticIsotropic", 1, 1.0e7, 0.25, 2.0)
        for i, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], 1):
            ops.node(i, float(x), float(y))
            if y == 0:
                ops.fix(i, 1, 1)
        ops.element("LadrunoQuad", 1, 1, 2, 3, 4, 1, "-type", "PlaneStrain", "-thick", 1.0,
                    "-formulation", "std", "-bulkViscosity", 0.5, 1.2)
        if attach_before_steps:
            ops.parameter(1, "element", 1, "E")
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(3, 0.0, -500.0)
        ops.load(4, 0.0, -500.0)
        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("Diagonal")
        ops.test("NormDispIncr", 1e-12, 1)
        ops.algorithm("Linear")
        ops.integrator("CentralDifferenceLadruno")
        ops.analysis("Transient")
        out = []
        for _ in range(20):
            assert ops.analyze(1, 5e-5) == 0
            out.append(ops.nodeDisp(3, 2))
        if not attach_before_steps:
            ops.parameter(1, "element", 1, "E")
        ops.updateParameter(1, 4.0e7)
        for _ in range(30):
            assert ops.analyze(1, 5e-5) == 0
            out.append(ops.nodeDisp(3, 2))
        return out

    u_pre = run(True)
    u_mid = run(False)
    _assert_bit_identical(u_pre, u_mid, "LadrunoQuad bvCdCache mid-run setParameter")
    assert max(abs(x) for x in u_pre) > 0.0


@pytest.mark.parametrize("elem_kind", ["brick", "quad"])
def test_bv_cd_cache_recvself_preserves_parameter_invalidation(elem_kind):
    """Finding 1 regression gate: in-place database restore (Domain::recvSelf) on an
    element with an attached Parameter must NOT reset bvCdCacheDisabled to false.
    Scenario: -bulkViscosity on, parameter on E, database File + save, analyze,
    restore, analyze (would populate stale cache if reset), updateParameter(E),
    analyze — must be 0-ULP bit-identical to the cache-disabled oracle (re-attaching
    a parameter immediately after restore)."""
    def run(reattach_after_restore, do_update=True):
        with tempfile.TemporaryDirectory() as tmp:
            db = os.path.join(tmp, f"bv_chk_{elem_kind}")
            ops.wipe()
            if elem_kind == "brick":
                ops.model("basic", "-ndm", 3, "-ndf", 3)
                ops.nDMaterial("ElasticIsotropic", 1, 1.0e7, 0.25, 2.0)
                coords = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                          (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
                for i, (x, y, z) in enumerate(coords, 1):
                    ops.node(i, float(x), float(y), float(z))
                    if z == 0:
                        ops.fix(i, 1, 1, 1)
                ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1,
                            "-formulation", "std", "-lumped", "-bulkViscosity", 0.5, 1.2)
                probe_node, probe_dof = 5, 3
                loads = [(n, (0.0, 0.0, -500.0)) for n in (5, 6, 7, 8)]
            else:
                ops.model("basic", "-ndm", 2, "-ndf", 2)
                ops.nDMaterial("ElasticIsotropic", 1, 1.0e7, 0.25, 2.0)
                for i, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], 1):
                    ops.node(i, float(x), float(y))
                    if y == 0:
                        ops.fix(i, 1, 1)
                ops.element("LadrunoQuad", 1, 1, 2, 3, 4, 1, "-type", "PlaneStrain",
                            "-thick", 1.0, "-formulation", "std", "-bulkViscosity", 0.5, 1.2)
                probe_node, probe_dof = 3, 2
                loads = [(3, (0.0, -500.0)), (4, (0.0, -500.0))]

            ops.parameter(1, "element", 1, "E")
            ops.timeSeries("Constant", 1)
            ops.pattern("Plain", 1, 1)
            for n, fv in loads:
                ops.load(n, *fv)
            ops.constraints("Plain")
            ops.numberer("Plain")
            ops.system("Diagonal")
            ops.test("NormDispIncr", 1e-12, 1)
            ops.algorithm("Linear")
            ops.integrator("ExplicitBathe", 0.54)
            ops.analysis("Transient")
            ops.database("File", db)

            for _ in range(10):
                assert ops.analyze(1, 5e-5) == 0
            ops.save(1)
            for _ in range(5):
                assert ops.analyze(1, 5e-5) == 0
            ops.restore(1)
            if reattach_after_restore:
                ops.parameter(2, "element", 1, "E")
            # Run steps AFTER restore BEFORE updateParameter: if recvSelf cleared
            # bvCdCacheDisabled, these steps populate bvCdCache with stale E=1e7!
            for _ in range(5):
                assert ops.analyze(1, 5e-5) == 0
            if do_update:
                ops.updateParameter(1, 4.0e7)
            out = []
            for _ in range(30):
                assert ops.analyze(1, 5e-5) == 0
                out.append(ops.nodeDisp(probe_node, probe_dof))
            ops.wipe()
            return out

    u_restored = run(reattach_after_restore=False, do_update=True)
    u_disabled_oracle = run(reattach_after_restore=True, do_update=True)
    u_no_update = run(reattach_after_restore=True, do_update=False)

    _assert_bit_identical(
        u_restored, u_disabled_oracle,
        f"{elem_kind} recvSelf + attached Parameter bvCdCache invalidation",
    )
    peak = max(abs(x) for x in u_disabled_oracle)
    assert peak > 0.0
    assert max(abs(a - b) for a, b in zip(u_restored, u_no_update)) / peak > 1e-2, (
        "non-vacuity: updateParameter(E) after restore must change trajectory"
    )


# ----------------------------------------------------------------------------
# 10. Lumped pure-alphaM O(n) fast paths: 0-ULP oracle vs Element::getRayleighDampingForces
# ----------------------------------------------------------------------------
@pytest.mark.parametrize("elem_kind", ["brick", "brick20", "quad", "cpml"])
def test_lumped_alphaM_fast_path_0ulp_oracle(elem_kind):
    """Compare pure alphaM fast path (betaK=0) against Element::getRayleighDampingForces
    dense fallback triggered via betaK=1e-300 (contributes 0.0 in FP64 -> 0-ULP oracle)."""
    def run(beta_k):
        ops.wipe()
        if elem_kind == "quad":
            ops.model("basic", "-ndm", 2, "-ndf", 2)
            ops.nDMaterial("ElasticIsotropic", 1, 1.0e7, 0.25, 2.0)
            for i, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], 1):
                ops.node(i, float(x), float(y))
                if y == 0:
                    ops.fix(i, 1, 1)
            ops.element("LadrunoQuad", 1, 1, 2, 3, 4, 1, "-type", "PlaneStrain", "-thick", 1.0)
            load_node, dof, dt = 3, 1, 1e-4
            load_vec = (100.0, 0.0)
        elif elem_kind == "brick20":
            ops.model("basic", "-ndm", 3, "-ndf", 3)
            ops.nDMaterial("ElasticIsotropic", 1, 1.0e7, 0.25, 2.0)
            corners = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                       (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
            mids = [(0.5, 0, 0), (1, 0.5, 0), (0.5, 1, 0), (0, 0.5, 0),
                    (0.5, 0, 1), (1, 0.5, 1), (0.5, 1, 1), (0, 0.5, 1),
                    (0, 0, 0.5), (1, 0, 0.5), (1, 1, 0.5), (0, 1, 0.5)]
            for idx, (x, y, z) in enumerate(corners + mids, start=1):
                ops.node(idx, float(x), float(y), float(z))
                if z == 0:
                    ops.fix(idx, 1, 1, 1)
            ops.element("LadrunoBrick20", 1, *range(1, 21), 1, "-formulation", "std", "-lumped")
            load_node, dof, dt = 7, 1, 5e-5
            load_vec = (100.0, 0.0, 0.0)
        else:
            ops.model("basic", "-ndm", 3, "-ndf", 3)
            coords = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                      (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
            for i, (x, y, z) in enumerate(coords, 1):
                ops.node(i, float(x), float(y), float(z))
                if z == 0:
                    ops.fix(i, 1, 1, 1)
            if elem_kind == "brick":
                ops.nDMaterial("ElasticIsotropic", 1, 1.0e7, 0.25, 2.0)
                ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1,
                            "-formulation", "uri", "-lumped")
                dt = 1e-4
            else:
                ops.element(
                    "LadrunoCPML", 1, 1, 2, 3, 4, 5, 6, 7, 8,
                    "-E", 1.0e8, "-nu", 0.25, "-rho", 2000.0,
                    "-box", 0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 1.0, 1.0, 1.0,
                    "-m", 2.0, "-R0", 1e-4, "-kappaMax", 1.0, "-alphaMax", 0.0,
                    "-lumped",
                )
                dt = 5e-4
            load_node, dof = 7, 1
            load_vec = (1000.0, 0.0, 0.0)

        ops.rayleigh(25.0, beta_k, 0.0, 0.0)
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(load_node, *load_vec)
        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("Diagonal")
        ops.test("NormDispIncr", 1e-12, 1)
        ops.algorithm("Linear")
        ops.integrator("CentralDifferenceLadruno")
        ops.analysis("Transient")
        out = []
        for _ in range(40):
            assert ops.analyze(1, dt) == 0
            out.append(ops.nodeDisp(load_node, dof))
        return out

    u_fast = run(0.0)
    u_oracle = run(1.0e-300)
    _assert_bit_identical(u_fast, u_oracle, f"lumped alphaM fast path ({elem_kind})")
    assert max(abs(x) for x in u_fast) > 0.0


# ----------------------------------------------------------------------------
# 11. LadrunoCPML lazy tangents revertToLastCommit + Biot u-p Rayleigh 8T + restart
# ----------------------------------------------------------------------------
def test_cpml_lazy_tangents_revert_to_last_commit():
    """Forced failed step at a different dt triggers revertToLastCommit; retrying at
    the original dt must rematerialize K_eff_/C_eff_ bit-identically to a clean run."""
    def run(inject_fail):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        coords = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                  (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
        for i, (x, y, z) in enumerate(coords, 1):
            ops.node(i, float(x), float(y), float(z))
            if z == 0:
                ops.fix(i, 1, 1, 1)
        ops.element(
            "LadrunoCPML", 1, 1, 2, 3, 4, 5, 6, 7, 8,
            "-E", 1.0e8, "-nu", 0.25, "-rho", 2000.0,
            "-box", 0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 1.0, 1.0, 1.0,
            "-m", 2.0, "-R0", 1e-4, "-kappaMax", 1.0, "-alphaMax", 0.0,
        )
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(5, 1.0e4, 0.0, 0.0)
        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("BandGeneral")
        ops.test("NormDispIncr", 1e-12, 10)
        ops.algorithm("Newton")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient")
        out = []
        for _ in range(5):
            assert ops.analyze(1, 1e-3) == 0
            out.append(ops.nodeDisp(5, 1))
        if inject_fail:
            ops.test("NormDispIncr", 1e-30, 1)
            assert ops.analyze(1, 5e-3) != 0
            ops.test("NormDispIncr", 1e-12, 10)
        for _ in range(5):
            assert ops.analyze(1, 1e-3) == 0
            out.append(ops.nodeDisp(5, 1))
        return out

    _assert_bit_identical(run(True), run(False), "CPML lazy tangents revertToLastCommit")


def test_cpml_biot_up_rayleigh_threaded_bit_identity():
    """Biot u-p LadrunoCPML with Rayleigh damping must be 0-ULP bit-identical 1T vs 8T."""
    def run(nt):
        ops.wipe()
        ops.ladrunoThreads(nt)
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        nx, ny = 4, 3

        def nid(i, j):
            return 1 + i + (nx + 1) * j

        for j in range(ny + 1):
            for i in range(nx + 1):
                ops.node(nid(i, j), float(i), float(j))
                if i == 0 or i == nx or j == 0:
                    ops.fix(nid(i, j), 1, 1, 1)
        etag = 1
        for j in range(ny):
            for i in range(nx):
                ns = (nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1))
                ops.element(
                    "LadrunoCPML", etag, *ns,
                    "-E", 1.0e7, "-nu", 0.25, "-rho", 2.0,
                    "-box", 1.0, 3.0, 1.0, 3.0, "-Lx", 1.0, "-Ly", 1.0,
                    "-m", 2.0, "-R0", 1e-4, "-kappaMax", 1.0, "-alphaMax", 0.0,
                    "-up", "-Kf", 2.0e6, "-poro", 0.4, "-rhoF", 1.0,
                    "-alphaB", 1.0, "-Ks", 1.0e10,
                    "-permH", 1e-3, 1e-3, "-gammaW", 9.81,
                )
                etag += 1
        ops.rayleigh(0.5, 1e-4, 0.0, 0.0)
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(nid(2, ny), 100.0, -50.0, 0.0)
        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("BandGeneral")
        ops.test("NormDispIncr", 1e-10, 10)
        ops.algorithm("Newton")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient")
        out = []
        for _ in range(15):
            assert ops.analyze(1, 1e-3) == 0
            out.extend(ops.nodeDisp(nid(2, ny)))
        ops.ladrunoThreads(1)
        return out

    t1 = run(1)
    t8 = run(8)
    _assert_bit_identical(t1, t8, "CPML Biot u-p + Rayleigh 1T vs 8T")
    assert max(abs(x) for x in t1) > 0.0


def test_explicit_bathe_and_cpml_database_save_restore():
    """Verify ExplicitBathe::sendSelf/recvSelf cache invalidation and LadrunoCPMLKernel
    packState/unpackState lazy tangent materialization across database save/restore."""
    with tempfile.TemporaryDirectory() as tmp:
        db = os.path.join(tmp, "wp260b_chk")

        def build():
            ops.wipe()
            ops.model("basic", "-ndm", 3, "-ndf", 3)
            coords = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                      (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
            for i, (x, y, z) in enumerate(coords, 1):
                ops.node(i, float(x), float(y), float(z))
                if z == 0:
                    ops.fix(i, 1, 1, 1)
            ops.element(
                "LadrunoCPML", 1, 1, 2, 3, 4, 5, 6, 7, 8,
                "-E", 1.0e8, "-nu", 0.25, "-rho", 2000.0,
                "-box", 0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 1.0, 1.0, 1.0,
                "-m", 2.0, "-R0", 1e-4, "-kappaMax", 1.0, "-alphaMax", 0.0,
                "-lumped",
            )
            ops.timeSeries("Constant", 1)
            ops.pattern("Plain", 1, 1)
            ops.load(5, 1.0e4, 0.0, 0.0)
            ops.constraints("Plain")
            ops.numberer("Plain")
            ops.system("Diagonal")
            ops.test("NormDispIncr", 1e-12, 1)
            ops.algorithm("Linear")
            ops.integrator("ExplicitBathe", 0.54)
            ops.analysis("Transient")
            ops.database("File", db)

        build()
        for _ in range(10):
            assert ops.analyze(1, 5e-4) == 0
        ops.save(1)
        u_cont = []
        for _ in range(10):
            assert ops.analyze(1, 5e-4) == 0
            u_cont.append(ops.nodeDisp(5, 1))

        build()
        ops.restore(1)
        u_rest = []
        for _ in range(10):
            assert ops.analyze(1, 5e-4) == 0
            u_rest.append(ops.nodeDisp(5, 1))
        ops.wipe()  # release Windows file handles before TemporaryDirectory cleanup

    _assert_bit_identical(u_cont, u_rest, "ExplicitBathe + CPML database save/restore")
    assert max(abs(x) for x in u_cont) > 0.0
