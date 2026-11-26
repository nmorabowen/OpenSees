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


# ----------------------------------------------------------------------------
# 4. ADR-67 P1.1 lifetime-invariant gate: consistent-SMS + mid-run removeElement
# ----------------------------------------------------------------------------
@pytest.mark.parametrize("integrator_kind", ["cdl_sms_cons", "eb_sms_cons"])
def test_sms_consistent_remove_element_lifetime_invariant(integrator_kind):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, 2.0e7, 0.25, 2.0)

    # 3-element chain: middle element is short (throttling) -> scaled by consistent SMS,
    # plus a redundant side truss that we remove mid-run to trigger domainChanged.
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

    ops.uniaxialMaterial("Elastic", 2, 1.0e5)
    ops.element("Truss", 99, nid[(0, 0, 0)], nid[(3, 0, 0)], 0.01, 2, "-rho", 1.0)

    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(nid[(3, 1, 1)], 1000.0, 0.0, 0.0)

    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    dt_target = 2.5e-4
    if integrator_kind == "cdl_sms_cons":
        ops.integrator("CentralDifferenceSMSConsistent", dt_target, "-maxAddedMass", 1.0)
    else:
        ops.integrator("ExplicitBathe", 0.54, "-sms", dt_target, "-consistent", "-maxAddedMass", 1.0)
    ops.analysis("Transient")

    dt_run = 0.8 * dt_target
    for _ in range(40):
        assert ops.analyze(1, dt_run) == 0
    u_before = ops.nodeDisp(nid[(3, 1, 1)], 1)
    assert math.isfinite(u_before) and abs(u_before) > 0.0

    # Mid-run topology change -> domainChanged -> ConsistentPCGScratch & blocks rebuild.
    ops.remove("element", 99)
    for _ in range(60):
        assert ops.analyze(1, dt_run) == 0
    u_after = ops.nodeDisp(nid[(3, 1, 1)], 1)
    assert math.isfinite(u_after) and abs(u_after) > 0.0


# ----------------------------------------------------------------------------
# 5. ADR-67 P1.2' massless-element short-circuit in CriticalTimeStep
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
    # Truss element pencil: K = (EA/L)[[1,-1],[-1,1]], M_lump = diag(m/2, m/2)
    # => lambda_max = 4*EA/(L*m) = 4*E/(rho*L^2) => dt_cr = 2/sqrt(lambda_max) = L/sqrt(E/rho)
    expected = L / math.sqrt(E / rho)
    assert abs(dt_cr - expected) / expected < 1e-10, f"dt_cr={dt_cr} vs expected={expected}"


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
# 7. LadrunoConstraintProjector pre-factored Cholesky L (nRet > 1 rigidLink)
# ----------------------------------------------------------------------------
def test_projector_cholesky_rigid_link_manifold():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 2.0, 0.0, 0.0)
    ops.node(3, 2.0, 1.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    # Give all 6 DOFs positive lumped mass (required by LadrunoProjection)
    ops.mass(2, 5.0, 5.0, 5.0, 0.5, 0.5, 0.5)
    ops.mass(3, 3.0, 3.0, 3.0, 0.3, 0.3, 0.3)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.element("elasticBeamColumn", 1, 1, 2, 1.0, 1.0e5, 4.0e4, 1.0, 0.1, 0.1, 1)
    # Rigid beam link couples 6 retained DOFs of node 2 -> exercises nRet=6 Cholesky L!
    ops.rigidLink("beam", 2, 3)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 10.0, 20.0, 5.0, 0.0, 0.0, 0.0)
    ops.constraints("LadrunoProjection")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno")
    ops.analysis("Transient")
    for _ in range(100):
        assert ops.analyze(1, 1.0e-4) == 0
    u2 = ops.nodeDisp(2)
    u3 = ops.nodeDisp(3)
    # Rigid link with r = (0, 1, 0): u3_x = u2_x - r_y * rz2, u3_y = u2_y, u3_z = u2_z + r_y * rx2
    assert abs(u3[0] - (u2[0] - 1.0 * u2[5])) < 1e-12
    assert abs(u3[1] - u2[1]) < 1e-12
    assert abs(u3[2] - (u2[2] + 1.0 * u2[3])) < 1e-12
    assert max(abs(x) for x in u2) > 0.0


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
