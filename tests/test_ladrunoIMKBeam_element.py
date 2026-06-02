"""LadrunoIMKBeam (Ladruno concentrated-plasticity IMK beam, classTag 33003) —
Zone-A element battery.

The element is a 3D beam-column with UNCOUPLED moment-rotation rotational hinges
at the ends (elastic interior in series with an IMK uniaxial law per bending
axis). v1 correctness anchors:

  * **elastic-degenerate**: with NO hinge materials the element is a pure elastic
    beam and must reproduce ``elasticBeamColumn`` to ~1e-9 (same A,E,G,Jx,Iy,Iz,
    same CrdTransf) on a 6-dof mixed load. This proves the basic-system assembly
    and the basic->global transform.

  * **single-hinge series**: a strong-axis hinge at end j, driven under rotation
    control, must reproduce an INDEPENDENT python model of "elastic rotational
    spring a = 4EIz/L in series with the IMK backbone M(theta)". This proves the
    exact series F = F_elastic + F_hinge and the 2x2 internal Newton (no
    n-factor). The check is non-tautological: python solves a(theta - thj) =
    M(thj) for the hinge rotation thj and compares BOTH thj and the moment to the
    element's reported hingeRotation / hingeMoment.

v2 adds the ADR-§6 real-world correctness anchors:

  * **diaphragm axial immunity** (THE headline guarantee): a ``rigidDiaphragm``
    injects a large spurious axial force into the beam; the moment-rotation
    response is bit-identical to the no-axial run, the element converges, and an
    equivalent fiber-hinge ``forceBeamColumn``+``HingeRadau`` does NOT (the fake
    axial drives the fibers to fracture -> state determination fails).

  * **column-face offset**: with a ``geomTransf -jntOffset`` rigid offset the
    hinge moment (clear-span end) and the joint moment (reaction) differ by
    exactly ``V*offset`` (closed form).

  * **Lignos-Krawinkler cyclic**: a steel ``Bilin`` (modified-IMK) end hinge
    under a standard cyclic protocol -- the elastic-interior series identity
    ``a*(theta_node - theta_h) == M`` holds exactly, and the reported hinge
    moment reproduces the bare uniaxial material driven by the same
    hinge-rotation history (cyclic deterioration faithfully delegated), with real
    hysteretic dissipation.

Plan: Ladruno_implementation/14_ladruno_imk_beam.md.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# member elastic properties (consistent SI-ish units)
_A = 1.0e-2
_E = 2.0e11
_G = 8.0e10
_Jx = 1.0e-5
_Iy = 1.0e-5
_Iz = 2.0e-5
_L = 3.0


def _new_model():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, _L, 0.0, 0.0)
    # local x = global x; vecxz = global z -> local y = global y
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)


# ---------------------------------------------------------------------------
# Test 1: elastic-degenerate == elasticBeamColumn
# ---------------------------------------------------------------------------
def _solve_cantilever_tip(make_element):
    """Fix node 1, build one element via make_element(), apply a mixed 6-dof
    load at node 2, linear-solve, return the node-2 displacement (6)."""
    _new_model()
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    make_element()

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    # mixed load so the whole 6x6 basic stiffness is exercised
    ops.load(2, 1.0e4, 2.0e4, -1.5e4, 3.0e3, -4.0e3, 5.0e3)

    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-12, 10)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    return [ops.nodeDisp(2, dof) for dof in range(1, 7)]


def test_elastic_degenerate_matches_elasticBeamColumn():
    def make_imk():
        ops.element("LadrunoIMKBeam", 1, 1, 2,
                    _A, _E, _G, _Jx, _Iy, _Iz, 1)  # no -matZ/-matY -> elastic

    def make_elastic():
        ops.element("elasticBeamColumn", 1, 1, 2,
                    _A, _E, _G, _Jx, _Iy, _Iz, 1)

    d_imk = _solve_cantilever_tip(make_imk)
    d_ref = _solve_cantilever_tip(make_elastic)

    for a, b in zip(d_imk, d_ref):
        assert abs(a - b) <= 1.0e-9 * (abs(b) + 1.0e-12), (a, b)

    # sanity: the tip actually moved (load really applied)
    assert any(abs(x) > 1.0e-9 for x in d_ref)


# ---------------------------------------------------------------------------
# Test 2: single strong-axis hinge reproduces the elastic-series + IMK backbone
# ---------------------------------------------------------------------------
# Steel01 monotonic backbone (pure bilinear kinematic): M(theta).
_My = 1000.0     # yield moment
_k0 = 1.0e6      # initial rotational stiffness of the hinge
_bh = 0.02       # strain-hardening ratio
_theta_y = _My / _k0


def _imk_backbone(theta):
    if theta <= _theta_y:
        return _k0 * theta
    return _My + _bh * _k0 * (theta - _theta_y)


def _series_predict(theta_imposed, a):
    """Solve a*(theta - thj) = M(thj) for the hinge rotation thj (monotonic,
    thj >= 0), return (thj, moment)."""
    thj = 0.0
    for _ in range(100):
        if thj <= _theta_y:
            M = _k0 * thj
            dM = _k0
        else:
            M = _My + _bh * _k0 * (thj - _theta_y)
            dM = _bh * _k0
        g = a * (theta_imposed - thj) - M
        if abs(g) <= 1.0e-12 * (a * theta_imposed + 1.0):
            break
        # dg/dthj = -a - dM
        thj -= g / (-a - dM)
    return thj, _imk_backbone(thj)


def test_single_hinge_series_matches_independent_model():
    # Iz chosen so a = 4 E Iz / L is a round number is unnecessary; use _Iz.
    a = 4.0 * _E * _Iz / _L

    _new_model()
    ops.uniaxialMaterial("Steel01", 10, _My, _k0, _bh)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    # node 2: only the strong-axis rotation (dof 6) is free
    ops.fix(2, 1, 1, 1, 1, 1, 0)

    ops.element("LadrunoIMKBeam", 1, 1, 2,
                _A, _E, _G, _Jx, _Iy, _Iz, 1,
                "-hinge", "j", "-matZ", 10)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0)  # reference moment on dof 6

    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    dtheta = 1.0e-4
    ops.integrator("DisplacementControl", 2, 6, dtheta)
    ops.test("NormDispIncr", 1.0e-12, 50)
    ops.algorithm("Newton")
    ops.analysis("Static")

    n_steps = 45  # spans elastic (theta_y=1e-3) well into hardening (~4.5e-3)
    for step in range(1, n_steps + 1):
        assert ops.analyze(1) == 0, f"step {step} failed"

        theta_imposed = ops.nodeDisp(2, 6)
        hingeRot = ops.eleResponse(1, "hingeRotation")    # [z_i,z_j,y_i,y_j]
        hingeMom = ops.eleResponse(1, "hingeMoment")      # [Mz_i,Mz_j,My_i,My_j]
        thj_ele = hingeRot[1]
        M_ele = hingeMom[1]

        thj_pred, M_pred = _series_predict(theta_imposed, a)

        assert abs(thj_ele - thj_pred) <= 1.0e-7 * (abs(thj_pred) + 1.0e-9), (
            f"step {step}: hinge rot {thj_ele} vs {thj_pred}")
        assert abs(M_ele - M_pred) <= 1.0e-6 * (abs(M_pred) + 1.0e-6), (
            f"step {step}: hinge moment {M_ele} vs {M_pred}")

    # we must have crossed yield for the test to be meaningful
    assert ops.nodeDisp(2, 6) > _theta_y * 2.0


def test_no_hinge_when_material_omitted_is_pure_elastic():
    """A LadrunoIMKBeam with -hinge but no -matZ/-matY must still be elastic
    (the hinge exists only where a material is supplied)."""
    def make_imk_hinge_no_mat():
        ops.element("LadrunoIMKBeam", 1, 1, 2,
                    _A, _E, _G, _Jx, _Iy, _Iz, 1, "-hinge", "both")

    def make_elastic():
        ops.element("elasticBeamColumn", 1, 1, 2,
                    _A, _E, _G, _Jx, _Iy, _Iz, 1)

    d_imk = _solve_cantilever_tip(make_imk_hinge_no_mat)
    d_ref = _solve_cantilever_tip(make_elastic)
    for a, b in zip(d_imk, d_ref):
        assert abs(a - b) <= 1.0e-9 * (abs(b) + 1.0e-12), (a, b)


# ---------------------------------------------------------------------------
# Test 4: rigid-diaphragm spurious axial -- moment-rotation immunity, AND the
#         fiber-hinge forceBeamColumn fails to converge where this element does
# ---------------------------------------------------------------------------
# The motivating failure mode (ADR section 1): a rigid floor diaphragm forces the
# beam ends to a common in-plane motion, injecting a spurious axial force P. A
# fiber hinge couples P-M, contaminating the moment capacity and breaking the
# solver; the uncoupled moment-rotation hinge is immune.
#
# Setup -- beam along global X; the xy-plane diaphragm (perpDirn=3) ties the
# beam's axial dof (uX), uY and rZ of node 2 to a master node 3 that is pulled in
# X, so the diaphragm IS the axial-force source. Bending is about global Y
# (rotation rY, the WEAK axis -> -matY hinge), which the xy-diaphragm leaves free;
# vecxz = global Z aligns local y with global Y. node 2's rY carries the bending
# demand (moment load about Y at node 2).
_diaP = 1.0e6     # axial force injected through the diaphragm (>> fiber capacity)
_diaM = 2500.0    # bending moment demand about Y (yields the My=1000 hinge)


def _diaphragm_imk(P, n_steps=40):
    """Run the LadrunoIMKBeam diaphragm model; return per-step
    (hingeRotation_Yj, hingeMoment_Yj, axialN) or None if it failed."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)            # i-end, fixed
    ops.node(2, _L, 0.0, 0.0)             # j-end (hinge); bends about free rY
    ops.node(3, 0.0, 0.0, 0.0)            # diaphragm master, pulled in X
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)

    ops.uniaxialMaterial("Steel01", 10, 1000.0, 1.0e6, 0.03)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.fix(2, 0, 0, 1, 1, 0, 0)          # uX,uY,rZ via diaphragm; uZ,rX fixed; rY free
    ops.fix(3, 0, 1, 1, 1, 1, 1)          # master: only uX free (axial driver)
    ops.rigidDiaphragm(3, 3, 2)           # perpDirn=3 (xy plane): ties uX,uY,rZ of 2 to 3

    ops.element("LadrunoIMKBeam", 1, 1, 2, _A, _E, _G, _Jx, _Iy, _Iz, 1,
                "-hinge", "j", "-matY", 10)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0, 0, 0, 0, _diaM, 0)     # bending moment about Y at j
    if P:
        ops.load(3, P, 0, 0, 0, 0, 0)     # axial pull at the diaphragm master
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 1.0 / n_steps)
    ops.test("NormDispIncr", 1.0e-10, 60, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")

    hist = []
    for _ in range(n_steps):
        if ops.analyze(1) != 0:
            return None
        thh = ops.eleResponse(1, "hingeRotation")[3]   # [Zi,Zj,Yi,Yj] -> Yj
        M = ops.eleResponse(1, "hingeMoment")[3]
        N = ops.eleResponse(1, "basicForce")[0]
        hist.append((thh, M, N))
    return hist


def _diaphragm_fiber(P, n_steps=40):
    """Same model/demand but an equivalent fiber-hinge forceBeamColumn+HingeRadau.
    Returns the step index at which it failed, or n_steps if it completed."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, _L, 0.0, 0.0)
    ops.node(3, 0.0, 0.0, 0.0)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.fix(2, 0, 0, 1, 1, 0, 0)
    ops.fix(3, 0, 1, 1, 1, 1, 1)
    ops.rigidDiaphragm(3, 3, 2)

    # member-scale steel fiber section with a fracture strain: the spurious axial
    # drives the fibers past it (capacity loss) -> genuine convergence failure.
    Fy = 8.0e6                                  # N_cap = A*Fy = 8e4 << _diaP
    ops.uniaxialMaterial("Steel01", 21, Fy, _E, 0.01)
    ops.uniaxialMaterial("MinMax", 20, 21, "-min", -0.02, "-max", 0.02)
    ops.section("Fiber", 99, "-GJ", _G * _Jx)
    ops.patch("quad", 20, 8, 8, -0.05, -0.05, 0.05, -0.05, 0.05, 0.05, -0.05, 0.05)
    ops.section("Elastic", 98, _E, _A, _Iz, _Iy, _G, _Jx)
    ops.beamIntegration("HingeRadau", 1, 99, 0.1, 99, 0.1, 98)
    ops.element("forceBeamColumn", 1, 1, 2, 1, 1)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0, 0, 0, 0, _diaM, 0)
    if P:
        ops.load(3, P, 0, 0, 0, 0, 0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 1.0 / n_steps)
    ops.test("NormDispIncr", 1.0e-10, 60, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    for step in range(n_steps):
        if ops.analyze(1) != 0:
            return step
    return n_steps


def test_diaphragm_axial_immunity_and_fiber_divergence():
    n = 40
    with_axial = _diaphragm_imk(_diaP, n)
    no_axial = _diaphragm_imk(0.0, n)

    # the LadrunoIMKBeam converges with and without the spurious axial
    assert with_axial is not None and len(with_axial) == n
    assert no_axial is not None and len(no_axial) == n

    # the diaphragm really injected a large axial force (and none without it)
    assert abs(with_axial[-1][2] - _diaP) <= 1.0e-3 * _diaP
    assert abs(no_axial[-1][2]) <= 1.0e-3 * _diaP

    # the hinge actually yielded (so the test exercises plasticity, not elastic)
    assert max(abs(M) for _, M, _ in no_axial) > 1.2 * 1000.0

    # THE GUARANTEE: moment-rotation is unchanged by the spurious axial force
    for (thA, MA, _), (thB, MB, _) in zip(with_axial, no_axial):
        assert abs(thA - thB) <= 1.0e-9 * (abs(thB) + 1.0e-9), (thA, thB)
        assert abs(MA - MB) <= 1.0e-7 * (abs(MB) + 1.0e-3), (MA, MB)

    # the equivalent fiber hinge converges WITHOUT the axial but FAILS with it
    assert _diaphragm_fiber(0.0, n) == n, "fiber model should run with no axial"
    fiber_failed_at = _diaphragm_fiber(_diaP, n)
    assert fiber_failed_at < n, (
        "fiber-hinge forceBeamColumn should fail to converge under the spurious "
        f"axial (failed at {fiber_failed_at} of {n})")


# ---------------------------------------------------------------------------
# Test 5: column-face offset -- hinge moment vs joint moment differ by V*offset
# ---------------------------------------------------------------------------
def _offset_run(off, V=1.0e3):
    """Cantilever along X with a member-axis rigid offset at end i (via the
    geomTransf -jntOffset). Returns (V, hingeMoment_Zi, jointReaction_Mz)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, _L, 0.0, 0.0)
    # rigid offset at the i-end along +X (member axis); offsets are GLOBAL coords
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0,
                   "-jntOffset", off, 0.0, 0.0, 0.0, 0.0, 0.0)
    ops.uniaxialMaterial("Steel01", 10, 1.0e9, 1.0e8, 0.01)   # high yield -> elastic
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.element("LadrunoIMKBeam", 1, 1, 2, _A, _E, _G, _Jx, _Iy, _Iz, 1,
                "-hinge", "i", "-matZ", 10)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, V, 0.0, 0.0, 0.0, 0.0)   # transverse Y -> strong-axis bending
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-12, 20, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    ops.reactions()
    Mz_i = ops.eleResponse(1, "hingeMoment")[0]   # Mz at the hinge (clear-span i-end)
    Rz1 = ops.nodeReaction(1, 6)                   # Mz reaction at the joint node
    return V, Mz_i, Rz1


def test_column_face_offset_moment_vs_joint_moment():
    """The rigid joint offset puts the hinge at the column face; the moment there
    is the joint moment minus V*offset (the shear transported across the rigid
    link). Closed form across several offsets."""
    # zero offset: hinge moment == joint moment == V*L (a cantilever)
    V, Mz0, Rz0 = _offset_run(0.0)
    assert abs(abs(Mz0) - V * _L) <= 1.0e-7 * V * _L
    assert abs(abs(Rz0) - V * _L) <= 1.0e-7 * V * _L

    for off in (0.3, 0.6, 0.9):
        V, Mz_i, Rz1 = _offset_run(off)
        # joint reaction is the full cantilever moment V*L (offset-independent)
        assert abs(abs(Rz1) - V * _L) <= 1.0e-7 * V * _L, (off, Rz1)
        # hinge moment is on the clear span: V*(L - off)
        assert abs(abs(Mz_i) - V * (_L - off)) <= 1.0e-7 * V * _L, (off, Mz_i)
        # the closed-form difference IS V*offset (shear x rigid-link length)
        assert abs((abs(Rz1) - abs(Mz_i)) - V * off) <= 1.0e-6 * V * (_L + off), off


# ---------------------------------------------------------------------------
# Test 6: Lignos-Krawinkler cyclic -- steel Bilin (modified-IMK) end hinge
# ---------------------------------------------------------------------------
_cyc_My, _cyc_Ke, _cyc_b = 1000.0, 1.0e6, 0.03


def _build_bilin(tag):
    """Lignos-Krawinkler modified-IMK bilinear steel-beam hinge (the classic
    deteriorating IMK moment-rotation law)."""
    ops.uniaxialMaterial("Bilin", tag,
                         _cyc_Ke, _cyc_b, _cyc_b, _cyc_My, -_cyc_My,
                         1.0, 1.0, 1.0, 1.0,        # LamdaS, LamdaD, LamdaA, LamdaK
                         1.0, 1.0, 1.0, 1.0,        # Cs, Cd, Ca, Ck
                         0.03, 0.03,                # Thetap_pos, Thetap_neg
                         0.20, 0.20,                # Thetapc_pos, Thetapc_neg
                         0.4, 0.4,                  # KPos, KNeg (residual ratio)
                         0.40, 0.40,                # Thetau_pos, Thetau_neg
                         1.0, 1.0)                  # PDPlus, PDNeg


def _cyclic_targets(amps=(0.004, 0.008, 0.016, 0.028), ncyc=2):
    targets = []
    for amp in amps:
        for _ in range(ncyc):
            targets += [amp, -amp]
    targets += [0.0]
    return targets


def test_lignos_krawinkler_cyclic_matches_bare_material():
    dth = 4.0e-4
    # --- drive the element through a standard symmetric cyclic protocol ---
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, _L, 0.0, 0.0)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    _build_bilin(10)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.fix(2, 1, 1, 1, 1, 1, 0)        # only strong-axis rotation free (like test 2)
    ops.element("LadrunoIMKBeam", 1, 1, 2, _A, _E, _G, _Jx, _Iy, _Iz, 1,
                "-hinge", "j", "-matZ", 10)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0, 0, 0, 0, 0, 1.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-12, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("DisplacementControl", 2, 6, dth)
    ops.analysis("Static")

    hist = []           # (theta_node, theta_hinge, moment)
    cur = 0.0
    for tgt in _cyclic_targets():
        step = dth if tgt >= cur else -dth
        nstep = max(1, int(round(abs(tgt - cur) / dth)))
        ops.integrator("DisplacementControl", 2, 6, step)
        for _ in range(nstep):
            assert ops.analyze(1) == 0, "cyclic element analysis failed"
            th2 = ops.nodeDisp(2, 6)
            thh = ops.eleResponse(1, "hingeRotation")[1]
            M = ops.eleResponse(1, "hingeMoment")[1]
            hist.append((th2, thh, M))
        cur = tgt

    Mscale = max(abs(M) for _, _, M in hist)
    thh_max = max(abs(thh) for _, thh, _ in hist)
    assert Mscale > 1.2 * _cyc_My          # well past yield (deterioration active)
    assert thh_max > 5.0 * (_cyc_My / _cyc_Ke)

    # (i) elastic-interior series identity holds cyclically and exactly:
    #     a*(theta_node - theta_hinge) == hinge moment, a = 4 E Iz / L
    a = 4.0 * _E * _Iz / _L
    series_err = max(abs(a * (th2 - thh) - M) for th2, thh, M in hist)
    assert series_err <= 1.0e-6 * Mscale, series_err

    # (ii) the element faithfully delegates the FULL cyclic path (incl. cyclic
    #      deterioration) to the hinge material: driving the bare Bilin material
    #      with the element's reported hinge-rotation history reproduces the
    #      reported hinge moment at every step.
    ops.wipe()
    _build_bilin(99)
    ops.testUniaxialMaterial(99)
    for (_, thh, M) in hist:
        ops.setStrain(thh)              # auto-commits each step
        bare = ops.getStress()
        assert abs(M - bare) <= 1.0e-6 * (Mscale) + 1.0e-6, (M, bare)

    # (iii) real hysteretic dissipation occurred (closed integral of M d(theta_h))
    energy = 0.0
    for k in range(1, len(hist)):
        energy += 0.5 * (hist[k][2] + hist[k - 1][2]) * (hist[k][1] - hist[k - 1][1])
    assert energy > 0.1 * _cyc_My * thh_max     # substantial enclosed hysteresis area
