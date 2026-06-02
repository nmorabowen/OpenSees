"""LadrunoIMKBeam2d (Ladruno concentrated-plasticity IMK beam, 2D, classTag
33004) — Zone-A element battery.

Planar counterpart of LadrunoIMKBeam (3D, classTag 33003). A 2D beam-column with
UNCOUPLED moment-rotation rotational hinges at the ends (elastic interior in
series with an IMK uniaxial law for the single strong axis Mz). Both dimensions
share the per-axis hinge kernel LadrunoIMKHinge.h, so these checks also guard
that the refactor preserves the planar reduction.

Correctness anchors (mirroring the 3D battery, reduced to the planar signature):

  * **ndm-dispatch + elastic-degenerate**: the SAME command ``LadrunoIMKBeam``,
    in an ndm-2/ndf-3 model, builds the 2D class; with NO hinge material it is a
    pure elastic beam and reproduces ``elasticBeamColumn`` to ~1e-9.

  * **single-hinge series**: a strong-axis hinge at end j under rotation control
    reproduces an INDEPENDENT python model of "elastic rotational spring
    a = 4EIz/L in series with the IMK backbone M(theta)" -- both the hinge
    rotation and the moment. Proves the exact series and the 2x2 internal Newton.

  * **axial immunity**: a spurious axial force injected through an equalDOF tie
    leaves the moment-rotation response bit-identical to the no-axial run (the
    uncoupled hinge has no P-M interaction).

  * **column-face offset**: a geomTransf -jntOffset rigid offset makes the hinge
    moment (clear span) and the joint reaction differ by exactly V*offset.

  * **asymmetric end laws**: -matZi/-matZj place an independent law at each end.

  * **cyclic Bilin**: a modified-IMK steel hinge under a cyclic protocol -- the
    elastic-interior series identity holds exactly and the reported hinge moment
    reproduces the bare uniaxial material driven by the same rotation history
    (cyclic deterioration faithfully delegated), with real hysteresis.

Plan: Ladruno_implementation/14_ladruno_imk_beam.md.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# member elastic properties (consistent SI-ish units)
_A = 1.0e-2
_E = 2.0e11
_Iz = 2.0e-5
_L = 3.0


def _new_model():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, _L, 0.0)
    ops.geomTransf("Linear", 1)


# ---------------------------------------------------------------------------
# Test 1: ndm-dispatch + elastic-degenerate == elasticBeamColumn (2D)
# ---------------------------------------------------------------------------
def _solve_cantilever_tip(make_element):
    """Fix node 1, build one element via make_element(), apply a mixed 3-dof
    load at node 2, linear-solve, return the node-2 displacement (3)."""
    _new_model()
    ops.fix(1, 1, 1, 1)
    make_element()

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    # mixed load so the whole 3x3 basic stiffness is exercised
    ops.load(2, 1.0e4, 2.0e4, 5.0e3)

    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-12, 10)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    return [ops.nodeDisp(2, dof) for dof in range(1, 4)]


def test_dispatch_elastic_degenerate_matches_elasticBeamColumn():
    def make_imk():
        # SAME command name as 3D -> ndm-2 model dispatches to LadrunoIMKBeam2d
        ops.element("LadrunoIMKBeam", 1, 1, 2, _A, _E, _Iz, 1)  # no -matZ -> elastic

    def make_elastic():
        ops.element("elasticBeamColumn", 1, 1, 2, _A, _E, _Iz, 1)

    d_imk = _solve_cantilever_tip(make_imk)
    d_ref = _solve_cantilever_tip(make_elastic)

    for a, b in zip(d_imk, d_ref):
        assert abs(a - b) <= 1.0e-9 * (abs(b) + 1.0e-12), (a, b)

    # sanity: the tip actually moved (load really applied)
    assert any(abs(x) > 1.0e-9 for x in d_ref)


def test_no_hinge_when_material_omitted_is_pure_elastic():
    """-hinge with no -matZ must still be elastic (the hinge exists only where a
    material is supplied)."""
    def make_imk_hinge_no_mat():
        ops.element("LadrunoIMKBeam", 1, 1, 2, _A, _E, _Iz, 1, "-hinge", "both")

    def make_elastic():
        ops.element("elasticBeamColumn", 1, 1, 2, _A, _E, _Iz, 1)

    d_imk = _solve_cantilever_tip(make_imk_hinge_no_mat)
    d_ref = _solve_cantilever_tip(make_elastic)
    for a, b in zip(d_imk, d_ref):
        assert abs(a - b) <= 1.0e-9 * (abs(b) + 1.0e-12), (a, b)


# ---------------------------------------------------------------------------
# Test 2: single strong-axis hinge reproduces elastic-series + IMK backbone
# ---------------------------------------------------------------------------
_My = 1000.0     # yield moment
_k0 = 1.0e6      # initial rotational stiffness of the hinge
_bh = 0.02       # strain-hardening ratio
_theta_y = _My / _k0


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
        thj -= g / (-a - dM)
    if thj <= _theta_y:
        return thj, _k0 * thj
    return thj, _My + _bh * _k0 * (thj - _theta_y)


def test_single_hinge_series_matches_independent_model():
    a = 4.0 * _E * _Iz / _L

    _new_model()
    ops.uniaxialMaterial("Steel01", 10, _My, _k0, _bh)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 0)              # only the strong-axis rotation (dof 3) free

    ops.element("LadrunoIMKBeam", 1, 1, 2, _A, _E, _Iz, 1, "-hinge", "j", "-matZ", 10)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, 0.0, 1.0)      # reference moment on dof 3

    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    dtheta = 1.0e-4
    ops.integrator("DisplacementControl", 2, 3, dtheta)
    ops.test("NormDispIncr", 1.0e-12, 50)
    ops.algorithm("Newton")
    ops.analysis("Static")

    n_steps = 45  # spans elastic (theta_y=1e-3) well into hardening (~4.5e-3)
    for step in range(1, n_steps + 1):
        assert ops.analyze(1) == 0, f"step {step} failed"

        theta_imposed = ops.nodeDisp(2, 3)
        hingeRot = ops.eleResponse(1, "hingeRotation")    # [Zi, Zj]
        hingeMom = ops.eleResponse(1, "hingeMoment")      # [Mz_i, Mz_j]
        thj_ele = hingeRot[1]
        M_ele = hingeMom[1]

        thj_pred, M_pred = _series_predict(theta_imposed, a)

        assert abs(thj_ele - thj_pred) <= 1.0e-7 * (abs(thj_pred) + 1.0e-9), (
            f"step {step}: hinge rot {thj_ele} vs {thj_pred}")
        assert abs(M_ele - M_pred) <= 1.0e-6 * (abs(M_pred) + 1.0e-6), (
            f"step {step}: hinge moment {M_ele} vs {M_pred}")

    assert ops.nodeDisp(2, 3) > _theta_y * 2.0


# ---------------------------------------------------------------------------
# Test 3: axial immunity -- a spurious axial force leaves M-theta unchanged
# ---------------------------------------------------------------------------
_axN = 1.0e6     # spurious axial force injected through an equalDOF tie
_axM = 2500.0    # bending moment demand (yields the My=1000 hinge)


def _axial_run(inject, n_steps=40):
    """Cantilever bending about z at end j; optionally inject a large axial force
    through an equalDOF tie of node-2 ux to a driver node pulled in x. Returns
    per-step (hingeRotation_Zj, hingeMoment_Zj, axialN) or None on failure."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)            # i-end, fixed
    ops.node(2, _L, 0.0)             # j-end (hinge), rz free
    ops.node(3, 0.0, 0.0)            # axial driver, pulled in x
    ops.geomTransf("Linear", 1)

    ops.uniaxialMaterial("Steel01", 10, 1000.0, 1.0e6, 0.03)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 0, 1, 0)              # ux tied to driver; uy fixed; rz free
    ops.fix(3, 0, 1, 1)             # driver: only ux free (axial source)
    ops.equalDOF(3, 2, 1)           # tie node-2 ux to driver node-3 ux

    ops.element("LadrunoIMKBeam", 1, 1, 2, _A, _E, _Iz, 1, "-hinge", "j", "-matZ", 10)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, 0.0, _axM)     # bending moment about z at j
    if inject:
        ops.load(3, _axN, 0.0, 0.0)  # axial pull at the driver
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
        thh = ops.eleResponse(1, "hingeRotation")[1]   # Zj
        M = ops.eleResponse(1, "hingeMoment")[1]
        N = ops.eleResponse(1, "basicForce")[0]
        hist.append((thh, M, N))
    return hist


def test_axial_immunity_moment_rotation_unaffected():
    n = 40
    with_axial = _axial_run(True, n)
    no_axial = _axial_run(False, n)

    assert with_axial is not None and len(with_axial) == n
    assert no_axial is not None and len(no_axial) == n

    # the tie really injected a large axial force (and ~none without it)
    assert abs(with_axial[-1][2] - _axN) <= 1.0e-3 * _axN
    assert abs(no_axial[-1][2]) <= 1.0e-3 * _axN

    # the hinge actually yielded (so plasticity, not elastic, is exercised)
    assert max(abs(M) for _, M, _ in no_axial) > 1.2 * 1000.0

    # THE GUARANTEE: moment-rotation is unchanged by the spurious axial force
    for (thA, MA, _), (thB, MB, _) in zip(with_axial, no_axial):
        assert abs(thA - thB) <= 1.0e-9 * (abs(thB) + 1.0e-9), (thA, thB)
        assert abs(MA - MB) <= 1.0e-7 * (abs(MB) + 1.0e-3), (MA, MB)


# ---------------------------------------------------------------------------
# Test 4: column-face offset -- hinge moment vs joint moment differ by V*offset
# ---------------------------------------------------------------------------
def _offset_run(off, V=1.0e3):
    """Cantilever along X with a member-axis rigid offset at end i (via the 2D
    geomTransf -jntOffset). Returns (V, hingeMoment_Zi, jointReaction_Mz)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, _L, 0.0)
    # rigid offset at the i-end along +X (member axis); 2D jntOffset: dXi,dYi,dXj,dYj
    ops.geomTransf("Linear", 1, "-jntOffset", off, 0.0, 0.0, 0.0)
    ops.uniaxialMaterial("Steel01", 10, 1.0e9, 1.0e8, 0.01)   # high yield -> elastic
    ops.fix(1, 1, 1, 1)
    ops.element("LadrunoIMKBeam", 1, 1, 2, _A, _E, _Iz, 1, "-hinge", "i", "-matZ", 10)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, V, 0.0)        # transverse Y -> strong-axis bending
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
    Rz1 = ops.nodeReaction(1, 3)                   # Mz reaction at the joint node
    return V, Mz_i, Rz1


def test_column_face_offset_moment_vs_joint_moment():
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
# Test 5: asymmetric end laws -- -matZi / -matZj put an independent law at each end
# ---------------------------------------------------------------------------
def test_asymmetric_end_laws_each_end_uses_its_own_material():
    MyA, MyB, Ke, b = 800.0, 1600.0, 1.0e6, 0.02   # distinct yield moments

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, _L, 0.0)
    ops.geomTransf("Linear", 1)
    ops.uniaxialMaterial("Steel01", 31, MyA, Ke, b)   # law A -> end i
    ops.uniaxialMaterial("Steel01", 32, MyB, Ke, b)   # law B -> end j
    ops.fix(1, 1, 1, 0)             # rz free + moment-driven at both ends
    ops.fix(2, 1, 1, 0)
    ops.element("LadrunoIMKBeam", 1, 1, 2, _A, _E, _Iz, 1, "-matZi", 31, "-matZj", 32)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0, 0, 1300.0)       # yields A (My=800)
    ops.load(2, 0, 0, 2400.0)       # yields B (My=1600)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 1.0 / 30.0)
    ops.test("NormDispIncr", 1.0e-12, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")

    thi_hist, thj_hist, Mi_hist, Mj_hist = [], [], [], []
    for _ in range(30):
        assert ops.analyze(1) == 0
        hr = ops.eleResponse(1, "hingeRotation")   # [Zi, Zj]
        hm = ops.eleResponse(1, "hingeMoment")
        thi_hist.append(hr[0]); thj_hist.append(hr[1])
        Mi_hist.append(hm[0]);  Mj_hist.append(hm[1])

    # both ends must have yielded (so the distinct laws are actually exercised)
    assert max(abs(t) for t in thi_hist) > 2.0 * (MyA / Ke)
    assert max(abs(t) for t in thj_hist) > 2.0 * (MyB / Ke)
    # and the two ends really carry different moments (distinct laws, not a copy)
    assert abs(max(abs(m) for m in Mi_hist) - max(abs(m) for m in Mj_hist)) > 100.0

    # end i must follow law A (tag 31); end j must follow law B (tag 32)
    ops.wipe()
    ops.uniaxialMaterial("Steel01", 31, MyA, Ke, b)
    ops.uniaxialMaterial("Steel01", 32, MyB, Ke, b)
    ops.testUniaxialMaterial(31)
    for th, M in zip(thi_hist, Mi_hist):
        ops.setStrain(th)
        assert abs(M - ops.getStress()) <= 1.0e-6 * abs(MyB) + 1.0e-6, ("i", M)
    ops.testUniaxialMaterial(32)
    for th, M in zip(thj_hist, Mj_hist):
        ops.setStrain(th)
        assert abs(M - ops.getStress()) <= 1.0e-6 * abs(MyB) + 1.0e-6, ("j", M)


# ---------------------------------------------------------------------------
# Test 6: Lignos-Krawinkler cyclic -- steel Bilin (modified-IMK) end hinge
# ---------------------------------------------------------------------------
_cyc_My, _cyc_Ke, _cyc_b = 1000.0, 1.0e6, 0.03


def _build_bilin(tag):
    ops.uniaxialMaterial("Bilin", tag,
                         _cyc_Ke, _cyc_b, _cyc_b, _cyc_My, -_cyc_My,
                         1.0, 1.0, 1.0, 1.0,
                         1.0, 1.0, 1.0, 1.0,
                         0.03, 0.03,
                         0.20, 0.20,
                         0.4, 0.4,
                         0.40, 0.40,
                         1.0, 1.0)


def _cyclic_targets(amps=(0.004, 0.008, 0.016, 0.028), ncyc=2):
    targets = []
    for amp in amps:
        for _ in range(ncyc):
            targets += [amp, -amp]
    targets += [0.0]
    return targets


def test_lignos_krawinkler_cyclic_matches_bare_material():
    dth = 4.0e-4
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, _L, 0.0)
    ops.geomTransf("Linear", 1)
    _build_bilin(10)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 0)             # only strong-axis rotation free
    ops.element("LadrunoIMKBeam", 1, 1, 2, _A, _E, _Iz, 1, "-hinge", "j", "-matZ", 10)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0, 0, 1.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-12, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("DisplacementControl", 2, 3, dth)
    ops.analysis("Static")

    hist = []           # (theta_node, theta_hinge, moment)
    cur = 0.0
    for tgt in _cyclic_targets():
        step = dth if tgt >= cur else -dth
        nstep = max(1, int(round(abs(tgt - cur) / dth)))
        ops.integrator("DisplacementControl", 2, 3, step)
        for _ in range(nstep):
            assert ops.analyze(1) == 0, "cyclic element analysis failed"
            th2 = ops.nodeDisp(2, 3)
            thh = ops.eleResponse(1, "hingeRotation")[1]
            M = ops.eleResponse(1, "hingeMoment")[1]
            hist.append((th2, thh, M))
        cur = tgt

    Mscale = max(abs(M) for _, _, M in hist)
    thh_max = max(abs(thh) for _, thh, _ in hist)
    assert Mscale > 1.2 * _cyc_My          # well past yield (deterioration active)
    assert thh_max > 5.0 * (_cyc_My / _cyc_Ke)

    # (i) elastic-interior series identity holds cyclically and exactly
    a = 4.0 * _E * _Iz / _L
    series_err = max(abs(a * (th2 - thh) - M) for th2, thh, M in hist)
    assert series_err <= 1.0e-6 * Mscale, series_err

    # (ii) the element faithfully delegates the FULL cyclic path to the hinge
    #      material: driving the bare Bilin with the reported hinge-rotation
    #      history reproduces the reported hinge moment at every step.
    ops.wipe()
    _build_bilin(99)
    ops.testUniaxialMaterial(99)
    for (_, thh, M) in hist:
        ops.setStrain(thh)
        bare = ops.getStress()
        assert abs(M - bare) <= 1.0e-6 * (Mscale) + 1.0e-6, (M, bare)

    # (iii) real hysteretic dissipation occurred
    energy = 0.0
    for k in range(1, len(hist)):
        energy += 0.5 * (hist[k][2] + hist[k - 1][2]) * (hist[k][1] - hist[k - 1][1])
    assert energy > 0.1 * _cyc_My * thh_max
