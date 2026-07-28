"""LadrunoBrick ``-geom hypo`` (v4, ADR 79) — hypoelastic rate-form UL battery.

``hypo`` integrates the Hughes-Winget midpoint objective strain increment in
the UNROTATED (Green-Naghdi) material frame and drives an UNCHANGED
small-strain material via accumulate-and-setTrialStrain — large strain with
the existing material library (the route for rate-form soil laws). Assembly
pushes the material stress/tangent forward with R = polar(F) on the current
configuration.

Gates (ADR 79 P1):
  1. RIGID ROTATION => ZERO STRESS through a large multi-step rotation — the
     HW midpoint makes the per-step increment EXACTLY zero for a rigid
     increment of any magnitude, so this holds to solver precision.
  2. hypo reduces to ``linear`` at infinitesimal strain.
  3. Load-driven cantilever agrees with -geom corot AND -geom finite in the
     large-rotation / small-strain regime (drives Newton through the hypo
     tangent).
  4. Uniaxial-strain compression to 25%%: exact vs the analytic hypoelastic
     target sigma = (lam+2mu)*ln(lambda) (Cauchy-rate convention), and the
     measured gap to -geom finite + LogStrain equals the 1/J work-conjugacy
     factor (ADR 79 §3.1) — the discrepancy asserted as a verified property.
  5. Elastic simple shear to gamma = 2 vs the numpy oracle integrating the
     SAME algorithm (element-vs-oracle exact; monotone GN signature — no
     Jaumann oscillation).
  6. J2 kinematic hardening objectivity in the hypo path's OWN convention
     (fixed unrotated material frame): deform-then-rigid-rotate leaves the
     material-frame stress EXACTLY frozen; simultaneous rotation+stretch
     agrees to the documented O(dtheta*deps) incremental-objectivity error.
  7. Database round-trip: the per-GP committed feed strain (the guarded
     Vector(48) in sendSelf) survives save/restore — the next step's stress
     after a restore matches the uninterrupted run.
"""
import math
import os
import tempfile

import numpy as np
import pytest

from _testbed import ops
import hypo_reference as href

pytestmark = [pytest.mark.zone_a]

_E = 200.0
_NU = 0.3
_LAM = _E * _NU / ((1.0 + _NU) * (1.0 - 2.0 * _NU))
_G = _E / (2.0 * (1.0 + _NU))

# distorted reference hex (the corot battery geometry) for the frame gates
_NODES_DIST = {
    1: (0.00, 0.00, 0.00),
    2: (1.00, 0.10, 0.05),
    3: (1.10, 1.00, 0.00),
    4: (0.05, 0.95, 0.10),
    5: (0.00, 0.05, 1.00),
    6: (1.00, 0.00, 1.05),
    7: (1.05, 1.00, 1.10),
    8: (0.00, 1.00, 0.95),
}
# unit cube for the homogeneous large-strain gates
_NODES_CUBE = {
    1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0), 3: (1.0, 1.0, 0.0), 4: (0.0, 1.0, 0.0),
    5: (0.0, 0.0, 1.0), 6: (1.0, 0.0, 1.0), 7: (1.0, 1.0, 1.0), 8: (0.0, 1.0, 1.0),
}
_CONN = [1, 2, 3, 4, 5, 6, 7, 8]


def _rot(thz, thy):
    cz, sz = math.cos(thz), math.sin(thz)
    Rz = np.array([[cz, -sz, 0.0], [sz, cz, 0.0], [0.0, 0.0, 1.0]])
    cy, sy = math.cos(thy), math.sin(thy)
    Ry = np.array([[cy, 0.0, sy], [0.0, 1.0, 0.0], [-sy, 0.0, cy]])
    return Ry @ Rz


def _affine_disp(nodes, A):
    u = np.zeros(24)
    I = np.eye(3)
    for tag, (x, y, z) in nodes.items():
        X = np.array([x, y, z])
        u[(tag - 1) * 3:(tag - 1) * 3 + 3] = (np.asarray(A) - I) @ X
    return u


def _build(nodes, geom="hypo", material="elastic"):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in nodes.items():
        ops.node(tag, x, y, z)
    if material == "elastic":
        ops.nDMaterial("ElasticIsotropic", 1, _E, _NU)
        mtag = 1
    elif material == "logstrain":
        ops.nDMaterial("ElasticIsotropic", 1, _E, _NU)
        ops.nDMaterial("LogStrain", 2, 1)
        mtag = 2
    elif material == "j2kin":
        K = _E / (3.0 * (1.0 - 2.0 * _NU))
        ops.nDMaterial("LadrunoJ2", 1, K, _G, "-iso", "voce", 0.30, 0.10, 20.0, 5.0,
                       "-kin", 2, 8.0, 50.0, 4.0, 20.0)
        mtag = 1
    else:
        raise ValueError(material)
    ops.element("LadrunoBrick", 1, *_CONN, mtag, "-formulation", "std", "-geom", geom)
    ops.constraints("Lagrange")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-11, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    ops.timeSeries("Constant", 1)   # factor == 1 => sp is the absolute target


def _step_to(nodes, A, k):
    """Impose the affine configuration A as a COMMITTED step (per-step pattern,
    the corot-battery idiom: committed material state persists across steps)."""
    u = _affine_disp(nodes, A)
    pat = 100 + k
    ops.pattern("Plain", pat, 1)
    for tag in nodes:
        base = (tag - 1) * 3
        for d in range(3):
            ops.sp(tag, d + 1, float(u[base + d]))
    ok = ops.analyze(1)
    assert ok == 0, f"hypo step {k} failed to converge"
    ops.remove("loadPattern", pat)
    return np.array(ops.eleResponse(1, "stresses"), dtype=float)


# --------------------------------------------------------------------------- #
#  Gate 1: multi-step large rigid rotation => zero stress, zero force.        #
# --------------------------------------------------------------------------- #
def test_hypo_rigid_rotation_is_stress_free():
    _build(_NODES_DIST, "hypo")
    # 5 committed steps to 2.5 rad about a skew axis (each increment 0.5 rad —
    # far beyond any small-increment assumption; the HW midpoint zeroes it
    # exactly per step).
    for k, th in enumerate((0.5, 1.0, 1.5, 2.0, 2.5)):
        s = _step_to(_NODES_DIST, _rot(th, 0.4 * th), k)
        smax = np.abs(s).max()
        assert smax <= 1.0e-7 * _E, (
            f"rigid rotation induced stress {smax:.3e} at step {k}")
    f = np.array(ops.eleForce(1), dtype=float)
    assert np.abs(f).max() <= 1.0e-7 * _E, (
        f"rigid rotation induced force {np.abs(f).max():.3e}")


# --------------------------------------------------------------------------- #
#  Gate 2: hypo reduces to linear at infinitesimal strain.                    #
# --------------------------------------------------------------------------- #
def test_hypo_reduces_to_linear_at_small_deformation():
    eps = 1.0e-5
    A = [[1.0 + eps, 0.5 * eps, 0.0],
         [0.5 * eps, 1.0 - 0.3 * eps, 0.0],
         [0.0, 0.0, 1.0 + 0.2 * eps]]
    forces = {}
    for geom in ("hypo", "linear"):
        _build(_NODES_DIST, geom)
        _step_to(_NODES_DIST, A, 0)
        forces[geom] = np.array(ops.eleForce(1), dtype=float)
    scale = max(np.abs(forces["linear"]).max(), 1.0e-30)
    assert np.abs(forces["hypo"] - forces["linear"]).max() <= 1.0e-3 * scale, (
        "hypo did not reduce to linear at small deformation: max diff "
        f"{np.abs(forces['hypo'] - forces['linear']).max():.3e} vs {scale:.3e}")


# --------------------------------------------------------------------------- #
#  Gate 3: load-driven cantilever — Newton through the hypo tangent, agreement #
#  with corot and finite in the large-rotation / small-strain regime.          #
# --------------------------------------------------------------------------- #
def _solve_cantilever(geom, nz=6, P=0.6, nsteps=20):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    corners = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    for k in range(nz + 1):
        for c, (x, y) in enumerate(corners):
            ops.node(4 * k + c + 1, x, y, float(k))
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU)
    if geom == "finite":
        ops.nDMaterial("LogStrain", 2, 1)
        mtag = 2
    else:
        mtag = 1
    for k in range(nz):
        bs = 4 * k
        conn = [bs + 1, bs + 2, bs + 3, bs + 4, bs + 5, bs + 6, bs + 7, bs + 8]
        ops.element("LadrunoBrick", k + 1, *conn, mtag, "-formulation", "std",
                    "-geom", geom)
    for c in range(4):
        ops.fix(c + 1, 1, 1, 1)
    top = [4 * nz + c + 1 for c in range(4)]
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in top:
        ops.load(t, P / 4.0, 0.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1.0e-8, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    for _ in range(nsteps):
        if ops.analyze(1) != 0:
            return None
    return float(np.mean([ops.nodeDisp(t, 1) for t in top]))


def test_hypo_cantilever_matches_corot_and_finite():
    ux_hypo = _solve_cantilever("hypo")
    assert ux_hypo is not None, "hypo cantilever Newton failed to converge"
    assert ux_hypo > 0.0 and math.isfinite(ux_hypo)

    ux_cor = _solve_cantilever("corot")
    assert ux_cor is not None
    ux_fin = _solve_cantilever("finite")
    assert ux_fin is not None

    rel_c = abs(ux_hypo - ux_cor) / abs(ux_cor)
    rel_f = abs(ux_hypo - ux_fin) / abs(ux_fin)
    print(f"cantilever tip ux: hypo={ux_hypo:.5e} corot={ux_cor:.5e} "
          f"finite={ux_fin:.5e} rel_c={rel_c:.3e} rel_f={rel_f:.3e}")
    assert rel_c <= 2.0e-2, f"hypo vs corot tip disp rel {rel_c:.3e}"
    assert rel_f <= 2.0e-2, f"hypo vs finite tip disp rel {rel_f:.3e}"


# --------------------------------------------------------------------------- #
#  Gate 4: uniaxial-strain compression to 25% — analytic hypoelastic target,  #
#  and the 1/J gap vs -geom finite + LogStrain asserted as a property.        #
# --------------------------------------------------------------------------- #
def _uniaxial_strain_run(geom, material, lam_final=0.75, nsteps=50):
    _build(_NODES_CUBE, geom, material)
    for k in range(1, nsteps + 1):
        lam = 1.0 + (lam_final - 1.0) * k / nsteps
        A = np.diag([lam, 1.0, 1.0])
        s = _step_to(_NODES_CUBE, A, k)
    return s  # last committed GP stresses (48)


def test_hypo_uniaxial_strain_analytic_and_J_factor():
    lam = 0.75
    s_hypo = _uniaxial_strain_run("hypo", "elastic", lam)
    # homogeneous: all 8 GPs identical
    sxx = s_hypo[0::6]
    assert np.abs(sxx - sxx[0]).max() <= 1.0e-9 * abs(sxx[0])
    target = (_LAM + 2.0 * _G) * math.log(lam)     # Cauchy-rate hypoelastic, exact
    assert abs(sxx[0] - target) <= 5.0e-4 * abs(target), (
        f"hypo uniaxial-strain sigma_xx {sxx[0]:.6e} != (lam+2mu)*ln(lambda) "
        f"{target:.6e}")

    # the one case finite+LogStrain can also do: same path, Kirchhoff-based —
    # sigma_finite = sigma_hypo / J with J = lambda (lateral fixed). The 1/J
    # work-conjugacy gap is MEASURED and asserted, per ADR 79 §3.1.
    s_fin = _uniaxial_strain_run("finite", "logstrain", lam)
    sxx_f = s_fin[0]
    ratio = sxx[0] / sxx_f
    assert abs(ratio - lam) <= 1.0e-3, (
        f"hypo/finite stress ratio {ratio:.5f} != J = {lam} — the Cauchy-vs-"
        f"Kirchhoff rate-law gap is not the expected 1/J factor")


# --------------------------------------------------------------------------- #
#  Gate 5: elastic simple shear to gamma = 2 — element == oracle, GN          #
#  signature (monotone shear stress, no Jaumann oscillation).                 #
# --------------------------------------------------------------------------- #
def _elastic_D():
    D = np.zeros((6, 6))
    for i in range(3):
        for j in range(3):
            D[i, j] = _LAM + (2.0 * _G if i == j else 0.0)
    for i in range(3, 6):
        D[i, i] = _G
    return D


def test_hypo_simple_shear_vs_oracle():
    N = 40
    gmax = 2.0
    _build(_NODES_CUBE, "hypo", "elastic")

    # oracle: the SAME algorithm on the homogeneous affine path (any linear-
    # complete gradient set gives the same H for an affine field — use the
    # oracle's frozen centroid gradients).
    D = _elastic_D()
    X = href.CORNERS
    eps_feed = np.zeros(6)
    s_hist_el, s_hist_or = [], []
    A_prev = np.eye(3)
    for k in range(1, N + 1):
        gam = gmax * k / N
        A = np.array([[1.0, gam, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        # oracle step
        un = (X @ A_prev.T) - X
        u1 = (X @ A.T) - X
        deps6, R1, J1 = href.hypo_increment(href.DNDX_CENTROID, un, u1)
        eps_feed += deps6
        s_hist_or.append(D @ eps_feed)          # material (unrotated-frame) stress
        # element step
        s = _step_to(_NODES_CUBE, A, k)
        s_hist_el.append(s[:6])
        A_prev = A

    s_el = np.array(s_hist_el)
    s_or = np.array(s_hist_or)
    scale = np.abs(s_or).max()
    assert np.abs(s_el - s_or).max() <= 1.0e-8 * scale, (
        "hypo element diverged from the oracle integrating the same algorithm: "
        f"max diff {np.abs(s_el - s_or).max():.3e} vs scale {scale:.3e}")

    # GN signature: the material-frame shear response is monotone up to
    # gamma = 2 (Jaumann would oscillate at large shear; GN does not).
    sxy = s_el[:, 3]
    assert np.all(np.diff(np.abs(sxy)) > 0.0), (
        "hypo shear stress not monotone — unexpected oscillation (Jaumann-like)")


# --------------------------------------------------------------------------- #
#  Gate 6: J2 kinematic hardening objectivity — the hypo path's OWN           #
#  convention (fixed unrotated material frame).                               #
# --------------------------------------------------------------------------- #
def _dev_stretch(a):
    return np.diag([1.0 + a, 1.0 / math.sqrt(1.0 + a), 1.0 / math.sqrt(1.0 + a)])


def test_hypo_j2_kin_deform_then_rigid_rotate_freezes_state():
    # yield the material (backstress develops), then rigidly rotate in steps:
    # every rigid increment gives deps == 0 EXACTLY, so the material-frame
    # stress must stay frozen to solver precision.
    _build(_NODES_DIST, "hypo", "j2kin")
    k = 0
    for a in (0.03, 0.06, 0.09):
        s_deformed = _step_to(_NODES_DIST, _dev_stretch(a), k)
        k += 1
    U = _dev_stretch(0.09)
    scale = np.abs(s_deformed).max()
    assert scale > 0.30, f"j2kin path did not yield (max |sigma| {scale:.3f}); vacuous"

    for th in (0.4, 0.8, 1.2):
        s_rot = _step_to(_NODES_DIST, _rot(th, -0.5 * th) @ U, k)
        k += 1
        assert np.abs(s_rot - s_deformed).max() <= 1.0e-6 * scale, (
            f"material-frame stress moved under a rigid increment (th={th}): "
            f"max diff {np.abs(s_rot - s_deformed).max():.3e} vs {scale:.3e}")


def test_hypo_j2_kin_simultaneous_rotation_objectivity():
    # simultaneous rotation + stretch: incremental objectivity holds to
    # O(dtheta * deps) per step (ADR 79 §3 — NOT an exact identity like the
    # rigid case above). 20 steps to 1.0 rad + 9% deviatoric stretch; the
    # documented error scale is ~sum(dtheta*deps) ~ 5e-3, gate at 2e-2.
    N = 20
    hist = {}
    for rotate in (False, True):
        _build(_NODES_DIST, "hypo", "j2kin")
        for k in range(1, N + 1):
            a = 0.09 * k / N
            Q = _rot(1.0 * k / N, -0.4 * k / N) if rotate else np.eye(3)
            s = _step_to(_NODES_DIST, Q @ _dev_stretch(a), k)
        hist[rotate] = s
    scale = np.abs(hist[False]).max()
    assert scale > 0.30, "unrotated j2kin path did not yield; vacuous"
    rel = np.abs(hist[True] - hist[False]).max() / scale
    print(f"hypo j2kin simultaneous-rotation objectivity rel err = {rel:.3e}")
    assert rel <= 2.0e-2, (
        f"material-frame stress under simultaneous rotation differs by {rel:.3e} "
        "— beyond the documented incremental-objectivity error scale")


# --------------------------------------------------------------------------- #
#  Gate 7: database round-trip preserves the committed feed strain.           #
#  Elastic material => stress is a direct function of eps_feed, so a lost /   #
#  zeroed Vector(48) makes the post-restore step wildly wrong.                #
# --------------------------------------------------------------------------- #
def test_hypo_database_roundtrip_preserves_feed_strain():
    steps = [(0.02, 0.2), (0.04, 0.4), (0.06, 0.6)]

    def drive(k0, subset):
        s = None
        for k, (a, th) in enumerate(subset, start=k0):
            s = _step_to(_NODES_DIST, _rot(th, -0.3 * th) @ _dev_stretch(a), k)
        return s

    with tempfile.TemporaryDirectory(prefix="ladruno_hypo_rt_",
                                     ignore_cleanup_errors=True) as td:
        dbpath = os.path.join(td, "hypo_rt")

        _build(_NODES_DIST, "hypo")
        drive(0, steps[:2])
        try:
            ops.database("File", dbpath)
        except Exception as exc:  # noqa: BLE001 - build without FE_Datastore
            pytest.skip(f"database() unsupported in this build: {exc}")
        saved = ops.save(1)
        if saved is not None and saved < 0:
            pytest.skip("database save returned failure on this build")
        s_ref = drive(2, steps[2:])            # uninterrupted continuation

        ops.wipe()
        _build(_NODES_DIST, "hypo")            # fresh skeleton, virgin state
        ops.database("File", dbpath)
        ops.restore(1)
        s_rt = drive(2, steps[2:])             # same step from restored state
        ops.wipe()                             # release FE_Datastore handles

    scale = max(np.abs(s_ref).max(), 1.0e-30)
    assert np.abs(s_rt - s_ref).max() <= 1.0e-9 * scale, (
        "post-restore step stress differs from the uninterrupted run — the "
        "hypo committed feed strain did not survive sendSelf/recvSelf "
        f"(max diff {np.abs(s_rt - s_ref).max():.3e} vs scale {scale:.3e})")
