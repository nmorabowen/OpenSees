"""LadrunoBrick20 dynamics — HRZ ``-lumped`` mass + explicit proof-of-life.
Zone-A P3 battery (ADR 72 §6 P3 row; plan task 3.2).

The P3 contract, in test order:

  1. EIGEN BRACKETING (fixed-free axial rod, exact 1-D analytic f1 = c/4L):
     consistent-mass frequency ABOVE analytic, HRZ BELOW, both CONVERGING
     (asserted on a coarse/fine mesh pair). Bending cantilever: consistent and
     HRZ bracket each other (HRZ below), both within 0.5% of the
     Euler-Bernoulli f1 and HRZ below it. NB the consistent-mass bending
     frequency sits ~2e-4 BELOW the EB value — the discrete model carries
     shear deformation EB theory omits, a physical (Timoshenko-direction)
     offset, not a mass-matrix artifact; the from-above bracket is exact
     against the shear-free axial analytic.

  2. ELEMENT-PATH HRZ VECTOR: on the unit cube, the ``-lumped`` getMass
     diagonal equals the P0 pinned fractions EXACTLY (corners 7/248*M,
     mid-edges 2/31*M, <=1e-14 rel; ADR 72 §3.5 / Cook §13.3 method),
     off-diagonals are EXACTLY zero, every entry positive, and the trace
     conserves 3*M_total. Extraction: two-dt Newmark-tangent subtraction
     (A(dt1)-A(dt2))/(c3(dt1)-c3(dt2)) — K cancels, M remains.

  3. criticalTimeStep(): equals the numpy 60-DOF pencil (max generalized eig
     of the oracle K vs the oracle HRZ-lumped M, dt_cr = 2/omega_max) to
     <=1e-10 rel, on the unit cube AND the P1 distorted hex, std AND uri.
     NOTE (deviation adjudicated in-file): the SIBLING LadrunoBrick has NO
     criticalTimeStep element method either — the fork architecture computes
     the per-element pencil exactly in CriticalTimeStep.cpp (ADR 65), and
     Element::getExplicitCriticalTimeStep's contract FORBIDS overriding when
     the element pencil is expressible (it REPLACES, not min-folds). This
     test pins the queryable ops.criticalTimeStep() value against the
     independent numpy pencil — the ADR §6 gate as stated.

  4. EXPLICIT WAVE BAR: 10-element H20 bar under CentralDifferenceLadruno,
     ``-lumped``, ~0.9x pencil dt, force-pulse excitation (NOT velocity ICs —
     the ADR-69 velocity-IC gates pin ERR~-100), 2500 steps: stable and
     bounded. The measured dt_cr ratio vs an equal-node-spacing LadrunoBrick
     H8 mesh is REPORTED (loose sanity band only): measured ~0.50 — BELOW
     the 1-D rod-theory 2/sqrt(6)~0.82 ballpark because the 3-D HRZ corner
     masses (7/248 vs the rod's 1/6-class shares) raise omega_max further.

  5. ENERGY BALANCE (ADR-69 recorder) on the wave bar: KE/IE closure — RES
     stays small vs the peak energy and the total drifts <1% over the run;
     the -v2 channel layout declares NO hourglass channel (E_hg absent):
     this element has no hourglass machinery AT ALL, so the channel must not
     even exist, which is stronger than a zero column.

  6. betaK-RAYLEIGH CLOBBER GATE (the ADR-66/#562 family bug):
     getResistingForceIncInertia under rayleigh(0, betaK, 0, 0) must keep
     M*a and -Q. Self-validating pattern from the plane-family x4 gate
     (test_ladrunoplane_dynamics): the undamped step-load run shows the ~2x
     dynamic overshoot vs static, then a tiny betaK moves the peak <5%.
     Under the bug the damped peak collapses toward 1x static. Parametrized
     {std,uri} x {consistent,-lumped}.

  7. massType SERIALIZATION: database round-trip of a ``-lumped`` element
     preserves the lumped mass through recvSelf (eigenvalues as the
     mass-sensitive probe; guarded non-vacuous against the consistent-mass
     eigenvalues).

Quirks honored: force excitation not velocity ICs (energy gates); dt held
CONSTANT all run (a dt change needs a revertToLastStep reseed); tables
imported from tests/hex20_reference.py (debt c — never re-transcribe).
"""
import math
import os

import pytest

from _testbed import ops

np = pytest.importorskip("numpy")
pytest.importorskip("sympy")             # hex20_reference derives symbolically
import hex20_reference as ref            # noqa: E402  (P0 oracle)

pytestmark = [pytest.mark.zone_a]

E, NU = 1000.0, 0.3
RHO = 2.0e-3

_CONN = list(range(1, 21))
_CUBE = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
         (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
_DIST = [(0.00, 0.00, 0.00), (1.00, 0.10, 0.00), (1.10, 1.00, 0.10), (0.05, 0.95, 0.00),
         (0.00, 0.05, 1.00), (1.00, 0.00, 1.05), (1.05, 1.00, 1.10), (0.00, 1.00, 0.95)]


# --------------------------------------------------------------------------
# fixtures — straight-edged hex20s from the oracle's tables (debt c)
# --------------------------------------------------------------------------
def _hex20_nodes(corners):
    """20 node coords (list, brcshl order) from 8 corners; mid-edges at exact
    midpoints via the oracle EDGES table (0-based corner pairs)."""
    nodes = [tuple(map(float, corners[i])) for i in range(8)]
    for a, b in ref.EDGES:
        ca, cb = corners[a], corners[b]
        nodes.append(tuple(0.5 * (ca[d] + cb[d]) for d in range(3)))
    return nodes


def _box(x0, x1, y0, y1, z0, z1):
    return [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
            (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)]


class _Mesh:
    """Coordinate-dedup node factory for multi-element hex20/hex8 bars."""

    def __init__(self):
        self.tag_of = {}
        self.coords = {}
        self.next = 1

    def add(self, c):
        key = tuple(round(v, 9) for v in c)
        if key not in self.tag_of:
            ops.node(self.next, *c)
            self.tag_of[key] = self.next
            self.coords[self.next] = c
            self.next += 1
        return self.tag_of[key]


def _h20_bar(nx, lumped, e_mod=E, nu=0.0, rho=1.0, sy=1.0, sz=1.0, hx=1.0,
             formulation="std"):
    """nx-element H20 bar along x, section sy x sz; returns the mesh."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, e_mod, nu, rho)
    m = _Mesh()
    extra = ["-formulation", formulation] + (["-lumped"] if lumped else [])
    for e in range(nx):
        conn = [m.add(c) for c in _hex20_nodes(_box(e * hx, (e + 1) * hx, 0, sy, 0, sz))]
        ops.element("LadrunoBrick20", e + 1, *conn, 1, *extra)
    return m


def _clamp_root(m, axial_only=False):
    """Fix the x=0 face; optionally pin y,z everywhere else (pure 1-D axial)."""
    for t, c in m.coords.items():
        if abs(c[0]) < 1e-9:
            ops.fix(t, 1, 1, 1)
        elif axial_only:
            ops.fix(t, 0, 1, 1)


# ==========================================================================
# 1. eigen bracketing — consistent ABOVE analytic, HRZ BELOW, both converging
# ==========================================================================
def _axial_f1(nx, lumped):
    m = _h20_bar(nx, lumped)                     # E=1000, nu=0, rho=1
    _clamp_root(m, axial_only=True)
    lam = ops.eigen("-fullGenLapack", 1)
    return math.sqrt(lam[0]) / (2.0 * math.pi)


def test_eigen_bracketing_axial_rod():
    """Fixed-free rod f1 = c/4L (exact 1-D analytic, nu=0): the consistent
    mass converges FROM ABOVE (Hughes §7.3.2), HRZ from BELOW, and both
    errors shrink under refinement. Measured (this build): nx=2 cons/analytic
    1.00026, HRZ 0.95358; nx=6 cons 1.00000(+), HRZ 0.99466."""
    c = math.sqrt(1000.0 / 1.0)
    results = {}
    for nx in (2, 6):
        L = float(nx)
        f1 = c / (4.0 * L)
        f_cons = _axial_f1(nx, lumped=False)
        f_hrz = _axial_f1(nx, lumped=True)
        # the bracket
        assert f_cons >= f1 * (1.0 - 1e-9), (
            f"nx={nx}: consistent f={f_cons} must sit ABOVE analytic {f1}"
        )
        assert f_hrz < f1, (
            f"nx={nx}: HRZ f={f_hrz} must sit BELOW analytic {f1}"
        )
        results[nx] = (f_cons / f1, f_hrz / f1)
    # convergence: both errors shrink from nx=2 to nx=6
    for i, name in ((0, "consistent"), (1, "HRZ")):
        e2 = abs(results[2][i] - 1.0)
        e6 = abs(results[6][i] - 1.0)
        assert e6 < e2, (
            f"{name} axial eigenfrequency not converging: |err| {e2} -> {e6}"
        )
    print(f"\n[P3 eigen axial] cons/analytic nx=2 {results[2][0]:.6f} nx=6 "
          f"{results[6][0]:.6f}; HRZ/analytic nx=2 {results[2][1]:.6f} nx=6 "
          f"{results[6][1]:.6f}")


def test_eigen_bracketing_bending_cantilever():
    """Slender cantilever (10 x 0.4 x 0.4, one honest mesh): HRZ f1 below the
    consistent f1 (the mass-matrix bracket), both within 0.5% of the
    Euler-Bernoulli analytic and HRZ BELOW it. The consistent value sits
    ~2e-4 below EB — shear deformation the EB reference omits (see module
    docstring); the strict from-above assertion lives on the shear-free
    axial gate above. Measured: cons/EB 0.9998, HRZ/EB 0.9966."""
    nx, sy, sz = 10, 0.4, 0.4
    L, A, I = float(nx), sy * sz, sy * sz ** 3 / 12.0
    f_eb = (1.875104 ** 2 / (2.0 * math.pi)) * math.sqrt(1000.0 * I / (1.0 * A * L ** 4))

    fs = {}
    for lumped in (False, True):
        m = _h20_bar(nx, lumped, sy=sy, sz=sz)
        _clamp_root(m)
        lam = ops.eigen("-fullGenLapack", 2)
        fs[lumped] = min(math.sqrt(x) / (2.0 * math.pi) for x in lam if x > 1e-12)

    f_cons, f_hrz = fs[False], fs[True]
    assert f_hrz < f_cons, (
        f"HRZ f1 {f_hrz} must sit below the consistent f1 {f_cons}"
    )
    assert f_hrz < f_eb, f"HRZ f1 {f_hrz} must sit below EB analytic {f_eb}"
    for name, f in (("consistent", f_cons), ("HRZ", f_hrz)):
        assert abs(f - f_eb) / f_eb < 5e-3, (
            f"{name} bending f1 {f} vs EB {f_eb}: off by "
            f"{abs(f - f_eb) / f_eb:.2%} (> 0.5%)"
        )
    print(f"\n[P3 eigen bending] cons/EB {f_cons / f_eb:.6f}  "
          f"HRZ/EB {f_hrz / f_eb:.6f}  (EB f1 = {f_eb:.6f})")


# ==========================================================================
# 2. element-path HRZ vector == the P0 cube fractions (7/248, 2/31)
# ==========================================================================
def _element_mass_matrix(corners, extra):
    """Full 60x60 element mass via two-dt Newmark-tangent subtraction on a
    FREE single element: A(dt) = K + M/(beta*dt^2), so
    (A1 - A2)/(c3_1 - c3_2) = M exactly (K cancels)."""
    nodes = _hex20_nodes(corners)

    def tangent(dt):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        for i, c in enumerate(nodes):
            ops.node(i + 1, *c)
        ops.nDMaterial("ElasticIsotropic", 1, E, NU, RHO)
        ops.element("LadrunoBrick20", 1, *_CONN, 1, *extra)
        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("FullGeneral")
        ops.algorithm("Linear")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient")
        assert ops.analyze(1, dt) == 0
        a = np.array(ops.printA("-ret"), dtype=float)
        n = int(round(len(a) ** 0.5))
        return a.reshape(n, n)

    dt1, dt2, beta = 1.0e-2, 2.0e-2, 0.25
    c3_1, c3_2 = 1.0 / (beta * dt1 * dt1), 1.0 / (beta * dt2 * dt2)
    return (tangent(dt1) - tangent(dt2)) / (c3_1 - c3_2)


@pytest.mark.parametrize("extra", [["-lumped"], ["-formulation", "uri", "-lumped"]],
                         ids=["std", "uri"])
def test_hrz_vector_matches_p0_cube_fractions(extra):
    """Unit-cube -lumped element mass: EXACTLY diagonal, every entry positive,
    corners = 7/248*M and mid-edges = 2/31*M to <=1e-14 rel (the P0 pinned
    fractions, ADR 72 §3.5), trace = 3*M_total. Mass is formulation-
    independent, so std and uri produce the identical vector."""
    M = _element_mass_matrix(_CUBE, extra)
    diag = np.diag(M).copy()
    off = M - np.diag(diag)
    m_tot = RHO * 1.0                              # unit cube volume

    assert np.max(np.abs(off)) == 0.0, (
        f"lumped mass must be EXACTLY diagonal; max offdiag {np.max(np.abs(off))}"
    )
    assert (diag > 0.0).all(), "HRZ produced a non-positive diagonal entry"

    f_corner = 7.0 / 248.0 * m_tot
    f_mid = 2.0 / 31.0 * m_tot
    for a in range(8):                              # corner nodes
        for d in range(3):
            v = diag[3 * a + d]
            assert abs(v - f_corner) <= 1e-14 * f_corner, (
                f"corner node {a} dof {d}: {v} != 7/248*M = {f_corner} "
                f"(rel {abs(v - f_corner) / f_corner:.2e})"
            )
    for a in range(8, 20):                          # mid-edge nodes
        for d in range(3):
            v = diag[3 * a + d]
            assert abs(v - f_mid) <= 1e-14 * f_mid, (
                f"mid-edge node {a} dof {d}: {v} != 2/31*M = {f_mid}"
            )
    assert abs(diag.sum() - 3.0 * m_tot) <= 1e-13 * 3.0 * m_tot, (
        f"HRZ trace {diag.sum()} does not conserve 3*M_total {3.0 * m_tot}"
    )


def test_consistent_mass_unchanged_when_lumped_absent():
    """FP-honesty guard: without -lumped the element mass is the (non-diagonal)
    27-pt consistent block — the P3 change must not perturb the std path."""
    M = _element_mass_matrix(_CUBE, [])
    off = M - np.diag(np.diag(M))
    assert np.max(np.abs(off)) > 1e-8, (
        "consistent mass unexpectedly diagonal — did -lumped leak into the default?"
    )


# ==========================================================================
# 3. criticalTimeStep == the numpy 60-DOF pencil
# ==========================================================================
def _oracle_dtcr(corners, formulation):
    """dt_cr = 2/omega_max of K v = lam * M_hrz v from the sympy-derived
    oracle tables (K at the formulation rule; M = HRZ lump of the 27-pt
    consistent mass — the mass actually in use under -lumped)."""
    X = ref.straight_edge_nodes(np.array(corners, dtype=float))
    D = ref.iso_D(E, NU)
    gp = ref.GP8_F if formulation == "uri" else ref.GP27_F
    K = ref.stiffness(X, D, gp)
    M = ref.consistent_mass(X, RHO, ref.GP27_F)
    md = ref.hrz_lump(M)
    s = np.sqrt(md)
    lam_max = np.linalg.eigvalsh(K / np.outer(s, s))[-1]
    return 2.0 / math.sqrt(lam_max)


def _ops_dtcr(corners, formulation, lumped=True):
    nodes = _hex20_nodes(corners)
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, c in enumerate(nodes):
        ops.node(i + 1, *c)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU, RHO)
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1)
    extra = ["-formulation", formulation] + (["-lumped"] if lumped else [])
    ops.element("LadrunoBrick20", 1, *_CONN, 1, *extra)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(7, 0.0, 0.0, 1.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator("ExplicitBathe", 0.54, "-cfl")
    ops.analysis("Transient")
    ops.analyze(1, 1.0e-9)
    return ops.criticalTimeStep()


@pytest.mark.parametrize("formulation", ["std", "uri"])
@pytest.mark.parametrize("mesh", ["cube", "dist"], ids=["cube", "distorted"])
def test_critical_timestep_matches_numpy_pencil(formulation, mesh):
    """ops.criticalTimeStep() (the ADR-65 per-element algebraic pencil, exact
    for any element) equals the independent numpy 60-DOF generalized-eig
    value to <=1e-10 rel — K per formulation, M = the HRZ lump in use.
    Measured rel ~3e-16 on all four combos."""
    corners = _CUBE if mesh == "cube" else _DIST
    d_ref = _oracle_dtcr(corners, formulation)
    d_ops = _ops_dtcr(corners, formulation)
    rel = abs(d_ops - d_ref) / d_ref
    assert rel <= 1e-10, (
        f"{formulation}/{mesh}: criticalTimeStep {d_ops} vs numpy pencil "
        f"{d_ref} (rel {rel:.2e})"
    )
    print(f"\n[P3 dt_cr] {formulation}/{mesh}: {d_ops:.8e} (pencil rel {rel:.1e})")


# ==========================================================================
# 4 + 5. explicit wave bar (CDL) — stability, dt-ratio report, energy closure
# ==========================================================================
def _run_wave_bar(efile, nsteps=2500, safety=0.9, v2=False):
    """10-element H20 wave bar, CentralDifferenceLadruno, -lumped, force-pulse
    excitation (Trig ramp; NO velocity ICs — the ADR-69 quirk), constant dt =
    safety * dt_cr (never changed mid-run — CDL dt changes need a
    revertToLastStep reseed). Returns (dt_cr, tip displacement history)."""
    m = _h20_bar(10, lumped=True)                 # E=1000, nu=0, rho=1
    xmax = max(c[0] for c in m.coords.values())
    tip = None
    for t, c in m.coords.items():
        if abs(c[0] - xmax) < 1e-9 and abs(c[1]) < 1e-9 and abs(c[2]) < 1e-9:
            tip = t
    _clamp_root(m, axial_only=True)

    ops.timeSeries("Trig", 1, 0.0, 0.04, 0.08, "-factor", 5.0)   # smooth pulse
    ops.pattern("Plain", 1, 1)
    ops.load(tip, 1.0, 0.0, 0.0)

    rec = ["EnergyBalance", "-file", efile, "-time"] + (["-v2"] if v2 else [])
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno")
    ops.recorder(*rec)                            # integrator BEFORE recorder (ADR-69)
    ops.analysis("Transient")

    ops.analyze(1, 1.0e-9)                        # prime; enables criticalTimeStep()
    dt_cr = ops.criticalTimeStep()
    assert dt_cr > 0.0 and math.isfinite(dt_cr)

    dt = safety * dt_cr
    disp = []
    for i in range(nsteps):
        assert ops.analyze(1, dt) == 0, f"explicit wave bar failed at step {i}"
        d = ops.nodeDisp(tip, 1)
        assert math.isfinite(d), f"non-finite tip displacement at step {i}"
        disp.append(d)
    return dt_cr, disp


def _h8_equalnode_dtcr():
    """Equal-node-spacing LadrunoBrick H8 bar (20 elements of h=0.5 — the H20
    bar's 0.5 nodal spacing along x), -lumped; returns its dt_cr."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.0, 1.0)
    m = _Mesh()
    for e in range(20):
        conn = [m.add(c) for c in _box(e * 0.5, (e + 1) * 0.5, 0, 1, 0, 1)]
        ops.element("LadrunoBrick", e + 1, *conn, 1, "-formulation", "std", "-lumped")
    _clamp_root(m, axial_only=True)
    tipn = max(m.coords, key=lambda t: m.coords[t][0])
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(tipn, 1.0, 0.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator("ExplicitBathe", 0.54, "-cfl")
    ops.analysis("Transient")
    ops.analyze(1, 1.0e-9)
    return ops.criticalTimeStep()


def test_explicit_wave_bar_stable_and_dt_ratio(tmp_path):
    """The wave bar runs 2500 steps at 0.9x the pencil dt: stable, finite,
    bounded. The dt_cr ratio vs the equal-node H8 mesh is REPORTED (loose
    sanity band only, per the ADR §6 gate wording): measured ~0.50 — below
    the 1-D 2/sqrt(6)~0.82 rod ballpark because the 3-D HRZ corner masses
    are far smaller than the rod-lumping shares (see module docstring)."""
    efile = str(tmp_path / "wave_energy.txt")
    dt20, disp = _run_wave_bar(efile)
    max_u = max(abs(d) for d in disp)
    assert max_u < 0.5, f"wave-bar response unbounded (max |u| = {max_u})"
    assert max_u > 1e-6, "wave bar did not move — vacuous run"

    dt8 = _h8_equalnode_dtcr()
    ratio = dt20 / dt8
    assert 0.2 < ratio < 1.0, (
        f"H20/H8 dt_cr ratio {ratio:.3f} outside the loose sanity band"
    )
    print(f"\n[P3 explicit] H20 dt_cr={dt20:.5e}  equal-node H8 dt_cr={dt8:.5e}  "
          f"ratio={ratio:.4f} (1-D rod theory 2/sqrt(6)={2 / math.sqrt(6):.4f})")


def test_energy_balance_closure_no_hourglass_channel(tmp_path, capfd):
    """ADR-69 energy closure on the wave bar (KE IE DW ULW RES ERR% columns):
    after the pulse ends the total energy drifts <1% of peak and |RES| stays
    <8% of peak (measured: 0.14% drift, 4.0% RES at 0.9x dt). The -v2
    channel-aware layout must NOT declare an E_hg (hourglass) column at all —
    this element has no hourglass machinery, so the channel does not exist
    (stronger than an all-zero column)."""
    efile = str(tmp_path / "wave_energy_v2.txt")
    capfd.readouterr()
    _run_wave_bar(efile, v2=True)
    echo = capfd.readouterr()
    text = echo.out + echo.err
    ops.wipe()                                    # flush/close the recorder

    assert "E_hg" not in text, (
        f"EnergyBalance -v2 declared a hourglass channel for a hourglass-free "
        f"element: {text[:400]!r}"
    )

    d = np.atleast_2d(np.loadtxt(efile))
    # column names from the echoed header — never parse by position
    import re
    mcols = re.findall(r"columns =( time)? \[model: ([^\]]+)\]", text)
    assert mcols, "no EnergyBalanceRecorder column echo captured"
    has_time, names = mcols[-1]
    cols = (["time"] if has_time else []) + names.split()

    t = d[:, cols.index("time")]
    # -v2 splits KE/DW into element/nodal buckets (KE_ele KE_nod ...)
    ke = d[:, cols.index("KE_ele")] + d[:, cols.index("KE_nod")]
    ie = d[:, cols.index("IE")]
    res = d[:, cols.index("RES")]
    tot = ke + ie
    peak = tot.max()
    assert peak > 0.0
    steady = t > 0.2                              # pulse (0.08 s) well over
    assert steady.any()
    res_rel = np.abs(res[steady]).max() / peak
    drift_rel = abs(tot[-1] - tot[steady][0]) / peak
    assert res_rel < 0.08, f"energy residual {res_rel:.1%} of peak (want <8%)"
    assert drift_rel < 0.01, f"total-energy drift {drift_rel:.2%} of peak (want <1%)"
    print(f"\n[P3 energy] peakE={peak:.4e}  steady max|RES|/peak={res_rel:.3%}  "
          f"drift={drift_rel:.3%}  (no E_hg channel)")


# ==========================================================================
# 6. betaK-Rayleigh clobber gate (ADR-66/#562 family bug), x4
# ==========================================================================
def _bar_rig(formulation, lumped):
    """3-element H20 cantilever bar for the clobber gate (small, fast)."""
    m = _h20_bar(3, lumped, e_mod=E, nu=0.0, rho=1.0, formulation=formulation)
    _clamp_root(m, axial_only=True)
    tip = max(m.coords, key=lambda t: m.coords[t][0])
    return m, tip


def _tip_static(formulation, lumped, load):
    m, tip = _bar_rig(formulation, lumped)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(tip, load, 0.0, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    return abs(ops.nodeDisp(tip, 1))


def _tip_dynamic_peak(formulation, lumped, load, betaK):
    """Step load at t=0 (Constant series), Newmark transient over ~1.6
    fundamental periods; returns the peak tip displacement."""
    m, tip = _bar_rig(formulation, lumped)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(tip, load, 0.0, 0.0)
    if betaK != 0.0:
        ops.rayleigh(0.0, betaK, 0.0, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")

    # fundamental axial period of the 3-long bar: T1 = 4L/c
    T1 = 4.0 * 3.0 / math.sqrt(1000.0 / 1.0)
    nsteps, dt = 120, 1.6 * T1 / 120
    peak = 0.0
    for _ in range(nsteps):
        assert ops.analyze(1, dt) == 0
        peak = max(peak, abs(ops.nodeDisp(tip, 1)))
    return peak


@pytest.mark.parametrize("lumped", [False, True], ids=["consistent", "lumped"])
@pytest.mark.parametrize("formulation", ["std", "uri"])
def test_betak_rayleigh_preserves_inertia(formulation, lumped):
    """getResistingForceIncInertia under rayleigh(0,betaK,0,0) must keep M*a
    and -Q (the ADR-66/#562 resid-clobber family bug). Self-validating: the
    undamped step response first PROVES the rig resolves dynamics (~2x static
    overshoot), then a tiny betaK must move the peak <5%. Under the bug the
    damped peak collapses toward 1x static (>40% shift)."""
    load = 1.0e-3
    static = _tip_static(formulation, lumped, load)
    peak0 = _tip_dynamic_peak(formulation, lumped, load, 0.0)
    # vacuity threshold 1.3x (not the SDOF 2x): the tip is a corner node of a
    # quadratic face, so the corner point-load carries a local static
    # flexibility share that the traveling-wave overshoot does not double —
    # measured undamped overshoot ~1.39x static, comfortably dynamic.
    assert peak0 > 1.3 * static, (
        f"[{formulation}/{'lumped' if lumped else 'consistent'}] undamped step "
        f"overshoot only {peak0 / static:.2f}x static — rig not resolving "
        "dynamics; gate would be vacuous"
    )
    peak_d = _tip_dynamic_peak(formulation, lumped, load, 1.0e-5)
    rel = abs(peak_d - peak0) / peak0
    assert rel < 0.05, (
        f"[{formulation}/{'lumped' if lumped else 'consistent'}] tiny betaK "
        f"moved the dynamic peak by {rel:.1%} ({peak0:.6g} -> {peak_d:.6g}, "
        f"static {static:.6g}) — inertia dropped in getResistingForceIncInertia"
    )


# ==========================================================================
# 7. massType serialization — the -lumped wire bit round-trips
# ==========================================================================
def test_masstype_survives_database_roundtrip():
    """recvSelf must restore massType=1: the round-tripped element's
    eigenvalues (mass-sensitive probe) match the pre-save lumped ones
    bit-for-bit, and differ from a consistent-mass build (non-vacuity)."""
    from _testbed.roundtrip import database_roundtrip

    nodes = _hex20_nodes(_DIST)

    def build(extra):
        def _b():
            ops.wipe()
            ops.model("basic", "-ndm", 3, "-ndf", 3)
            for i, c in enumerate(nodes):
                ops.node(i + 1, *c)
            ops.nDMaterial("ElasticIsotropic", 1, E, NU, RHO)
            for n in (1, 2, 3, 4, 9, 10, 11, 12):
                ops.fix(n, 1, 1, 1)
            ops.element("LadrunoBrick20", 1, *_CONN, 1, *extra)
            ops.timeSeries("Linear", 1)
            ops.pattern("Plain", 1, 1)
            ops.load(7, 0.0, 0.0, -1.0)
            ops.system("FullGeneral")
            ops.numberer("Plain")
            ops.constraints("Plain")
            ops.integrator("LoadControl", 1.0)
            ops.algorithm("Linear")
            ops.analysis("Static")
            assert ops.analyze(1) == 0
        return _b

    def eig_probe():
        return list(ops.eigen("-fullGenLapack", 6))

    # non-vacuity reference: consistent vs lumped eigenvalues differ
    build(["-lumped"])()
    eig_lumped = eig_probe()
    build([])()
    eig_cons = eig_probe()
    diff = max(abs(a - b) / abs(b) for a, b in zip(eig_lumped, eig_cons))
    assert diff > 1e-3, (
        "lumped and consistent eigenvalues identical — the probe cannot "
        f"detect a lost massType (rel diff {diff:.2e})"
    )

    # the round-trip itself: probe_fn compares eigenvalues before/after restore
    database_roundtrip(build(["-lumped"]), probe_nodes=[5, 6, 7, 8], ndf=3,
                       dbname="lb20_dyn_rt", probe_fn=eig_probe)
