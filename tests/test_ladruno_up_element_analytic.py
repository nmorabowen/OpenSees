"""LadrunoUP (ELE 33017) — ADR-71 P1 ANALYTIC Zone-B battery (WP1.F(i)).

Analytic / closed-form gates for the unified Biot u-p element under the
honest-p contract (ADR-71 §3.2): the nodal pressure DOF *is* p, blocks split
by the time-derivative they multiply —

    getTangentStiff  [ K , -Q ; 0 , H ]      getDamp  [ C_Ray , 0 ; Q^T , S+Ht ]
    getMass          [ M , 0  ; 0 , 0 ]      resisting force carries NO rate terms

The battery (one test per WP1.F(i) bullet):

  1. Terzaghi 1D consolidation vs the numpy-summed series (Q4 1x10 column)
     + settlement-vs-sqrt(t) behaviour;
  2. B2 ZS84 Example-1 column (ADR §7.1) at the ~1e-3 gate;
  3. B3 Boone-Ingraffea hard numbers u(0+)=0.254 mm / p(0+)=0.410 MPa
     (+ u(inf)=0.079 mm);
  4. static drained + steady-seepage NEW capability (does not exist upstream:
     the UP-ucsd/UWelements rate-DOF contract makes `Static` structurally
     singular — this element's H block lives in getTangentStiff, so a static
     seepage solve is well-posed);
  5. all-impervious static singularity smoke (structural singularity pinned
     by SVD + the solver-dependent arbitrary-p symptom; the ADR's "loud
     failure" expectation is xfail-documented — see the test docstrings);
  6. wrong-solver divergence (ProfileSPD silently drops one Q block,
     ADR §3.2 <FW-F1>) — documented + LEDGER_quirks row;
  7. Rayleigh contract: with betaK-only damping the p-rows/cols of the
     assembled effective tangent are EXACTLY c2*Q^T / H + c2*(S+Ht) — no
     betaK*(-Q)/betaK*H contamination (+ alphaM free-vibration decay sanity);
  8. FD consistency: static tangent vs central-difference residual, and the
     damp/mass blocks on the transient (IncInertia) path vs FD w.r.t. nodal
     vel/accel, cross-checked against an independent numpy block oracle.

MATRIX/RESIDUAL EXTRACTION MECHANISMS (established here, reusable):
  * assembled effective tangent: ``ops.printA('-ret')`` under
    system('FullGeneral') + numberer('Plain') + constraints('Plain') on a
    1-element model — NOTE the returned buffer is COLUMN-major (LAPACK
    storage): reshape then TRANSPOSE.  With constraints('Plain') and no fix,
    equation numbering == node-major DOF order, so element blocks are read
    off directly.
  * residual R(x): prescribe every DOF with ``ops.sp`` under
    constraints('Penalty') (NEVER 'Transformation' — a fully-prescribed rig
    leaves zero free equations and FullGenLinSOE dies, a known quirk), then
    ``ops.reactions()``; nodeReaction == assembled resisting force.
  * transient residual (getResistingForceIncInertia): after the static
    penalty solve, impose trial rates with ``ops.setNodeVel/-Accel
    (..., '-commit')`` and call ``ops.reactions('-dynamic')``.

DELIBERATE PINS (measured while building this battery):
  * '-dynSeepage off' in every accuracy-gated consolidation run: the default
    'on' evaluates the seepage drive with TRIAL accelerations, and the
    Newmark-recovered accel noise of numerically-damped wave modes pollutes
    the p-field WITHOUT BOUND as dt shrinks (ZS84 column, gamma=0.6: max rel
    err 1.8e-2 at dt=0.08 growing to 8.7e-1 at dt=0.005; with 'off' the same
    sweep sits at 1.3e-3..7e-4 and CONVERGES).  Reported upward — element
    source is not touched here.
  * gamma >= 0.6 Newmark everywhere transient: gamma=1/2 leaves the p-row
    parasitic root -(1-gamma)/gamma on the unit circle (ADR §3.2) — a
    gamma=0.5 consolidation march measurably DIVERGES (~1e10 by Tv=0.5).

Run:  py -3.12 -m pytest tests/test_ladruno_up_element_analytic.py -x -q
(the worktree dist/bin opensees.pyd is CPython 3.12 — plain `python` is 3.11
and will not import it).
"""
import os
import sys
from pathlib import Path

import numpy as np
import pytest

# --------------------------------------------------------------------------
# bootstrap: THIS worktree's engine (evict any boot-.pth preloaded stale pyd)
# --------------------------------------------------------------------------
_DIST = str(Path(__file__).resolve().parents[1] / "dist" / "bin")
if not os.path.isfile(os.path.join(_DIST, "opensees.pyd")):
    pytest.skip(f"worktree engine not built: {_DIST}", allow_module_level=True)

os.environ["PATH"] = _DIST + os.pathsep + os.environ.get("PATH", "")
try:
    os.add_dll_directory(_DIST)
except (FileNotFoundError, OSError):
    pass
if _DIST not in sys.path:
    sys.path.insert(0, _DIST)
# the installer's ladruno_opensees.pth boot module PRELOADS another build's
# opensees.pyd at interpreter start — evict it so THIS worktree's pyd wins.
for _m in ("opensees", "openseespy", "openseespy.opensees"):
    sys.modules.pop(_m, None)
import opensees as ops  # noqa: E402

assert os.path.normcase(os.path.dirname(ops.__file__)) == os.path.normcase(_DIST), (
    f"wrong opensees.pyd imported: {ops.__file__} (want {_DIST})"
)

pytestmark = [pytest.mark.zone_b]

_GENERAL = ("UmfPack",)   # the honest-p tangent is unsymmetric (ADR §3.2)


# --------------------------------------------------------------------------
# numpy oracles
# --------------------------------------------------------------------------
def terzaghi_series(Z, Tv, nterms=300):
    """Excess-pore-pressure ratio p/p0 of the Terzaghi series.

    Z = z/H measured from the DRAINED face; single-drainage path H."""
    m = np.arange(nterms)
    M = (2 * m + 1) * np.pi / 2.0
    return np.sum((2.0 / M) * np.sin(M * Z) * np.exp(-M * M * Tv))


def terzaghi_U(Tv, nterms=300):
    """Average degree of consolidation U(Tv) of the Terzaghi series."""
    m = np.arange(nterms)
    M = (2 * m + 1) * np.pi / 2.0
    return 1.0 - np.sum((2.0 / M**2) * np.exp(-M * M * Tv))


def q4_unit_square_blocks(kbar, one_over_Qbar, biot_alpha=1.0, rho=None,
                          rhoF=None):
    """Independent numpy oracle: Q (8x4), H, S, Ht-unit (4x4), (optionally)
    the consistent solid mass M (8x8) and the dynamic-seepage acceleration
    coupling G[b, a*2+i] = int dNp_{b,i} kbar_i rhoF Nu_a dV (4x8) of the
    UNIT SQUARE Q4, 2x2 Gauss, thickness 1 — the same integrals
    LadrunoUPKernel.h assembles."""
    XY = np.array([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])
    g = 1.0 / np.sqrt(3.0)
    Q = np.zeros((8, 4))
    H = np.zeros((4, 4))
    S = np.zeros((4, 4))
    Ht = np.zeros((4, 4))
    M = np.zeros((8, 8))
    G = np.zeros((4, 8))
    for xi, eta in [(-g, -g), (g, -g), (g, g), (-g, g)]:
        dN = np.array([[-(1 - eta), -(1 - xi)],
                       [(1 - eta), -(1 + xi)],
                       [(1 + eta), (1 + xi)],
                       [-(1 + eta), (1 - xi)]]) / 4.0
        N = np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                      (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)]) / 4.0
        J = dN.T @ XY
        detJ = np.linalg.det(J)
        dNx = dN @ np.linalg.inv(J).T
        w = detJ  # unit GP weight product folded (2x2 rule weights = 1)
        for a in range(4):
            for i in range(2):
                Q[2 * a + i, :] += dNx[a, i] * biot_alpha * N * w
        H += (dNx * kbar) @ dNx.T * w      # diagonal kbar
        S += np.outer(N, N) * one_over_Qbar * w
        Ht += dNx @ dNx.T * w              # unit-alpha stabilization Laplacian
        if rho is not None:
            for i in range(2):
                M[i::2, i::2] += rho * np.outer(N, N) * w
        if rhoF is not None:
            for b in range(4):
                for a in range(4):
                    for i in range(2):
                        G[b, 2 * a + i] += dNx[b, i] * kbar * rhoF * N[a] * w
    return Q, H, S, Ht, M, G


# element-node-major DOF indexing of the 12x12 single-Q4 system [ux,uy,p]/node
_UIDX = [0, 1, 3, 4, 6, 7, 9, 10]
_PIDX = [2, 5, 8, 11]


# --------------------------------------------------------------------------
# model builders
# --------------------------------------------------------------------------
def _column(nely, L, E, nu, rho, Kf, poro, rhoF, kbar, extra=(),
            drained_top=True, fix_p_top_val=None):
    """2-node-wide Q4 consolidation column: ux=0 everywhere, uy=0 at base,
    p=0 (or prescribed) on the top row when drained_top. Returns nid map."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    nid = {}
    k = 1
    for j in range(nely + 1):
        for i in range(2):
            ops.node(k, float(i), L * j / nely)
            nid[(i, j)] = k
            k += 1
    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)
    for j in range(nely):
        ops.element("LadrunoUP", j + 1,
                    nid[(0, j)], nid[(1, j)], nid[(1, j + 1)], nid[(0, j + 1)],
                    1, "-Kf", Kf, "-poro", poro, "-rhoF", rhoF,
                    "-perm", kbar, kbar, *extra)
    for j in range(nely + 1):
        for i in range(2):
            ops.fix(nid[(i, j)], 1, (1 if j == 0 else 0), 0)
    if drained_top:
        for i in range(2):
            ops.fix(nid[(i, nely)], 0, 0, 1)
    return nid


def _analysis_transient(gamma=0.6, beta=0.3025):
    ops.system(*_GENERAL)
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.algorithm("Linear")
    ops.integrator("Newmark", gamma, beta)
    ops.analysis("Transient")


# ==========================================================================
# 1. Terzaghi 1D consolidation vs series
# ==========================================================================
@pytest.mark.t1
def test_terzaghi_1d_consolidation_vs_series():
    """Q4 1x10 column, drained top / sealed sides+bottom, instantaneous load.

    p(z,t) vs the numpy Terzaghi series at Tv in {0.05, 0.2, 0.5} and
    settlement-vs-sqrt(t) via the degree of consolidation U(Tv).

    TOLERANCE RATIONALE (discretization-aware, mesh h = H/10): the spatial
    error of the bilinear pair is O(h^2) once the pressure front spans a few
    elements — measured 0.59% / 0.16% / 0.03% of p0 at Tv = 0.1 / 0.2 / 0.5
    (gated 2% / 1% / 0.5%, ~3x margin) — but at Tv=0.05 the drainage
    boundary layer sqrt(Tv)*H = 2.2 m spans only ~2 elements and the damped
    step-load transient still rings, so Tv=0.05 is gated loose at 15%
    (measured 11.0%; a mesh-refinement leg is B2's job, not this one).
    gamma=0.6/beta=0.3025 (the production rider), dt = 0.02 s (well above
    the equal-order small-dt oscillation floor for this cv).
    """
    E, nu, rho = 1.0e4, 0.2, 2.0          # kPa, Mg/m3
    Kf, poro, rhoF, kbar = 4444.0, 0.4, 1.0, 1.0e-4
    L = 10.0
    q = 10.0
    nely = 10
    Eoed = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu))
    ooQ = poro / Kf
    Qbar = 1.0 / ooQ
    cv = kbar / (1.0 / Eoed + ooQ)
    p0 = q * Qbar / (Qbar + Eoed)          # undrained (compressible-fluid) p
    s0 = q * L / (Eoed + Qbar)             # undrained settlement
    sinf = q * L / Eoed                    # drained settlement

    nid = _column(nely, L, E, nu, rho, Kf, poro, rhoF, kbar,
                  extra=("-stab", "off", "-dynSeepage", "off"))
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for i in range(2):
        ops.load(nid[(i, nely)], 0.0, -q / 2.0, 0.0)
    _analysis_transient()

    assert ops.analyze(1, 1.0e-3) == 0
    t = 1.0e-3
    dt = 0.02
    H = L
    settle = {}
    perr = {}
    for Tv in (0.05, 0.1, 0.2, 0.5):
        tt = Tv * H * H / cv
        nst = int(round((tt - t) / dt))
        assert ops.analyze(nst, dt) == 0
        t += nst * dt
        TvA = t * cv / (H * H)
        errs = [abs(ops.nodeDisp(nid[(0, j)], 3)
                    - p0 * terzaghi_series((L - L * j / nely) / H, TvA)) / p0
                for j in range(1, nely)]
        perr[Tv] = max(errs)
        settle[Tv] = -ops.nodeDisp(nid[(0, nely)], 2)

    assert perr[0.05] < 0.15, f"Tv=0.05 max rel err {perr[0.05]:.3%}"
    assert perr[0.1] < 0.02, f"Tv=0.1 max rel err {perr[0.1]:.3%}"
    assert perr[0.2] < 0.01, f"Tv=0.2 max rel err {perr[0.2]:.3%}"
    assert perr[0.5] < 0.005, f"Tv=0.5 max rel err {perr[0.5]:.3%}"

    # settlement vs sqrt(t): U_fem(Tv) tracks the series U(Tv) (itself
    # ~sqrt(4Tv/pi) early), and the early-time sqrt(t) proportionality holds
    # (Tv=0.05 excluded: the damped step-load transient still contaminates
    # the settlement there — measured dU = +0.07)
    U_fem = {Tv: (settle[Tv] - s0) / (sinf - s0) for Tv in settle}
    for Tv in (0.1, 0.2, 0.5):
        assert abs(U_fem[Tv] - terzaghi_U(Tv)) < 0.02, (
            f"U(Tv={Tv}) fem {U_fem[Tv]:.4f} vs series {terzaghi_U(Tv):.4f}")
    ratio = U_fem[0.2] / U_fem[0.1]
    assert abs(ratio - np.sqrt(2.0)) < 0.07, (
        f"early-time settlement is not sqrt(t): U(0.2)/U(0.1) = {ratio:.3f}")


# ==========================================================================
# 2. B2 — ZS84 Example-1 column (ADR §7.1)
# ==========================================================================
@pytest.mark.t1
def test_b2_zs84_column_consolidation_gate():
    """Zienkiewicz-Shiomi 1984 Example-1 soil column (ADR §7.1 B2 pins):
    30 m column, q = 1 kN/m^2 ramped by sin(pi*t/2tp), tp = 0.1 s;
    E = 30 MN/m^2, nu = 0.2, rho = 2, rho_f = 1 Mg/m^3, K_f = 100 MN/m^2,
    K_s = inf, k' = 1e-2 m/s (kbar = k'/(rho_f*g), g = 9.81), n = 0.3;
    top drained, base impermeable/fixed.

    Gate: p(z,t) vs the compressible-storage consolidation series at the
    PINNED discretization nely = 30 (h = 1 m), dt = 0.02 s, gamma = 0.6 /
    beta = 0.3025, '-dynSeepage off' (see module docstring):
      Tv = 0.5 : max rel err (of q) <= 2e-3  (measured 1.05e-3 — the ADR's
                 ~1e-3 class; the dt-sweep 0.08..0.0025 is flat at ~1e-3 =
                 the O(h^2) spatial floor);
      Tv = 0.2 : <= 5e-3 (measured 3.55e-3) — the excess over the spatial
                 floor is the TAIL OF THE PHYSICAL compressible-wave ringing
                 (ADR §3.6 artifact (a)): it does NOT shrink with mesh
                 refinement (3.38e-3 at nely=60) and GROWS with smaller dt
                 (2.1e-2 at dt=0.01, where gamma=0.6 damps resolved wave
                 modes more slowly per unit time), exactly the ZS84 Fig 2/5/7
                 oscillation behaviour the paper itself reports.
    """
    L, E, nu, rho, rhoF = 30.0, 3.0e4, 0.2, 2.0, 1.0
    Kf, n, g = 1.0e5, 0.3, 9.81
    kbar = 1.0e-2 / (rhoF * g)
    q, tp = 1.0, 0.1
    nely, dt = 30, 0.02

    Eoed = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu))
    ooQ = n / Kf
    cv = kbar / (1.0 / Eoed + ooQ)
    p0 = q / (1.0 + ooQ * Eoed)

    nid = _column(nely, L, E, nu, rho, Kf, n, rhoF, kbar,
                  extra=("-stab", "off", "-dynSeepage", "off"))
    ts = list(np.linspace(0.0, tp, 21))
    vals = [float(np.sin(np.pi * t / (2 * tp))) for t in ts]
    ops.timeSeries("Path", 1, "-time", *(ts + [1.0e12]), "-values", *(vals + [1.0]))
    ops.pattern("Plain", 1, 1)
    for i in range(2):
        ops.load(nid[(i, nely)], 0.0, -q / 2.0, 0.0)
    _analysis_transient()

    t = 0.0
    gates = {0.2: 5.0e-3, 0.5: 2.0e-3}
    for Tv in (0.2, 0.5):
        tt = Tv * L * L / cv
        nst = int(round((tt - t) / dt))
        assert ops.analyze(nst, dt) == 0
        t += nst * dt
        TvA = t * cv / (L * L)
        err = max(abs(ops.nodeDisp(nid[(0, j)], 3)
                      - p0 * terzaghi_series((L - L * j / nely) / L, TvA)) / q
                  for j in range(1, nely))
        assert err < gates[Tv], (
            f"ZS84 Tv={Tv}: max rel err {err:.2e} (gate {gates[Tv]:.0e})")


# ==========================================================================
# 3. B3 — Boone-Ingraffea hard numbers
# ==========================================================================
@pytest.mark.t1
def test_b3_boone_ingraffea_hard_numbers():
    """B3 (ADR §7.1, Book §6.4.1): 6 m poroelastic column, instantaneous
    q = 1 MPa (+ p = 1 MPa at the drained top for the long-time limit);
    G = 6000, nu = 0.2 (E = 14400), K_s = 36000, K_f = 3000 MPa, n = 0.19,
    alpha = 0.79, k = 2e-5 m^2/(MPa*s). The ONLY gate exercising alpha != 1.

    Pinned targets: u(0+) = 0.254 mm, p(0+) = 0.410 MPa, u(inf) = 0.079 mm.

    Protocol: the undrained response is captured on the SEALED column (no
    p BC anywhere — exactly undrained, no drainage boundary layer) with a
    quasi-static ramp (rho = 2e-3 in MPa/Mg-consistent units => wave transit
    1.7e-3 s << tramp = 0.05 s << consolidation time 257 s); then the top
    p = 1 MPa head + drainage is switched on and marched to steady state.

    Tolerance 2%: the pinned numbers derive from nu_u = 0.33 (rounded) —
    first-principles alpha/Q̄ algebra gives u0 = 0.2521 mm / p0 = 0.4149 MPa,
    i.e. the pins themselves sit ~1% from the exact parameter set.
    """
    G, nu = 6000.0, 0.2
    E = 2 * G * (1 + nu)
    alpha, Kf, Ks, n = 0.79, 3000.0, 36000.0, 0.19
    k = 2.0e-5
    L, sigma = 6.0, 1.0
    rho = 2.0e-3          # 2 Mg/m^3 in MPa-consistent mass units
    nely = 12

    nid = _column(nely, L, E, nu, rho, Kf, n, 1.0e-3, k,
                  extra=("-alpha", alpha, "-Ks", Ks, "-stab", "off",
                         "-dynSeepage", "off"),
                  drained_top=False)      # SEALED: exactly undrained
    tramp = 0.05
    ts = [0.0, tramp, 1.0e12]
    ops.timeSeries("Path", 1, "-time", *ts, "-values", 0.0, 1.0, 1.0)
    ops.pattern("Plain", 1, 1)
    for i in range(2):
        ops.load(nid[(i, nely)], 0.0, -sigma / 2.0, 0.0)
    _analysis_transient()
    nramp = 25
    assert ops.analyze(nramp, tramp / nramp) == 0

    u0 = -ops.nodeDisp(nid[(0, nely)], 2) * 1000.0   # mm
    p0 = ops.nodeDisp(nid[(0, nely // 2)], 3)        # MPa, mid-column
    assert abs(u0 - 0.254) / 0.254 < 0.02, f"u(0+) = {u0:.4f} mm (want 0.254)"
    assert abs(p0 - 0.410) / 0.410 < 0.02, f"p(0+) = {p0:.4f} MPa (want 0.410)"

    # long-time limit: open the drained top with p = 1 MPa and march.
    # Adding SPs mid-analysis changes the equation count — UmfPack's
    # symbolic setSize fails (-8) on the stale handle, so rebuild the
    # analysis objects (standard staged-analysis wipeAnalysis idiom).
    # MEASURED TRAP (reported upward, not fixed here): a NONZERO p-DOF sp
    # STAGED onto a running transient under constraints('Transformation')
    # converges to a WRONG steady state (interior p uniform at the wrong
    # level, H*p != 0 at the constrained boundary — minimal repro: 4-element
    # sealed column, load 25 steps, then sp p=1 at top + wipeAnalysis:
    # Transformation -> p_int = -73.2 with p_top = 1; Penalty/Lagrange ->
    # the correct p = 1 everywhere; the same sp applied from step 1 is fine
    # under all three). Phase B therefore runs under Penalty.
    ops.timeSeries("Constant", 2)
    ops.pattern("Plain", 2, 2)
    for i in range(2):
        ops.sp(nid[(i, nely)], 3, 1.0)
    ops.wipeAnalysis()
    ops.system(*_GENERAL)
    ops.numberer("RCM")
    ops.constraints("Penalty", 1e12, 1e12)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")
    t, cur = tramp, 0.05
    while t < 2.5e3:                        # ~10x the consolidation time
        assert ops.analyze(1, cur) == 0
        t += cur
        cur = min(cur * 1.4, 25.0)
    uinf = -ops.nodeDisp(nid[(0, nely)], 2) * 1000.0
    assert abs(uinf - 0.079) / 0.079 < 0.02, f"u(inf) = {uinf:.4f} mm (want 0.079)"


# ==========================================================================
# 4. static steady-seepage NEW capability
# ==========================================================================
@pytest.mark.t1
def test_static_steady_seepage_new_capability():
    """STATIC steady seepage under a prescribed head difference — a genuinely
    NEW capability: it does not exist upstream. Every upstream u-p element
    (UP-ucsd, UWelements) parks H in getDamp under the rate-DOF trick, so a
    `Static` integrator sees structurally-zero pressure rows and the solve is
    singular; LadrunoUP's honest-p contract puts H in getTangentStiff and the
    drained-solid + steady-seepage problem becomes well-posed (ADR §3.2).

    Gate: linear p(z) profile through interior nodes (exact for the P1 basis)
    and eleResponse 'flux' = kbar*dh/L per GP, plus 'porePressure' GP
    consistency — all to 1e-9 (this discretization is exact for linear p).
    """
    L, nely = 4.0, 4
    kbar, dP = 1.0e-4, 12.0
    nid = _column(nely, L, 1.0e4, 0.2, 2.0, 5000.0, 0.4, 1.0, kbar,
                  extra=("-stab", "off"), drained_top=False)
    # pure seepage rig: fix ALL u, prescribe head at both ends
    # (_column already fixed ux everywhere + uy at the base — add the rest)
    for j in range(1, nely + 1):
        for i in range(2):
            ops.fix(nid[(i, j)], 0, 1, 0)
    for i in range(2):
        ops.fix(nid[(i, 0)], 0, 0, 1)          # base p = 0
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for i in range(2):
        ops.sp(nid[(i, nely)], 3, dP)          # top p = dP
    ops.system(*_GENERAL)
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0                  # a STATIC u-p solve succeeds

    for j in range(nely + 1):
        p_expect = dP * j / nely
        for i in range(2):
            assert abs(ops.nodeDisp(nid[(i, j)], 3) - p_expect) < 1e-9 * dP

    q_expect = kbar * dP / L                    # Darcy |q_y| = kbar*dh/L
    for e in range(1, nely + 1):
        flux = np.array(ops.eleResponse(e, "flux")).reshape(-1, 2)
        assert np.max(np.abs(flux[:, 0])) < 1e-12          # no lateral flow
        assert np.max(np.abs(flux[:, 1] + q_expect)) < 1e-9 * q_expect, (
            f"ele {e} flux {flux[:,1]} vs -kbar*dh/L = {-q_expect}")
        # GP pore pressure consistency (interpolated nodal p)
        gp_p = np.array(ops.eleResponse(e, "porePressure"))
        lo, hi = dP * (e - 1) / nely, dP * e / nely
        assert np.all(gp_p > lo - 1e-9) and np.all(gp_p < hi + 1e-9)


# ==========================================================================
# 5. all-impervious static — singularity smoke
# ==========================================================================
def _sealed_static(system):
    """4-element sealed column (k > 0, NO p fixed anywhere), static solve.
    Returns (rc, p-profile)."""
    nid = _column(4, 4.0, 1.0e4, 0.2, 2.0, 4444.0, 0.4, 1.0, 1.0e-4,
                  extra=("-stab", "off"), drained_top=False)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for i in range(2):
        ops.load(nid[(i, 4)], 0.0, -5.0, 0.0)
    ops.system(*system)
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    rc = ops.analyze(1)
    return rc, np.array([ops.nodeDisp(nid[(0, j)], 3) for j in range(5)])


@pytest.mark.t0
def test_all_impervious_static_tangent_is_singular():
    """A fully sealed static model (>= 1 drained node per region VIOLATED,
    ADR §3.2 <UP-6/FW-F7>) assembles a structurally singular tangent: H is a
    pure-Neumann Laplacian (H*1 = 0) and det = det(K)*det(H) = 0.  Pinned by
    SVD of the assembled A (sigma_min/sigma_max < 1e-12) AND by the practical
    symptom: the solved pressure LEVEL is solver-dependent garbage (each
    factorization amplifies whatever round-off lands on the zero pivot)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    for i, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], 1):
        ops.node(i, float(x), float(y))
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e4, 0.2, 2.0)
    ops.element("LadrunoUP", 1, 1, 2, 3, 4, 1, "-Kf", 5000.0, "-poro", 0.4,
                "-rhoF", 1.0, "-perm", 1e-4, 1e-4, "-stab", "off")
    ops.fix(1, 1, 1, 0)
    ops.fix(2, 1, 1, 0)                    # u restrained; p sealed everywhere
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 1.0, -2.0, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    ops.analyze(1)
    A = np.array(ops.printA("-ret"))
    N = int(round(len(A) ** 0.5))
    A = A.reshape(N, N).T                  # printA('-ret') is column-major
    sv = np.linalg.svd(A, compute_uv=False)
    assert sv[-1] / sv[0] < 1e-12, (
        f"sealed static tangent should be singular; sigma ratio {sv[-1]/sv[0]:.2e}")

    # practical symptom: the p level is an arbitrary solver-dependent constant
    rc_u, p_u = _sealed_static(("UmfPack",))
    rc_f, p_f = _sealed_static(("FullGeneral",))
    assert np.max(np.abs(p_u - p_f)) > 1.0, (
        "expected solver-dependent arbitrary pressure levels on the singular "
        f"sealed model; got UmfPack {p_u} vs FullGeneral {p_f}")


@pytest.mark.t0
@pytest.mark.xfail(strict=True, reason=(
    "ADR-71 §3.2/§7 expects the all-impervious static to fail LOUDLY "
    "(analyze rc != 0). MEASURED 2026-07-11: every serial general solver "
    "(UmfPack/FullGeneral/BandGeneral/SparseGeneral) factorizes the singular "
    "system through round-off and returns rc = 0 with an arbitrary, "
    "solver-dependent pressure level (p-RHS is consistent, so no solver "
    "reports the rank deficiency). The singularity itself is pinned by the "
    "companion test; this xfail documents the refuted loud-failure claim "
    "for the ADR/guide to reword."))
def test_all_impervious_static_fails_loudly():
    rc, _ = _sealed_static(("UmfPack",))
    assert rc != 0, "sealed static returned rc == 0 (silent, not loud)"


# ==========================================================================
# 6. wrong-solver divergence (ProfileSPD vs UmfPack)
# ==========================================================================
@pytest.mark.t0
def test_wrong_solver_divergence_profilespd_vs_umfpack():
    """DOCUMENTED SILENT-WRONG-ANSWER TRAP (ADR §3.2 <FW-F1>, LEDGER_quirks):
    the honest-p effective tangent is unsymmetric (-Q in getTangentStiff,
    +Q^T in getDamp). ProfileSPD — the no-`system`-command DEFAULT — stores
    only the upper-triangle-in-profile, silently discarding one coupling
    block; analyze() returns 0 every step and the pore pressures are garbage.

    Gate: identical Terzaghi columns under ProfileSPD vs UmfPack; ProfileSPD
    must report SUCCESS (that is the trap) while its pressures differ from
    UmfPack's by orders of magnitude (measured ~1e88 vs ~5 kPa)."""
    def run(system):
        nid = _column(10, 10.0, 1.0e4, 0.2, 2.0, 4444.0, 0.4, 1.0, 1.0e-4)
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        for i in range(2):
            ops.load(nid[(i, 10)], 0.0, -5.0, 0.0)
        ops.system(*system)
        ops.numberer("RCM")
        ops.constraints("Transformation")
        ops.algorithm("Linear")
        ops.integrator("Newmark", 0.6, 0.3025)
        ops.analysis("Transient")
        rcs = [ops.analyze(1, 0.1) for _ in range(60)]
        ps = np.array([ops.nodeDisp(nid[(0, j)], 3) for j in range(11)])
        return rcs, ps

    rc_u, p_u = run(("UmfPack",))
    rc_p, p_p = run(("ProfileSPD",))
    assert all(r == 0 for r in rc_u)
    assert all(r == 0 for r in rc_p), (
        "ProfileSPD reported failure — the trap documented here is that it "
        "does NOT; update the quirks row if the framework starts failing loudly")
    scale = np.max(np.abs(p_u))
    assert np.isfinite(scale) and scale > 0
    diff = np.max(np.abs(p_u - p_p)) / scale
    assert diff > 1.0e3, (
        f"expected materially wrong ProfileSPD pressures (measured ~1e88); "
        f"relative divergence only {diff:.2e}")


# ==========================================================================
# 7. Rayleigh contract
# ==========================================================================
def _effective_tangent(betaK=0.0, dt=0.1, gamma=0.6, beta=0.3025,
                       transient=True, stab_off=True):
    """Assembled effective tangent of the free-floating unit-square Q4 via
    printA under Plain numbering (equation == node-major DOF order)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    for i, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], 1):
        ops.node(i, float(x), float(y))
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e4, 0.2, 2.0)
    extra = ("-stab", "off") if stab_off else ()
    ops.element("LadrunoUP", 1, 1, 2, 3, 4, 1, "-Kf", 5000.0, "-poro", 0.4,
                "-rhoF", 1.0, "-perm", 1e-4, 1e-4, *extra)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 0.0, 0.0, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    if transient:
        ops.integrator("Newmark", gamma, beta)
        if betaK != 0.0:
            ops.rayleigh(0.0, betaK, 0.0, 0.0)
        ops.algorithm("Linear")
        ops.analysis("Transient")
        assert ops.analyze(1, dt) == 0
    else:
        ops.integrator("LoadControl", 1.0)
        ops.algorithm("Linear")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
    A = np.array(ops.printA("-ret"))
    N = int(round(len(A) ** 0.5))
    return A.reshape(N, N).T               # column-major -> true row-major


@pytest.mark.t0
def test_rayleigh_contract_beta_k_never_touches_p_blocks():
    """ADR §3.2 obligation 2 contract test: with betaK-only damping, the
    p-rows/cols of the assembled damp must be EXACTLY Q^T / S (+Ht) — no
    betaK*(-Q) / betaK*H contamination (the base-class helpers would compose
    betaK*getTangentStiff(), smearing -Q and H into the damping).

    On the 1-element model under Newmark(gamma,beta,dt) the assembled
    effective tangent is A = K_tot + c2*C + c3*M (c2 = gamma/(beta*dt),
    c3 = 1/(beta*dt^2)), so the check is EXACT against the independent numpy
    block oracle:
      A[p,u] = c2*Q^T,  A[u,p] = -Q,  A[p,p] = H + c2*S  — identical for
      betaK = 0 and betaK = 0.05;  A(0.05) - A(0) nonzero ONLY in the u-u
      block (= c2*betaK*K_solid)."""
    kbar, ooQ = 1e-4, 0.4 / 5000.0
    Q, H, S, Ht, _, _ = q4_unit_square_blocks(kbar, ooQ)
    dt, gamma, beta = 0.1, 0.6, 0.3025
    c2 = gamma / (beta * dt)

    A0 = _effective_tangent(betaK=0.0, dt=dt)
    A1 = _effective_tangent(betaK=0.05, dt=dt)

    scal = np.max(np.abs(A0))
    for A, tag in ((A0, "betaK=0"), (A1, "betaK=0.05")):
        pu = A[np.ix_(_PIDX, _UIDX)]
        up = A[np.ix_(_UIDX, _PIDX)]
        pp = A[np.ix_(_PIDX, _PIDX)]
        assert np.max(np.abs(pu - c2 * Q.T)) < 1e-9 * scal, f"{tag}: pu != c2*Q^T"
        assert np.max(np.abs(up + Q)) < 1e-9 * scal, f"{tag}: up != -Q"
        assert np.max(np.abs(pp - (H + c2 * S))) < 1e-9 * scal, f"{tag}: pp != H+c2*S"

    D = A1 - A0                             # only c2*betaK*K_solid may differ
    assert np.max(np.abs(D[np.ix_(_PIDX, _UIDX)])) == 0.0
    assert np.max(np.abs(D[np.ix_(_UIDX, _PIDX)])) == 0.0
    assert np.max(np.abs(D[np.ix_(_PIDX, _PIDX)])) == 0.0
    assert np.max(np.abs(D[np.ix_(_UIDX, _UIDX)])) > 1e-6 * scal, (
        "betaK had no effect on the u-u block — Rayleigh not wired?")


@pytest.mark.t0
def test_rayleigh_alpha_m_free_vibration_decays():
    """Sanity: alphaM-only Rayleigh damps a drained free vibration. Two
    identical runs (gamma = 0.6 algorithmic damping in BOTH) differing only
    by rayleigh(alphaM = 0.8): the damped run's late-time envelope must fall
    well below the undamped run's."""
    def run(alphaM):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        for i, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], 1):
            ops.node(i, float(x), float(y))
        ops.nDMaterial("ElasticIsotropic", 1, 1.0e4, 0.2, 2.0)
        ops.element("LadrunoUP", 1, 1, 2, 3, 4, 1, "-Kf", 5000.0, "-poro", 0.4,
                    "-rhoF", 1.0, "-perm", 1e-4, 1e-4, "-stab", "off")
        ops.fix(1, 1, 1, 1)
        ops.fix(2, 1, 1, 1)
        ops.fix(3, 0, 0, 1)                # drained: p pinned everywhere
        ops.fix(4, 0, 0, 1)
        # hold a preload for 25 steps (resolved by the march), then release
        # into free vibration (fundamental period ~0.06 s, dt = 2e-3)
        ops.timeSeries("Path", 1, "-time", 0.0, 0.05, 0.0501, 1.0e12,
                       "-values", 1.0, 1.0, 0.0, 0.0)
        ops.pattern("Plain", 1, 1)
        ops.load(3, 2.0, 0.0, 0.0)
        ops.load(4, 2.0, 0.0, 0.0)
        ops.system(*_GENERAL)
        ops.numberer("Plain")
        ops.constraints("Transformation")
        if alphaM != 0.0:
            ops.rayleigh(alphaM, 0.0, 0.0, 0.0)
        ops.algorithm("Linear")
        ops.integrator("Newmark", 0.6, 0.3025)
        ops.analysis("Transient")
        dt = 2.0e-3
        u = []
        for _ in range(425):
            assert ops.analyze(1, dt) == 0
            u.append(ops.nodeDisp(3, 1))
        u = np.array(u[25:])                 # free-vibration window
        early = np.max(np.abs(u[:100]))
        late = np.max(np.abs(u[300:]))
        assert early > 1e-8, "no free vibration excited"
        return late / early

    r_damped = run(20.0)                     # zeta = alphaM/(2w) ~ 0.1
    r_free = run(0.0)
    assert r_damped < 0.5 * r_free, (
        f"alphaM damping ineffective: late/early envelope {r_damped:.4f} "
        f"(damped) vs {r_free:.4f} (undamped)")


# ==========================================================================
# 8. FD consistency (static tangent + IncInertia damp/mass blocks)
# ==========================================================================
_FD_NODES = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
_FD_MAT = ("ElasticIsotropic", 1, 1.0e4, 0.2, 2.0)
_FD_ELE = ("LadrunoUP", 1, 1, 2, 3, 4, 1, "-Kf", 5000.0, "-poro", 0.4,
           "-rhoF", 1.0, "-perm", 1e-4, 1e-4, "-stab", "off")


def _fd_model_prescribed(xvec):
    """Build the 1-element model with EVERY DOF penalty-prescribed to xvec
    and solve one static step (constraints('Penalty') — Transformation on a
    fully-prescribed rig leaves zero free equations and FullGenLinSOE dies)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    for i, (x, y) in enumerate(_FD_NODES, 1):
        ops.node(i, x, y)
    ops.nDMaterial(*_FD_MAT)
    ops.element(*_FD_ELE)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    k = 0
    for nd in (1, 2, 3, 4):
        for d in (1, 2, 3):
            ops.sp(nd, d, float(xvec[k]))
            k += 1
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Penalty", 1e12, 1e12)
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0


def _gather_reactions():
    R = np.zeros(12)
    k = 0
    for nd in (1, 2, 3, 4):
        rr = ops.nodeReaction(nd)
        for d in range(3):
            R[k] = rr[d]
            k += 1
    return R


def _static_resid(xvec):
    _fd_model_prescribed(xvec)
    ops.reactions()
    return _gather_reactions()


def _dyn_resid(xvec, vvec, avec):
    """getResistingForceIncInertia residual: static state from xvec, then
    trial rates imposed and reactions('-dynamic') assembled."""
    _fd_model_prescribed(xvec)
    k = 0
    for nd in (1, 2, 3, 4):
        for d in (1, 2, 3):
            ops.setNodeVel(nd, d, float(vvec[k]), "-commit")
            ops.setNodeAccel(nd, d, float(avec[k]), "-commit")
            k += 1
    ops.reactions("-dynamic")
    return _gather_reactions()


_X0 = np.array([0.0, 0.0, 0.5, 0.01, -0.005, 2.0,
                0.008, 0.012, 3.0, -0.003, 0.009, 1.5])
_V0 = np.array([0.001, -0.002, 0.03, 0.002, 0.001, -0.02,
                -0.001, 0.003, 0.05, 0.0, -0.001, 0.01])
_A0 = np.array([0.01, -0.02, 0.3, 0.02, 0.01, -0.2,
                -0.01, 0.03, 0.5, 0.0, -0.01, 0.1])


@pytest.mark.t0
def test_fd_static_tangent_consistency():
    """Central-difference d(resisting force)/d(u,p) at a random-ish state vs
    the assembled static tangent [K, -Q; 0, H] from printA. Tolerance 1e-6
    relative (measured 9.4e-9 — linear elastic, so FD error is pure
    penalty/solver round-off)."""
    h = 1e-5
    Kfd = np.zeros((12, 12))
    for j in range(12):
        xp = _X0.copy()
        xp[j] += h
        xm = _X0.copy()
        xm[j] -= h
        Kfd[:, j] = (_static_resid(xp) - _static_resid(xm)) / (2 * h)
    Ka = _effective_tangent(transient=False)
    scal = np.max(np.abs(Ka))
    assert np.max(np.abs(Kfd - Ka)) < 1e-6 * scal, (
        f"FD static tangent mismatch: {np.max(np.abs(Kfd - Ka))/scal:.2e} rel")
    # and the honest-p block pattern itself
    Q, H, _, _, _, _ = q4_unit_square_blocks(1e-4, 0.4 / 5000.0)
    assert np.max(np.abs(Ka[np.ix_(_PIDX, _UIDX)])) == 0.0   # static pu = 0
    assert np.max(np.abs(Ka[np.ix_(_UIDX, _PIDX)] + Q)) < 1e-9 * scal
    assert np.max(np.abs(Ka[np.ix_(_PIDX, _PIDX)] - H)) < 1e-9 * scal


@pytest.mark.t0
def test_fd_inc_inertia_damp_and_mass_blocks():
    """FD of the TRANSIENT residual (getResistingForceIncInertia via
    reactions('-dynamic')) w.r.t. nodal velocities == getDamp, and w.r.t.
    accelerations == getMass PLUS the dynamic-seepage acceleration coupling —
    cross-checked against the independent numpy oracle:
      dR/dv = C = [0, 0; Q^T, S]     (no Rayleigh set, stab off)
      dR/da = [M_solid, 0; +G, 0],   G = int dNp^T kbar rhoF Nu dV
    The +G p-row block is the '-dynSeepage on' (default) RESIDUAL term of
    ADR §3.1 — deliberately kept out of the TANGENT (test 7 shows the
    assembled A[p,u] is exactly c2*Q^T with no c3*G) — so this FD leg both
    verifies M and MEASURES the omitted-tangent term the ADR §7 P0 row asks
    to keep quantified. Tolerance 1e-6 relative (measured ~4e-15/1e-9)."""
    kbar, ooQ, rho, rhoF = 1e-4, 0.4 / 5000.0, 2.0, 1.0
    Q, H, S, Ht, Msol, G = q4_unit_square_blocks(kbar, ooQ, rho=rho, rhoF=rhoF)
    h = 1e-5

    Cfd = np.zeros((12, 12))
    for j in range(12):
        vp = _V0.copy()
        vp[j] += h
        vm = _V0.copy()
        vm[j] -= h
        Cfd[:, j] = (_dyn_resid(_X0, vp, _A0) - _dyn_resid(_X0, vm, _A0)) / (2 * h)
    cs = max(np.max(np.abs(Q)), np.max(np.abs(S)))
    assert np.max(np.abs(Cfd[np.ix_(_UIDX, _UIDX)])) < 1e-6 * cs   # no Rayleigh
    assert np.max(np.abs(Cfd[np.ix_(_UIDX, _PIDX)])) < 1e-6 * cs   # damp up = 0
    assert np.max(np.abs(Cfd[np.ix_(_PIDX, _UIDX)] - Q.T)) < 1e-6 * cs
    assert np.max(np.abs(Cfd[np.ix_(_PIDX, _PIDX)] - S)) < 1e-6 * cs

    Mfd = np.zeros((12, 12))
    for j in range(12):
        ap = _A0.copy()
        ap[j] += h
        am = _A0.copy()
        am[j] -= h
        Mfd[:, j] = (_dyn_resid(_X0, _V0, ap) - _dyn_resid(_X0, _V0, am)) / (2 * h)
    ms = np.max(np.abs(Msol))
    assert np.max(np.abs(Mfd[np.ix_(_UIDX, _UIDX)] - Msol)) < 1e-6 * ms
    assert np.max(np.abs(Mfd[:, _PIDX])) < 1e-9 * ms               # p-cols zero
    # p-rows: NOT zero — they carry +G, the dynamic-seepage residual coupling
    gs = np.max(np.abs(G))
    assert gs > 0
    assert np.max(np.abs(Mfd[np.ix_(_PIDX, _UIDX)] - G)) < 1e-6 * gs, (
        "dynamic-seepage acceleration coupling in the IncInertia residual "
        "does not match the analytic block int dNp^T kbar rhoF Nu dV")
