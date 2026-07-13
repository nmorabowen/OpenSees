r"""LadrunoUP (ELE 33017) — ADR-71 P3 B1 ZCB80 phasor gate (WP3.C-B1, OPUS).

The Taylor-Hood BT6 lane (quadratic u on the 6-node Bernstein triangle, linear
vertex p, heterogeneous ndf) driven by the Zienkiewicz-Chang-Bettess 1980
periodic-load poroelastic column (ADR-71 §7.1 B1; Book Fig 5.25). OpenSees is
real-valued / time-domain, so the frequency-domain benchmark is realised as a
transient sinusoidal drive rung to steady state, with per-node phasors
reconstructed by a least-squares A·sin+B·cos fit over the last cycle and
compared to an INDEPENDENT numpy closed form.

===========================================================================
CLOSED FORM — derivation route + independent numpy validation
===========================================================================
ZCB80 eqs 12-30 are the frequency-domain transfer functions of a 1D u-p Biot
column. Rather than transcribe the paper's algebra, the closed form here is
DERIVED from the u-p governing equations directly (the ADR treats the two as
equivalent). 1D column, vertical z (base z=0, top z=L=30 m), harmonic ansatz
field(z,t) = Im{field_hat(z)·e^{iωt}} (the Im convention matches a sin drive):

  momentum (mixture, tension-positive σ'=D u', total σ=σ'-α p):
      D u'' + ρ ω² u - α p' = 0                                        (1)
  fluid mass balance (Darcy + storage + coupling; dynamic-seepage term
      k̄ρ_f ω² u' verified negligible below and DROPPED, matching the FE run
      under -dynSeepage off):
      -k̄ p'' + (iω/Q̄) p + iωα u' = 0                                  (2)

Modal ansatz [u,p]=[U,P]e^{λz} gives a quadratic in μ=λ²:
      k̄D μ² + [k̄ρω² - iωD/Q̄ - iωα²] μ - iω ρ ω² /Q̄ = 0
→ 4 roots λ=±√μ₁,±√μ₂; mode ratio r=P/U=(Dλ²+ρω²)/(αλ). Four boundary
conditions close the system:
      u(0)=0        (rigid base)         p'(0)=0     (impermeable base)
      p(L)=0        (drained top)        D u'(L)=-q  (total-stress load, q=100)
The two waves are wildly disparate (a ~5 cm consolidation boundary layer, |λ|L
~ 600, and a ~km elastic wave), so the modal basis is written with a per-mode
reference point exp(λ(z-z_ref)), z_ref∈{0,L}, to keep every term ≤1 in
magnitude (naive exp(λL) overflows).

INTERNAL VALIDATION (test_zcb80_closed_form_self_validation), all in numpy:
  (i)   ODE residual of the analytic solution < 1e-5 (both legs);
  (ii)  independent 2nd-order FD BVP (2×2 block-tridiagonal Thomas, no scipy)
        converges to the analytic solution (< 3% at h≈1.9 mm, decreasing);
  (iii) undrained bulk pressure |p̂| ≈ q·Q̄/(D+α²Q̄) (92.3 Pa for Q̄=1e4 MPa,
        100 Pa for 1e9) — the physical anchor;
  (iv)  ω→0 gives the drained limit (|p̂|→0);
  (v)   dynamic-seepage on-vs-off changes the phasors < 1e-3 (k̄ρ_f ω²α ~1e-7
        vs iωα² ~3.4 — utterly negligible here; justifies -dynSeepage off).

===========================================================================
WHERE THE BENCHMARK LIVES ON THE ZCB80 CHART (why a graded mesh)
===========================================================================
With k'=1e-7 m/s the drainage is tiny: Π₁ = k'V_c²/(gβωL²) ≈ 3.6e-5 — DEEP in
ZCB80's undrained strip. The consolidation coefficient cv = k̄/(1/D+1/Q̄) ≈
8e-3 m²/s gives a periodic drainage boundary layer δ = √(cv/ω) ≈ 5 cm at the
drained top of a 30 m column. The interior pore pressure is therefore the
uniform UNDRAINED value; all the spatial action is a 5 cm layer under the top.

A uniform mesh cannot resolve this (a 3.75 m cell over a 5 cm layer produces the
classic Vermeer-Verruijt drainage-boundary overshoot: with 8 uniform cells the
node just below the top reads |p̂|≈207 vs the true 92 — a 45% profile L2 error).
The physically correct realisation, and "the cleanest triangulation" the pin
asks the agent to pick, is a strip mesh GEOMETRICALLY GRADED toward the drained
top (top cell ≈ few mm, growth ratio ≈1.5, ~18-20 cells spanning 30 m). Uniform
bisection of that graded mesh preserves the grading and serves as the
refinement ladder. Every cell is a straight-sided rectangle split into two BT6
triangles by its SW-NE diagonal (mid-edge nodes exactly at edge midpoints, so
the affine-map / straight-side precondition of the Bézier provider holds).

===========================================================================
MESH & LOADING details
===========================================================================
Strip mesh, 1 m wide × 30 m, one rectangle per depth level split into two BT6:
  T_lower vertices (SW,SE,NE), mid-edges (SW-SE bottom, SE-NE right, NE-SW diag)
  T_upper vertices (SW,NE,NW), mid-edges (SW-NE diag[shared], NE-NW top, NW-SW left)
Both windings are CCW (detJ>0, the BT6 element rejects detJ<0). Modeling dance:
vertices declared under model(ndf=3) (u_x,u_y,p carriers), mid-edge nodes under
model(ndf=2) (u_x,u_y only) — the heterogeneous-ndf pattern P3 introduces.
Element node order = 3 vertices then the 3 mid-edges on edges (v1,v2),(v2,v3),
(v3,v1); '-pOrder linear' selects Taylor-Hood; '-stab' is FATAL on TH.

BCs: ux=0 on ALL nodes (exact 1D confined column = ZCB80's assumption; the pin's
"side rollers ux fixed" taken to its 1D limit); uy=0 at the base; p=0 at the top
vertices (drained); sides and base carry NO p BC (impervious = natural).

LOAD — the quadratic-Bernstein edge subtlety: the top edge interpolates u with
the quadratic Bernstein basis {(1-t)², 2t(1-t), t²}. EACH Bernstein basis
function of degree n integrates to 1/(n+1) over the edge, so the consistent
nodal load for a uniform traction q on the 3-node quadratic edge is q·ℓ/3 at
EVERY node (the two vertices AND the mid-edge) — NOT the Lagrange 1/6-2/3-1/6
split. Here ℓ=1 m (width), thickness=1 m, q=100 Pa → q·ℓ/3 = 33.33 N per node,
scaled by a Trig time series sin(ωt). The 0.06% displacement agreement below
confirms the thirds weighting reproduces a uniform 100 Pa top traction.

===========================================================================
GATES & measured numbers (this build, UmfPack, Newmark γ=0.6/β=0.3025)
===========================================================================
Δt = T/128 ≈ 0.0145 s (T=2π/ω=1.859 s, so Δt≤T/100 ✓); ring up until per-cycle
amplitude drift < 0.5% (settles in 4 cycles, cap 15).

(a) Q̄=1e4 MPa leg (compressible), base graded mesh N=20 (top cell 4 mm, r=1.5):
      |p̂(z)| L2-rel = 1.47%  (gate ≤ 2%)
      |û(z)| L2-rel = 0.06%  (gate ≤ 2% — the pin's headline tolerance)
      phase p̂ L2 = 0.011 rad over the significant-|p| nodes (gate ≤ 0.1)
    The phase agreement (≈0, not π/2) confirms the Im-convention phasor sign.

(b) Q̄=1e9 MPa leg (undrained-limit checkerboard stressor), ladder N=16→32→64
    (bisected graded mesh):
      no-checkerboard: CB_fe/CB_cf = 0.98, 1.00, 1.00 — the FE vertex-p profile
        carries NO alternating (checkerboard) content beyond the true solution's
        own (both fall as the mesh refines). CB = |⟨p,χ⟩|/(‖p‖‖χ‖), χ_j=(-1)^j
        over interior vertices (the 1D analog of the B4 lattice metric).
      pressure L2 mesh-refinement convergence: 2.00% → 0.85% → 0.38% (monotone,
        rates 1.23, 1.18). The rate is ~1.2, NOT the O(h³) of quadratic-u,
        because the coupled 5 cm pressure boundary layer bottlenecks BOTH fields
        to ~1st order on the graded ladder (see the u-rate note).

(c) quadratic-u preserved — measured via a clean companion (test_..._quadratic_u):
    the pinned "u-convergence rate > 2" is UNACHIEVABLE in this near-undrained B1
    — global (and even bulk) u-convergence is bottlenecked to ~1.0-1.2 by the
    coupled pressure boundary layer (u = ∫ε dz carries the boundary-layer error
    through the momentum coupling; bulk-u rate measured 1.02/0.98). REPORTED to
    MAIN as a refutation. Quadratic-u preservation is instead demonstrated
    rigorously and mesh-independently: a static drained self-weight column has
    the EXACT parabolic displacement u(z)=(ρg/2D)(z²-2Lz); the quadratic-u BT6
    reproduces it to 1.5e-14 on any mesh (a linear element cannot) — the direct
    proof that the u space is genuinely quadratic.

(d) Δt adequacy — one Δt-halving (T/128 vs T/256, Q̄=1e4 base mesh): the phasor
    amplitudes move |p̂| 0.002%, |û| 0.000% (gate < 0.5%) — γ=0.6 numerical
    amplitude damping is already Δt-converged at T/128.

Runtime: ~15 s for the whole battery (well under the 120 s budget); ring-up
dominates but the graded base meshes are small (≤64 cells × 2 triangles).

Run: py -3.12 -m pytest tests/test_ladruno_up_element_th_b1.py -v
(the worktree dist/bin opensees.pyd is CPython 3.12; plain `python` is 3.11 and
will not import it. Solver = UmfPack — the honest-p tangent is unsymmetric.)
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
if not (os.path.isfile(os.path.join(_DIST, "opensees.pyd"))
        or os.path.isfile(os.path.join(_DIST, "opensees.so"))):
    pytest.skip(f"worktree engine not built: {_DIST}", allow_module_level=True)

os.environ["PATH"] = _DIST + os.pathsep + os.environ.get("PATH", "")
try:
    os.add_dll_directory(_DIST)
except (FileNotFoundError, OSError, AttributeError):
    pass
if _DIST not in sys.path:
    sys.path.insert(0, _DIST)
for _m in ("opensees", "openseespy", "openseespy.opensees"):
    sys.modules.pop(_m, None)
import opensees as ops  # noqa: E402

assert os.path.normcase(os.path.dirname(ops.__file__)) == os.path.normcase(_DIST), (
    f"wrong opensees.pyd imported: {ops.__file__} (want {_DIST}); run with "
    "py -3.12 so the boot .pth cannot preload another build")

pytestmark = [pytest.mark.zone_b]

_UMF = ("UmfPack",)   # honest-p tangent is unsymmetric (ADR §3.2)

# --------------------------------------------------------------------------
# B1 parameters (ADR-71 §7.1). ρ_s = 2000 is the SATURATED MIXTURE density
# (ZCB80 convention): the LadrunoUP mass block uses the material getRho()
# DIRECTLY as the mixture density (LadrunoUP.cpp:784,1002 — "SATURATED mixture,
# the pinned convention"), so the material rho passed below IS ρ used in ρü, and
# the closed form uses the SAME ρ=2000 → the two are self-consistent regardless
# of the grain-vs-mixture reading of "ρ_s".
_E, _NU = 7.492e8, 0.2
_D = _E * (1 - _NU) / ((1 + _NU) * (1 - 2 * _NU))      # constrained modulus
_RHO, _RHOF, _N = 2000.0, 1000.0, 0.333
_KBAR = 1.0e-7 / (_RHOF * 9.81)                        # k̄ = k'/(ρ_w g)
_OMEGA, _ALPHA, _L = 3.379, 1.0, 30.0
_T = 2.0 * np.pi / _OMEGA
_Q = 100.0                                             # top traction amplitude


# ==========================================================================
# INDEPENDENT numpy closed form (ZCB80 / u-p BVP)
# ==========================================================================
def zcb80_phasor(z, Qbar_Pa, omega=_OMEGA, dyn_seepage=False):
    """Complex phasors û(z), p̂(z), convention field(t)=Im{field_hat e^{iωt}}
    (matches a sin drive), forcing D û'(L)=-q. Numerically stable per-mode exp."""
    ooQ = 1.0 / Qbar_Pa
    ds = 1.0 if dyn_seepage else 0.0
    a = _KBAR * _D
    b = (_KBAR * _RHO * omega**2 - 1j * omega * _D * ooQ - 1j * omega * _ALPHA**2
         - ds * _KBAR * _RHOF * omega**2 * _ALPHA)
    c = -1j * omega * _RHO * omega**2 * ooQ
    disc = np.sqrt(b * b - 4 * a * c)
    mu = np.array([(-b + disc) / (2 * a), (-b - disc) / (2 * a)])
    lam = np.array([np.sqrt(mu[0]), -np.sqrt(mu[0]),
                    np.sqrt(mu[1]), -np.sqrt(mu[1])], dtype=complex)
    r = (_D * lam**2 + _RHO * omega**2) / (_ALPHA * lam)
    zr = np.where(lam.real > 0, _L, 0.0)               # keep exp arg ≤ 0

    def refexp(zbc):
        return np.exp(lam * (zbc - zr))
    M = np.array([refexp(0.0), r * lam * refexp(0.0),
                  r * refexp(_L), _D * lam * refexp(_L)], dtype=complex)
    amp = np.linalg.solve(M, np.array([0, 0, 0, -_Q], dtype=complex))
    z = np.atleast_1d(np.asarray(z, dtype=float))
    Ez = np.exp(np.outer(z, lam) - lam * zr)
    return Ez @ amp, Ez @ (amp * r), (lam, amp, r, zr)


def _fd_bvp(Qbar, N, omega=_OMEGA):
    """Independent 2nd-order FD of the same BVP via 2×2 block-tridiagonal
    Thomas (no scipy). Validates the analytic eigen-solution."""
    ooQ = 1.0 / Qbar
    h = _L / N
    nn = N + 1
    A = np.zeros((nn, 2, 2), dtype=complex)
    B = np.zeros((nn, 2, 2), dtype=complex)
    C = np.zeros((nn, 2, 2), dtype=complex)
    d = np.zeros((nn, 2), dtype=complex)
    for k in range(1, N):
        A[k, 0, 0] = _D / h**2
        B[k, 0, 0] = -2 * _D / h**2 + _RHO * omega**2
        C[k, 0, 0] = _D / h**2
        C[k, 0, 1] = -_ALPHA / (2 * h)
        A[k, 0, 1] = _ALPHA / (2 * h)
        A[k, 1, 1] = -_KBAR / h**2
        B[k, 1, 1] = 2 * _KBAR / h**2 + 1j * omega * ooQ
        C[k, 1, 1] = -_KBAR / h**2
        C[k, 1, 0] = 1j * omega * _ALPHA / (2 * h)
        A[k, 1, 0] = -1j * omega * _ALPHA / (2 * h)
    B[0, 0, 0] = 1.0                                   # u(0)=0
    B[0, 1, 1] = -1.0 / h; C[0, 1, 1] = 1.0 / h        # p'(0)=0
    B[N, 0, 0] = _D / h; A[N, 0, 0] = -_D / h; d[N, 0] = -_Q   # D u'(L)=-q
    B[N, 1, 1] = 1.0                                   # p(L)=0
    Cp = np.zeros((nn, 2, 2), dtype=complex)
    dp = np.zeros((nn, 2), dtype=complex)
    inv = np.linalg.inv(B[0])
    Cp[0] = inv @ C[0]; dp[0] = inv @ d[0]
    for k in range(1, nn):
        inv = np.linalg.inv(B[k] - A[k] @ Cp[k - 1])
        Cp[k] = inv @ C[k]
        dp[k] = inv @ (d[k] - A[k] @ dp[k - 1])
    x = np.zeros((nn, 2), dtype=complex)
    x[nn - 1] = dp[nn - 1]
    for k in range(nn - 2, -1, -1):
        x[k] = dp[k] - Cp[k] @ x[k + 1]
    return np.linspace(0, _L, nn), x[:, 0], x[:, 1]


# --------------------------------------------------------------------------
# small numpy helpers
# --------------------------------------------------------------------------
def _l2rel(a, b):
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-300))


def _cb(pv):
    """1D checkerboard index on the interior vertex-p profile: |<p,chi>|/
    (||p|| ||chi||), chi_j=(-1)^j (the 1D analog of the ADR B4 lattice mode)."""
    pi = pv[1:-1]
    chi = np.array([(-1) ** j for j in range(len(pi))], float)
    return abs(np.vdot(chi, pi)) / (np.linalg.norm(pi) * np.linalg.norm(chi) + 1e-30)


def _phase_l2(p_fe, p_cf, wthr=0.2):
    """RMS phase error (rad) over nodes where |p_cf| is significant."""
    m = np.abs(p_cf) > wthr * np.abs(p_cf).max()
    dphi = np.angle(p_fe[m]) - np.angle(p_cf[m])
    dphi = (dphi + np.pi) % (2 * np.pi) - np.pi
    return float(np.linalg.norm(dphi) / np.sqrt(max(m.sum(), 1)))


# --------------------------------------------------------------------------
# graded straight-sided BT6-TH strip mesh
# --------------------------------------------------------------------------
def _graded_levels(Ncell, top_cell, ratio):
    """y-levels 0..L, cells growing from `top_cell` at the drained top by
    `ratio` downward, normalised to span L (fine near the drainage layer)."""
    hs = np.array([top_cell * ratio ** k for k in range(Ncell)])
    hs = hs * _L / hs.sum()
    y = np.concatenate([[0.0], np.cumsum(hs[::-1])])
    y[-1] = _L
    return y


def _bisect(y):
    out = [y[0]]
    for j in range(len(y) - 1):
        out.append(0.5 * (y[j] + y[j + 1]))
        out.append(y[j + 1])
    return np.array(out)


def _build(Ncell, Qbar, ylevels):
    """1 m-wide column, Ncell rectangles (graded), each split into two straight-
    sided BT6 by the SW-NE diagonal. Vertices ndf=3, mid-edges ndf=2 (the TH
    modeling dance). Returns (vert{(i,j)->tag}, mid, Kf)."""
    ops.wipe()
    Kf = _N * Qbar                                     # 1/Q̄ = n/K_f (α=1, K_s→∞)
    ops.model("basic", "-ndm", 2, "-ndf", 3)           # vertex carriers
    tag = 1
    vert = {}
    for j in range(Ncell + 1):
        for i in range(2):
            ops.node(tag, float(i), float(ylevels[j])); vert[(i, j)] = tag; tag += 1
    ops.model("basic", "-ndm", 2, "-ndf", 2)           # mid-edge nodes (no p)
    mid = {}
    for j in range(Ncell + 1):                         # shared horizontal mids
        ops.node(tag, 0.5, float(ylevels[j])); mid[('h', j)] = tag; tag += 1
    for j in range(Ncell):                             # per-cell l/r/diag mids
        yc = 0.5 * (ylevels[j] + ylevels[j + 1])
        ops.node(tag, 0.0, yc); mid[('l', j)] = tag; tag += 1
        ops.node(tag, 1.0, yc); mid[('r', j)] = tag; tag += 1
        ops.node(tag, 0.5, yc); mid[('d', j)] = tag; tag += 1
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU, _RHO)
    e = 1
    for j in range(Ncell):
        SW, SE = vert[(0, j)], vert[(1, j)]
        NE, NW = vert[(1, j + 1)], vert[(0, j + 1)]
        mb, mt = mid[('h', j)], mid[('h', j + 1)]
        ml, mr, md = mid[('l', j)], mid[('r', j)], mid[('d', j)]
        # T_lower verts (SW,SE,NE): mids (SW,SE)=mb,(SE,NE)=mr,(NE,SW)=md ; CCW
        ops.element("LadrunoUP", e, SW, SE, NE, mb, mr, md, 1,
                    "-thick", 1.0, "-Kf", Kf, "-poro", _N, "-rhoF", _RHOF,
                    "-perm", _KBAR, _KBAR, "-pOrder", "linear",
                    "-dynSeepage", "off"); e += 1
        # T_upper verts (SW,NE,NW): mids (SW,NE)=md,(NE,NW)=mt,(NW,SW)=ml ; CCW
        ops.element("LadrunoUP", e, SW, NE, NW, md, mt, ml, 1,
                    "-thick", 1.0, "-Kf", Kf, "-poro", _N, "-rhoF", _RHOF,
                    "-perm", _KBAR, _KBAR, "-pOrder", "linear",
                    "-dynSeepage", "off"); e += 1
    # BCs: ux=0 everywhere (1D confined); uy=0 base; p=0 top vertices (drained)
    for (i, j), nd in vert.items():
        ops.fix(nd, 1, 1 if j == 0 else 0, 1 if j == Ncell else 0)
    for (typ, j), nd in mid.items():
        ops.fix(nd, 1, 1 if (typ == 'h' and j == 0) else 0)
    return vert, mid


def _run_phasor(Ncell, Qbar, ylevels, nppc=128, ncyc=15, drift_tol=5e-3):
    """Ring the sinusoidal drive to steady state, LS-fit A·sin+B·cos per vertex
    over the last cycle. Returns (z, p_hat, u_hat, cycles_used)."""
    vert, mid = _build(Ncell, Qbar, ylevels)
    dt = _T / nppc
    w = _Q * 1.0 / 3.0                                 # Bernstein thirds, ℓ=t=1
    ops.timeSeries("Trig", 1, 0.0, 1.0e9, _T, "-factor", 1.0)   # sin(ωt)
    ops.pattern("Plain", 1, 1)
    ops.load(vert[(0, Ncell)], 0.0, -w, 0.0)
    ops.load(vert[(1, Ncell)], 0.0, -w, 0.0)
    ops.load(mid[('h', Ncell)], 0.0, -w)               # mid-edge node (ndf=2)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system(*_UMF)
    ops.test("NormDispIncr", 1e-10, 20)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")

    probe = vert[(0, Ncell // 2)]
    prevA = None
    cyc = 0
    for cyc in range(ncyc):
        cp = []
        for _ in range(nppc):
            assert ops.analyze(1, dt) == 0, f"ring-up step failed (cyc {cyc})"
            cp.append(ops.nodeDisp(probe, 3))
        A = 0.5 * (max(cp) - min(cp))
        drift = abs(A - prevA) / max(abs(A), 1e-30) if prevA is not None else 1.0
        prevA = A
        if cyc >= 3 and drift < drift_tol:
            break

    vnodes = [vert[(0, j)] for j in range(Ncell + 1)]
    rec = {nd: [] for nd in vnodes}
    recu = {nd: [] for nd in vnodes}
    tt = []
    t = 0.0
    for _ in range(nppc):
        assert ops.analyze(1, dt) == 0
        t += dt
        tt.append(t)
        for nd in vnodes:
            rec[nd].append(ops.nodeDisp(nd, 3))
            recu[nd].append(ops.nodeDisp(nd, 2))
    tt = np.array(tt)
    G = np.column_stack([np.sin(_OMEGA * tt), np.cos(_OMEGA * tt)])

    def fit(series):                                   # field_hat = A + iB
        coef, *_ = np.linalg.lstsq(G, np.array(series), rcond=None)
        return coef[0] + 1j * coef[1]
    z = np.array(ylevels)
    p = np.array([fit(rec[nd]) for nd in vnodes])
    u = np.array([fit(recu[nd]) for nd in vnodes])
    return z, p, u, cyc + 1


# base meshes (graded toward the drained top)
_BASE_ACC = (20, 0.004, 1.5)      # Q̄=1e4 accuracy leg — top cell 4 mm
_BASE_LADDER = (16, 0.006, 1.6)   # Q̄=1e9 convergence ladder base — top cell 6 mm


# ==========================================================================
# 0. closed-form self-validation (numpy only — the reference must be right)
# ==========================================================================
@pytest.mark.t0
def test_zcb80_closed_form_self_validation():
    """Independent checks that the numpy ZCB80 closed form is correct BEFORE it
    gates the FE: ODE residual, an independent FD-BVP solve, the undrained bulk
    anchor, the drained ω→0 limit, dynamic-seepage negligibility, and the
    Bernstein-thirds edge-load weight."""
    # (i) ODE residual of the analytic solution (both legs)
    for Qbar in (1e10, 1e15):
        _, _, (lam, amp, r, zr) = zcb80_phasor(np.array([0.0]), Qbar)
        zt = np.array([5.0, 15.0, 25.0])
        Ez = np.exp(np.outer(zt, lam) - lam * zr)
        u = Ez @ amp; up = Ez @ (amp * lam); upp = Ez @ (amp * lam**2)
        p = Ez @ (amp * r); pp = Ez @ (amp * r * lam); ppp = Ez @ (amp * r * lam**2)
        res1 = _D * upp + _RHO * _OMEGA**2 * u - _ALPHA * pp
        res2 = -_KBAR * ppp + (1j * _OMEGA / Qbar) * p + 1j * _OMEGA * _ALPHA * up
        s1 = abs(_D * upp).max()
        s2 = max(abs(_KBAR * ppp).max(), abs(1j * _OMEGA * _ALPHA * up).max())
        assert np.abs(res1).max() / s1 < 1e-5, f"ODE1 residual {Qbar:.0e}"
        assert np.abs(res2).max() / s2 < 1e-5, f"ODE2 residual {Qbar:.0e}"

    # (ii) independent FD-BVP converges to the analytic solution
    for Qbar in (1e10, 1e15):
        prev = None
        for Nf in (4000, 16000):
            zf, uf, pf = _fd_bvp(Qbar, Nf)
            ua, pa, _ = zcb80_phasor(zf, Qbar)
            ep = _l2rel(pf, pa)
            if prev is not None:
                assert ep < prev, f"FD not converging ({Qbar:.0e})"
            prev = ep
        assert ep < 0.03, f"FD-vs-analytic {Qbar:.0e}: {ep:.2%} (want < 3%)"

    # (iii) undrained bulk anchor |p̂| ≈ q·Q̄/(D+α²Q̄)
    for Qbar in (1e10, 1e15):
        _, pa, _ = zcb80_phasor(np.array([_L / 2]), Qbar)
        p_und = _Q * Qbar / (_D + _ALPHA**2 * Qbar)
        assert abs(abs(pa[0]) - p_und) / p_und < 0.01, (
            f"bulk |p̂| {abs(pa[0]):.2f} vs undrained {p_und:.2f} ({Qbar:.0e})")

    # (iv) ω→0 drained limit: excess pressure vanishes. The drainage length
    # √(cv/ω) must swamp L (30 m); at ω=1e-8 it is ~890 m so |p̂| → 0 (a monotone
    # sweep 1e-6→1e-9 gives 5.3, 0.53, 0.053, 0.0053 — clean 1st-order-in-ω decay)
    _, p0, _ = zcb80_phasor(np.linspace(0, _L, 7), 1e10, omega=1e-8)
    assert np.abs(p0).max() < 0.005 * _Q, "ω→0 should be drained (p→0)"

    # (v) dynamic-seepage on-vs-off is negligible (justifies -dynSeepage off)
    zc = np.linspace(0, _L, 9)
    for Qbar in (1e10, 1e15):
        u0, p0, _ = zcb80_phasor(zc, Qbar, dyn_seepage=False)
        u1, p1, _ = zcb80_phasor(zc, Qbar, dyn_seepage=True)
        assert _l2rel(p1, p0) < 1e-3 and _l2rel(u1, u0) < 1e-3, (
            f"dynamic seepage not negligible ({Qbar:.0e})")

    # (vi) Bernstein quadratic edge basis integrates to 1/3 each (the load rule
    # — the consistent nodal load is q·ℓ/3 per node, NOT Lagrange 1/6-2/3-1/6)
    t = np.linspace(0, 1, 20001)
    dt = t[1] - t[0]
    for Bi in ((1 - t)**2, 2 * t * (1 - t), t**2):
        integ = (0.5 * Bi[0] + Bi[1:-1].sum() + 0.5 * Bi[-1]) * dt   # trapezoid
        assert abs(integ - 1.0 / 3.0) < 1e-4


# ==========================================================================
# 1. B1 Q̄=1e4 leg — phasor accuracy gate (the pin's ≤2% headline)
# ==========================================================================
@pytest.mark.t1
def test_b1_q1e4_phasor_accuracy():
    """Q̄=1e4 MPa (compressible) leg: reconstruct |p̂|, arg p̂, |û| on the base
    graded mesh and compare to the closed form.

    Gates (measured, this build): |p̂| L2-rel 1.47% (≤ 2%), |û| L2-rel 0.06%
    (≤ 2%), phase-p L2 0.011 rad (≤ 0.1 — the ≈0 value, not π/2, confirms the
    Im-convention phasor sign is correct). Ring-up settles in 4 cycles."""
    Nc, tc, rr = _BASE_ACC
    yl = _graded_levels(Nc, tc, rr)
    z, p_fe, u_fe, cyc = _run_phasor(Nc, 1e10, yl)
    u_cf, p_cf, _ = zcb80_phasor(z, 1e10)

    assert cyc <= 15, f"ring-up did not settle within 15 cycles ({cyc})"
    ep = _l2rel(np.abs(p_fe), np.abs(p_cf))
    eu = _l2rel(np.abs(u_fe[1:]), np.abs(u_cf[1:]))   # z=0 is fixed (u=0)
    eph = _phase_l2(p_fe, p_cf)
    assert ep <= 0.02, f"|p̂| L2-rel {ep:.3%} > 2% (Q̄=1e4)"
    assert eu <= 0.02, f"|û| L2-rel {eu:.3%} > 2% (Q̄=1e4)"
    assert eph <= 0.1, f"phase p̂ L2 {eph:.4f} rad > 0.1 (sign/convention?)"
    # sanity: the interior really is at the undrained level (a dead field would
    # also pass a small L2), and the top vertex is the drained p=0 node
    p_und = _Q * 1e10 / (_D + 1e10)
    assert abs(abs(p_fe[len(z) // 2]) - p_und) / p_und < 0.02
    assert abs(p_fe[-1]) < 1e-6 * p_und


# ==========================================================================
# 2. B1 Q̄=1e9 leg — no checkerboard + pressure mesh-convergence
# ==========================================================================
@pytest.mark.t1
def test_b1_q1e9_no_checkerboard_and_convergence():
    """Q̄=1e9 MPa (undrained-limit checkerboard stressor) leg, ladder N=16→32→64
    (bisected graded mesh):

      * no-checkerboard: the FE vertex-p checkerboard index CB stays COMPARABLE
        to the closed form's own on the same nodes (CB_fe/CB_cf = 0.98, 1.00,
        1.00) — the Taylor-Hood pair introduces NO spurious alternating mode at
        the near-incompressible undrained limit (the whole point of P3 TH);
      * pressure L2 mesh-refinement convergence: 2.00% → 0.85% → 0.38%,
        monotone, positive rate (~1.2 — bottlenecked by the 5 cm coupled
        pressure boundary layer, not the quadratic-u order; see the module
        docstring gate (b)/(c))."""
    Nc, tc, rr = _BASE_LADDER
    yl = _graded_levels(Nc, tc, rr)
    eps = []
    for _lvl in range(3):
        N = len(yl) - 1
        z, p_fe, u_fe, cyc = _run_phasor(N, 1e15, yl)
        u_cf, p_cf, _ = zcb80_phasor(z, 1e15)
        cb_fe, cb_cf = _cb(np.abs(p_fe)), _cb(np.abs(p_cf))
        # FE carries no MORE alternating content than the true solution sampled
        # on the same (graded) nodes — no spurious checkerboard mode
        assert cb_fe <= 1.3 * cb_cf + 1e-3, (
            f"N={N}: CB_fe {cb_fe:.4f} >> CB_cf {cb_cf:.4f} — checkerboard?")
        eps.append(_l2rel(np.abs(p_fe), np.abs(p_cf)))
        yl = _bisect(yl)
    # monotone mesh-refinement convergence of the pressure profile
    assert eps[1] < eps[0] and eps[2] < eps[1], (
        f"pressure L2 not monotonically converging: {eps}")
    r1 = np.log(eps[0] / eps[1]) / np.log(2.0)
    r2 = np.log(eps[1] / eps[2]) / np.log(2.0)
    assert r1 > 0.8 and r2 > 0.8, f"pressure convergence rates too low: {r1:.2f},{r2:.2f}"
    assert eps[2] < 0.006, f"finest-mesh pressure L2 {eps[2]:.3%} (want < 0.6%)"


# ==========================================================================
# 3. quadratic-u preserved (companion exactness) + B1 u-rate refutation
# ==========================================================================
def _selfweight_parabola(Ncell, ylevels):
    """Static DRAINED (p=0 everywhere) confined column under a self-weight body
    force (ctor -body): pure quadratic-u elastostatics with the EXACT solution
    u(z)=(ρg/2D)(z²-2Lz). Returns (rc, z, u_fe, u_exact)."""
    g = 9.81
    ops.wipe()
    Kf = _N * 1e10
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    tag = 1
    vert = {}
    for j in range(Ncell + 1):
        for i in range(2):
            ops.node(tag, float(i), float(ylevels[j])); vert[(i, j)] = tag; tag += 1
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    mid = {}
    for j in range(Ncell + 1):
        ops.node(tag, 0.5, float(ylevels[j])); mid[('h', j)] = tag; tag += 1
    for j in range(Ncell):
        yc = 0.5 * (ylevels[j] + ylevels[j + 1])
        ops.node(tag, 0.0, yc); mid[('l', j)] = tag; tag += 1
        ops.node(tag, 1.0, yc); mid[('r', j)] = tag; tag += 1
        ops.node(tag, 0.5, yc); mid[('d', j)] = tag; tag += 1
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU, _RHO)
    e = 1
    for j in range(Ncell):
        SW, SE = vert[(0, j)], vert[(1, j)]
        NE, NW = vert[(1, j + 1)], vert[(0, j + 1)]
        mb, mt = mid[('h', j)], mid[('h', j + 1)]
        ml, mr, md = mid[('l', j)], mid[('r', j)], mid[('d', j)]
        for tri in ((SW, SE, NE, mb, mr, md), (SW, NE, NW, md, mt, ml)):
            ops.element("LadrunoUP", e, *tri, 1, "-thick", 1.0, "-Kf", Kf,
                        "-poro", _N, "-rhoF", _RHOF, "-perm", _KBAR, _KBAR,
                        "-pOrder", "linear", "-dynSeepage", "off",
                        "-body", 0.0, -g); e += 1
    for (i, j), nd in vert.items():
        ops.fix(nd, 1, 1 if j == 0 else 0, 1)          # drained (p=0), base uy=0
    for (typ, j), nd in mid.items():
        ops.fix(nd, 1, 1 if (typ == 'h' and j == 0) else 0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.system(*_UMF)
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    rc = ops.analyze(1)
    z = np.array(ylevels)
    u = np.array([ops.nodeDisp(vert[(0, j)], 2) for j in range(Ncell + 1)])
    uex = (_RHO * g / (2 * _D)) * (z**2 - 2 * _L * z)
    return rc, z, u, uex


@pytest.mark.t1
def test_b1_quadratic_u_preserved():
    """Quadratic-u preservation.

    The pinned 'u-convergence rate > 2' is UNACHIEVABLE in the near-undrained B1
    (REPORTED to MAIN): the coupled 5 cm pressure boundary layer bottlenecks
    global — and even bulk — u-convergence to ~1.0-1.2 (u=∫ε dz inherits the
    pressure rate through the momentum coupling). Instead, quadratic-u is proven
    directly and mesh-independently: a static drained self-weight column has the
    EXACT parabola u(z)=(ρg/2D)(z²-2Lz), which the quadratic-u BT6 space
    reproduces to machine precision on ANY mesh (a linear u space cannot).
    Measured: 1.5e-14 (uniform 8) / 4.8e-13 (graded 16)."""
    for Ncell, levels in ((8, _graded_levels(8, 0.5, 1.3)),
                          (16, _graded_levels(16, 0.006, 1.6))):
        rc, z, u, uex = _selfweight_parabola(Ncell, levels)
        assert rc == 0, f"static drained self-weight solve failed (N={Ncell})"
        err = np.max(np.abs(u - uex)) / max(np.max(np.abs(uex)), 1e-30)
        assert err < 1e-8, (
            f"quadratic-u parabola not exact (N={Ncell}): {err:.2e} — the BT6 "
            "u space should represent u(z)=(ρg/2D)(z²-2Lz) exactly")


# ==========================================================================
# 4. Δt adequacy — one Δt-halving check (pins Δt ≤ T/100)
# ==========================================================================
@pytest.mark.t1
def test_b1_dt_halving_adequacy():
    """One Δt-halving on the Q̄=1e4 base mesh (T/128 vs T/256): the reconstructed
    phasor amplitudes move < 0.5% (measured |p̂| 0.002%, |û| 0.000%). γ=0.6
    numerical amplitude damping is already Δt-converged at Δt=T/128 ≈ 0.0145 s
    (≤ T/100), so the steady-state amplitudes are not Δt-biased."""
    Nc, tc, rr = _BASE_ACC
    yl = _graded_levels(Nc, tc, rr)
    _, p1, u1, _ = _run_phasor(Nc, 1e10, yl, nppc=128)
    _, p2, u2, _ = _run_phasor(Nc, 1e10, yl, nppc=256)
    dp = _l2rel(np.abs(p1), np.abs(p2))
    du = _l2rel(np.abs(u1), np.abs(u2))
    assert dp < 0.005, f"|p̂| moves {dp:.3%} on Δt halving (want < 0.5%)"
    assert du < 0.005, f"|û| moves {du:.3%} on Δt halving (want < 0.5%)"
