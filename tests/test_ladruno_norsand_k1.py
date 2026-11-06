"""WP-144 gate G2 (Zone A): the K1 CLOSED-FORM gates of the equation sheet THROUGH OpenSees (stdBrick).

Scope.  Ladruno_files/testbed/norsand_oracle/tests/test_g1_k1_closed_forms.py proves K1.1 - K1.10 on the two ORACLES.
This file proves the ones that matter for the shipped object on the real `nDMaterial LadrunoNorSand` inside a real
`stdBrick`, driven by the real analysis stack, with the expected values taken from the SHEET only (section 13), never
from the oracles and never from a run.  Zone A: openseespy + numpy, no scipy / sympy / oracle import.

  K1.2   closed non-coaxial elastic strain loop through a stdBrick: zero net work <= 1e-12, the state returns   (13.2)
  K1.7   undrained (isochoric) TXC critical-state endpoint p_cs in BOTH CSL modes, closed form                  (13.7)
  K1.8   drained critical-state asymptote: eta -> M, psi_i -> 0, pi_i -> p, D -> 0                              (13.8)
  K1.9   dissipation census D >= 0 on every step of the K1.7 / K1.8 paths (the `D` response)                    (13.9)
  K1.5   flow-rule ratio d eps^p_v / d eps^p_s = sqrt(3/2) beta F_p / Omega at every backward-Euler step        (13.5)

WHERE EVERY EXPECTED VALUE COMES FROM (none is harvested from the shell's own output).
  * sheet 13.2   W = oint sigma:d eps = 0 (<= 1e-12 of oint |sigma:d eps|) and the state returns to round-off for any
                 closed strain loop (sigma = dPsi/d eps^e).  The work integral along a straight strain segment is taken
                 by 12-point Gauss-Legendre in the segment parameter t: sigma(E_j + t D_j) is exp(linear) over a range of
                 ~0.4 in the exponent, for which GL-12 is exact to ~1e-15 (the G1 argument).  The strain is prescribed
                 DOF by DOF through a Path time series, so every segment is a straight line in strain space; the elastic
                 response is path independent, so the Gauss nodes are visited by variable load steps.  A non-vacuity
                 control: the work of the FIRST (open) segment equals the stored-energy difference
                 Psi(E1) - Psi(0) = -kappa_hat p0 (exp(-eps_v/kappa_hat) - 1) + mu0 e:e (sheet S.5, alpha0 = 0:
                 d Psi = p d eps_v + 2 mu0 e:de), 1e-10 relative; the vertex stresses equal (S.5).
  * sheet 13.7   isochoric => v = v0 exp(tr eps) = v0 (tr eps = 0); critical state <=> psi_i = 0, pi_i = p, D = 0,
                 eta = M (theta = pi/3): p_cs = -exp((v_c0 - v)/lambda_tilde) (paper) or
                 -p_a ((e0 - e)/lambda_c)^(1/xi), e = v - 1 (fork); q_cs = M |p_cs| / zeta(pi/3) = M |p_cs| in TXC.  v0 is
                 built by inverting that closed form for a target p_cs = -250 kPa (sheet 6, S.22).
                 Tolerances are the G1 ones, derived there from the sheet: the CS is a fixed point of the backward-Euler
                 map (D = 0 => eps^p_v increment 0 => p stationary) and an attractor, psi_i relaxing like exp(-|chi| eps^p_s)
                 (|chi| = 3.5): at |eps_ax| = 5 the residual is ~ e^-17..e^-8 ~ 1e-8, so p within 1e-6 relative and the
                 derived quantities within 1e-5, and the error at |eps_ax| = 2, 3.5, 5 non-increasing.
  * sheet 13.8   at large eps_s: psi_i -> 0, pi_i -> p, eta -> M (q/|p| -> M / zeta(theta)), D -> 0.  Tolerances (G1,
                 relaxation argument e^-10 at |eps_ax| = 3): eta/M - 1 within 2e-3, |psi_i| 1e-3, |pi_i/p - 1| and the
                 closed-form dilatancy D(eta) 5e-3, and the eta error at 3.0 at most half of that at 1.5 (not stalled).
  * sheet 13.9   D^p = Delta lambda sum_a sigma_a q_a >= 0 every plastic step, = 0 on elastic ones (round-off floor
                 1e-9 |p0| |step| as in G1).
  * sheet 13.5   for backward Euler Delta eps^p = Delta lambda q(sigma_{n+1}, pi_{i,n+1}) EXACTLY, so over ONE
                 unsubstepped step Delta eps^p_v / Delta eps^p_s = sqrt(3/2) beta F_p / Omega at the converged state
                 (sheet: 2.5e-15 over a full TXC path; G1 gate 1e-8 relative).  F_p from (S.12), Omega from (S.20) with
                 y_a of (S.6) and the corner value zeta_bar_y(pi/3) = +sqrt6 zeta_bar''(pi/3)/9 (S.9; 0.8164965809 for
                 rho_bar = 0.8, the number printed in the sheet 4.2 table); 2-invariant (rho = rho_bar = 1): Omega =
                 sqrt(3/2) exactly and D = (eta - M)/(1 - N_bar).  Identity per UNSUBSTEPPED step only: over a
                 substepped step the ratio is a weighted mean (sheet 13.5), so the check is gated on stepInfo
                 substeps == 1 (a premise count guards it).

DRIVERS.  Undrained: a unit stdBrick with every nodal DOF prescribed (the isochoric strain diag(1/2, 1/2, -1) |eps_ax|,
no Newton tolerance can enter).  Drained: the symmetric octant of a triaxial sample, one stdBrick; the three symmetry
planes are rollers, the top face z = 1 is displaced (axial strain, displacement control: a drained test must pass the
peak and reach the critical state, which load control cannot) and the x = 1 / y = 1 faces carry the constant -100 kPa
confining tractions (lumped, -25 per node); every Gauss point sees the uniform solution.  `system FullGeneral`
(NON-symmetric consistent tangent), Newton, NormDispIncr 1e-12.

WALL TIME (measured on the dev box, Windows, CPython 3.12, 2026-10-01, whole file): 89 s for 17 tests, of which the
five path runs (two undrained 148-step prescribed runs ~2 s each, three drained 191-step Newton runs ~25 s each) are
shared (cached) by the K1.7 / K1.8 / K1.9 / K1.5 tests; K1.2 takes 0.2 s.  Nothing is marked slow.  Measured at write time on
the pre-decision (linear v-update) binary: every closed form here holds on either v law (the K1 identities involve v only
through psi_i at the shell's own v, and the isochoric paths have tr eps = 0), K1.5 to 2.6e-15.

Every test names, in a `KILLS:` line, the mutant (a silent regression of the shell or the kernel) it fails on.
"""
import math
import os

import numpy as np
import pytest

# ---- engine binding (same idiom as tests/test_ladruno_norsand.py) -----------
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
_DIST_BIN = os.path.join(_ROOT, "dist", "bin")
if os.path.isfile(os.path.join(_DIST_BIN, "opensees.pyd")):
    from _engine import bind_worktree_engine
    ops = bind_worktree_engine(_DIST_BIN)
else:
    from _testbed import ops

pytestmark = [pytest.mark.zone_a]

SQ32 = math.sqrt(1.5)
SQ6 = math.sqrt(6.0)
I3 = np.eye(3)

# ===========================================================================
#  Parameter sets (sheet 14 / 15, the G1 sets): K2 energy + surface, paper or fork CSL
# ===========================================================================
_P0, _KH, _MU = -100.0, 0.01, 5400.0
_M, _N, _NB, _CHI, _H = 1.2, 0.4, 0.2, -3.5, 280.0
_PAPER = dict(rho=0.7, rho_bar=0.8, csl="paper", lambda_tilde=0.0135, v_c0=1.81)
_FORK = dict(rho=0.7, rho_bar=0.8, csl="fork", lambda_tilde=None, v_c0=None, e0=0.83, lambda_c=0.027, xi=0.45,
             p_a=101.325)
_PAPER_2INV = dict(_PAPER, rho=1.0, rho_bar=1.0)
_P_CS_TARGET = -250.0                      # the undrained CS pressure the initial v0 is built to give


def _mat_args(kw, v0, pi0, **over):
    d = {"p0": _P0, "kappa_hat": _KH, "mu0": _MU, "M": _M, "N": _N, "N_bar": _NB, "chi": _CHI, "h": _H}
    d.update(kw)
    d.update(over)
    d["v0"], d["pi0"] = v0, pi0
    out = []
    for k, v in d.items():
        if v is None:
            continue
        out += ["-" + k, v]
    return out


# ---- sheet closed forms, written here from the sheet (not from the oracles or the shell) ----------------
def _v0_for_psi(kw, pi0, psi0):
    """Specific volume giving image state parameter psi0 at pi_i0 (sheet 6, S.22)."""
    if kw["csl"] == "paper":
        return psi0 + kw["v_c0"] - kw["lambda_tilde"] * math.log(-pi0)
    return 1.0 + kw["e0"] + psi0 - kw["lambda_c"] * (-pi0 / kw["p_a"]) ** kw["xi"]


def _v0_for_pcs(kw, pcs):
    """Invert the undrained CS closed form (sheet 13.7) for v0."""
    if kw["csl"] == "paper":
        return kw["v_c0"] - kw["lambda_tilde"] * math.log(-pcs)
    return 1.0 + kw["e0"] - kw["lambda_c"] * (-pcs / kw["p_a"]) ** kw["xi"]


def _pcs_closed(kw, v):
    """Sheet 13.7: isochoric => v constant, CS <=> psi_i = 0 with pi_i = p."""
    if kw["csl"] == "paper":
        return -math.exp((kw["v_c0"] - v) / kw["lambda_tilde"])
    e = v - 1.0
    return -kw["p_a"] * ((kw["e0"] - e) / kw["lambda_c"]) ** (1.0 / kw["xi"])


def _psi(kw, v, pi):
    """psi_i = v - v_c(pi_i) (S.22)."""
    if kw["csl"] == "paper":
        return v - kw["v_c0"] + kw["lambda_tilde"] * math.log(-pi)
    return (v - 1.0) - kw["e0"] + kw["lambda_c"] * (-pi / kw["p_a"]) ** kw["xi"]


def _eta_S12(p, pi, M=_M, N=_N):
    return (M / N) * (1.0 - (1.0 - N) * (p / pi) ** (N / (1.0 - N)))


def _pq_theta(sig6):
    """p, q, principal values of a Voigt stress (S.1); the brick paths here are axisymmetric, so theta is a corner."""
    S = np.array([[sig6[0], sig6[3], sig6[5]], [sig6[3], sig6[1], sig6[4]], [sig6[5], sig6[4], sig6[2]]])
    w = np.linalg.eigvalsh(0.5 * (S + S.T))
    p = float(w.sum()) / 3.0
    xi = w - p
    return p, SQ32 * float(np.linalg.norm(xi)), w


def _zeta_bar_corner(rho_bar, corner):
    """(zeta_bar, zeta_bar_y) of S.9 at an axisymmetric corner, Willam-Warnke: 'C' = theta pi/3 (TXC), 'E' = theta 0
    (TXE).  zeta_y(pi/3) = +sqrt6 zeta''(pi/3)/9, zeta_y(0) = -sqrt6 zeta''(0)/9 (S.9), with the zeta'' of the sheet 4.2
    table at rho_bar = 0.8 (zeta''(pi/3) = 3.0, zeta''(0) = -0.75: zeta_y = 0.8164965809 / 0.2041241452, the numbers
    printed there), and everything zero at rho_bar = 1 (zeta = 1)."""
    if rho_bar == 1.0:
        return 1.0, 0.0
    if rho_bar == 0.8:
        return (1.0, +0.8164965809) if corner == "C" else (1.0 / 0.8, +0.2041241452)
    raise ValueError("only the sheet-4.2 table rows rho_bar = 0.8 and 1.0 are used here")


def _omega(rho_bar, sig6, corner):
    """Omega of (S.20) at an axisymmetric corner stress: Omega^2 = (3/2) zeta_b^2 + (zeta_b_y q)^2 sum_a y_a^2, with y_a
    of (S.6) (delta_a = 1)."""
    p, q, w = _pq_theta(sig6)
    zb, zby = _zeta_bar_corner(rho_bar, corner)
    xi = w - p
    R = float(np.linalg.norm(xi))
    ya = 3.0 * xi ** 2 / R ** 3 - 3.0 * float((xi ** 3).sum()) * xi / R ** 5 - 1.0 / R
    return math.sqrt(1.5 * zb ** 2 + (zby * q) ** 2 * float((ya ** 2).sum()))


def _dilatancy_closed(kw, sig6, pi, corner):
    """D = sqrt(3/2) beta F_p / Omega, beta F_p = (eta - M)/(1 - N_bar), eta from (S.12) at the state (sheet 5.3)."""
    p, _, _ = _pq_theta(sig6)
    return SQ32 * (_eta_S12(p, pi) - _M) / ((1.0 - _NB) * _omega(kw["rho_bar"], sig6, corner))


# ===========================================================================
#  Decks
# ===========================================================================
_CUBE = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]


def _unit_cube(mat_tag):
    for i, (x, y, z) in enumerate(_CUBE, 1):
        ops.node(i, float(x), float(y), float(z))
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, mat_tag)


def _stack(test_tol=1e-12, test_iter=50):
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")                  # the consistent tangent is NON-symmetric
    ops.test("NormDispIncr", test_tol, test_iter, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


def _deck_prescribed(kw, v0, pi0, times, vertices, **over):
    """Unit stdBrick, EVERY nodal DOF prescribed: u(t) = E(t) x with E(t) piecewise linear through `vertices` (3x3 strain
    tensors) at pseudo-times `times`.  One Path series per DOF (24), so each segment is a straight line in strain space.
    The series is held at its last value for one more unit of pseudo-time: a Path series is ZERO beyond its last time, and
    the accumulated load-control time can land 1e-16 past the last vertex (measured: it silently unloaded the last
    steps of the first version of the undrained run)."""
    times = [float(t) for t in times] + [float(times[-1]) + 1.0]
    vertices = list(vertices) + [vertices[-1]]
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("LadrunoNorSand", 1, *_mat_args(kw, v0, pi0, **over))
    _unit_cube(1)
    tag = 0
    for nid, xyz in enumerate(_CUBE, 1):
        for d in range(3):
            vals = [float(sum(E[d][k] * xyz[k] for k in range(3))) for E in vertices]
            tag += 1
            ops.timeSeries("Path", tag, "-time", *times, "-values", *vals)
            ops.pattern("Plain", tag, tag)
            ops.sp(nid, d + 1, 1.0)
    _stack()


def _deck_drained(kw, v0, pi0, ax_total):
    """Symmetric octant of a drained triaxial sample.  Rollers on x = 0, y = 0, z = 0; the top face (z = 1) is displaced by
    ax_total * lambda (compression < 0); constant -25 per node on the x = 1 and y = 1 faces (sigma = -100 kPa over a unit
    face).  Free DOFs: x of nodes 2, 3, 6, 7 and y of nodes 3, 4, 7, 8."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("LadrunoNorSand", 1, *_mat_args(kw, v0, pi0))
    _unit_cube(1)
    for nid, (x, y, z) in enumerate(_CUBE, 1):
        fx, fy, fz = int(x == 0), int(y == 0), int(z == 0)
        if fx or fy or fz:
            ops.fix(nid, fx, fy, fz)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for nid, (x, y, z) in enumerate(_CUBE, 1):
        fx = -25.0 if x == 1 else 0.0
        fy = -25.0 if y == 1 else 0.0
        if fx or fy:
            ops.load(nid, fx, fy, 0.0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for nid, (x, y, z) in enumerate(_CUBE, 1):
        if z == 1:
            ops.sp(nid, 3, float(ax_total))
    _stack()


def _gp(name):
    return [float(x) for x in ops.eleResponse(1, "material", 1, name)]


def _read():
    return dict(stress=_gp("stress"), state=_gp("state"), info=_gp("stepInfo"), D=_gp("D")[0])


# stage schedules: (|axial strain of the stage|, number of steps), the G1 ones (steps 2e-3, 1e-2 / 2e-2, 5e-2 / 2.5e-2)
U_SCHED = [(0.04, 20), (0.96, 48), (4.0, 80)]               # undrained: |eps_ax| = 5.0
TX_SCHED = [(0.04, 20), (0.46, 46), (2.5, 125)]             # drained:   |eps_ax| = 3.0 (finer last stage: Newton)
TX_SCHED_SHORT = [(0.04, 20), (0.46, 46)]                   # drained 2-invariant flow-rule case: |eps_ax| = 0.5


def _steps_of(sched):
    """[(axial strain step size, cumulative |axial strain| after the step)]."""
    out, c = [], 0.0
    for tot, n in sched:
        for _ in range(n):
            c += tot / n
            out.append((tot / n, c))
    return out


# ===========================================================================
#  Path runs (cached): each returns the per-step records read after every analyze()
# ===========================================================================
_UND_CASES = {
    "UND_paper": dict(kw=_PAPER, pi0=-60.0),
    "UND_fork": dict(kw=_FORK, pi0=-60.0),
}
_TX_CASES = {   # kw, sign (-1 compression), pi0, psi0, schedule, corner
    "TXC_paper": dict(kw=_PAPER, sign=-1, pi0=-80.0, psi0=-0.05, sched=TX_SCHED, corner="C"),
    "TXC_fork": dict(kw=_FORK, sign=-1, pi0=-80.0, psi0=-0.05, sched=TX_SCHED, corner="C"),
    "TXE_paper": dict(kw=_PAPER, sign=+1, pi0=-80.0, psi0=-0.05, sched=TX_SCHED, corner="E"),
    "TXC_2inv": dict(kw=_PAPER_2INV, sign=-1, pi0=-80.0, psi0=-0.05, sched=TX_SCHED_SHORT, corner="C"),
}
_RUNS = {}


def _run_undrained(case):
    key = ("U", case)
    if key in _RUNS:
        return _RUNS[key]
    spec = _UND_CASES[case]
    kw = spec["kw"]
    v0 = _v0_for_pcs(kw, _P_CS_TARGET)
    ax_total = 5.0
    E = np.diag([0.5 * ax_total, 0.5 * ax_total, -ax_total])           # isochoric TXC, |eps_ax| = 5
    _deck_prescribed(kw, v0, spec["pi0"], [0.0, 1.0], [np.zeros((3, 3)), E])
    recs = []
    for step, cum in _steps_of(U_SCHED):
        ops.integrator("LoadControl", step / ax_total)
        assert ops.analyze(1) == 0, f"{case}: analyze failed at |eps_ax| = {cum:.4f}"
        r = _read()
        r.update(cum=cum, step=step)
        recs.append(r)
    out = dict(kw=kw, v0=v0, pi0=spec["pi0"], recs=recs, corner="C")
    _RUNS[key] = out
    return out


def _run_drained(case):
    key = ("D", case)
    if key in _RUNS:
        return _RUNS[key]
    spec = _TX_CASES[case]
    kw = spec["kw"]
    v0 = _v0_for_psi(kw, spec["pi0"], spec["psi0"])
    ax_total = sum(t for t, _ in spec["sched"])
    _deck_drained(kw, v0, spec["pi0"], spec["sign"] * ax_total)
    recs = []
    for step, cum in _steps_of(spec["sched"]):
        ops.integrator("LoadControl", step / ax_total)
        assert ops.analyze(1) == 0, f"{case}: Newton did not converge at |eps_ax| = {cum:.4f}"
        r = _read()
        r.update(cum=cum, step=step)
        recs.append(r)
    out = dict(kw=kw, v0=v0, pi0=spec["pi0"], recs=recs, corner=spec["corner"])
    _RUNS[key] = out
    return out


def _run(case):
    return _run_undrained(case) if case in _UND_CASES else _run_drained(case)


# ===========================================================================
#  K1.2  closed non-coaxial elastic strain loop through a stdBrick (sheet 13.2)
# ===========================================================================
# vertices of a closed loop in total-strain space (the G1 loop): every segment changes the principal axes
_LOOP = [
    np.zeros((3, 3)),
    np.array([[0.0020, 0.0015, 0.0], [0.0015, -0.0010, 0.0010], [0.0, 0.0010, 0.0005]]),
    np.array([[-0.0010, -0.0020, 0.0015], [-0.0020, 0.0015, 0.0], [0.0015, 0.0, -0.0015]]),
    np.array([[0.0005, 0.0, -0.0020], [0.0, 0.0020, 0.0015], [-0.0020, 0.0015, -0.0005]]),
    np.zeros((3, 3)),
]


def _voigt_work(sig6, D):
    """sigma:D with the Voigt (11 22 33 12 23 13) stress and the TENSOR strain increment D."""
    return (sig6[0] * D[0, 0] + sig6[1] * D[1, 1] + sig6[2] * D[2, 2]
            + 2.0 * (sig6[3] * D[0, 1] + sig6[4] * D[1, 2] + sig6[5] * D[0, 2]))


def _hyperelastic_stress(E):
    """(S.5), alpha0 = 0: sigma = p 1 + 2 mu0 dev(E), p = p0 exp(-tr E / kappa_hat); Voigt (11 22 33 12 23 13)."""
    ev = float(np.trace(E))
    p = _P0 * math.exp(-ev / _KH)
    S = p * I3 + 2.0 * _MU * (E - ev / 3.0 * I3)
    return [S[0, 0], S[1, 1], S[2, 2], S[0, 1], S[1, 2], S[0, 2]]


@pytest.mark.parametrize("alpha0", [0.0, 0.2])
def test_k1_2_closed_nonconaxial_elastic_loop_through_a_stdbrick(alpha0):
    """K1.2 through OpenSees.  A closed strain loop 0 -> E1 -> E2 -> E3 -> 0 (every segment non-coaxial, shear on all
    pairs) is prescribed on a stdBrick; the stress is sampled at the 12 Gauss-Legendre nodes of every segment and
    W = sum_j sum_i w_i sigma(E_j + t_i D_j) : D_j is formed.
    EXPECTED (sheet 13.2): |W| <= 1e-12 sum |sigma:D| (zero net work, sigma = dPsi/d eps^e, for alpha0 = 0 and the
    coupled alpha0 = 0.2), and the state returns: stress = -100 I to 1e-12 |p0|, elastic strain 0 to 1e-12 (absolute),
    pi_i and the plastic strains unchanged (every sample elastic: stepInfo plastic == 0).  Non-vacuity (alpha0 = 0, S.5):
    the vertex stresses equal p 1 + 2 mu0 dev(E) at 1e-10 |p|, and the work of the first, OPEN segment equals the stored
    energy Psi(E1) - Psi(0) = -kappa_hat p0 (exp(-eps_v/kappa_hat) - 1) + mu0 e:e at 1e-10 relative: the quadrature is
    sensitive, only the closed loop cancels.
    KILLS: a stress that is not the gradient of a potential (a plastic leak on an elastic path, a non-conservative
    coupling term in the alpha0 stress, an asymmetric tangent-built update), a state that does not return (committed
    strain / pi_i / eps_p drift), a wrong exponent sign or kappa_hat in p(eps_v), a Voigt shear factor in the stress
    <-> tensor map (the vertex closed form and the work both see it)."""
    pi0, v0 = -250.0, 1.7                               # apex pi_c = -538 kPa: the whole loop stays elastic (asserted)
    _deck_prescribed(_PAPER_2INV, v0, pi0, [0.0, 1.0, 2.0, 3.0, 4.0], _LOOP, alpha0=alpha0)
    xg, wg = np.polynomial.legendre.leggauss(12)
    tg, wg = 0.5 * (xg + 1.0), 0.5 * wg
    cur = 0.0

    def advance(to):
        nonlocal cur
        ops.integrator("LoadControl", to - cur)
        assert ops.analyze(1) == 0, f"analyze failed at pseudo-time {to}"
        cur = to
        r = _read()
        assert r["info"][1] == 0.0 and r["info"][0] == 0.0, f"loop left the elastic domain at pseudo-time {to}: {r['info']}"
        return r

    W = Wabs = W_open = 0.0
    for j in range(4):
        D = _LOOP[j + 1] - _LOOP[j]
        for t, w in zip(tg, wg):
            r = advance(j + float(t))
            d = _voigt_work(r["stress"], D)
            W += w * d
            Wabs += w * abs(d)
            if j == 0:
                W_open += w * d
        r = advance(float(j + 1))                                           # the vertex
        if alpha0 == 0.0:
            ref = _hyperelastic_stress(_LOOP[j + 1])
            p_ref = abs(ref[0] + ref[1] + ref[2]) / 3.0
            err = max(abs(a - b) for a, b in zip(r["stress"], ref))
            assert err <= 1e-10 * p_ref, f"vertex {j + 1}: stress {r['stress']} vs closed form {ref}"
    assert Wabs > 0.1, f"the loop does little work ({Wabs:.3e}): the gate would prove nothing"
    assert abs(W) <= 1e-12 * Wabs, f"net work {W:.3e}, sum|sigma:D| {Wabs:.3e}, ratio {abs(W) / Wabs:.2e}"

    # the state returns
    end = _read()
    assert max(abs(a - b) for a, b in zip(end["stress"], [-100.0] * 3 + [0.0] * 3)) <= 1e-12 * abs(_P0), end["stress"]
    assert max(abs(x) for x in _gp("elasticStrain")) <= 1e-12, _gp("elasticStrain")
    assert end["state"][0] == pi0 and end["state"][4] == 0.0 and end["state"][5] == 0.0, end["state"]

    if alpha0 == 0.0:
        E1 = _LOOP[1]
        ev = float(np.trace(E1))
        e_dev = E1 - ev / 3.0 * I3
        psi = -_KH * _P0 * (math.exp(-ev / _KH) - 1.0) + _MU * float(np.sum(e_dev * e_dev))
        assert abs(W_open - psi) <= 1e-10 * abs(psi), (W_open, psi)


# ===========================================================================
#  K1.7  undrained critical-state endpoint, both CSL modes (sheet 13.7)
# ===========================================================================
@pytest.mark.parametrize("case", ["UND_paper", "UND_fork"])
def test_k1_7_undrained_txc_ends_on_the_csl_closed_form(case):
    """Isochoric triaxial compression through a stdBrick (all nodal DOFs prescribed, |eps_ax| = 5) ends at
    p_cs = -exp((v_c0 - v)/lambda_tilde) (paper) or -p_a((e0 - e)/lambda_c)^(1/xi), e = v - 1 (fork), with
    q_cs = M |p_cs| (theta = pi/3, zeta = 1), pi_i = p, psi_i = 0 (S.22), D = 0 (closed form of 5.3).
    EXPECTED (sheet 13.7): v = v0 exp(tr eps) = v0 on an isochoric path (1e-12 relative); the closed form built from v0
    reproduces the target p_cs = -250 kPa the v0 was inverted from (1e-9); p within 1e-6 at |eps_ax| = 5 and the error
    non-increasing at 2, 3.5, 5; q, pi_i/p - 1, psi_i, D within 1e-5 (derivation in the module docstring).
    KILLS: a CSL wired to the wrong branch (paper vs fork) or with a wrong v_c0 / lambda / xi / p_a; a specific volume that
    drifts on an isochoric path (any v-update that does not preserve v on tr d_eps = 0); a hardening law that does not
    drive psi_i to 0 (h, chi, pi_i* wrong); an image pressure that does not approach p; a dilatancy that does not vanish."""
    run = _run(case)
    kw, v0, recs = run["kw"], run["v0"], run["recs"]
    pcs = _pcs_closed(kw, v0)
    assert abs(pcs - _P_CS_TARGET) <= 1e-9 * abs(_P_CS_TARGET), "v0 inversion of the closed form is inconsistent"
    assert abs(recs[-1]["cum"] - 5.0) <= 1e-9
    assert max(abs(r["state"][2] - v0) for r in recs) <= 1e-12 * v0, "specific volume changed on an isochoric path"

    def p_err(at):
        k = int(np.argmin([abs(r["cum"] - at) for r in recs]))
        return abs(_pq_theta(recs[k]["stress"])[0] - pcs) / abs(pcs)

    e2, e35, e5 = p_err(2.0), p_err(3.5), p_err(5.0)
    assert e35 <= e2 + 1e-12 and e5 <= e35 + 1e-12, (e2, e35, e5)
    assert e5 <= 1e-6, f"p_cs relative error at eps_ax = 5: {e5:.3e}"
    end = recs[-1]
    p, q, _ = _pq_theta(end["stress"])
    pi, psi_state, v = end["state"][0], end["state"][1], end["state"][2]
    assert abs(q - _M * abs(pcs)) <= 1e-5 * _M * abs(pcs), (q, _M * abs(pcs))
    assert abs(pi / p - 1.0) <= 1e-5, pi / p
    assert abs(_psi(kw, v, pi)) <= 1e-5
    assert abs(psi_state - _psi(kw, v, pi)) <= 1e-12                   # the state response slot is psi_i (S.22)
    assert abs(_dilatancy_closed(kw, end["stress"], pi, "C")) <= 1e-5


# ===========================================================================
#  K1.8  drained critical-state asymptote (sheet 13.8)
# ===========================================================================
@pytest.mark.parametrize("case", ["TXC_paper", "TXC_fork", "TXE_paper"])
def test_k1_8_drained_critical_state_asymptote_through_a_stdbrick(case):
    """Drained triaxial (symmetric octant, lateral traction -100 kPa, axial displacement control) to |eps_ax| = 3:
    psi_i -> 0, pi_i -> p, eta -> M, i.e. -zeta(theta) q/p -> M so q/|p| -> M/zeta(theta) (M in TXC, rho M in TXE),
    D -> 0 (closed-form dilatancy at the state).  EXPECTED (sheet 13.8) tolerances in the module docstring (G1,
    relaxation e^-10); premise: the lateral stress held at -100 kPa (the drained condition, 1e-6 kPa) at every step; the
    eta error at 3.0 at most half of that at 1.5.
    KILLS: a critical state line not reached (wrong CSL, h or chi); eta saturating below M (a wrong M, a Lode shape
    zeta evaluated at the wrong corner: rho M vs M in TXE); a hardening law that stalls; a drained path that is not
    drained (lateral stress drifting: a wrong shear / volumetric split in the tangent breaks Newton or the constraint)."""
    run = _run(case)
    kw, recs, corner = run["kw"], run["recs"], run["corner"]
    zc = 1.0 if corner == "C" else 1.0 / kw["rho"]
    assert abs(recs[-1]["cum"] - 3.0) <= 1e-9
    assert max(abs(r["stress"][0] - _P0) for r in recs) <= 1e-6, "premise: the lateral stress must stay at -100 kPa"
    assert max(abs(r["stress"][0] - r["stress"][1]) for r in recs) <= 1e-9, "premise: axisymmetric"

    def eta_err(r):
        p, q, _ = _pq_theta(r["stress"])
        return abs(-zc * q / p / _M - 1.0)

    kmid = int(np.argmin([abs(r["cum"] - 1.5) for r in recs]))
    end = recs[-1]
    p, q, _ = _pq_theta(end["stress"])
    pi, v = end["state"][0], end["state"][2]
    assert eta_err(end) <= 2e-3, (eta_err(end), q / abs(p), _M / zc)
    assert eta_err(end) <= 0.5 * eta_err(recs[kmid]), (eta_err(end), eta_err(recs[kmid]))
    assert abs(_psi(kw, v, pi)) <= 1e-3
    assert abs(end["state"][1] - _psi(kw, v, pi)) <= 1e-12             # psi_i slot = (S.22) at the shell's own v, pi_i
    assert abs(pi / p - 1.0) <= 5e-3
    assert abs(_dilatancy_closed(kw, end["stress"], pi, corner)) <= 5e-3


# ===========================================================================
#  K1.9  dissipation census, the `D` response (sheet 13.9)
# ===========================================================================
@pytest.mark.parametrize("case", ["UND_paper", "UND_fork", "TXC_paper", "TXC_fork", "TXE_paper"])
def test_k1_9_dissipation_census_every_step_through_a_stdbrick(case):
    """D^p = Delta lambda sum_a sigma_a q_a >= 0 at every plastic step and = 0 at every elastic step (sheet 13.9, S.38 /
    S.39; condition A holds for these sets), read from the `D` response after every analyze().  Floor: 1e-9 |p0| x the
    axial step (D is stress x strain) on plastic steps, 1e-15 |p0| on elastic ones (G1).  At least 10 plastic steps.
    KILLS: a dissipation with the wrong sign or a missing term in the shell's D response (stress:flow-vector product), a
    plastic step reported as elastic (D == 0 while stepInfo says plastic), a non-conservative flow rule (N_bar > N or
    rho/rho_bar < beta slipping past the parser would give D < 0 near eta ~ M/N)."""
    run = _run(case)
    n_pl = n_el = 0
    for k, r in enumerate(run["recs"]):
        tol = 1e-9 * abs(_P0) * r["step"]
        if r["info"][1] == 1.0:
            n_pl += 1
            assert r["D"] >= -tol, f"{case} step {k + 1} (|eps_ax| = {r['cum']:.4f}): D = {r['D']:.3e}"
        else:
            n_el += 1
            assert abs(r["D"]) <= 1e-15 * abs(_P0), f"{case} elastic step {k + 1}: D = {r['D']:.3e}"
    print(f"\n[{case}] steps {len(run['recs'])}: plastic {n_pl}, elastic {n_el}; "
          f"sum D = {sum(r['D'] for r in run['recs']):.4f}")
    assert n_pl >= 10, "path is mostly elastic: not a dissipation census"


# ===========================================================================
#  K1.5  flow-rule ratio per backward-Euler step (sheet 13.5)
# ===========================================================================
@pytest.mark.parametrize("case", ["TXC_2inv", "TXC_paper", "TXC_fork", "UND_paper", "UND_fork"])
def test_k1_5_flow_rule_ratio_at_every_backward_euler_step(case):
    """Delta eps^p_v / Delta eps^p_s = sqrt(3/2) beta F_p / Omega = sqrt(3/2) (eta - M) / ((1 - N_bar) Omega) at the
    CONVERGED state of the step (AB06 39, sheet 5.3 / 13.5): exact for backward Euler, Delta eps^p = Delta lambda
    q(sigma_{n+1}, pi_{i,n+1}).  The left side is the difference of the shell's committed eps^p_v, eps^p_s
    (`state` response) over the step; the right side is recomputed here from the committed stress and pi_i: eta from
    (S.12), Omega from (S.20) (y_a of S.6, zeta_bar_y(pi/3) the sheet-4.2 number; 2-invariant case: Omega = sqrt(3/2)
    exactly, D = (eta - M)/(1 - N_bar)).  Gate 1e-8 relative to max(1, |D|) (G1; sheet: 2.5e-15 on a full path), every
    step that is plastic, has q > 0, deviatoric plastic increment > 1e-12 and is NOT substepped (stepInfo substeps ==
    1: over a substepped step the ratio is a weighted mean, sheet 13.5).  At least 20 such steps.
    KILLS: a flow rule that is not the sheet's (beta, N_bar, Omega, the corner zeta_bar_y), eps^p_v / eps^p_s accumulated
    with a wrong factor (sqrt(2/3), a missing Omega), plastic strain recovered from the wrong state (trial instead of
    converged), a dilatancy sign error."""
    run = _run(case)
    kw, recs, corner = run["kw"], run["recs"], run["corner"]
    prev = [0.0, 0.0]
    n_checked = n_subst = n_plastic = 0
    worst = 0.0
    for k, r in enumerate(recs):
        epv, eps_s = r["state"][4], r["state"][5]
        dv, ds = epv - prev[0], eps_s - prev[1]
        prev = [epv, eps_s]
        p, q, _ = _pq_theta(r["stress"])
        if r["info"][1] != 1.0 or q <= 0.0:
            continue
        n_plastic += 1
        if r["info"][6] != 1.0:
            n_subst += 1
            continue
        if ds <= 1e-12:
            continue
        d_cf = _dilatancy_closed(kw, r["stress"], r["state"][0], corner)
        err = abs(dv / ds - d_cf) / max(1.0, abs(d_cf))
        worst = max(worst, err)
        n_checked += 1
        assert err <= 1e-8, f"{case} step {k + 1} (|eps_ax| = {r['cum']:.4f}): ratio {dv / ds!r} vs closed form {d_cf!r}"
    print(f"\n[{case}] plastic steps {n_plastic}, substepped (skipped) {n_subst}, checked {n_checked}; worst {worst:.2e}")
    assert n_checked >= 20, f"only {n_checked} unsubstepped plastic steps checked ({n_subst} substepped of {n_plastic})"
