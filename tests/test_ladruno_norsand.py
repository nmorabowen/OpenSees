"""WP-144 gate G2 (Zone A): the MATERIAL-LEVEL contract of `nDMaterial LadrunoNorSand`.

Scope.  This file is the shell contract only -- construction, refusals, state
plumbing, bounded work and the commit latch -- driven through the PUBLIC
interface (openseespy, plus ONE classic-Tcl subprocess smoke).  The mathematics
(return map, tangent vs the oracles, K1/K2 benchmarks) is Zone B:
`Ladruno_files/testbed/norsand_oracle/` (O1 rate oracle, O2 algorithmic oracle,
kernel parity).  Nothing here imports an oracle, scipy or sympy: numpy is not
even needed.

WHERE EVERY EXPECTED VALUE COMES FROM (none is harvested from the shell).
  * Equation sheet `Ladruno_implementation/144a_norsand_equation_sheet.md`:
      (S.5)  BA06 energy, alpha0 = 0:   p = p0 exp(-(eps_v - eps_v0)/kappa_hat),
             s = 2 mu0 e (so sigma12 = mu0 gamma12), K = -p/kappa_hat, G = mu0;
      (S.12) F(p, pi_i) = p eta, eta = (M/N)[1 - (1-N)(p/pi_i)^(N/(1-N))] on the
             hydrostatic axis; apex pi_c = pi_i/(1-N)^((1-N)/N), i.e. the
             surface through an isotropic p is at pi_i = p (1-N)^((1-N)/N);
      (S.22) psi_i = v - v_c0 + lambda_tilde ln(-pi_i)   (paper CSL),
             psi_i = (v-1) - e0 + lambda_c (-pi_i/p_a)^xi  (fork CSL);
      sec.1.2  v = v0 exp(tr eps) (G2 owner decision 2026-10-01, exponential update,
             dv/deps = v; was v0 (1 + tr eps)): v0 is the INITIAL specific volume, a separate
             datum; eps = eps^e + eps^p (additivity); eps^p_s is the accumulated
             sqrt(2/3)||eps^p_dev|| (S.20) -- an EQUALITY on a monotone
             axisymmetric path (the deviatoric flow direction never turns);
      sec.4.1/4.2, 11 (S.39): refusal ranges (GA [7/9,1], WW (1/2,1], both rho
             and rho_bar), dissipation N_bar <= N and rho/rho_bar >= (1-N)/(1-N_bar).
  * Plan sec.2.2 / 2.4 / 2.8 / 2.9 (owner decisions): refuse rho = 1/2, hard
    refusals + warning-only rho > rho_bar, -pi0 REQUIRED, fork CSL needs -p_a,
    the WP-99 commit latch, v0 as separate state (G1 critic N3).
  * AB06 section 6.1 (published): the K2 parameter set used as the base deck.
  * Closed-form arithmetic on those (the numbers are computed in this file by
    numpy-free python, never typed in from a run).

Elements used, and why.  `stdBrick` DISCARDS a material's setTrialStrain return
code (the "discarding element" of the plan); `LadrunoBrick` FORWARDS it (WP-86b).
Both are driven on a ZERO-FREE-DOF unit cube (every node DOF prescribed by `sp`,
u = f(t) E x for a homogeneous strain E), so the answer is the material's own
response to a prescribed strain and no Newton tolerance can enter.  The solver
is `system FullGeneral` throughout: the consistent tangent is NON-symmetric.
`ops.NDTest('SetStrain'|'CommitState'|'GetStress'|'GetResponse', tag, ...)` drives
the prototype material point directly (the setTrialStrain path, no element).

HONEST LIMITS OF THIS FILE (stated, not hidden).
  * NDTest ignores setTrialStrain's return code, so a trial refusal is observed
    through the `refusal` / `stepInfo` responses and the WARNING, and its RETURN
    CODE through a forwarding element (LadrunoBrick -> analyze() != 0).
  * `Domain::addElement` runs ONE `element->update()` at creation, so a clone
    taken from a stepped prototype is immediately integrated to zero strain by
    the element.  Plastic history (pi_i, eps^p) is untouched by that elastic
    unloading and is what the clone tests read; v and the stress are not.
  * sendSelf/recvSelf is reached only through `database File` on a brick deck.
    A trial that differs from the committed state is reached by the clone route
    described above (test_database_roundtrip_trial_differs_from_committed).
  * The tangent is checked against its CLOSED FORM in the elastic regime only;
    the plastic tangent vs finite differences and vs O2 is Zone B (G1 + parity).

WALL TIME (measured on the dev box, Windows, 2026-10-01, CPython 3.12): the whole
file is 9.6-10.4 s for 22 tests; the five wild-increment tests cost 1.8-2.1 s each
(a refused 2^8-substep trial at 8 Gauss points x Newton iterations), the other 17
together ~1 s.  No test is marked slow.
The bounded-work child runs the wild-increment legs under a hard subprocess
timeout so a mutant that removes the work caps fails instead of hanging the
battery.

Every test below names, in a `KILLS:` line, the mutant (a silent regression of
the shell) that would make it fail.
"""
import json
import math
import os
import re
import subprocess
import tempfile

import pytest

# ---- engine binding ---------------------------------------------------------
# Dev box: bind THIS worktree's dist/bin/opensees.pyd (os.add_dll_directory +
# sys.path insert, Ladruno_internal/BUILD_GOTCHAS.md sec.0/4) BEFORE anything
# imports `_testbed`, whose __init__ imports opensees eagerly.  CI: the module is
# copied into tests/ and `_testbed.ops` finds it.
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
_DIST_BIN = os.path.join(_ROOT, "dist", "bin")
if os.path.isfile(os.path.join(_DIST_BIN, "opensees.pyd")):
    from _engine import bind_worktree_engine
    ops = bind_worktree_engine(_DIST_BIN)
else:
    from _testbed import ops
from _testbed.subprocess_run import run_python_script  # noqa: E402

pytestmark = [pytest.mark.zone_a]

_OpsErr = getattr(ops, "OpenSeesError", Exception)


# ===========================================================================
#  The K2 parameter set (AB06 sec.6.1, published) and the deck builders
# ===========================================================================
# kappa_hat 0.01, p0 -100 kPa (eps_v0 0), mu0 5400, alpha0 0, lambda_tilde 0.0135,
# M 1.2, N 0.4, N_bar 0.2, h 280, v 1.59, v_c0 1.81; rho 0.7 / rho_bar 0.8 with
# Willam-Warnke.  chi = -3.5 (BA06/AB06 "alpha ~ -3.5 for sands").  pi_i0 = -60.4
# is a value INSIDE the surface through the isotropic start (checked below against
# the closed form, not assumed).
_P0, _KH, _MU, _M, _N, _NB = -100.0, 0.01, 5400.0, 1.2, 0.4, 0.2
_RHO, _RHOB, _CHI, _H = 0.7, 0.8, -3.5, 280.0
_LT, _VC0, _V0, _PI0 = 0.0135, 1.81, 1.59, -60.4

_BASE = {
    "p0": _P0, "kappa_hat": _KH, "mu0": _MU, "M": _M, "N": _N, "N_bar": _NB,
    "rho": _RHO, "rho_bar": _RHOB, "chi": _CHI, "h": _H,
    "csl": "paper", "lambda_tilde": _LT, "v_c0": _VC0, "v0": _V0, "pi0": _PI0,
}


def _args(**over):
    """Flat `-flag value` list for nDMaterial.  A value of None drops the flag."""
    d = dict(_BASE)
    d.update(over)
    out = []
    for k, v in d.items():
        if v is None:
            continue
        out += ["-" + k, v]
    return out


def _make(tag, **over):
    ops.nDMaterial("LadrunoNorSand", tag, *_args(**over))


def _fresh():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)


# ---- closed forms (python only, from the sheet) ----------------------------
def _F_hydro(p, pi, M=_M, N=_N):
    """(S.12) on the hydrostatic axis (q = 0): F = p * eta."""
    eta = (M / N) * (1.0 - (1.0 - N) * (p / pi) ** (N / (1.0 - N)))
    return p * eta


def _pi_on_surface(p=_P0, N=_N):
    """pi_i of the surface THROUGH an isotropic stress p (eta = 0, sheet 5.1)."""
    return p * (1.0 - N) ** ((1.0 - N) / N)


def _psi_paper(v, pi):
    return v - _VC0 + _LT * math.log(-pi)


def _tensor(e):
    """Voigt (11,22,33,g12,g23,g13; engineering shear) -> 3x3 tensor."""
    return [[e[0], e[3] / 2, e[5] / 2],
            [e[3] / 2, e[1], e[4] / 2],
            [e[5] / 2, e[4] / 2, e[2]]]


def _elastic_closed_form(e):
    """Stress (Voigt) and tangent (6x6, engineering shear) of the BA06 energy
    with alpha0 = 0 from the initial state (sheet S.5): sigma = p 1 + 2 mu0 e_dev,
    p = p0 exp(-eps_v/kappa_hat), tangent = K dd + 2 mu0 (I - dd/3), K = -p/kh."""
    E = _tensor(e)
    ev = E[0][0] + E[1][1] + E[2][2]
    p = _P0 * math.exp(-ev / _KH)
    sig = [[(p if i == j else 0.0) + 2 * _MU * (E[i][j] - (ev / 3 if i == j else 0.0))
            for j in range(3)] for i in range(3)]
    sv = [sig[0][0], sig[1][1], sig[2][2], sig[0][1], sig[1][2], sig[0][2]]
    K = -p / _KH
    T = [[0.0] * 6 for _ in range(6)]
    for i in range(3):
        for j in range(3):
            T[i][j] = K - 2 * _MU / 3 + (2 * _MU if i == j else 0.0)
    for i in range(3, 6):
        T[i][i] = _MU                       # engineering shear: sigma12 = mu0 * gamma12
    return sv, T, p


def _flat(T):
    return [T[i][j] for i in range(6) for j in range(6)]


def _maxabs(a, b):
    return max(abs(x - y) for x, y in zip(a, b))


# ---- deck: a unit cube, every DOF prescribed (zero free equations) ----------
_CUBE = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
         (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]


def _cube(elem, ele_tag, mat_tag, n0=0):
    ids = [n0 + i + 1 for i in range(8)]
    for nid, (x, y, z) in zip(ids, _CUBE):
        ops.node(nid, float(x), float(y), float(z))
    ops.element(elem, ele_tag, *ids, mat_tag)
    return ids


def _prescribe(ids, e, pat, ser):
    """Homogeneous strain: u = E x on every node, every DOF (sp), in pattern `pat`
    driven by time series `ser` (defined by the caller)."""
    E = _tensor(e)
    ops.pattern("Plain", pat, ser)
    for nid, xyz in zip(ids, _CUBE):
        for d in range(3):
            ops.sp(nid, d + 1, float(sum(E[d][k] * xyz[k] for k in range(3))))


def _fix_all(ids):
    for nid in ids:
        ops.fix(nid, 1, 1, 1)


def _analysis(dlam):
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")          # UNSYMMETRIC: the tangent is non-symmetric
    ops.test("NormDispIncr", 1.0e-13, 25, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", dlam)
    ops.analysis("Static")


def _point_deck(elem, e, mat_tag=1, **over):
    """wipe + material + ONE cube + a Linear ramp to the homogeneous strain `e`
    at pseudo-time 1.  Returns the node ids.  Analysis NOT set up."""
    _fresh()
    _make(mat_tag, **over)
    ids = _cube(elem, 1, mat_tag)
    ops.timeSeries("Linear", 1)
    _prescribe(ids, e, 1, 1)
    return ids


def _gp(name, ele=1):
    return list(ops.eleResponse(ele, "material", 1, name))


def _text(capfd):
    c = capfd.readouterr()
    return c.out + c.err


# the shell's refusal vocabulary (LadrunoNorSand.cpp refusalName / subName), as
# the SPEC the warnings are held to: kernel enum Refusal 1..8 and EvalErr 1..8.
_REFUSAL_TEXT = {
    1: "local Newton did not converge",
    2: "local line search failed",
    3: "pi_i root not bracketed",
    4: "pi_i solve did not converge",
    5: "pi_i* guard B <= 0",
    6: "p or pi_i not negative",
    7: "negative plastic multiplier",
    8: "substeps exhausted (2^8)",
}
_SUB_TEXT = {1: "p_or_pi_nonneg", 2: "pi_nonneg", 3: "B_nonpos", 4: "pi_fold",
             5: "pi_nobracket", 6: "pi_noconv", 7: "nonfinite", 8: "singular_J"}


# ===========================================================================
#  1. CONSTRUCTION ECHO AND EVERY REFUSAL
# ===========================================================================
def test_construction_echo_paper_csl(capfd):
    """The material SAYS what it runs, once per construction, with numbers that
    match the sheet's closed forms.

    EXPECTED (sheet S.12, S.22, hand arithmetic): psi_i0 = v0 - v_c0 +
    lambda_tilde ln(-pi_i0); F0 = p0 eta(p0, pi_i0) < 0 (inside).  The echo prints
    6 significant digits, so the numbers are compared to 5e-6 absolute / 5e-5.
    KILLS: an echo that drops a parameter or swaps rho / rho_bar (the line is
    asserted verbatim); a psi_i0 evaluated at p instead of pi_i or with the wrong
    CSL branch; an F0 computed with the wrong eta (N vs N_bar) or sign; an
    "inside/ON/OUTSIDE" classification that is wrong for this start; the
    classTag 33023 echo (the base class tag).
    """
    capfd.readouterr()
    _fresh()
    _make(1)
    out = _text(capfd)
    assert "LadrunoNorSand tag 1" in out and "classTag 33023" in out, out
    assert "p0=-100 kappa_hat=0.01 eps_v0=0 mu0=5400 alpha0=0" in out, out
    assert "M=1.2 N=0.4 N_bar=0.2 rho=0.7 rho_bar=0.8" in out, out
    assert "chi=-3.5 h=280" in out, out
    assert "CSL: paper" in out and "lambda_tilde=0.0135 v_c0=1.81" in out, out
    assert "Lode shape zeta: WW (Willam-Warnke); Q-cap: none" in out, out
    m = re.search(r"pi_i0=(\S+), psi_i0=(\S+), F0=(\S+) \(([^)]*)\)", out)
    assert m, ("no initial-state line in the echo", out)
    pi0, psi0, F0 = float(m.group(1)), float(m.group(2)), float(m.group(3))
    assert pi0 == _PI0
    assert abs(psi0 - _psi_paper(_V0, _PI0)) <= 5e-6, (psi0, _psi_paper(_V0, _PI0))
    F_ref = _F_hydro(_P0, _PI0)
    assert F_ref < 0.0                                    # the start is inside, by the closed form
    assert abs(F0 - F_ref) <= 5e-5 * abs(F_ref), (F0, F_ref)
    assert m.group(4) == "inside the surface", m.group(4)
    # negative control: a baseline start emits no rho > rho_bar warning
    assert "rho > rho_bar" not in out


def test_construction_echo_fork_csl(capfd):
    """Fork CSL echo and the fork psi_i0 closed form (sheet S.22, p_a is a stress).

    EXPECTED: psi_i0 = (v0-1) - e0 + lambda_c (-pi_i0/p_a)^xi.
    KILLS: an echo that prints the paper CSL for the fork deck; a fork psi_i0
    that forgets the -1 (v vs e), drops p_a from the ratio, or swaps xi/lambda_c.
    """
    e0, lc, xi, pa = 0.83, 0.027, 0.45, 101.325
    capfd.readouterr()
    _fresh()
    _make(1, csl="fork", lambda_tilde=None, v_c0=None, e0=e0, lambda_c=lc, xi=xi, p_a=pa)
    out = _text(capfd)
    assert "CSL: fork" in out and "e0=0.83 lambda_c=0.027 xi=0.45 p_a=101.325" in out, out
    m = re.search(r"psi_i0=(\S+), F0=", out)
    assert m, out
    ref = (_V0 - 1.0) - e0 + lc * (-_PI0 / pa) ** xi
    assert abs(float(m.group(1)) - ref) <= 5e-6, (m.group(1), ref)


def test_missing_pi0_is_refused(capfd):
    """-pi0 is REQUIRED (plan 2.8): there is no silent on-the-surface default.

    The refusal must be an error (not a constructed material) and must say which
    flag is missing.  A deck that supplies -pi0 is constructed (positive control,
    same tag, so the first call really did not register one).
    KILLS: a parser that defaults pi_i0 onto the yield surface (the first loading
    step would then be plastic by construction); a refusal that is printed but
    does not abort the command.
    """
    capfd.readouterr()
    _fresh()
    with pytest.raises(_OpsErr):
        _make(1, pi0=None)
    out = _text(capfd)
    assert "missing required -pi0" in out, out
    _make(1)                                   # the same tag is free: nothing was registered
    assert "LadrunoNorSand tag 1" in _text(capfd)


def test_start_outside_the_surface_refused_and_printed_remedy_is_accepted(capfd):
    """A start OUTSIDE the surface is refused; the %.17g remedy it prints, pasted
    back as -pi0, is ACCEPTED (on the surface, with a warning).

    EXPECTED (sheet S.12, apex rule): the surface through an isotropic p0 is at
    pi_i = p0 (1-N)^((1-N)/N); pi_i0 = -20 gives F = p eta > 0 (closed form below).
    Band (one constant, 1e-6 |p0| = 1e-4 here): pi_on*(1 - 1e-3) is outside
    (F ~ +0.2), pi_on*(1 + 1e-3) is inside (F ~ -0.2, NOT 'on').
    KILLS: an outside start accepted (sign of the F > tol test flipped, or the
    check dropped); a remedy computed from the wrong apex rule (compared to the
    closed form to 1e-12); a remedy printed with fewer than 15 significant digits
    (the paste-back would then land off the surface); the on-surface band
    mis-set so the pasted value is refused; 'inside' silently classed as 'ON'.
    """
    pi_on = _pi_on_surface()
    assert _F_hydro(_P0, -20.0) > 1.0e-4                  # outside by the closed form
    assert _F_hydro(_P0, pi_on * (1 - 1e-3)) > 1.0e-4     # just outside the band
    assert _F_hydro(_P0, pi_on * (1 + 1e-3)) < -1.0e-4    # just inside the band

    capfd.readouterr()
    _fresh()
    with pytest.raises(_OpsErr):
        _make(1, pi0=-20.0)
    out = _text(capfd)
    assert "OUTSIDE the yield surface" in out, out
    m = re.search(r"the surface through sigma0 is at pi_i = (\S+):", out)
    assert m, ("the refusal does not print the remedy", out)
    digits = re.sub(r"[^0-9]", "", m.group(1).lstrip("-0.")).rstrip("0")
    assert len(digits) >= 15, ("remedy is not printed with %.17g", m.group(1))
    remedy = float(m.group(1))
    assert abs(remedy - pi_on) <= 1e-12 * abs(pi_on), (remedy, pi_on)

    # pasted back: accepted, on the surface, with the warning
    capfd.readouterr()
    _make(1, pi0=remedy)
    out = _text(capfd)
    assert "ON the yield surface" in out, out
    assert "(ON the surface)" in out, out
    f0 = float(re.search(r"F0=(\S+)", out).group(1))
    assert abs(f0) <= 1.0e-4, f0

    # the band edges
    capfd.readouterr()
    _fresh()
    with pytest.raises(_OpsErr):
        _make(2, pi0=pi_on * (1 - 1e-3))
    assert "OUTSIDE the yield surface" in _text(capfd)
    _make(3, pi0=pi_on * (1 + 1e-3))
    out = _text(capfd)
    assert "(inside the surface)" in out and "ON the yield surface" not in out, out


def test_dissipation_refusal_code_13(capfd):
    """Dissipation rule (sheet S.39, plan 2.2): refuse unless N_bar <= N AND
    rho/rho_bar >= (1-N)/(1-N_bar).  Code 13 is the ratio rule, code 12 N_bar > N.

    EXPECTED: the counterexample N_bar = N = 0.4, rho = 0.7, rho_bar = 0.8 passes
    the paper's rho <= rho_bar yet has ratio 0.875 < beta = 1 -> code 13.  The
    K2 set (beta = 0.75, ratio 0.875) is accepted.  Exact boundary: beta = 0.75,
    rho 0.75 / rho_bar 1.0 (ratio 0.75) accepted; rho 0.74 (ratio 0.74) code 13.
    N_bar = 0.4 > N = 0.2 -> code 12.
    KILLS: a refusal that only checks rho <= rho_bar (the paper's test); the
    inequality flipped or off by one boundary (>= vs >); the ratio compared with
    beta^-1 or with (1-N_bar)/(1-N); N_bar > N not refused; a refusal that does
    not carry its code (the message must say 'code 13').
    """
    capfd.readouterr()
    _fresh()
    with pytest.raises(_OpsErr):
        _make(1, N_bar=0.4, rho=0.7, rho_bar=0.8)
    out = _text(capfd)
    assert "REFUSED (code 13)" in out and "dissipation refusal" in out, out
    assert "rho/rho_bar=0.875" in out, out                # 0.7/0.8, closed form

    _make(2)                                              # the published K2 set
    _make(3, rho=0.75, rho_bar=1.0)                       # ratio == beta exactly: accepted
    capfd.readouterr()
    with pytest.raises(_OpsErr):
        _make(4, rho=0.74, rho_bar=1.0)                   # just below the boundary
    assert "(code 13)" in _text(capfd)

    with pytest.raises(_OpsErr):
        _make(5, N=0.2, N_bar=0.4, rho=0.8, rho_bar=0.8)  # N_bar > N
    out = _text(capfd)
    assert "(code 12)" in out and "N_bar" in out, out


def test_rho_one_half_refused_code_11_willam_warnke(capfd):
    """Willam-Warnke: rho = 1/2 EXACTLY is refused (owner decision 2026-10-01, the
    compression corner becomes a vertex: sheet 4.2), for rho AND rho_bar; the
    range is (1/2, 1].  Code 11.

    EXPECTED (sheet 4.2): 0.5 and 0.4 refused, 0.51 accepted ("admissible but a
    poor choice"), 1.0 accepted, 1.0001 refused.
    KILLS: the open end closed (>= 0.5); the check applied to rho only (the
    sheet: 'every zeta range applies to rho_bar exactly as to rho'); the upper
    end not enforced; a wrong code.
    """
    capfd.readouterr()
    _fresh()
    with pytest.raises(_OpsErr):
        _make(1, rho=0.5, rho_bar=None)               # rho_bar defaults to rho = 0.5
    out = _text(capfd)
    assert "(code 11)" in out and "rho=0.5" in out, out
    with pytest.raises(_OpsErr):
        _make(2, rho=0.7, rho_bar=0.5)                # rho_bar alone
    out = _text(capfd)
    assert "(code 11)" in out and "rho_bar=0.5" in out, out
    with pytest.raises(_OpsErr):
        _make(3, rho=0.4, rho_bar=0.4)
    assert "(code 11)" in _text(capfd)
    with pytest.raises(_OpsErr):
        _make(4, rho=1.0001, rho_bar=1.0001)
    assert "(code 11)" in _text(capfd)
    _make(5, rho=0.51, rho_bar=0.51)
    _make(6, rho=1.0, rho_bar=1.0)
    capfd.readouterr()


def test_gudehus_argyris_rho_below_7_9_refused(capfd):
    """Gudehus-Argyris is convex only for 7/9 <= rho <= 1 (sheet 4.1): refused
    below, for rho and for rho_bar, code 11; 7/9 itself is accepted.  The SAME
    rho = 0.7 is fine under WW (TIMs' c = 0.71 is the reason WW is the default).

    KILLS: GA given the WW range (0.7 accepted); the 7/9 boundary excluded; the
    GA check applied to rho only; WW and GA ranges swapped.
    """
    capfd.readouterr()
    _fresh()
    with pytest.raises(_OpsErr):
        _make(1, zeta="GA", rho=0.7, rho_bar=0.9)             # rho < 7/9
    out = _text(capfd)
    assert "(code 11)" in out and "7/9" in out, out
    with pytest.raises(_OpsErr):
        _make(2, zeta="GA", rho=0.77, rho_bar=0.9)
    assert "(code 11)" in _text(capfd)
    with pytest.raises(_OpsErr):
        _make(3, zeta="GA", rho=0.9, rho_bar=0.7)             # rho_bar < 7/9
    out = _text(capfd)
    assert "(code 11)" in out and "rho_bar=0.69999999" in out, out     # printed with %.17g
    _make(4, zeta="GA", rho=7.0 / 9.0, rho_bar=1.0)           # boundary accepted
    out = _text(capfd)
    assert "Lode shape zeta: GA" in out, out
    _make(5, zeta="WW", rho=0.7, rho_bar=0.9)                 # control: WW takes 0.7
    capfd.readouterr()


def test_rho_above_rho_bar_only_warns(capfd):
    """rho > rho_bar violates only AB06's psi_c <= phi_c reading: a WARNING, not a
    refusal (plan 2.2: dissipation under reading A still holds, ratio 1.125 >=
    beta = 0.75).  The material is constructed.

    KILLS: a hard refusal for rho > rho_bar (the paper's reading re-imposed); a
    warning that is silently dropped.
    """
    capfd.readouterr()
    _fresh()
    _make(1, rho=0.9, rho_bar=0.8)                         # must NOT raise
    out = _text(capfd)
    assert "rho > rho_bar" in out and "only a warning" in out, out
    assert "REFUSED" not in out, out
    assert "LadrunoNorSand tag 1" in out                   # and it echoed: it was built


def test_fork_csl_without_p_a_refused(capfd):
    """The fork CSL multiplies a stress (p_a): no unit-blind default, so -p_a is
    REQUIRED for -csl fork (and not for the paper CSL).

    KILLS: a silent default p_a (e.g. 101.325 or 100) that is wrong in any other
    stress unit; the requirement leaking into the paper CSL.
    """
    capfd.readouterr()
    _fresh()
    fork = dict(csl="fork", lambda_tilde=None, v_c0=None, e0=0.83, lambda_c=0.027, xi=0.45)
    with pytest.raises(_OpsErr):
        _make(1, **fork)                                   # no -p_a
    out = _text(capfd)
    assert "-csl fork needs -p_a" in out, out
    _make(2, p_a=101.325, **fork)                          # with it: accepted
    _make(3)                                               # paper CSL needs no p_a
    capfd.readouterr()


# ===========================================================================
#  2. STATE PLUMBING
# ===========================================================================
def _nd(tag, e):
    ops.NDTest("SetStrain", tag, *e)


def test_elastic_response_and_initial_tangent_match_the_closed_form():
    """Hyperelastic response of the BA06 energy with alpha0 = 0 (sheet S.5), at a
    non-coaxial strain, through the engineering-shear Voigt shell.

    EXPECTED: sigma = p 1 + 2 mu0 e_dev with p = p0 exp(-eps_v/kappa_hat); tangent
    = K dd + 2 mu0 (I - dd/3) (K = -p/kappa_hat) with SHEAR ENTRIES mu0 (engineering
    gamma); initial tangent at the start: C11 = K + 4 mu/3 = 17200, C12 = 6400,
    C44 = 5400.  The strain is elastic by construction (F ~ -30 < 0, asserted by
    stepInfo 'plastic' == 0).  Tolerances: stress 1e-11 |p| and tangent 1e-9 K
    (spectral difference quotient: eps_mach |sigma|/|d eps| ~ 1e-10 kPa).
    KILLS: a factor 2 in the stress<->tensor shear map or in the tangent's shear
    columns; K = -p/kappa_hat replaced by a constant or by p0; a wrong exponent
    sign in p(eps_v); the initial tangent evaluated at the wrong state.
    """
    _fresh()
    _make(1)
    sv0, T0, _ = _elastic_closed_form([0.0] * 6)
    assert _maxabs(ops.NDTest("GetTangentStiffness", 1), _flat(T0)) <= 1e-9 * 17200.0
    assert _flat(T0)[0] == 17200.0 and T0[0][1] == 6400.0 and T0[3][3] == 5400.0
    assert _maxabs(ops.NDTest("GetStress", 1), sv0) == 0.0

    e = [-4.0e-4, -3.0e-4, -2.0e-4, 2.0e-4, -1.5e-4, 1.0e-4]
    _nd(1, e)
    sv, T, p = _elastic_closed_form(e)
    pi_c = _PI0 / (1.0 - _N) ** ((1.0 - _N) / _N)             # compression apex, sheet 5.1: -129.95
    assert pi_c < p < _P0                                     # hand check: p = -109.4, below the apex
    assert ops.NDTest("GetResponse", 1, "stepInfo")[1] == 0.0  # stepInfo.plastic: an elastic step
    assert _maxabs(ops.NDTest("GetStress", 1), sv) <= 1e-11 * abs(p)
    K = -p / _KH
    assert _maxabs(ops.NDTest("GetTangentStiffness", 1), _flat(T)) <= 1e-9 * K


def test_state_response_carries_v_and_v0_separately():
    """v0 != v through the PUBLIC route (plan 2.8, G1 critic N3).

    A monotone axisymmetric (triaxial-compression) path with a volume change is
    driven on the material point (setTrialStrain + commit, two plastic steps);
    the 'state' response [pi_i, psi_i, v, v0, eps_p_v, eps_p_s] must then carry
    v and v0 SEPARATELY.
    EXPECTED (all from the sheet): v0 == 1.59 exactly (it is an input, never
    updated); v = v0 exp(tr eps) (sec.1.2 / 13.10, to 1e-13) and v != v0 by 5e-3;
    psi_i = v - v_c0 + lambda_tilde ln(-pi_i) (S.22, 1e-13); additivity
    eps_p_v = tr eps - eps^e_v with eps^e_v = -kappa_hat ln(p/p0) (S.5, 1e-10);
    on this axisymmetric monotone path the plastic deviatoric flow direction never
    turns, so eps_p_s = sqrt(2/3) ||e - s/(2 mu0)|| exactly (1e-10); pi_i has
    moved off pi_i0 (hardening is live) and the step really was plastic.
    KILLS: v0 following v (v0 := v at commit); the v0 slot returning v (or 0);
    v updated linearly (v += v0 tr de, or v0 (1 + tr eps): off by v0 x^2/2 = 2e-5 here,
    x = tr eps = -5e-3, against the 1e-13 gate); psi_i evaluated with the
    wrong v or pi; eps_p_v dropped / double counted / wrong sign; eps_p_s missing
    its sqrt(2/3) or the Omega factor; pi_i frozen (no hardening).
    """
    _fresh()
    _make(1)
    steps = [[-0.006, 0.0015, 0.0015, 0, 0, 0], [-0.012, 0.0035, 0.0035, 0, 0, 0]]
    for e in steps:
        _nd(1, e)
        assert ops.NDTest("GetResponse", 1, "refusal")[1] == 0.0       # no refused trial
        assert ops.NDTest("GetResponse", 1, "stepInfo")[1] == 1.0      # plastic
        ops.NDTest("CommitState", 1)
    pi, psi, v, v0, epv, eps = ops.NDTest("GetResponse", 1, "state")
    sig = ops.NDTest("GetStress", 1)

    E = _tensor(steps[-1])
    ev = E[0][0] + E[1][1] + E[2][2]
    assert v0 == _V0
    assert abs(v - _V0 * math.exp(ev)) <= 1e-13, (v, _V0 * math.exp(ev))
    assert abs(_V0 * (1.0 + ev) - _V0 * math.exp(ev)) >= 1e-5      # the linear rule is far outside the gate: not vacuous
    assert abs(v - v0) >= 5e-3
    assert abs(psi - _psi_paper(v, pi)) <= 1e-13
    assert pi < _PI0 - 1e-3 and pi < 0.0                              # hardened: |pi_i| grew

    p = (sig[0] + sig[1] + sig[2]) / 3.0
    assert abs(epv - (ev + _KH * math.log(p / _P0))) <= 1e-10, (epv, ev + _KH * math.log(p / _P0))
    s_d = [sig[0] - p, sig[1] - p, sig[2] - p]
    e_d = [E[i][i] - ev / 3.0 for i in range(3)]
    ep_d = [e_d[i] - s_d[i] / (2.0 * _MU) for i in range(3)]
    eps_ref = math.sqrt(2.0 / 3.0) * math.sqrt(sum(x * x for x in ep_d))
    assert eps > 1e-4 and abs(eps - eps_ref) <= 1e-10, (eps, eps_ref)


def _drive_brick(elem, e, n):
    """n LoadControl steps to the homogeneous strain e on a fresh single cube."""
    _point_deck(elem, e)
    _analysis(1.0 / n)
    for k in range(n):
        assert ops.analyze(1) == 0, f"step {k + 1} failed"


_E_PATH = [-0.006, 0.0015, 0.0015, 0.0012, 0.0, 0.0]      # compression + a shear: off the corner


def test_revert_to_start_restores_the_initial_state_including_v0():
    """revertToStart returns to the INITIAL state (sigma0, v0, pi_i0): stress, strain,
    state [pi_i0, psi_i0, v0, v0, 0, 0], the hyperelastic initial tangent, the
    counters -- and a REPLAY of the same path is then bit-identical.

    The deck is the zero-free-DOF stdBrick cube; ops.reset() is revertToStart.
    Before the reset the state must have moved (v != v0, eps_p > 0, pi_i hardened)
    so the post-reset equality is not vacuous; v0 must be 1.59 on both sides.
    KILLS: revertToStart that leaves v (or v0 := v), pi_i, eps_p, the trial strain
    or the refusal counters behind; a revert to the LAST COMMIT instead of the
    start; a stale trial tangent after the revert; non-reproducible replay
    (state left over from the first run).
    """
    _drive_brick("stdBrick", _E_PATH, 10)
    s_run = _gp("stress")
    st = _gp("state")
    assert st[3] == _V0 and abs(st[2] - st[3]) >= 1e-3         # v moved, v0 did not
    assert st[5] > 1e-4 and st[0] < _PI0 - 1e-3                # plastic and hardened
    assert _gp("refusal") == [0.0, 0.0, 0.0, 0.0, 0.0]

    ops.reset()                                                # Domain::revertToStart
    sv0, T0, _ = _elastic_closed_form([0.0] * 6)
    assert _gp("stress") == [-100.0, -100.0, -100.0, 0.0, 0.0, 0.0]
    assert _gp("strain") == [0.0] * 6
    st0 = _gp("state")
    assert st0[0] == _PI0 and st0[2] == _V0 and st0[3] == _V0
    assert st0[4] == 0.0 and st0[5] == 0.0
    assert abs(st0[1] - _psi_paper(_V0, _PI0)) <= 1e-13
    assert _maxabs(_gp("tangent"), _flat(T0)) <= 1e-9 * 17200.0
    assert _gp("refusal") == [0.0, 0.0, 0.0, 0.0, 0.0]

    for k in range(10):                                        # replay
        assert ops.analyze(1) == 0, f"replay step {k + 1} failed"
    assert _gp("stress") == s_run, "replay after revertToStart is not bit-identical"
    assert _gp("state") == st


def _asymmetry(T):
    n = max(abs(x) for x in T)
    return max(abs(T[6 * i + j] - T[6 * j + i]) for i in range(6) for j in range(6)) / n


def _roundtrip_deck(elem, e, prime=None):
    """Fresh cube deck; `prime` (strain or None) first steps the PROTOTYPE with
    NDTest and commits it, so the element's Gauss-point clones are made AFTER it
    stepped."""
    _fresh()
    _make(1)
    if prime is not None:
        _nd(1, prime)
        ops.NDTest("CommitState", 1)
    ids = _cube(elem, 1, 1)
    ops.timeSeries("Linear", 1)
    _prescribe(ids, e, 1, 1)
    _analysis(1.0 / 10.0)


def _snapshot():
    return {k: _gp(k) for k in ("stress", "tangent", "state", "strain", "elasticStrain")}


def _restored_tangent_is_a_tangent_of_the_restored_state(T_after, T_before, p_after):
    """After a restore the tangent is either the saved one bit-exactly (the wire's
    trial block survived) or the hyperelastic closed form at the restored p
    (the element's post-restore update() re-integrated a zero increment)."""
    K = -p_after / _KH
    Tc = [[0.0] * 6 for _ in range(6)]
    for i in range(3):
        for j in range(3):
            Tc[i][j] = K - 2 * _MU / 3 + (2 * _MU if i == j else 0.0)
    for i in range(3, 6):
        Tc[i][i] = _MU
    return T_after == T_before or _maxabs(T_after, _flat(Tc)) <= 1e-9 * K


# WHAT THE DATABASE ROUND TRIP CAN AND CANNOT SHOW (measured, and why).
# `Domain::recvSelf` calls `element->update()` on every element after its recvSelf
# (Domain.cpp ~4273; `Domain::addElement` does the same at :497), which
# re-integrates every Gauss point from its RESTORED committed state to the strain
# of the restored nodal displacements.  For a restored material that is a ZERO
# increment (the strain is bit-identical), so the committed state, the trial
# STATE and every counter are untouched, but the trial STRESS and the trial
# TANGENT are RECOMPUTED (stress within a last bit of the saved one, tangent the
# elastic one rather than the saved consistent plastic one).  Consequently the 42
# trial stress/tangent doubles that sendSelf writes are DEAD DATA from the public
# route: no openseespy call can tell a wire that carries them from one that does
# not.  The tests below therefore hold the COMMITTED state (bit-exact where the
# update cannot move it, and through a continuation that re-integrates from it)
# and the LATCH (which the update does not recompute).  This is reported as a
# finding, not hidden.
def test_database_roundtrip_carries_the_committed_state_with_v0_distinct():
    """sendSelf / recvSelf through `database File` (the SANISAND DB pattern): a
    material restored into a fresh skeleton reports the SAME state (incl. v0 != v
    and a hardened pi_i) and strain BIT-EXACTLY, an elastic strain and a stress
    within rounding of the saved ones (1e-14, see the note above), a tangent that is
    a legitimate tangent of the restored state (the saved CTO bit-exactly, or the
    hyperelastic closed form at the restored p), and then CONTINUES the path to the
    same answer as the never-saved reference (1e-12).

    The tangent at the save is non-symmetric (asserted > 1e-3) and the state has
    v != v0 (>= 1e-3) and eps_p_s > 1e-5, so a v0 written as v, a dropped eps_p, a
    pi_i taken from the initial state or a transposed committed tangent cannot
    pass.  The continuation (4 more plastic steps) re-integrates from the
    RESTORED committed state, which is what the wire is for.
    KILLS: a wire layout that drops or reorders a committed-state slot (eps_e,
    pi_i, v, v0, eps_p_v, eps_p_s, total strain); v0 stored as v; committed state
    not restored (continuation differs from the reference).
    """
    n1, n2 = 6, 4

    def run_to(n):
        _roundtrip_deck("stdBrick", _E_PATH)
        for k in range(n):
            assert ops.analyze(1) == 0

    run_to(n1)                                                  # reference, never saved
    for k in range(n2):
        assert ops.analyze(1) == 0
    ref = _snapshot()

    run_to(n1)
    before = _snapshot()
    assert _asymmetry(before["tangent"]) > 1e-3, "tangent not asymmetric here: nothing to transpose"
    assert before["state"][3] == _V0 and abs(before["state"][2] - _V0) >= 1e-3
    assert before["state"][5] > 1e-5

    with tempfile.TemporaryDirectory(prefix="ladruno_norsand_", ignore_cleanup_errors=True) as td:
        db = os.path.join(td, "norsand_rt")
        ops.database("File", db)
        ops.save(1)
        _roundtrip_deck("stdBrick", _E_PATH)                    # fresh, uncommitted skeleton
        assert _gp("state") != before["state"]
        ops.database("File", db)
        ops.restore(1)
        after = _snapshot()
        for k in ("state", "strain"):
            assert after[k] == before[k], (f"restored '{k}' differs from the saved one", before[k], after[k])
        # the post-restore update re-composes the elastic strain and the stress through
        # the spectral decomposition: equal to rounding, not to the bit
        assert _maxabs(after["elasticStrain"], before["elasticStrain"]) <= 1e-14 * max(abs(x) for x in before["elasticStrain"])
        assert _maxabs(after["stress"], before["stress"]) <= 1e-14 * max(abs(x) for x in before["stress"])
        p_after = sum(after["stress"][:3]) / 3.0
        assert _restored_tangent_is_a_tangent_of_the_restored_state(after["tangent"], before["tangent"], p_after)

        ops.wipeAnalysis()
        _analysis(1.0 / 10.0)
        for k in range(n2):
            assert ops.analyze(1) == 0, f"post-restore step {k + 1} failed"
        cont = _snapshot()
        ops.wipe()
    assert _maxabs(cont["stress"], ref["stress"]) <= 1e-12 * max(abs(x) for x in ref["stress"]), (cont["stress"], ref["stress"])
    assert _maxabs(cont["state"], ref["state"]) <= 1e-12 * max(abs(x) for x in ref["state"])


def test_database_roundtrip_with_committed_differing_from_trial():
    """The same round trip where the clone's COMMITTED state is far from its trial.

    Route: the PROTOTYPE is stepped (setTrialStrain at E1, committed), then the
    element is created; Domain::addElement runs one update() which integrates the
    clone from its copied committed state (strain E1, v = v0 exp(tr E1)) to ZERO
    strain, so the clone's trial (stress ~ -100, v = v0) differs from its
    committed state (E1) -- precondition asserted against the prototype's own
    committed stress.  After save / restore the responses are bit-exact (the
    post-restore update re-integrates the same zero step from the restored
    committed state, so equality here requires that committed state to be
    restored bit-exactly, v included), and the continuation (one plastic step to
    E2, from the restored COMMITTED state, whose v differs from the trial v by
    v0 (1 - exp(tr E1)) = 4.8e-3 in psi_i, 3 % of psi_i) equals the never-saved reference to
    1e-12.
    KILLS: a committed v lost or replaced by the trial v; a committed eps_e,
    pi_i or eps_p not restored; sC and sT confused on the wire.
    """
    e1 = [-0.006, 0.0015, 0.0015, 0.0012, 0.0, 0.0]            # tr = -0.003
    e2 = [-0.012, 0.0035, 0.0035, 0.0010, 0.0, 0.0]

    def deck():
        _roundtrip_deck("stdBrick", e2, prime=e1)

    deck()
    sig_committed = list(ops.NDTest("GetStress", 1))            # prototype committed at E1
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    ref = _snapshot()

    deck()
    before = _snapshot()
    rel = _maxabs(before["stress"], sig_committed) / max(abs(x) for x in sig_committed)
    assert rel > 1e-3, ("trial == committed here: the gate would prove nothing", rel)
    assert before["state"][3] == _V0 and abs(before["state"][2] - _V0) <= 1e-12   # trial v = v0 (zero strain)
    with tempfile.TemporaryDirectory(prefix="ladruno_norsand_", ignore_cleanup_errors=True) as td:
        db = os.path.join(td, "norsand_rt2")
        ops.database("File", db)
        ops.save(1)
        _fresh()                                               # skeleton: element + pattern, NO priming
        _make(1)
        ids = _cube("stdBrick", 1, 1)
        ops.timeSeries("Linear", 1)
        _prescribe(ids, e2, 1, 1)
        _analysis(1.0)
        ops.database("File", db)
        ops.restore(1)
        after = _snapshot()
        for k in before:
            assert after[k] == before[k], (f"restored '{k}' differs", before[k], after[k])
        ops.wipeAnalysis()
        _analysis(1.0)
        assert ops.analyze(1) == 0
        cont = _snapshot()
        ops.wipe()
    assert _maxabs(cont["stress"], ref["stress"]) <= 1e-12 * max(abs(x) for x in ref["stress"]), (cont["stress"], ref["stress"])
    assert _maxabs(cont["state"], ref["state"]) <= 1e-12 * max(abs(x) for x in ref["state"]), (cont["state"], ref["state"])


def _hist(state):
    """The HISTORY part of the state response: pi_i, eps_p_v, eps_p_s."""
    return [state[0], state[4], state[5]]


def _primed_two_clones(e1):
    """Prototype (tag 1) stepped to e1 and committed; THEN two cubes (elements 1, 2)
    are built from it, each on its own nodes.  Returns (ids1, ids2)."""
    _fresh()
    _make(1)
    _nd(1, e1)
    ops.NDTest("CommitState", 1)
    ids1 = _cube("stdBrick", 1, 1, n0=0)
    ids2 = _cube("stdBrick", 2, 1, n0=8)
    return ids1, ids2


def test_getcopy_of_a_stepped_source_carries_its_history_by_value():
    """DOCUMENTED CONTRACT (LadrunoNorSand.h, copyFrom: 'copy EVERYTHING ... the
    clone shares no storage with the source'): a clone taken AFTER the source has
    stepped starts from the source's history, by VALUE.  This is deliberately the
    opposite of LadrunoSANISAND::getCopy (fresh from scalars); it is asserted
    separately so a change of that contract fails HERE and not inside the
    isolation test.

    EXPECTED: the source's (pi_i, eps_p_v, eps_p_s) after one plastic step,
    bit-exactly, in BOTH element clones; the history is non-trivial (eps_p_s > 1e-5).
    KILLS: a getCopy that builds a fresh point (history silently dropped, so an
    element made from a stepped material restarts); a clone that forgets part of
    the history (pi_i reset, eps_p lost).
    """
    e1 = [-0.006, 0.0015, 0.0015, 0.0012, 0.0, 0.0]
    _primed_two_clones(e1)
    h_src = _hist(ops.NDTest("GetResponse", 1, "state"))
    assert h_src[2] > 1e-5 and h_src[0] != _PI0
    assert _hist(_gp("state", 1)) == h_src
    assert _hist(_gp("state", 2)) == h_src


def test_getcopy_clones_are_history_isolated():
    """Clones made AFTER the source stepped share NO history with it or with each
    other (the SANISAND IMPL-EX pattern, done right: the source really has
    history to leak, and each direction is driven).

    Sequence: prototype stepped to E1 (committed) -> elements 1 and 2 cloned from
    it -> (a) the prototype is driven further (a different plastic step,
    committed): both clones' history must still be the E1 history; (b) element 1
    is driven plastically by analysis while element 2 is left at rest (all nodes
    fixed: its Gauss points only see elastic unloading, which cannot change
    history): element 1's history grows, element 2's and the prototype's do not.
    Comparison is exact on pi_i, eps_p_v, eps_p_s (elastic unloading is not
    allowed to move them) and on the prototype's whole state vector.
    KILLS: any static / shared member (a shared State, a shared strain buffer, a
    getCopy that aliases rather than copies): a step of one object would show in
    another.
    """
    e1 = [-0.006, 0.0015, 0.0015, 0.0012, 0.0, 0.0]
    e_proto = [-0.010, 0.0030, 0.0030, 0.0, 0.0, 0.0]
    e2 = [-0.011, 0.0030, 0.0030, 0.0020, 0.0, 0.0]
    ids1, ids2 = _primed_two_clones(e1)
    h1 = _hist(ops.NDTest("GetResponse", 1, "state"))

    # (a) drive the PROTOTYPE after the clones exist
    _nd(1, e_proto)
    ops.NDTest("CommitState", 1)
    state_proto = ops.NDTest("GetResponse", 1, "state")
    assert _hist(state_proto) != h1, "the prototype did not move: (a) proves nothing"
    assert _hist(_gp("state", 1)) == h1
    assert _hist(_gp("state", 2)) == h1

    # (b) drive ELEMENT 1 by analysis; element 2 is fixed (rest)
    ops.timeSeries("Linear", 1)
    _prescribe(ids1, e2, 1, 1)
    _fix_all(ids2)
    _analysis(1.0)
    assert ops.analyze(1) == 0
    h_el1 = _hist(_gp("state", 1))
    assert h_el1[2] > h1[2] + 1e-5, ("element 1 was not driven plastically", h_el1, h1)
    assert _hist(_gp("state", 2)) == h1, "a sibling clone's history changed when element 1 was driven"
    assert ops.NDTest("GetResponse", 1, "state") == state_proto, "the prototype changed when a clone was driven"


# ===========================================================================
#  3. REFUSAL, THE LATCH, AND BOUNDED WORK
# ===========================================================================
_E_WILD = [-0.5, 0.2, 0.2, 0.0, 0.0, 0.0]                  # a finite trial the return map cannot do
_E_BENIGN = 1.0e-3                                         # pseudo-time increment of a benign (elastic) step


def _wild_deck(elem):
    """Cube driven by Linear u = f(t) E_WILD x; analysis at a benign increment."""
    _point_deck(elem, _E_WILD)
    _analysis(_E_BENIGN)


_CHILD = r'''
import json, os, sys, time
d = sys.argv[1]
if d and os.path.isdir(d):
    os.environ["PATH"] = d + os.pathsep + os.environ.get("PATH", "")
    _add = getattr(os, "add_dll_directory", None)
    if _add is not None:
        try:
            _add(d)
        except OSError:
            pass
    sys.path.insert(0, d)
try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops
print("ENGINE=%s" % os.path.abspath(ops.__file__))
args = json.loads(sys.argv[2])
wild = json.loads(sys.argv[3])

def emit(k, v):
    print("%s=%s" % (k, json.dumps(v)))

ops.wipe()
ops.model("basic", "-ndm", 3, "-ndf", 3)
ops.nDMaterial("LadrunoNorSand", 1, *args)
t0 = time.perf_counter()
ops.NDTest("SetStrain", 1, *wild)
emit("POINT_ELAPSED", time.perf_counter() - t0)
emit("POINT_REFUSAL", ops.NDTest("GetResponse", 1, "refusal"))
emit("POINT_STEPINFO", ops.NDTest("GetResponse", 1, "stepInfo"))
emit("POINT_STRESS", ops.NDTest("GetStress", 1))
emit("POINT_STATE", ops.NDTest("GetResponse", 1, "state"))

# the same wild increment through a FORWARDING element (LadrunoBrick), 8 Gauss points
CUBE = [(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]
ops.wipe()
ops.model("basic", "-ndm", 3, "-ndf", 3)
ops.nDMaterial("LadrunoNorSand", 1, *args)
ids = [i + 1 for i in range(8)]
for n, (x, y, z) in zip(ids, CUBE):
    ops.node(n, float(x), float(y), float(z))
ops.element("LadrunoBrick", 1, *ids, 1)
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
E = [[wild[0], wild[3]/2, wild[5]/2], [wild[3]/2, wild[1], wild[4]/2], [wild[5]/2, wild[4]/2, wild[2]]]
for n, xyz in zip(ids, CUBE):
    for dd in range(3):
        ops.sp(n, dd + 1, float(sum(E[dd][k] * xyz[k] for k in range(3))))
ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
ops.test("NormDispIncr", 1.0e-13, 25, 0); ops.algorithm("Newton")
ops.integrator("LoadControl", 1.0); ops.analysis("Static")
t0 = time.perf_counter()
rc = ops.analyze(1)
emit("BRICK_ELAPSED", time.perf_counter() - t0)
emit("BRICK_RC", rc)
emit("BRICK_REFUSAL", list(ops.eleResponse(1, "material", 1, "refusal")))
'''


def test_wild_trial_is_refused_within_a_wall_clock_bound():
    """A wild finite trial increment is REFUSED, and quickly; the committed state
    is untouched; a forwarding element cuts the step.

    Run in a CHILD process under a hard timeout, so a shell whose work caps are
    removed FAILS here (timeout) instead of hanging the whole battery; the child
    prints its engine path (a stale binary would test yesterday's shell).
    EXPECTED:
      * material point (NDTest): refusal response [code, n, latched, finest,
        sub] with code 8 (SUBSTEPS_EXHAUSTED: every refusal went down the 2^8
        ladder, kernel docs), n == 1, latched == 0 (a TRIAL refusal never
        latches), finest in 1..7 (the finest-level cause, never the ladder
        wrapper and never 0); stepInfo.substeps == 256 (2^8);
        stress == the initial isotropic stress and state == the initial state
        BIT-EXACTLY (the committed state is untouched); wall time < 10 s (kernel
        bound: 511 solves x 331 residual evaluations; measured ~0.1 s);
      * LadrunoBrick forwarding: analyze() != 0 (the step is cut), no latch,
        wall time < 60 s (measured ~2 s: 8 Gauss points x Newton iterations).
    KILLS: a missing/raised substep or local-iteration cap (hang -> timeout); a
    refusal that force-accepts; a refused trial that moves the committed state;
    a trial refusal that latches; the sentinel dropped on the way to the element.
    """
    eng = os.path.dirname(os.path.abspath(ops.__file__))
    rc, text = run_python_script(
        _CHILD, argv=(eng, json.dumps(_args()), json.dumps(_E_WILD)),
        merge_stderr=True, timeout=300)
    assert rc == 0, ("the child failed (harness or crash)", text[-3000:])
    eng_line = next((ln for ln in text.splitlines() if ln.startswith("ENGINE=")), None)
    assert eng_line is not None, text[-3000:]
    assert os.path.normcase(os.path.dirname(eng_line.split("=", 1)[1])) == os.path.normcase(eng), (
        "the child loaded a different opensees than this session", eng_line, eng)

    def get(name):
        for ln in text.splitlines():
            if ln.startswith(name + "="):
                return json.loads(ln.split("=", 1)[1])
        raise AssertionError(f"marker {name} missing\n" + text[-3000:])

    ref = get("POINT_REFUSAL")
    assert ref[0] == 8.0 and ref[1] == 1.0 and ref[2] == 0.0, ref
    assert 1.0 <= ref[3] <= 7.0, ref
    info = get("POINT_STEPINFO")
    assert info[6] == 256.0, info                               # 2^8 sub-increments tried
    assert get("POINT_STRESS") == [-100.0, -100.0, -100.0, 0.0, 0.0, 0.0] or \
        _maxabs(get("POINT_STRESS"), [-100.0, -100.0, -100.0, 0.0, 0.0, 0.0]) == 0.0
    st = get("POINT_STATE")
    assert st[0] == _PI0 and st[2] == _V0 and st[3] == _V0 and st[4] == 0.0 and st[5] == 0.0, st
    assert get("POINT_ELAPSED") < 10.0, get("POINT_ELAPSED")

    assert get("BRICK_RC") != 0, "the forwarding element did not cut the step"
    assert get("BRICK_REFUSAL")[2] == 0.0, "a TRIAL refusal latched"
    assert get("BRICK_ELAPSED") < 60.0, get("BRICK_ELAPSED")


def test_forwarding_element_cuts_the_step_and_the_retry_is_untouched(capfd):
    """On a FORWARDING element (LadrunoBrick) a wild trial fails analyze(), WITHOUT
    a latch, and a smaller step then converges to EXACTLY what a run that never
    saw the wild trial gives (the refused trial did not touch the committed
    state).  The WARNING names the return code and the finest cause.

    KILLS: a refused trial that corrupts the committed state (the retry would
    differ); a trial refusal that latches (the retry would fail); the sentinel
    -33086 not in the message; a refusal reported with no finest-level cause.
    """
    # reference: two benign steps, never wild
    _wild_deck("LadrunoBrick")
    for _ in range(2):
        assert ops.analyze(1) == 0
    ref = _gp("stress")

    capfd.readouterr()
    _wild_deck("LadrunoBrick")
    assert ops.analyze(1) == 0                                  # benign step 1 (t = 1e-3)
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) != 0                                  # the wild step is cut
    out = _text(capfd)
    r = _gp("refusal")
    assert r[0] == 8.0 and r[1] >= 1.0 and r[2] == 0.0, r       # refused, counted, NOT latched
    assert "REFUSED this trial strain" in out and "code -33086" in out, out
    fin = int(r[3])
    assert 1 <= fin <= 7
    assert "finest level: " + _REFUSAL_TEXT[fin] in out, (fin, out)
    ops.integrator("LoadControl", _E_BENIGN)
    assert ops.analyze(1) == 0, "the step was not retryable after a trial refusal"
    assert _gp("stress") == ref, "a refused trial moved the committed state"


def test_discarding_element_aborts_the_commit_and_latches_until_revert(capfd):
    """On a DISCARDING element (stdBrick ignores setTrialStrain's return code) a
    refused trial reaches commitState: the commit is ABORTED (analyze() != 0, time
    not advanced), the point LATCHES (every later trial and commit refused, even a
    benign one), and revertToStart clears it -- after which the SAME benign step
    reproduces the pre-wild stress bit-exactly.

    EXPECTED: before the wild step: latched == 0, time = 1e-3; after: analyze !=
    0, latched == 1, time still 1e-3 (the failed step was reverted), stress equal
    to the pre-wild committed stress; a benign retry still fails (latch) and the
    latch keeps the SAME refusal cause; ops.reset() -> refusal response all zero,
    stress = initial; first benign step == the original first benign step.
    KILLS: a commitState that commits the frozen refused state (analyze() == 0 and
    time advances: the silent plateau of the TIMs strip); a missing
    ladrunoNoteCommitRefusal (the Domain does not abort); a latch that does not
    stick (benign retry succeeds); a latch that overwrites the recorded cause;
    revertToStart leaving the latch or the counters.
    """
    capfd.readouterr()
    _wild_deck("stdBrick")
    assert ops.analyze(1) == 0
    s_benign = _gp("stress")
    assert _gp("refusal")[2] == 0.0 and abs(ops.getTime() - _E_BENIGN) <= 1e-15

    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) != 0, "the wild step was committed by a discarding element"
    latched = _gp("refusal")
    assert latched[0] == 8.0 and latched[2] == 1.0, latched
    assert abs(ops.getTime() - _E_BENIGN) <= 1e-15, "the failed step advanced the domain time"
    assert _gp("stress") == s_benign, "the aborted commit changed the stress"

    ops.integrator("LoadControl", _E_BENIGN)
    assert ops.analyze(1) != 0, "a latched point accepted a benign step"
    again = _gp("refusal")
    assert again[2] == 1.0 and again[0] == latched[0] and again[3] == latched[3] and again[4] == latched[4], (
        "the latch changed the recorded refusal cause", latched, again)
    assert again[1] > latched[1]                                # still counted

    ops.reset()
    assert _gp("refusal") == [0.0, 0.0, 0.0, 0.0, 0.0]
    assert _gp("stress") == [-100.0, -100.0, -100.0, 0.0, 0.0, 0.0]
    assert ops.analyze(1) == 0, "revertToStart did not clear the latch"
    assert _gp("stress") == s_benign


def test_latch_survives_a_database_round_trip(capfd):
    """A latched point stays latched, with its recorded cause, through sendSelf /
    recvSelf, and revertToStart then clears it.  (Unlike the trial stress/tangent,
    the latch and the refusal record are NOT recomputed by the post-restore
    update(), so they ARE observable through the wire.)

    EXPECTED: saved [code 8, n, latched 1, finest, sub] -> restored [8, >= n, 1,
    same finest, same sub] (the update() re-integrates a latched point: it is
    refused again and counted, never un-latched and never given a new cause); a
    benign step on the restored point is still refused; ops.reset() clears it and
    the same benign step then succeeds.
    KILLS: a wire that drops the latch flag or the recorded cause (a restored
    latched point would silently resume integrating: the plateau of the TIMs
    strip); recvSelf re-initialising the counters/cause.
    """
    capfd.readouterr()
    _wild_deck("stdBrick")
    assert ops.analyze(1) == 0
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) != 0
    saved = _gp("refusal")
    assert saved[0] == 8.0 and saved[2] == 1.0 and 1.0 <= saved[3] <= 7.0, saved
    with tempfile.TemporaryDirectory(prefix="ladruno_norsand_", ignore_cleanup_errors=True) as td:
        db = os.path.join(td, "norsand_latch")
        ops.database("File", db)
        ops.save(1)
        _wild_deck("stdBrick")                                  # fresh skeleton, never latched
        assert _gp("refusal") == [0.0, 0.0, 0.0, 0.0, 0.0]
        ops.database("File", db)
        ops.restore(1)
        got = _gp("refusal")
        assert got[0] == saved[0] and got[2] == 1.0 and got[3] == saved[3] and got[4] == saved[4], (saved, got)
        assert got[1] >= saved[1], (saved, got)
        ops.wipeAnalysis()
        _analysis(_E_BENIGN)
        assert ops.analyze(1) != 0, "a restored latched point accepted a step"
        ops.reset()
        assert _gp("refusal") == [0.0, 0.0, 0.0, 0.0, 0.0]
        assert ops.analyze(1) == 0, "revertToStart did not clear the restored latch"
        ops.wipe()
    capfd.readouterr()


def test_latch_and_warnings_carry_the_finest_refusal_reason(capfd):
    """The trial WARNING, the latch WARNING and the `refusal` response all carry the
    FINEST-level cause, consistently (plan 2.8 / header: 'both reach the one-time
    WARNING and the refusal / stepInfo responses').

    EXPECTED (shell vocabulary, kernel enums): the refusal response [code,
    count, latched, finest, finest_sub] has code 8 and finest in 1..7; the trial
    WARNING contains 'substeps exhausted (2^8); finest level: <name of finest>'
    (+ ' (<sub>)' when finest_sub != 0); the latch WARNING ('reached commitState',
    'LATCHES') names the same finest cause; stepInfo's finest/finest_sub equal the
    refusal response's.
    KILLS: a warning that prints only the ladder wrapper ('substeps exhausted')
    and loses the cause; a latch warning without the cause; finest / finest_sub
    swapped or zeroed in the response; the latched repeat overwriting the cause.
    """
    capfd.readouterr()
    _wild_deck("stdBrick")
    assert ops.analyze(1) == 0
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) != 0
    out = _text(capfd)
    code, n, latched, fin, sub = _gp("refusal")
    assert code == 8.0 and latched == 1.0 and 1.0 <= fin <= 7.0
    fin, sub = int(fin), int(sub)
    cause = "finest level: " + _REFUSAL_TEXT[fin] + ((" (" + _SUB_TEXT[sub] + ")") if sub else "")
    trial = [ln for ln in out.splitlines() if "REFUSED this trial strain" in ln]
    latch = [ln for ln in out.splitlines() if "reached commitState" in ln and "LATCHES" in ln]
    assert trial, out
    assert latch, out
    assert all("substeps exhausted (2^8)" in ln and cause in ln for ln in trial), (cause, trial[0])
    assert all(cause in ln for ln in latch), (cause, latch[0])
    info = _gp("stepInfo")
    assert int(info[7]) == fin and int(info[8]) == sub, (info, fin, sub)


# ===========================================================================
#  4. Classic Tcl: one deck through OpenSees.exe
# ===========================================================================
def _opensees_exe():
    for root in (_DIST_BIN, os.path.join(_ROOT, "build", "Release")):
        for name in ("OpenSees.exe", "OpenSees"):
            cand = os.path.join(root, name)
            if os.path.isfile(cand):
                return cand
    return None


_TCL_INIT_MISSING = "Can't find a usable init.tcl"


def _find_tcl_library():
    roots = [os.path.join(os.path.expanduser("~"), ".conan2"), os.path.join(_ROOT, "build")]
    for root in roots:
        if not os.path.isdir(root):
            continue
        for dirpath, dirnames, filenames in os.walk(root):
            if "init.tcl" in filenames and "tcl8.6" in dirpath.replace("\\", "/"):
                return dirpath
            dirnames[:] = [d for d in dirnames if "tcl" in d.lower() or d in
                           ("p", "b", "lib", "library", "share", "Release", "res")]
    return None


_TCL_K2 = ("-p0 -100 -kappa_hat 0.01 -mu0 5400 -M 1.2 -N 0.4 -N_bar 0.2 -chi -3.5 -h 280 "
           "-csl paper -lambda_tilde 0.0135 -v_c0 1.81 -v0 1.59")
# Every refusal marker is written AFTER the refusal's message on the same console stream (OpenSees routes both
# opserr and puts there), so the text between two markers is exactly the message of that command.
_TCL_REFUSALS = ("RC_MISSING_PI0", "RC_OUTSIDE", "RC_CODE13", "RC_CODE12", "RC_RHO_HALF", "RC_FORK_NO_PA")
_TCL_DECK = """\
wipe
model basic -ndm 3 -ndf 3
nDMaterial LadrunoNorSand 1 {k2} -rho 0.7 -rho_bar 0.8 -pi0 -60.4
puts "RC_MISSING_PI0=[catch {{nDMaterial LadrunoNorSand 2 {k2} -rho 0.7 -rho_bar 0.8}}]"
puts "RC_OUTSIDE=[catch {{nDMaterial LadrunoNorSand 3 {k2} -rho 0.7 -rho_bar 0.8 -pi0 -20}}]"
puts "RC_CODE13=[catch {{nDMaterial LadrunoNorSand 4 {k2b} -rho 0.7 -rho_bar 0.8 -pi0 -60.4}}]"
puts "RC_CODE12=[catch {{nDMaterial LadrunoNorSand 7 {k2c} -rho 0.8 -rho_bar 0.8 -pi0 -60.4}}]"
puts "RC_RHO_HALF=[catch {{nDMaterial LadrunoNorSand 5 {k2} -rho 0.5 -rho_bar 0.5 -pi0 -60.4}}]"
puts "RC_FORK_NO_PA=[catch {{nDMaterial LadrunoNorSand 6 -p0 -100 -kappa_hat 0.01 -mu0 5400 -M 1.2 -N 0.4 -chi -3.5 -h 280 -csl fork -e0 0.83 -lambda_c 0.027 -xi 0.45 -v0 1.59 -rho 0.7 -pi0 -60.4}}]"
set n 1
foreach {{x y z}} {{0 0 0 1 0 0 1 1 0 0 1 0 0 0 1 1 0 1 1 1 1 0 1 1}} {{ node $n $x $y $z; incr n }}
element stdBrick 1 1 2 3 4 5 6 7 8 1
pattern Plain 1 Linear {{
  set n 1
  foreach {{x y z}} {{0 0 0 1 0 0 1 1 0 0 1 0 0 0 1 1 0 1 1 1 1 0 1 1}} {{
    sp $n 1 [expr {{-1.0e-4*$x}}]
    sp $n 2 [expr {{-1.0e-4*$y}}]
    sp $n 3 [expr {{-1.0e-4*$z}}]
    incr n
  }}
}}
constraints Transformation
numberer Plain
system FullGeneral
test NormDispIncr 1e-13 25 0
algorithm Newton
integrator LoadControl 1.0
analysis Static
puts "RC_ANALYZE=[analyze 1]"
puts "STRESS=[eleResponse 1 material 1 stress]"
puts "STATE=[eleResponse 1 material 1 state]"
puts "LADRUNO-NORSAND-TCL-OK"
"""


def _refusal_segments(out):
    """{marker: (catch result, the console text printed by that command)} from the ordered console output."""
    seg, pos = {}, 0
    for name in _TCL_REFUSALS:
        m = re.search(r"^" + name + r"=(\d+)[ \t]*\r?$", out[pos:], re.M)
        assert m, (f"{name}: marker missing from the console output (or out of order)", out)
        seg[name] = (int(m.group(1)), out[pos:pos + m.start()])
        pos += m.end()
    return seg


def test_tcl_subprocess_smoke(tmp_path):
    """Registration in this fork is FIVE sites and two gates: a material wired into
    openseespy only fails in OpenSees.exe.  One short deck through the classic-Tcl
    binary: construct, SIX refusals (each must be a Tcl error, catch == 1, AND print
    its own reason), a stdBrick cube under an isotropic elastic compression, and the
    closed form.

    EXPECTED refusal messages (the shell vocabulary, sheet S.12 / S.39, 4.2; each
    asserted on the text printed by THAT command, not merely on catch == 1):
      missing -pi0   'missing required -pi0'
      outside start  'OUTSIDE the yield surface of -pi0 = -20', F(sigma0, pi_i0) =
                     p eta(p, pi_i) = 226.323 (closed form S.12 at p = -100, N = 0.4,
                     M = 1.2, printed to 6 digits: 1e-5), and the remedy pi_i =
                     p0 (1-N)^((1-N)/N) = -46.4758001544890 (S.12 apex rule, 1e-12)
      code 13        'REFUSED (code 13)', 'dissipation refusal',
                     rho/rho_bar = 0.875 < (1-N)/(1-N_bar) = 1 (N_bar = N = 0.4)
      code 12        'REFUSED (code 12)', N_bar = 0.4 > N = 0.2
      code 11        'REFUSED (code 11)', rho=0.5 outside the admissible range
                     (1/2, 1] of zeta='WW'
      fork, no p_a   '-csl fork needs -p_a'
    and no refusal prints another's code (a code-13 deck must not say code 11 / 12).
    Elastic closed form (sheet S.5): isotropic eps = -1e-4 per axis (eps_v = -3e-4):
    sigma_ii = p0 exp(0.03) = -103.0455 (deviator 0), 1e-10 relative; v = v0 exp(tr
    eps) = 1.59 exp(-3e-4) (exponential update, sheet 1.2), v0 = 1.59; the echo
    reaches the Tcl console with classTag 33023.
    KILLS: the material registered for openseespy but not for Tcl (unknown
    nDMaterial); the Tcl parser dropping a flag (-pi0, -rho_bar, -csl) while
    accepting the command (the echo and the closed-form stress would be off); a
    refusal that does not surface as a Tcl error; a refusal that surfaces with the
    wrong reason code or a missing message (e.g. every parameter refusal reporting
    code 11, or the code-12 / code-13 rules swapped).
    """
    exe = _opensees_exe()
    if exe is None:
        pytest.skip(f"no OpenSees executable in {_DIST_BIN} or build/Release")
    deck = tmp_path / "norsand_smoke.tcl"
    k2 = _TCL_K2
    deck.write_text(_TCL_DECK.format(k2=k2, k2b=k2.replace("-N_bar 0.2", "-N_bar 0.4"),
                                     k2c=k2.replace("-N 0.4 -N_bar 0.2", "-N 0.2 -N_bar 0.4")))

    def _run(env):
        p = subprocess.run([exe, str(deck)], cwd=str(tmp_path), env=env,
                           stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=180)
        return p, p.stdout + p.stderr

    proc, out = _run(None)
    if _TCL_INIT_MISSING in out:
        tcl_lib = _find_tcl_library()
        if tcl_lib is None:
            pytest.skip("the Tcl runtime is missing for this binary (no init.tcl found); an "
                        "environment gap, see Ladruno_internal/WORKFLOW_GOTCHAS.md sec.7")
        proc, out = _run(dict(os.environ, TCL_LIBRARY=tcl_lib))
        assert _TCL_INIT_MISSING not in out, out

    assert "LADRUNO-NORSAND-TCL-OK" in out, ("the deck did not reach its last line", out)
    assert "unknown nDMaterial" not in out, out
    assert "LadrunoNorSand tag 1" in out and "classTag 33023" in out, out
    assert "psi_i0=-0.164637" in out and "(inside the surface)" in out, out

    # ---- every refusal: a Tcl error AND its own reason text
    seg = _refusal_segments(out)
    for name, (rc, _txt) in seg.items():
        assert rc == 1, (f"{name}: the refusal is not a Tcl error", out)
    t = {k: v[1] for k, v in seg.items()}

    assert "missing required -pi0" in t["RC_MISSING_PI0"], t["RC_MISSING_PI0"]
    assert "REFUSED (code" not in t["RC_MISSING_PI0"], t["RC_MISSING_PI0"]

    o = t["RC_OUTSIDE"]
    assert "OUTSIDE the yield surface of -pi0 = -20" in o, o
    assert "tag 3" in o and "REFUSED (code" not in o, o
    f_print = float(re.search(r"F\(sigma0, pi_i0\) = (\S+) >", o).group(1))
    f_ref = _F_hydro(_P0, -20.0)                                    # closed form S.12 (226.32...)
    assert abs(f_print - f_ref) <= 1e-5 * abs(f_ref), (f_print, f_ref)
    remedy = float(re.search(r"the surface through sigma0 is at pi_i = (\S+):", o).group(1))
    assert abs(remedy - _pi_on_surface()) <= 1e-12 * abs(_pi_on_surface()), (remedy, _pi_on_surface())

    c13 = t["RC_CODE13"]
    assert "tag 4" in c13 and "parameter set REFUSED (code 13)" in c13, c13
    assert "dissipation refusal" in c13 and "rho/rho_bar=0.875" in c13, c13
    assert "(1-N)/(1-N_bar)=1" in c13, c13                           # beta = 0.6/0.6 = 1 for N_bar = N = 0.4
    assert "(code 11)" not in c13 and "(code 12)" not in c13, c13

    c12 = t["RC_CODE12"]
    assert "tag 7" in c12 and "parameter set REFUSED (code 12)" in c12, c12
    assert re.search(r"dissipation refusal: N_bar=0\.4\d* > N=0\.2\d*", c12), c12
    assert "(code 11)" not in c12 and "(code 13)" not in c12, c12

    c11 = t["RC_RHO_HALF"]
    assert "tag 5" in c11 and "parameter set REFUSED (code 11)" in c11, c11
    assert "rho=0.5 outside the admissible range (1/2, 1] of zeta='WW'" in c11, c11
    assert "(code 12)" not in c11 and "(code 13)" not in c11, c11

    fk = t["RC_FORK_NO_PA"]
    assert "-csl fork needs -p_a" in fk and "REFUSED (code" not in fk, fk

    assert "RC_ANALYZE=0" in out, out

    def vec(name):
        m = re.search(name + r"=(.+)", out)
        assert m, (name, out)
        return [float(x) for x in m.group(1).split()]

    p_ref = _P0 * math.exp(3.0e-4 / _KH)
    stress = vec("STRESS")
    assert _maxabs(stress, [p_ref] * 3 + [0.0] * 3) <= 1e-10 * abs(p_ref), (stress, p_ref)
    state = vec("STATE")
    assert state[3] == _V0 and abs(state[2] - _V0 * math.exp(-3.0e-4)) <= 1e-13, state
    assert state[0] == _PI0 and state[4] == 0.0 and state[5] == 0.0, state   # elastic: history untouched


# ===========================================================================
#  4. WP-144 G2 WIRING DECISIONS: getInitialTangent after plastic history, and
#     getCopy of a LATCHED source (appended; the contracts are in
#     LadrunoNorSand.cpp / LEDGER_quirks "getCopy of a latched LadrunoNorSand")
# ===========================================================================
def _linear_step_displacement(prime, initial):
    """ONE `algorithm Linear [-initial]` step on a stdBrick cube (bottom face fixed,
    the four top nodes free), loaded with the element's OWN resisting force R0 at
    its starting state plus a small known increment d: the unbalance of the step is
    exactly d (P - R0), so the displacement u = K^-1 d is the response of whichever
    stiffness K the algorithm assembled (the initial one under -initial, the current
    one otherwise).  Returns u of the four top nodes (12 numbers).  `prime` (a strain,
    or None) first steps the PROTOTYPE with NDTest and commits it, so the element's
    Gauss-point clones are taken from a point that carries plastic history.

    (printA is no use here: it re-forms the CURRENT tangent itself.)"""
    d = [1.0e-6 * (1.0 + 0.37 * i) * (1 if i % 2 == 0 else -1) for i in range(12)]
    _fresh()
    _make(1)
    if prime is not None:
        _nd(1, prime)
        assert ops.NDTest("GetResponse", 1, "stepInfo")[1] == 1.0, "the priming step was not plastic"
        ops.NDTest("CommitState", 1)
    ids = _cube("stdBrick", 1, 1)
    for n in ids[:4]:
        ops.fix(n, 1, 1, 1)
    r0 = list(ops.eleResponse(1, "forces"))[12:]               # resisting force of the 4 free nodes
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for k, n in enumerate(ids[4:]):
        ops.load(n, *[r0[3 * k + c] + d[3 * k + c] for c in range(3)])
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 10, 0)
    ops.algorithm("Linear", "-initial") if initial else ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    return [x for n in ids[4:] for x in ops.nodeDisp(n)]


def test_initial_tangent_is_the_initial_elastic_tangent_after_plastic_steps():
    """`-initial` must stay a genuine initial-stiffness iteration: getInitialTangent
    is the hyperelastic tangent at the INITIAL state (sigma0, v0, pi_i0), however
    much plastic history the point carries (the base default returns getTangent()).

    Observed through the public route, with no oracle: the same known unbalance d
    is solved by one `Linear -initial` step on (a) a cube made from a never-stepped
    material and (b) a cube made from a PLASTICALLY PRIMED prototype (committed
    plastic history in every Gauss point).  Both solve K_init u = d with the same
    K_init, so u must agree; and the same step by plain `Linear` (the CURRENT
    tangent) on the primed cube must NOT agree with it (the history moved the
    elastic tangent: K = -p/kappa_hat at a plastically shifted p), so the equality
    is not vacuous.
    EXPECTED: |u_primed_init - u_fresh| <= 1e-6 max|u_fresh| (the two unbalances differ
    only by the 1e-9 round-off of r0 + d - r0); |u_primed_cur - u_fresh| >= 1e-3
    max|u_fresh| (measured 1.1e-2); at a never-stepped point both arms agree (1e-9).
    KILLS: getInitialTangent returning getTangent() (the primed `-initial` arm would
    then equal the current-tangent arm, off by 1e-2); an initial tangent taken from
    the committed or trial State instead of s0; s0 overwritten by a step, a commit
    or a revert.  (A factor-2 error in the initial tangent's shear columns is the
    closed-form test's job, above.)
    """
    e1 = [-0.006, 0.0015, 0.0015, 0.0012, 0.0, 0.0]            # plastic (hardens pi_i)
    u_fresh = _linear_step_displacement(None, True)
    scale = max(abs(x) for x in u_fresh)
    assert scale > 1e-12
    assert _maxabs(_linear_step_displacement(None, False), u_fresh) <= 1e-9 * scale, (
        "a never-stepped point: the current tangent must BE the initial tangent")
    u_cur = _linear_step_displacement(e1, False)
    assert _maxabs(u_cur, u_fresh) >= 1e-3 * scale, (
        "the plastic history did not move the current tangent: the test would prove nothing")
    u_init = _linear_step_displacement(e1, True)
    assert _maxabs(u_init, u_fresh) <= 1e-6 * scale, (
        "-initial on a plastically primed point is NOT the initial elastic tangent")


def test_getcopy_of_a_latched_source_is_not_latched():
    """DOCUMENTED DECISION (LadrunoNorSand.cpp copyFrom): the refusal latch is
    per integration point and is NOT inherited by a clone.  A fresh element made
    from a latched material must integrate and commit, not start refused.

    Route (public): NDTest drives the PROTOTYPE with a wild increment (refused
    trial) and commits it, which latches the prototype (NDTest ignores the
    return codes, so the commit is the "host that discards the code").  A stdBrick
    cube is then built from it: the element's Gauss points are getCopy clones of
    the latched prototype.  The cube is driven to a mild ELASTIC strain.
    EXPECTED: prototype still latched with its cause kept ([8, 1, 1, finest, sub]);
    the element's `refusal` response all zeros (no refusal, no latch, no count);
    analyze() == 0 for both steps; the Gauss-point stress is the closed-form
    elastic stress of the strain (sheet S.5, 1e-9 |p|), i.e. the clone integrated
    from the source's committed (initial) state.
    KILLS: copyFrom copying `latched` / `trialRefused` / the refusal record (the
    clone would refuse from its first update(), analyze() != 0 and the response
    would show latched 1); a clone that inherits the source's refused trial
    strain; clearing the SOURCE's latch while copying.
    """
    _fresh()
    _make(1)
    _nd(1, _E_WILD)                                            # refused trial (counted, not latched)
    r = list(ops.NDTest("GetResponse", 1, "refusal"))
    assert r[0] == 8.0 and r[1] == 1.0 and r[2] == 0.0, r
    ops.NDTest("CommitState", 1)                               # a refused trial reaches commit: LATCH
    latched = list(ops.NDTest("GetResponse", 1, "refusal"))
    assert latched[0] == 8.0 and latched[2] == 1.0, latched

    e = [-1.2e-4, -1.0e-4, -0.8e-4, 1.0e-4, 0.0, 0.0]          # elastic (F ~ -41 < 0): p = -103, with a shear
    ids = _cube("stdBrick", 1, 1)
    ops.timeSeries("Linear", 1)
    _prescribe(ids, e, 1, 1)
    _analysis(0.5)
    for k in range(2):
        assert ops.analyze(1) == 0, f"step {k + 1}: the clone of a latched material started refused"
    assert _gp("refusal") == [0.0, 0.0, 0.0, 0.0, 0.0], _gp("refusal")
    assert list(ops.NDTest("GetResponse", 1, "refusal")) == latched, "the source's latch/record changed when it was cloned"
    sv, _, p = _elastic_closed_form(e)
    assert _gp("stepInfo")[1] == 0.0, "the mild step was plastic: not the closed-form regime"
    assert _maxabs(_gp("stress"), sv) <= 1e-9 * abs(p), (_gp("stress"), sv)
