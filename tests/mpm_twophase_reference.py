"""WP-232l (232l-V): closed forms and FE (LadrunoUP) references for two-phase particles -- never MPM.

Never reads the MPM C++ (SRC/domain/pattern/ladrunoMPM/, SRC/element/ladrunoMPMPoint/). Every
number comes from a closed form or from the fork's own LadrunoUP (ELE 33017) on the grid, time
step and Newmark pair the gate names (wp-briefs.md "232l", rows T2-T6). Units kN, t, m, s, kPa.
Sign: tension-positive total stress; p compression-positive (LadrunoUP.cpp:2235-2236).

ROUND 1 (this file as committed in Round 1) holds the closed forms and the FE legs that fix the
max(2 x FE, 0.02) bands before any two-phase MPM number. The numpy oracle (T0 i-iii) arrives in
Round 2, written from tests/mpm_twophase_interface.md alone.

    python3.12 tests/mpm_twophase_reference.py closed      # closed-form self-checks (no binary)
    python3.12 tests/mpm_twophase_reference.py terzaghi    # T2 FE: 1x10 and 2x10 Q4, step load
    python3.12 tests/mpm_twophase_reference.py mandel      # T3 FE: 10x4 half strip, tied plate
    python3.12 tests/mpm_twophase_reference.py hydro       # T4 FE: self-weight column, steady
    python3.12 tests/mpm_twophase_reference.py b4          # T5 FE: McGann B4 CB indices
    python3.12 tests/mpm_twophase_reference.py ramp        # T6 FE: T2/T3 cases with the T6 ramp

FE legs need a built fork binary: LADRUNO_BIN=<dir with opensees.pyd>, else this worktree's
dist/bin. Each leg writes tests/data/mpm/twophase/<leg>.json with ops.ladrunoBuild().

TERZAGHI is the case of tests/test_ladruno_up_element_analytic.py:213-224 (E 1e4 kPa, nu 0.2,
rho 2 t/m3, Kf 4444 kPa, n 0.4, rhoF 1, kbar 1e-4 m4/(kN s), H = L = 10 m, q = 10 kPa, 10 cells,
Newmark 0.6 / 0.3025, dt 0.02 s after a 1e-3 s first step, -dynSeepage off). terzaghi_series and
terzaghi_U are that file's :89-102 verbatim; that module binds an engine at import (:74-79), so
Round 1 keeps a copy and tests/test_mpm_twophase_reference.py pins the copy against it.

MANDEL (plane strain, rigid frictionless plate, lateral drainage): the series of Mandel (1953)
in the generalized form of Cheng & Detournay (1988) and Abousleiman et al. (1996), with the root
equation tan a = (1 - nu)/(nu_u - nu) a. It is cross-checked here against an independent
Laplace-domain solution derived from the plane-strain field equations (mandel_p_laplace, fixed
Talbot inversion), so the series does not rest on recalled equation numbers alone.
"""
from __future__ import annotations

import datetime as _dt
import json
import math
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
DATA = HERE / "data" / "mpm" / "twophase"


# ==============================================================================================
# closed forms
# ==============================================================================================
def terzaghi_series(Z, Tv, nterms=300):
    """p/p0 of the Terzaghi series; Z = z/H from the DRAINED face (copy of
    test_ladruno_up_element_analytic.py:89-95)."""
    m = np.arange(nterms)
    M = (2 * m + 1) * np.pi / 2.0
    return np.sum((2.0 / M) * np.sin(M * Z) * np.exp(-M * M * Tv))


def terzaghi_U(Tv, nterms=300):
    """Average degree of consolidation U(Tv) (copy of :98-102)."""
    m = np.arange(nterms)
    M = (2 * m + 1) * np.pi / 2.0
    return 1.0 - np.sum((2.0 / M**2) * np.exp(-M * M * Tv))


def terzaghi_ramp(Z, Tv, Tvr, nterms=300):
    """p/p0 under a linear load ramp 0 -> 1 over Tv in [0, Tvr], then held (Duhamel on the
    series; p0 is the step's undrained pressure per full load). Tvr = 0 is the step."""
    if Tvr <= 0.0:
        return terzaghi_series(Z, Tv, nterms)
    m = np.arange(nterms)
    M = (2 * m + 1) * np.pi / 2.0
    t1 = min(Tv, Tvr)
    g = (np.exp(-M * M * (Tv - t1)) - np.exp(-M * M * Tv)) / (M * M)
    return np.sum((2.0 / M) * np.sin(M * Z) * g) / Tvr


def terzaghi_settle_ramp(Tv, Tvr, s0, sinf, nterms=300):
    """Settlement under the ramp: Duhamel on s_step = s0 + (sinf - s0) U(Tv)."""
    if Tvr <= 0.0:
        return s0 + (sinf - s0) * terzaghi_U(Tv, nterms)
    m = np.arange(nterms)
    M = (2 * m + 1) * np.pi / 2.0
    t1 = min(Tv, Tvr)
    # int_0^t1 [s0 + (sinf-s0)(1 - sum 2/M^2 e^{-M^2 (Tv-tau)})] dtau / Tvr
    g = (np.exp(-M * M * (Tv - t1)) - np.exp(-M * M * Tv)) / (M * M)
    return (sinf * t1 - (sinf - s0) * np.sum((2.0 / M**2) * g)) / Tvr


def terzaghi_constants(E, nu, Kf, poro, kbar, L, q, alpha=1.0):
    """The constants test_ladruno_up_element_analytic.py:218-224 forms (alpha = 1, Ks = inf)."""
    Eoed = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu))
    ooQ = poro / Kf
    Qbar = 1.0 / ooQ
    cv = kbar / (1.0 / Eoed + ooQ)
    p0 = q * Qbar / (Qbar + Eoed)
    s0 = q * L / (Eoed + Qbar)
    sinf = q * L / Eoed
    return dict(Eoed=Eoed, Qbar=Qbar, cv=cv, p0=p0, s0=s0, sinf=sinf)


# ---------------------------------------------------------------------------------- Mandel
def mandel_constants(E, nu, Kf, poro, kbar, alpha=1.0, Ks=None):
    """Biot constants for Mandel. Storage S = n/Kf + (alpha - n)/Ks (LadrunoUP's 1/Qbar,
    ADR-71 §3.1); Ks = None means incompressible grains (LadrunoUP's default -Ks -1)."""
    K = E / (3.0 * (1.0 - 2.0 * nu))
    G = E / (2.0 * (1.0 + nu))
    S = poro / Kf + (0.0 if Ks is None else (alpha - poro) / Ks)
    M = 1.0 / S
    Ku = K + alpha * alpha * M
    B = alpha * M / Ku
    nu_u = (3.0 * Ku - 2.0 * G) / (2.0 * (3.0 * Ku + G))
    c = kbar * M * (K + 4.0 * G / 3.0) / (Ku + 4.0 * G / 3.0)
    c_cd = (2.0 * kbar * B * B * G * (1 - nu) * (1 + nu_u) ** 2
            / (9.0 * (1 - nu_u) * (nu_u - nu)))       # the papers' form; equals c when alpha = 1
    return dict(K=K, G=G, S=S, M=M, Ku=Ku, B=B, nu=nu, nu_u=nu_u, c=c, c_cd=c_cd,
                alpha=alpha, kbar=kbar)


def mandel_roots(nu, nu_u, n=200):
    """First n positive roots of tan(a) = (1 - nu)/(nu_u - nu) * a; root k lies in
    (k pi, k pi + pi/2) (k = 0, 1, ...) because the slope ratio exceeds 1."""
    from scipy.optimize import brentq
    r = (1.0 - nu) / (nu_u - nu)
    f = lambda a: math.sin(a) - r * a * math.cos(a)   # noqa: E731  (no poles)
    out = []
    eps = 1e-12
    for k in range(n):
        lo, hi = k * math.pi + eps, k * math.pi + math.pi / 2.0 - eps
        if k == 0:
            lo = 1e-9
        out.append(brentq(f, lo, hi, xtol=1e-15, rtol=1e-15, maxiter=200))
    return np.array(out)


def mandel_p(x, t, a, sig0, mc, nterms=200, t_ramp=0.0, roots=None):
    """Mandel pore pressure p(x, t), plate force 2F = 2 a sig0 on the full width 2a
    (half model: force a*sig0 on 0 <= x <= a), drained at |x| = a. Series:
        p = (2 F B (1+nu_u) / (3a)) sum sin a_n/(a_n - sin a_n cos a_n)
                                       (cos(a_n x/a) - cos a_n) exp(-a_n^2 c t / a^2)
    with F = a sig0. t_ramp > 0: linear load ramp over [0, t_ramp] then held (Duhamel)."""
    al = mandel_roots(mc["nu"], mc["nu_u"], nterms) if roots is None else roots
    F = a * sig0
    pref = 2.0 * F * mc["B"] * (1.0 + mc["nu_u"]) / (3.0 * a)
    coef = np.sin(al) / (al - np.sin(al) * np.cos(al)) * (np.cos(al * x / a) - np.cos(al))
    lam = al * al * mc["c"] / (a * a)
    if t_ramp <= 0.0:
        return pref * float(np.sum(coef * np.exp(-lam * t)))
    t1 = min(t, t_ramp)
    g = (np.exp(-lam * (t - t1)) - np.exp(-lam * t)) / lam
    return pref * float(np.sum(coef * g)) / t_ramp


def mandel_p0(sig0, mc):
    """Undrained (t = 0+) pressure: Skempton under plane strain, B (1 + nu_u) sig0 / 3."""
    return mc["B"] * (1.0 + mc["nu_u"]) * sig0 / 3.0


def mandel_p_laplace(x, t, a, sig0, mc, nT=32):
    """INDEPENDENT check of mandel_p (step load), derived here from the plane-strain field
    equations, not from the papers. sigma_xx = 0 everywhere (lateral faces free, sigma_xy = 0),
    eps_yy = eps(t) uniform (rigid frictionless plate), p(+-a) = 0:
        eps_xx = (alpha p - lam eps)/(lam + 2G)
        c1 p_t - kbar p_xx = -c2 eps_t,  c1 = S + alpha^2/(lam+2G),  c2 = 2 G alpha/(lam+2G)
        2a [4G(lam+G)/(lam+2G)] eps - c2 int p dx = -2 F      (plate equilibrium)
    Laplace: p^ = B^ (1 - cosh(q x)/cosh(q a)), q^2 = c1 s/kbar,
        B^ [1 + c2 alpha (1 - tanh(qa)/(qa)) / (2 c1 (lam+G))] = (c2/c1) sig0 (lam+2G)/(4G(lam+G) s)
    inverted by the fixed Talbot rule (Abate & Valko 2004), nT nodes."""
    G, K, S, al = mc["G"], mc["K"], mc["S"], mc["alpha"]
    lam = K - 2.0 * G / 3.0
    c1 = S + al * al / (lam + 2 * G)
    c2 = 2 * G * al / (lam + 2 * G)

    def phat(s):
        q = np.sqrt(c1 * s / mc["kbar"])
        qa = q * a
        T = np.tanh(qa) / qa
        Bh = (c2 / c1) * sig0 * (lam + 2 * G) / (4 * G * (lam + G) * s)
        Bh = Bh / (1.0 + c2 * al * (1.0 - T) / (2.0 * c1 * (lam + G)))
        # 1 - cosh(qx)/cosh(qa), overflow-safe
        ratio = np.exp(q * (x - a)) * (1 + np.exp(-2 * q * x)) / (1 + np.exp(-2 * qa))
        return Bh * (1.0 - ratio)

    r = 2.0 * nT / (5.0 * t)
    th = np.arange(1, nT) * np.pi / nT
    cot = np.cos(th) / np.sin(th)
    s = r * th * (cot + 1j)
    sig = th + (th * cot - 1.0) * cot
    tot = 0.5 * phat(complex(r, 0.0)) * np.exp(r * t)
    tot = tot + np.sum(np.exp(t * s) * np.array([phat(si) for si in s]) * (1 + 1j * sig))
    return float(np.real(tot) * r / nT)


def mandel_peak(mc, a=1.0, nt=4000):
    """Closed-form Mandel-Cryer peak of p(0,t)/p0 and its time factor c t/a^2."""
    al = mandel_roots(mc["nu"], mc["nu_u"], 200)
    p0 = mandel_p0(1.0, mc)
    tv = np.geomspace(1e-5, 1.0, nt)
    pv = np.array([mandel_p(0.0, T * a * a / mc["c"], a, 1.0, mc, roots=al) for T in tv]) / p0
    k = int(np.argmax(pv))
    return float(pv[k]), float(tv[k])


# ------------------------------------------------------------------------------ hydrostatics
def hydrostatic(Y, H, rho_sat, rho_f, g=9.81):
    """Steady self-weight column, drained top at Y = H, p = 0 there:
    p = rho_f g (H - Y); sigma'_yy = -(rho_sat - rho_f) g (H - Y) (tension-positive)."""
    return rho_f * g * (H - Y), -(rho_sat - rho_f) * g * (H - Y)


def closed_selfcheck():
    """Closed-form self-checks: Mandel series vs the Laplace solution; p(x, 0+) uniform;
    terzaghi_ramp -> step as Tvr -> 0. Returns a dict of measured differences."""
    out = {}
    mc = mandel_constants(**MANDEL_SOIL)
    out["mandel_c_vs_cd_form"] = abs(mc["c"] - mc["c_cd"]) / mc["c"]
    a = MANDEL["a"]
    p0 = mandel_p0(MANDEL["sig0"], mc)
    worst = 0.0
    for tv in (1e-3, 0.01, 0.05, 0.1, 0.3, 1.0):
        t = tv * a * a / mc["c"]
        for x in (0.0, 0.3 * a, 0.7 * a, 0.95 * a):
            ps = mandel_p(x, t, a, MANDEL["sig0"], mc)
            pl = mandel_p_laplace(x, t, a, MANDEL["sig0"], mc)
            worst = max(worst, abs(ps - pl) / p0)
    out["mandel_series_vs_laplace_max_over_p0"] = worst
    t = 1e-9 * a * a / mc["c"]
    out["mandel_p_0plus_uniform_err"] = max(
        abs(mandel_p(x, t, a, MANDEL["sig0"], mc, nterms=4000) / p0 - 1.0) for x in (0.0, 0.5 * a))
    out["terzaghi_ramp_to_step"] = abs(terzaghi_ramp(0.4, 0.2, 1e-6) - terzaghi_series(0.4, 0.2))
    return out


# ==============================================================================================
# the cases (Round 1 proposals; the orchestrator freezes them)
# ==============================================================================================
TERZ = dict(E=1.0e4, nu=0.2, rho=2.0, Kf=4444.0, poro=0.4, rhoF=1.0, kbar=1.0e-4,
            L=10.0, q=10.0, nely=10, dt=0.02, dt0=1.0e-3, gamma=0.6, beta=0.3025)
TV_T2 = (0.05, 0.1, 0.2, 0.5)

# Mandel: Terzaghi's soil with nu 0.1 (was 0.2): closed-form peak 1.0791 p0 (step), 1.0781 p0
# under the 5 s ramp. nu 0.2 gives only 1.059. A STEP load is not gated: the tied strip rings
# (undrained lateral transit a/V_u = 0.097 s; gamma 0.6 barely damps the slow modes at dt 0.02)
# and the FE misses p(0,t) by 0.05-1.1 p0 before tv 0.07. The ramp (5 s = 52 lateral transits,
# tv_r = 0.0266) is the gated load for T3 and T6 alike.
MANDEL_SOIL = dict(E=1.0e4, nu=0.1, Kf=4444.0, poro=0.4, kbar=1.0e-4)
MANDEL = dict(a=10.0, b=4.0, nx=10, ny=4, rho=2.0, rhoF=1.0, sig0=10.0, t_ramp=5.0,
              dt=0.02, dt0=1.0e-3, gamma=0.6, beta=0.3025)
MANDEL_TVS = (0.04, 0.078, 0.12, 0.2, 0.4)     # 0.078 = the ramp closed form's peak (tv 0.0781)
MANDEL_NOTE = ("Series of Mandel (1953) in the form of Cheng & Detournay (1988) and Abousleiman et "
               "al. (1996), root equation tan a = (1 - nu)/(nu_u - nu) a. The papers' EQUATION "
               "NUMBERS ARE UNVERIFIED (the papers were not opened). The warrant is the independent "
               "Laplace-domain solution derived in mandel_p_laplace (fixed Talbot inversion), which "
               "the series matches to selfcheck.mandel_series_vs_laplace_max_over_p0.")
TERZ_T6_RAMP = 5.0   # s: 53 undrained transits of H (0.0949 s each); 2 s (21) leaves 0.0175 at Tv 0.1, dt 0.002

HYDRO = dict(E=1.0e4, nu=0.2, rho=2.0, Kf=4444.0, poro=0.4, rhoF=1.0, kbar=1.0e-4,
             L=10.0, nely=10, g=9.81)

B4 = dict(E=25000.0, nu=0.3, rho=2.67, poro=0.4, Kf=2.2e12, rhoF=1.0,
          perm=1.0e-7, Ncell=10, L=30.0, q=0.1, strip=3.0, dt=0.01, T=1.0)   # equiv.py:332-334


# ==============================================================================================
# FE engine
# ==============================================================================================
_OPS = None


def fe_engine():
    global _OPS
    if _OPS is not None:
        return _OPS
    b = os.environ.get("LADRUNO_BIN") or str(REPO / "dist" / "bin")
    if not (Path(b) / "opensees.pyd").exists() and not list(Path(b).glob("opensees*.so")):
        raise SystemExit(f"no fork binary at {b} (set LADRUNO_BIN); FE legs need one")
    os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")
    os.environ["PATH"] = b + os.pathsep + os.environ.get("PATH", "")
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(b)
    sys.path.insert(0, b)
    import opensees as ops  # noqa: E402
    if os.path.normcase(os.path.dirname(os.path.abspath(ops.__file__))) != os.path.normcase(b):
        raise SystemExit(f"WRONG ENGINE: imported {ops.__file__}, expected {b}")
    _OPS = ops
    return ops


def _provenance(leg: str, deck: str) -> dict:
    ops = fe_engine()
    return dict(leg=leg, build=ops.ladrunoBuild(), engine=os.path.abspath(ops.__file__),
                deck=f"tests/mpm_twophase_reference.py::{deck}",
                date=_dt.date.today().isoformat(),
                units="kN, t, m, s, kPa; sigma tension-positive, p compression-positive")


def _save(name: str, payload: dict) -> Path:
    DATA.mkdir(parents=True, exist_ok=True)
    p = DATA / f"{name}.json"
    payload["output"] = str(p.relative_to(REPO)).replace("\\", "/")
    p.write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    return p


def _stab_args(stab):
    if stab == "off":
        return ["-stab", "off"]
    if isinstance(stab, tuple) and stab[0] == "auto":
        return ["-stab", "auto", float(stab[1])]
    return ["-stab", float(stab)]


def _transient(gamma, beta, constraints=("Transformation",)):
    ops = fe_engine()
    ops.system("UmfPack")
    ops.numberer("RCM")
    ops.constraints(*constraints)
    ops.algorithm("Linear")
    ops.integrator("Newmark", gamma, beta)
    ops.analysis("Transient")


def _quad_p(pn, xi, eta):
    """Bilinear interpolation of the 4 corner values (ccw from (-1,-1)) at (xi, eta)."""
    N = np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                  (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)]) / 4.0
    return float(N @ pn)


# ---------------------------------------------------------------------------------------- T2
def fe_terzaghi(nx=1, stab="off", t_ramp=0.0, first_small=True, dt=None, tvs=TV_T2,
                plate=True):
    """LadrunoUP Q4 nx x 10 column (h = 1 m): ux rollers on x = 0 and x = nx (every node when
    nx = 1), uy = 0 at the base, p = 0 on the top row, impervious elsewhere; top load q through
    the frozen interface's rigid plate (P10: top-row u_y tied by equalDOF, force q nx h on the
    retained node) or, plate=False, as consistent nodal forces (identical on this 1D column). Step (Constant series) or
    a linear ramp over t_ramp. Errors at Tv: p at nodes (interior j = 1..9, as the FE test; and
    j = 0..9), p at the GIMP ppc-2 particle points (y = (j + 1/4) h, (j + 3/4) h) and at the
    2x2 Gauss points (y = (j + 1/2 +- 1/(2 sqrt 3)) h), and |s - s_exact|/(sinf - s0)."""
    ops = fe_engine()
    c = TERZ
    dt = c["dt"] if dt is None else dt
    k = terzaghi_constants(c["E"], c["nu"], c["Kf"], c["poro"], c["kbar"], c["L"], c["q"])
    L, ny = c["L"], c["nely"]
    h = L / ny
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    nid = {}
    n = 1
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.node(n, i * h, j * h)
            nid[(i, j)] = n
            n += 1
    ops.nDMaterial("ElasticIsotropic", 1, c["E"], c["nu"], c["rho"])
    e = 1
    for j in range(ny):
        for i in range(nx):
            ops.element("LadrunoUP", e, nid[(i, j)], nid[(i + 1, j)], nid[(i + 1, j + 1)],
                        nid[(i, j + 1)], 1, "-Kf", c["Kf"], "-poro", c["poro"],
                        "-rhoF", c["rhoF"], "-perm", c["kbar"], c["kbar"],
                        "-dynSeepage", "off", *_stab_args(stab))
            e += 1
    for j in range(ny + 1):
        for i in range(nx + 1):
            fx = 1 if i in (0, nx) else 0
            fy = 1 if j == 0 else 0
            fp = 1 if j == ny else 0
            if fx or fy or fp:
                ops.fix(nid[(i, j)], fx, fy, fp)
    if t_ramp > 0.0:
        ops.timeSeries("Path", 1, "-time", 0.0, t_ramp, 1.0e12, "-values", 0.0, 1.0, 1.0)
    else:
        ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    if plate:   # rigid frictionless plate: tie the top row's u_y, total force q * nx * h
        for i in range(1, nx + 1):
            ops.equalDOF(nid[(0, ny)], nid[(i, ny)], 2)
        ops.load(nid[(0, ny)], 0.0, -c["q"] * h * nx, 0.0)
    else:
        for i in range(nx + 1):
            w = 0.5 if i in (0, nx) else 1.0
            ops.load(nid[(i, ny)], 0.0, -c["q"] * h * w, 0.0)
    _transient(c["gamma"], c["beta"])
    t = 0.0
    if first_small:
        assert ops.analyze(1, c["dt0"]) == 0
        t = c["dt0"]
    Tvr = t_ramp * k["cv"] / L**2
    rows = []
    for Tv in tvs:
        tt = Tv * L * L / k["cv"]
        nst = int(round((tt - t) / dt))
        if nst > 0:
            assert ops.analyze(nst, dt) == 0
            t += nst * dt
        TvA = t * k["cv"] / (L * L)
        pn = np.array([ops.nodeDisp(nid[(0, j)], 3) for j in range(ny + 1)])
        pex = lambda y: k["p0"] * terzaghi_ramp((L - y) / L, TvA, Tvr)   # noqa: E731
        node_err = [abs(pn[j] - pex(j * h)) / k["p0"] for j in range(ny + 1)]
        part, gaus = [], []
        for j in range(ny):
            for f in (0.25, 0.75):
                y = (j + f) * h
                pp = pn[j] * (1 - f) + pn[j + 1] * f
                part.append(abs(pp - pex(y)) / k["p0"])
            for f in (0.5 - 0.5 / math.sqrt(3), 0.5 + 0.5 / math.sqrt(3)):
                y = (j + f) * h
                pp = pn[j] * (1 - f) + pn[j + 1] * f
                gaus.append(abs(pp - pex(y)) / k["p0"])
        s = -ops.nodeDisp(nid[(0, ny)], 2)
        s_ex = terzaghi_settle_ramp(TvA, Tvr, k["s0"], k["sinf"])
        rows.append(dict(Tv=Tv, Tv_actual=TvA, t=t,
                         p_err_nodes_interior=max(node_err[1:ny]),
                         p_err_nodes_all=max(node_err[:ny]),
                         p_err_at=int(np.argmax(node_err[:ny])),
                         p_err_gimp_points=max(part), p_err_gauss_points=max(gaus),
                         U_err=abs(s - s_ex) / (k["sinf"] - k["s0"]),
                         settle=s, settle_exact=s_ex,
                         U_fem=(s - k["s0"]) / (k["sinf"] - k["s0"]),
                         U_exact=(s_ex - k["s0"]) / (k["sinf"] - k["s0"])))
    return dict(nx=nx, stab=stab if isinstance(stab, str) else list(stab), t_ramp=t_ramp,
                first_small_step=first_small, dt=dt, plate=plate, constants=k, rows=rows)


# ---------------------------------------------------------------------------------------- T3
def fe_mandel(stab="off", t_ramp=0.0, tvs=None, nx=None, ny=None, dt=None, first_small=True):
    """LadrunoUP Q4 half strip 0 <= x <= a, 0 <= y <= b (h = 1 m): ux = 0 and impervious on x = 0
    (symmetry), uy = 0 and impervious on y = 0 (symmetry), p = 0 and traction-free on x = a
    (drained face), the top row's uy tied to its left node by equalDOF (the rigid frictionless
    plate) and the plate force a*sig0 on that node. Observable p(0, t)/p0 at the left column."""
    ops = fe_engine()
    c = MANDEL
    nx = c["nx"] if nx is None else nx
    ny = c["ny"] if ny is None else ny
    dt = c["dt"] if dt is None else dt
    mc = mandel_constants(**MANDEL_SOIL)
    a, b = c["a"], c["b"]
    hx, hy = a / nx, b / ny
    p0 = mandel_p0(c["sig0"], mc)
    al = mandel_roots(mc["nu"], mc["nu_u"], 200)
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    nid = {}
    n = 1
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.node(n, i * hx, j * hy)
            nid[(i, j)] = n
            n += 1
    s = MANDEL_SOIL
    ops.nDMaterial("ElasticIsotropic", 1, s["E"], s["nu"], c["rho"])
    e = 1
    for j in range(ny):
        for i in range(nx):
            ops.element("LadrunoUP", e, nid[(i, j)], nid[(i + 1, j)], nid[(i + 1, j + 1)],
                        nid[(i, j + 1)], 1, "-Kf", s["Kf"], "-poro", s["poro"],
                        "-rhoF", c["rhoF"], "-perm", s["kbar"], s["kbar"],
                        "-dynSeepage", "off", *_stab_args(stab))
            e += 1
    for j in range(ny + 1):
        for i in range(nx + 1):
            fx = 1 if i == 0 else 0
            fy = 1 if j == 0 else 0
            fp = 1 if i == nx else 0
            if fx or fy or fp:
                ops.fix(nid[(i, j)], fx, fy, fp)
    master = nid[(0, ny)]
    for i in range(1, nx + 1):
        ops.equalDOF(master, nid[(i, ny)], 2)
    if t_ramp > 0.0:
        ops.timeSeries("Path", 1, "-time", 0.0, t_ramp, 1.0e12, "-values", 0.0, 1.0, 1.0)
    else:
        ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(master, 0.0, -a * c["sig0"], 0.0)
    _transient(c["gamma"], c["beta"])
    t = 0.0
    if first_small:
        assert ops.analyze(1, c["dt0"]) == 0
        t = c["dt0"]
    rows = []
    for tv in tvs:
        tt = tv * a * a / mc["c"]
        nst = int(round((tt - t) / dt))
        if nst > 0:
            assert ops.analyze(nst, dt) == 0
            t += nst * dt
        tvA = t * mc["c"] / (a * a)
        pc = np.mean([ops.nodeDisp(nid[(0, j)], 3) for j in range(ny + 1)])
        pc_spread = np.ptp([ops.nodeDisp(nid[(0, j)], 3) for j in range(ny + 1)])
        pex = mandel_p(0.0, t, a, c["sig0"], mc, t_ramp=t_ramp, roots=al)
        # the MPM observable: the first particle column (ppc 2: x = hx/4), FE interpolated
        p1 = np.mean([ops.nodeDisp(nid[(1, j)], 3) for j in range(ny + 1)])
        pq = 0.75 * pc + 0.25 * p1
        pexq = mandel_p(0.25 * hx, t, a, c["sig0"], mc, t_ramp=t_ramp, roots=al)
        rows.append(dict(tv=tv, tv_actual=tvA, t=t, p_fe_over_p0=pc / p0,
                         p_exact_over_p0=pex / p0, err_over_p0=abs(pc - pex) / p0,
                         spread_over_p0=pc_spread / p0,
                         p_fe_quarter_over_p0=pq / p0, p_exact_quarter_over_p0=pexq / p0,
                         err_quarter_over_p0=abs(pq - pexq) / p0))
    return dict(stab=stab if isinstance(stab, str) else list(stab), t_ramp=t_ramp, nx=nx, ny=ny,
                dt=dt, p0=p0, constants=mc, rows=rows)


# ---------------------------------------------------------------------------------------- T4
def fe_hydro(nx=1, stab="off", t_end=None, dt=1.0, nely=None):
    """Self-weight saturated column (Terzaghi's soil), drained top p = 0, gravity on from t = 0
    (LadrunoUP -body is an ACCELERATION, x rho_mix in the element: quirk), run to steady state.
    Errors: p at nodes and particle points vs rho_f g (H - Y) over rho_f g H; sigma'_yy at the
    Gauss points and at the GIMP ppc-2 points (element-constant on this column) over rho' g H."""
    ops = fe_engine()
    c = HYDRO
    L, ny = c["L"], (c["nely"] if nely is None else nely)
    h = L / ny
    k = terzaghi_constants(c["E"], c["nu"], c["Kf"], c["poro"], c["kbar"], L, 1.0)
    t_end = 5.0 * L * L / k["cv"] if t_end is None else t_end
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    nid = {}
    n = 1
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.node(n, i * h, j * h)
            nid[(i, j)] = n
            n += 1
    ops.nDMaterial("ElasticIsotropic", 1, c["E"], c["nu"], c["rho"])
    e = 1
    for j in range(ny):
        for i in range(nx):
            ops.element("LadrunoUP", e, nid[(i, j)], nid[(i + 1, j)], nid[(i + 1, j + 1)],
                        nid[(i, j + 1)], 1, "-Kf", c["Kf"], "-poro", c["poro"],
                        "-rhoF", c["rhoF"], "-perm", c["kbar"], c["kbar"],
                        "-body", 0.0, -c["g"], "-dynSeepage", "off", *_stab_args(stab))
            e += 1
    for j in range(ny + 1):
        for i in range(nx + 1):
            fx = 1 if i in (0, nx) else 0
            fy = 1 if j == 0 else 0
            fp = 1 if j == ny else 0
            if fx or fy or fp:
                ops.fix(nid[(i, j)], fx, fy, fp)
    _transient(0.6, 0.3025)
    nst = int(round(t_end / dt))
    assert ops.analyze(nst, dt) == 0
    gw = c["rhoF"] * c["g"] * L
    gp = (c["rho"] - c["rhoF"]) * c["g"] * L
    perr = [abs(ops.nodeDisp(nid[(0, j)], 3) - hydrostatic(j * h, L, c["rho"], c["rhoF"], c["g"])[0]) / gw
            for j in range(ny + 1)]
    serr_g, serr_p = [], []
    for j in range(ny):
        st = np.array(ops.eleResponse(j * nx + 1, "stresses")).reshape(4, -1)
        # GP order of LadrunoUP Q4: assume the element's own GP coordinates; use the mean and the
        # element-constant check (this column's sigma'_yy is constant per element)
        syy = st[:, 1]
        for f in (0.5 - 0.5 / math.sqrt(3), 0.5 + 0.5 / math.sqrt(3)):
            serr_g.append(abs(np.mean(syy) - hydrostatic((j + f) * h, L, c["rho"], c["rhoF"], c["g"])[1]) / gp)
        for f in (0.25, 0.75):
            serr_p.append(abs(np.mean(syy) - hydrostatic((j + f) * h, L, c["rho"], c["rhoF"], c["g"])[1]) / gp)
        spread = float(np.ptp(syy))
    return dict(nx=nx, nely=ny, stab=stab if isinstance(stab, str) else list(stab), t_end=t_end,
                dt=dt, p_err_over_gwH=max(perr), sig_err_gauss_over_rhopgH=max(serr_g),
                sig_err_gimp_over_rhopgH=max(serr_p), last_elem_gp_spread=spread,
                note="sigma' is element-constant on this column; the error at a point is "
                     "rho' g |y - y_mid|, i.e. h/(2 sqrt3) at Gauss points, h/4 at GIMP points")


# ---------------------------------------------------------------------------------------- T5
def fe_b4(stab):
    """B4 McGann footing (tests/test_ladruno_up_element_equiv.py:332-411, re-run here): CB_lat,
    CB_lap at t = 1 s over the interior nodes."""
    ops = fe_engine()
    b = B4
    N = b["Ncell"]
    h = b["L"] / N
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    nid = {}
    n = 1
    for j in range(N + 1):
        for i in range(N + 1):
            ops.node(n, i * h, j * h)
            nid[(i, j)] = n
            n += 1
    ops.nDMaterial("ElasticIsotropic", 1, b["E"], b["nu"], b["rho"])
    e = 1
    for j in range(N):
        for i in range(N):
            ops.element("LadrunoUP", e, nid[(i, j)], nid[(i + 1, j)], nid[(i + 1, j + 1)],
                        nid[(i, j + 1)], 1, "-thick", 1.0, "-Kf", b["Kf"], "-poro", b["poro"],
                        "-rhoF", b["rhoF"], "-perm", b["perm"], b["perm"],
                        "-dynSeepage", "off", *_stab_args(stab))
            e += 1
    for j in range(N + 1):
        for i in range(N + 1):
            fx = 1 if (i == 0 or i == N) else 0
            fy = 1 if j == 0 else 0
            fp = 1 if j == N else 0
            if fx or fy or fp:
                ops.fix(nid[(i, j)], fx, fy, fp)
    ops.timeSeries("Path", 1, "-time", 0.0, 0.1, 1.0e6, "-values", 0.0, 1.0, 1.0)
    ops.pattern("Plain", 1, 1)
    ops.load(nid[(0, N)], 0.0, -b["q"] * b["strip"] / 2.0, 0.0)
    ops.load(nid[(1, N)], 0.0, -b["q"] * b["strip"] / 2.0, 0.0)
    ops.rayleigh(0.05, 0.02, 0.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormUnbalance", 1e-8, 50)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")
    for _ in range(int(round(b["T"] / b["dt"]))):
        assert ops.analyze(1, b["dt"]) == 0
    P = {(i, j): ops.nodeDisp(nid[(i, j)], 3) for i in range(N + 1) for j in range(N + 1)}
    inter = [(i, j) for i in range(1, N) for j in range(1, N)]
    p = np.array([P[q] for q in inter])
    chi = np.array([(-1.0) ** (i + j) for (i, j) in inter])
    cb_lat = abs(float(p @ chi)) / max(float(np.linalg.norm(p) * np.linalg.norm(chi)), 1e-30)
    lap = np.array([P[(i, j)] - 0.25 * (P[(i + 1, j)] + P[(i - 1, j)] + P[(i, j + 1)]
                                        + P[(i, j - 1)]) for (i, j) in inter])
    cb_lap = float(np.linalg.norm(lap)) / max(float(np.linalg.norm(p)), 1e-30)
    return dict(stab=stab if isinstance(stab, str) else list(stab), CB_lat=cb_lat, CB_lap=cb_lap,
                pmax=float(np.max(np.abs(p))))


# ==============================================================================================
def main(argv):
    leg = argv[1] if len(argv) > 1 else "closed"
    if leg == "closed":
        mc = mandel_constants(**MANDEL_SOIL)
        pk, tpk = mandel_peak(mc)
        out = dict(selfcheck=closed_selfcheck(), mandel_constants=mc, mandel_peak=pk,
                   mandel_peak_tv=tpk)
        print(json.dumps(out, indent=1))
        return
    if leg == "terzaghi":
        rec = _provenance("terzaghi_fe", "fe_terzaghi")
        rec["case"] = TERZ
        rec["runs"] = [fe_terzaghi(nx, st, first_small=fs)
                       for nx in (1, 2) for st in ("off", ("auto", 0.25)) for fs in (True, False)]
        print(_save("t2_terzaghi_fe", rec))
        return
    if leg == "mandel":
        mc = mandel_constants(**MANDEL_SOIL)
        rec = _provenance("mandel_fe", "fe_mandel")
        rec["case"] = dict(MANDEL, soil=MANDEL_SOIL)
        pk, tpk = mandel_peak(mc)
        rec["closed_form_peak"] = dict(p_over_p0=pk, tv=tpk)
        al = mandel_roots(mc["nu"], mc["nu_u"])
        tv = np.linspace(0.03, 0.3, 5401)
        pv = np.array([mandel_p(0.0, T * MANDEL["a"] ** 2 / mc["c"], MANDEL["a"], 1.0, mc,
                                t_ramp=MANDEL["t_ramp"], roots=al) for T in tv]) / mandel_p0(1.0, mc)
        rec["closed_form_peak_ramp"] = dict(p_over_p0=float(pv.max()), tv=float(tv[pv.argmax()]))
        rec["selfcheck"] = closed_selfcheck()
        rec["tvs"] = list(MANDEL_TVS)
        rec["undrained_lateral_transit_s"] = MANDEL["a"] / math.sqrt(
            (mc["Ku"] + 4 * mc["G"] / 3) / MANDEL["rho"])
        rec["closed_form_note"] = MANDEL_NOTE
        rec["peak_gate"] = 1.0 + 0.5 * (rec["closed_form_peak_ramp"]["p_over_p0"] - 1.0)
        # GATED (ruling 5820012e5): 5 s linear ramp from t = 0, uniform dt 0.02, both stab legs
        rec["gated"] = [fe_mandel(st, tvs=MANDEL_TVS, t_ramp=MANDEL["t_ramp"], first_small=False)
                        for st in ("off", ("auto", 0.25))]
        for run in rec["gated"]:
            run["band"] = [max(2.0 * max(r["err_over_p0"], r["err_quarter_over_p0"]), 0.02)
                           for r in run["rows"]]
        # REPORTED: the step load (unreachable by the FE: it rings), a 2 s ramp, dt 0.002
        rec["reported"] = [fe_mandel(st, tvs=MANDEL_TVS, t_ramp=tr)
                           for tr in (0.0, 2.0) for st in ("off", ("auto", 0.25))]
        rec["reported"] += [fe_mandel(st, tvs=MANDEL_TVS, t_ramp=MANDEL["t_ramp"], dt=0.002,
                                      first_small=False) for st in ("off", ("auto", 0.25))]
        print(_save("t3_mandel_fe", rec))
        return
    if leg == "t6fe":
        # T6 (ruling): the FE reference at route S's OWN dt, 5 s ramp, both cases, both stab legs
        dt = float(argv[2])
        rec = _provenance(f"t6_fe_dt{dt:g}", "fe_terzaghi+fe_mandel")
        rec["dt"] = dt
        rec["t_ramp"] = dict(terzaghi=TERZ_T6_RAMP, mandel=MANDEL["t_ramp"])
        rec["closed_form_note"] = MANDEL_NOTE
        rec["terzaghi"] = [fe_terzaghi(1, st, t_ramp=TERZ_T6_RAMP, first_small=False, dt=dt,
                                       tvs=TV_T2) for st in ("off", ("auto", 0.25))]
        rec["mandel"] = [fe_mandel(st, tvs=MANDEL_TVS, t_ramp=MANDEL["t_ramp"], dt=dt,
                                   first_small=False) for st in ("off", ("auto", 0.25))]
        for run in rec["terzaghi"]:
            run["band_p"] = [max(2.0 * r["p_err_gimp_points"], 0.02) for r in run["rows"]]
            run["band_U"] = [max(2.0 * r["U_err"], 0.02) for r in run["rows"]]
        for run in rec["mandel"]:
            run["band"] = [max(2.0 * max(r["err_over_p0"], r["err_quarter_over_p0"]), 0.02)
                           for r in run["rows"]]
        print(_save(f"t6_fe_dt{dt:g}", rec))
        return
    if leg == "ramp":
        rec = _provenance("terzaghi_ramp_fe", "fe_terzaghi")
        rec["case"] = dict(TERZ, t_ramp=TERZ_T6_RAMP)
        k = terzaghi_constants(TERZ["E"], TERZ["nu"], TERZ["Kf"], TERZ["poro"], TERZ["kbar"],
                               TERZ["L"], TERZ["q"])
        rec["undrained_transit_s"] = TERZ["L"] / math.sqrt((k["Eoed"] + k["Qbar"]) / TERZ["rho"])
        rec["runs"] = [fe_terzaghi(nx, st, t_ramp=TERZ_T6_RAMP, first_small=False, dt=dt)
                       for nx in (1, 2) for st in ("off", ("auto", 0.25)) for dt in (0.02, 0.002)]
        print(_save("t6_terzaghi_ramp_fe", rec))
        return
    if leg == "hydro":
        rec = _provenance("hydro_fe", "fe_hydro")
        rec["case"] = dict(HYDRO, nely=20, ruling="h = H/20 gated, band 0.03; H/10 reported")
        rec["closed_form"] = ("p = rho_f g (H - Y); sigma'_yy = -(rho_sat - rho_f) g (H - Y); "
                              "bars 0.01 rho_f g H (p) and 0.03 (rho_sat - rho_f) g H (sigma')")
        rec["gated"] = [fe_hydro(nx, st, nely=20) for nx in (1, 2) for st in ("off", ("auto", 0.25))]
        rec["reported"] = [fe_hydro(nx, st, nely=10) for nx in (1, 2)
                           for st in ("off", ("auto", 0.25))]
        print(_save("t4_hydro_fe", rec))
        return
    if leg == "dtsweep":
        # the FE error under the T6 ramp is Delta-t dependent: gamma 0.6 damps the slow
        # undrained modes per STEP, so a small (route-S-like) Delta-t leaves them ringing
        rec = _provenance("dtsweep_fe", "fe_terzaghi+fe_mandel")
        rec["t_ramp"] = dict(terzaghi=TERZ_T6_RAMP, mandel=MANDEL["t_ramp"])
        rec["terzaghi"] = [fe_terzaghi(1, "off", t_ramp=TERZ_T6_RAMP, first_small=False, dt=dt)
                           for dt in (0.05, 0.02, 0.01, 0.005, 0.002, 0.001)]
        rec["mandel"] = [fe_mandel("off", tvs=MANDEL_TVS, t_ramp=MANDEL["t_ramp"], dt=dt,
                                   first_small=False)
                         for dt in (0.05, 0.02, 0.01, 0.005, 0.002, 0.001)]
        print(_save("t6_dtsweep_fe", rec))
        return
    if leg == "plate":
        # Terzaghi under the interface's `mpm plate` mechanism: top-row u_y tied (equalDOF)
        rec = _provenance("terzaghi_plate_fe", "fe_terzaghi(plate=True)")
        rec["runs"] = [fe_terzaghi(nx, "off", plate=True) for nx in (1, 2)]
        rec["runs_nodal"] = [fe_terzaghi(nx, "off", plate=False) for nx in (1, 2)]
        print(_save("t2_terzaghi_plate_fe", rec))
        return
    if leg == "b4":
        rec = _provenance("b4_fe", "fe_b4")
        rec["case"] = B4
        rec["runs"] = [fe_b4(st) for st in ("off", ("auto", 0.25))]
        print(_save("t5_b4_fe", rec))
        return
    raise SystemExit(f"unknown leg {leg}")


if __name__ == "__main__":
    main(sys.argv)
