"""WP-232l (232l-V): closed forms and FE (LadrunoUP) references for two-phase particles -- never MPM.

Never reads the MPM C++ (SRC/domain/pattern/ladrunoMPM/, SRC/element/ladrunoMPMPoint/). Every
number comes from a closed form or from the fork's own LadrunoUP (ELE 33017) on the grid, time
step and Newmark pair the gate names (wp-briefs.md "232l", rows T2-T6). Units kN, t, m, s, kPa.
Sign: tension-positive total stress; p compression-positive (LadrunoUP.cpp:2235-2236).

Holds the closed forms, the FE legs that fix the T2-T6 bands (rulings 5820012e5), T1's LadrunoUP
references, and the numpy oracle for T0 i-iii, written from tests/mpm_twophase_interface.md §1,
§7a and the rulings 9064acca3 / 3d3699a50 alone. The oracle's particle element is self-checked
against LadrunoUP's own A and b (oracle_vs_fe_t1, <= 8e-16 on all 12 T1 legs).

    python3.12 tests/mpm_twophase_reference.py closed      # closed-form self-checks (no binary)
    python3.12 tests/mpm_twophase_reference.py terzaghi    # T2 FE: 1x10 and 2x10 Q4, step load
    python3.12 tests/mpm_twophase_reference.py mandel      # T3 FE: 10x4 half strip, tied plate
    python3.12 tests/mpm_twophase_reference.py hydro       # T4 FE: self-weight column, steady
    python3.12 tests/mpm_twophase_reference.py b4          # T5 FE: McGann B4 CB indices
    python3.12 tests/mpm_twophase_reference.py ramp        # T6 FE: T2/T3 cases with the T6 ramp
    python3.12 tests/mpm_twophase_reference.py t6fe <dt>   # T6 FE at the MPM run's own dt (bands)
    python3.12 tests/mpm_twophase_reference.py dtsweep     # T6 FE error against dt
    python3.12 tests/mpm_twophase_reference.py plate       # T2 plate (equalDOF) vs nodal loads
    python3.12 tests/mpm_twophase_reference.py t1          # T1 LadrunoUP A, b, first step
    python3.12 tests/mpm_twophase_reference.py t0          # T0 i-iii oracle (numpy only)

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



# ---------------------------------------------------------------------------------- R: TIMs
def fe_first_settlement(stab, dt=0.02, nsteps=10, rho=None, Kf=None):
    """The R row's first-settlement diagnostic (TIMs' C0-up: 1.71 stabilized, 0.94 not) on T2's
    column (1 x 10, plate, step load from t = 0, Newmark 0.6 / 0.3025): the plate settlement after
    each of the first nsteps steps of dt, over the closed-form undrained settlement
    s0 = q L/(Eoed + Qbar) (compressible fluid, alpha 1, Ks inf). rho: the mixture density
    (default T2's 2.0; a tiny rho removes inertia, so the first step is the quasi-static u(0+))."""
    ops = fe_engine()
    c = TERZ
    rho = c["rho"] if rho is None else rho
    Kf = c["Kf"] if Kf is None else Kf
    k = terzaghi_constants(c["E"], c["nu"], Kf, c["poro"], c["kbar"], c["L"], c["q"])
    L, ny = c["L"], c["nely"]
    h = L / ny
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    nid = {}
    n = 1
    for j in range(ny + 1):
        for i in range(2):
            ops.node(n, i * h, j * h)
            nid[(i, j)] = n
            n += 1
    ops.nDMaterial("ElasticIsotropic", 1, c["E"], c["nu"], rho)
    for j in range(ny):
        ops.element("LadrunoUP", j + 1, nid[(0, j)], nid[(1, j)], nid[(1, j + 1)],
                    nid[(0, j + 1)], 1, "-Kf", Kf, "-poro", c["poro"], "-rhoF", c["rhoF"],
                    "-perm", c["kbar"], c["kbar"], "-dynSeepage", "off", *_stab_args(stab))
    for j in range(ny + 1):
        for i in range(2):
            ops.fix(nid[(i, j)], 1, 1 if j == 0 else 0, 1 if j == ny else 0)
    ops.equalDOF(nid[(0, ny)], nid[(1, ny)], 2)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(nid[(0, ny)], 0.0, -c["q"] * h, 0.0)
    _transient(c["gamma"], c["beta"])
    hist = []
    for s_ in range(nsteps):
        assert ops.analyze(1, dt) == 0
        u = -ops.nodeDisp(nid[(0, ny)], 2)
        hist.append({"t": (s_ + 1) * dt, "settle": u, "ratio": u / k["s0"],
                     "p_mid_over_p0": ops.nodeDisp(nid[(0, ny // 2)], 3) / k["p0"]})
    return {"stab": stab if isinstance(stab, str) else list(stab), "dt": dt, "rho": rho, "Kf": Kf,
            "alpha_st": (0.0 if stab == "off" else
                         stab_alpha(h, c["E"] / (3 * (1 - 2 * c["nu"])),
                                    c["E"] / (2 * (1 + c["nu"])), stab[1])),
            "u_undrained": k["s0"], "u_drained": k["sinf"],
            "ratio_first_step": hist[0]["ratio"], "history": hist}


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


# ---------------------------------------------------------------------------------------- T1
# FE identity (route E): linear basis, particles at the 2x2 (2x2x2) Gauss points, V_p = Gauss
# weight x detJ = (h/2)^ndm, small strain. The common state is interface §7a Q5's re-seed: the
# particles carry sigma' = 0, F = I and the values of t1_particle_values; the FE nodal state is
# their projection, P5a for p, pd, pdd (u = 0, v = v^, a = a^, p = p^, pd = pd^, pdd = pdd^).
T1 = {"E": 1.0e4, "nu": 0.2, "rho": 2.0, "Kf": 4444.0, "poro": 0.4, "rhoF": 1.0, "kbar": 1.0e-2,
      "g": 9.81, "dt": 0.01, "gamma": 0.6, "beta": 0.3025, "h": 1.0, "renorm_cut": 0.05}
T1_MESHES = {"q4_1": (2, 1), "q4_2x2": (2, 2), "h8_1": (3, 1)}     # name -> (ndm, cells per axis)


def t1_dp_args():
    """DruckerPrager (UW) with the elastic leg's K, G; a high cone (sigma_y 100 kPa) so the
    first step stays elastic; non-associated (rho 0.2, rho_bar 0.1); hardening off."""
    c = T1
    K = c["E"] / (3 * (1 - 2 * c["nu"]))
    G = c["E"] / (2 * (1 + c["nu"]))
    return [K, G, 100.0, 0.2, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, c["rho"]]


def t1_particle_values(x):
    """The particles' prescribed state (z = 0 in 2D):
        v_i   = 0.10 sin(1.1 x + 0.7 y + 0.5 z + 0.3 i + 0.2)
        a_i   = 0.50 cos(0.9 x - 0.6 y + 0.4 z + 0.5 i + 0.1)
        p     = 5.0 + 2.0 sin(0.8 x + 1.2 y + 0.3 z)
        pd    = 0.5 cos(1.3 x - 0.4 y + 0.2 z)
        pdd   = 0.2 sin(0.6 x + 0.9 y - 0.7 z + 0.4)          (i = 0-based component)"""
    X = list(x) + [0.0] * (3 - len(x))
    xx, yy, zz = X
    nd = len(x)
    v = [0.10 * math.sin(1.1 * xx + 0.7 * yy + 0.5 * zz + 0.3 * i + 0.2) for i in range(nd)]
    a = [0.50 * math.cos(0.9 * xx - 0.6 * yy + 0.4 * zz + 0.5 * i + 0.1) for i in range(nd)]
    p = 5.0 + 2.0 * math.sin(0.8 * xx + 1.2 * yy + 0.3 * zz)
    pd = 0.5 * math.cos(1.3 * xx - 0.4 * yy + 0.2 * zz)
    pdd = 0.2 * math.sin(0.6 * xx + 0.9 * yy - 0.7 * zz + 0.4)
    return {"v": v, "a": a, "p": p, "pd": pd, "pdd": pdd}


def t1_particles(ndm, n, h=1.0, rho=2.0):
    """Particles at the 2^ndm Gauss points of each cell of an n^ndm grid on [0, n h]^ndm; order:
    cells x fastest, then y, z; inside a cell the Gauss points x fastest (-, +)."""
    g = 0.5 / math.sqrt(3.0)
    V = (h / 2.0) ** ndm
    out = []
    for cell in np.ndindex(*([n] * ndm)):
        cell = cell[::-1]                       # x fastest
        for gp in np.ndindex(*([2] * ndm)):
            gp = gp[::-1]
            x = [(cell[d] + 0.5 + (g if gp[d] else -g)) * h for d in range(ndm)]
            out.append(dict(x=x, V=V, m=rho * V, **t1_particle_values(x)))
    return out


def linear_support(x, n, h=1.0):
    """Linear hat (mpm_interface.md §4.2) on the real grid [0, n h]^ndm: per axis the clamped
    cell's two nodes; tensor product, x fastest. Returns [(idx tuple, S, dS)]."""
    ndm = len(x)
    ax = []
    for d in range(ndm):
        c = min(int(math.floor(x[d] / h)), n - 1)
        xi = (x[d] - c * h) / h
        ax.append([(c, 1.0 - xi, -1.0 / h), (c + 1, xi, 1.0 / h)])
    out = []
    for combo in np.ndindex(*([2] * ndm)):
        combo = combo[::-1]
        idx = tuple(ax[d][combo[d]][0] for d in range(ndm))
        S = 1.0
        for d in range(ndm):
            S *= ax[d][combo[d]][1]
        dS = []
        for d in range(ndm):
            v = ax[d][combo[d]][2]
            for e in range(ndm):
                if e != d:
                    v *= ax[e][combo[e]][1]
            dS.append(v)
        out.append((idx, S, np.array(dS)))
    return out


def renormalize(supports, masses, cut):
    """Interface §7a Q3 steps 2-4: nodal mass over every support node, keep m_a > cut max m_p,
    rescale S' = S/sigma, dS' = (dS - S' gamma)/sigma per particle over its kept nodes."""
    mnode = {}
    for sup, m in zip(supports, masses):
        for idx, S, _ in sup:
            mnode[idx] = mnode.get(idx, 0.0) + S * m
    thr = cut * max(masses)
    out = []
    for sup in supports:
        kept = [(i, S, dS) for (i, S, dS) in sup if mnode[i] > thr]
        sig = sum(S for _, S, _ in kept)
        gam = sum(dS for _, _, dS in kept)
        out.append([(i, S / sig, (dS - (S / sig) * gam) / sig) for (i, S, dS) in kept])
    return out, mnode


def p5_projection(parts, supports, two_phase=None, consistent=True):
    """Interface §7a Q5 with §7b P5a: m_a = sum S m; v^, a^ mass-weighted; p^, pd^, pdd^ by the
    consistent L2 projection over the two-phase particles (consistent=False: P5's lumped,
    volume-weighted form, kept for the record). Returns (nodal dict, M's smallest eigenvalue)."""
    acc = {}
    for k, (P, sup) in enumerate(zip(parts, supports)):
        tp = True if two_phase is None else two_phase[k]
        for idx, S, _ in sup:
            r = acc.setdefault(idx, {"m": 0.0, "mv": 0.0, "ma": 0.0, "V": 0.0, "Vp": 0.0,
                                     "Vpd": 0.0, "Vpdd": 0.0})
            r["m"] += S * P["m"]
            r["mv"] = r["mv"] + S * P["m"] * np.array(P["v"])
            r["ma"] = r["ma"] + S * P["m"] * np.array(P["a"])
            if tp:
                r["V"] += S * P["V"]
                r["Vp"] += S * P["V"] * P["p"]
                r["Vpd"] += S * P["V"] * P["pd"]
                r["Vpdd"] += S * P["V"] * P["pdd"]
    out = {}
    for idx, r in acc.items():
        V = r["V"]
        out[idx] = {"m": r["m"], "v": (r["mv"] / r["m"]).tolist(),
                    "a": (r["ma"] / r["m"]).tolist(),
                    "p": r["Vp"] / V if V > 0 else 0.0, "pd": r["Vpd"] / V if V > 0 else 0.0,
                    "pdd": r["Vpdd"] / V if V > 0 else 0.0}
    if not consistent:
        return out, None
    # P5a (interface §7b): p, pd, pdd by the consistent L2 projection over the p-active nodes
    # (every node a two-phase particle touches, here), M x^ = P^T V x_part, solved densely
    tp = [True] * len(parts) if two_phase is None else list(two_phase)
    pn = sorted({i for k, sup in enumerate(supports) if tp[k] for i, _, _ in sup},
                key=lambda t: t[::-1])
    X, eigmin = consistent_l2(
        [[(i, S) for i, S, _ in sup] for k, sup in enumerate(supports) if tp[k]],
        [P["V"] for k, P in enumerate(parts) if tp[k]],
        [[P["p"], P["pd"], P["pdd"]] for k, P in enumerate(parts) if tp[k]], pn)
    for a, idx in enumerate(pn):
        out[idx]["p"], out[idx]["pd"], out[idx]["pdd"] = (float(v) for v in X[a])
    return out, eigmin


def consistent_l2(weights, V, values, nodes):
    """P5a: solve M x^ = P^T V x with M_ab = sum_p V_p S_pa S_pb over `nodes` (dense; the oracle
    rule of the P5a freeze). weights[k] = [(node, S^p)], values[k] = the particle's row of fields.
    Returns (x^ [len(nodes), nfields], the smallest eigenvalue of M)."""
    ix = {n: a for a, n in enumerate(nodes)}
    nn = len(nodes)
    vals = np.atleast_2d(np.asarray(values, float))
    if vals.shape[0] != len(weights):
        vals = vals.T
    M = np.zeros((nn, nn))
    R = np.zeros((nn, vals.shape[1]))
    for w, Vp, x in zip(weights, V, vals):
        for i, Si in w:
            R[ix[i]] += Vp * Si * x
            for j, Sj in w:
                M[ix[i], ix[j]] += Vp * Si * Sj
    return np.linalg.solve(M, R), float(np.linalg.eigvalsh(M).min())


def stab_alpha(h, K, G, a0):
    """LadrunoUPKernel.h stabAlphaAuto: a0 h^2/(K + 4G/3)."""
    return a0 * h * h / (K + 4.0 * G / 3.0)


def fe_t1(mesh, mat, stab):
    """LadrunoUP on the T1 mesh at the common state: the Newmark first-iteration tangent A
    (printA), right side b (printB) and the first step (Linear algorithm: one exact solve).
    No boundary conditions (the effective tangent is regular: c3 M on u rows, H + c2 S* on p).
    Gravity g on the mixture (-body) and as the fluid drive (the -fluidBody default)."""
    ops = fe_engine()
    c = T1
    ndm, n = T1_MESHES[mesh]
    h = c["h"]
    ndf = ndm + 1
    parts = t1_particles(ndm, n, h, c["rho"])
    sups = [linear_support(P["x"], n, h) for P in parts]
    sups, _ = renormalize(sups, [P["m"] for P in parts], c["renorm_cut"])
    proj, eigmin = p5_projection(parts, sups)
    ops.wipe()
    ops.model("basic", "-ndm", ndm, "-ndf", ndf)
    tag = {}
    for k, idx in enumerate(sorted(proj, key=lambda t: t[::-1]), start=1):
        ops.node(k, *[float(i * h) for i in idx])
        tag[idx] = k
    if mat == "elastic":
        ops.nDMaterial("ElasticIsotropic", 1, c["E"], c["nu"], c["rho"])
    else:
        ops.nDMaterial("DruckerPrager", 1, *t1_dp_args())
    grav = [0.0] * ndm
    grav[1] = -c["g"]
    if ndm == 2:
        corners = [(0, 0), (1, 0), (1, 1), (0, 1)]
    else:
        corners = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                   (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    e = 1
    for cell in np.ndindex(*([n] * ndm)):
        cell = cell[::-1]
        conn = [tag[tuple(cell[d] + cc[d] for d in range(ndm))] for cc in corners]
        ops.element("LadrunoUP", e, *conn, 1, "-Kf", c["Kf"], "-poro", c["poro"],
                    "-rhoF", c["rhoF"], "-perm", *([c["kbar"]] * ndm), "-body", *grav,
                    "-lumped", "-geom", "linear", "-dynSeepage", "off", *_stab_args(stab))
        e += 1
    for idx, k in tag.items():
        st = proj[idx]
        for d in range(ndm):
            ops.setNodeDisp(k, d + 1, 0.0, "-commit")
            ops.setNodeVel(k, d + 1, st["v"][d], "-commit")
            ops.setNodeAccel(k, d + 1, st["a"][d], "-commit")
        ops.setNodeDisp(k, ndf, st["p"], "-commit")
        ops.setNodeVel(k, ndf, st["pd"], "-commit")
        ops.setNodeAccel(k, ndf, st["pdd"], "-commit")
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.algorithm("Linear")
    ops.integrator("Newmark", c["gamma"], c["beta"])
    ops.analysis("Transient")
    U0 = {k: [ops.nodeDisp(k, d + 1) for d in range(ndf)] for k in tag.values()}
    assert ops.analyze(1, c["dt"]) == 0
    A = np.array(ops.printA("-ret"))
    N = int(round(len(A) ** 0.5))
    A = A.reshape(N, N).T                           # column-major -> row-major
    # printB after a Linear step returns the post-solve SOE right side (~1e-14), not b. The step
    # IS one exact LU solve of A dU = b (FullGeneral), so b = A dU, dU = U1 - U0 (Newmark's
    # displacement predictor leaves U at U0 on every DOF, p included).
    b_after = np.array(ops.printB("-ret"))
    labels = [None] * N
    dU = np.zeros(N)
    nodes_out = []
    for idx, k in sorted(tag.items(), key=lambda kv: kv[1]):
        eqs = ops.nodeDOFs(k)
        for d in range(ndf):
            labels[eqs[d]] = [k] + [float(i * h) for i in idx] + [d]
            dU[eqs[d]] = ops.nodeDisp(k, d + 1) - U0[k][d]
        st = proj[idx]
        nodes_out.append({"tag": k, "x": [float(i * h) for i in idx], "m": st["m"],
                          "v0": st["v"], "a0": st["a"], "p0": st["p"], "pd0": st["pd"],
                          "pdd0": st["pdd"],
                          "u1": [ops.nodeDisp(k, d + 1) for d in range(ndm)],
                          "p1": ops.nodeDisp(k, ndf),
                          "v1": [ops.nodeVel(k, d + 1) for d in range(ndm)],
                          "pd1": ops.nodeVel(k, ndf),
                          "a1": [ops.nodeAccel(k, d + 1) for d in range(ndm)],
                          "pdd1": ops.nodeAccel(k, ndf)})
    b = A @ dU
    K0 = c["E"] / (3 * (1 - 2 * c["nu"]))
    G0 = c["E"] / (2 * (1 + c["nu"]))
    if stab == "off":
        alpha_st = 0.0
    elif isinstance(stab, tuple):
        alpha_st = stab_alpha(h, K0, G0, stab[1])
    else:
        alpha_st = float(stab)
    return {"mesh": mesh, "ndm": ndm, "cells": n, "material": mat,
            "stab": stab if isinstance(stab, str) else list(stab), "alpha_st": alpha_st,
            "projection": "P5a consistent L2 (interface §7b)", "M_eigmin": eigmin,
            "particles": parts, "nodes": nodes_out, "dofs": labels,
            "A": A.tolist(), "b": b.tolist(), "dU": dU.tolist(),
            "printB_after_step_norm": float(np.linalg.norm(b_after)),
            "A_fro": float(np.linalg.norm(A)), "b_norm": float(np.linalg.norm(b))}


# ==============================================================================================
# THE NUMPY ORACLE (T0 i-iii), written from tests/mpm_twophase_interface.md §1 and §7a, and the
# frozen interfaces it cites (mpm_interface.md §2, §4, §5, §6, §7; mpm_implicit_interface.md
# D1.5-D1.9, §3). Never from the MPM C++. Units and signs as above; Voigt plane strain
# [xx, yy, xy], 3D [xx, yy, zz, xy, yz, zx], engineering shear.
# ==============================================================================================
def gimp1(r, h, lp):
    """uGIMP 1D factor (mpm_interface.md §4.1; interface §7a Q3 step 1)."""
    a = abs(r)
    s = 1.0 if r > 0 else (-1.0 if r < 0 else 0.0)
    if a < lp:
        return 1.0 - (r * r + lp * lp) / (2.0 * h * lp), -r / (h * lp)
    if a < h - lp:
        return 1.0 - a / h, -s / h
    if a < h + lp:
        return (h + lp - a) ** 2 / (4.0 * h * lp), -s * (h + lp - a) / (2.0 * h * lp)
    return 0.0, 0.0


class OGrid:
    """Uniform grid with one ghost layer (mpm_interface.md §2). Node index tuples (i, j[, k]),
    -1 <= i <= n_x + 1; the C++ node id is reported for convenience."""

    def __init__(self, x0, h, n):
        self.x0 = [float(v) for v in x0]
        self.h = float(h)
        self.n = list(n)
        self.ndm = len(n)

    def node_id(self, idx):
        if self.ndm == 2:
            return (idx[0] + 1) + (self.n[0] + 3) * (idx[1] + 1)
        return (idx[0] + 1) + (self.n[0] + 3) * ((idx[1] + 1) + (self.n[1] + 3) * (idx[2] + 1))

    def coord(self, idx):
        return [self.x0[d] + idx[d] * self.h for d in range(self.ndm)]

    def face_nodes(self, face):
        """A face's node set: its face layer and its ghost layer (mpm_interface.md §2)."""
        ax = "xyz".index(face[0])
        lay = (-1, 0) if face[1:] == "min" else (self.n[ax], self.n[ax] + 1)
        return ax, lay

    def on_face(self, idx, face):
        ax, lay = self.face_nodes(face)
        return idx[ax] in lay

    def support(self, x, lp, basis="gimp"):
        """[(idx, S, dS)] in the C++ order: per-axis scan c-1..c+2 ascending, entries with
        S = dS = 0 dropped, tensor product x fastest (interface §7a Q1)."""
        h = self.h
        ax = []
        for d in range(self.ndm):
            c = int(math.floor((x[d] - self.x0[d]) / h))
            c = min(max(c, 0), self.n[d] - 1)
            lst = []
            if basis == "gimp":
                for i in range(c - 1, c + 3):
                    if i < -1 or i > self.n[d] + 1:
                        continue
                    S, dS = gimp1(x[d] - (self.x0[d] + i * h), h, lp[d])
                    if S == 0.0 and dS == 0.0:
                        continue
                    lst.append((i, S, dS))
            else:
                xi = (x[d] - (self.x0[d] + c * h)) / h
                lst = [(c, 1.0 - xi, -1.0 / h), (c + 1, xi, 1.0 / h)]
            ax.append(lst)
        out = []
        for combo in np.ndindex(*[len(a) for a in ax][::-1]):
            combo = combo[::-1]
            idx = tuple(ax[d][combo[d]][0] for d in range(self.ndm))
            S = 1.0
            for d in range(self.ndm):
                S *= ax[d][combo[d]][1]
            dS = np.empty(self.ndm)
            for d in range(self.ndm):
                v = ax[d][combo[d]][2]
                for e in range(self.ndm):
                    if e != d:
                        v *= ax[e][combo[e]][1]
                dS[d] = v
            out.append((idx, S, dS))
        return out


def elastic_D(E, nu, ndm):
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    G = E / (2 * (1 + nu))
    if ndm == 2:
        return np.array([[lam + 2 * G, lam, 0.0], [lam, lam + 2 * G, 0.0], [0.0, 0.0, G]])
    D = np.zeros((6, 6))
    D[:3, :3] = lam
    for i in range(3):
        D[i, i] = lam + 2 * G
        D[3 + i, 3 + i] = G
    return D


def bmat(dS, ndm):
    """Strain-displacement rows of one node (Voigt, engineering shear)."""
    if ndm == 2:
        x, y = dS
        return np.array([[x, 0.0], [0.0, y], [y, x]])
    x, y, z = dS
    return np.array([[x, 0, 0], [0, y, 0], [0, 0, z], [y, x, 0], [0, z, y], [z, 0, x]], float)


def voigt_tensor(s, ndm):
    if ndm == 2:
        return np.array([[s[0], s[2]], [s[2], s[1]]])
    return np.array([[s[0], s[3], s[5]], [s[3], s[1], s[4]], [s[5], s[4], s[2]]])


def stab_value(fl, D0, h, ndm):
    """alpha_st (interface §1): off 0, a raw value, or a0 h^2/(K + 4G/3) with G = D0(2,2)
    [3D: D0(3,3)], lam = D0(0,1), K = lam + 2G/3."""
    st = fl.get("stab", ("auto", 0.25))
    if st == "off":
        return 0.0
    if isinstance(st, (tuple, list)) and st[0] == "auto":
        G = D0[2, 2] if ndm == 2 else D0[3, 3]
        lam = D0[0, 1]
        return stab_alpha(h, lam + 2.0 * G / 3.0, G, float(st[1]))
    return float(st)


def particle_blocks(sup, V, m, D, sig_n, fl, alpha_st, grav, lam, tr, rayleigh=(0.0, 0.0)):
    """One two-phase particle element (interface §1, §7a Q2) at the trial nodal values
    tr = {u: [nen][ndm], p: [nen], v, pd, a, pdd}. sup = [(idx, S, dS)] (S^p = S on GIMP and
    linear). Returns dict K, C, M, R, R_static (node-major [u..., p]), sig (trial sigma')."""
    nen = len(sup)
    ndm = len(sup[0][2])
    ndf = ndm + 1
    nd = nen * ndf
    S = np.array([s for _, s, _ in sup])
    dS = np.array([d for _, _, d in sup])
    al = fl.get("alpha", 1.0)
    Ks = fl.get("Ks", -1.0)
    ooQ = fl["poro"] / fl["Kf"] + ((al - fl["poro"]) / Ks if Ks > 0 else 0.0)
    kb = np.array(fl["perm"], float)
    u = np.array(tr["u"], float)
    deps = sum(bmat(dS[k], ndm) @ u[k] for k in range(nen))
    sig = np.array(sig_n, float) + D @ deps
    st = voigt_tensor(sig, ndm)
    Ku = np.zeros((nen * ndm, nen * ndm))
    for a in range(nen):
        Ba = bmat(dS[a], ndm)
        for b in range(nen):
            Ku[a * ndm:(a + 1) * ndm, b * ndm:(b + 1) * ndm] = V * Ba.T @ D @ bmat(dS[b], ndm)
    Q = np.zeros((nen * ndm, nen))
    for a in range(nen):
        for i in range(ndm):
            Q[a * ndm + i, :] = al * V * dS[a, i] * S
    H = V * (dS * kb) @ dS.T
    Sm = V * ooQ * np.outer(S, S)
    Ht = V * dS @ dS.T
    Sst = Sm + alpha_st * Ht
    Mu = np.diag(np.repeat(S * m, ndm))
    Cu = rayleigh[0] * Mu + rayleigh[1] * Ku
    uI = [k * ndf + i for k in range(nen) for i in range(ndm)]
    pI = [k * ndf + ndm for k in range(nen)]
    K = np.zeros((nd, nd))
    C = np.zeros((nd, nd))
    M = np.zeros((nd, nd))
    K[np.ix_(uI, uI)] = Ku
    K[np.ix_(uI, pI)] = -Q
    K[np.ix_(pI, pI)] = H
    C[np.ix_(uI, uI)] = Cu
    C[np.ix_(pI, uI)] = Q.T
    C[np.ix_(pI, pI)] = Sst
    M[np.ix_(uI, uI)] = Mu
    p = np.array(tr["p"], float)
    g = np.array(grav, float) * lam
    fprime = np.array([V * st @ dS[k] for k in range(nen)]).ravel()
    fbody = np.array([S[k] * m * g for k in range(nen)]).ravel()
    fseep = V * fl["rhoF"] * dS @ (kb * g)
    Rs = np.zeros(nd)
    Rs[uI] = fprime - Q @ p - fbody
    Rs[pI] = H @ p - fseep
    xv = np.zeros(nd)
    xa = np.zeros(nd)
    xv[uI] = np.array(tr["v"], float).ravel()
    xv[pI] = tr["pd"]
    xa[uI] = np.array(tr["a"], float).ravel()
    xa[pI] = tr["pdd"]
    R = Rs + M @ xa + C @ xv
    return {"K": K, "C": C, "M": M, "R": R, "R_static": Rs, "sig": sig, "deps": deps,
            "nodes": [i for i, _, _ in sup], "S": S, "dS": dS, "alpha_st": alpha_st}


def pattern_values(nen, ndm, s):
    """Interface §7a Q1: value_j = s_kind sin(1.7 j + 0.3), j the element DOF index."""
    ndf = ndm + 1
    out = {k: [] for k in ("u", "p", "v", "pd", "a", "pdd")}
    for k in range(nen):
        u, v, a = [], [], []
        for c in range(ndf):
            w = math.sin(1.7 * (ndf * k + c) + 0.3)
            if c < ndm:
                u.append(s["u"] * w)
                v.append(s["v"] * w)
                a.append(s["a"] * w)
            else:
                out["p"].append(s["p"] * w)
                out["pd"].append(s["pd"] * w)
                out["pdd"].append(s["pdd"] * w)
        out["u"].append(u)
        out["v"].append(v)
        out["a"].append(a)
    return out


# ---------------------------------------------------------------- particle sets and state
def lattice(grid, lo, hi, ppc, rho):
    """mpm_interface.md §5: sub-lattice x0 + ((a + 0.5) h)/n, a_x fastest; box closed."""
    ndm = grid.ndm
    axes = []
    for d in range(ndm):
        pts = [grid.x0[d] + ((a + 0.5) * grid.h) / ppc for a in range(ppc * grid.n[d])]
        axes.append([x for x in pts if lo[d] <= x <= hi[d]])
    X = []
    for combo in np.ndindex(*[len(a) for a in axes][::-1]):
        combo = combo[::-1]
        X.append([axes[d][combo[d]] for d in range(ndm)])
    V0 = (grid.h / ppc) ** ndm
    return np.array(X), V0, rho * V0, grid.h / (2.0 * ppc)


def t0_initial(X, case):
    """The T0 initial particle fields (formulas quoted in the cases file)."""
    N, ndm = X.shape
    st = {"x": X.copy(), "v": np.zeros((N, ndm)), "a": np.zeros((N, ndm)),
          "p": np.zeros(N), "pd": np.zeros(N), "pdd": np.zeros(N),
          "sig": np.zeros((N, 3 if ndm == 2 else 6))}
    for k, (x, y) in enumerate(X):
        st["p"][k] = 5.0 + 2.0 * x + 3.0 * (case["H"] - y)
        if case.get("modeA", False):
            st["v"][k] = [0.05 * math.sin(math.pi * y / 2), 0.03 * math.sin(math.pi * x / 2)]
            st["a"][k] = [0.1 * math.cos(x), -0.05 * math.sin(y)]
            st["pd"][k] = 0.5 * math.sin(x + y)
            st["pdd"][k] = 0.1 * math.cos(y)
    return st


T0A = {"x0": (0.0, 0.0), "h": 1.0, "cells": (2, 2), "box": ((0.0, 0.0), (2.0, 2.0)), "ppc": 2,
       "H": 2.0, "rho": 2.0, "E": 1.0e4, "nu": 0.2,
       "fluid": {"Kf": 4444.0, "poro": 0.4, "rhoF": 1.0, "perm": (1.0e-2, 1.0e-2),
                 "stab": ("auto", 0.25)},
       "gravity": (0.0, -9.81), "fix": [("ymin", (1,))], "pfix": [("ymax", 0.0)],
       "rayleigh": (0.2, 0.001), "newmark": (0.6, 0.3025), "dt": 0.01, "steps": 3,
       "flip": 1.0, "pblend": 1.0, "renorm_cut": 0.05, "modeA": True,
       "dump": {"particles": (0, 5, 15),
                "s": {"u": 1.0e-3, "p": 2.0, "v": 0.05, "pd": 0.3, "a": 0.4, "pdd": 0.1}}}

T0S = {"x0": (0.0, 0.0), "h": 0.5, "cells": (4, 4), "box": ((0.0, 0.0), (2.0, 2.0)), "ppc": 2,
       "H": 2.0, "rho": 2.0, "E": 1.0e4, "nu": 0.2,
       "fluid": {"Kf": 4444.0, "poro": 0.4, "rhoF": 1.0, "perm": (1.0e-2, 1.0e-2),
                 "stab": ("auto", 0.25)},
       "gravity": (0.0, -9.81), "fix": [("xmin", (0,)), ("xmax", (0,)), ("ymin", (1,))],
       "pfix": [("ymax", 0.0)], "dt": 1.0e-3, "steps": 10, "flip": 1.0, "pblend": 1.0,
       "masscut": 1.0e-10, "pcut": 0.05, "modeA": False}


def _fixed_comps(grid, idx, case):
    out = set()
    for face, comps in case["fix"]:
        if grid.on_face(idx, face):
            out.update(comps)
    return out


def _pfix(grid, idx, case):
    for face, p0 in case["pfix"]:
        if grid.on_face(idx, face):
            return p0
    return None


# --------------------------------------------------------------------------- T0 (i)
def t0_particle_dump(case=T0A):
    """T0 (i): `mpm debug particleBlocks` at the initial committed state (before step 1),
    for each listed particle, with the case's pattern scales. λ = 1 (series none)."""
    grid = OGrid(case["x0"], case["h"], case["cells"])
    X, V0, m, lpv = lattice(grid, *case["box"], case["ppc"], case["rho"])
    ndm = grid.ndm
    lp = [lpv] * ndm
    D = elastic_D(case["E"], case["nu"], ndm)
    a_st = stab_value(case["fluid"], D, grid.h, ndm)
    sups = [grid.support(x, lp) for x in X]
    sups, _ = renormalize(sups, [m] * len(X), case["renorm_cut"])
    sig0 = np.zeros(3 if ndm == 2 else 6)
    out = []
    for pid in case["dump"]["particles"]:
        sup = sups[pid]
        tr = pattern_values(len(sup), ndm, case["dump"]["s"])
        blk = particle_blocks(sup, V0, m, D, sig0, case["fluid"], a_st, case["gravity"], 1.0,
                              tr, case["rayleigh"])
        out.append({"particle": pid, "x": X[pid].tolist(), "V": V0, "m": m, "alpha_st": a_st,
                    "nodes": [grid.node_id(i) for i in blk["nodes"]],
                    "node_idx": [list(i) for i in blk["nodes"]],
                    "coords": [grid.coord(i) for i in blk["nodes"]],
                    "S": blk["S"].tolist(), "dS": blk["dS"].tolist(), "trial": tr,
                    "K": blk["K"].tolist(), "C": blk["C"].tolist(), "M": blk["M"].tolist(),
                    "R": blk["R"].tolist(), "R_static": blk["R_static"].tolist()})
    return out


# --------------------------------------------------------------------------- T0 (ii)
def modeA_step(st, grid, case, V0, m, lp, D, a_st, lam=1.0):
    """One mode-A Newmark step, algorithm Linear (one exact solve), `kinematics small`
    (interface §1 route E, §7a Q3/Q5; mpm_implicit_interface.md §3, D1.9). Returns the new
    particle state and the grid solution."""
    ndm = grid.ndm
    ndf = ndm + 1
    gam, bet = case["newmark"]
    dt = case["dt"]
    N = len(st["x"])
    sups = [grid.support(x, lp) for x in st["x"]]
    sups, _ = renormalize(sups, [m] * N, case["renorm_cut"])
    parts = [{"m": m, "V": V0, "v": st["v"][k], "a": st["a"][k], "p": st["p"][k],
              "pd": st["pd"][k], "pdd": st["pdd"][k]} for k in range(N)]
    proj, eigmin = p5_projection(parts, sups)
    nodes = sorted(proj, key=lambda t: t[::-1])
    eq = {idx: k * ndf for k, idx in enumerate(nodes)}
    neq = len(nodes) * ndf
    U0, V0n, A0 = np.zeros(neq), np.zeros(neq), np.zeros(neq)
    fixed = np.zeros(neq, bool)
    for idx in nodes:
        e0 = eq[idx]
        fc = _fixed_comps(grid, idx, case)
        for d in range(ndm):
            if d in fc:
                fixed[e0 + d] = True
            else:
                V0n[e0 + d] = proj[idx]["v"][d]
                A0[e0 + d] = proj[idx]["a"][d]
        p0 = _pfix(grid, idx, case)
        if p0 is not None:
            fixed[e0 + ndm] = True
            U0[e0 + ndm] = p0
        else:
            U0[e0 + ndm] = proj[idx]["p"]
            V0n[e0 + ndm] = proj[idx]["pd"]
            A0[e0 + ndm] = proj[idx]["pdd"]
    c2 = gam / (bet * dt)
    c3 = 1.0 / (bet * dt * dt)
    Up = U0.copy()
    Vp = (1 - gam / bet) * V0n + dt * (1 - gam / (2 * bet)) * A0
    Ap = -V0n / (bet * dt) + (1 - 1 / (2 * bet)) * A0
    Vp[fixed] = 0.0
    Ap[fixed] = 0.0
    Amat = np.zeros((neq, neq))
    b = np.zeros(neq)
    for k in range(N):
        sup = sups[k]
        dofs = [eq[i] + c for i, _, _ in sup for c in range(ndf)]
        tr = {"u": [[Up[eq[i] + d] for d in range(ndm)] for i, _, _ in sup],
              "p": [Up[eq[i] + ndm] for i, _, _ in sup],
              "v": [[Vp[eq[i] + d] for d in range(ndm)] for i, _, _ in sup],
              "pd": [Vp[eq[i] + ndm] for i, _, _ in sup],
              "a": [[Ap[eq[i] + d] for d in range(ndm)] for i, _, _ in sup],
              "pdd": [Ap[eq[i] + ndm] for i, _, _ in sup]}
        blk = particle_blocks(sup, V0, m, D, st["sig"][k], case["fluid"], a_st,
                              case["gravity"], lam, tr, case.get("rayleigh", (0.0, 0.0)))
        Amat[np.ix_(dofs, dofs)] += blk["K"] + c2 * blk["C"] + c3 * blk["M"]
        b[dofs] -= blk["R"]
    free = ~fixed
    dU = np.zeros(neq)
    dU[free] = np.linalg.solve(Amat[np.ix_(free, free)], b[free])
    U1 = Up + dU
    V1 = Vp + c2 * dU
    A1 = Ap + c3 * dU
    V1[fixed] = 0.0
    A1[fixed] = 0.0
    new = {k: (v.copy() if hasattr(v, "copy") else v) for k, v in st.items()}
    aF, aP = case["flip"], case["pblend"]
    for k in range(N):
        sup = sups[k]
        S = np.array([s for _, s, _ in sup])
        rows = [eq[i] for i, _, _ in sup]
        u = np.array([U1[r:r + ndm] for r in rows])
        vg = np.array([V1[r:r + ndm] for r in rows])
        ag = np.array([A1[r:r + ndm] for r in rows])
        vh = np.array([proj[i]["v"] for i, _, _ in sup])
        ah = np.array([proj[i]["a"] for i, _, _ in sup])
        pg, pdg, pddg = (np.array([U1[r + ndm] for r in rows]), np.array([V1[r + ndm] for r in rows]),
                         np.array([A1[r + ndm] for r in rows]))
        ph = np.array([proj[i]["p"] for i, _, _ in sup])
        pdh = np.array([proj[i]["pd"] for i, _, _ in sup])
        pddh = np.array([proj[i]["pdd"] for i, _, _ in sup])
        new["x"][k] = st["x"][k] + S @ u
        new["v"][k] = aF * (st["v"][k] + S @ (vg - vh)) + (1 - aF) * (S @ vg)
        new["a"][k] = aF * (st["a"][k] + S @ (ag - ah)) + (1 - aF) * (S @ ag)
        new["p"][k] = aP * (st["p"][k] + S @ (pg - ph)) + (1 - aP) * (S @ pg)
        new["pd"][k] = aP * (st["pd"][k] + S @ (pdg - pdh)) + (1 - aP) * (S @ pdg)
        new["pdd"][k] = aP * (st["pdd"][k] + S @ (pddg - pddh)) + (1 - aP) * (S @ pddg)
        deps = sum(bmat(sup[j][2], ndm) @ u[j] for j in range(len(sup)))
        new["sig"][k] = st["sig"][k] + D @ deps
    grid_out = {"nodes": [grid.node_id(i) for i in nodes], "node_idx": [list(i) for i in nodes],
                "U": U1.tolist(), "V": V1.tolist(), "A": A1.tolist(), "M_eigmin": eigmin,
                "phat": [proj[i]["p"] for i in nodes], "pdhat": [proj[i]["pd"] for i in nodes],
                "pddhat": [proj[i]["pdd"] for i in nodes]}
    return new, grid_out


def _state_out(st):
    return {k: np.asarray(v).tolist() for k, v in st.items()}


def t0_modeA(case=T0A):
    """T0 (ii): three mode-A Newmark steps from the T0 initial state; the particle fields after
    each step (and the grid solution, for diagnosis)."""
    grid = OGrid(case["x0"], case["h"], case["cells"])
    X, V0, m, lpv = lattice(grid, *case["box"], case["ppc"], case["rho"])
    lp = [lpv] * grid.ndm
    D = elastic_D(case["E"], case["nu"], grid.ndm)
    a_st = stab_value(case["fluid"], D, grid.h, grid.ndm)
    st = t0_initial(X, case)
    hist = [{"step": 0, "particles": _state_out(st)}]
    for n in range(case["steps"]):
        st, g = modeA_step(st, grid, case, V0, m, lp, D, a_st)
        hist.append({"step": n + 1, "t": (n + 1) * case["dt"], "particles": _state_out(st),
                     "grid": g})
    return {"V0": V0, "m": m, "lp": lpv, "alpha_st": a_st, "history": hist}


# --------------------------------------------------------------------------- T0 (iii)
def routeS_step(st, grid, case, V0, m, lp, D, a_st, lam=1.0):
    """One route-S step (interface §1 'Route S step', §7a Q4; mpm_interface.md §6 MUSL),
    `kinematics small`: Δε = sym(Δt Σ v** ⊗ ∇S) from the conditioned second projection, V = V0.
    p-active and the S^p renormalisation: ruling 3d3699a50 (see the pressure block)."""
    ndm = grid.ndm
    dt = case["dt"]
    fl = case["fluid"]
    N = len(st["x"])
    sups = [grid.support(x, lp) for x in st["x"]]
    al = fl.get("alpha", 1.0)
    ooQ = fl["poro"] / fl["Kf"]
    kb = np.array(fl["perm"], float)
    g = np.array(case["gravity"], float) * lam
    nodes = sorted({i for sup in sups for i, _, _ in sup}, key=lambda t: t[::-1])
    nid = {idx: k for k, idx in enumerate(nodes)}
    nn = len(nodes)
    # raw nodal masses (fix 3 and the u-activity test run on raw weights)
    mass = np.zeros(nn)
    Vw = np.zeros(nn)
    for k in range(N):
        for idx, S, _ in sups[k]:
            mass[nid[idx]] += S * m
            Vw[nid[idx]] += S * V0
    # §7c (1) corrected (9e20c87e0): a node is u-active iff m_a > max(massCutoff, pressureCutoff)
    # m_max (raw masses). Each particle's kept entries (on u-active nodes) are renormalised ALWAYS,
    # mode A's formula: S' = S/sigma, dS' = (dS - S' gamma)/sigma over the kept entries; entries on
    # other nodes are absent everywhere. The p weights are NOT touched by this pass (fix 3 below
    # runs on the raw weights).
    c_u = max(case["masscut"], case["pcut"])
    active = mass > c_u * m
    dropped_u = [a for a in range(nn) if mass[a] > 0.0 and not active[a]]
    usups = []
    for k in range(N):
        kept = [(i, S, dS) for i, S, dS in sups[k] if active[nid[i]]]
        sg = sum(S for _, S, _ in kept)
        gm = sum(dS for _, _, dS in kept)
        usups.append([(i, S / sg, (dS - (S / sg) * gm) / sg) for i, S, dS in kept])
    # P2G on the renormalised u weights (mass, momentum, forces)
    massu = np.zeros(nn)
    mom = np.zeros((nn, ndm))
    f = np.zeros((nn, ndm))
    for k in range(N):
        st_k = voigt_tensor(st["sig"][k], ndm)
        for idx, S, dS in usups[k]:
            a = nid[idx]
            massu[a] += S * m
            mom[a] += S * m * st["v"][k]
            f[a] += S * m * g - V0 * st_k @ dS
    vhat = np.zeros((nn, ndm))
    vhat[active] = mom[active] / massu[active, None]
    vstar = np.zeros((nn, ndm))
    vstar[active] = vhat[active] + dt * f[active] / massu[active, None]
    Minv = np.zeros((nn, ndm))
    Minv[active] = 1.0 / massu[active, None]
    for idx in nodes:
        for c in _fixed_comps(grid, idx, case):
            vstar[nid[idx], c] = 0.0
            Minv[nid[idx], c] = 0.0
    # p-activity (ruling 3d3699a50): m_a > c_m m_max (mode A's cutoff, default 5e-2; m_max =
    # max_p m_p, interface §7a Q3 step 3) AND S^p > 0 from a two-phase particle. Each particle's
    # S^p, grad S^p are renormalised over its p-active nodes (renormalize's formula); dropped
    # nodes are absent from the pressure system, the projection and G2P. Q keeps the
    # unrenormalised mode-F grad S on the u side.
    m_max = m
    pact = [a for a in range(nn) if Vw[a] > 0 and mass[a] > case["pcut"] * m_max]
    pidx = {a: k for k, a in enumerate(pact)}
    npn = len(pact)
    spr = []
    for k in range(N):
        kept = [(nid[i], S, dS) for i, S, dS in sups[k] if nid[i] in pidx]
        if not any(S > 0.0 and nid[i] not in pidx for i, S, _ in sups[k]):
            spr.append(kept)        # nothing with S^p > 0 dropped: no renormalisation (§7b fix 3)
            continue
        sg = sum(S for _, S, _ in kept)          # unsigned raw weights (GIMP: W = S)
        gm = sum(dS for _, _, dS in kept)
        spr.append([(a, S / sg, (dS - (S / sg) * gm) / sg) for a, S, dS in kept])
    # P5a: p^ is the consistent L2 projection over the p-active nodes (pinned ones included)
    Xh, eigmin = consistent_l2([[(a, S) for a, S, _ in spr[k]] for k in range(N)], [V0] * N,
                               [[st["p"][k]] for k in range(N)], pact)
    phat = np.zeros(nn)
    for k_, a in enumerate(pact):
        phat[a] = Xh[k_, 0]
    Sst = np.zeros((npn, npn))
    Hm = np.zeros((npn, npn))
    Qm = np.zeros((nn * ndm, npn))
    fse = np.zeros(npn)
    for k in range(N):
        sup = [(nid[i], S, dS) for i, S, dS in usups[k]]
        ps = [(pidx[a], S, dS) for a, S, dS in spr[k]]
        for b, Sb, dSb in ps:
            fse[b] += V0 * fl["rhoF"] * dSb @ (kb * g)
            for c_, Sc, dSc in ps:
                Sst[b, c_] += V0 * (ooQ * Sb * Sc + a_st * dSb @ dSc)
                Hm[b, c_] += V0 * dSb @ (kb * dSc)
            for a, Sa, dSa in sup:
                for i in range(ndm):
                    Qm[a * ndm + i, b] += al * V0 * dSa[i] * Sb
    Mi = Minv.ravel()
    Amat = Sst + dt * Hm + dt * dt * Qm.T @ (Mi[:, None] * Qm)
    ph = np.array([phat[a] for a in pact])
    rhs = Sst @ ph - dt * Qm.T @ vstar.ravel() + dt * fse
    drained = np.zeros(npn, bool)
    pval = np.zeros(npn)
    for k, a in enumerate(pact):
        p0 = _pfix(grid, nodes[a], case)
        if p0 is not None:
            drained[k] = True
            pval[k] = p0
    fr = ~drained
    rhs = rhs - Amat[:, drained] @ pval[drained]
    pval[fr] = np.linalg.solve(Amat[np.ix_(fr, fr)], rhs[fr])
    v1 = vstar + dt * (Mi * (Qm @ pval)).reshape(nn, ndm)
    acc = np.zeros((nn, ndm))
    acc[active] = (v1[active] - vhat[active]) / dt
    pg = np.zeros(nn)
    for k, a in enumerate(pact):
        pg[a] = pval[k]
    new = {k: (v.copy() if hasattr(v, "copy") else v) for k, v in st.items()}
    aF, aP = case["flip"], case["pblend"]
    for k in range(N):
        S = np.array([s for _, s, _ in usups[k]])
        rows = [nid[i] for i, _, _ in usups[k]]
        new["v"][k] = aF * (st["v"][k] + dt * S @ acc[rows]) + (1 - aF) * (S @ v1[rows])
        new["x"][k] = st["x"][k] + dt * S @ v1[rows]
        prow = [a for a, _, _ in spr[k]]
        Sp = np.array([s for _, s, _ in spr[k]])
        new["p"][k] = aP * (st["p"][k] + Sp @ (pg[prow] - phat[prow])) + (1 - aP) * (Sp @ pg[prow])
    # §7c (2): a two-phase particle's skeleton strain comes from the grid velocity v1 (USL), the
    # conditioned velocity after the correction (0 on fixed and dropped components)
    for k in range(N):
        L = sum(np.outer(v1[nid[i]], dS) for i, _, dS in usups[k])
        E = 0.5 * dt * (L + L.T)
        deps = np.array([E[0, 0], E[1, 1], 2 * E[0, 1]])
        new["sig"][k] = st["sig"][k] + D @ deps
    return new, {"nodes": [grid.node_id(nodes[a]) for a in pact],
                 "node_idx": [list(nodes[a]) for a in pact], "p": pval.tolist(),
                 "phat": ph.tolist(), "M_eigmin": eigmin,
                 "u_dropped": [list(nodes[a]) for a in dropped_u]}


def drained_pencil(st, grid, case, V0, m, lp, D):
    """§7c (3): dt_d = 2/sqrt(lambda_max) of M_L^-1 K_d at the step-start state -- lumped mass on
    the u-active free components (two-phase rule), K_d = sum_p V Bt D B on the raw supports; fixed
    and dropped components removed. Dense symmetric eigenproblem (the C++ runs 200 power
    iterations toward the same lambda_max)."""
    ndm = grid.ndm
    sups = [grid.support(x, lp) for x in st["x"]]
    nodes = sorted({i for sup in sups for i, _, _ in sup}, key=lambda t: t[::-1])
    nid = {idx: k for k, idx in enumerate(nodes)}
    nn = len(nodes)
    mass = np.zeros(nn)
    K = np.zeros((nn * ndm, nn * ndm))
    for sup in sups:
        for idx, S, _ in sup:
            mass[nid[idx]] += S * m
        rows = [nid[i] * ndm + c for i, _, _ in sup for c in range(ndm)]
        B = np.hstack([bmat(dS, ndm) for _, _, dS in sup])
        K[np.ix_(rows, rows)] += V0 * B.T @ D @ B
    act = mass > max(case["masscut"], case["pcut"]) * m
    keep = []
    for idx in nodes:
        a = nid[idx]
        fc = _fixed_comps(grid, idx, case)
        for c in range(ndm):
            if act[a] and c not in fc:
                keep.append(a * ndm + c)
    Kk = K[np.ix_(keep, keep)]
    sm = 1.0 / np.sqrt(np.array([mass[r // ndm] for r in keep]))
    lam = float(np.max(np.linalg.eigvalsh(sm[:, None] * Kk * sm[None, :])))
    return 2.0 / math.sqrt(lam)


def t0_routeS(case=T0S):
    grid = OGrid(case["x0"], case["h"], case["cells"])
    X, V0, m, lpv = lattice(grid, *case["box"], case["ppc"], case["rho"])
    lp = [lpv] * grid.ndm
    D = elastic_D(case["E"], case["nu"], grid.ndm)
    a_st = stab_value(case["fluid"], D, grid.h, grid.ndm)
    st = t0_initial(X, case)
    keep = ("x", "v", "p", "sig")
    hist = [{"step": 0, "particles": {k: np.asarray(st[k]).tolist() for k in keep}}]
    dtd = []
    for n in range(case["steps"]):
        dtd.append(drained_pencil(st, grid, case, V0, m, lp, D))
        st, g = routeS_step(st, grid, case, V0, m, lp, D, a_st)
        hist.append({"step": n + 1, "t": (n + 1) * case["dt"],
                     "particles": {k: np.asarray(st[k]).tolist() for k in keep}, "grid": g})
    return {"V0": V0, "m": m, "lp": lpv, "alpha_st": a_st, "history": hist,
            "dt_drained_pencil": dtd, "safety": 0.5,
            "rules": "interface §1, §7a Q4, §7b P5a/fix 3, §7c (u-activity, USL strain, pencil)"}


T13 = dict(T0S, h=1.0, cells=(6, 4), box=((1.0, 1.0), (3.0, 3.0)), gravity=(0.0, 0.0), fix=[],
           pfix=[], dt=5.0e-3, steps=200)


def t13_translation(case=T13, v0=(1.0, 0.0)):
    """T13's simplest leg (§7c (1) corrected), numpy expectation from the text: a free 2 x 2-cell
    gimp two-phase block in rigid translation, v = (1, 0), no loads, no conditions, p = 0, route S,
    dt 5e-3, 200 steps. Returns the max deviations of v (from v0), sigma' and p (from 0), and of
    x (from X + v0 t), over all steps and particles."""
    grid = OGrid(case["x0"], case["h"], case["cells"])
    X, V0, m, lpv = lattice(grid, *case["box"], case["ppc"], case["rho"])
    lp = [lpv] * grid.ndm
    D = elastic_D(case["E"], case["nu"], grid.ndm)
    a_st = stab_value(case["fluid"], D, grid.h, grid.ndm)
    N = len(X)
    st = {"x": X.copy(), "v": np.tile(np.array(v0, float), (N, 1)), "p": np.zeros(N),
          "sig": np.zeros((N, 3))}
    dv = dsig = dp = dx = 0.0
    for n in range(case["steps"]):
        st, _ = routeS_step(st, grid, case, V0, m, lp, D, a_st)
        t = (n + 1) * case["dt"]
        dv = max(dv, float(np.max(np.abs(st["v"] - np.array(v0)))))
        dsig = max(dsig, float(np.max(np.abs(st["sig"]))))
        dp = max(dp, float(np.max(np.abs(st["p"]))))
        dx = max(dx, float(np.max(np.abs(st["x"] - (X + np.array(v0) * t)))))
    return {"v_dev": dv, "sig_dev": dsig, "p_dev": dp, "x_dev": dx, "particles": N,
            "final_x_mean": st["x"].mean(axis=0).tolist()}


# --------------------------------------------------------------------------- self-checks
def oracle_vs_fe_t1(leg):
    """Self-check: the oracle's particle blocks at the T1 Gauss-point particles, assembled
    under Newmark's predictor from the P5 state, against LadrunoUP's A and b (t1_fe.json)."""
    c = T1
    ndm, n = leg["ndm"], leg["cells"]
    ndf = ndm + 1
    parts = leg["particles"]
    sups = [linear_support(P["x"], n, c["h"]) for P in parts]
    sups, _ = renormalize(sups, [P["m"] for P in parts], c["renorm_cut"])
    proj, eigmin = p5_projection(parts, sups)
    lab = {(tuple(int(round(v / c["h"])) for v in L[1:1 + ndm]), L[-1]): e
           for e, L in enumerate(leg["dofs"])}
    neq = len(leg["dofs"])
    U0, V0n, A0 = np.zeros(neq), np.zeros(neq), np.zeros(neq)
    for idx, s in proj.items():
        for d in range(ndm):
            V0n[lab[(idx, d)]] = s["v"][d]
            A0[lab[(idx, d)]] = s["a"][d]
        U0[lab[(idx, ndm)]] = s["p"]
        V0n[lab[(idx, ndm)]] = s["pd"]
        A0[lab[(idx, ndm)]] = s["pdd"]
    gam, bet, dt = c["gamma"], c["beta"], c["dt"]
    c2, c3 = gam / (bet * dt), 1.0 / (bet * dt * dt)
    Vp = (1 - gam / bet) * V0n + dt * (1 - gam / (2 * bet)) * A0
    Ap = -V0n / (bet * dt) + (1 - 1 / (2 * bet)) * A0
    D = elastic_D(c["E"], c["nu"], ndm)
    fl = {"Kf": c["Kf"], "poro": c["poro"], "rhoF": c["rhoF"], "perm": [c["kbar"]] * ndm}
    grav = [0.0] * ndm
    grav[1] = -c["g"]
    A = np.zeros((neq, neq))
    b = np.zeros(neq)
    for P, sup in zip(parts, sups):
        dofs = [lab[(i, cc)] for i, _, _ in sup for cc in range(ndf)]
        tr = {"u": [[0.0] * ndm for _ in sup], "p": [U0[lab[(i, ndm)]] for i, _, _ in sup],
              "v": [[Vp[lab[(i, d)]] for d in range(ndm)] for i, _, _ in sup],
              "pd": [Vp[lab[(i, ndm)]] for i, _, _ in sup],
              "a": [[Ap[lab[(i, d)]] for d in range(ndm)] for i, _, _ in sup],
              "pdd": [Ap[lab[(i, ndm)]] for i, _, _ in sup]}
        blk = particle_blocks(sup, P["V"], P["m"], D, np.zeros(3 if ndm == 2 else 6), fl,
                              leg["alpha_st"], grav, 1.0, tr)
        A[np.ix_(dofs, dofs)] += blk["K"] + c2 * blk["C"] + c3 * blk["M"]
        b[dofs] -= blk["R"]
    Afe, bfe = np.array(leg["A"]), np.array(leg["b"])
    return (float(np.linalg.norm(A - Afe) / np.linalg.norm(Afe)),
            float(np.linalg.norm(b - bfe) / np.linalg.norm(bfe)))


# --------------------------------------------------------------------- §7b numpy checks
def sealed_column_check(steps=3):
    """The P5a sealed-column invariant (§7b): motionless (every u fixed), sealed (no pfix),
    kbar = 0, pd = pdd = 0, gravity 0, sigma' = 0, and a particle p field the grid CANNOT
    represent (9 particles, 4 nodes). Mode-A steps must keep every particle p bit for bit,
    whatever x^ is. Returns (max |p_n - p_0| over the steps, the invisible residual norm)."""
    case = dict(T0A, cells=(1, 1), box=((0.0, 0.0), (1.0, 1.0)), ppc=3, gravity=(0.0, 0.0),
                fix=[("xmin", (0, 1)), ("xmax", (0, 1)), ("ymin", (0, 1)), ("ymax", (0, 1))],
                pfix=[], rayleigh=(0.0, 0.0),
                fluid=dict(T0A["fluid"], perm=(0.0, 0.0)))
    grid = OGrid(case["x0"], case["h"], case["cells"])
    X, V0, m, lpv = lattice(grid, *case["box"], case["ppc"], case["rho"])
    lp = [lpv] * 2
    D = elastic_D(case["E"], case["nu"], 2)
    a_st = stab_value(case["fluid"], D, grid.h, 2)
    N = len(X)
    st = {"x": X.copy(), "v": np.zeros((N, 2)), "a": np.zeros((N, 2)),
          "p": np.array([5.0 + 3.0 * math.sin(7.0 * x) * math.cos(5.0 * y) for x, y in X]),
          "pd": np.zeros(N), "pdd": np.zeros(N), "sig": np.zeros((N, 3))}
    p0 = st["p"].copy()
    worst = 0.0
    for _ in range(steps):
        st, g = modeA_step(st, grid, case, V0, m, lp, D, a_st)
        worst = max(worst, float(np.max(np.abs(st["p"] - p0))))
    sups, _ = renormalize([grid.support(x, lp) for x in X], [m] * N, case["renorm_cut"])
    interp = np.array([sum(S * g["U"][3 * g["node_idx"].index(list(i)) + 2] for i, S, _ in sup)
                       for sup in sups])
    return worst, float(np.linalg.norm(p0 - interp)), st


def bspline2_1d(r, h):
    """mpm_bspline_interface.md B1, quadratic."""
    xi = r / h
    a = abs(xi)
    s = 1.0 if xi > 0 else (-1.0 if xi < 0 else 0.0)
    if a < 0.5:
        return 0.75 - xi * xi, -2.0 * xi / h
    if a < 1.5:
        return 0.5 * (1.5 - a) ** 2, -s * (1.5 - a) / h
    return 0.0, 0.0


def fix2_hydrostatic_check(H=2.0, W=1.0, h=0.5, ppc=2, rhoF=1.0, g=9.81, lam=1.0,
                           affine=True):
    """§7b fix 2, the static algebra, on a bspline2 column x in [0, W], y in [0, H]: xmin/xmax
    impervious with a tangential drive (p images even), ymin impervious with the normal drive
    b_f = (0, -g) (the AFFINE even image c_g = c_m + rhoF lam b_f,y (Y_g - Y_m)), ymax drained
    p0 = 0 (odd image). With the nodal coefficients c_a = rhoF g (H - Y_a) on the real nodes:
      (a) sum_a S'_a c_a + pi_p = rhoF g (H - y_p) at every particle (values and gradients),
      (b) the p-row flux integrand grad(p_h) - rhoF b_f vanishes at every particle (zero flux:
          the hydrostatic field is an exact steady solution of the p rows),
      (c) the P5a projection of x_part - pi returns the c_a, and G2P (adding pi back) returns
          the particle p exactly.
    affine=False is the negative control (P9's plain even image). Returns the three max errors
    (over gamma_w H; (b) over gamma_w)."""
    nx, ny = int(round(W / h)), int(round(H / h))
    b = np.array([0.0, -g]) * lam
    phys = lambda y: rhoF * g * (H - y)          # noqa: E731
    pts = []
    for j in range(ppc * ny):
        for i in range(ppc * nx):
            pts.append(((i + 0.5) * h / ppc, (j + 0.5) * h / ppc))
    errs = [0.0, 0.0, 0.0]
    weights, values, nodes = [], [], set()
    for (x, y) in pts:
        ax = []
        for d, (xp, n) in enumerate(((x, nx), (y, ny))):
            c = min(int(math.floor(xp / h)), n - 1)
            lst = []
            for i in range(c - 1, c + 3):
                S, dS = bspline2_1d(xp - i * h, h)
                if S == 0.0 and dS == 0.0:
                    continue
                lst.append((i, S, dS))
            ax.append(lst)
        Sp = {}
        pi_v, pi_g = 0.0, np.zeros(2)
        for (ix, Sx, dSx) in ax[0]:
            for (iy, Sy, dSy) in ax[1]:
                S, dS = Sx * Sy, np.array([dSx * Sy, Sx * dSy])
                jx, sx = ix, 1.0
                if ix < 0:
                    jx = -ix                         # xmin image, even
                elif ix > nx:
                    jx = 2 * nx - ix                 # xmax image, even
                jy, sy, aff = iy, 1.0, 0.0
                if iy < 0:
                    jy = -iy                         # ymin image, affine even
                    aff = rhoF * b[1] * (iy * h - jy * h) if affine else 0.0
                elif iy > ny:
                    jy, sy = 2 * ny - iy, -1.0       # ymax image, odd (p0 = 0)
                key = (jx, jy)
                w = sx * sy
                if key in Sp:
                    Sp[key] = (Sp[key][0] + w * S, Sp[key][1] + w * dS)
                else:
                    Sp[key] = (w * S, w * dS)
                pi_v += S * aff * sx
                pi_g += dS * aff * sx
        ph = sum(S * phys(k[1] * h) for k, (S, _) in Sp.items()) + pi_v
        gr = sum(dS * phys(k[1] * h) for k, (_, dS) in Sp.items()) + pi_g
        errs[0] = max(errs[0], abs(ph - phys(y)) / (rhoF * g * H))
        errs[1] = max(errs[1], float(np.max(np.abs(gr - rhoF * b))) / (rhoF * g))
        weights.append([(k, S) for k, (S, _) in Sp.items()])
        values.append([phys(y) - pi_v])
        nodes.update(Sp)
    nodes = sorted(nodes, key=lambda t: t[::-1])
    X, _ = consistent_l2(weights, [(h / ppc) ** 2] * len(pts), values, nodes)
    cmap = {n: X[a, 0] for a, n in enumerate(nodes)}
    for n, cval in cmap.items():
        errs[2] = max(errs[2], abs(cval - phys(n[1] * h)) / (rhoF * g * H))
    return errs


# ==============================================================================================
# R ROWS, route E (mode A) -- reported, never gated. The MPM is driven only through the public
# `mpm` commands of tests/mpm_twophase_interface.md (§3-§7c) and the cases file; the MPM C++ is
# never read. LadrunoUP legs come from the functions above, in the same process (same build).
# ==============================================================================================
def _mpm_deck(ops, basis, cells, box, E=1.0e4, nu=0.2, rho=2.0, Kf=4444.0, n=0.4, rhoF=1.0,
              kbar=1.0e-4, stab="off", fix=(), pfix=(), plate=None, gravity=None, linear=True,
              blend=None):
    """One mode-A two-phase deck (ndm 2), Newmark 0.6 / 0.3025, UmfPack, Transformation."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)
    h = (box[1][1] - box[0][1]) / cells[1]
    ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", h, "-cells", cells[0], cells[1],
            "-basis", basis)
    ops.mpm("mode", "implicit")
    ops.mpm("body", 1, "box", box[0][0], box[0][1], box[1][0], box[1][1], "-ppc", 2, "-mat", 1,
            "-rho", rho)
    sa = ["-stab", "off"] if stab == "off" else ["-stab", "auto", float(stab[1])]
    ops.mpm("twoPhase", 1, "-Kf", Kf, "-poro", n, "-rhoF", rhoF, "-perm", kbar, kbar, *sa)
    if gravity:
        ops.mpm("gravity", 0.0, -gravity)
    for face, comps in fix:
        ops.mpm("fix", face, comps)
    for face, p0 in pfix:
        ops.mpm("pfix", face, p0)
    if plate:
        ops.mpm("plate", plate[0], plate[1])
    if blend is not None:
        ops.mpm("pressureBlend", blend)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    if linear:
        ops.algorithm("Linear")
    else:
        ops.test("NormDispIncr", 1e-10, 25)
        ops.algorithm("Newton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")


def _flat(v):
    return [float(x) for x in v] if isinstance(v, (list, tuple)) else [float(v)]


def _plate_node(ops):
    grid = set(int(t) for t in _flat(ops.mpm("info", "gridNodeTags")))
    rest = [t for t in ops.getNodeTags() if t not in grid]
    return rest[0] if len(rest) == 1 else None


def r1_mpm_first_settlement(basis, Kf, stab, dt, nsteps=4, rho=2.0e-6):
    """Row 1: T2's column under the plate (step), quasi-static (rho 2e-6), the plate settlement
    after each step over s0 = q L/(Eoed + Kf/n). The plate node carries the step's increment
    (the grid is rebuilt every step), so U accumulates."""
    ops = fe_engine()
    c = TERZ
    w = 2.0 if basis.startswith("bspline") else 1.0
    k = terzaghi_constants(c["E"], c["nu"], Kf, c["poro"], c["kbar"], c["L"], c["q"])
    _mpm_deck(ops, basis, (int(w), 10), ((0.0, 0.0), (w, c["L"])), rho=rho, Kf=Kf, stab=stab,
              fix=[("xmin", "x"), ("xmax", "x"), ("ymin", "y")], pfix=[("ymax", 0.0)],
              plate=("ymax", c["q"] * w))
    U, pt, hist = 0.0, None, []
    for s_ in range(nsteps):
        ops.mpm("analyze", 1, dt)
        if pt is None:
            pt = _plate_node(ops)
        U += -ops.nodeDisp(pt, 2)
        hist.append({"t": (s_ + 1) * dt, "settle": U, "ratio": U / k["s0"]})
    return {"basis": basis, "Kf": Kf, "stab": stab if isinstance(stab, str) else list(stab),
            "dt": dt, "rho": rho, "u_undrained": k["s0"], "ratio_first_step": hist[0]["ratio"],
            "history": hist}


def _field_at_particles(ops, tmpdir):
    """The grid field sum_a S^p_a(x_p) p_a (+ the dump's pi) at every particle, from
    `mpm debug particleBlocks` (zero pattern scales); nodal p read BEFORE the dumps."""
    nodal = {int(t): ops.nodeDisp(int(t), 3) for t in _flat(ops.mpm("info", "gridNodeTags"))}
    n = len(_flat(ops.mpm("particleData", "m")))
    out = []
    for p in range(n):
        path = str(Path(tmpdir) / ("r2_%d.json" % p))
        ops.mpm("debug", "particleBlocks", p, path, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        d = json.loads(Path(path).read_text())
        out.append(sum(d["Sp"][j] * nodal.get(d["nodes"][j], 0.0) for j in range(d["nen"]))
                   + d.get("pi", 0.0))
    return out


def r2_projection_diffusion(basis, blend, dt=0.02, t_end=10.0, samples=(0.5, 1.0, 2.0, 5.0, 10.0)):
    """Row 2: a motionless, sealed column (every u fixed on gimp's 1-cell column; bspline2's
    interior x = 1 column cannot be fixed by faces, so E = 1e9 there), kbar 1e-4, no gravity,
    stab off, p_init = 5 + 2 cos(pi Y/L). Closed form (Neumann diffusion, c = kbar Kf/n with the
    skeleton held): p = 5 + 2 cos(pi y/L) exp(-c pi^2 t/L^2). Reported, in units of the amplitude
    2: particle p vs closed form, grid field at the particles vs closed form, particle vs field."""
    import tempfile
    ops = fe_engine()
    L, n, Kf, kb = 10.0, 0.4, 4444.0, 1.0e-4
    w = 2.0 if basis.startswith("bspline") else 1.0
    E = 1.0e9 if basis.startswith("bspline") else 1.0e4
    _mpm_deck(ops, basis, (int(w), 10), ((0.0, 0.0), (w, L)), E=E, Kf=Kf, kbar=kb, stab="off",
              fix=[("xmin", "xy"), ("xmax", "xy"), ("ymin", "xy"), ("ymax", "xy")],
              blend=blend)
    X = [_flat(ops.mpm("particleData", "x"))[i:i + 2]
         for i in range(0, len(_flat(ops.mpm("particleData", "x"))), 2)]
    p0 = [5.0 + 2.0 * math.cos(math.pi * x[1] / L) for x in X]
    ops.mpm("particlePressure", *p0)
    cdiff = kb * Kf / n
    rows, t, done = [], 0.0, 0
    with tempfile.TemporaryDirectory() as tmp:
        for ts in samples:
            st = int(round(ts / dt))
            ops.mpm("analyze", st - done, dt)
            done = st
            t = st * dt
            pp = _flat(ops.mpm("particleData", "p"))
            xs = _flat(ops.mpm("particleData", "x"))
            ex = [5.0 + 2.0 * math.cos(math.pi * xs[2 * k + 1] / L)
                  * math.exp(-cdiff * math.pi ** 2 * t / L ** 2) for k in range(len(pp))]
            fld = _field_at_particles(ops, tmp)
            rows.append({"t": t, "amp_exact": 2.0 * math.exp(-cdiff * math.pi ** 2 * t / L ** 2),
                         "particle_vs_exact": max(abs(a - b) for a, b in zip(pp, ex)) / 2.0,
                         "field_vs_exact": max(abs(a - b) for a, b in zip(fld, ex)) / 2.0,
                         "particle_vs_field": max(abs(a - b) for a, b in zip(pp, fld)) / 2.0,
                         "mean_p": sum(pp) / len(pp)})
    return {"basis": basis, "pressureBlend": blend, "E": E, "dt": dt, "c": cdiff, "rows": rows,
            "unit": "the initial amplitude 2 kPa"}


def r4_cost(nsteps=900, dt=0.02):
    """Row 4: wall time per step on T2's column (stab auto 0.25, Newmark 0.6/0.3025, dt 0.02), route
    E gimp (1 x 10) and bspline2 (2 x 10) against LadrunoUP Q4 1 x 10, Linear and Newton."""
    import time
    ops = fe_engine()
    c = TERZ
    out = []
    for algo in ("Linear", "Newton"):
        for basis in ("gimp", "bspline2"):
            w = 2.0 if basis.startswith("bspline") else 1.0
            _mpm_deck(ops, basis, (int(w), 10), ((0.0, 0.0), (w, c["L"])), stab=("auto", 0.25),
                      fix=[("xmin", "x"), ("xmax", "x"), ("ymin", "y")], pfix=[("ymax", 0.0)],
                      plate=("ymax", c["q"] * w), linear=(algo == "Linear"))
            t0 = time.perf_counter()
            ops.mpm("analyze", nsteps, dt)
            out.append({"code": "route E", "basis": basis, "algorithm": algo,
                        "particles": len(_flat(ops.mpm("particleData", "m"))),
                        "ms_per_step": 1e3 * (time.perf_counter() - t0) / nsteps})
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        for j in range(11):
            for i in range(2):
                ops.node(1 + 2 * j + i, float(i), float(j))
        ops.nDMaterial("ElasticIsotropic", 1, c["E"], c["nu"], c["rho"])
        for j in range(10):
            ops.element("LadrunoUP", j + 1, 1 + 2 * j, 2 + 2 * j, 4 + 2 * j, 3 + 2 * j, 1,
                        "-Kf", c["Kf"], "-poro", c["poro"], "-rhoF", c["rhoF"],
                        "-perm", c["kbar"], c["kbar"], "-dynSeepage", "off",
                        "-stab", "auto", 0.25)
        for j in range(11):
            for i in range(2):
                ops.fix(1 + 2 * j + i, 1, 1 if j == 0 else 0, 1 if j == 10 else 0)
        ops.equalDOF(21, 22, 2)
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(21, 0.0, -c["q"], 0.0)
        ops.system("UmfPack")
        ops.numberer("RCM")
        ops.constraints("Transformation")
        if algo == "Linear":
            ops.algorithm("Linear")
        else:
            ops.test("NormDispIncr", 1e-10, 25)
            ops.algorithm("Newton")
        ops.integrator("Newmark", 0.6, 0.3025)
        ops.analysis("Transient")
        t0 = time.perf_counter()
        assert ops.analyze(nsteps, dt) == 0
        out.append({"code": "LadrunoUP Q4", "basis": "1 x 10", "algorithm": algo,
                    "ms_per_step": 1e3 * (time.perf_counter() - t0) / nsteps})
    return out


R3 = {"H": 10.0, "h": 0.5, "g": 9.81, "rho": 2.0, "rhoF": 1.0, "E": 1.0e4, "nu": 0.2,
      "Kf": 4444.0, "poro": 0.4, "kbar": 1.0e-4, "D": 2.0, "V": 0.05, "Tr": 2.0,
      "settle1": (900, 1.0), "translate_dt": 0.1, "settle2": (100, 1.0)}


def r3_cell_crossing(basis, translate=True, stab=("auto", 0.25)):
    """Row 3 (orchestrator's protocol): T4's column, SEALED (every face impervious, no pfix), particle
    p initialised hydrostatic gamma_w (H - Y) by `mpm particlePressure`, the skeleton loaded by its
    own gravity from t = 0 and settled (900 s, dt 1 s: Tv 5), then translated rigidly down 4 cells
    (D = 2 m) by K0's floor, `gridvel floor box 0 -1 W 0 v free -V profile t0 Tr t1 translate y`
    (mpm_soil_interface.md D17), V = 0.05 m/s, ramps Tr = 2 s: |a| = 0.025 m/s^2, inertial pressure
    <= rho_sat a H = 0.5 kPa = 0.51 % gamma_w H; dt 0.1 s through the translation; then settled
    (100 s, dt 1 s). Twin: v = 0, the same timeline. Metric: the offset-free RMS of
    (p - mean p) - (gamma_w (H - Y0) - mean), Y0 the initial particle Y, over gamma_w H."""
    ops = fe_engine()
    c = R3
    H, h = c["H"], c["h"]
    W = 2 * h if basis.startswith("bspline") else h
    gw = c["rhoF"] * c["g"]
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, c["E"], c["nu"], c["rho"])
    ops.mpm("grid", 232, "-origin", 0.0, -5.0, "-h", h, "-cells", int(round(W / h)), 31,
            "-basis", basis)
    ops.mpm("mode", "implicit")
    ops.mpm("body", 1, "box", 0.0, 0.0, W, H, "-ppc", 2, "-mat", 1, "-rho", c["rho"])
    # 232k's D19a carrier: a bonded single-phase elastic block y in [-1.5, 0] under the column; the
    # translating gridvel prescribes the rows 2 to 4 cells below the soil, inside the block, so no
    # prescribed row enters or leaves a soil particle's support (the bare-nodes floor of D19 drops
    # the column ~0.3 m when its top row leaves the moving box: measured, then replaced)
    ops.nDMaterial("ElasticIsotropic", 2, c["E"], c["nu"], c["rho"])
    ops.mpm("body", 2, "box", 0.0, -1.5, W, 0.0, "-ppc", 2, "-mat", 2, "-rho", c["rho"])
    sa = ["-stab", "off"] if stab == "off" else ["-stab", "auto", float(stab[1])]
    ops.mpm("twoPhase", 1, "-Kf", c["Kf"], "-poro", c["poro"], "-rhoF", c["rhoF"],
            "-perm", c["kbar"], c["kbar"], *sa)
    ops.mpm("gravity", 0.0, -c["g"])
    ops.mpm("fix", "xmin", "x")
    ops.mpm("fix", "xmax", "x")
    n1, dt1 = c["settle1"]
    t0 = n1 * dt1
    t1 = t0 + c["D"] / c["V"] + c["Tr"]
    v = -c["V"] if translate else 0.0
    ops.mpm("gridvel", "floor", "box", 0.0, -2.0, W, -1.0, "v", "free", v,
            "profile", t0, c["Tr"], t1, "translate", "y")
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")
    xs = _flat(ops.mpm("particleData", "x"))
    Yall = np.array(xs[1::2])
    soil = Yall > 0.0                      # body 1 (the block's particles lie below y = 0)
    Y0 = Yall[soil]
    ops.mpm("particlePressure", *[gw * (H - y) if y > 0.0 else 0.0 for y in Yall])
    ref = gw * (H - Y0)

    def metric():
        p = np.array(_flat(ops.mpm("particleData", "p")))[soil]
        e = (p - p.mean()) - (ref - ref.mean())
        return float(np.sqrt(np.mean(e * e)) / (gw * H)), float(p.mean() - ref.mean())

    out = {"basis": basis, "translate": translate, "stab": stab if isinstance(stab, str) else list(stab),
           "W": W, "t0": t0, "t1": t1, "protocol": dict(c), "particles": len(Y0),
           "floor": "D19a block: elastic carrier y in [-1.5, 0], gridvel rows y in [-2, -1]"}
    out["metric_initial"] = metric()[0]
    ops.mpm("analyze", n1, dt1)
    out["after_settle1"], out["offset_after_settle1"] = metric()
    nt = int(round((t1 - t0) / c["translate_dt"]))
    worst = 0.0
    hist = []
    for k in range(nt):
        ops.mpm("analyze", 1, c["translate_dt"])
        if k % 10 == 9:
            m_, off = metric()
            worst = max(worst, m_)
            yy = np.array(_flat(ops.mpm("particleData", "x"))[1::2])[soil]
            hist.append({"t": t0 + (k + 1) * c["translate_dt"], "metric": m_, "offset": off,
                         "mean_dy": float(np.mean(yy - Y0))})
    out["max_during_translation"] = worst
    out["translation_history"] = hist
    out["after_translation"] = metric()[0]
    n2, dt2 = c["settle2"]
    ops.mpm("analyze", n2, dt2)
    out["after_settle2"], out["offset_after_settle2"] = metric()
    y = np.array(_flat(ops.mpm("particleData", "x"))[1::2])[soil]
    out["mean_displacement_y"] = float(np.mean(y - Y0))
    out["displacement_spread"] = float(np.ptp(y - Y0))
    return out


# ---------------------------------------------------------------------- R rows, route S
RS_DEFINITION = (
    "u(0+) for the dynamic decks (chosen before running): the plate settlement at the end of a "
    "linear load ramp of T_r = 20 undrained wave transits of the column, T_r = 20 L/c_u, c_u = "
    "sqrt((Eoed + Kf/n)/rho), rho = 2 (T2's): 1.897 s for T2's soil, 0.1204 s for Kf 2.2e6. "
    "ratio = u(T_r)/s0, s0 = q L/(Eoed + Kf/n). Beside it the closed form's own ratio under the same "
    "ramp, s_ramp(T_r)/s0 (Duhamel; the real drainage during the ramp), and LadrunoUP and route E "
    "(mode A) on the same ramp at dt = T_r/50. Route S: mode F at dt = 0.49 x its drained pencil "
    "(just under the default cap s 0.5); settlement read as the particle quadrature "
    "-sum_p (F_yy - 1) V0_p / width (mode F has no plate node). Quasi-static rho = 2e-6 is not "
    "reachable explicitly (the drained pencil scales with sqrt(rho)).")


def _rs_deck(ops, basis, Kf, stab, t_ramp, mode="explicit"):
    """T2's column under the plate, load ramped 0 -> 1 over t_ramp (Path series), route S (mode F)
    or route E (mode A, Newmark 0.6/0.3025, Linear)."""
    c = TERZ
    w = 2.0 if basis.startswith("bspline") else 1.0
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, c["E"], c["nu"], c["rho"])
    ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 1.0, "-cells", int(w), 10, "-basis", basis)
    ops.mpm("mode", mode)
    ops.mpm("body", 1, "box", 0.0, 0.0, w, c["L"], "-ppc", 2, "-mat", 1, "-rho", c["rho"])
    sa = ["-stab", "off"] if stab == "off" else ["-stab", "auto", float(stab[1])]
    ops.mpm("twoPhase", 1, "-Kf", Kf, "-poro", c["poro"], "-rhoF", c["rhoF"],
            "-perm", c["kbar"], c["kbar"], *sa)
    ops.mpm("fix", "xmin", "x")
    ops.mpm("fix", "xmax", "x")
    ops.mpm("fix", "ymin", "y")
    ops.mpm("pfix", "ymax", 0.0)
    ops.mpm("plate", "ymax", c["q"] * w)
    ops.timeSeries("Path", 77, "-time", 0.0, t_ramp, 1.0e6, "-values", 0.0, 1.0, 1.0)
    ops.mpm("timeSeries", 77)
    if mode != "explicit":
        ops.constraints("Transformation")
        ops.numberer("RCM")
        ops.system("UmfPack")
        ops.algorithm("Linear")
        ops.integrator("Newmark", 0.6, 0.3025)
        ops.analysis("Transient")
    return w


def _ramp_consts(Kf):
    c = TERZ
    k = terzaghi_constants(c["E"], c["nu"], Kf, c["poro"], c["kbar"], c["L"], c["q"])
    cu = math.sqrt((k["Eoed"] + k["Qbar"]) / c["rho"])
    tr = 20.0 * c["L"] / cu
    s_cf = terzaghi_settle_ramp(k["cv"] * tr / c["L"] ** 2, k["cv"] * tr / c["L"] ** 2,
                                k["s0"], k["sinf"])
    return k, cu, tr, s_cf


def rs_first_settlement(route, basis, Kf, stab):
    """Row 1 (route S or route E) on the ramp of RS_DEFINITION. Returns the ratio at T_r."""
    import time
    ops = fe_engine()
    k, cu, tr, s_cf = _ramp_consts(Kf)
    if route == "S":
        w = _rs_deck(ops, basis, Kf, stab, tr)
        dtd = _flat(ops.mpm("criticalTimeStep", "drainedPencil"))[0]
        dt = 0.49 * dtd
    else:
        w = _rs_deck(ops, basis, Kf, stab, tr, mode="implicit")
        dt, dtd = tr / 50.0, None
    V0 = np.array(_flat(ops.mpm("particleData", "V")))
    nst = int(math.ceil(tr / dt - 1e-9))
    dt = tr / nst
    t0 = time.perf_counter()
    ops.mpm("analyze", nst, dt)
    wall = time.perf_counter() - t0
    if route == "S":
        F = np.array(_flat(ops.mpm("particleData", "F"))).reshape(-1, 4)
        U = -float(np.sum((F[:, 3] - 1.0) * V0)) / w
    else:
        F = np.array(_flat(ops.mpm("particleData", "F"))).reshape(-1, 4)
        U = -float(np.sum((F[:, 3] - 1.0) * V0)) / w
    return {"route": route, "basis": basis, "Kf": Kf,
            "stab": stab if isinstance(stab, str) else list(stab), "t_ramp": tr, "dt": dt,
            "steps": nst, "dt_drained_pencil": dtd, "u": U, "ratio": U / k["s0"],
            "ratio_closed_form": s_cf / k["s0"], "u_undrained": k["s0"],
            "ms_per_step": 1e3 * wall / nst}


def fe_ramp_settlement(Kf, stab):
    """LadrunoUP Q4 1 x 10, plate, the same ramp, dt = T_r/50: u(T_r)/s0."""
    k, cu, tr, s_cf = _ramp_consts(Kf)
    old = TERZ["Kf"]
    TERZ["Kf"] = Kf
    try:
        r = fe_terzaghi(1, stab, t_ramp=tr, first_small=False, dt=tr / 50.0,
                        tvs=(k["cv"] * tr / TERZ["L"] ** 2,))
    finally:
        TERZ["Kf"] = old
    row = r["rows"][0]
    return {"route": "LadrunoUP", "Kf": Kf, "stab": stab if isinstance(stab, str) else list(stab),
            "t_ramp": tr, "dt": tr / 50.0, "u": row["settle"], "ratio": row["settle"] / k["s0"],
            "ratio_closed_form": row["settle_exact"] / k["s0"]}


def rs_cost_and_dt():
    """Row 2: on T2's column (T2's soil, stab auto 0.25): route S's pencils and the admitted dt at
    s 0.5; wall time per step for route S at 0.49 dt_d, route E and LadrunoUP at dt 0.02 (Linear),
    and per simulated second."""
    import time
    ops = fe_engine()
    out = []
    for basis in ("gimp", "bspline2"):
        _rs_deck(ops, basis, TERZ["Kf"], ("auto", 0.25), 1.0)
        dtd = _flat(ops.mpm("criticalTimeStep", "drainedPencil"))[0]
        dtu = _flat(ops.mpm("criticalTimeStep", "undrainedPencil"))[0]
        dte = _flat(ops.mpm("criticalTimeStep"))[0]
        dt = 0.49 * dtd
        n = 2000
        t0 = time.perf_counter()
        ops.mpm("analyze", n, dt)
        ms = 1e3 * (time.perf_counter() - t0) / n
        out.append({"route": "S", "basis": basis, "dt_drained_pencil": dtd,
                    "dt_undrained_pencil": dtu, "dt_estimator_h_over_c": dte,
                    "dt_admitted_s0.5": 0.5 * dtd, "dt_run": dt, "ms_per_step": ms,
                    "ms_per_simulated_s": ms / dt})
        _rs_deck(ops, basis, TERZ["Kf"], ("auto", 0.25), 1.0, mode="implicit")
        n = 900
        t0 = time.perf_counter()
        ops.mpm("analyze", n, 0.02)
        ms = 1e3 * (time.perf_counter() - t0) / n
        out.append({"route": "E", "basis": basis, "dt_run": 0.02,
                    "dt_admitted": "unconditional (Newmark 0.6/0.3025)", "ms_per_step": ms,
                    "ms_per_simulated_s": ms / 0.02})
    fe = [r for r in r4_cost(nsteps=900, dt=0.02) if r["code"] == "LadrunoUP Q4"
          and r["algorithm"] == "Linear"][0]
    out.append({"route": "LadrunoUP Q4 1 x 10", "dt_run": 0.02, "dt_admitted": "unconditional",
                "ms_per_step": fe["ms_per_step"], "ms_per_simulated_s": fe["ms_per_step"] / 0.02})
    return out


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
        rec["newmark"] = dict(terzaghi=[TERZ["gamma"], TERZ["beta"]],
                              mandel=[MANDEL["gamma"], MANDEL["beta"]])
        rec["rule"] = "cases t6_ref_ramp: band per leg = max(2 x FE error at the run's dt, 0.02)"
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
    if leg == "t1":
        rec = _provenance("t1_fe", "fe_t1")
        rec["case"] = dict(T1, dp_args=t1_dp_args(), meshes=T1_MESHES,
                           particle_values=t1_particle_values.__doc__)
        rec["legs"] = {f"{m}_{mat}_{'off' if st == 'off' else 'auto'}": fe_t1(m, mat, st)
                       for m in T1_MESHES for mat in ("elastic", "dp")
                       for st in ("off", ("auto", 0.25))}
        print(_save("t1_fe", rec))
        return
    if leg == "t0":
        # the oracle: numpy only, no binary (provenance: the git head, not ladrunoBuild)
        rec = {"leg": "t0_oracle", "deck": "tests/mpm_twophase_reference.py::t0_*",
               "date": _dt.date.today().isoformat(),
               "units": "kN, t, m, s, kPa; sigma tension-positive, p compression-positive",
               "source": "tests/mpm_twophase_interface.md §1, §7a, §7b (+ rulings 9064acca3, 3d3699a50, P5a freeze 2ecebf6ae)",
               "case_modeA": T0A, "case_routeS": T0S}
        rec["t0_i"] = t0_particle_dump(T0A)
        rec["t0_ii"] = t0_modeA(T0A)
        rec["t0_iii"] = t0_routeS(T0S)
        rec["t13_translation_gimp_free"] = dict(t13_translation(), case=T13,
                                                expectation="v, sigma', p exact to round-off")
        print(_save("t0_oracle", rec))
        return
    if leg == "settle0":
        rec = _provenance("first_settlement_fe", "fe_first_settlement")
        rec["case"] = dict(TERZ, load="plate step from t = 0")
        rec["definition"] = ("ratio = plate settlement after the first step / s0, s0 = q L/(Eoed + "
                             "Qbar) = the closed-form undrained (compressible-fluid) settlement")
        stabs = ("off", ("auto", 0.1), ("auto", 0.25))
        # quasi-static (rho 2e-6: no inertia) -- the diagnostic's own meaning
        rec["legs"] = [fe_first_settlement(st, dt=dt, rho=2.0e-6) for st in stabs
                       for dt in (1.0e-3, 0.02, 0.2)]
        # REPORTED context: a near-incompressible fluid (Kf 2.2e6, water), where S -> 0 and the
        # stabilization dominates the storage block (TIMs' 1.71 regime)
        rec["legs_stiff_fluid"] = [fe_first_settlement(st, dt=dt, rho=2.0e-6, Kf=2.2e6)
                                   for st in stabs for dt in (1.0e-3, 0.02, 0.2)]
        # T2's dynamic deck (rho 2): the first steps are the compression-wave transient (transit
        # 0.095 s), so the "first settlement" there is inertia, not the undrained response
        rec["legs_dynamic"] = [fe_first_settlement(st, dt=dt) for st in stabs
                               for dt in (0.02, 1.0e-3)]
        print(_save("r_first_settlement_fe", rec))
        return
    if leg == "rrows_e":
        # R rows for route E (reported, never gated): rows 1, 2, 4; written to r_rows_route_e.json
        which = argv[2:] or ["r1", "r2", "r4"]
        p_ = DATA / "r_rows_route_e.json"
        rec = json.loads(p_.read_text(encoding="utf-8")) if p_.exists() else {}
        rec.update(_provenance("r_rows_route_e", "r1_mpm_first_settlement, r2_projection_diffusion, "
                               "r4_cost (+ fe_first_settlement beside)"))
        rec["reported_never_gated"] = True
        stabs = ("off", ("auto", 0.1), ("auto", 0.25))
        if "r1" in which:
            rows = []
            for Kf in (4444.0, 2.2e6):
                for st in stabs:
                    for dt in (1.0e-3, 0.02, 0.2):
                        fe = fe_first_settlement(st, dt=dt, nsteps=4, rho=2.0e-6, Kf=Kf)
                        for basis in ("gimp", "bspline2"):
                            r = r1_mpm_first_settlement(basis, Kf, st, dt)
                            r["ladrunoup_ratio_first_step"] = fe["ratio_first_step"]
                            r["ladrunoup_history"] = [h["ratio"] for h in fe["history"]]
                            rows.append(r)
            rec["r1_first_settlement"] = rows
        if "r2" in which:
            rec["r2_projection_diffusion"] = [r2_projection_diffusion(b, bl)
                                              for b in ("gimp", "bspline2") for bl in (1.0, 0.0)]
        if "r4" in which:
            rec["r4_cost"] = r4_cost()
        if "r3" in which:
            rec["r3_cell_crossing"] = [r3_cell_crossing(b, tr) for b in ("gimp", "bspline2")
                                       for tr in (True, False)]
        if "r5" in which:
            rec["r5_undrained_simple_shear"] = {
                "status": "not run: no drainage switch at particles",
                "reason": ("tests/mpm_twophase_interface.md (P8 `mpm pfix <face> <p0>`; §1-§7c) offers "
                           "no way to release a pfix or change the permeability mid-run, which the "
                           "protocol's drained confinement -> sealed shear needs; per the "
                           "orchestrator's protocol none is invented")}
        print(_save("r_rows_route_e", rec))
        return
    if leg == "rrows_s":
        # R rows for route S (reported, never gated), and the route comparison
        rec = _provenance("r_rows_route_s", "rs_first_settlement, fe_ramp_settlement, "
                          "rs_cost_and_dt")
        rec["reported_never_gated"] = True
        rec["definition_u0"] = RS_DEFINITION
        stabs = ("off", ("auto", 0.1), ("auto", 0.25))
        rows = []
        for Kf in (4444.0, 2.2e6):
            for st in stabs:
                rows.append(fe_ramp_settlement(Kf, st))
                for route in ("S", "E"):
                    for basis in ("gimp", "bspline2"):
                        rows.append(rs_first_settlement(route, basis, Kf, st))
        rec["r1_first_settlement"] = rows
        rec["r2_cost_dt"] = rs_cost_and_dt()

        def pick(route, basis, Kf, stab):
            for r in rows:
                if (r["route"] == route and r.get("basis", basis) == basis and r["Kf"] == Kf
                        and r["stab"] == stab):
                    return r["ratio"]
            return None
        cost = {(r["route"], r.get("basis")): r for r in rec["r2_cost_dt"]}
        rec["r3_route_comparison"] = {
            "first_settlement_T2soil_auto0.25_ramp": {
                "closed_form": [r for r in rows if r["Kf"] == 4444.0][0]["ratio_closed_form"],
                "LadrunoUP": pick("LadrunoUP", None, 4444.0, ["auto", 0.25]),
                "route_E_gimp": pick("E", "gimp", 4444.0, ["auto", 0.25]),
                "route_E_bspline2": pick("E", "bspline2", 4444.0, ["auto", 0.25]),
                "route_S_gimp": pick("S", "gimp", 4444.0, ["auto", 0.25]),
                "route_S_bspline2": pick("S", "bspline2", 4444.0, ["auto", 0.25])},
            "CB_lap_B4 (quoted, orchestrator)": {"route_E": {"auto0.25": 0.134, "off": 0.621},
                                                 "route_S": {"auto0.25": 0.137, "off": 0.614},
                                                 "LadrunoUP": {"auto0.25": 0.1259, "off": 0.6030}},
            "cost_and_dt_T2column": {
                f"{k[0]}{'_' + k[1] if k[1] else ''}": {
                    "ms_per_step": v["ms_per_step"], "ms_per_simulated_s": v["ms_per_simulated_s"],
                    "dt": v.get("dt_admitted_s0.5", v.get("dt_run"))} for k, v in cost.items()}}
        print(_save("r_rows_route_s", rec))
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
