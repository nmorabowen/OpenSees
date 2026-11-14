"""WP-232f (232f-V): the references behind every expected number of tests/mpm_implicit_cases.txt.

Never reads the mode-A C++ (SRC/domain/pattern/ladrunoMPM/ or the particle element). Each number
comes from a closed form, an independent numpy solve, or a fork FE run -- a different code path
from mode A (LadrunoBrick, quad, DruckerPrager through an FE element).

    python3.12 tests/mpm_implicit_reference.py closed          # closed forms (numpy only)
    python3.12 tests/mpm_implicit_reference.py elastica        # cantilever elastica (numpy)
    python3.12 tests/mpm_implicit_reference.py bar             # discrete bar (numpy)
    python3.12 tests/mpm_implicit_reference.py seed            # plate seeding checksums
    python3.12 tests/mpm_implicit_reference.py fe_cantilever_qs
    python3.12 tests/mpm_implicit_reference.py fe_cantilever_dyn
    python3.12 tests/mpm_implicit_reference.py fe_footing
    python3.12 tests/mpm_implicit_reference.py fe_plate
    python3.12 tests/mpm_implicit_reference.py contract1 <out.h5>  # the plate particle set

FE legs need a built fork binary: dist/bin of this worktree, or LADRUNO_BIN=<dir with
opensees.pyd>. Each writes tests/data/mpm/implicit/<leg>.json with the build hash
(ops.ladrunoBuild()), the deck (this file + the leg's function), the date and the numbers;
tests/test_mpm_implicit_reference.py checks those files (no binary needed).
"""
from __future__ import annotations

import datetime as _dt
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
DATA = HERE / "data" / "mpm" / "implicit"
CASES = HERE / "mpm_implicit_cases.txt"

# ----------------------------------------------------------------------------------------------
# The cases file (the same reader shape as tests/mpm_reference.py parse_cases)
# ----------------------------------------------------------------------------------------------


def parse_cases(path: Path = CASES) -> dict:
    """{name: {'lines': [(key, tokens)], 'kv': {key: [tokens, ...]}}} in file order."""
    out, cur, name = {}, None, None
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        tok = line.split()
        if tok[0] == "case":
            if cur is not None:
                raise ValueError(f"case {tok[1]} opens inside {name}")
            name, cur = tok[1], {"lines": [], "kv": {}}
            continue
        if tok[0] == "end":
            if cur is None:
                raise ValueError("'end' outside a case")
            if name in out:
                raise ValueError(f"duplicate case {name}")
            out[name] = cur
            cur = name = None
            continue
        if cur is None:
            raise ValueError(f"line outside a case: {line}")
        cur["lines"].append((tok[0], tok[1:]))
        cur["kv"].setdefault(tok[0], []).append(tok[1:])
    if cur is not None:
        raise ValueError(f"case {name} not closed")
    return out


# ----------------------------------------------------------------------------------------------
# Closed forms
# ----------------------------------------------------------------------------------------------


def alpha_from_phi_txc(phi_deg: float) -> float:
    """Drucker-Prager cone slope (sqrt(J2) + alpha I1 = k) through the triaxial-compression
    meridian of Mohr-Coulomb (Chen & Han); dp_collapse.py's convention."""
    s = math.sin(math.radians(phi_deg))
    return 2.0 * s / (math.sqrt(3.0) * (3.0 - s))


def phi_ps_from_alpha(alpha: float) -> float:
    """Plane-strain Mohr-Coulomb angle matching the cone: tan phi = 3 alpha / sqrt(1 - 12 a^2)."""
    return math.degrees(math.atan(3.0 * alpha / math.sqrt(1.0 - 12.0 * alpha ** 2)))


def prandtl_nq(phi_deg: float) -> float:
    """Prandtl-Reissner N_q = exp(pi tan phi) tan^2(pi/4 + phi/2)."""
    p = math.radians(phi_deg)
    return math.exp(math.pi * math.tan(p)) * math.tan(math.pi / 4 + p / 2) ** 2


def davis_phi(phi_deg: float, psi_deg: float = 0.0) -> float:
    """Davis (1968) reduced angle for non-associated flow: tan phi* = cos psi sin phi / (1 -
    sin phi sin psi) -- context for the non-associated leg (it is NOT the gate)."""
    p, s = math.radians(phi_deg), math.radians(psi_deg)
    return math.degrees(math.atan(math.cos(s) * math.sin(p) / (1.0 - math.sin(p) * math.sin(s))))


def footing_closed() -> dict:
    a = alpha_from_phi_txc(20.0)
    phi = phi_ps_from_alpha(a)
    nq = prandtl_nq(phi)
    q0 = 10.0
    # the cohesion the sigma_y = 0.2 kPa apex regularizer implies in the plane-strain match:
    # sqrt(J2) + alpha I1 = k, k = sigma_y / sqrt(3) (UW DruckerPrager: ||s|| + rho I1 =
    # sqrt(2/3) sigma_y, rho = sqrt2 alpha); c = k sqrt(9 + 12 tan^2 phi) / 3 (Chen & Han)
    k = 0.2 / math.sqrt(3.0)
    c = k * math.sqrt(9.0 + 12.0 * math.tan(math.radians(phi)) ** 2) / 3.0
    nc = (nq - 1.0) / math.tan(math.radians(phi))
    return dict(alpha=a, rho_dp=math.sqrt(2.0) * a, phi_ps=phi, Nq=nq, q_exact=q0 * nq,
                band=[0.9 * q0 * nq, 1.1 * q0 * nq], c_regularizer=c, cNc=c * nc,
                cNc_fraction=c * nc / (q0 * nq), phi_davis=davis_phi(phi),
                q_davis=q0 * prandtl_nq(davis_phi(phi)))


# Howland (1930), Phil. Trans. R. Soc. A 229:49-86: sigma_max / sigma_gross for a strip of width W
# with a central circular hole of diameter d under tension, as tabulated in Pilkey & Pilkey,
# Peterson's Stress Concentration Factors, 3rd ed. (2008), Chart 4.1. K_tn = K_tg (1 - d/W).
HOWLAND_KTG = {0.0: 3.00, 0.1: 3.03, 0.2: 3.14, 0.3: 3.36, 0.4: 3.74, 0.5: 4.32}


def howland_ktn(d_over_w: float) -> float:
    return HOWLAND_KTG[d_over_w] * (1.0 - d_over_w)


def heywood_ktn(d_over_w: float) -> float:
    """Heywood's empirical fit K_tn = 2 + (1 - d/W)^3 -- a second, independent source."""
    return 2.0 + (1.0 - d_over_w) ** 3


def bar_closed() -> dict:
    L, E, rho = 25.0, 100.0, 1.0
    c = math.sqrt(E / rho)
    return dict(L=L, c=c, period=4.0 * L / c)


def bar_discrete(n_el=25, dt=0.01, lumped=True) -> dict:
    """An independent numpy solve of the bar's discretization gap: linear 2-node bar elements,
    h = 1, fixed at x = 0, free at L = 25; first-mode period of the semi-discrete system, then
    the Newmark average-acceleration period elongation T_n = T * (w dt / 2) / atan(w dt / 2)."""
    L, E, rho = 25.0, 100.0, 1.0
    h = L / n_el
    K = np.zeros((n_el + 1, n_el + 1))
    M = np.zeros_like(K)
    ke = E / h * np.array([[1, -1], [-1, 1]])
    me = rho * h * (np.diag([0.5, 0.5]) if lumped else np.array([[1 / 3, 1 / 6], [1 / 6, 1 / 3]]))
    for e in range(n_el):
        K[e:e + 2, e:e + 2] += ke
        M[e:e + 2, e:e + 2] += me
    K, M = K[1:, 1:], M[1:, 1:]
    w = math.sqrt(min(np.linalg.eigvals(np.linalg.solve(M, K)).real))
    T = 2 * math.pi / w
    x = w * dt / 2
    return dict(n_el=n_el, lumped=lumped, period_semidiscrete=T,
                period_newmark=T * x / math.atan(x), period_exact=4 * L / math.sqrt(E / rho))


def cantilever_linear(E: float, rho=1000.0, g=10.0, L=10.0, H=1.0) -> dict:
    q = rho * g * H
    EI = E * H ** 3 / 12.0
    mu = rho * H
    w1 = 1.8751040687119611 ** 2 * math.sqrt(EI / (mu * L ** 4))
    return dict(qL3_EI=q * L ** 3 / EI, tip_static_linear=q * L ** 4 / (8 * EI),
                period1_linear=2 * math.pi / w1)


# ----------------------------------------------------------------------------------------------
# Elastica: inextensible Euler-Bernoulli cantilever under a distributed weight q per unit length
# ----------------------------------------------------------------------------------------------


def elastica_tip(qL3_EI: float, n: int = 4000) -> dict:
    """theta = downward slope angle, s = arc length from the root, z downward:
    x' = cos theta, z' = sin theta, M(s) = int_s^L q (x(xi) - x(s)) dxi (hogging),
    theta' = M / EI, so theta'' = -(q/EI) (L - s) cos theta; theta(0) = 0, theta'(L) = 0.
    Nondimensional (L = 1, lam = q L^3 / EI). Shooting on theta'(0) (bisection; theta'(1) grows
    with theta'(0)) with RK4 over n steps. Returns the tip deflection / L (downward positive)."""
    lam = qL3_EI

    def f(s, y):
        th = y[0]
        return np.array([y[1], -lam * (1.0 - s) * math.cos(th), math.cos(th), math.sin(th)])

    def shoot(k0):
        h = 1.0 / n
        y = np.array([0.0, k0, 0.0, 0.0])          # theta, theta', x, z
        s = 0.0
        for _ in range(n):
            k1 = f(s, y)
            k2 = f(s + h / 2, y + h / 2 * k1)
            k3 = f(s + h / 2, y + h / 2 * k2)
            k4 = f(s + h, y + h * k3)
            y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            s += h
        return y

    lo, hi = 0.0, lam
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if shoot(mid)[1] > 0.0:
            hi = mid
        else:
            lo = mid
    y = shoot(0.5 * (lo + hi))
    return dict(tip_down_over_L=float(y[3]), tip_x_over_L=float(y[2]), theta_tip=float(y[0]),
                kappa0=0.5 * (lo + hi))


# ----------------------------------------------------------------------------------------------
# The fork FE engine (a different code path from mode A)
# ----------------------------------------------------------------------------------------------
_OPS = None


def fe_engine():
    """Import opensees.pyd from LADRUNO_BIN or this worktree's dist/bin (never an installed one
    picked up through a site .pth: the path is checked)."""
    global _OPS
    if _OPS is not None:
        return _OPS
    b = os.environ.get("LADRUNO_BIN") or str(REPO / "dist" / "bin")
    if not (Path(b) / "opensees.pyd").exists():
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
                deck=f"tests/mpm_implicit_reference.py::{deck}",
                date=_dt.date.today().isoformat())


def _save(name: str, payload: dict) -> Path:
    DATA.mkdir(parents=True, exist_ok=True)
    p = DATA / f"{name}.json"
    payload["output"] = str(p.relative_to(REPO)).replace("\\", "/")
    p.write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    return p


def _slab(nx, ny, L, H):
    """A plane-strain slab: (nx x ny) hexes over [0,L] x [0,H], one element (thickness H/ny)
    out of plane, z fixed everywhere. Node id = 1 + i + (nx+1)(j + (ny+1) k)."""
    xs, ys = np.linspace(0, L, nx + 1), np.linspace(0, H, ny + 1)
    t = H / ny

    def nid(i, j, k):
        return 1 + i + (nx + 1) * (j + (ny + 1) * k)

    nodes = {nid(i, j, k): (xs[i], ys[j], k * t)
             for k in (0, 1) for j in range(ny + 1) for i in range(nx + 1)}
    hexes = []
    for j in range(ny):
        for i in range(nx):
            hexes.append([nid(i, j, 0), nid(i + 1, j, 0), nid(i + 1, j + 1, 0), nid(i, j + 1, 0),
                          nid(i, j, 1), nid(i + 1, j, 1), nid(i + 1, j + 1, 1), nid(i, j + 1, 1)])
    return xs, ys, t, nid, nodes, hexes


def _sample_u(ops, xs, ys, nid, pts, dof=2):
    """Mean over pts of the bilinear interpolation of nodal u_dof (layer k = 0)."""
    out = []
    for X, Y in pts:
        i = min(max(int(np.searchsorted(xs, X) - 1), 0), len(xs) - 2)
        j = min(max(int(np.searchsorted(ys, Y) - 1), 0), len(ys) - 2)
        r = (X - xs[i]) / (xs[i + 1] - xs[i])
        s = (Y - ys[j]) / (ys[j + 1] - ys[j])
        u = [ops.nodeDisp(nid(a, b, 0), dof)
             for a, b in ((i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1))]
        out.append((1 - r) * (1 - s) * u[0] + r * (1 - s) * u[1] + r * s * u[2]
                   + (1 - r) * s * u[3])
    return float(np.mean(out))


CANT = dict(L=10.0, H=1.0, rho=1000.0, g=10.0)
MPM_H = {"r1": 0.25, "r2": 0.125, "r3": 0.0625}


def tip_points(h, H=1.0, L=10.0, ppc=2):
    """The MPM tip particle column: x = L - h/(2 ppc), y = (k + 1/2) h/ppc (lattice rule)."""
    n = int(round(H / (h / ppc)))
    return [(L - h / (2 * ppc), (k + 0.5) * h / ppc) for k in range(n)]


def _cantilever_model(ops, ny, E, geom, form, lumped=False):
    L, H, rho, g = CANT["L"], CANT["H"], CANT["rho"], CANT["g"]
    nx = int(round(ny * L / H))
    xs, ys, t, nid, nodes, hexes = _slab(nx, ny, L, H)
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for n, (x, y, z) in nodes.items():
        ops.node(n, float(x), float(y), float(z))
    ops.nDMaterial("ElasticIsotropic", 1, E, 0.0, rho)
    extra = ["-geom", geom]
    if form != "std":
        extra += ["-formulation", form]
    if lumped:
        extra += ["-lumped"]
    for e, c in enumerate(hexes, start=1):
        ops.element("LadrunoBrick", e, *c, 1, *extra, "-b", 0.0, -rho * g, 0.0)
    for n, (x, y, z) in nodes.items():
        ops.fix(n, *((1, 1, 1) if abs(x) < 1e-12 else (0, 0, 1)))
    return xs, ys, nid, len(hexes)


def fe_cantilever_qs(nys=(4, 8, 16, 32), geom="hypo", form="std"):
    """Quasi-static gravity ramp in 20 LoadControl steps (as the MPM case); the body load goes
    through eleLoad -brickW so the ramp scales it (LadrunoBrick applies -b unconditionally
    otherwise). Output: the tip-column mean u_y at each MPM resolution's particle points."""
    ops = fe_engine()
    E = 6.0e7
    res = []
    for ny in nys:
        t0 = time.time()
        xs, ys, nid, ne = _cantilever_model(ops, ny, E, geom, form)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.eleLoad("-ele", *range(1, ne + 1), "-type", "-brickW")
        ops.constraints("Plain")
        ops.numberer("RCM")
        ops.system("UmfPack")
        ops.test("NormDispIncr", 1e-10, 50, 0)
        ops.algorithm("Newton")
        ops.integrator("LoadControl", 0.05)
        ops.analysis("Static")
        if ops.analyze(20) != 0:
            raise RuntimeError(f"cantilever QS ny={ny} {geom}/{form} failed")
        tips = {r: _sample_u(ops, xs, ys, nid, tip_points(h)) for r, h in MPM_H.items()}
        tipL = _sample_u(ops, xs, ys, nid, [(10.0, y) for y in np.linspace(0, 1, 2 * ny + 1)])
        res.append(dict(ny=ny, nx=int(round(10 * ny)), hexes=ne, tip=tips, tip_at_L=tipL,
                        wall_s=round(time.time() - t0, 1)))
        print(f"QS {geom}/{form} ny={ny}: tip(L) {tipL:.6f}  "
              + " ".join(f"{r} {v:.6f}" for r, v in tips.items()), flush=True)
    return res


def extrema_parabola(t, u):
    """Local minima of a sampled history, each refined by the parabola through the minimum
    sample and its two neighbours. Returns [(t*, u*), ...]."""
    out = []
    for k in range(1, len(u) - 1):
        if u[k] < u[k - 1] and u[k] <= u[k + 1]:
            a, b, c = u[k - 1], u[k], u[k + 1]
            den = a - 2 * b + c
            d = 0.5 * (a - c) / den if den != 0 else 0.0
            out.append((t[k] + d * (t[k + 1] - t[k]), b - 0.25 * (a - c) * d))
    return out


def fe_cantilever_dyn(nys=(4, 8, 16, 32), geom="hypo", form="std", lumped=False, dt=0.01,
                      nsteps=400):
    """Suddenly applied gravity (constant -b, no pattern), Newmark 0.5 0.25, from rest."""
    ops = fe_engine()
    E = 1.2e8
    res = []
    for ny in nys:
        t0 = time.time()
        xs, ys, nid, ne = _cantilever_model(ops, ny, E, geom, form, lumped)
        ops.constraints("Plain")
        ops.numberer("RCM")
        ops.system("UmfPack")
        ops.test("NormDispIncr", 1e-10, 50, 0)
        ops.algorithm("Newton")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient")
        pts = {r: tip_points(h) for r, h in MPM_H.items()}
        hist = {r: [0.0] for r in MPM_H}
        tt = [0.0]
        for k in range(nsteps):
            if ops.analyze(1, dt) != 0:
                raise RuntimeError(f"cantilever dyn ny={ny} failed at step {k + 1}")
            tt.append(ops.getTime())
            for r in MPM_H:
                hist[r].append(_sample_u(ops, xs, ys, nid, pts[r]))
        out = {}
        for r in MPM_H:
            mins = extrema_parabola(np.array(tt), np.array(hist[r]))
            out[r] = dict(first_peak=mins[0][1], first_period=mins[1][0] - mins[0][0],
                          t_min=[m[0] for m in mins[:3]])
        res.append(dict(ny=ny, hexes=ne, lumped=lumped, dt=dt, tip=out,
                        wall_s=round(time.time() - t0, 1)))
        print(f"DYN {geom}/{form} lumped={lumped} ny={ny}: " + " ".join(
            f"{r} peak {v['first_peak']:.6f} T {v['first_period']:.6f}" for r, v in out.items()),
            flush=True)
    return res


# ----------------------------------------------------------------------------------------------
# Footing: the ADR-79 mild-cone control on the MPM case's own domain and staging
# ----------------------------------------------------------------------------------------------
FOOT = dict(B=2.0, X=12.0, D=6.0, q0=10.0, K=1.5e5, nu=0.45, sy=0.2, ds=0.002, s_max=0.2)


def fe_footing(n_across_b=4, rough=True, geom="linear", form="bbar", ds=None, s_max=None):
    """Half model, x in [0, 12], y in [-6, 0], uniform h = B / n_across_b, plane-strain slab of
    LadrunoBrick (H8, B-bar, small-strain geometry: the ADR-79 control's element) with the case's
    DruckerPrager. Stage 1: the surcharge ramps on every top node (tributary) in 10 steps.
    Stage 2: the footing nodes (top, x <= B/2) are pushed down 2 mm per step to s = 0.2 m; rough =
    u_x held (as the MPM gridvel), else free (the deck's smooth footing). q = (-sum R_y + the
    surcharge lumped on the footing nodes) / (B/2 t): dp_strip.py's r_corr convention, the MPM
    case's cond_force. Each step: a retry ladder (KrylovNewton, NewtonLineSearch, relaxed
    KrylovNewton), then halving."""
    ops = fe_engine()
    f = FOOT
    ds = ds or f["ds"]
    s_max = s_max or f["s_max"]
    h = f["B"] / n_across_b
    nx, ny = int(round(f["X"] / h)), int(round(f["D"] / h))
    xs, ys, t, nid, nodes, hexes = _slab(nx, ny, f["X"], f["D"])
    ys = ys - f["D"]
    a = alpha_from_phi_txc(20.0)
    G = 3 * f["K"] * (1 - 2 * f["nu"]) / (2 * (1 + f["nu"]))
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for n, (x, y, z) in nodes.items():
        ops.node(n, float(x), float(y - f["D"]), float(z))
    ops.nDMaterial("DruckerPrager", 1, f["K"], G, f["sy"], math.sqrt(2) * a, 0.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    extra = ["-geom", geom] + (["-formulation", form] if form != "std" else [])
    for e, c in enumerate(hexes, start=1):
        ops.element("LadrunoBrick", e, *c, 1, *extra)
    for n, (x, y, z) in nodes.items():
        y = y - f["D"]
        if abs(y + f["D"]) < 1e-9:
            ops.fix(n, 1, 1, 1)
        else:
            ops.fix(n, 1 if (abs(x) < 1e-9 or abs(x - f["X"]) < 1e-9) else 0, 0, 1)
    top = [(i, k) for k in (0, 1) for i in range(nx + 1)]
    trib = {nid(i, ny, k): (h * t / 4.0) * (1 if i in (0, nx) else 2) for i, k in top}
    foot = [nid(i, ny, k) for i, k in top if xs[i] <= f["B"] / 2 + 1e-9]
    r_corr = f["q0"] * sum(trib[n] for n in foot)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n, A in trib.items():
        ops.load(n, 0.0, -f["q0"] * A, 0.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    tol = 1e-6
    ops.test("NormUnbalance", tol, 40, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.1)
    ops.analysis("Static")
    if ops.analyze(10) != 0:
        raise RuntimeError("footing surcharge stage failed")
    ops.loadConst("-time", 0.0)
    uy0 = {n: ops.nodeDisp(n, 2) for n in foot}
    ux0 = {n: ops.nodeDisp(n, 1) for n in foot}
    for k, n in enumerate(foot):
        ops.timeSeries("Path", 100 + k, "-time", 0.0, 1.0, "-values", uy0[n], uy0[n] - s_max,
                       "-useLast")
        ops.pattern("Plain", 100 + k, 100 + k)
        ops.sp(n, 2, 1.0)
        if rough:
            ops.sp(n, 1, ux0[n])
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.analysis("Static")
    ladder = [("KrylovNewton", tol, 40), ("NewtonLineSearch", tol, 60),
              ("KrylovNewton", 10 * tol, 80)]
    rows, dsc, nfail, nrelax, t0 = [], ds, 0, 0, time.time()
    while True:
        s_now = ops.getTime() * s_max
        if s_now >= s_max - 1e-12:
            break
        dsc = min(dsc, s_max - s_now)
        ops.integrator("LoadControl", dsc / s_max)
        ok = False
        for i, (alg, tl, it) in enumerate(ladder):
            ops.test("NormUnbalance", tl, it, 0)
            ops.algorithm(alg)
            if ops.analyze(1) == 0:
                ok, nrelax = True, nrelax + (i == 2)
                break
            nfail += 1
        if not ok:
            dsc *= 0.5
            if dsc < 1e-6:
                break
            continue
        dsc = min(2 * dsc, ds)
        ops.reactions()
        R = -sum(ops.nodeReaction(n, 2) for n in foot) + r_corr
        rows.append((ops.getTime() * s_max, R / (f["B"] / 2 * t)))
    s = np.array([r[0] for r in rows])
    q = np.array([r[1] for r in rows])
    out = dict(n_across_b=n_across_b, h=h, hexes=len(hexes), rough=rough, geom=geom, form=form,
               q_max=float(q.max()), s_at_qmax=float(s[q.argmax()]), s_end=float(s[-1]),
               q_plateau=float(q[s >= 0.17 - 1e-12].mean()) if s[-1] >= 0.17 else None,
               q_at={f"{sb:.3f}": float(np.interp(sb * f["B"], s, q))
                     for sb in (0.01, 0.02, 0.05, 0.075, 0.1) if sb * f["B"] <= s[-1]},
               steps=len(rows), failed_attempts=nfail, relaxed=nrelax,
               wall_s=round(time.time() - t0, 1))
    print(f"FOOT n={n_across_b} rough={rough}: q_max {out['q_max']:.3f} at s {out['s_at_qmax']:.4f}"
          f" plateau {out['q_plateau']} end s {out['s_end']:.4f} ({out['wall_s']} s)", flush=True)
    return out


# ----------------------------------------------------------------------------------------------
# Plate with a hole, d/W = 0.5: the O-grid seeding rule of the case, and the FE reference
# ----------------------------------------------------------------------------------------------
PLATE = dict(a=1.0, half_w=2.0, top=4.0, h=0.03125, n_theta=256, n_r=120, E=2.0e5, nu=0.25,
             strain=1.0e-4)
_GL8 = np.polynomial.legendre.leggauss(8)


def _rho_o(th, half_w):
    return half_w / np.maximum(np.cos(th), np.sin(th))


def holeplate_particles(a=1.0, half_w=2.0, n_theta=256, n_r=120):
    """The case's body 'holeplate': particle centres, exact cell areas, GIMP half-lengths.
    Ids i + n_theta j (theta fastest)."""
    if n_theta % 2:
        raise ValueError("n_theta must be even (no cell may straddle pi/4)")
    dth, ds = 0.5 * math.pi / n_theta, 1.0 / n_r
    th = (np.arange(n_theta) + 0.5) * dth
    s = (np.arange(n_r) + 0.5) * ds
    T, S = np.meshgrid(th, s, indexing="xy")          # shape (n_r, n_theta): theta fastest
    R = a + S * (_rho_o(T, half_w) - a)
    pos = np.stack([R * np.cos(T), R * np.sin(T)], axis=-1).reshape(-1, 2)
    # exact area of [th +- dth/2] x [s +- ds/2]: int dth (r(s+)^2 - r(s-)^2)/2, Gauss in theta
    xg, wg = _GL8
    vol = np.zeros((n_r, n_theta))
    for x, w in zip(xg, wg):
        tq = th + 0.5 * dth * x                         # (n_theta,)
        ro = _rho_o(tq, half_w)
        rp = a + (s[:, None] + 0.5 * ds) * (ro[None, :] - a)
        rm = a + (s[:, None] - 0.5 * ds) * (ro[None, :] - a)
        vol += 0.5 * w * 0.5 * dth * (rp ** 2 - rm ** 2)
    vol = vol.reshape(-1)
    # lp = sqrt(V)/2, clipped so no GIMP domain leaves the square [0, half_w]^2 (Contract 1,
    # C1-POS): only the cells next to the symmetry axes and the outer edges are clipped
    l0 = 0.5 * np.sqrt(vol)
    lp = np.stack([np.minimum.reduce([l0, pos[:, 0], half_w - pos[:, 0]]),
                   np.minimum.reduce([l0, pos[:, 1], half_w - pos[:, 1]])], axis=1)
    return pos, vol, lp


def lattice_box(lo, hi, origin, h, ppc=2):
    """tests/mpm_interface.md §5: sub-lattice points x0 + ((a + 0.5) * h) / n inside the closed
    box, x fastest. 2D."""
    pts = []
    na = [int(math.floor((hi[d] - origin[d]) / (h / ppc))) + 2 for d in range(2)]
    ys = [origin[1] + ((j + 0.5) * h) / ppc for j in range(na[1])]
    xs = [origin[0] + ((i + 0.5) * h) / ppc for i in range(na[0])]
    for y in ys:
        if lo[1] <= y <= hi[1]:
            for x in xs:
                if lo[0] <= x <= hi[0]:
                    pts.append((x, y))
    return np.array(pts)


def plate_particle_set(p=PLATE):
    pos1, vol1, lp1 = holeplate_particles(p["a"], p["half_w"], p["n_theta"], p["n_r"])
    pos2 = lattice_box((0.0, p["half_w"]), (p["half_w"], p["top"]), (0.0, 0.0), p["h"])
    vol2 = np.full(len(pos2), (p["h"] / 2) ** 2)
    lp2 = np.full((len(pos2), 2), p["h"] / 4)
    pos = np.vstack([pos1, pos2])
    vol = np.r_[vol1, vol2]
    lp = np.vstack([lp1, lp2])
    body = np.r_[np.full(len(pos1), 1), np.full(len(pos2), 2)]
    return pos, vol, lp, body


def write_plate_contract1(path, p=PLATE) -> None:
    """The plate's particle set as ONE Contract-1 file, through the writer of
    tests/mpm_contract_fixtures.py (the way apeGmsh writes it). Units N, m, s; zero initial
    stress; density 1 (quasi-static)."""
    import mpm_contract_fixtures as F  # noqa: E402  (h5py)
    import mpm_contracts as MC  # noqa: E402
    pos, vol, lp, body = plate_particle_set(p)
    n = len(pos)
    cells = [int(round(p["half_w"] / p["h"])), int(round(p["top"] / p["h"])) + 1]
    model = {
        "attrs": {
            "FORMAT": MC.FORMAT, "FORMAT_VERSION": np.int32(1), "GENERATOR": "wp232f-reference",
            "SET_ID": "232f0000-0000-4000-8000-0000000000d5", "CREATED_UTC": "2026-11-14T00:00:00Z",
            "NDM": np.int32(2), "NUM_PARTICLES": np.int64(n), "PLANE_STRAIN_THICKNESS": np.float64(1.0),
            "UNIT_LENGTH": "m", "UNIT_FORCE": "N", "UNIT_TIME": "s", "UNIT_SYSTEM": "N_m_s",
            "MASS_SOURCE": "density", "STRESS_SIGN": "tension-positive", "VOIGT_ORDER": MC.VOIGT[2],
        },
        "grid": {"ORIGIN": np.array([0.0, 0.0]), "CELL_SIZE": np.array([p["h"], p["h"]]),
                 "NUM_CELLS": np.array(cells, dtype=np.int32), "BASIS": "gimp"},
        "bodies": [(1, "holeplate"), (2, "upper")],
        "provenance": {"apegmsh_version": "none", "apegmsh_commit": "unknown",
                       "session_id": "232f-V plate_hole_dw050"},
        "seeding": {"METHOD": "holeplate-ogrid+lattice",
                    "PARTICLES_PER_CELL": np.array([2, 2], dtype=np.int32)},
        "particles": {
            "ID": np.arange(n, dtype="<i8"),
            "POSITION": pos.astype("<f8"),
            "VOLUME": vol.astype("<f8"),
            "DENSITY": np.ones(n),
            "GIMP_HALF_LENGTH": lp.astype("<f8"),
            "STRESS": np.zeros((n, 3)),
            "STRESS33": np.zeros(n),
            "BODY_ID": body.astype("<i4"),
            "MATERIAL_TAG": np.ones(n, dtype="<i4"),
        },
    }
    F.write_contract1(path, model)


def plate_seed_checksums(p=PLATE) -> dict:
    pos, vol, lp, body = plate_particle_set(p)
    return dict(count=int(len(pos)), volume_sum=float(vol.sum()),
                exact_area=p["half_w"] * p["top"] - math.pi * p["a"] ** 2 / 4,
                moment=[float((vol * pos[:, 0]).sum()), float((vol * pos[:, 1]).sum())],
                max_lp_over_h=float(lp.max() / p["h"]), min_volume=float(vol.min()))


def fe_plate(n_theta=256, n_r=120, n_up=None, top=4.0, load="disp", a=1.0, half_w=2.0,
             nu=0.25):
    """Fork FE (quad, PlaneStrain ElasticIsotropic, small strain) on the O-grid mesh whose CELLS
    are the case's particle cells (corner nodes at the theta / s edges), plus a block above with
    the O-grid's top-edge x stations. load = 'disp': top edge u_y prescribed (as the MPM gridvel);
    'traction': uniform end traction (Howland's condition). Returns K_tn = peak sigma_yy on the
    ligament at the hole edge / sigma_net, the peak from a quadratic fit of the Gauss-point
    sigma_yy of the innermost element row nearest theta = 0, extrapolated to r = a."""
    ops = fe_engine()
    t0 = time.time()
    dth, ds = 0.5 * math.pi / n_theta, 1.0 / n_r
    the = np.arange(n_theta + 1) * dth
    se = np.arange(n_r + 1) * ds
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    tag = {}

    def node(key, x, y):
        if key not in tag:
            tag[key] = len(tag) + 1
            ops.node(tag[key], float(x), float(y))
        return tag[key]

    for j in range(n_r + 1):
        for i in range(n_theta + 1):
            th = the[i]
            ro = half_w / max(math.cos(th), math.sin(th))
            if i == n_theta:
                ro = half_w
            r = a + se[j] * (ro - a)
            x, y = r * math.cos(th), r * math.sin(th)
            if i == n_theta:
                x = 0.0
            if i == 0:
                y = 0.0
            if j == n_r and th >= math.pi / 4:
                y = half_w
            if j == n_r and th <= math.pi / 4:
                x = half_w
            node(("o", i, j), x, y)
    ops.nDMaterial("ElasticIsotropic", 1, PLATE["E"], nu, 1.0)
    e = 0
    elems = {}
    for j in range(n_r):
        for i in range(n_theta):
            e += 1
            elems[e] = (i, j)
            # counter-clockwise: +s (outward) before +theta
            ops.element("quad", e, tag[("o", i, j)], tag[("o", i, j + 1)], tag[("o", i + 1, j + 1)],
                        tag[("o", i + 1, j)], 1.0, "PlaneStrain", 1)
    # upper block: x stations from the O-grid top edge (theta in [pi/4, pi/2]), reversed to
    # increasing x; y uniform with about the same spacing
    top_i = list(range(n_theta, n_theta // 2 - 1, -1))
    xst = [ops.nodeCoord(tag[("o", i, n_r)])[0] for i in top_i]
    n_up = n_up or int(round((top - half_w) / (half_w / (n_theta // 2))))
    yst = np.linspace(half_w, top, n_up + 1)
    for jj, y in enumerate(yst):
        for ii, i in enumerate(top_i):
            if jj == 0:
                tag[("u", ii, 0)] = tag[("o", i, n_r)]
            else:
                node(("u", ii, jj), xst[ii], y)
    for jj in range(n_up):
        for ii in range(len(top_i) - 1):
            e += 1
            ops.element("quad", e, tag[("u", ii, jj)], tag[("u", ii + 1, jj)],
                        tag[("u", ii + 1, jj + 1)], tag[("u", ii, jj + 1)], 1.0, "PlaneStrain", 1)
    for n in sorted(set(tag.values())):          # ("u", ii, 0) aliases the O-grid top nodes
        x, y = ops.nodeCoord(n)
        fx = 1 if abs(x) < 1e-12 else 0
        fy = 1 if abs(y) < 1e-12 else 0
        if fx or fy:
            ops.fix(n, fx, fy)
    topn = [tag[("u", ii, n_up)] for ii in range(len(top_i))]
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    if load == "disp":
        for n in topn:
            ops.sp(n, 2, PLATE["strain"] * top)
    else:
        xs = np.array(xst)
        for ii, n in enumerate(topn):
            w = 0.0
            if ii > 0:
                w += 0.5 * (xs[ii] - xs[ii - 1])
            if ii < len(xs) - 1:
                w += 0.5 * (xs[ii + 1] - xs[ii])
            ops.load(n, 0.0, 1.0 * w)
    ops.constraints("Transformation" if load == "disp" else "Plain")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-12, 10, 0)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    if ops.analyze(1) != 0:
        raise RuntimeError("plate FE failed")
    ops.reactions()
    if load == "disp":
        Ry = sum(ops.nodeReaction(n, 2) for n in topn)   # the support pulls: positive
    else:
        Ry = half_w * 1.0
    s_net = Ry / (half_w - a)
    # Gauss-point sigma_yy of the innermost elements of the theta = 0 row, extrapolated to r = a
    rr, sy = [], []
    for j in range(4):
        etag = 1 + 0 + n_theta * j
        st = np.array(ops.eleResponse(etag, "stresses")).reshape(4, -1)
        crd = np.array([ops.nodeCoord(tag[("o", i2, j2)]) for i2, j2 in
                        ((0, j), (0, j + 1), (1, j + 1), (1, j))])
        g = 1 / math.sqrt(3)
        for (xi, eta), srow in zip(((-g, -g), (g, -g), (g, g), (-g, g)), st):
            N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                                 (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
            xg = N @ crd
            if xg[1] < 0.5 * (crd[:, 1].max()):            # the GP row nearest theta = 0
                rr.append(math.hypot(*xg))
                sy.append(srow[1])
    rr, sy = np.array(rr), np.array(sy)
    c = np.polyfit(rr - a, sy, 2)
    peak = float(np.polyval(c, 0.0))
    out = dict(n_theta=n_theta, n_r=n_r, n_up=n_up, top=top, load=load, nu=nu,
               elements=e, sigma_net=s_net, peak_syy=peak, Ktn=peak / s_net,
               # the sampling floor: sigma_yy / sigma_net at the innermost particle centre of
               # the case's seeding (n_theta 256, n_r 120): r = a + (rho_o(dth/2) - a) / 240
               Ktn_at_first_particle=float(np.polyval(c, (_rho_o(math.pi / 1024, half_w) - a)
                                                      / 240.0) / s_net),
               Ktn_innermost_gp=float(sy[np.argmin(rr)] / s_net), r_innermost_gp=float(rr.min()),
               wall_s=round(time.time() - t0, 1))
    print(f"PLATE n_theta={n_theta} n_r={n_r} top={top} {load}: K_tn {out['Ktn']:.5f} "
          f"(innermost GP {out['Ktn_innermost_gp']:.5f} at r {out['r_innermost_gp']:.5f}) "
          f"{e} quads {out['wall_s']} s", flush=True)
    return out


# ----------------------------------------------------------------------------------------------
# Assembling the references (each leg writes tests/data/mpm/implicit/<leg>.json)
# ----------------------------------------------------------------------------------------------


def richardson(fine: float, coarse: float, ratio: float = 2.0, p: float = 2.0) -> float:
    return fine + (fine - coarse) / (ratio ** p - 1.0)


def observed_order(a: float, b: float, c: float, ratio: float = 2.0) -> float:
    """Order from three successive values on meshes refined by `ratio`."""
    return math.log(abs((a - b) / (b - c))) / math.log(ratio)


def leg_cantilever_qs() -> dict:
    hypo = fe_cantilever_qs((4, 8, 16, 32, 64), "hypo", "std")
    corot = fe_cantilever_qs((8, 16, 32), "corot", "std")
    ref = {r: richardson(hypo[-1]["tip"][r], hypo[-2]["tip"][r]) for r in MPM_H}
    el = elastica_tip(2.0)
    return dict(**_provenance("fe_cantilever_qs", "leg_cantilever_qs"),
                element="LadrunoBrick -geom hypo -formulation std, plane-strain slab "
                        "(u_z = 0), ElasticIsotropic E = 6e7 nu = 0, square elements",
                hypo=hypo, corot=corot, elastica=el,
                elastica_tip=-el["tip_down_over_L"] * 10.0,
                order_hypo=observed_order(*[x["tip"]["r3"] for x in hypo[-3:]]),
                reference={r: ref[r] for r in MPM_H},
                reference_note="Richardson (p = 2, ratio 2) of the two finest hypo meshes, at "
                               "each MPM resolution's tip-column particle points")


def leg_cantilever_dyn(nys=(4, 8, 16, 32)) -> dict:
    cons = fe_cantilever_dyn(nys)
    lump = fe_cantilever_dyn((16,), lumped=True)
    half = fe_cantilever_dyn((16,), dt=0.005, nsteps=800)
    ref = {r: dict(first_peak=richardson(cons[-1]["tip"][r]["first_peak"],
                                         cons[-2]["tip"][r]["first_peak"]),
                   first_period=richardson(cons[-1]["tip"][r]["first_period"],
                                           cons[-2]["tip"][r]["first_period"])) for r in MPM_H}
    return dict(**_provenance("fe_cantilever_dyn", "leg_cantilever_dyn"),
                element="LadrunoBrick -geom hypo -formulation std, plane-strain slab, "
                        "ElasticIsotropic E = 1.2e8 nu = 0 rho = 1000, consistent mass; Newmark "
                        "0.5 0.25, dt 0.01, 400 steps, constant -b from rest",
                consistent=cons, lumped=lump, dt_half=half, reference=ref,
                linear=cantilever_linear(1.2e8),
                reference_note="Richardson (p = 2) of the two finest consistent-mass meshes")


def leg_footing() -> dict:
    rough = [fe_footing(n) for n in (4, 8, 16)]
    smooth = fe_footing(4, rough=False)
    # the locking probe: the same H8 WITHOUT B-bar (full integration), as a GIMP grid with
    # 2 x 2 particles per cell integrates -- mode A has no anti-locking until ADR §9 phase 4
    std = [fe_footing(n, form="std") for n in (4, 8, 16)]
    return dict(**_provenance("fe_footing", "leg_footing"),
                std_no_antilocking=std,
                element="LadrunoBrick -formulation bbar -geom linear (the ADR-79 control's "
                        "element), plane-strain slab, DruckerPrager (UW) mild cone",
                closed=footing_closed(), rough=rough, smooth_4=smooth,
                adr79_context=dict(
                    source="Ladruno_implementation/79_ladruno_geom_hypo_adr.md (nq20_nonassoc, "
                           "graded 30 m x 20 m box, smooth footing, H8 B-bar)",
                    q_over_exact=[1.0020, 0.9517, 0.9260], cells_across_b=[4, 8, 16]))


def leg_plate() -> dict:
    levels = [fe_plate(64, 30), fe_plate(128, 60), fe_plate(256, 120), fe_plate(512, 240)]
    long_traction = fe_plate(256, 120, top=16.0, load="traction")
    return dict(**_provenance("fe_plate", "leg_plate"),
                element="quad (4-node, PlaneStrain ElasticIsotropic E = 2e5 nu = 0.25), the "
                        "case's O-grid cells as elements; peak = quadratic Gauss-point fit "
                        "extrapolated to r = a on the ligament",
                levels=levels, long_traction=long_traction,
                howland_ktn=howland_ktn(0.5), heywood_ktn=heywood_ktn(0.5),
                seeding=plate_seed_checksums())


def main(argv=None) -> None:
    argv = list(sys.argv[1:] if argv is None else argv)
    what = argv[0] if argv else "closed"
    if what == "closed":
        print(json.dumps(dict(footing=footing_closed(), bar=bar_closed(),
                              bar_discrete=bar_discrete(), howland=howland_ktn(0.5),
                              heywood=heywood_ktn(0.5),
                              cantilever_qs=cantilever_linear(6e7),
                              cantilever_dyn=cantilever_linear(1.2e8)), indent=1))
    elif what == "elastica":
        print(json.dumps(elastica_tip(2.0), indent=1))
    elif what == "bar":
        print(json.dumps([bar_discrete(), bar_discrete(lumped=False)], indent=1))
    elif what == "seed":
        print(json.dumps(plate_seed_checksums(), indent=1))
    elif what == "contract1":
        write_plate_contract1(argv[1])
    elif what in ("fe_cantilever_qs", "fe_cantilever_dyn", "fe_footing", "fe_plate"):
        leg = {"fe_cantilever_qs": leg_cantilever_qs, "fe_cantilever_dyn": leg_cantilever_dyn,
               "fe_footing": leg_footing, "fe_plate": leg_plate}[what]
        print("wrote", _save(what, leg()))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
