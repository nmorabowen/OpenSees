"""Numpy oracle for the explicit MPM engine (WP-232b).

An independent reference implementation of the step frozen in tests/mpm_interface.md (uGIMP or
linear basis, P2G, mass cutoff, FLAC local damping, MUSL with a PIC/FLIP blend, grid faces and
the Coulomb wall, the critical time step, the incremental hypo-kernel entry). Written from the
interface and the literature cited there, never from SRC/. The fork's C++ worker (232c) must not
read this file: it codes against the interface and the case files only.

Running this file as a script fills every `pending-oracle` value in tests/mpm_cases.txt:

    python3.12 tests/mpm_reference.py            # fill pending values
    python3.12 tests/mpm_reference.py --report   # print every case's observables
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from hypo_incremental_reference import (  # noqa: E402
    hypo_increment_from_F,
    hypo_increment_from_F_plane_strain,
    push_stress,
    voigt3_to_tensor,
    voigt_to_tensor_stress,
)

CASES_PATH = os.path.join(HERE, "mpm_cases.txt")
PENDING = "pending-oracle"
FACES = ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")
AXIS = {"x": 0, "y": 1, "z": 2}


# --------------------------------------------------------------------------- #
#  Basis (interface section 4)                                                #
# --------------------------------------------------------------------------- #

def gimp_1d(r, h, lp):
    """uGIMP 1D factor S(r) and dS/dr, r = x_p - x_i (Steffen, Kirby & Berzins 2008)."""
    r = np.asarray(r, dtype=float)
    a = np.abs(r)
    s = np.sign(r)
    lp = np.broadcast_to(lp, r.shape)
    S = np.zeros_like(r)
    dS = np.zeros_like(r)
    m1 = a < lp
    m2 = (~m1) & (a < h - lp)
    m3 = (~m1) & (~m2) & (a < h + lp)
    S[m1] = 1.0 - (r[m1] ** 2 + lp[m1] ** 2) / (2.0 * h * lp[m1])
    dS[m1] = -r[m1] / (h * lp[m1])
    S[m2] = 1.0 - a[m2] / h
    dS[m2] = -s[m2] / h
    t = h + lp[m3] - a[m3]
    S[m3] = t * t / (4.0 * h * lp[m3])
    dS[m3] = -s[m3] * t / (2.0 * h * lp[m3])
    return S, dS


GHOST = 1   # ghost node layers beyond each face (interface section 2, Amendment 5)


class Grid:
    """Real grid [x0, x0 + n h] plus GHOST node layers beyond every face.

    Node indices run from -GHOST to n_d + GHOST per axis; the id of node (i, j, k) is
    (i+G) + (n_x+1+2G)*((j+G) + (n_y+1+2G)*(k+G)), x fastest (interface section 2).
    """

    def __init__(self, origin, h, cells):
        self.x0 = np.asarray(origin, dtype=float)
        self.h = float(h)
        self.n = np.asarray(cells, dtype=int)
        self.ndm = len(self.n)
        self.nn_axis = self.n + 1 + 2 * GHOST
        self.nnodes = int(np.prod(self.nn_axis))
        ids = np.arange(self.nnodes)
        self.ijk = np.zeros((self.nnodes, self.ndm), dtype=int)   # signed node indices
        rem = ids.copy()
        for a in range(self.ndm):
            self.ijk[:, a] = rem % self.nn_axis[a] - GHOST
            rem //= self.nn_axis[a]

    def node_id(self, idx):
        nid, mult = 0, 1
        for a in range(self.ndm):
            nid = nid + (np.asarray(idx[a]) + GHOST) * mult
            mult *= self.nn_axis[a]
        return nid

    def face_nodes(self, face):
        """The face layer and its ghost layers, with the face's inward normal."""
        a = AXIS[face[0]]
        if a >= self.ndm:
            raise ValueError(f"face {face} not in {self.ndm}D")
        normal = np.zeros(self.ndm)
        if face.endswith("min"):
            sel = self.ijk[:, a] <= 0
            normal[a] = 1.0
        else:
            sel = self.ijk[:, a] >= self.n[a]
            normal[a] = -1.0
        return np.nonzero(sel)[0], normal

    def shape(self, x, lp, basis):
        """Support of every particle: node ids (N, K), S (N, K), grad S (N, K, ndm)."""
        N = x.shape[0]
        Ss, dSs, idx = [], [], []
        for a in range(self.ndm):
            if basis == "gimp":
                c = np.floor((x[:, a] - self.x0[a]) / self.h).astype(int)
                cand = c[:, None] + np.arange(-1, 3)[None, :]
                r = x[:, a, None] - (self.x0[a] + cand * self.h)
                S, dS = gimp_1d(r, self.h, np.broadcast_to(lp[:, a, None], r.shape).copy())
                out = (cand < -GHOST) | (cand > self.n[a] + GHOST)
                # only zero-weight candidates may fall outside the ghost layers
                assert not np.any(S[out] != 0.0) and not np.any(dS[out] != 0.0)
                S[out] = 0.0
                dS[out] = 0.0
                cand = np.clip(cand, -GHOST, self.n[a] + GHOST)
            elif basis == "linear":
                c = np.floor((x[:, a] - self.x0[a]) / self.h).astype(int)
                c = np.clip(c, 0, self.n[a] - 1)
                xi = (x[:, a] - (self.x0[a] + c * self.h)) / self.h
                cand = np.stack([c, c + 1], axis=1)
                S = np.stack([1.0 - xi, xi], axis=1)
                dS = np.stack([np.full(N, -1.0 / self.h), np.full(N, 1.0 / self.h)], axis=1)
            else:
                raise ValueError(basis)
            Ss.append(S)
            dSs.append(dS)
            idx.append(cand + GHOST)
        S = np.ones((N, 1))
        G = [np.ones((N, 1)) for _ in range(self.ndm)]
        nid = np.zeros((N, 1), dtype=int)
        mult = 1
        for a in range(self.ndm):
            S = (S[:, :, None] * Ss[a][:, None, :]).reshape(N, -1)
            for b in range(self.ndm):
                fac = dSs[a] if a == b else Ss[a]
                G[b] = (G[b][:, :, None] * fac[:, None, :]).reshape(N, -1)
            nid = (nid[:, :, None] + mult * idx[a][:, None, :]).reshape(N, -1)
            mult *= self.nn_axis[a]
        return nid, S, np.stack(G, axis=-1)

    def admissible(self, x, lp, basis):
        """The particle centre inside the real grid (and 0 < lp <= h/2 for GIMP)."""
        lo = self.x0
        hi = self.x0 + self.n * self.h
        ok = np.all((x >= lo) & (x <= hi), axis=1)
        if basis == "gimp":
            ok &= np.all((lp > 0) & (lp <= 0.5 * self.h), axis=1)
        return ok


# --------------------------------------------------------------------------- #
#  Material (interface section 9)                                             #
# --------------------------------------------------------------------------- #

def lame(E, nu):
    E = np.asarray(E, dtype=float)
    nu = np.asarray(nu, dtype=float)
    return E * nu / ((1 + nu) * (1 - 2 * nu)), E / (2 * (1 + nu))


def elastic_D(lam, mu, ndm):
    """Batched isotropic elasticity matrices, Voigt with engineering shear strain."""
    lam = np.atleast_1d(lam)
    mu = np.atleast_1d(mu)
    N = lam.shape[0]
    if ndm == 3:
        D = np.zeros((N, 6, 6))
        for i in range(3):
            for j in range(3):
                D[:, i, j] = lam
            D[:, i, i] = lam + 2 * mu
            D[:, 3 + i, 3 + i] = mu
    else:
        D = np.zeros((N, 3, 3))
        D[:, 0, 0] = D[:, 1, 1] = lam + 2 * mu
        D[:, 0, 1] = D[:, 1, 0] = lam
        D[:, 2, 2] = mu
    return D


# --------------------------------------------------------------------------- #
#  Engine                                                                     #
# --------------------------------------------------------------------------- #

def velocity_gradient(v_grid, nid, G):
    """L_p[a][b] = sum_k v_grid[nid_pk, a] * dS_pk/dx_b (interface section 6.6)."""
    return np.einsum("pka,pkb->pab", v_grid[nid], G)


class StepRefused(Exception):
    pass


class Engine:
    """Particle state of interface section 3 plus the engine parameters of section 11."""

    def __init__(self, grid: Grid, basis="gimp"):
        self.grid = grid
        self.ndm = grid.ndm
        self.nv = 6 if self.ndm == 3 else 3
        self.basis = basis
        self.flip = 1.0
        self.damping = 0.0
        self.masscut = 1e-10
        self.safety = 0.5
        self.gravity = np.zeros(self.ndm)
        self.series = ("none",)
        self.fixes = []    # (face, [axes])
        self.walls = []    # (face, mu)
        self.t = 0.0
        self.parts = []

    # -- particles ----------------------------------------------------------
    def add_particles(self, ids, body, mat, x, lp, V0, m, v, s0, E, nu):
        self.parts.append(dict(id=np.asarray(ids, dtype=np.int64), body=np.asarray(body, dtype=np.int32),
                               mat=np.asarray(mat, dtype=np.int32), x=np.asarray(x, float),
                               lp=np.asarray(lp, float), V0=np.asarray(V0, float), m=np.asarray(m, float),
                               v=np.asarray(v, float), s0=np.asarray(s0, float),
                               E=np.asarray(E, float), nu=np.asarray(nu, float)))

    def finalize(self):
        P = {k: np.concatenate([p[k] for p in self.parts]) for k in self.parts[0]}
        order = np.argsort(P["id"], kind="stable")
        P = {k: a[order] for k, a in P.items()}
        N = P["id"].shape[0]
        self.id, self.body, self.mat = P["id"], P["body"], P["mat"]
        self.x, self.lp, self.V0, self.m, self.v = P["x"], P["lp"], P["V0"], P["m"], P["v"]
        self.V = self.V0.copy()
        self.F = np.broadcast_to(np.eye(self.ndm), (N, self.ndm, self.ndm)).copy()
        self.R = self.F.copy()
        self.eps_hat = np.zeros((N, self.nv))
        self.sig0_hat = P["s0"]
        self.sig_hat = self.sig0_hat.copy()
        self.sig = self.sig0_hat.copy()
        lam, mu = lame(P["E"], P["nu"])
        self.lam, self.mu = lam, mu
        self.D = elastic_D(lam, mu, self.ndm)
        self.m_ref = float(np.max(self.m))
        self.X0 = self.x.copy()

    # -- section 6.1 / 6.2 ---------------------------------------------------
    def critical_dt(self):
        # M_p from the COMMITTED tangent; for this hypoelastic law the tangent is constant D
        M = self.lam + 2 * self.mu
        rho = self.m / self.V
        c = np.sqrt(M / rho)
        speed = np.linalg.norm(self.v, axis=1)
        return self.grid.h / np.max(c + speed)

    def load_factor(self, t):
        if self.series[0] == "none":
            return 1.0
        return min(t / self.series[1], 1.0)

    # -- grid conditions (section 7) ----------------------------------------
    def _apply_conditions(self, v, active, full):
        for face in FACES:
            for (wface, mu) in self.walls:
                if wface != face:
                    continue
                nodes, n = self.grid.face_nodes(face)
                nodes = nodes[active[nodes]]
                vv = v[nodes]
                vn = vv @ n
                into = vn < 0.0
                if not np.any(into):
                    continue
                vt = vv - vn[:, None] * n[None, :]
                if full:
                    tn = np.linalg.norm(vt, axis=1)
                    red = np.minimum(mu * np.abs(vn), tn)
                    scale = np.where(tn > 0.0, 1.0 - red / np.where(tn > 0.0, tn, 1.0), 0.0)
                    new = vt * scale[:, None]
                else:
                    new = vt
                vv[into] = new[into]
                v[nodes] = vv
        for (face, axes) in self.fixes:
            nodes, _ = self.grid.face_nodes(face)
            for a in axes:
                v[nodes, a] = 0.0
        return v

    # -- one step (section 6) -------------------------------------------------
    def step(self, dt):
        g = self.grid
        ndm = self.ndm
        # 6.1 admissibility and the critical time step
        ok = g.admissible(self.x, self.lp, self.basis)
        if not np.all(ok):
            raise StepRefused(f"particle {int(self.id[np.argmin(ok)])} leaves the grid")
        dtc = self.critical_dt()
        if dt > self.safety * dtc:
            raise StepRefused(f"dt {dt!r} > {self.safety} * dt_crit {dtc!r}")
        # 6.2 load factor at t_{n+1}
        lam_t = self.load_factor(self.t + dt)
        b = lam_t * self.gravity
        # 6.3 P2G
        nid, S, G = g.shape(self.x, self.lp, self.basis)
        flat = nid.ravel()
        nn = g.nnodes
        mS = S * self.m[:, None]
        m_i = np.bincount(flat, weights=mS.ravel(), minlength=nn)
        p_i = np.stack([np.bincount(flat, weights=(mS * self.v[:, a, None]).ravel(), minlength=nn)
                        for a in range(ndm)], axis=1)
        fext = m_i[:, None] * b[None, :]
        sigT = self._cauchy_tensor(self.sig)
        fp = -self.V[:, None, None] * np.einsum("pab,pkb->pka", sigT, G)
        fint = np.stack([np.bincount(flat, weights=fp[:, :, a].ravel(), minlength=nn)
                         for a in range(ndm)], axis=1)
        f = fint + fext
        # 6.4 active nodes
        active = m_i > self.masscut * self.m_ref
        m_safe = np.where(active, m_i, 1.0)
        # 6.5 damping, momentum update, conditions
        vn_i = np.where(active[:, None], p_i / m_safe[:, None], 0.0)
        f = f - self.damping * np.abs(f) * np.sign(vn_i)
        v1_i = np.where(active[:, None], vn_i + dt * f / m_safe[:, None], 0.0)
        v1_i = self._apply_conditions(v1_i, active, full=True)
        a_i = np.where(active[:, None], (v1_i - vn_i) / dt, 0.0)
        # 6.6 particle update (MUSL)
        Sa = np.einsum("pk,pka->pa", S, a_i[nid])
        Sv1 = np.einsum("pk,pka->pa", S, v1_i[nid])
        v_new = self.flip * (self.v + dt * Sa) + (1.0 - self.flip) * Sv1
        x_new = self.x + dt * Sv1
        ps = np.stack([np.bincount(flat, weights=(mS * v_new[:, a, None]).ravel(), minlength=nn)
                       for a in range(ndm)], axis=1)
        vs_i = np.where(active[:, None], ps / m_safe[:, None], 0.0)
        vs_i = self._apply_conditions(vs_i, active, full=False)
        L = velocity_gradient(vs_i, nid, G)
        F_new = (np.eye(ndm)[None] + dt * L) @ self.F
        if ndm == 3:
            status, deps, R1, J1 = hypo_increment_from_F(self.F, F_new)
        else:
            status, deps, R1, J1 = hypo_increment_from_F_plane_strain(self.F, F_new)
        if np.any(status != 0):
            raise StepRefused("kernel status -1: the step fails before any material commits")
        eps_new = self.eps_hat + deps
        sig_hat_new = self.sig0_hat + np.einsum("pij,pj->pi", self.D, eps_new)
        sig_new = push_stress(R1, sig_hat_new)
        # 6.7 commit (elastic never refuses); 6.8 trial values become committed
        self.x, self.v, self.F, self.R = x_new, v_new, F_new, R1
        self.V = J1 * self.V0
        self.eps_hat, self.sig_hat, self.sig = eps_new, sig_hat_new, sig_new
        self.t += dt
        self._last_active = active
        return active

    def _cauchy_tensor(self, s):
        return voigt_to_tensor_stress(s) if self.ndm == 3 else voigt3_to_tensor(s)

    # -- unit-level hooks -----------------------------------------------------
    def active_nodes(self):
        nid, S, _ = self.grid.shape(self.x, self.lp, self.basis)
        m_i = np.bincount(nid.ravel(), weights=(S * self.m[:, None]).ravel(), minlength=self.grid.nnodes)
        return np.nonzero(m_i > self.masscut * self.m_ref)[0]

    # -- observables (section 10) ---------------------------------------------
    def kinetic(self):
        return 0.5 * float(np.sum(self.m * np.sum(self.v ** 2, axis=1)))

    def strain_energy(self):
        return 0.5 * float(np.sum(self.V * np.einsum("pi,pij,pj->p", self.eps_hat, self.D, self.eps_hat)))

    def com(self, body=None):
        sel = slice(None) if body is None else (self.body == body)
        return np.sum(self.m[sel, None] * self.x[sel], axis=0) / np.sum(self.m[sel])

    def com_v(self, body=None):
        sel = slice(None) if body is None else (self.body == body)
        return np.sum(self.m[sel, None] * self.v[sel], axis=0) / np.sum(self.m[sel])


# --------------------------------------------------------------------------- #
#  Case file (format in the header of tests/mpm_cases.txt)                    #
# --------------------------------------------------------------------------- #

def parse_cases(path: str = CASES_PATH):
    """Ordered dict name -> list of token lists (key first). 'end' closes a case."""
    cases, cur, name = {}, None, None
    with open(path, encoding="ascii") as fh:
        for line in fh:
            tok = line.split()
            if not tok or tok[0].startswith("#"):
                continue
            if tok[0] == "case":
                cur, name = [], tok[1]
            elif tok[0] == "end":
                cases[name] = cur
                cur = None
            else:
                cur.append(tok)
    return cases


def _get(lines, key, default=None):
    for t in lines:
        if t[0] == key:
            return t[1:]
    return default


def _all(lines, key):
    return [t[1:] for t in lines if t[0] == key]


def _floats(tok):
    return [float(t) for t in tok]


def _lattice(grid: Grid, body_tok):
    ndm = grid.ndm
    kind = body_tok[1]
    i = body_tok.index("ppc")
    ppc = int(body_tok[i + 1])
    mat = int(body_tok[body_tok.index("mat") + 1])
    rho = float(body_tok[body_tok.index("rho") + 1])
    axes = [grid.x0[a] + (np.arange(ppc * grid.n[a]) + 0.5) * grid.h / ppc for a in range(ndm)]
    # a_x fastest: meshgrid 'ij' then ravel in Fortran order
    mesh = np.meshgrid(*axes, indexing="ij")
    pts = np.stack([m.ravel(order="F") for m in mesh], axis=1)
    if kind == "box":
        lo = np.array(_floats(body_tok[2:2 + ndm]))
        hi = np.array(_floats(body_tok[2 + ndm:2 + 2 * ndm]))
        keep = np.all((pts >= lo) & (pts <= hi), axis=1)
    elif kind in ("disc", "ball"):
        c = np.array(_floats(body_tok[2:2 + ndm]))
        R = float(body_tok[2 + ndm])
        keep = np.sum((pts - c) ** 2, axis=1) < R * R
    else:
        raise ValueError(kind)
    pts = pts[keep]
    V0 = (grid.h / ppc) ** ndm
    lp = grid.h / (2 * ppc)
    return pts, V0, lp, mat, rho


def build_engine(lines):
    ndm = int(_get(lines, "ndm")[0])
    grid = Grid(_floats(_get(lines, "origin")), float(_get(lines, "h")[0]),
                [int(t) for t in _get(lines, "cells")])
    eng = Engine(grid, _get(lines, "basis", ["gimp"])[0])
    nv = eng.nv
    mats = {}
    for t in _all(lines, "material"):
        assert t[1] == "elastic"
        mats[int(t[0])] = (float(t[2]), float(t[3]))
    vel = {int(t[0]): t[1:] for t in _all(lines, "velocity")}
    stress = {int(t[0]): t[1:] for t in _all(lines, "stress")}
    shifts = {int(t[0]): t[1:] for t in _all(lines, "shift")}
    next_id = 0
    for bt in _all(lines, "body"):
        bid = int(bt[0])
        pts, V0, lp, mat, rho = _lattice(grid, bt)
        if bid in shifts:
            pts = pts + np.array(_floats(shifts[bid][:ndm]))
        N = pts.shape[0]
        v = np.zeros((N, ndm))
        if bid in vel:
            vt = vel[bid]
            if vt[0] == "uniform":
                v[:] = _floats(vt[1:1 + ndm])
            elif vt[0] == "sine":
                comp, ax = AXIS[vt[1]], AXIS[vt[2]]
                v0, beta, xref = _floats(vt[3:6])
                v[:, comp] = v0 * np.sin(beta * (pts[:, ax] - xref))
            else:
                raise ValueError(vt)
        s0 = np.zeros((N, nv))
        if bid in stress:
            st = stress[bid]
            if st[0] == "geostatic":
                assert ndm == 2
                ytop, rho_s, gg, K0 = _floats(st[1:5])
                sv = -rho_s * gg * (ytop - pts[:, 1])
                s0[:, 0] = K0 * sv
                s0[:, 1] = sv
            elif st[0] == "uniform":
                s0[:] = _floats(st[1:1 + nv])
            else:
                raise ValueError(st)
        E, nu = mats[mat]
        eng.add_particles(np.arange(next_id, next_id + N), np.full(N, bid), np.full(N, mat), pts,
                          np.full((N, ndm), lp), np.full(N, V0), np.full(N, rho * V0), v, s0,
                          np.full(N, E), np.full(N, nu))
        next_id += N
    for pt in _all(lines, "particle"):
        pid, bid, mat = int(pt[0]), int(pt[1]), int(pt[2])
        rho, V0 = float(pt[3]), float(pt[4])
        i = pt.index("x")
        x = _floats(pt[i + 1:i + 1 + ndm])
        lp = _floats(pt[i + 2 + ndm:i + 2 + 2 * ndm])
        v = _floats(pt[i + 3 + 2 * ndm:i + 3 + 3 * ndm])
        s = _floats(pt[i + 4 + 3 * ndm:i + 4 + 3 * ndm + nv])
        E, nu = mats[mat]
        eng.add_particles([pid], [bid], [mat], [x], [lp], [V0], [rho * V0], [v], [s], [E], [nu])
    for t in _all(lines, "fix"):
        eng.fixes.append((t[0], [AXIS[c] for c in t[1]]))
    for t in _all(lines, "wall"):
        eng.walls.append((t[0], float(t[1])))
    if _get(lines, "gravity"):
        eng.gravity = np.array(_floats(_get(lines, "gravity")))
    ser = _get(lines, "series", ["none"])
    eng.series = ("none",) if ser[0] == "none" else ("linear", float(ser[1]))
    eng.damping = float(_get(lines, "damping", ["0"])[0])
    eng.flip = float(_get(lines, "flip", ["1"])[0])
    eng.masscut = float(_get(lines, "masscut", ["1e-10"])[0])
    eng.safety = float(_get(lines, "safety", ["0.5"])[0])
    eng.finalize()
    return eng


def _period(t, vc):
    cross = []
    for k in range(len(vc) - 1):
        if vc[k] < 0.0 <= vc[k + 1]:
            cross.append(t[k] + (t[k + 1] - t[k]) * (0.0 - vc[k]) / (vc[k + 1] - vc[k]))
    if len(cross) < 2:
        return float("nan")
    return float(np.mean(np.diff(cross)))


def run_case(lines) -> dict:
    """Run one case and return every observable it may ask for, keyed by observable name."""
    kind = _get(lines, "kind")[0]
    res = {}
    if kind == "shape1d":
        return res  # evaluated per expect line (see shape1d_values)
    if kind == "shape":
        ndm = int(_get(lines, "ndm")[0])
        grid = Grid(_floats(_get(lines, "origin")), float(_get(lines, "h")[0]),
                    [int(t) for t in _get(lines, "cells")])
        pt = _get(lines, "point")
        x = np.array([_floats(pt[:ndm])])
        lp = np.array([_floats(pt[ndm + 1:2 * ndm + 1])])
        nid, S, G = grid.shape(x, lp, _get(lines, "basis", ["gimp"])[0])
        nodes = {}
        for k in range(nid.shape[1]):
            if S[0, k] != 0.0 or np.any(G[0, k] != 0.0):
                nodes[int(nid[0, k])] = np.concatenate([[S[0, k]], G[0, k]])
        res["node"] = nodes
        return res
    eng = build_engine(lines)
    if kind == "critdt":
        res["dt_crit"] = np.array([eng.critical_dt()])
        return res
    if kind == "masscut":
        res["active_nodes"] = eng.active_nodes().astype(float)
        return res
    dt = float(_get(lines, "dt")[0])
    steps = int(_get(lines, "steps")[0])
    E0 = eng.kinetic() + eng.strain_energy()
    com0 = {int(b): eng.com(int(b)) for b in np.unique(eng.body)}
    t_hist, vcx, vmax = [0.0], [eng.com_v()[0]], [float(np.max(np.linalg.norm(eng.v, axis=1)))]
    res["refused"] = False
    try:
        for _ in range(steps):
            eng.step(dt)
            t_hist.append(eng.t)
            vcx.append(eng.com_v()[0])
            vmax.append(float(np.max(np.linalg.norm(eng.v, axis=1))))
    except StepRefused as exc:
        res["refused"] = True
        res["refused_at"] = len(t_hist) - 1
        res["refusal"] = str(exc)
        return res
    res["engine"] = eng
    res["period_comx"] = np.array([_period(t_hist, vcx)])
    E1 = eng.kinetic() + eng.strain_energy()
    res["energy_ratio_end"] = np.array([E1 / E0 if E0 > 0.0 else float("nan")])
    res["momentum"] = np.sum(eng.m[:, None] * eng.v, axis=0)
    res["max_speed"] = np.array([max(vmax)])
    res["com_disp_x"] = {b: np.array([eng.com(b)[0] - com0[b][0]]) for b in com0}
    res["com_v_x"] = {b: np.array([eng.com_v(b)[0]]) for b in com0}
    res["com_v"] = {b: eng.com_v(b) for b in com0}
    res["com_disp"] = {b: eng.com(b) - com0[b] for b in com0}
    res["F_error"] = np.array([np.max(np.abs(eng.F - np.eye(eng.ndm)))])
    res["max_abs_stress"] = np.array([np.max(np.abs(eng.sig))])
    res["fields"] = particle_fields(eng)
    return res


def particle_fields(eng):
    """id -> concatenated [x, v, F, V, sig, sig_hat, eps_hat] (field-line order)."""
    out = {}
    for p in range(eng.id.shape[0]):
        out[int(eng.id[p])] = np.concatenate([eng.x[p], eng.v[p], eng.F[p].ravel(), [eng.V[p]],
                                              eng.sig[p], eng.sig_hat[p], eng.eps_hat[p]])
    return out


def selfweight_residual(eng, rho0, g, l0, ybase):
    Y = eng.X0[:, 1]
    exact = -rho0 * abs(g) * (l0 - (Y - ybase))
    return (eng.sig[:, 1] - exact) / (rho0 * abs(g) * l0)


def observable(res, obs, args):
    """The value of observable `obs` with its arguments, as a float array."""
    if obs in ("period_comx", "energy_ratio_end", "momentum", "max_speed", "F_error",
               "max_abs_stress", "dt_crit", "active_nodes"):
        return res[obs]
    if obs in ("com_disp_x", "com_v_x", "com_v", "com_disp"):
        return res[obs][int(args[0])]
    if obs == "translation_error":
        eng = res["engine"]
        v0 = np.array(_floats(args))
        return np.array([np.max(np.linalg.norm(eng.x - (eng.X0 + v0 * eng.t), axis=1))])
    if obs in ("selfweight_error", "cell_crossing_meter"):
        r = selfweight_residual(res["engine"], *_floats(args))
        if obs == "selfweight_error":
            return np.array([np.max(np.abs(r))])
        return np.array([math.sqrt(float(np.mean(r ** 2)))])
    raise KeyError(obs)


def split_expect(tok):
    """'expect <obs> [args] <op> <values...> [tol a r]' -> (obs, args, op, values, tol)."""
    obs = tok[0]
    k = next(i for i, t in enumerate(tok) if t in ("eq", "le", "ge"))
    args, op = tok[1:k], tok[k]
    rest = tok[k + 1:]
    tol = None
    if "tol" in rest:
        j = rest.index("tol")
        tol = (float(rest[j + 1]), float(rest[j + 2]))
        rest = rest[:j]
    return obs, args, op, rest, tol


def shape1d_values(tok):
    """'s1 h lp r eq S dS' -> (S, dS) from the oracle."""
    h, lp, r = float(tok[1]), float(tok[2]), float(tok[3])
    S, dS = gimp_1d(np.array([r]), h, np.array([lp]))
    return np.array([S[0], dS[0]])


# --------------------------------------------------------------------------- #
#  Filling pending values                                                     #
# --------------------------------------------------------------------------- #

def _fmt(a):
    return " ".join(f"{v:.17e}" for v in np.asarray(a, dtype=float).ravel())


def field_line(pid, vec, ndm):
    nv = 6 if ndm == 3 else 3
    i = 0
    parts = [f"field {pid}"]
    for key, n in (("x", ndm), ("v", ndm), ("F", ndm * ndm), ("V", 1), ("sig", nv), ("sighat", nv),
                   ("epshat", nv)):
        parts.append(f"{key} {_fmt(vec[i:i + n])}")
        i += n
    return " ".join(parts)


def fill_pending(path: str = CASES_PATH, verbose=True) -> int:
    with open(path, encoding="ascii") as fh:
        text = fh.read().split("\n")
    cases = parse_cases(path)
    out, name, n, cache = [], None, 0, {}
    for line in text:
        tok = line.split()
        if tok and tok[0] == "case":
            name = tok[1]
        elif tok and tok[0] == "end":
            name = None
        if name and tok and tok[0] == "expect" and PENDING in tok:
            if name not in cache:
                if verbose:
                    print(f"  running {name} ...", flush=True)
                cache[name] = run_case(cases[name])
            res = cache[name]
            obs, args, op, vals, tol = split_expect(tok[1:])
            assert op == "eq", line
            if obs == "particle_fields":
                out.append(line.replace(PENDING, "fields"))
                ndm = res["engine"].ndm
                for pid, vec in sorted(res["fields"].items()):
                    out.append(field_line(pid, vec, ndm))
            else:
                out.append(line.replace(PENDING, _fmt(observable(res, obs, args))))
            n += 1
            continue
        out.append(line)
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("\n".join(out))
    return n


def report(path: str = CASES_PATH):
    for name, lines in parse_cases(path).items():
        kind = _get(lines, "kind")[0]
        if kind in ("shape1d", "shape"):
            continue
        res = run_case(lines)
        if res.get("refused"):
            print(f"{name}: refused at step {res['refused_at']}: {res['refusal']}")
            continue
        if kind != "run":
            print(name, {k: v for k, v in res.items()})
            continue
        line = [f"{name}:"]
        for t in _all(lines, "expect"):
            if t[0] in ("refused", "accepted", "particle_fields"):
                continue
            obs, args, op, vals, tol = split_expect(t)
            line.append(f"{obs}{args if args else ''}={observable(res, obs, args)}")
        print(" ".join(line))


if __name__ == "__main__":
    if "--report" in sys.argv:
        report()
    else:
        print(f"filled {fill_pending()} pending values in {CASES_PATH}")
