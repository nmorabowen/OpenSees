"""Reference (numpy) of the Model C MCFT strip and its fibre section, for Ladruno WP-254.

An independent vectorised implementation of the strip law that LadrunoMCFTStrip
(SRC/material/nD/LadrunoMCFTKernel.h) implements in C++, written from the published law: rotating
MCFT (Vecchio & Collins 1986) with one transverse unknown per strip, Popovics/Mander compression with
Coleman-Spacone e20, beta softening, Collins-Mitchell tension stiffening, the 'sum' crack check,
Palermo-Vecchio unloading, and the bar/tie steel law; calibrated in the Model C beam on Sezen & Moehle
(2002) and Lynn et al. (1996). Units N, mm, MPa; tension positive.

It evaluates every strip of a section at once, scans all 81 root candidates and bisects them as
arrays: a different code path from the C++ kernel (which scans lazily, strip by strip), with the same
floating-point operations, so the two agree to round-off.

Used by tests/test_ladrunoMCFTStrip_cpp.py (strip and section paths against the g++ gate) and
tests/test_ladrunoMCFTStrip.py (the OpenSees classes).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

ES = 200e3


@dataclass
class Steel:
    fy: float
    fu: float
    Es: float = ES
    esh: float = 0.008
    esu: float = 0.12

    def env(self, a):
        a = np.asarray(a, float)
        ey = self.fy / self.Es
        r = np.clip((a - self.esh) / (self.esu - self.esh), 0.0, 1.0)
        hard = self.fy + (self.fu - self.fy) * (1.0 - (1.0 - r) ** 2)
        out = np.where(a <= ey, self.Es * a, np.where(a <= self.esh, self.fy, hard))
        return np.where(a > self.esu, 1e-3 * self.fy, out)

    def stress(self, e, ep):
        e = np.asarray(e, float)
        ey = self.fy / self.Es
        trial = self.Es * (e - ep)
        up = np.where(e > self.esu, 1e-3 * self.fy, self.env(np.maximum(e, ey)))
        lo = -np.where(-e > self.esu, 1e-3 * self.fy, self.env(np.maximum(-e, ey)))
        return np.clip(trial, lo, up)


@dataclass
class Concrete:
    fc: float
    fp: np.ndarray        # per strip peak stress
    ep: np.ndarray        # per strip peak strain
    e20: np.ndarray       # per strip strain at 0.2 fp
    smx: float
    smy: float
    agg: float = 19.0
    k_beta: float = 170.0
    c_vci: float = 0.18

    def __post_init__(self):
        self.Ec = 4700.0 * math.sqrt(self.fc)
        self.fcr = 0.33 * math.sqrt(self.fc)
        self.ecr = self.fcr / self.Ec

    def f_env(self, c):
        fp, ep, e20, Ec = self.fp, self.ep, self.e20, self.Ec
        x = c / ep
        r = Ec / (Ec - fp / ep)
        pop = fp * x * r / (r - 1.0 + x ** r)
        post = np.maximum(fp - 0.8 * fp * (c - ep) / (e20 - ep), 0.2 * fp)
        return np.where(c <= ep, pop, post)

    def ft(self, x):
        return np.where(x <= self.ecr, self.Ec * x, self.fcr / (1.0 + np.sqrt(500.0 * np.maximum(x - self.ecr, 0.0))))

    def membrane(self, ex, ey, g, rho, fyt, fst, hist):
        c0 = 0.5 * (ex + ey)
        R = np.sqrt((0.5 * (ex - ey)) ** 2 + (0.5 * g) ** 2) + 1e-30
        e1, e2 = c0 + R, c0 - R
        cos2, sin2 = (ex - ey) / (2 * R), g / (2 * R)
        beta = np.minimum(1.0, 1.0 / (0.8 + self.k_beta * np.maximum(e1, 0.0)))
        cc = np.maximum(-e2, 0.0)
        fcomp = self.f_env(cc)
        f1 = self.ft(e1)
        if hist is not None:
            em1, em2 = hist
            xx = em2 / self.ep
            p2 = np.minimum(self.ep * (0.166 * xx ** 2 + 0.132 * xx), 0.9 * em2)
            reld = self.f_env(em2) * np.clip((cc - p2) / np.maximum(em2 - p2, 1e-12), 0.0, 1.0)
            fcomp = np.where(cc < em2, reld, fcomp)
            sec = self.ft(em1) * e1 / np.maximum(em1, 1e-12)
            f1 = np.where((e1 > 0) & (e1 < em1), sec, f1)
        phi = 0.5 * np.arctan2(sin2, cos2)
        th = np.abs(0.5 * np.pi - np.abs(phi))
        s_th = 1.0 / (np.abs(np.sin(th)) / self.smx + np.abs(np.cos(th)) / self.smy + 1e-12)
        w = np.maximum(e1, 0.0) * s_th
        vci = self.c_vci * math.sqrt(self.fc) / (0.31 + 24.0 * w / (self.agg + 16.0))
        thc = np.clip(th, 1e-3, 0.5 * np.pi - 1e-3)
        cap = vci * np.tan(thc) + rho * np.maximum(fyt - fst, 0.0)
        f1 = np.where(e1 > 0, np.minimum(f1, cap), f1)
        f2 = -beta * fcomp
        f2 = np.where(e2 > 0, 0.0, f2)
        f1 = np.where(e1 < 0, -self.f_env(np.maximum(-e1, 0.0)), f1)
        sx = 0.5 * (f1 + f2) + 0.5 * (f1 - f2) * cos2
        sy = 0.5 * (f1 + f2) - 0.5 * (f1 - f2) * cos2
        txy = 0.5 * (f1 - f2) * sin2
        return sx, sy, txy, e1, e2, w


@dataclass
class State:
    ey: np.ndarray            # committed eps_y per strip (root reference)
    em1: np.ndarray
    em2: np.ndarray
    ept: np.ndarray           # tie plastic strain per strip
    epb: np.ndarray           # bar plastic strain per bar
    use_hist: bool = False


@dataclass
class Strips:
    """n independent strips, each driven at its own (eps_x, gamma_xy); the section builds on it."""
    conc: Concrete
    rho: np.ndarray
    tie_steel: Steel

    def new_state(self, nbar=0) -> State:
        n = len(self.rho)
        return State(np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(nbar))

    def resid(self, ex, ey, g, st):
        fst = self.tie_steel.stress(ey, st.ept)
        hist = (st.em1, st.em2) if st.use_hist else None
        sx, sy, txy, e1, e2, w = self.conc.membrane(ex, ey, g, self.rho, self.tie_steel.fy, fst, hist)
        return sy + self.rho * fst, (sx, txy, e1, e2, w, fst)

    def solve_ey(self, ex, g, st):
        if st.use_hist:
            off = np.geomspace(1e-7, 0.15, 40)
            cand = np.concatenate([-off[::-1], [0.0], off])
            Y = np.clip(st.ey[None, :] + cand[:, None], -0.01, 0.15)
            F, _ = self.resid(ex[None, :], Y, g[None, :], st)
            inc = (F[:-1] <= 0) & (F[1:] >= 0)
            mid = 0.5 * (Y[:-1] + Y[1:])
            dist = np.where(inc, np.abs(mid - st.ey[None, :]), np.inf)
            k = np.argmin(dist, axis=0)
            ok = np.isfinite(dist[k, np.arange(len(k))])
            lo = np.where(ok, Y[k, np.arange(len(k))], -0.004)
            hi = np.where(ok, Y[k + 1, np.arange(len(k))], 0.15)
        else:
            lo = np.full_like(ex, -0.004); hi = np.full_like(ex, 0.15)
        for _ in range(60):
            m = 0.5 * (lo + hi)
            f, _ = self.resid(ex, m, g, st)
            lo = np.where(f <= 0, m, lo); hi = np.where(f <= 0, hi, m)
        return 0.5 * (lo + hi)

    def stresses(self, ex, g, st):
        """(sx, txy, ey, e1, e2, w, fst) per strip at strains (ex, g), committed state st."""
        ex = np.asarray(ex, float); g = np.asarray(g, float)
        ey = self.solve_ey(ex, g, st)
        _, (sx, txy, e1, e2, w, fst) = self.resid(ex, ey, g, st)
        return dict(sx=sx, txy=txy, ey=ey, e1=e1, e2=e2, w=w, fst=fst)

    def commit(self, ex, g, st):
        o = self.stresses(ex, g, st)
        st.em1 = np.maximum(st.em1, o["e1"]) if st.use_hist else np.maximum(0.0, o["e1"])
        st.em2 = np.maximum(st.em2, np.maximum(-o["e2"], 0.0)) if st.use_hist else np.maximum(0.0, -o["e2"])
        st.ept = o["ey"] - self.tie_steel.stress(o["ey"], st.ept) / self.tie_steel.Es
        st.ey = o["ey"]
        st.use_hist = True
        return o


@dataclass
class Section:
    """Model C's section: n strips through h, bars at discrete y (area off the nearest strip), parabolic
    gamma profile, N = sum sx Ac + sum fs A, M = -sum (sx Ac y + fs A y), V = sum txy Ac (no dowel)."""
    b: float
    h: float
    conc: Concrete
    bars_y: np.ndarray
    bars_A: np.ndarray
    rho: np.ndarray
    bar_steel: Steel
    tie_steel: Steel
    nstrip: int = 40
    y: np.ndarray = field(init=False)
    Ac: np.ndarray = field(init=False)
    gy: np.ndarray = field(init=False)

    def __post_init__(self):
        dy = self.h / self.nstrip
        self.y = -self.h / 2 + dy * (np.arange(self.nstrip) + 0.5)
        Ac = np.full(self.nstrip, self.b * dy)
        for yb, A in zip(self.bars_y, self.bars_A):
            Ac[int(np.argmin(np.abs(self.y - yb)))] -= A
        self.Ac = Ac
        self.gy = 1.5 * (1 - (2 * self.y / self.h) ** 2)
        self.strips = Strips(self.conc, self.rho, self.tie_steel)

    def new_state(self) -> State:
        return self.strips.new_state(len(self.bars_y))

    def forces(self, e, st, out=False):
        e0, kap, gam = e
        ex = e0 - self.y * kap
        g = gam * self.gy
        o = self.strips.stresses(ex, g, st)
        eb = e0 - self.bars_y * kap
        fs = self.bar_steel.stress(eb, st.epb)
        sx, txy = o["sx"], o["txy"]
        N = np.sum(sx * self.Ac) + np.sum(fs * self.bars_A)
        M = -np.sum(sx * self.Ac * self.y) - np.sum(fs * self.bars_A * self.bars_y)
        V = np.sum(txy * self.Ac)
        s = np.array([N, M, V])
        if out:
            return s, dict(o, eb=eb, fs=fs, ex=ex, g=g)
        return s

    def commit(self, e, st):
        _, o = self.forces(e, st, out=True)
        st.epb = o["eb"] - self.bar_steel.stress(o["eb"], st.epb) / self.bar_steel.Es
        self.strips.commit(o["ex"], o["g"], st)
        return o


def mander_rect(fc, bc, dc, wi, s, dbt, Asy, Asz, fyh, rho_cc):
    """Mander et al. (1988) confined peak for a rectangular core."""
    s1 = max(s - dbt, 0.0)
    ke = max((1 - sum(w * w for w in wi) / (6 * bc * dc)) * (1 - s1 / (2 * bc)) * (1 - s1 / (2 * dc)) / (1 - rho_cc), 0.0)
    ry, rz = Asy / (s * dc), Asz / (s * bc)
    fl = math.sqrt(ke * ry * fyh * ke * rz * fyh)
    fcc = fc * (-1.254 + 2.254 * math.sqrt(1 + 7.94 * fl / fc) - 2 * fl / fc)
    ecc = 0.002 * (1 + 5 * (fcc / fc - 1))
    return fcc, ecc


def column_section(fc=25.6, fyl=331.0, ful=496.0, fyt=400.0, fut=600.0, rho_t=0.001, s=457.2, d=381.0,
                   dbl=32.3, Ab=819.0, b=457.2, h=457.2, lp_h=0.25, gf_mult=0.5, dbt=9.5, nstrip=40):
    """A Lynn-type 18 in column section (3CLH18-like): 8 bars 3-2-3, ties in the core strips only, e20 from
    the crushing energy with Lip = lp_h h. Returns the Section plus the per-strip data (y, core, fp, ep, e20)
    a fibre model needs."""
    yb = d - h / 2
    cov = h / 2 - yb - dbl / 2 - dbt
    Av = rho_t * b * s
    bc = h - 2 * (cov + dbt / 2)
    fcc, ecc = mander_rect(fc, bc, bc, [yb - dbl] * 8, s, dbt, Av, Av, fyt, 8 * Ab / bc ** 2)
    y = -h / 2 + h / nstrip * (np.arange(nstrip) + 0.5)
    core = np.abs(y) <= h / 2 - cov
    fp = np.where(core, fcc, fc); ep = np.where(core, ecc, 0.002)
    Ec = 4700 * math.sqrt(fc)
    Gf = gf_mult * np.where(core, 1.7 * 2 * fc, 2 * fc)
    e20 = Gf / (0.6 * fp * lp_h * h) - 0.8 * fp / Ec + ep
    conc = Concrete(fc, fp, ep, e20, smx=yb, smy=s)
    sec = Section(b, h, conc, np.array([-yb, 0.0, yb]), np.array([3, 2, 3]) * Ab, np.where(core, rho_t, 0.0),
                  Steel(fyl, ful), Steel(fyt, fut), nstrip=nstrip)
    return sec, dict(y=y, core=core, fp=fp, ep=ep, e20=e20, yb=yb, s=s, fyt=fyt, fut=fut, rho_t=rho_t,
                     fyl=fyl, ful=ful, Ab=Ab, b=b, h=h, fc=fc)
