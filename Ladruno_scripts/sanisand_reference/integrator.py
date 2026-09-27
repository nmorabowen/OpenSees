"""Stiff-ODE reference integration of the DM04 rate equations over one increment.

The rate equations are integrated in pseudo-time t in [0, 1] over a prescribed
(possibly mixed stress/strain controlled) increment with scipy's implicit Radau
IIA (or BDF) at a tight tolerance.  The right-hand side is discontinuous at
loading/unloading, at an alpha_in reseat and at the Macaulay kinks; each of those
is an EVENT that ends the current smooth segment, after which the mode is decided
again and integration restarts.  So every segment the solver sees is smooth.

Modes
  elastic : dsigma = E:deps, dalpha = dz = 0.  Ends when f reaches +ftol.
  plastic : DM04 loading index  L = N/H  with
              N = Q:E:deps = 2G n:de - K (n:r) deps_v            (numerator)
              H = Kp + Q:E:R = Kp + 2G n:R' - K D (n:r)          (denominator)
            dsigma = E:(deps - L R), dalpha = L (2/3) h b,
            dz = -cz <-L D> (z_max n + z),  de = -(1+e) deps_v.
            The consistency condition df = 0 holds identically on this branch
            (it is how L is derived), so f stays at its entry value up to the
            ODE tolerance; f at exit is reported.
            Ends when N -> 0 (unloading), when H changes sign, when
            (alpha - alpha_in):n -> 0 (alpha_in reseat), at a kink, or at p_floor.

Mode decision at a point ON the yield surface (|f| <= ftol):
  paper rule: if (alpha - alpha_in):n < 0 then alpha_in := alpha (DM04: a new
  loading process starts); then plastic iff N > 0 and H > 0; elastic iff N <= 0.
  N > 0 with H <= 0 has NO admissible rate solution in the continuous theory
  (the elastic branch violates f <= 0, the plastic branch has L < 0): the
  integrator STOPS with status "H_nonpositive" -- it never maps that case to
  "elastic".  (Mechanism F of WP-128/129 is exactly that mapping.)
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .model import (Options, Params, State, SQ23, I3, bounding_report,
                    elastic_mandel, elastic_moduli, m2t, mac, norm, plastic_weights,
                    quantities, strain_v2t, t2m, t2v, v2t, yield_f, ddot)

_ONE = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
_SQ2 = math.sqrt(2.0)


@dataclass
class Control:
    """Prescribed increment over the pseudo-time interval [0, 1].

    mask[i] True  -> component i is STRAIN controlled, value[i] = strain
                     increment (Voigt, engineering shear, compression positive);
    mask[i] False -> component i is STRESS controlled, value[i] = stress
                     increment (tensor component, kPa, compression positive).
    Order xx yy zz xy yz zx."""
    mask: tuple
    value: tuple

    @classmethod
    def strain(cls, deps):
        return cls(tuple([True] * 6), tuple(float(x) for x in deps))

    @property
    def pure_strain(self):
        return all(self.mask)

    def mandel(self):
        """(strain-known Mandel values, stress-known Mandel values)."""
        em = np.zeros(6)
        sm = np.zeros(6)
        for i in range(6):
            v = self.value[i]
            if i < 3:
                em[i] = v
                sm[i] = v
            else:
                em[i] = v / _SQ2      # engineering gamma -> sqrt2 * eps_ij
                sm[i] = v * _SQ2      # tensor sigma_ij -> sqrt2 * sigma_ij
        return em, sm


@dataclass
class Result:
    status: str                 # "ok" or the reason integration stopped early
    t_end: float                # pseudo-time reached (1.0 on success)
    state: State                # state at t_end
    start: dict                 # bounding report at the start
    end: dict                   # bounding report at t_end
    f_end: float
    max_rho_b: float            # max ||alpha|| / ||alpha^b_theta(n)|| along the path
    max_rho_alpha: float        # max ||alpha|| / bounding radius in alpha's own direction
    max_abs_f_plastic: float    # consistency drift: max |f| over plastic segments
    min_H_sign_margin: float    # min Hs/scale over plastic segments (<0 never taken)
    segments: list = field(default_factory=list)
    reseats: list = field(default_factory=list)
    notes: list = field(default_factory=list)
    nfev: int = 0
    path: dict = field(default_factory=dict)
    deps_total: list = field(default_factory=list)   # the realised strain increment
    uw_negative_h: bool = False   # plastic with (alpha-alpha_in):n < 0 (h < 0) anywhere
    min_a_plastic: float = float("inf")   # min (alpha-alpha_in):n over plastic samples

    def summary(self):
        s = self.state
        return dict(status=self.status, t_end=self.t_end, f_end=self.f_end,
                    max_rho_b=self.max_rho_b, rho_b_end=self.end.get("rho_b"),
                    max_rho_alpha=self.max_rho_alpha,
                    rho_alpha_end=self.end.get("rho_alpha"),
                    p_end=self.end.get("p"), eta_end=self.end.get("eta"),
                    e_end=s.e, n_segments=len(self.segments),
                    n_reseats=len(self.reseats), nfev=self.nfev,
                    max_abs_f_plastic=self.max_abs_f_plastic,
                    sigma=t2v(s.sigma).tolist(), alpha=t2v(s.alpha).tolist(),
                    z=t2v(s.z).tolist(), alpha_in=t2v(s.alpha_in).tolist(),
                    deps=self.deps_total, notes=self.notes,
                    uw_negative_h=self.uw_negative_h)


def _pack(sig, alpha, z, e, lcum, eps):
    return np.concatenate([t2v(sig), t2v(alpha), t2v(z), [e, lcum],
                           eps])


def _unpack(y):
    return v2t(y[0:6]), v2t(y[6:12]), v2t(y[12:18]), float(y[18]), float(y[19])


class _Increment:
    def __init__(self, state, control, P, O, rtol, atol_scale, method):
        self.P, self.O = P, O
        self.control = control
        self.em, self.sm = control.mandel()
        self.mask = np.array(control.mask, dtype=bool)
        self.S = np.where(self.mask)[0]
        self.T = np.where(~self.mask)[0]
        self.alpha_in = state.alpha_in.copy()
        self.rtol = rtol
        self.method = method
        self.frozen = None
        if O.elastic_moduli in ("frozen", "frozen_increment"):
            self.frozen = elastic_moduli(float(np.trace(state.sigma)) / 3.0,
                                         state.e, P, O)
        p0 = max(abs(float(np.trace(state.sigma)) / 3.0), 1.0e-3)
        emax = max(float(np.max(np.abs(self.em))) if len(self.S) else 0.0, 1.0e-12)
        self.emax = emax
        atol = np.empty(26)
        atol[0:6] = atol_scale * rtol * p0
        atol[6:12] = atol_scale * rtol
        atol[12:18] = atol_scale * rtol * max(1.0, P.z_max)
        atol[18] = atol_scale * rtol
        atol[19] = atol_scale * rtol * emax
        atol[20:26] = atol_scale * rtol * emax
        self.atol = atol
        typ = np.empty(26)
        typ[0:6] = p0
        typ[6:12] = 0.1
        typ[12:18] = 1.0
        typ[18] = 1.0
        typ[19:26] = emax
        self.jac_typ = typ
        self.nfev = 0
        self.uw_negative_h = False

    # ------------------------------------------------------------------
    def q_of(self, y, mode):
        sig, alpha, z, e, _ = _unpack(y)
        moduli = None
        if self.frozen is not None and (mode == "elastic"
                                        or self.O.elastic_moduli == "frozen_increment"):
            moduli = self.frozen
        return quantities(sig, alpha, z, e, self.alpha_in, self.P, self.O, moduli)

    def _solve_mixed(self, Emat):
        if len(self.T) == 0:
            return self.em.copy()
        d = self.em.copy()
        d[self.T] = 0.0
        rhs = self.sm[self.T] - Emat[np.ix_(self.T, self.S)] @ self.em[self.S]
        d[self.T] = np.linalg.solve(Emat[np.ix_(self.T, self.T)], rhs)
        return d

    def rates(self, y, mode, q=None):
        """(dy/dt, info) for the given mode."""
        self.nfev += 1
        if q is None:
            q = self.q_of(y, mode)
        P, O = self.P, self.O
        Ee = elastic_mandel(q.G, q.K)
        Qm = t2m(q.n - (q.nr / 3.0) * I3)
        EQ = Ee @ Qm
        info = dict(q=q)
        if mode == "elastic":
            dm = self._solve_mixed(Ee)
            N = float(EQ @ dm)
            dsm = Ee @ dm
            dalpha = np.zeros((3, 3))
            dz = np.zeros((3, 3))
            L = 0.0
            info.update(N=N, sg=0, w=0.0)
        else:
            w, hw, sg = plastic_weights(q, O)
            Rm = t2m(q.R)
            ER = Ee @ Rm
            Emat = Ee - w * np.outer(ER, EQ)
            dm = self._solve_mixed(Emat)
            N = float(EQ @ dm)
            L = w * N
            hL = hw * N
            dsm = Ee @ (dm - L * Rm)
            dalpha = hL * (2.0 / 3.0) * q.b
            dz = -P.cz * mac(-L * q.D) * (P.z_max * q.n + v2t(y[12:18]))
            info.update(N=N, sg=sg, w=w, hw=hw, L=L)
        deps = m2t(dm)
        dv = float(np.trace(deps))
        e = float(y[18])
        de = -(1.0 + (P.e_init if O.void_ratio_law == "initial" else e)) * dv
        dy = np.concatenate([t2v(m2t(dsm)), t2v(dalpha), t2v(dz), [de, L],
                             [deps[0, 0], deps[1, 1], deps[2, 2], 2 * deps[0, 1],
                              2 * deps[1, 2], 2 * deps[0, 2]]])
        return dy, info

    def jacobian(self, y, mode):
        """Forward-difference Jacobian with a FIXED relative step per column.

        scipy's own num_jac grows its step without bound on columns the
        right-hand side does not depend on (eps, the plastic-multiplier sum, z
        while the fabric is inactive) and overflows; this one does not.  Only
        the first 19 columns (sigma, alpha, z, e) can be non-zero: the rates do
        not depend on the accumulated strain or on sum(L)."""
        f0, _ = self.rates(y, mode)
        J = np.zeros((26, 26))
        typ = self.jac_typ
        for j in range(19):
            h = 1.0e-7 * max(abs(y[j]), typ[j])
            yp = y.copy()
            yp[j] += h
            fp, _ = self.rates(yp, mode)
            J[:, j] = (fp - f0) / h
        return J

    # ------------------------------------------------------------------
    def ftol(self, q):
        return self.O.ftol_rel * SQ23 * self.P.m * max(q.p, 0.0) + 1.0e-300

    def ntol(self, q):
        return 1.0e-12 * (2.0 * q.G + 3.0 * q.K) * self.emax * (1.0 + abs(q.nr))

    def decide(self, t, y, notes, reseats):
        """Return mode ('elastic'|'plastic') or ('stop', status)."""
        O = self.O
        q = self.q_of(y, "plastic")
        if not (q.p > O.p_floor):
            return ("stop", "p_floor")
        ft = self.ftol(q)
        if q.f < -ft:
            return "elastic"
        if q.f > ft:
            if t == 0.0:
                if q.f > max(ft, O.ftol_abs_start):
                    if O.start_outside == "stop":
                        return ("stop", "start_outside_yield")
                    notes.append(f"start outside the yield surface, f0 = {q.f:.3e}; "
                                 "integrated as on-surface (start_outside='plastic')")
            else:
                notes.append(f"t={t:.6g}: f = {q.f:.3e} > ftol {ft:.1e} at a segment "
                             "boundary (drift); treated as on-surface")
        # on the surface: DM04 alpha_in reseat at the onset of a loading process
        # ((alpha - alpha_in):n < 0).  Under the UW rule the only reseat is the
        # trial-direction test at t = 0 (in integrate()); a plastic onset with
        # (alpha - alpha_in):n < 0 inside the increment then runs with h < 0,
        # exactly what the C++ does (WP-128 mechanism G) -- noted, not repaired.
        if q.a < 0.0:
            if O.alpha_in_rule == "paper":
                self.alpha_in = v2t(y[6:12]).copy()
                reseats.append(dict(t=t, rule="paper",
                                    where="start" if t == 0.0 else "onset"))
                q = self.q_of(y, "plastic")
            else:
                self.uw_negative_h = True
                notes.append(f"t={t:.6g}: plastic onset with (alpha - alpha_in):n = "
                             f"{q.a:.3e} < 0 and no reseat (UW rule): h < 0")
        _, ip = self.rates(y, "plastic", q)
        nt = self.ntol(q)
        if ip["N"] > nt and ip["sg"] > 0:
            return "plastic"
        _, ie = self.rates(y, "elastic")
        if ie["N"] <= nt:
            return "elastic"
        if ip["N"] > nt and ip["sg"] <= 0:
            return ("stop", "H_nonpositive")
        return ("stop", "no_consistent_mode")

    # ------------------------------------------------------------------
    def events(self, mode, y0):
        P, O = self.P, self.O
        evs, names = [], []

        def add(fun, name, direction, scale):
            g0 = fun(0.0, y0)
            if abs(g0) <= 1.0e-13 * scale and direction == 0:
                return                          # starts on the switch: skip this segment
            fun.terminal = True
            fun.direction = direction
            evs.append(fun)
            names.append(name)

        def e_pfloor(t, y):
            return float(np.trace(v2t(y[0:6]))) / 3.0 + O.p_residual - O.p_floor
        add(e_pfloor, "p_floor", -1, 1.0)

        if mode == "elastic":
            def e_f(t, y):
                q = self.q_of(y, "elastic")
                return q.f - self.ftol(q)
            add(e_f, "yield", 1, 1.0)
        else:
            def e_N(t, y):
                _, info = self.rates(y, "plastic")
                return info["N"]
            add(e_N, "unload", -1, 1.0)

            def e_H(t, y):
                return self.q_of(y, "plastic").Hs
            add(e_H, "H_sign", 0, 1.0)

            def e_a(t, y):
                return self.q_of(y, "plastic").a
            add(e_a, "alpha_in_reversal", -1 if O.alpha_in_rule == "paper" else 0, 1.0)
            if O.kink_events:
                def e_zn(t, y):
                    q = self.q_of(y, "plastic")
                    return ddot(v2t(y[12:18]), q.n)
                add(e_zn, "kink_zn", 0, 1.0)

                def e_D(t, y):
                    return self.q_of(y, "plastic").D
                add(e_D, "kink_D", 0, 1.0)
        return evs, names


def integrate(state, control, P, O=None, rtol=1.0e-10, atol_scale=1.0e-2,
              method="Radau", max_segments=5000, record=True):
    """Integrate one prescribed increment from `state`.  Returns a Result."""
    from scipy.integrate import solve_ivp

    if O is None:
        O = Options()
    if not isinstance(control, Control):
        control = Control.strain(control)
    inc = _Increment(state, control, P, O, rtol, atol_scale, method)
    notes, reseats, segments = [], [], []

    # U6 (UW): the reversal test once per increment, on the ELASTIC trial
    # direction (alpha_n - alpha_in_n):(Ce:deps) < 0, before anything else.
    if O.alpha_in_rule == "uw":
        G0_, K0_ = inc.frozen if inc.frozen is not None else elastic_moduli(
            float(np.trace(state.sigma)) / 3.0, state.e, P, O)
        Ee = elastic_mandel(G0_, K0_)
        dm = inc._solve_mixed(Ee)
        trial = m2t(Ee @ dm)
        if ddot(state.alpha - state.alpha_in, trial) < 0.0:
            inc.alpha_in = state.alpha.copy()
            reseats.append(dict(t=0.0, rule="uw", where="start"))

    y = _pack(state.sigma, state.alpha, state.z, state.e, 0.0, np.zeros(6))
    start = bounding_report(state.sigma, state.alpha, state.z, state.e,
                            inc.alpha_in, P, O)
    t = 0.0
    status = "ok"
    max_rho = start["rho_b"]
    max_rho_a = start["rho_alpha"]
    max_abs_f_pl = 0.0
    min_hmargin = float("inf")
    min_a_pl = float("inf")
    tt, yy = [0.0], [y.copy()]
    stall = 0
    while t < 1.0:
        if len(segments) >= max_segments:
            status = "max_segments"
            break
        mode = inc.decide(t, y, notes, reseats)
        if isinstance(mode, tuple):
            status = mode[1]
            break
        evs, names = inc.events(mode, y)
        fun = (lambda tt_, yy_, m=mode: inc.rates(yy_, m)[0])
        jac = (lambda tt_, yy_, m=mode: inc.jacobian(yy_, m))
        sol = solve_ivp(fun, (t, 1.0), y, method=method, rtol=rtol, atol=inc.atol,
                        events=evs, jac=jac)
        if sol.status == -1:
            # keep the last accepted point: where the solver gave up is the finding
            t, y = float(sol.t[-1]), sol.y[:, -1].copy()
            pq = float(np.trace(v2t(y[0:6]))) / 3.0
            status = "solver_failed"
            qf = inc.q_of(y, mode)
            notes.append(f"solver failed at t={t:.6g}, p={pq:.4g} kPa, mode {mode}, "
                         f"(alpha-alpha_in):n={qf.a:.3e}, b:n={qf.bn:.3e}, "
                         f"Hs={qf.Hs:.3e}, rho_b={qf.rho_b:.3f}: {sol.message}")
            break
        fired = None
        if sol.status == 1:
            for k, te in enumerate(sol.t_events):
                if len(te) and abs(te[0] - sol.t[-1]) <= 1e-15 * max(1.0, sol.t[-1]):
                    fired = names[k]
                    break
            if fired is None:
                for k, te in enumerate(sol.t_events):
                    if len(te):
                        fired = names[k]
                        break
        t_new = float(sol.t[-1])
        y_new = sol.y[:, -1].copy()
        segments.append(dict(mode=mode, t0=t, t1=t_new, event=fired,
                             steps=int(len(sol.t) - 1)))
        if record:
            tt.extend(sol.t[1:].tolist())
            yy.extend([sol.y[:, k].copy() for k in range(1, sol.y.shape[1])])
        for k in range(1, sol.y.shape[1]):
            q = inc.q_of(sol.y[:, k], mode)
            if q.rho_b > max_rho:
                max_rho = q.rho_b
            max_rho_a = max(max_rho_a, q.rho_alpha)
            if mode == "plastic":
                max_abs_f_pl = max(max_abs_f_pl, abs(q.f))
                min_a_pl = min(min_a_pl, q.a)
                sc = max(abs(q.Hs), 1e-300)
                min_hmargin = min(min_hmargin, q.Hs / (abs(q.Hs) + abs(q.a * q.X) + 1e-300))
        stall = stall + 1 if (t_new - t) <= 1e-13 else 0
        if stall > 50:
            status = "chatter"
            t, y = t_new, y_new
            break
        t, y = t_new, y_new
        if fired == "p_floor":
            status = "p_floor"
            break
        if fired == "alpha_in_reversal" and O.alpha_in_rule == "paper":
            inc.alpha_in = v2t(y[6:12]).copy()
            reseats.append(dict(t=t, rule="paper", where="in-plastic"))
        if fired is None and t < 1.0:
            status = "solver_stopped"
            break

    sig, alpha, z, e, lcum = _unpack(y)
    end_state = State(sig, alpha, z, e, inc.alpha_in.copy())
    end = bounding_report(sig, alpha, z, e, inc.alpha_in, P, O)
    path = {}
    if record:
        Y = np.array(yy)
        path = dict(t=np.array(tt), y=Y)
    res = Result(status=status, t_end=t, state=end_state, start=start, end=end,
                 f_end=end["f"], max_rho_b=max(max_rho, end["rho_b"]),
                 max_rho_alpha=max(max_rho_a, end["rho_alpha"]),
                 max_abs_f_plastic=max_abs_f_pl,
                 min_H_sign_margin=min_hmargin, segments=segments,
                 reseats=reseats, notes=notes, nfev=inc.nfev, path=path,
                 deps_total=y[20:26].tolist(), min_a_plastic=min_a_pl,
                 uw_negative_h=inc.uw_negative_h or min_a_pl < -1.0e-12)
    return res


def path_table(res, P, O=None, every=1):
    """Per recorded point: t, p, q, eta, e, psi, rho_b, f, eps (6)."""
    if O is None:
        O = Options()
    out = []
    Y = res.path["y"]
    T = res.path["t"]
    for k in range(0, len(T), every):
        sig, alpha, z, e, lcum = _unpack(Y[k])
        p = float(np.trace(sig)) / 3.0
        s = sig - p * I3
        qd = math.sqrt(1.5) * norm(s)
        br = bounding_report(sig, alpha, z, e, res.state.alpha_in, P, O)
        out.append(dict(t=float(T[k]), p=p, q=qd, eta=qd / p if p > 0 else float("nan"),
                        e=e, psi=br["psi"], rho_b=br["rho_b"], f=br["f"],
                        eps=Y[k][20:26].tolist(), sigma=t2v(sig).tolist(),
                        lcum=lcum))
    return out
