"""O1 continuum rate oracle: the rate equations (S.40)-(S.41) integrated by Radau.

One increment = a linear-in-pseudo-time ramp t in [0, 1] of the prescribed strain (or, for
the driver helpers, mixed strain / zero-stress-rate control).  There is NO return map and NO
projection: in the plastic mode the closed-form multiplier (S.41) keeps F = 0 up to the ODE
tolerance and the residual |F| is reported (flags['max_F_rel']).

Modes and switches (every segment the solver sees is smooth):
  elastic : eps^e' = eps', pi_i' = 0.  Ends at the yield crossing F = 0 (terminal event,
            root-found by solve_ivp; this is the drift control at elastic -> plastic).
  plastic : lambda' = <f:a^e:eps'> / (f:a^e:q + H), eps^e' = eps' - lambda' q,
            pi_i' = sqrt(2/3) h lambda' (pi_i* - pi_i) Omega.
            Ends at unloading (numerator -> 0), at a loss of uniqueness (denominator -> 0:
            the run STOPS with status 'den_nonpositive', never mapped to elastic), at the
            B > 0 guard of (S.23), at a cap switch (planar: the corner eta = c1 M; smooth:
            the clamp points eta_1, eta_2), or p -> 0.
  decide  : F < -ftol -> elastic; |F| <= ftol -> plastic iff the plastic numerator > 0
            (and the denominator > 0), elastic iff the elastic numerator <= 0.
            F > ftol at a segment start -> status 'outside' (never projected).
H = 0 crossings are non-terminal events, recorded in flags['H_zero'] (K1.6 check).

With p_min > 0 (sheet 9.7) the floor modes are 'floor' and 'plastic_floor' (Koiter, S.55) and, since round 4,
'floor_rule': the rate form (S.60) of the split with the image-pressure rule (S.57), used on the floor and on F
wherever lam' > 0, lam_f' - lam_f^tr' > 0 and f:a^e:P < 0 hold (checked first at every floor decision); it is left
at the edges of that set (events rule_off / rule_off_fAP), and re-decided when P:a^e:eps' changes sign
(rule_branch).  floor_rule=False keeps (S.55) everywhere.

State vector y (19): eps^e (Mandel 6), pi_i, v, W = int sigma:deps, int |sigma:deps|,
Dp = int lambda' sigma:q, eps_p_v, eps_p_s, total strain (Mandel 6).

Specific volume (sheet 1.2 / S.40, G2 owner decision 2026-10-01): v is ALGEBRAIC from the total
strain, v = v_n exp(tr eps - tr eps_n) = v0 exp(tr eps) (the exponential update; v' = v tr eps').
Every constitutive call uses that algebraic value (`_Increment.vof`); y[7] still carries the ODE
v' = v tr eps' as a diagnostic only, and the committed State.v is the algebraic value, so the
identity v = v0 exp(tr eps) holds to round-off (sheet 13.10).  It is the same rule for
kin='small' and kin='log' (under log strain tr eps = ln J, so v = v0 J; sheet 1.4).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

from .model import (I3, PFLOOR, R_TOL_REL, SQ23, c2m, elastic_strain_scale, energy, m2t, plastic, t2m)
from .params import Params

FTOL_REL = 1.0e-8        # |F| / (M |p|) accepted as "on the surface"
OUTSIDE_REL = 1.0e-6     # F / (M |p|) above this at a segment start -> status 'outside'
FLOOR_TOL = 1.0e-8       # (p + p_min) / p_min >= -FLOOR_TOL: "on the floor" (sheet 9.7, S.55)
FLOOR_OUT = 1.0e-6       # (p + p_min) / p_min above this at a segment start -> status 'floor_outside'
PM = t2m(PFLOOR)         # Mandel form of the floor direction P = 1/3 (S.55)
# 'floor_rule' = the rate form (S.60) of the split with the image-pressure rule (S.57) (sheet 9.7, round 4)
FLOOR_MODES = ("floor", "plastic_floor", "floor_rule")
PLASTIC_MODES = ("plastic", "plastic_floor", "floor_rule")
FAP_TOL = 1.0e-12        # f:a^e:P < -FAP_TOL |f||a^e||P|: the (S.60) condition 'the floor pushes the stress out of F'


@dataclass
class State:
    sigma: np.ndarray
    eps_e: np.ndarray
    pi_i: float
    v: float
    D: float = 0.0            # plastic dissipation over the last increment, int sigma:deps^p
    eps_p_v: float = 0.0      # accumulated
    eps_p_s: float = 0.0      # accumulated
    flags: dict = field(default_factory=dict)


def _pack(st: State, floor=False):
    eps = st.flags.get("eps_total", np.zeros((3, 3)))
    y = np.concatenate([t2m(st.eps_e), [st.pi_i, st.v, st.flags.get("W", 0.0),
                                        st.flags.get("W_abs", 0.0),
                                        st.flags.get("Dp_total", 0.0),
                                        st.eps_p_v, st.eps_p_s], t2m(eps)])
    if floor:                 # y[19] = eps^f_v = int lam_f' (the floor mechanism, S.55); only when p_min > 0
        y = np.concatenate([y, [st.flags.get("eps_f_v", 0.0)]])
    return y


class _Increment:
    def __init__(self, P: Params, st: State, deps, smask, dsig, rtol, kin, floor_rule=True):
        self.P = P
        self.kin = kin
        self.v0 = st.flags.get("v0", st.v)
        # exponential v-update anchored at the increment start (exp of a sum: exact composition)
        self.v_start = float(st.v)
        self.tr_start = float(np.trace(st.flags.get("eps_total", np.zeros((3, 3)))))
        self.dm = t2m(np.asarray(deps, dtype=float))
        self.smask = np.zeros(6, dtype=bool) if smask is None else np.asarray(smask, dtype=bool)
        self.dsm = np.zeros(6) if dsig is None else t2m(np.asarray(dsig, dtype=float))
        self.S = np.where(~self.smask)[0]
        self.T = np.where(self.smask)[0]
        self.pmin = P.pmin
        self.floor = self.pmin > 0.0
        self.rule = self.floor and bool(floor_rule)   # (S.60) on the floor where it applies; False = (S.55) only
        self.s60_b = 0          # (S.60): branch of the trial floor rate <P:a:eps'> (1 = active), frozen per segment
        self.ptr_tol = 0.0
        ny = 20 if self.floor else 19
        emax = max(float(np.max(np.abs(self.dm))), 1e-9)
        pabs = max(abs(float(np.trace(st.sigma))) / 3.0, 1e-3)
        e_sc = elastic_strain_scale(P, pabs)
        at = np.empty(ny)
        at[0:6] = 1e-1 * rtol * e_sc
        at[6] = 1e-2 * rtol * max(abs(st.pi_i), 1e-3)
        at[7] = 1e-2 * rtol
        at[8:11] = 1e-2 * rtol * pabs * emax
        at[11:13] = 1e-2 * rtol * emax
        at[13:19] = 1e-2 * rtol * emax
        if self.floor:
            at[19] = 1e-2 * rtol * emax
        self.atol = at
        self.ny = ny
        self.typ = np.concatenate([np.full(6, 1e-3), [max(abs(st.pi_i), 1e-3), st.v],
                                   np.zeros(5), np.full(6, 1e-3)] + ([np.zeros(1)] if self.floor else []))
        self.nfev = 0
        self.w_force = None      # planar cap: the side (w = 1 or 0) chosen for this segment
        self.cap_side = 0

    # ------------------------------------------------------------------
    def vof(self, y):
        """Specific volume v = v_n exp(tr eps - tr eps_n) from the total strain in y[13:19]."""
        return self.v_start * math.exp(float(y[13] + y[14] + y[15]) - self.tr_start)

    def _solve_mixed(self, Cm):
        d = self.dm.copy()
        if len(self.T):
            d[self.T] = 0.0
            rhs = self.dsm[self.T] - Cm[np.ix_(self.T, self.S)] @ self.dm[self.S]
            d[self.T] = np.linalg.solve(Cm[np.ix_(self.T, self.T)], rhs)
        return d

    def rates(self, y, mode):
        self.nfev += 1
        P = self.P
        ee = m2t(y[0:6])
        pi, v = float(y[6]), self.vof(y)
        el = energy(ee, P)
        Am = c2m(el.a4)
        pq = plastic(el.sig, pi, v, P, w_override=self.w_force)
        fm, qm = t2m(pq.f), t2m(pq.qflow)
        sm = t2m(el.sig)
        info = dict(pq=pq, sig=el.sig)
        lam_f = 0.0
        if mode == "elastic":
            d = self._solve_mixed(Am)
            N = float(fm @ Am @ d)
            lam = 0.0
            info.update(N=N, den=float("nan"))
        elif mode == "plastic":
            Aq = Am @ qm
            fA = fm @ Am
            den = float(fm @ Aq) + pq.H
            Cm = Am - np.outer(Aq, fA) / den
            d = self._solve_mixed(Cm)
            N = float(fA @ d)
            lam = N / den
            info.update(N=N, den=den)
        elif mode == "floor":                   # (S.55), floor only: lam_f = P:a^e:eps' / P:a^e:P
            AP = Am @ PM
            PA = PM @ Am
            PAP = float(PM @ AP)
            Cm = Am - np.outer(AP, PA) / PAP
            d = self._solve_mixed(Cm)
            lam = 0.0
            lam_f = float(PA @ d) / PAP
            N = float(fm @ Am @ (d - lam_f * PM))        # f:sigma' (pi_i frozen): the yield-approach rate
            info.update(N=N, den=float("nan"), PAP=PAP)
        elif mode == "plastic_floor":           # (S.55), both mechanisms (Koiter)
            Aq = Am @ qm
            fA = fm @ Am
            AP = Am @ PM
            PA = PM @ Am
            G = np.array([[float(fA @ qm) + pq.H, float(fA @ PM)], [float(PA @ qm), float(PA @ PM)]])
            Gi = np.linalg.inv(G)
            Cm = Am - (Gi[0, 0] * np.outer(Aq, fA) + Gi[0, 1] * np.outer(Aq, PA)
                       + Gi[1, 0] * np.outer(AP, fA) + Gi[1, 1] * np.outer(AP, PA))
            d = self._solve_mixed(Cm)
            lam, lam_f = (float(x) for x in Gi @ np.array([fA @ d, PA @ d]))
            den = float(G[0, 0] - G[0, 1] * G[1, 0] / G[1, 1])     # Schur complement: the plastic denominator
            N = lam * den
            info.update(N=N, den=den, PAP=float(G[1, 1]), G=G)
        elif mode == "floor_rule":              # (S.60), sheet 9.7 round 4: the split with the rule (S.57)
            # lam_f^tr = <P:a:eps'>/(P:a:P)           (step 1f; branch s60_b frozen for the segment)
            # lam      = f:a:(eps' - lam_f^tr P)/(f:a:q + H)
            # lam_f    = P:a:(eps' - lam q)/(P:a:P)   (p' = 0)
            # pi_i'    = sqrt(2/3) h lam (pi_i* - pi_i) Omega + (f:a:P)(lam_f - lam_f^tr)/F_pi
            # every multiplier is linear in eps' on a branch: lam = l.eps', lam_f = m.eps'
            Aq = Am @ qm
            fA = fm @ Am
            AP = Am @ PM
            PA = PM @ Am
            PAP = float(PM @ AP)
            fAP = float(fA @ PM)
            PAq = float(PA @ qm)
            den = float(fA @ qm) + pq.H
            b = float(self.s60_b)
            lv = (fA - b * (fAP / PAP) * PA) / den
            mv = (PA - PAq * lv) / PAP
            Cm = Am - np.outer(Aq, lv) - np.outer(AP, mv)
            d = self._solve_mixed(Cm)
            lam = float(lv @ d)
            lam_f = float(mv @ d)
            ptr = float(PA @ d)                  # P:a:eps' (the trial floor rate times P:a:P)
            lam_ftr = b * ptr / PAP
            N = lam * den                        # f:a:(eps' - lam_f^tr P)
            info.update(N=N, den=den, PAP=PAP, fAP=fAP, lam_ftr=lam_ftr, ptr=ptr,
                        rule=(lam_f - lam_ftr) * PAP)
        else:
            raise ValueError(mode)
        dee = d - lam * qm
        if lam_f != 0.0:
            dee = dee - lam_f * PM
        dpi = SQ23 * P.h * lam * (pq.pistar - pi) * pq.Omega
        if mode == "floor_rule":                 # the reset term of (S.60): F' = 0 with p' = 0
            dpi = dpi + info["fAP"] * (lam_f - info["lam_ftr"]) / pq.Fpi
        trd = float(d[0] + d[1] + d[2])
        dv = v * trd                      # (S.40) v' = v tr eps' (diagnostic copy in y[7])
        sw = float(sm @ d)
        dD = lam * float(sm @ qm)
        info.update(lam=lam, dD=dD, d=d, dee=dee, dpi=dpi, a4=el.a4)
        dy = np.concatenate([dee, [dpi, dv, sw, abs(sw), dD, lam * float(np.trace(pq.qflow)),
                                   lam * SQ23 * pq.Omega], d])
        if self.floor:
            info.update(lam_f=lam_f, pdot=float(PM @ Am @ dee), p=el.p)
            dy = np.concatenate([dy, [lam_f]])
        return dy, info

    def jac(self, y, mode):
        f0, _ = self.rates(y, mode)
        J = np.zeros((self.ny, self.ny))
        for j in (0, 1, 2, 3, 4, 5, 6, 13, 14, 15):   # v depends on y only via tr eps (y[13:16])
            hj = 1e-7 * max(abs(y[j]), self.typ[j])
            yp = y.copy()
            yp[j] += hj
            fp, _ = self.rates(yp, mode)
            J[:, j] = (fp - f0) / hj
        return J

    # ------------------------------------------------------------------
    def Frel(self, y):
        el = energy(m2t(y[0:6]), self.P, tangent=False)
        pq = plastic(el.sig, float(y[6]), self.vof(y), self.P, w_override=self.w_force)
        return pq.F / (self.P.M * abs(pq.p)), pq

    def _eta_rate(self, y, w):
        """d(eta)/dt on the plastic branch with the cap weight forced to w (planar cap)."""
        old = self.w_force
        self.w_force = w
        try:
            _, ip = self.rates(y, "plastic")
        finally:
            self.w_force = old
        pq = ip["pq"]
        dsig = np.einsum("ijkl,kl->ij", ip["a4"], m2t(ip["dee"]))
        pdot = float(np.trace(dsig)) / 3.0
        return ((pq.Fp - pq.eta) / pq.p) * pdot + (pq.Fpi / pq.p) * ip["dpi"], ip

    def decide(self, y):
        self.w_force = None
        self.cap_side = 0
        Fr, pq = self.Frel(y)
        if self.P.cap == "planar" and Fr >= -FTOL_REL and                 abs(pq.eta - self.P.c1 * self.P.M) <= 1e-10 * self.P.M:
            # on the Q-corner eta = chi_cap M: pick the side whose flow moves AWAY from it
            e1, i1 = self._eta_rate(y, 1.0)
            e0, i0 = self._eta_rate(y, 0.0)
            if i1["N"] > 0 and i0["N"] > 0:
                if e1 < 0.0 < e0:
                    return ("stop", "cap_sliding"), Fr
                if e1 >= 0.0:
                    self.w_force, self.cap_side = 1.0, 1
                elif e0 <= 0.0:
                    self.w_force, self.cap_side = 0.0, -1
                Fr, pq = self.Frel(y)
        if self.P.cap == "planar" and self.cap_side == 0:
            # off the corner: freeze the side for the whole segment so that the RHS the
            # solver sees is smooth; the crossing is caught by the 'cap_switch' event
            self.w_force = 1.0 if pq.eta >= self.P.c1 * self.P.M else 0.0
            self.cap_side = 1 if self.w_force == 1.0 else -1
            Fr, pq = self.Frel(y)
        if self.floor:
            pf = (pq.p + self.pmin) / self.pmin
            if pf > FLOOR_OUT:
                return ("stop", "floor_outside"), Fr
            if pf >= -FLOOR_TOL:
                return self._decide_floor(y, Fr, pq), Fr
        if Fr < -FTOL_REL:
            return "elastic", Fr
        if Fr > OUTSIDE_REL:
            return ("stop", "outside"), Fr
        if self.P.N > 0 and pq.B <= 0.0:
            return ("stop", "B_guard"), Fr
        _, ip = self.rates(y, "plastic")
        el = energy(m2t(y[0:6]), self.P)
        scale = 1e-12 * np.linalg.norm(t2m(pq.f)) * np.linalg.norm(c2m(el.a4)) * \
            max(np.linalg.norm(self.dm), 1e-30)
        if ip["N"] > scale and ip["den"] > 0.0:
            if pq.vertex:
                dd = m2t(ip["d"])
                dev = dd - np.trace(dd) / 3.0 * np.eye(3)
                if np.linalg.norm(dev) > 1e-12 * max(np.linalg.norm(dd), 1e-300):
                    # vertex rule (sheet 3.2) has no deviatoric gradient: F would leave 0
                    return ("stop", "vertex_nonisotropic"), Fr
            return "plastic", Fr
        _, ie = self.rates(y, "elastic")
        if ie["N"] <= scale:
            return "elastic", Fr
        if ip["den"] <= 0.0:
            return ("stop", "den_nonpositive"), Fr
        return ("stop", "no_consistent_mode"), Fr

    def _decide_floor(self, y, Fr, pq):
        """Mode at a state ON the floor p = -p_min (S.55): the KKT pair lam, lam_f >= 0 with F' <= 0 / p' <= 0
        for the inactive mechanism.  Tried in the order plastic+floor, plastic, floor, elastic; an active
        multiplier must exceed +tol*scale, an inactive constraint rate must not exceed +tol*scale (sheet 12
        tie-break, extended to the second mechanism)."""
        P = self.P
        if Fr > OUTSIDE_REL:
            return ("stop", "outside")
        el = energy(m2t(y[0:6]), P)
        Am = c2m(el.a4)
        dn = max(np.linalg.norm(self.dm), 1e-30)
        sc_f = 1e-12 * np.linalg.norm(t2m(pq.f)) * np.linalg.norm(Am) * dn
        sc_p = 1e-12 * np.linalg.norm(PM) * np.linalg.norm(Am) * dn
        on_yield = Fr >= -FTOL_REL
        if on_yield:
            if P.N > 0 and pq.B <= 0.0:
                return ("stop", "B_guard")
            if self.rule:
                b = self._s60_branch(y, Am, pq, sc_f, sc_p)
                if b is not None:
                    self.s60_b = b
                    return "floor_rule"
            _, ib = self.rates(y, "plastic_floor")
            if ib["lam"] * ib["den"] > sc_f and ib["lam_f"] * ib["PAP"] > sc_p and ib["den"] > 0.0:
                return "plastic_floor"
            _, ip = self.rates(y, "plastic")
            if ip["N"] > sc_f and ip["den"] > 0.0 and ip["pdot"] <= sc_p:
                return "plastic"
        _, ifl = self.rates(y, "floor")
        if ifl["lam_f"] * ifl["PAP"] > sc_p and (not on_yield or ifl["N"] <= sc_f):
            return "floor"
        _, ie = self.rates(y, "elastic")
        if ie["pdot"] <= sc_p and (not on_yield or ie["N"] <= sc_f):
            return "elastic"
        return ("stop", "no_consistent_mode_floor")

    def _s60_branch(self, y, Am, pq, sc_f, sc_p):
        """(S.60) validity at a state on the floor and on F (sheet 9.7, round 4): lam' > 0, lam_f' - lam_f^tr' > 0 and
        f:a^e:P < 0 (and the plastic denominator > 0).  Returns the branch b of <P:a:eps'> (1: the trial floor rate is
        active) for which (S.60) holds self-consistently, or None (then (S.55) applies, unchanged).  Under strain
        control exactly one branch is consistent; under mixed control b = 1 is tried first.  The probe does not count
        in nfev, so a floor path on which (S.60) never applies is bit-identical to the round-3b O1."""
        fAP = float(t2m(pq.f) @ Am @ PM)
        if not fAP < -FAP_TOL * np.linalg.norm(t2m(pq.f)) * np.linalg.norm(Am) * np.linalg.norm(PM):
            return None
        n0, old = self.nfev, self.s60_b
        try:
            for b in (1, 0):
                self.s60_b = b
                _, ir = self.rates(y, "floor_rule")
                consistent = ir["ptr"] > sc_p if b == 1 else ir["ptr"] <= sc_p
                if not consistent:
                    continue
                if ir["den"] > 0.0 and ir["N"] > sc_f and ir["rule"] > sc_p:
                    self.ptr_tol = sc_p
                    return b
                return None
            return None
        finally:
            self.nfev, self.s60_b = n0, old

    def events(self, mode, y0):
        P = self.P
        evs, names = [], []

        def add(fun, name, direction, terminal=True, skip_zero=True):
            if skip_zero:
                g0 = fun(0.0, y0)
                if abs(g0) <= 1e-14:
                    return
            fun.terminal = terminal
            fun.direction = direction
            evs.append(fun)
            names.append(name)

        if mode not in FLOOR_MODES:
            def e_p(t, y):
                el = energy(m2t(y[0:6]), P, tangent=False)
                return -el.p - 1e-6 * P.p_ref
            add(e_p, "p_to_zero", -1)
            if self.floor:                      # the floor is reached: p + p_min crosses 0 upward (S.55)
                pm = self.pmin

                def e_fl(t, y):
                    return (energy(m2t(y[0:6]), P, tangent=False).p + pm) / pm
                add(e_fl, "floor_hit", 1)
        else:
            def e_rel(t, y):                    # the floor mechanism's multiplier turns negative: release
                ii = self.rates(y, mode)[1]
                return ii["lam_f"] * ii["PAP"]
            add(e_rel, "floor_release", -1)
            if mode == "floor_rule":            # the edges of the (S.60) set, and the kink of <P:a:eps'>
                def e_rule(t, y):               # lam_f' - lam_f^tr' -> 0: the return no longer relaxes p
                    return self.rates(y, mode)[1]["rule"]
                add(e_rule, "rule_off", -1)

                def e_fap(t, y):                # f:a^e:P -> 0: the floor no longer pushes the stress out of F
                    return self.rates(y, mode)[1]["fAP"]
                add(e_fap, "rule_off_fAP", 1)

                # P:a:eps' changes sign: the other branch of lam_f^tr (a hysteresis of ptr_tol: a pure deviatoric
                # rate under BA06 alpha0 = 0 has P:a:eps' = 0 identically and must not chatter on its round-off)
                if self.s60_b == 1:
                    def e_tr(t, y):
                        return self.rates(y, mode)[1]["ptr"]
                    add(e_tr, "rule_branch", -1)
                else:
                    tol = self.ptr_tol

                    def e_tr(t, y):
                        return self.rates(y, mode)[1]["ptr"] - tol
                    add(e_tr, "rule_branch", 1)

        if mode in ("elastic", "floor"):
            F0, _ = self.Frel(y0)
            thr = max(F0, 0.0)

            def e_F(t, y):
                return self.Frel(y)[0] - thr
            add(e_F, "yield", 1, skip_zero=False)
            return evs, names

        def e_N(t, y):
            return self.rates(y, mode)[1]["N"]
        add(e_N, "unload", -1)

        def e_den(t, y):
            return self.rates(y, mode)[1]["den"]
        add(e_den, "den_nonpositive", -1)

        def e_H(t, y):
            pq = self.Frel(y)[1]
            return (pq.pistar - float(y[6])) / abs(float(y[6]))
        add(e_H, "H_zero", 0, terminal=False)

        if P.N > 0:
            def e_B(t, y):
                return self.Frel(y)[1].B - 1e-8
            add(e_B, "B_guard", -1)

        y0pq = self.Frel(y0)[1]
        if not y0pq.vertex:
            def e_vx(t, y):
                pq = self.Frel(y)[1]
                return pq.R / abs(pq.p) - 10.0 * R_TOL_REL
            add(e_vx, "vertex_reached", -1)

        if P.cap == "planar":
            side = self.cap_side
            if side != 0:
                def e_c(t, y):
                    return side * (self.Frel(y)[1].eta - P.c1 * P.M)
                add(e_c, "cap_switch", -1, skip_zero=False)
            else:
                def e_c(t, y):
                    return self.Frel(y)[1].eta - P.c1 * P.M
                add(e_c, "cap_switch", 0)
        elif P.cap == "smooth":
            def e_c1(t, y):
                return self.Frel(y)[1].eta - P.c1 * P.M
            add(e_c1, "cap_eta1", 0)

            def e_c2(t, y):
                return self.Frel(y)[1].eta - P.c2 * P.M
            add(e_c2, "cap_eta2", 0)
        return evs, names


def integrate_increment(P: Params, st: State, deps, smask=None, dsig=None, rtol=1e-10,
                        kin="small", max_segments=500, floor_rule=True) -> State:
    """Integrate one increment over pseudo-time [0, 1] from state `st`.

    deps  : 3x3 strain increment (strain-controlled components; tensor strains).
    smask : optional Mandel-order bool[6] (xx,yy,zz,xy,yz,xz); True = stress-controlled
            component with stress increment dsig (default 0).
    kin   : 'small' or 'log' (eps = log strain).  Both use v = v_n exp(tr deps) (sheet 1.2/1.4;
            the pre-G2 small-strain rule v' = v0 tr eps' is superseded); kin is recorded only.
    floor_rule : True (default, sheet round 4): on the floor and on F, the rate form (S.60) wherever its
            validity set holds, the Koiter form (S.55) elsewhere.  False: (S.55) only (the round-3b O1;
            used to measure the gap between the two rate forms).  No effect when p_min = 0.
    """
    inc = _Increment(P, st, deps, smask, dsig, rtol, kin, floor_rule)
    y = _pack(st, inc.floor)
    t = 0.0
    status = "ok"
    segments = []
    Hz = list(st.flags.get("H_zero", []))
    notes = []
    max_F = 0.0
    min_rate = float("inf")
    min_den = float("inf")
    max_floor_drift = 0.0
    min_lam_f = float("inf")
    inf = float("inf")
    rule_diag = dict(min_lam=inf, min_lam_f=inf, max_lam_ftr=-inf, min_rule=inf, max_fAP=-inf, min_den=inf,
                     max_F=0.0, max_floor_dp=0.0)
    stall = 0
    mode = None
    while t < 1.0:
        if len(segments) >= max_segments:
            status = "max_segments"
            break
        mode, Fr = inc.decide(y)
        if isinstance(mode, tuple):
            status = mode[1]
            break
        evs, names = inc.events(mode, y)
        try:
            sol = solve_ivp(lambda tt, yy, m=mode: inc.rates(yy, m)[0], (t, 1.0), y,
                            method="Radau", rtol=rtol, atol=inc.atol, events=evs,
                            jac=lambda tt, yy, m=mode: inc.jac(yy, m))
        except FloatingPointError as exc:
            # a Radau stage evaluated a state outside the model domain (p >= 0; e.g. a HAR trial past the edge of dom Psi with
            # p_min = 0, sheet 2.3 / 13.13): the p_to_zero stop of the event machinery, reached by a stage before its event.
            # The state returned is the last accepted one (t, y unchanged).  Round 3b verify fix.
            status = "p_to_zero"
            notes.append("domain: " + str(exc))
            break
        if sol.status == -1:
            t, y = float(sol.t[-1]), sol.y[:, -1].copy()
            status = "solver_failed"
            notes.append(sol.message)
            break
        fired = None
        if sol.status == 1:
            tend = sol.t[-1]
            for k, te in enumerate(sol.t_events):
                if evs[k].terminal and len(te) and abs(te[-1] - tend) <= 1e-14 * max(1.0, tend):
                    fired = names[k]
                    break
        for k, nm in enumerate(names):
            if nm == "H_zero":
                for te, ye in zip(sol.t_events[k], sol.y_events[k]):
                    el = energy(m2t(ye[0:6]), P, tangent=False)
                    pq = plastic(el.sig, float(ye[6]), float(ye[7]), P)
                    Hz.append(dict(eps_p_s=float(ye[12]), D=pq.dilatancy,
                                   chi_psi=P.chi * pq.psi, p=pq.p, q=pq.q, eta_ratio=pq.q / -pq.p))
        # diagnostics at the accepted output points
        for k in range(1, sol.y.shape[1]):
            yk = sol.y[:, k]
            dyk, info = inc.rates(yk, mode)
            pq = info["pq"]
            if mode in PLASTIC_MODES:
                max_F = max(max_F, abs(pq.F) / (P.M * abs(pq.p)))
                min_den = min(min_den, info["den"])
                rate_scale = abs(pq.p) * max(np.linalg.norm(inc.dm), 1e-30)
                min_rate = min(min_rate, info["dD"] / rate_scale)
            if mode in FLOOR_MODES:
                max_floor_drift = max(max_floor_drift, abs(pq.p + inc.pmin) / inc.pmin)
                min_lam_f = min(min_lam_f, info["lam_f"] / max(np.linalg.norm(inc.dm), 1e-30))
            if mode == "floor_rule":
                dn = max(np.linalg.norm(inc.dm), 1e-30)
                rd = rule_diag
                rd["min_lam"] = min(rd["min_lam"], info["lam"] / dn)
                rd["min_lam_f"] = min(rd["min_lam_f"], info["lam_f"] / dn)
                rd["max_lam_ftr"] = max(rd["max_lam_ftr"], info["lam_ftr"] / dn)
                rd["min_rule"] = min(rd["min_rule"], info["rule"] / (info["PAP"] * dn))
                rd["max_fAP"] = max(rd["max_fAP"], info["fAP"])
                rd["min_den"] = min(rd["min_den"], info["den"])
                rd["max_F"] = max(rd["max_F"], abs(pq.F))
                rd["max_floor_dp"] = max(rd["max_floor_dp"], abs(pq.p + inc.pmin))
        segments.append(dict(mode=mode, t0=t, t1=float(sol.t[-1]), event=fired,
                             steps=int(len(sol.t) - 1)))
        if mode == "floor_rule":
            segments[-1]["s60_b"] = inc.s60_b
        t_new, y_new = float(sol.t[-1]), sol.y[:, -1].copy()
        stall = stall + 1 if t_new - t <= 1e-13 else 0
        t, y = t_new, y_new
        if stall > 50:
            status = "chatter"
            break
        if fired in ("p_to_zero", "den_nonpositive", "B_guard", "vertex_reached"):
            status = fired
            break
        if fired is None and t < 1.0:
            status = "solver_stopped"
            break

    ee = m2t(y[0:6])
    el = energy(ee, P, tangent=False)
    Fr = float("nan")
    try:
        pq = plastic(el.sig, float(y[6]), inc.vof(y), P)
        Fr = pq.F / (P.M * abs(pq.p))
    except FloatingPointError:
        pq = None
    flags = dict(st.flags)
    last_mode = segments[-1]["mode"] if segments else (mode if isinstance(mode, str) else None)
    if inc.floor:
        efv = float(y[19])
        nseg_f = sum(1 for sg in segments if sg["mode"] in FLOOR_MODES)
        flags.update(eps_f_v=efv, W_f=inc.pmin * efv, floor=(last_mode in FLOOR_MODES),
                     at_floor=bool(el.p > -inc.pmin * (1.0 + 1e-10)), floor_segments=nseg_f,
                     n_floor_increments=st.flags.get("n_floor_increments", 0) + (1 if nseg_f else 0),
                     max_floor_drift=max_floor_drift, min_lam_f_rel=min_lam_f, p_min=inc.pmin)
        # (S.60), sheet 9.7 round 4: segments integrated in the rule's rate form (the split's pattern -Pf*)
        nseg_r = sum(1 for sg in segments if sg["mode"] == "floor_rule")
        flags.update(rule_segments=nseg_r, floor_rule=(last_mode == "floor_rule"), rule_branch=(segments[-1].get("s60_b") if segments else None),
                     n_rule_increments=st.flags.get("n_rule_increments", 0) + (1 if nseg_r else 0),
                     rule_diag=(rule_diag if nseg_r else None))
    flags.update(
        status=status, plastic=(last_mode in PLASTIC_MODES), mode=last_mode,
        F_rel=Fr, max_F_rel=max_F, min_den=min_den, min_Dp_rate_rel=min_rate,
        segments=segments, n_segments=len(segments), nfev=inc.nfev, notes=notes,
        v0=inc.v0, W=float(y[8]), W_abs=float(y[9]), Dp_total=float(y[10]),
        eps_total=m2t(y[13:19]), H_zero=Hz, kin=kin, v_ode=float(y[7]),
        psi=(pq.psi if pq else float("nan")), pistar=(pq.pistar if pq else float("nan")),
        H=(pq.H if pq else float("nan")))
    return State(sigma=el.sig, eps_e=ee, pi_i=float(y[6]), v=inc.vof(y),
                 D=float(y[10]) - st.flags.get("Dp_total", 0.0),
                 eps_p_v=float(y[11]), eps_p_s=float(y[12]), flags=flags)
