"""ADR-92 P0 -- the single-Gauss-point SANISAND / IMPL-EX oracle.

A numpy transcription of `ManzariDafalias`'s stress update at ONE Gauss point, exact
to the C++ at `ladruno` 3f003d110 (`SRC/material/nD/UWmaterials/ManzariDafalias.cpp`),
INCLUDING the code's departures from Dafalias & Manzari (2004):

  * `GetElasticModuli` uses `m_e_init`, never the current void ratio
    (`mUseCurrentVoidRatioInG` is false in every constructor);
  * the `D_factor` dilatancy sigmoid below `p = 0.05*P_atm`;
  * `ForwardEuler`'s block-scoped `r` -- identically zero -- reproduced VERBATIM and
    only in the `forward_euler` path;
  * `ModifiedEuler`'s `TolE = 1e-4` (mTolR ignored at `-honorTolR 0`) and its
    FORCE-ACCEPT at `dT == dT_min`;
  * `Stress_Correction`, which is ON by default (`mStressCorrectionInUse = true` in all
    four constructors, `ManzariDafalias.cpp:217/303/368/429`) -- the source-extraction
    memo said "inert", it is not;
  * `BackwardEuler_CPPM`'s low-`p` branch, whose Newton is disabled with a literal
    `errFlag = 0` (`:2264`) so it ALWAYS falls through to `explicit_integrator`.

Sign convention: everything inside is the material's own (compression POSITIVE), i.e.
`LadrunoSANISAND3D::setTrialStrain` negates before calling. CSV columns from
`probe_binary_triaxial.py` are in the element's convention and are negated on read.

Run with python3.12 + numpy.  `python3.12 sanisand_implex_oracle.py --gate all`
"""
from __future__ import annotations

import argparse
import csv
import glob
import math
import os
from dataclasses import dataclass, field

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

ONE3 = 1.0 / 3.0
TWO3 = 2.0 / 3.0
ROOT23 = math.sqrt(2.0 / 3.0)
SMALL = 1.0e-10
I1 = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])

# ---------------------------------------------------------------- tensor algebra
# ManzariDafalias.cpp:5077-5300.  Storage is (11,22,33,12,23,13); "contravariant"
# means stress-like, "covariant" means strain-like (engineering shear).


def trace(v):
    return v[0] + v[1] + v[2]


def dev(v):
    r = v.copy()
    p = ONE3 * trace(v)
    r[0] -= p
    r[1] -= p
    r[2] -= p
    return r


def single_dot(v1, v2):
    r = np.empty(6)
    r[0] = v1[0] * v2[0] + v1[3] * v2[3] + v1[5] * v2[5]
    r[1] = v1[3] * v2[3] + v1[1] * v2[1] + v1[4] * v2[4]
    r[2] = v1[5] * v2[5] + v1[4] * v2[4] + v1[2] * v2[2]
    r[3] = 0.5 * (v1[0] * v2[3] + v1[3] * v2[0] + v1[3] * v2[1]
                  + v1[1] * v2[3] + v1[5] * v2[4] + v1[4] * v2[5])
    r[4] = 0.5 * (v1[3] * v2[5] + v1[5] * v2[3] + v1[1] * v2[4]
                  + v1[4] * v2[1] + v1[4] * v2[2] + v1[2] * v2[4])
    r[5] = 0.5 * (v1[0] * v2[5] + v1[5] * v2[0] + v1[3] * v2[4]
                  + v1[4] * v2[3] + v1[5] * v2[2] + v1[2] * v2[5])
    return r


def dd_contr(v1, v2):
    return float(v1 @ v2 + v1[3:] @ v2[3:])


def dd_mixed(v1, v2):
    return float(v1 @ v2)


def norm_contr(v):
    return math.sqrt(dd_contr(v, v))


def to_cov(v):
    r = v.copy()
    r[3:] *= 2.0
    return r


def to_contra(v):
    r = v.copy()
    r[3:] *= 0.5
    return r


def stiffness(K, G):
    a = K + 4.0 * ONE3 * G
    b = K - 2.0 * ONE3 * G
    C = np.zeros((6, 6))
    C[0, 0] = C[1, 1] = C[2, 2] = a
    C[3, 3] = C[4, 4] = C[5, 5] = G
    C[0, 1] = C[0, 2] = C[1, 2] = b
    C[1, 0] = C[2, 0] = C[2, 1] = b
    return C


def compliance(K, G):
    a = 1.0 / (9 * K) + 1.0 / (3 * G)
    b = 1.0 / (9 * K) - 1.0 / (6 * G)
    c = 1.0 / G
    D = np.zeros((6, 6))
    D[0, 0] = D[1, 1] = D[2, 2] = a
    D[3, 3] = D[4, 4] = D[5, 5] = c
    D[0, 1] = D[0, 2] = D[1, 2] = b
    D[1, 0] = D[2, 0] = D[2, 1] = b
    return D


def macauley(x):
    return x if x > 0.0 else 0.0


# ---------------------------------------------------------------- the material

CONSTS = dict(G0=264.32, nu=0.3129, e_init=0.6944, Mc=1.3309, c=0.71,
              lambda_c=0.027, e0=0.83, ksi=0.45, P_atm=101.0, m=0.005,
              h0=1.3, ch=0.968, nb=3.5, A0=0.05, nd=5.75, z_max=12.5, cz=1100.0)


@dataclass
class StateDep:
    n: np.ndarray
    d: np.ndarray
    b: np.ndarray
    cos3t: float
    h: float
    psi: float
    aB: float
    aD: float
    b0: float
    A: float
    D: float
    B: float
    C: float
    R: np.ndarray


class Abandoned(Exception):
    """ModifiedEuler's bare `return` at p < m_Presidual with dT == dT_min."""


@dataclass
class Counters:
    steps: int = 0
    me_calls: int = 0
    me_substeps: int = 0
    me_force_accept: int = 0
    me_abandon: int = 0
    me_lowp_clamp: int = 0
    me_neg_dgamma: int = 0
    sc_iterations: int = 0
    be_calls: int = 0
    be_lowp_trial: int = 0
    be_newton_fail: int = 0
    be_recursive: int = 0
    be_explicit_fallback: int = 0
    be_maxsubstep: int = 0

    def add(self, other):
        for k in self.__dataclass_fields__:
            setattr(self, k, getattr(self, k) + getattr(other, k))


class Sanisand:
    """`LadrunoSANISAND` at one Gauss point, in the material's own sign convention."""

    def __init__(self, consts=None, scheme=1, Pmin=0.0101, Presidual=0.0,
                 TolF=1e-10, TolR=1e-10, honor_tolR=False, elast_flag=1,
                 stress_correction=True):
        c = dict(CONSTS)
        if consts:
            c.update(consts)
        for k, v in c.items():
            setattr(self, "m_" + k, v)
        self.scheme = scheme
        self.m_Pmin = Pmin
        self.m_Presidual = Presidual
        self.mTolF = TolF
        self.mTolR = TolR
        self.mHonorTolRInME = honor_tolR
        self.mElastFlag = elast_flag
        self.mStressCorrectionInUse = stress_correction
        self.cnt = Counters()
        self._z6 = np.zeros(6)
        # trial
        self.eps = np.zeros(6)
        self.sig = np.zeros(6)
        self.epsE = np.zeros(6)
        self.alpha = np.zeros(6)
        self.alpha_in = np.zeros(6)
        self.fabric = np.zeros(6)
        self.dGamma = 0.0
        # committed
        self.eps_n = np.zeros(6)
        self.sig_n = np.zeros(6)
        self.epsE_n = np.zeros(6)
        self.alpha_n = np.zeros(6)
        self.alpha_in_n = np.zeros(6)
        self.fabric_n = np.zeros(6)
        self.dGamma_n = 0.0
        self.void_ratio = self.m_e_init
        self.mK, self.mG = 0.0, 0.0
        self._refresh_moduli()

    # ---- committed-state helpers
    def _refresh_moduli(self):
        self.mK, self.mG = self.elastic_moduli(self.sig_n, self.void_ratio)

    def set_committed(self, eps, sig, epsE, alpha, fabric, alpha_in, dGamma=0.0):
        self.eps_n = np.array(eps, float)
        self.sig_n = np.array(sig, float)
        self.epsE_n = np.array(epsE, float)
        self.alpha_n = np.array(alpha, float)
        self.fabric_n = np.array(fabric, float)
        self.alpha_in_n = np.array(alpha_in, float)
        self.dGamma_n = dGamma
        self.eps = self.eps_n.copy()
        self.sig = self.sig_n.copy()
        self.epsE = self.epsE_n.copy()
        self.alpha = self.alpha_n.copy()
        self.fabric = self.fabric_n.copy()
        self.alpha_in = self.alpha_in_n.copy()
        self.dGamma = dGamma
        self.void_ratio = self.m_e_init - (1 + self.m_e_init) * trace(self.eps_n)
        self._refresh_moduli()

    def commit(self):
        """ManzariDafalias::commitState (:467-511)."""
        self.alpha_in_n = self.alpha_in.copy()
        self.sig_n = self.sig.copy()
        self.eps_n = self.eps.copy()
        self.epsE_n = self.epsE.copy()
        self.alpha_n = self.alpha.copy()
        self.fabric_n = self.fabric.copy()
        self.dGamma_n = self.dGamma
        self.void_ratio = self.m_e_init - (1 + self.m_e_init) * trace(self.eps)
        self._refresh_moduli()

    def eps_p_n(self):
        return self.eps_n - self.epsE_n

    # ---- kernel
    def elastic_moduli(self, sigma, e):
        pn = ONE3 * trace(sigma)
        if pn <= self.m_Pmin:
            pn = self.m_Pmin
        eG = self.m_e_init                      # mUseCurrentVoidRatioInG is false
        G = self.m_G0 * self.m_P_atm * (2.97 - eG) ** 2 / (1 + eG)
        if self.mElastFlag != 0:
            G *= math.sqrt(pn / self.m_P_atm)
        K = TWO3 * (1 + self.m_nu) / (1 - 2 * self.m_nu) * G
        return K, G

    def get_F(self, sigma, alpha):
        s = dev(sigma)
        p = ONE3 * trace(sigma) + self.m_Presidual
        s = s - p * alpha
        return norm_contr(s) - ROOT23 * self.m_m * p

    def get_psi(self, e, p):
        return e - (self.m_e0 - self.m_lambda_c * (p / self.m_P_atm) ** self.m_ksi)

    def normal_to_yield(self, sigma, alpha):
        p = ONE3 * trace(sigma) + self.m_Presidual
        if abs(p) < SMALL:
            return np.zeros(6)
        n = dev(sigma) - p * alpha
        nn = norm_contr(n)
        if nn < SMALL:
            nn = 1.0
        return n / nn

    @staticmethod
    def lode(n):
        c = math.sqrt(6.0) * trace(single_dot(n, single_dot(n, n)))
        return min(1.0, max(-1.0, c))

    def g(self, cos3t):
        c = self.m_c
        return 2 * c / ((1 + c) - (1 - c) * cos3t)

    def state_dependent(self, sigma, alpha, fabric, e, alpha_in) -> StateDep:
        p = ONE3 * trace(sigma) + self.m_Presidual
        if p < SMALL:
            p = SMALL
        n = self.normal_to_yield(sigma, alpha)
        aain = dd_contr(alpha - alpha_in, n)
        psi = self.get_psi(e, p)
        cos3t = self.lode(n)
        gc = self.g(cos3t)
        try:
            aB = gc * self.m_Mc * math.exp(-self.m_nb * psi) - self.m_m
            aD = gc * self.m_Mc * math.exp(self.m_nd * psi) - self.m_m
        except OverflowError:
            # psi has run away (the step is far past what the substepper can hold,
            # e.g. d(eps_z) = 2e-3 at p0 = 100, where the BINARY's own global solve
            # stalls too -- probe tx_p100 n100 ez0.20). Report it as an abandoned
            # step, which is what the C++ would print as garbage, not crash on.
            raise Abandoned(f"exp overflow in GetStateDependent: psi = {psi:.3e}")
        b0 = self.m_G0 * self.m_h0 * (1.0 - self.m_ch * e) / math.sqrt(p / self.m_P_atm)
        d = ROOT23 * aD * n - alpha
        b = ROOT23 * aB * n - alpha
        h = 1.0e10 if abs(aain) < SMALL else b0 / aain
        A = self.m_A0 * (1 + macauley(dd_contr(fabric, n)))
        D = A * dd_contr(d, n)
        if p < 0.05 * self.m_P_atm:
            D *= 1.0 / (1.0 + math.exp(7.6349 - 7.2713 * 101.0 / self.m_P_atm * p))
        B = 1.0 + 1.5 * (1 - self.m_c) / self.m_c * gc * cos3t
        C = 3.0 * math.sqrt(1.5) * (1 - self.m_c) / self.m_c * gc
        R = B * n - C * (single_dot(n, n) - ONE3 * I1) + ONE3 * D * I1
        return StateDep(n, d, b, cos3t, h, psi, aB, aD, b0, A, D, B, C, R)

    # ---- integrate()  (:957-994)
    def integrate(self, eps_trial):
        self.eps = np.array(eps_trial, float)
        Ce = stiffness(self.mK, self.mG)
        trial_dir = Ce @ (self.eps - self.eps_n)
        if dd_contr(self.alpha_n - self.alpha_in_n, trial_dir) < 0.0:
            self.alpha_in = self.alpha_n.copy()
        else:
            self.alpha_in = self.alpha_in_n.copy()
        self.cnt.steps += 1
        if self.mElastFlag == 0:
            out = self.elastic_integrator()
        elif self.scheme == 2:
            out = self.backward_euler(self.sig_n, self.eps_n, self.epsE_n,
                                      self.alpha_n, self.fabric_n, self.alpha_in,
                                      self.eps, self.mG, self.mK, 0)
        else:
            out = self.explicit_integrator(self.sig_n, self.eps_n, self.epsE_n,
                                           self.alpha_n, self.fabric_n,
                                           self.alpha_in, self.eps,
                                           self.mG, self.mK)
        (self.epsE, self.sig, self.alpha, self.fabric, self.dGamma,
         self.void_ratio) = out[:6]
        return self.sig

    def elastic_integrator(self):
        d = self.eps - self.eps_n
        vr = self.m_e_init - (1 + self.m_e_init) * trace(self.eps)
        epsE = self.epsE_n + d
        K, G = self.elastic_moduli(self.sig_n, vr)
        sig = self.sig_n + stiffness(K, G) @ d
        p = ONE3 * trace(sig) + self.m_Presidual
        alpha = dev(sig) / p if p > SMALL else self.alpha_n.copy()
        return epsE, sig, alpha, self.fabric_n.copy(), 0.0, vr

    # ---- explicit_integrator (:1031-1147): elastic / plastic split
    def explicit_integrator(self, cS, cE, cEE, cA, cF, ain, nE, G, K):
        vr = self.m_e_init - (1 + self.m_e_init) * trace(nE)
        dStrain = nE - cE
        epsE = cEE + dStrain
        aC = stiffness(K, G)
        dSigma = aC @ dStrain
        nStress = cS + dSigma
        f = self.get_F(nStress, cA)
        p = ONE3 * trace(nStress) + self.m_Presidual
        p_tr_pos = not (p < self.m_Presidual)
        if p_tr_pos and f <= self.mTolF:
            return epsE, nStress, cA.copy(), cF.copy(), 0.0, vr
        fn = self.get_F(cS, cA)
        pn = ONE3 * trace(cS) + self.m_Presidual
        if pn < self.m_Presidual:
            return epsE, self.m_Pmin * I1, np.zeros(6), cF.copy(), 0.0, vr
        if fn > self.mTolF:
            return self._exp_scheme(cS, cE, cEE, cA, cF, ain, nE, G, K, vr)
        if fn < -self.mTolF:
            a = self.intersection_factor(cS, cE, nE, cA, 0.0, 1.0)
            de = dStrain * a
            ds = aC @ de
            return self._exp_scheme(cS + ds, cE + de, cEE + de, cA, cF, ain, nE, G, K, vr)
        # |fn| < TolF
        nd = norm_contr(dSigma)
        ratio = dd_contr(self.normal_to_yield(cS, cA), dSigma) / (1.0 if nd == 0 else nd)
        if ratio > -math.sqrt(self.mTolF):
            return self._exp_scheme(cS, cE, cEE, cA, cF, ain, nE, G, K, vr)
        a = self.intersection_factor_unloading(cS, cE, nE, cA)
        de = dStrain * a
        ds = aC @ de
        return self._exp_scheme(cS + ds, cE + de, cEE + de, cA, cF, ain, nE, G, K, vr)

    def _exp_scheme(self, cS, cE, cEE, cA, cF, ain, nE, G, K, vr):
        if self.scheme == 5:
            fn = self.forward_euler
        else:                                # default: ModifiedEuler (:1060-1061)
            fn = self.modified_euler
        return fn(cS, cE, cEE, cA, cF, ain, nE, G, K, vr)

    def intersection_factor(self, cS, cE, nE, cA, a0, a1):
        si = nE - cE
        vr = self.m_e_init - (1 + self.m_e_init) * trace(cE + a0 * si)
        K, G = self.elastic_moduli(cS, vr)
        f0 = self.get_F(cS + a0 * (stiffness(K, G) @ si), cA)
        vr = self.m_e_init - (1 + self.m_e_init) * trace(cE + a1 * si)
        K, G = self.elastic_moduli(cS, vr)
        f1 = self.get_F(cS + a1 * (stiffness(K, G) @ si), cA)
        a = a0
        for i in range(1, 11):
            a = a1 - f1 * (a1 - a0) / (f1 - f0)
            f = self.get_F(cS + a * (stiffness(K, G) @ si), cA)
            if abs(f) < self.mTolF:
                break
            if f * f0 < 0:
                a1, f1 = a, f
            else:
                f1 = f1 * f0 / (f0 + f)
                a0, f0 = a, f
            if i == 10:
                a = 0.0
                break
        if a > 1 - SMALL:
            a = 1.0
        if a < SMALL:
            a = 0.0
        return a

    def intersection_factor_unloading(self, cS, cE, nE, cA):
        a0, a1, nSub = 0.0, 1.0, 20
        si = nE - cE
        vr = self.m_e_init - (1 + self.m_e_init) * trace(cE)
        K, G = self.elastic_moduli(cS, vr)
        dS = stiffness(K, G) @ si
        for i in range(1, nSub):
            a = a1 - (a1 - a0) / 2.0
            f = self.get_F(cS + a * dS, cA)
            if f > self.mTolF:
                a1 = a
            elif f < -self.mTolF:
                a0 = a
                break
            else:
                return a
        return self.intersection_factor(cS, cE, nE, cA, a0, a1)

    # ---- ForwardEuler (:1313-1362) -- with the zero-`r` defect VERBATIM
    def forward_euler(self, cS, cE, cEE, cA, cF, ain, nE, G, K, vr_ignored):
        cvr = self.m_e_init - (1 + self.m_e_init) * trace(cE)
        vr = self.m_e_init - (1 + self.m_e_init) * trace(nE)
        sd = self.state_dependent(cS, cA, cF, cvr, ain)
        dvol = trace(nE - cE)
        ddev = dev(nE - cE)
        p = ONE3 * trace(cS) + self.m_Presidual
        r = np.zeros(6)                       # <-- ManzariDafalias.cpp:1330-1332
        Kp = TWO3 * p * sd.h * dd_contr(sd.b, sd.n)
        t4 = (Kp + 2.0 * G * (sd.B - sd.C * trace(single_dot(sd.n, single_dot(sd.n, sd.n))))
              - K * sd.D * dd_contr(sd.n, r))
        if abs(t4) < SMALL:
            t4 = SMALL
        dg = (2.0 * G * dd_mixed(sd.n, ddev) - K * dvol * dd_contr(sd.n, r)) / t4
        dSig = (2.0 * G * to_contra(ddev) + K * dvol * I1 - macauley(dg)
                * (2.0 * G * (sd.B * sd.n - sd.C * (single_dot(sd.n, sd.n) - ONE3 * I1))
                   + K * sd.D * I1))
        dAl = macauley(dg) * TWO3 * sd.h * sd.b
        dFa = -macauley(dg) * self.m_cz * macauley(-sd.D) * (self.m_z_max * sd.n + cF)
        dPs = dg * to_cov(sd.R)
        return (cEE + (nE - cE) - dPs, cS + dSig, cA + dAl, cF + dFa, dg, vr)

    # ---- ModifiedEuler (:1365-1692) -- the deck default
    def modified_euler(self, cS, cE, cEE, cA, cF, ain, nE, G, K, vr_ignored):
        self.cnt.me_calls += 1
        T, dT, dT_min = 0.0, 1.0, 1e-6
        TolE = self.mTolR if self.mHonorTolRInME else 1e-4
        dStrain = nE - cE
        epsE = cEE + dStrain
        nStress, nAlpha, nFabric = cS.copy(), cA.copy(), cF.copy()
        dg = 0.0
        vr = self.m_e_init - (1 + self.m_e_init) * trace(nE)

        p = ONE3 * trace(nStress) + self.m_Presidual
        if p < self.m_Pmin + self.m_Presidual:
            self.cnt.me_lowp_clamp += 1
            nStress = dev(nStress) + self.m_Pmin * I1

        while T < 1.0:
            self.cnt.me_substeps += 1
            vr = self.m_e_init - (1 + self.m_e_init) * trace(cE + T * dStrain)
            dvol = dT * trace(dStrain)
            ddev = dev(dStrain) * dT

            # --- stage 1
            p = ONE3 * trace(nStress) + self.m_Presidual
            sd = self.state_dependent(nStress, nAlpha, nFabric, vr, ain)
            r = dev(nStress) / p
            Kp = TWO3 * p * sd.h * dd_contr(sd.b, sd.n)
            t4 = (Kp + 2.0 * G * (sd.B - sd.C * trace(single_dot(sd.n, single_dot(sd.n, sd.n))))
                  - K * sd.D * dd_contr(sd.n, r))
            if abs(t4) < SMALL:
                dS1 = np.zeros(6); dA1 = np.zeros(6); dF1 = np.zeros(6)
                dP1 = ddev + dvol * I1
            else:
                dg = (2.0 * G * dd_mixed(sd.n, ddev) - K * dvol * dd_contr(sd.n, r)) / t4
                if dg < -SMALL:
                    self.cnt.me_neg_dgamma += 1
                    dg = 0.0
                    dS1 = 2.0 * G * to_contra(ddev) + K * dvol * I1
                    dA1 = 3.0 * (dev(nStress + dS1) / trace(nStress + dS1)
                                 - dev(nStress) / trace(nStress))
                    dF1 = np.zeros(6); dP1 = np.zeros(6)
                else:
                    dS1 = (2.0 * G * to_contra(ddev) + K * dvol * I1 - macauley(dg)
                           * (2.0 * G * (sd.B * sd.n - sd.C * (single_dot(sd.n, sd.n) - ONE3 * I1))
                              + K * sd.D * I1))
                    dA1 = macauley(dg) * TWO3 * sd.h * sd.b
                    dF1 = (-macauley(dg) * self.m_cz * macauley(-sd.D)
                           * (self.m_z_max * sd.n + nFabric))
                    dP1 = dg * to_cov(sd.R)

            # --- stage 2
            s1 = nStress + dS1
            p = ONE3 * trace(s1) + self.m_Presidual
            if p < self.m_Presidual:
                if dT == dT_min:
                    # ManzariDafalias.cpp bare `return`: the caller keeps whatever
                    # the partially-completed substep loop already wrote.
                    self.cnt.me_abandon += 1
                    return epsE, nStress, nAlpha, nFabric, dg, vr
                dT = max(0.1 * dT, dT_min)
                continue
            sd = self.state_dependent(s1, nAlpha + dA1, nFabric + dF1, vr, ain)
            r = dev(s1) / p
            Kp = TWO3 * p * sd.h * dd_contr(sd.b, sd.n)
            t4 = (Kp + 2.0 * G * (sd.B - sd.C * trace(single_dot(sd.n, single_dot(sd.n, sd.n))))
                  - K * sd.D * dd_contr(sd.n, r))
            if abs(t4) < SMALL:
                dS2 = np.zeros(6); dA2 = np.zeros(6); dF2 = np.zeros(6)
                dP2 = ddev + dvol * I1
            else:
                dg = (2.0 * G * dd_mixed(sd.n, ddev) - K * dvol * dd_contr(sd.n, r)) / t4
                if dg < 0.0:
                    self.cnt.me_neg_dgamma += 1
                    dg = 0.0
                    dS2 = 2.0 * G * to_contra(ddev) + K * dvol * I1
                    dA2 = 3.0 * (dev(nStress + dS2) / trace(nStress + dS2)
                                 - dev(nStress) / trace(nStress))
                    dF2 = np.zeros(6); dP2 = np.zeros(6)
                else:
                    dS2 = (2.0 * G * to_contra(ddev) + K * dvol * I1 - macauley(dg)
                           * (2.0 * G * (sd.B * sd.n - sd.C * (single_dot(sd.n, sd.n) - ONE3 * I1))
                              + K * sd.D * I1))
                    dA2 = macauley(dg) * TWO3 * sd.h * sd.b
                    dF2 = (-macauley(dg) * self.m_cz * macauley(-sd.D)
                           * (self.m_z_max * sd.n + nFabric + dF1))
                    dP2 = dg * to_cov(sd.R)

            sN = nStress + 0.5 * (dS1 + dS2)
            fN = nFabric + 0.5 * (dF1 + dF2)
            aN = nAlpha + 0.5 * (dA1 + dA2)
            p = ONE3 * trace(sN) + self.m_Presidual
            if p < self.m_Presidual:
                if dT == dT_min:
                    # ManzariDafalias.cpp bare `return`: the caller keeps whatever
                    # the partially-completed substep loop already wrote.
                    self.cnt.me_abandon += 1
                    return epsE, nStress, nAlpha, nFabric, dg, vr
                dT = max(0.1 * dT, dT_min)
                continue

            sNorm = norm_contr(nStress)
            e_ = norm_contr(dS2 - dS1)
            err = e_ if sNorm < 0.5 else e_ / (2 * sNorm)

            if err > TolE:
                q = max(0.8 * math.sqrt(TolE / err), 0.1)
                if dT == dT_min:
                    # FORCE-ACCEPT (:1649-1662): the update always "succeeds".
                    self.cnt.me_force_accept += 1
                    epsE = epsE - 0.5 * (dP1 + dP2)
                    nStress = sN
                    eta = math.sqrt(13.5) * norm_contr(dev(nStress)) / trace(nStress)
                    if eta > self.m_Mc:
                        nStress = (ONE3 * trace(nStress) * I1
                                   + self.m_Mc / eta * dev(nStress))
                    nAlpha = cA + 3.0 * (dev(nStress) / trace(nStress)
                                         - dev(cS) / trace(cS))
                    T += dT
                dT = max(q * dT, dT_min)
            else:
                epsE = epsE - 0.5 * (dP1 + dP2)
                nStress, nAlpha, nFabric = sN, aN, fN
                nStress, nAlpha, nFabric, epsE = self.stress_correction(
                    cS, cEE, cA, cF, ain, nStress, nAlpha, nFabric, epsE, vr, G, K)
                T += dT
                # err == 0 gives inf in IEEE (the C++ path); guard Python's raise.
                q = (float("inf") if err == 0.0
                     else max(0.8 * math.sqrt(TolE / err), 0.5))
                dT = max(q * dT, dT_min)
                dT = min(dT, 1 - T)
        return epsE, nStress, nAlpha, nFabric, dg, vr

    # ---- Stress_Correction (:2562-2755) -- ACTIVE (mStressCorrectionInUse = true)
    def stress_correction(self, cS, cEE, cA, cF, ain, nS, nA, nF, epsE, vr, G, K):
        if not self.mStressCorrectionInUse:
            return nS, nA, nF, epsE
        maxIter = 50
        p = ONE3 * trace(nS) + self.m_Presidual
        if p < self.m_Pmin + self.m_Presidual:
            p = self.m_Pmin + self.m_Presidual
            return p * I1, np.zeros(6), nF, epsE
        fr = self.get_F(nS, nA)
        if abs(fr) < self.mTolF:
            return nS, nA, nF, epsE
        sN, aN = nS.copy(), nA.copy()
        for i in range(1, maxIter + 1):
            self.cnt.sc_iterations += 1
            devS = dev(sN)
            aC = stiffness(K, G)
            sd = self.state_dependent(sN, aN, nF, vr, ain)
            dSigmaP = aC @ to_cov(sd.R)
            aBar = TWO3 * sd.h * sd.b
            r = devS / p
            dfdS = sd.n - ONE3 * dd_contr(sd.n, r) * I1
            dfdA = -p * sd.n
            lam = fr / (dd_contr(dfdS, dSigmaP) - dd_contr(dfdA, aBar))
            if abs(self.get_F(sN - lam * dSigmaP, aN + lam * aBar)) < abs(fr):
                sN = sN - lam * dSigmaP
                aN = aN + lam * aBar
            else:
                lam = fr / dd_contr(dfdS, dfdS)
                if abs(self.get_F(sN - lam * dfdS, aN)) < abs(fr):
                    sN = sN - lam * dfdS
                else:
                    return nS, nA, nF, epsE
            fr = self.get_F(sN, aN)
            if abs(fr) < self.mTolF:
                nS, nA = sN, aN
                break
            if i == maxIter:
                if self.get_F(cS, nA) < self.mTolF:
                    dS = nS - cS
                    a_up, a_mid, a_dn = 1.0, 0.5, 0.0
                    fr_old = self.get_F(cS + a_mid * dS, nA)
                    for _ in range(maxIter):
                        if fr_old < 0.0:
                            a_dn = a_mid
                            a_mid = 0.5 * (a_up + a_mid)
                        else:
                            a_up = a_mid
                            a_mid = 0.5 * (a_dn + a_mid)
                        fr_old = self.get_F(cS + a_mid * dS, nA)
                        if abs(fr_old) < self.mTolF:
                            nS = cS + a_mid * dS
                            break
                else:
                    nS, nA, nF = cS.copy(), cA.copy(), cF.copy()
            p = ONE3 * trace(nS) + self.m_Presidual
        epsE = cEE + compliance(K, G) @ (nS - cS)
        return nS, nA, nF, epsE

    # ---- BackwardEuler_CPPM (:2189-2461)
    def _newton_res(self, x, inv):
        stress, alpha, fabric, dg = x[0:6], x[6:12], x[12:18], x[18]
        strain, cStrain, cStress, cEE, cA, cF = (inv["strain"], inv["cStrain"],
                                                 inv["cStress"], inv["cEE"],
                                                 inv["cA"], inv["cF"])
        vr, ain = inv["vr"], inv["ain"]
        trialE = cEE + (strain - cStrain)
        aD = compliance(self.mK, self.mG)
        sd = self.state_dependent(stress, alpha, fabric, vr, ain)
        aBar = TWO3 * sd.h * sd.b
        zBar = -self.m_cz * macauley(-sd.D) * (self.m_z_max * sd.n + fabric)
        eStrain = cEE + aD @ (stress - cStress)
        res = np.empty(19)
        res[0:6] = eStrain - trialE + dg * to_cov(sd.R)
        res[6:12] = alpha - cA - dg * aBar
        res[12:18] = fabric - cF - dg * zBar
        res[18] = self.get_F(stress, alpha)
        return res

    def _newton_iter2(self, x0, inv):
        """NewtonIter2 (:2876).  The C++ builds an analytical Jacobian by hand in
        NewtonSol; the oracle differentiates the SAME residual numerically. The
        residual, tolerance (`||R0||*TolR + TolR`) and iteration cap (30) are exact;
        the iterate path is not, so failure counts below are indicative, not exact."""
        x = np.array(x0, float)
        res = self._newton_res(x, inv)
        tol = np.linalg.norm(res) * self.mTolR + self.mTolR
        for _ in range(1, 31):
            if np.linalg.norm(res) < tol:
                return 1, x
            J = np.empty((19, 19))
            hstep = math.sqrt(2.0 * 2.220446049250313e-16)
            for i in range(19):
                xt = x.copy()
                t = xt[i]
                hh = hstep * max(1.0, abs(t))
                xt[i] = t + hh
                hh = xt[i] - t
                J[:, i] = (self._newton_res(xt, inv) - res) / hh
            try:
                dx = np.linalg.solve(J, -res)
            except np.linalg.LinAlgError:
                return -1, x
            if not np.all(np.isfinite(dx)):
                return -1, x
            x = x + dx
            res = self._newton_res(x, inv)
        return 0, x

    def _check(self, trial, stress, cA, nA):
        if dd_contr(self.normal_to_yield(stress, cA),
                    self.normal_to_yield(trial, cA)) < 0:
            return -4
        return 1

    def backward_euler(self, cS, cE, cEE, cA, cF, ain, nE, G, K, level):
        self.cnt.be_calls += 1
        if level > 10:
            self.cnt.be_maxsubstep += 1
            return (None, None, None, None, None, None, -3)
        cvr = self.m_e_init - (1 + self.m_e_init) * trace(cE)
        vr = self.m_e_init - (1 + self.m_e_init) * trace(nE)
        epsE = cEE + (nE - cE)
        nAlpha, nFabric, dg = cA.copy(), cF.copy(), 0.0
        self.mK, self.mG = self.elastic_moduli(cS, cvr)
        K, G = self.mK, self.mG
        aC = stiffness(K, G)
        trial = cS + aC @ (epsE - cEE)
        nStress = trial.copy()
        NextF = self.get_F(nStress, nAlpha)
        p = ONE3 * trace(nStress)
        inv = dict(strain=nE, cStrain=cE, cStress=cS, cEE=cEE, cA=cA, cF=cF,
                   vr=vr, ain=ain)

        if p < self.m_Pmin:
            # :2234-2286.  The tension-cutoff Newton is disabled by a literal
            # `errFlag = 0` (:2264), so this ALWAYS falls through to explicit.
            self.cnt.be_lowp_trial += 1
            self.cnt.be_explicit_fallback += 1
            r = self.explicit_integrator(cS, cE, cEE, cA, cF, ain, nE, G, K)
            return (*r, 1)

        if NextF > self.mTolF:
            errFlag, x = self._newton_iter2(
                np.concatenate([nStress, nAlpha, nFabric, [dg]]), inv)
            if errFlag == 1:
                nStress, nAlpha, nFabric, dg = x[0:6], x[6:12], x[12:18], x[18]
                errFlag = self._check(trial, nStress, cA, nAlpha)
            if errFlag != 1:
                self.cnt.be_newton_fail += 1
            SchemeControl = 2
            guard = 0
            while errFlag != 1:
                guard += 1
                if guard > 6:
                    break
                if errFlag == -1:
                    SchemeControl = 3
                if errFlag == -2:
                    SchemeControl = 2
                si = nE - cE
                if SchemeControl == 2:
                    self.cnt.be_recursive += 1
                    lvl = level + 1
                    r1 = self.backward_euler(cS, cE, cEE, cA, cF, ain,
                                             cE + si / 2, G, K, lvl)
                    if r1[-1] == -3:
                        SchemeControl += 1
                        continue
                    hE, hS, hA, hF, hg, hvr, ef = r1
                    r2 = self.backward_euler(hS, cE + si / 2, hE, hA, hF, ain,
                                             nE, G, K, lvl)
                    if r2[-1] == -3:
                        SchemeControl += 1
                        continue
                    if r2[-1] == 1:
                        _, nStress, nAlpha, nFabric, dg, _, _ = (None,) + r2[1:]
                        nStress, nAlpha, nFabric, dg = r2[1], r2[2], r2[3], r2[4]
                        errFlag = 1
                    else:
                        SchemeControl += 1
                        continue
                else:
                    self.cnt.be_explicit_fallback += 1
                    r = self.explicit_integrator(cS, cE, cEE, cA, cF, ain, nE, G, K)
                    epsE, nStress, nAlpha, nFabric, dg, vr = r
                    errFlag = 1
            sd = self.state_dependent(nStress, nAlpha, nFabric, vr, ain)
            # NOTE: this overwrites the elastic strain even when the explicit
            # fallback produced the stress, using that path's LAST-substep dGamma.
            epsE = cEE + (nE - cE) - dg * to_cov(sd.R)
        return (epsE, nStress, nAlpha, nFabric, dg, vr, 1)


# ---------------------------------------------------------------- IMPL-EX

class Implex:
    """IMPL-EX wrapper: extrapolated response every pass, true return at commit."""

    def __init__(self, mat: Sanisand, form="A", alpha=1.0):
        assert form in ("A", "B")
        self.mat = mat
        self.form = form
        self.alpha = alpha
        self.d_eps_p = np.zeros(6)
        self.dGamma_n = 0.0
        self.R_n = np.zeros(6)
        self.sig_tilde = mat.sig_n.copy()
        self.sig_implicit = mat.sig_n.copy()
        self.errors = []
        self._refresh_Rn()

    def _refresh_Rn(self):
        m = self.mat
        sd = m.state_dependent(m.sig_n, m.alpha_n, m.fabric_n,
                               m.void_ratio, m.alpha_in_n)
        self.R_n = sd.R

    def Ce(self):
        return stiffness(self.mat.mK, self.mat.mG)

    def extrapolate(self, eps_trial, f=None):
        """The response the element sees. Incremental form -- see the memo."""
        m = self.mat
        f = self.alpha if f is None else f
        deps = np.array(eps_trial, float) - m.eps_n
        if self.form == "A":
            dep = f * self.d_eps_p
        else:
            dep = (f * self.dGamma_n) * to_cov(self.R_n)
        self.sig_tilde = m.sig_n + self.Ce() @ (deps - dep)
        return self.sig_tilde

    def commit(self, eps_trial):
        m = self.mat
        epsp_old = m.eps_p_n()
        m.integrate(eps_trial)
        self.sig_implicit = m.sig.copy()
        m.commit()
        self.d_eps_p = m.eps_p_n() - epsp_old
        self.dGamma_n = m.dGamma_n
        self._refresh_Rn()
        den = norm_contr(self.sig_implicit) + m.m_P_atm * norm_contr(m.eps_n)
        self.errors.append(norm_contr(self.sig_tilde - self.sig_implicit) / den)


# ---------------------------------------------------------------- CSV / paths

def read_probe_csv(path):
    meta, rows = {}, []
    with open(path, newline="", encoding="utf-8") as fh:
        lines = fh.readlines()
    body = []
    for ln in lines:
        if ln.startswith("#"):
            k, _, v = ln[2:].strip().partition(": ")
            meta[k] = v
        else:
            body.append(ln)
    rd = csv.DictReader(body)
    for r in rd:
        rows.append({k: float(v) for k, v in r.items()})
    return meta, rows


def committed_from_row(row):
    """Row -> the material's committed state (internal sign)."""
    eps = -np.array([row[f"eps{i}"] for i in range(6)])
    sig = -np.array([row[f"sig{i}"] for i in range(6)])
    st = [row[f"state{i}"] for i in range(26)]
    epsE = np.array(st[0:6])
    alpha = np.array(st[6:12])
    fabric = np.array(st[12:18])
    alpha_in = np.array(st[18:24])
    return eps, sig, epsE, alpha, fabric, alpha_in, st[24], st[25]


# ---------------------------------------------------------------- drivers

def make_material(consts=None, scheme=1, **kw):
    return Sanisand(consts=consts, scheme=scheme, **kw)


def seed_from_csv(mat, row):
    eps, sig, epsE, alpha, fabric, alpha_in, e, dg = committed_from_row(row)
    mat.set_committed(eps, sig, epsE, alpha, fabric, alpha_in, dg)
    return e


def solve_lateral(step_fn, p0, de_zz, eps_n, tol=1e-11, itmax=80):
    """Mixed control: prescribe d(eps_zz), find d(eps_xx)=d(eps_yy) with sig_xx = p0.

    `step_fn(eps_trial) -> sigma`. Secant iteration; returns (eps_trial, sigma, its)."""
    def make(dl):
        e = eps_n.copy()
        e[0] += dl
        e[1] += dl
        e[2] += de_zz
        return e

    x0 = 0.0
    s0 = step_fn(make(x0))
    g0 = s0[0] - p0
    x1 = -0.3 * de_zz if de_zz != 0 else 1e-9
    s1 = step_fn(make(x1))
    g1 = s1[0] - p0
    for it in range(itmax):
        if abs(g1) < tol * max(1.0, abs(p0)):
            return make(x1), s1, it
        if g1 == g0:
            break
        x2 = x1 - g1 * (x1 - x0) / (g1 - g0)
        x0, g0 = x1, g1
        x1 = x2
        s1 = step_fn(make(x1))
        g1 = s1[0] - p0
    return make(x1), s1, itmax


@dataclass
class Run:
    ez: list = field(default_factory=list)
    sig: list = field(default_factory=list)
    e: list = field(default_factory=list)
    eta: list = field(default_factory=list)
    p: list = field(default_factory=list)
    q: list = field(default_factory=list)
    implex_err: list = field(default_factory=list)
    cnt: Counters = field(default_factory=Counters)
    abandoned: int = 0


def _invariants(sig):
    p = ONE3 * trace(sig)
    s = dev(sig)
    q = math.sqrt(1.5 * dd_contr(s, s))
    return p, q


def drive_triaxial(seed_row, p0, nstep, ez_max, kind="implicit", scheme=1,
                   companion_scheme=None, consts=None, sample=None, **kw):
    """Drained triaxial compression at one Gauss point.

    kind: 'implicit' | 'implex_A' | 'implex_B'.  'reference' = 'implicit' at a tiny
    increment (the caller chooses nstep)."""
    mat = make_material(consts, scheme=scheme, **kw)
    seed_from_csv(mat, seed_row)
    eps0 = mat.eps_n.copy()
    run = Run()
    ix = None
    if kind.startswith("implex"):
        cs = companion_scheme if companion_scheme is not None else scheme
        mat.scheme = cs
        ix = Implex(mat, form=kind[-1])
    de = ez_max / nstep
    sample = set(sample) if sample else None

    def record(k):
        sig = mat.sig_n if ix is None else ix.sig_tilde
        p, q = _invariants(sig)
        run.ez.append(mat.eps_n[2] - eps0[2])
        run.sig.append(sig.copy())
        run.e.append(mat.m_e_init - (1 + mat.m_e_init) * trace(mat.eps_n))
        run.p.append(p)
        run.q.append(q)
        run.eta.append(q / p if p != 0 else float("nan"))

    record(0)
    for k in range(1, nstep + 1):
        if ix is None:
            def step_fn(e):
                try:
                    return mat.integrate(e)
                except Abandoned:
                    run.abandoned += 1
                    return mat.sig.copy()
            eps_t, sig, _ = solve_lateral(step_fn, p0, de, mat.eps_n)
            try:
                mat.integrate(eps_t)
            except Abandoned:
                run.abandoned += 1
            mat.commit()
        else:
            eps_t, sig, _ = solve_lateral(ix.extrapolate, p0, de, mat.eps_n)
            ix.extrapolate(eps_t)
            try:
                ix.commit(eps_t)
            except Abandoned:
                run.abandoned += 1
                mat.commit()
        if sample is None or k in sample:
            record(k)
    run.cnt = mat.cnt
    if ix is not None:
        run.implex_err = ix.errors
    return run, mat, ix


def drive_prescribed(seed_row, path_fn, nstep, kind="implicit", scheme=1,
                     companion_scheme=None, consts=None, **kw):
    """Prescribed full 6-component strain path: eps(t), t in [0,1]."""
    mat = make_material(consts, scheme=scheme, **kw)
    seed_from_csv(mat, seed_row)
    run = Run()
    ix = None
    if kind.startswith("implex"):
        cs = companion_scheme if companion_scheme is not None else scheme
        mat.scheme = cs
        ix = Implex(mat, form=kind[-1])
    eps0 = mat.eps_n.copy()

    def record():
        sig = mat.sig_n if ix is None else ix.sig_tilde
        p, q = _invariants(sig)
        run.ez.append(mat.eps_n[2] - eps0[2])
        run.sig.append(sig.copy())
        run.e.append(mat.m_e_init - (1 + mat.m_e_init) * trace(mat.eps_n))
        run.p.append(p)
        run.q.append(q)
        run.eta.append(q / p if p != 0 else float("nan"))

    record()
    for k in range(1, nstep + 1):
        eps_t = eps0 + path_fn(k / nstep)
        if ix is None:
            try:
                mat.integrate(eps_t)
            except Abandoned:
                run.abandoned += 1
            mat.commit()
        else:
            ix.extrapolate(eps_t)
            try:
                ix.commit(eps_t)
            except Abandoned:
                run.abandoned += 1
                mat.commit()
        record()
    run.cnt = mat.cnt
    if ix is not None:
        run.implex_err = ix.errors
    return run, mat, ix


# ---------------------------------------------------------------- gates

def _fmt(x, w=12, p=4):
    return f"{x:>{w}.{p}e}"


def _mat_from_meta(meta, scheme=None):
    consts = dict(CONSTS)
    consts["e_init"] = float(meta["e_init"])
    return make_material(consts,
                         scheme=int(meta["scheme"]) if scheme is None else scheme,
                         Pmin=float(meta["Pmin"]),
                         Presidual=float(meta["Presidual"]),
                         TolF=float(meta["tol"]), TolR=float(meta["tol"]),
                         honor_tolR=bool(int(meta["honorTolR"])))


def _replay(meta, rows, seed_perturb=0.0):
    """Replay the recorded strain sequence; return the per-step committed states."""
    mat = _mat_from_meta(meta, rows[0])
    mat = _mat_from_meta(meta)
    seed_from_csv(mat, rows[0])
    if seed_perturb:
        mat.sig_n[2] *= (1.0 + seed_perturb)
        mat._refresh_moduli()
    out = []
    for r in rows[1:]:
        eps = committed_from_row(r)[0]
        mat.integrate(eps)
        mat.commit()
        out.append((mat.sig_n.copy(), mat.void_ratio, mat.alpha_n.copy(),
                    mat.fabric_n.copy(), mat.epsE_n.copy()))
    return mat, out


def _conditioning_twin(meta, rows, perturb=1e-15):
    """CONTROL A: how far does the ORACLE diverge from ITSELF when the seed is moved
    by `perturb` relative?  This bounds what any transcription can reproduce."""
    _, a = _replay(meta, rows)
    _, b = _replay(meta, rows, seed_perturb=perturb)
    worst = 0.0
    for x, y in zip(a, b):
        worst = max(worst,
                    norm_contr(x[0] - y[0]) / norm_contr(x[0]),
                    abs(x[1] - y[1]) / abs(x[1]),
                    norm_contr(x[2] - y[2]) / max(norm_contr(x[2]), 1e-30),
                    (norm_contr(x[3] - y[3]) / norm_contr(x[3])
                     if norm_contr(x[3]) > 1e-30 else 0.0),
                    norm_contr(x[4] - y[4]) / norm_contr(x[4]))
    return worst


def _residual_audit(meta, rows):
    """CONTROL B (scheme 2 only): evaluate the CODE'S OWN 19-residual at the binary's
    answer and at the oracle's answer, from the same committed state.  If the binary's
    residual is no smaller than the oracle's, both are converged to the same declared
    tolerance and the difference is Newton-stopping ambiguity, not a formula error."""
    mat = _mat_from_meta(meta)
    seed_from_csv(mat, rows[0])
    wb = wo = 0.0
    for r in rows[1:]:
        eps, sig, epsE, alpha, fab, ain, e_, dg = committed_from_row(r)
        cS, cE, cEE, cA, cF = (mat.sig_n, mat.eps_n, mat.epsE_n,
                               mat.alpha_n, mat.fabric_n)
        Ce = stiffness(mat.mK, mat.mG)
        td = Ce @ (eps - cE)
        ain_use = (cA.copy() if dd_contr(cA - mat.alpha_in_n, td) < 0
                   else mat.alpha_in_n.copy())
        vr = mat.m_e_init - (1 + mat.m_e_init) * trace(eps)
        cvr = mat.m_e_init - (1 + mat.m_e_init) * trace(cE)
        mat.mK, mat.mG = mat.elastic_moduli(cS, cvr)
        inv = dict(strain=eps, cStrain=cE, cStress=cS, cEE=cEE, cA=cA, cF=cF,
                   vr=vr, ain=ain_use)
        wb = max(wb, float(np.linalg.norm(
            mat._newton_res(np.concatenate([sig, alpha, fab, [dg]]), inv))))
        mat.integrate(eps)
        wo = max(wo, float(np.linalg.norm(mat._newton_res(
            np.concatenate([mat.sig, mat.alpha, mat.fabric, [mat.dGamma]]), inv))))
        mat.commit()
    return wb, wo


def gate_G0(verbose=True):
    """Replay the binary's recorded strain sequence through `implicit`."""
    print("=" * 78)
    print("G0 -- oracle `implicit` vs the fork binary on the binary's own strain path")
    print("=" * 78)
    files = sorted(glob.glob(os.path.join(DATA, "tx_*.csv")))
    ok = True
    table = []
    for f in files:
        meta, rows = read_probe_csv(f)
        if len(rows) < 3:
            print(f"  {os.path.basename(f)}: only {len(rows)} rows -- SKIPPED "
                  f"(the binary's own global solve stalled)")
            continue
        mat, states = _replay(meta, rows)
        worst = dict(sig=0.0, e=0.0, alpha=0.0, z=0.0, epsE=0.0, step=0)
        for r, st in zip(rows[1:], states):
            _, sig, epsE, alpha, fabric, _, e_ref, _ = committed_from_row(r)
            oS, oe, oA, oZ, oE = st
            ds = norm_contr(oS - sig) / norm_contr(sig)
            de = abs(oe - e_ref) / abs(e_ref)
            da = norm_contr(oA - alpha) / max(norm_contr(alpha), 1e-30)
            nz = norm_contr(fabric)
            dz = (norm_contr(oZ - fabric) / nz) if nz > 1e-30 else norm_contr(oZ - fabric)
            dE = norm_contr(oE - epsE) / norm_contr(epsE)
            if ds > worst["sig"]:
                worst.update(sig=ds, step=int(r["step"]))
            worst["e"] = max(worst["e"], de)
            worst["alpha"] = max(worst["alpha"], da)
            worst["z"] = max(worst["z"], dz)
            worst["epsE"] = max(worst["epsE"], dE)
        raw = max(worst["sig"], worst["e"], worst["alpha"], worst["z"], worst["epsE"])
        table.append((os.path.basename(f), len(rows) - 1, worst, raw,
                      int(meta["scheme"]), meta, rows))
    print(f"{'file':<32}{'n':>4}{'sigma':>13}{'e':>13}{'alpha':>13}{'z':>13}"
          f"{'epsE':>13}  raw <= 1e-8")
    for name, n, w, raw, sch, _, _ in table:
        print(f"{name:<32}{n:>4}{_fmt(w['sig'])}{_fmt(w['e'])}{_fmt(w['alpha'])}"
              f"{_fmt(w['z'])}{_fmt(w['epsE'])}  {'yes' if raw <= 1e-8 else 'NO'}")

    print("\n  Controls for every row that misses 1e-8 raw:")
    verdicts = []
    for name, n, w, raw, sch, meta, rows in table:
        if raw <= 1e-8:
            verdicts.append((name, "PASS", "raw <= 1e-8"))
            continue
        if sch == 2:
            wb, wo = _residual_audit(meta, rows)
            good = wo <= wb
            print(f"    {name}: CONTROL B (scheme 2, FD Jacobian). worst "
                  f"||R|| at the BINARY's answer = {wb:.3e}, at the ORACLE's = "
                  f"{wo:.3e}  -> the oracle satisfies the code's own residual "
                  f"{'at least as well' if good else 'WORSE'}")
            verdicts.append((name, "PASS" if good else "FAIL",
                             f"Newton-stopping ambiguity: ||R||_bin={wb:.2e} "
                             f">= ||R||_ora={wo:.2e}"))
        else:
            cond = _conditioning_twin(meta, rows)
            good = raw <= cond
            print(f"    {name}: CONTROL A (path conditioning). A 1e-15 relative seed "
                  f"perturbation moves the ORACLE from ITSELF by {cond:.3e}; the "
                  f"oracle-vs-binary gap is {raw:.3e}  -> "
                  f"{'within' if good else 'ABOVE'} the path's own round-off "
                  f"amplification")
            verdicts.append((name, "PASS" if good else "FAIL",
                             f"conditioning bound {cond:.2e} >= gap {raw:.2e}"))
    ok = all(v[1] == "PASS" for v in verdicts)
    print(f"\nG0 verdict: {'PASS' if ok else 'FAIL'}")
    for nm, v, why in verdicts:
        print(f"    {nm:<32} {v}   ({why})")
    return ok, table, verdicts


def _err_vs_ref(run, ref, key="sig"):
    """Error at the common terminal state."""
    if key == "sig":
        a, b = run.sig[-1], ref.sig[-1]
        return norm_contr(a - b) / norm_contr(b)
    a = getattr(run, key)[-1]
    b = getattr(ref, key)[-1]
    return abs(a - b) / abs(b)


def _paths():
    """The seed rows available, keyed by label."""
    out = {}
    for f in sorted(glob.glob(os.path.join(DATA, "tx_*.csv"))):
        meta, rows = read_probe_csv(f)
        if not rows:
            continue
        out[os.path.basename(f)] = (meta, rows[0])
    return out


def _seed(label_sub):
    for k, v in _paths().items():
        if label_sub in k:
            return v
    raise SystemExit(f"no probe CSV matching {label_sub!r} in {DATA}")


def gate_G1(ez_max=0.004, ns=(25, 50, 100, 200, 400), ref_mult=16):
    """Two errors are reported and they answer different questions.

      TOTAL      = || X(N) - reference || / || reference ||
                   -- everything wrong with the answer at increment N.
      EXTRAP     = || X(N) - implicit(N) || / || implicit(N) ||
                   -- the IMPL-EX perturbation ALONE, at the same increment.
                   THIS is the quantity ADR sec.2 claims is first order in dt.
    """
    print("=" * 78)
    print("G1 -- convergence order under increment halving")
    print("=" * 78)
    results = {}
    for label, sub, p0 in (("T1 p0=100", "tx_p100_e0.6944_s1_n40", 100.0),
                           ("T2 p0=5", "tx_p5_", 5.0)):
        meta, row = _seed(sub)
        consts = dict(CONSTS); consts["e_init"] = float(meta["e_init"])
        refn = ns[-1] * ref_mult
        ref, _, _ = drive_triaxial(row, p0, refn, ez_max, "implicit", consts=consts)
        print(f"\n{label}   ez_max={ez_max}  reference = implicit at "
              f"d(eps_z)={ez_max/refn:.3e} ({refn} steps)")
        print(f"{'N':>5}{'d eps_z':>12} | {'TOTAL A':>11}{'TOT B(s1)':>11}"
              f"{'TOT B(s2)':>11}{'TOT impl':>11} | {'EXTRAP A':>11}"
              f"{'EXT B(s1)':>11}{'EXT B(s2)':>11} | {'ieA max':>12}"
              f"{'ieA mean':>12}{'ieA last':>12}")
        rows_ = []
        for n in ns:
            de = ez_max / n
            rA, _, ixA = drive_triaxial(row, p0, n, ez_max, "implex_A", consts=consts)
            rB1, _, ixB1 = drive_triaxial(row, p0, n, ez_max, "implex_B", consts=consts)
            rB2, _, ixB2 = drive_triaxial(row, p0, n, ez_max, "implex_B",
                                          companion_scheme=2, consts=consts)
            rI, _, _ = drive_triaxial(row, p0, n, ez_max, "implicit", consts=consts)
            tot = (_err_vs_ref(rA, ref), _err_vs_ref(rB1, ref),
                   _err_vs_ref(rB2, ref), _err_vs_ref(rI, ref))
            ext = (_err_vs_ref(rA, rI), _err_vs_ref(rB1, rI), _err_vs_ref(rB2, rI))
            ie = max(ixA.errors)
            iem = float(np.mean(ixA.errors[1:]))
            iel = ixA.errors[-1]
            rows_.append((n, de) + tot + ext + (ie, iem, iel))
            print(f"{n:>5}{de:>12.3e} | {tot[0]:>11.4e}{tot[1]:>11.4e}"
                  f"{tot[2]:>11.4e}{tot[3]:>11.4e} | {ext[0]:>11.4e}"
                  f"{ext[1]:>11.4e}{ext[2]:>11.4e} | {ie:>12.4e}"
                  f"{iem:>12.4e}{iel:>12.4e}")
        h = np.array([r[1] for r in rows_])
        slopes = {}
        names = ["TOT A", "TOT B1", "TOT B2", "TOT impl",
                 "EXT A", "EXT B1", "EXT B2",
                 "implexErr A max", "implexErr A mean", "implexErr A last"]
        line = f"{'slope':>5}{'':>12} | "
        for j, nm in enumerate(names, start=2):
            y = np.array([r[j] for r in rows_])
            m_ = float(np.polyfit(np.log(h), np.log(np.maximum(y, 1e-300)), 1)[0])
            slopes[nm] = m_
            line += f"{m_:>11.3f}"
            if j in (5, 8):
                line += " | "
        print(line)
        # slope over the finest three increments only (asymptotic regime)
        print(f"{'  (fine-3)':>17} | ", end="")
        for j, nm in enumerate(names, start=2):
            y = np.array([r[j] for r in rows_[-3:]])
            m_ = float(np.polyfit(np.log(h[-3:]), np.log(np.maximum(y, 1e-300)), 1)[0])
            slopes[nm + "_fine"] = m_
            print(f"{m_:>11.3f}", end=" | " if j in (5, 8) else "")
        print()
        results[label] = (rows_, slopes)
    return results


def gate_G2(ez_max=0.02, ref_de=1.0e-5):
    """The error as a function of the increment, spanning the campaign's range.

    `1e-4 m / h` is the NOMINAL element strain of the act's platen step (1.0e-4 at
    h = 1.0 m, 2.0e-4 at h = 0.5 m).  README section 1: that is a LOWER bound -- the
    Gauss point at the footing corner sees a strain concentration of one to two orders
    more -- so the increment is swept, not quoted at a single value.

    TOTAL  = vs the implicit reference at d(eps_z) = ref_de.
    EXTRAP = vs the implicit run at the SAME increment: the IMPL-EX perturbation alone.
    """
    print("=" * 78)
    print("G2 -- error vs increment, over the campaign's range")
    print("=" * 78)
    des = [1.0e-4, 2.0e-4, 5.0e-4, 1.0e-3, 2.0e-3, 5.0e-3, 1.0e-2]
    tags = {1.0e-4: "  <- h=1.0 m nominal", 2.0e-4: "  <- h=0.5 m nominal"}
    out = {}
    for label, sub, p0 in (("T1 p0=100", "tx_p100_e0.6944_s1_n40", 100.0),
                           ("T2 p0=5", "tx_p5_", 5.0),
                           ("U dense e=0.60", "tx_p100_e0.6_", 100.0)):
        meta, row = _seed(sub)
        consts = dict(CONSTS); consts["e_init"] = float(meta["e_init"])
        refn = int(round(ez_max / ref_de))
        ref, _, _ = drive_triaxial(row, p0, refn, ez_max, "implicit", consts=consts)
        print(f"\n  {label}   to eps_z = {ez_max}   reference = implicit at "
              f"d(eps_z) = {ref_de:.1e} ({refn} steps); "
              f"reference terminal q/p = {ref.eta[-1]:.5f}")
        print(f"    {'d eps_z':>9}{'N':>5} | {'TOT impl':>10}{'TOT A':>10}"
              f"{'TOT B':>10} | {'EXTRAP A':>10}{'EXTRAP B':>10} | "
              f"{'de/e A':>9}{'deta/eta A':>11}{'ieA mean':>10} | {'q/p A':>8}"
              f"{'q/p impl':>9}  FA")
        for de in des:
            n = int(round(ez_max / de))
            rI, _, _ = drive_triaxial(row, p0, n, ez_max, "implicit", consts=consts)
            rA, _, ixA = drive_triaxial(row, p0, n, ez_max, "implex_A", consts=consts)
            rB, _, ixB = drive_triaxial(row, p0, n, ez_max, "implex_B",
                                        companion_scheme=1, consts=consts)
            row_ = dict(
                de=de, n=n,
                tot_I=_err_vs_ref(rI, ref), tot_A=_err_vs_ref(rA, ref),
                tot_B=_err_vs_ref(rB, ref),
                ext_A=_err_vs_ref(rA, rI), ext_B=_err_vs_ref(rB, rI),
                de_e=_err_vs_ref(rA, ref, "e"), deta=_err_vs_ref(rA, ref, "eta"),
                ie=float(np.mean(ixA.errors[1:])) if len(ixA.errors) > 1 else
                float(ixA.errors[0]),
                eta_A=rA.eta[-1], eta_I=rI.eta[-1],
                fa=rI.cnt.me_force_accept)
            out[(label, de)] = row_
            print(f"    {de:>9.1e}{n:>5} | {row_['tot_I']:>10.3e}"
                  f"{row_['tot_A']:>10.3e}{row_['tot_B']:>10.3e} | "
                  f"{row_['ext_A']:>10.3e}{row_['ext_B']:>10.3e} | "
                  f"{row_['de_e']:>9.2e}{row_['deta']:>11.3e}{row_['ie']:>10.3e}"
                  f" | {row_['eta_A']:>8.4f}{row_['eta_I']:>9.4f}"
                  f"{row_['fa']:>4}" + tags.get(de, ""))
    return out


def gate_G3():
    print("=" * 78)
    print("G3 -- the corner: a path driven onto the Pmin = 0.0101 kPa floor")
    print("=" * 78)
    meta, row = _seed("tx_p5_")
    consts = dict(CONSTS); consts["e_init"] = float(meta["e_init"])

    # volumetric extension (internal sign: negative volumetric strain) with a
    # superposed shear, so p collapses onto the floor and the point keeps flowing.
    def path(t):
        e = np.zeros(6)
        ev = -3.0e-3 * min(t / 0.5, 1.0)          # unload p to the floor by t=0.5
        e[0] = e[1] = e[2] = ev / 3.0
        e[2] += 6.0e-3 * t                        # keep shearing throughout
        e[0] -= 3.0e-3 * t
        e[1] -= 3.0e-3 * t
        return e

    out = {}
    refn = 160 * 64
    ref, _, _ = drive_prescribed(row, path, refn, "implicit", consts=consts)
    print(f"  reference = implicit at {refn} prescribed steps")
    for n in (40, 80, 160):
        rI, mI, _ = drive_prescribed(row, path, n, "implicit", consts=consts)
        rA, mA, ixA = drive_prescribed(row, path, n, "implex_A", consts=consts)
        rB, mB, ixB = drive_prescribed(row, path, n, "implex_B",
                                       companion_scheme=2, consts=consts)
        rS2, mS2, _ = drive_prescribed(row, path, n, "implicit", scheme=2,
                                       consts=consts)
        print(f"\n  n = {n}   terminal p: ref={ref.p[-1]:.5f}  implicit={rI.p[-1]:.5f}"
              f"  A={rA.p[-1]:.5f}  B={rB.p[-1]:.5f}  scheme2={rS2.p[-1]:.5f} kPa")
        print(f"    min p along path: ref={min(ref.p):.5f}  implicit={min(rI.p):.5f}"
              f"  A={min(rA.p):.5f}  B={min(rB.p):.5f}")
        for nm, r in (("implicit", rI), ("implex_A", rA), ("implex_B", rB)):
            es = _err_vs_ref(r, ref, "sig")
            bounded = np.all(np.isfinite(np.array(r.sig)))
            print(f"    {nm:<10} |dsig|/|sig|={es:.4e}  finite={bounded}"
                  f"  abandoned={r.abandoned}")
        print(f"    ModifiedEuler force-accept at dT_min: implicit={rI.cnt.me_force_accept}"
              f" / {rI.cnt.me_substeps} substeps"
              f"   low-p clamps={rI.cnt.me_lowp_clamp}"
              f"   abandons={rI.cnt.me_abandon}")
        print(f"    scheme-2 companion: calls={rS2.cnt.be_calls}"
              f"  low-p trial branch={rS2.cnt.be_lowp_trial}"
              f"  newton failures={rS2.cnt.be_newton_fail}"
              f"  recursive halvings={rS2.cnt.be_recursive}"
              f"  EXPLICIT fallbacks={rS2.cnt.be_explicit_fallback}")
        out[n] = dict(ref=ref, I=rI, A=rA, B=rB, S2=rS2)
    return out


def gate_G4():
    print("=" * 78)
    print("G4 -- implex_A (plastic strain) vs implex_B (dGamma), head to head")
    print("=" * 78)
    out = {}
    for label, sub, p0 in (("T1 p0=100", "tx_p100_e0.6944_s1_n40", 100.0),
                           ("T2 p0=5", "tx_p5_", 5.0),
                           ("U dense e=0.60", "tx_p100_e0.6_", 100.0)):
        meta, row = _seed(sub)
        consts = dict(CONSTS); consts["e_init"] = float(meta["e_init"])
        ez_max, n = 0.004, 40
        ref, _, _ = drive_triaxial(row, p0, 50 * n, ez_max, "implicit", consts=consts)
        variants = {
            "A  (eps_p, companion s1)": dict(kind="implex_A", companion_scheme=1),
            "A  (eps_p, companion s2)": dict(kind="implex_A", companion_scheme=2),
            "B  (dGamma, companion s1)": dict(kind="implex_B", companion_scheme=1),
            "B  (dGamma, companion s2)": dict(kind="implex_B", companion_scheme=2),
        }
        print(f"\n  {label}  d(eps_z)={ez_max/n:.1e}  ({n} steps to eps_z={ez_max})")
        for nm, kw in variants.items():
            r, m, ix = drive_triaxial(row, p0, n, ez_max, consts=consts, **kw)
            es = _err_vs_ref(r, ref, "sig")
            et = _err_vs_ref(r, ref, "eta")
            print(f"    {nm:<28} |dsig|/|sig|={es:.4e}  deta/eta={et:.4e}"
                  f"  implexError max={max(ix.errors):.3e}")
            out[(label, nm)] = (es, et, max(ix.errors))
    # T4: compression then reversal
    print("\n  T4 -- compression then reversal (p0=100, +-0.004 in eps_z)")
    meta, row = _seed("tx_p100_e0.6944_s1_n40")
    consts = dict(CONSTS); consts["e_init"] = float(meta["e_init"])

    def path(t):
        e = np.zeros(6)
        ez = 0.004 * (2 * t if t <= 0.5 else 2 - 2 * t)
        e[2] = ez
        e[0] = e[1] = -0.45 * ez
        return e

    ref, _, _ = drive_prescribed(row, path, 80 * 50, "implicit", consts=consts)
    for nm, kw in (("implicit", dict(kind="implicit")),
                   ("implex_A", dict(kind="implex_A", companion_scheme=1)),
                   ("implex_B", dict(kind="implex_B", companion_scheme=2))):
        r, m, ix = drive_prescribed(row, path, 80, consts=consts, **kw)
        es = _err_vs_ref(r, ref, "sig")
        extra = f"  implexError max={max(ix.errors):.3e}" if ix else ""
        print(f"    {nm:<10} |dsig|/|sig|={es:.4e}{extra}")
        out[("T4", nm)] = es
    return out


def gate_G5():
    print("=" * 78)
    print("G5 -- tangent identity  d(sigma~)/d(eps) == Ce(p_n)")
    print("=" * 78)
    meta, row = _seed("tx_p100_e0.6944_s1_n40")
    consts = dict(CONSTS); consts["e_init"] = float(meta["e_init"])
    worst = 0.0
    for form in ("A", "B"):
        mat = make_material(consts, scheme=1)
        seed_from_csv(mat, row)
        ix = Implex(mat, form=form)
        # take a few real steps first so d_eps_p / dGamma are non-trivial
        eps = mat.eps_n.copy()
        for _ in range(5):
            eps = mat.eps_n.copy()
            eps[2] += 1e-4
            eps[0] -= 0.45e-4
            eps[1] -= 0.45e-4
            ix.extrapolate(eps)
            ix.commit(eps)
        Ce = ix.Ce()
        base = ix.extrapolate(eps).copy()
        num = np.empty((6, 6))
        for j in range(6):
            hh = 1e-9
            e2 = eps.copy()
            e2[j] += hh
            num[:, j] = (ix.extrapolate(e2) - base) / hh
        ix.extrapolate(eps)
        rel = np.max(np.abs(num - Ce)) / np.max(np.abs(Ce))
        worst = max(worst, rel)
        print(f"  implex_{form}: max|d(sig~)/d(eps) - Ce| / max|Ce| = {rel:.3e}")
    print(f"\nG5 verdict: {'PASS' if worst < 1e-6 else 'FAIL'}  "
          f"(worst {worst:.3e}; FD step 1e-9)")
    return worst


def control_total_vs_incremental():
    """ADR sec.2 writes sigma~ = Ce(p_n):(eps_{n+1} - eps_p~) in TOTAL form.  With a
    pressure-dependent Ce that does not reproduce sigma_n at a zero increment.
    Quantify the defect."""
    print("=" * 78)
    print("CONTROL -- ADR sec.2's TOTAL-strain form vs the incremental form")
    print("=" * 78)
    for sub, p0 in (("tx_p100_e0.6944_s1_n40", 100.0), ("tx_p5_", 5.0)):
        meta, row = _seed(sub)
        consts = dict(CONSTS); consts["e_init"] = float(meta["e_init"])
        mat = make_material(consts, scheme=1)
        seed_from_csv(mat, row)
        Ce = stiffness(mat.mK, mat.mG)
        eps_p = mat.eps_p_n()
        sig_total = Ce @ (mat.eps_n - eps_p)      # zero increment, total form
        rel = norm_contr(sig_total - mat.sig_n) / norm_contr(mat.sig_n)
        print(f"  p0={p0:>5.0f} kPa: at a ZERO increment the total form returns "
              f"p={ONE3*trace(sig_total):.4f} kPa against the committed "
              f"p={ONE3*trace(mat.sig_n):.4f} kPa  (rel err {rel:.3e})")
    print("  => the oracle uses the incremental form "
          "sigma~ = sigma_n + Ce(p_n):(d eps - f*d eps_p).")


def control_forward_euler_r_bug():
    """The zero-`r` defect is reproduced only on the ForwardEuler path (scheme 5)."""
    print("=" * 78)
    print("CONTROL -- the ForwardEuler zero-`r` defect (scheme 5 only)")
    print("=" * 78)
    meta, row = _seed("tx_p100_e0.6944_s1_n40")
    consts = dict(CONSTS); consts["e_init"] = float(meta["e_init"])
    for scheme, nm in ((1, "ModifiedEuler (r correct)"), (5, "ForwardEuler (r == 0)")):
        r, m, _ = drive_triaxial(row, 100.0, 40, 0.004, "implicit",
                                 scheme=scheme, consts=consts)
        print(f"  scheme {scheme} {nm:<28} terminal q/p = {r.eta[-1]:.6f}  "
              f"p = {r.p[-1]:.4f} kPa")


def gate_instr():
    """The two instrumented counts the ADR asked for, on their own.

    (1) ModifiedEuler's FORCE-ACCEPT at dT_min -- how often, and what it does.
    (2) BackwardEuler_CPPM's silent explicit fallback -- how often the retry ladder
        (and the low-p trial branch, which is an UNCONDITIONAL fallback) reaches it.
    """
    print("=" * 78)
    print("INSTRUMENTATION -- ModifiedEuler force-accept and BackwardEuler fallback")
    print("=" * 78)
    out = {}
    for label, sub, p0 in (("T1 p0=100", "tx_p100_e0.6944_s1_n40", 100.0),
                           ("T2 p0=5", "tx_p5_", 5.0),
                           ("T2b p0=1", "tx_p1_e0.6944_s1_n400", 1.0)):
        meta, row = _seed(sub)
        consts = dict(CONSTS); consts["e_init"] = float(meta["e_init"])
        print(f"\n  {label} -- prescribed triaxial-shaped strain ramp to eps_z = 0.02,"
              f" one step of size d(eps_z)")

        def path(t, EZ=0.02):
            e = np.zeros(6)
            e[2] = EZ * t
            e[0] = e[1] = -0.45 * EZ * t
            return e

        print(f"    {'d eps_z':>10}{'steps':>7}{'ME substeps':>13}"
              f"{'FORCE-ACCEPT':>14}{'abandon':>9}{'lowp clamp':>12}"
              f"{'terminal q/p':>14}")
        for n in (200, 100, 40, 20, 10, 5, 2, 1):
            r, m, _ = drive_prescribed(row, path, n, "implicit", scheme=1,
                                       consts=consts)
            print(f"    {0.02/n:>10.2e}{n:>7}{m.cnt.me_substeps:>13}"
                  f"{m.cnt.me_force_accept:>14}{m.cnt.me_abandon:>9}"
                  f"{m.cnt.me_lowp_clamp:>12}{r.eta[-1]:>14.5f}")
            out[(label, "s1", n)] = m.cnt
        print(f"    {'d eps_z':>10}{'steps':>7}{'BE calls':>10}{'low-p trial':>12}"
              f"{'newton fail':>12}{'recursive':>11}{'EXPLICIT fb':>13}"
              f"{'terminal q/p':>14}")
        for n in (200, 100, 40, 20, 10, 5, 2, 1):
            r, m, _ = drive_prescribed(row, path, n, "implicit", scheme=2,
                                       consts=consts)
            print(f"    {0.02/n:>10.2e}{n:>7}{m.cnt.be_calls:>10}"
                  f"{m.cnt.be_lowp_trial:>12}{m.cnt.be_newton_fail:>12}"
                  f"{m.cnt.be_recursive:>11}{m.cnt.be_explicit_fallback:>13}"
                  f"{r.eta[-1]:>14.5f}")
            out[(label, "s2", n)] = m.cnt
    return out


def control_A_equals_B_under_scheme2():
    """With a BackwardEuler companion the two D1 candidates are the SAME algorithm.

    BackwardEuler_CPPM sets  epsE = cEE + dEps - dGamma*ToCovariant(R)  at the converged
    state (`:2453-2455`), so d_eps_p == dGamma_n * ToCovariant(R_n) identically.  With a
    substepped explicit companion d_eps_p is the SUM of the substeps' plastic strains
    while mDGamma is only the LAST substep's multiplier, and the two diverge."""
    print("=" * 78)
    print("CONTROL -- implex_A == implex_B when the companion is BackwardEuler")
    print("=" * 78)
    meta, row = _seed("tx_p100_e0.6944_s1_n40")
    consts = dict(CONSTS); consts["e_init"] = float(meta["e_init"])
    for cs in (2, 1):
        rA, _, _ = drive_triaxial(row, 100.0, 40, 0.004, "implex_A",
                                  companion_scheme=cs, consts=consts)
        rB, _, _ = drive_triaxial(row, 100.0, 40, 0.004, "implex_B",
                                  companion_scheme=cs, consts=consts)
        d = max(norm_contr(a - b) for a, b in zip(rA.sig, rB.sig))
        print(f"  companion scheme {cs}: max |sigma~_A - sigma~_B| = {d:.3e} kPa "
              f"(||sigma|| ~ {norm_contr(rA.sig[-1]):.1f} kPa)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gate", default="all",
                    choices=["all", "G0", "G1", "G2", "G3", "G4", "G5",
                             "controls", "instr"])
    a = ap.parse_args()
    np.set_printoptions(precision=6, suppress=False, linewidth=140)
    if a.gate in ("all", "G0"):
        ok = gate_G0()[0]
        if not ok and a.gate == "all":
            print("\n*** G0 FAILED -- P0 stops here by the spec. "
                  "No IMPL-EX number below is admissible. ***")
            return
    if a.gate in ("all", "controls"):
        control_total_vs_incremental()
        control_forward_euler_r_bug()
        control_A_equals_B_under_scheme2()
    if a.gate in ("all", "G5"):
        gate_G5()
    if a.gate in ("all", "G1"):
        gate_G1()
    if a.gate in ("all", "G2"):
        gate_G2()
    if a.gate in ("all", "G3"):
        gate_G3()
    if a.gate in ("all", "G4"):
        gate_G4()
    if a.gate in ("all", "instr"):
        gate_instr()


if __name__ == "__main__":
    main()
