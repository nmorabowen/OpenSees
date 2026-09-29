"""(b) Calibrated behaviour on the TIMs campaign set, every R1 variant vs DM04.

Element tests at a material point, integrated exactly (the oracle's Radau with
events) along the whole path.  Axial = x (Voigt 0), compression positive.

Monotonic (one integrate() per test):
  TCd_p{25,100,400}   drained triaxial compression, e0 0.6944, eps_a to 10 %
  TEd_p100            drained triaxial extension, eps_a to -10 %
  TCu_p100_e{0.6944,0.80}  undrained triaxial compression, eps_a to 3 %
  PSd_p100            drained plane-strain compression (eps_zz = 0), eps_xx to 10 %
  SSd_p100            drained simple shear (sigma_yy const), gamma to 10 %
  SSu_p100            undrained (constant-volume) simple shear, gamma to 5 %
Cyclic:
  CTXu_*              stress-controlled undrained cyclic triaxial, q = +-q_cyc, to
                      eps_a double amplitude 5 % or N_max cycles
  CSSu_*              stress-controlled constant-volume cyclic simple shear,
                      tau = +-tau_cyc, to gamma DA 7.5 % or N_max
  CSSd_g*             strain-controlled drained cyclic simple shear gamma = +-g,
                      10 cycles: secant G, damping, eps_v per cycle
  CTXd_e*             strain-controlled drained cyclic triaxial, small amplitude
"""
from __future__ import annotations

import math

import numpy as np

import r1common as C
from sanisand_r1 import driver
from sanisand_r1.integrator import Control, integrate, path_table
from sanisand_r1.model import I3

E0 = 0.6944
RTOL_MONO = 1e-10
RTOL_CYC = 1e-9


def _tab_arrays(tab):
    return dict(
        eps=[r["eps"] for r in tab], sigma=[r["sigma"] for r in tab],
        p=[r["p"] for r in tab], q=[r["q"] for r in tab], e=[r["e"] for r in tab],
        rho_b=[r["rho_b"] for r in tab])


def _iso(p0, e0=E0):
    return driver.isotropic_state(p0, e0)


def _mono(state, ctl, O, rtol=RTOL_MONO):
    res = integrate(state, ctl, C.P, O, rtol=rtol)
    tab = path_table(res, C.P, O)
    d = _tab_arrays(tab)
    d.update(status=res.status, t_end=res.t_end, reseats=len(res.reseats),
             min_hx=res.min_h_over_x, n_cap=res.n_cap, segments=len(res.segments))
    return d


def TCd(p0, O):
    return _mono(_iso(p0), Control((True, False, False, True, True, True),
                                   (0.10, 0, 0, 0, 0, 0)), O)


def TEd(p0, O):
    return _mono(_iso(p0), Control((True, False, False, True, True, True),
                                   (-0.10, 0, 0, 0, 0, 0)), O)


def TCu(p0, e0, O):
    return _mono(_iso(p0, e0), Control.strain([0.03, -0.015, -0.015, 0, 0, 0]), O)


def PSd(p0, O):
    return _mono(_iso(p0), Control((True, False, True, True, True, True),
                                   (0.10, 0, 0, 0, 0, 0)), O)


def SSd(p0, O):
    return _mono(_iso(p0), Control((True, False, True, True, True, True),
                                   (0, 0, 0, 0.10, 0, 0)), O)


def SSu(p0, O):
    return _mono(_iso(p0), Control.strain([0, 0, 0, 0.05, 0, 0]), O)


# ---------------------------------------------------------------------------
# cyclic, stress-controlled: one integrate() per half cycle, stopped by an event
# ---------------------------------------------------------------------------
def _qtx(y):
    return y[0] - y[1]           # sigma_xx - sigma_yy (compression positive)


def _tau(y):
    return y[3]                  # sigma_xy


def cyclic_stress(kind, p0, e0, amp, O, n_max=30, stop_da=None, span=0.02):
    """kind 'CTXu' (q = +-amp) or 'CSSu' (tau = +-amp).  Undrained (isochoric)."""
    st = _iso(p0, e0)
    if kind == "CTXu":
        mk = lambda s: Control.strain([s, -0.5 * s, -0.5 * s, 0, 0, 0])
        meas, strain_i, stop_da = _qtx, 0, (stop_da or 0.05)
    else:
        mk = lambda s: Control.strain([0, 0, 0, s, 0, 0])
        meas, strain_i, stop_da = _tau, 3, (stop_da or 0.075)
    eps_acc = np.zeros(6)
    hist = dict(eps=[], sigma=[], p=[], q=[], half=[])
    halves, status, reseats, min_hx, caps = 0, "ok", 0, float("inf"), [0, 0, 0]
    emin, emax = 0.0, 0.0
    n_liq = None
    for k in range(2 * n_max):
        target = amp if k % 2 == 0 else -amp
        sgn = 1.0 if target > 0 else -1.0
        g = (lambda q_, y_, tg=target: meas(y_) - tg)
        res = integrate(st, mk(sgn * span), C.P, O, rtol=RTOL_CYC,
                        extra_events=[("target", g, 1 if sgn > 0 else -1)])
        reseats += len(res.reseats)
        min_hx = min(min_hx, res.min_h_over_x)
        caps = [a + b for a, b in zip(caps, res.n_cap)]
        tab = path_table(res, C.P, O)
        for r in tab[1:]:
            e6 = eps_acc + np.array(r["eps"])
            hist["eps"].append(e6.tolist()); hist["sigma"].append(r["sigma"])
            hist["p"].append(r["p"]); hist["q"].append(r["q"]); hist["half"].append(k)
            emin, emax = min(emin, e6[strain_i]), max(emax, e6[strain_i])
        eps_acc = eps_acc + np.array(res.deps_total)
        st = res.state
        halves = k + 1
        if res.status not in ("ok", "event:target"):
            status = res.status
            break
        if res.status == "ok":          # the target was never reached: runaway strain
            status = "runaway"
            n_liq = (k + 1) / 2.0
            break
        if emax - emin >= stop_da:
            n_liq = (k + 1) / 2.0
            status = "DA_reached"
            break
    return dict(hist=hist, status=status, halves=halves, n_liq=n_liq, reseats=reseats,
                min_hx=min_hx, n_cap=caps, da=emax - emin,
                p_end=float(np.trace(st.sigma)) / 3.0)


def cyclic_strain(kind, p0, e0, amp, O, n_cyc=10):
    """Drained strain-controlled cycles: kind 'CSSd' (gamma = +-amp, sigma_yy const)
    or 'CTXd' (eps_a = +-amp, lateral stresses const).  Returns per-cycle secant
    stiffness, damping ratio (loop area / (4 pi W)) and eps_v at the end of each
    cycle."""
    st = _iso(p0, e0)
    if kind == "CSSd":
        mk = lambda s: Control((True, False, True, True, True, True), (0, 0, 0, s, 0, 0))
        si, ti = 3, 3
    else:
        mk = lambda s: Control((True, False, False, True, True, True), (s, 0, 0, 0, 0, 0))
        si, ti = 0, None
    # first quarter to +amp, then halves of 2 amp
    legs = [amp] + [(-2 * amp if k % 2 == 0 else 2 * amp) for k in range(2 * n_cyc)]
    eps_acc = np.zeros(6)
    E, S = [0.0], [0.0 if ti is not None else 0.0]
    EV = []
    hist = dict(eps=[], sigma=[], p=[], q=[], leg=[])
    status, reseats, min_hx, caps = "ok", 0, float("inf"), [0, 0, 0]
    for k, d in enumerate(legs):
        res = integrate(st, mk(d), C.P, O, rtol=RTOL_CYC)
        reseats += len(res.reseats)
        min_hx = min(min_hx, res.min_h_over_x)
        caps = [a + b for a, b in zip(caps, res.n_cap)]
        tab = path_table(res, C.P, O)
        for r in tab[1:]:
            e6 = eps_acc + np.array(r["eps"])
            hist["eps"].append(e6.tolist()); hist["sigma"].append(r["sigma"])
            hist["p"].append(r["p"]); hist["q"].append(r["q"]); hist["leg"].append(k)
        eps_acc = eps_acc + np.array(res.deps_total)
        st = res.state
        EV.append(float(eps_acc[0] + eps_acc[1] + eps_acc[2]))
        if res.status != "ok":
            status = res.status
            break
    # per-cycle loops (legs 1..: each full cycle = two legs)
    eps = np.array(hist["eps"])[:, si]
    sig = np.array(hist["sigma"])
    tau = sig[:, 3] if kind == "CSSd" else sig[:, 0] - sig[:, 1]
    leg = np.array(hist["leg"])
    cyc = []
    for c in range(n_cyc):
        m = (leg == 1 + 2 * c) | (leg == 2 + 2 * c)
        if m.sum() < 4:
            break
        x, y = eps[m], tau[m]
        area = 0.5 * abs(np.sum(x[:-1] * y[1:] - x[1:] * y[:-1]))    # loop (shoelace)
        ks = (y.max() - y.min()) / (x.max() - x.min())
        W = 0.5 * ks * (0.5 * (x.max() - x.min())) ** 2
        cyc.append(dict(k_sec=ks, damping=area / (4 * math.pi * W) if W > 0 else float("nan"),
                        epsv_end=EV[min(2 + 2 * c, len(EV) - 1)]))
    return dict(hist=hist, status=status, cycles=cyc, reseats=reseats, min_hx=min_hx,
                n_cap=caps)


TESTS = {
    "TCd_p25": lambda O: TCd(25.0, O),
    "TCd_p100": lambda O: TCd(100.0, O),
    "TCd_p400": lambda O: TCd(400.0, O),
    "TEd_p100": lambda O: TEd(100.0, O),
    "TCu_p100_e0.6944": lambda O: TCu(100.0, E0, O),
    "TCu_p100_e0.80": lambda O: TCu(100.0, 0.80, O),
    "PSd_p100": lambda O: PSd(100.0, O),
    "SSd_p100": lambda O: SSd(100.0, O),
    "SSu_p100": lambda O: SSu(100.0, O),
}
CYCLIC = {
    # amplitudes chosen from a DM04 pilot (cyc_pilot.py) so that DM04 liquefies in
    # a moderate number of cycles; set in cyc_amps.json
}
