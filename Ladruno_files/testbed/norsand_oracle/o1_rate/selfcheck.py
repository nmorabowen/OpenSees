"""O1 self-checks (NOT the gate suite; P0d writes that).  Prints the numbers.

    cd Ladruno_files/testbed/norsand_oracle
    python -m o1_rate.selfcheck            # everything (~3-5 min)
    python -m o1_rate.selfcheck quick      # rtol check without 1e-12
    python -m o1_rate.selfcheck k2 rtol    # named groups only (see GROUPS)
"""
from __future__ import annotations

import math
import sys
import time
import warnings

import numpy as np

from . import Params, initial_state, run_path, tangent, triaxial, k2_path, acoustic_min_det
from .integrator import integrate_increment
from .localization import finite_spatial_tangent
from .model import (ISYM, SQ6, c2m, elastic_tangent_spectral, energy, plastic, zeta_fun, t2m, m2t)

QUICK = "quick" in sys.argv[1:]
K2 = dict()          # Params() defaults are the AB06 6.1 set (paper CSL, WW, no cap)
SIG0 = -100.0 * np.eye(3)


def hdr(s):
    print("\n== " + s)


def sym(A):
    return 0.5 * (A + A.T)


def check_zeta():
    hdr("K1.3 zeta corners (WW and GA)")
    for kind, rhos in (("WW", (0.5, 0.7, 0.8, 1.0)), ("GA", (7 / 9, 0.8, 1.0))):
        for r in rhos:
            z0, _ = zeta_fun(0.0, 1 / SQ6, r, kind)
            z60, _ = zeta_fun(math.pi / 3, -1 / SQ6, r, kind)
            print(f"  {kind} rho={r:.4f}: zeta(0)-1/rho = {z0 - 1 / r:+.1e}, zeta(pi/3)-1 = {z60 - 1:+.1e}")


def check_elastic_tangent():
    hdr("a^e: coordinate-free tensor form vs (S.3)+(S.33) spectral, and vs central FD of sigma(eps^e)")
    rng = np.random.default_rng(1)
    for P in (Params(), Params(alpha0=5.0, check=False)):
        for _ in range(3):
            ee = sym(rng.normal(size=(3, 3))) * 2e-3 - 0.004 * np.eye(3)
            C = energy(ee, P).a4
            Cs = elastic_tangent_spectral(ee, P)
            h = 1e-7
            Cfd = np.zeros_like(C)
            for k in range(3):
                for l in range(3):
                    E = np.zeros((3, 3))
                    E[k, l] += 0.5 * h
                    E[l, k] += 0.5 * h
                    Cfd[:, :, k, l] = (energy(ee + E, P, False).sig - energy(ee - E, P, False).sig) / (2 * h)
            Cref = np.einsum("ijmn,mnkl->ijkl", C, ISYM)
            print(f"  alpha0={P.alpha0}: |C - C_spectral|/|C| = {np.linalg.norm(C - Cs) / np.linalg.norm(C):.1e}"
                  f"   |C - C_FD|/|C| = {np.linalg.norm(Cref - Cfd) / np.linalg.norm(C):.1e}")


def check_iso():
    hdr("K1.1 isotropic compression, eps_v = -0.01 (elastic, pi_i0 = -300)")
    P = Params()
    s0 = initial_state(P, SIG0, 1.59, -300.0)
    s = run_path(P, s0, [np.eye(3) * (-0.01 / 3)])[-1]
    p = np.trace(s.sigma) / 3
    ex = -100.0 * math.e
    print(f"  p = {p:.10f}  closed form {ex:.10f}  rel err {abs(p - ex) / abs(ex):.1e}  status {s.flags['status']}"
          f" plastic={s.flags['plastic']}")


def check_loop():
    hdr("K1.2 closed non-coaxial elastic loop (pi_i0 = -300)")
    P = Params(rho=0.7, rho_bar=0.8)
    s0 = initial_state(P, SIG0 + np.diag([-20.0, 5.0, 0.0]), 1.59, -300.0)
    a = np.array([[-1e-3, 6e-4, 0.0], [6e-4, 2e-4, 0.0], [0.0, 0.0, 3e-4]])
    b = np.array([[2e-4, 0.0, -3e-4], [0.0, -8e-4, 5e-4], [-3e-4, 5e-4, 4e-4]])
    out = run_path(P, s0, np.array([a, b, -a, -b]))
    s = out[-1]
    W, Wa = s.flags["W"], s.flags["W_abs"]
    print(f"  all elastic: {not any(o.flags['plastic'] for o in out)}; W = {W:.3e}, int|dW| = {Wa:.3e}, "
          f"W/int|dW| = {abs(W) / Wa:.1e}")
    print(f"  |sigma_end - sigma_0|/|sigma_0| = {np.linalg.norm(s.sigma - s0.sigma) / np.linalg.norm(s0.sigma):.1e}, "
          f"|eps^e_end - eps^e_0| = {np.linalg.norm(s.eps_e - s0.eps_e):.1e}, pi_i change {s.pi_i - s0.pi_i:.1e}")


def check_dissipation_and_peak():
    hdr("K1.9 D >= 0 on drained TXC / TXE and undrained TXC (rho 0.7 / rho_bar 0.8); K1.6 at H = 0")
    P = Params(rho=0.7, rho_bar=0.8)
    for kind, ax, v0 in (("drained", -0.2, 1.59), ("drained", +0.05, 1.59), ("undrained", -0.2, 1.59),
                         ("drained", -0.2, 1.75)):
        s0 = initial_state(P, SIG0, v0, -60.4)
        out = triaxial(P, s0, kind, ax, 40)
        Ds = [o.D for o in out]
        rates = [o.flags["min_Dp_rate_rel"] for o in out if o.flags["min_Dp_rate_rel"] < float("inf")]
        Hz = out[-1].flags["H_zero"]
        nplast = sum(o.flags["plastic"] or o.D > 0 for o in out)
        print(f"  {kind:9s} ax={ax:+.2f} v0={v0}: status {out[-1].flags['status']}, plastic steps {nplast}, "
              f"min D_step = {min(Ds):.3e}, min pointwise sigma:eps_p rate/(|p||deps|) = "
              f"{min(rates) if rates else float('nan'):.3e}, max |F|/(M|p|) = "
              f"{max(o.flags['max_F_rel'] for o in out):.1e}")
        for hz in Hz:
            print(f"      H = 0 at eps_p_s = {hz['eps_p_s']:.5f}: D = {hz['D']:.12f}, chi psi_i = {hz['chi_psi']:.12f},"
                  f" |D - chi psi_i| = {abs(hz['D'] - hz['chi_psi']):.1e}")


def check_undrained_cs():
    hdr("K1.7 undrained critical state endpoint (TXC, rho 0.7 / 0.8), axial strain -8 in 80 increments")
    for csl, v0 in (("paper", 1.76), ("paper", 1.74), ("fork", 1.80)):
        P = Params(rho=0.7, rho_bar=0.8, csl_mode=csl)
        s0 = initial_state(P, SIG0, v0, -60.4)
        out = triaxial(P, s0, "undrained", -8.0, 80)
        s = out[-1]
        pq = plastic(s.sigma, s.pi_i, s.v, P)
        e = v0 - 1.0
        pcs = (-math.exp((P.v_c0 - v0) / P.lambda_tilde) if csl == "paper"
               else -P.p_a * ((P.e0 - e) / P.lambda_c) ** (1.0 / P.xi))
        print(f"  {csl} v0={v0}: p = {pq.p:.10f} vs p_cs = {pcs:.10f} (rel {abs(pq.p / pcs - 1):.1e}); "
              f"q = {pq.q:.8f} vs M|p_cs| = {P.M * abs(pcs):.8f} (rel {abs(pq.q / (P.M * abs(pcs)) - 1):.1e}); "
              f"psi_i = {pq.psi:.1e}, D = {pq.dilatancy:.1e}, H = {pq.H:.1e}, |v - v0| = {abs(s.v - v0):.1e}")


def check_drained_cs():
    hdr("K1.8 drained critical-state asymptote (TXC, rho 0.7 / 0.8), axial strain -6 in 60 increments; K1.10 v identity")
    for csl, v0 in (("paper", 1.59), ("fork", 1.70)):
        P = Params(rho=0.7, rho_bar=0.8, csl_mode=csl)
        s0 = initial_state(P, SIG0, v0, -60.4)
        out = triaxial(P, s0, "drained", -6.0, 60)
        for s in (out[9], out[29], out[-1]):
            pq = plastic(s.sigma, s.pi_i, s.v, P)
            print(f"  {csl} v0={v0} eps_a={s.flags['eps_total'][0, 0]:+.1f}: psi_i = {pq.psi:+.2e}, "
                  f"-zeta q/p - M = {-pq.zeta * pq.q / pq.p - P.M:+.2e}, pi_i/p - 1 = {s.pi_i / pq.p - 1:+.2e}, "
                  f"D = {pq.dilatancy:+.2e}, H = {pq.H:+.2e}, status {s.flags['status']}")
        # K1.10 (sheet 13.10, G2): v = v0 exp(tr eps) at every committed state, to round-off
        idr = max(abs(o.v / (v0 * math.exp(np.trace(o.flags["eps_total"]))) - 1.0) for o in out)
        ode = max(abs(o.flags["v_ode"] / o.v - 1.0) for o in out)
        x = float(np.trace(out[-1].flags["eps_total"]))
        print(f"  {csl} K1.10: max |v / (v0 exp(tr eps)) - 1| = {idr:.1e} over {len(out)} states; "
              f"ODE copy v' = v tr eps' vs algebraic {ode:.1e}; tr eps_end = {x:+.4f}, "
              f"linear-rule gap v0(1+x-e^x) = {v0 * (1 + x - math.exp(x)):+.3e}")


def check_rtol():
    hdr("rtol convergence, endpoint stress (drained TXC to -0.1 and undrained TXC to -0.1, 20 increments)")
    P = Params(rho=0.7, rho_bar=0.8)
    for kind in ("drained", "undrained"):
        res = {}
        for rt in ((1e-8, 1e-10) if QUICK else (1e-8, 1e-10, 1e-12)):
            t0 = time.time()
            s0 = initial_state(P, SIG0, 1.59, -60.4)
            res[rt] = triaxial(P, s0, kind, -0.1, 20, rtol=rt)[-1].sigma
            print(f"  {kind} rtol {rt:.0e}: {time.time() - t0:.1f} s")
        keys = sorted(res, reverse=True)
        for k1, k2 in zip(keys[:-1], keys[1:]):
            d = np.linalg.norm(res[k1] - res[k2]) / np.linalg.norm(res[k2])
            print(f"  {kind}: |sigma(rtol {k1:.0e}) - sigma(rtol {k2:.0e})| / |sigma| = {d:.1e}")


def check_refusal():
    hdr("K1.9 refusal counterexample N_bar = N = 0.4, rho 0.7, rho_bar 0.8")
    try:
        Params(N=0.4, N_bar=0.4, rho=0.7, rho_bar=0.8)
        print("  NOT refused  <-- FAIL")
    except ValueError as e:
        print(f"  refused: {e}")
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        Params(rho=0.8, rho_bar=0.7)
        print(f"  rho 0.8 > rho_bar 0.7 (condition A holds): warned = {len(w) > 0}")
    for kw in (dict(zeta="GA", rho=0.7, rho_bar=0.8), dict(rho=0.45, rho_bar=0.8), dict(N_bar=0.5),
               dict(rho=0.5, rho_bar=0.5), dict(rho=0.5, rho_bar=0.8), dict(rho=0.7, rho_bar=0.5)):
        try:
            Params(**kw)
            print(f"  {kw}: NOT refused  <-- FAIL")
        except ValueError:
            print(f"  {kw}: refused")
    # forced: a state on the yield surface at the extension corner (theta = 0), p/pi_i = 0.1
    P = Params(N=0.4, N_bar=0.4, rho=0.7, rho_bar=0.8, check=False)
    pi, p = -100.0, -10.0
    eta = (P.M / P.N) * (1 - (1 - P.N) * (p / pi) ** (P.N / (1 - P.N)))
    q = -p * eta * P.rho                                  # zeta(0) = 1/rho
    sig = np.diag([p - q / 3, p - q / 3, p + 2 * q / 3])  # TXE: sigma1 = sigma2 < sigma3
    s0 = initial_state(P, sig, 1.59, pi)
    pq = plastic(s0.sigma, s0.pi_i, s0.v, P)
    C = energy(s0.eps_e, P).a4
    d = m2t(np.linalg.solve(c2m(C), t2m(pq.f)))
    d *= 1e-6 / np.linalg.norm(d)
    s1 = integrate_increment(P, s0, d)
    print(f"  forced (validation bypassed): eta = {eta:.4f}, theta = {pq.theta:.2e}, F_rel = "
          f"{pq.F / (P.M * abs(pq.p)):.1e}; one step along a^e^-1:f -> status {s1.flags['status']}, "
          f"plastic {s1.flags['plastic']}, D_step = {s1.D:.4e} (< 0 expected)")
    # also along a drained TXE path from isotropic
    s0 = initial_state(P, SIG0, 1.59, -60.4)
    out = triaxial(P, s0, "drained", +0.2, 80)
    Dmin = min(o.D for o in out)
    print(f"  forced drained TXE path from isotropic -100: min D_step = {Dmin:.3e} "
          f"(status {out[-1].flags['status']}, {len(out)} steps)")


def check_finite_tangent():
    hdr("(S.34)/(S.44) finite-strain spatial tangent vs FD of tau(f b^e f^T) (elastic, non-coaxial g)")
    P = Params(rho=0.7, rho_bar=0.8)
    s0 = initial_state(P, np.diag([-110.0, -150.0, -90.0]), 1.59, -400.0)
    ae = tangent(P, s0, plastic_branch=False)
    A = finite_spatial_tangent(s0, ae)
    rng = np.random.default_rng(3)
    be = None
    w, V = np.linalg.eigh(s0.eps_e)
    be = V @ np.diag(np.exp(2 * w)) @ V.T

    def tau_of(g, h):
        f = np.eye(3) + h * g
        b = f @ be @ f.T
        ww, VV = np.linalg.eigh(b)
        ee = VV @ np.diag(0.5 * np.log(ww)) @ VV.T
        return energy(ee, P, False).sig
    for _ in range(3):
        g = rng.normal(size=(3, 3))
        h = 1e-6
        dtau = (tau_of(g, h) - tau_of(g, -h)) / (2 * h)
        fd = dtau - s0.sigma @ g.T
        an = np.einsum("ijkl,kl->ij", A, g)
        print(f"  |a:g - (dtau - tau g^T)| / |a:g| = {np.linalg.norm(an - fd) / np.linalg.norm(fd):.1e}")


def check_cont_tangent_fd():
    hdr("continuum a^ep (S.42) vs the integrated response (one-sided, loading directions, h -> 0)")
    P = Params(rho=0.7, rho_bar=0.8)
    s0 = initial_state(P, SIG0, 1.59, -60.4)
    s = triaxial(P, s0, "drained", -0.02, 4)[-1]
    C = tangent(P, s)
    rng = np.random.default_rng(5)
    base = np.diag([-1.0, 0.5, 0.5])
    for _ in range(3):
        d = base + 0.3 * sym(rng.normal(size=(3, 3)))
        d /= np.linalg.norm(d)
        errs = []
        for h in (1e-5, 1e-6, 1e-7):
            s1 = integrate_increment(P, s, h * d, rtol=1e-12)
            fd = (s1.sigma - s.sigma) / h
            an = np.einsum("ijkl,kl->ij", C, d)
            errs.append(np.linalg.norm(fd - an) / np.linalg.norm(an))
        print("  rel err at h = 1e-5, 1e-6, 1e-7: " + ", ".join(f"{e:.1e}" for e in errs)
              + f"   (plastic: {s1.flags['plastic']})")


def check_caps():
    hdr("cap modes x (WW 0.7/0.8, GA 0.8/0.8): near-isotropic compression diag(-1,-.9,-.9)e-3 x 30, and hydrostatic")
    for cap in ("none", "planar", "smooth"):
        for zk, r, rb in (("WW", 0.7, 0.8), ("GA", 0.8, 0.8)):
            kw = dict(c1=0.1, c2=0.1) if cap == "planar" else {}
            P = Params(rho=r, rho_bar=rb, cap=cap, zeta=zk, **kw)
            s0 = initial_state(P, SIG0, 1.59, -60.4)
            out = run_path(P, s0, np.array([np.diag([-1.0, -0.9, -0.9]) * 1e-3] * 30))
            s = out[-1]
            pq = plastic(s.sigma, s.pi_i, s.v, P)
            print(f"  {cap:6s} {zk}: {len(out)} steps, status {s.flags['status']}, eta = {pq.eta:.4f}, w = {pq.w:.3f}, "
                  f"min D_step = {min(o.D for o in out):.2e}, max |F|/(M|p|) = {max(o.flags['max_F_rel'] for o in out):.1e}")
            out = run_path(P, s0, np.array([np.eye(3) * -1e-3] * 20))
            pc = s0.pi_i / (1 - P.N) ** ((1 - P.N) / P.N)
            print(f"         hydrostatic x20: p_end = {np.trace(out[-1].sigma) / 3:.10f} vs pi_c = {pc:.10f}, "
                  f"pi_i frozen: {out[-1].pi_i == s0.pi_i}, status {out[-1].flags['status']}")


def _k2_case(args):
    pi0, chi, vc, rho, rb = args
    P = Params(rho=rho, rho_bar=rb, chi=chi, v_c0=vc)
    r = k2_path(P, initial_state(P, SIG0, 1.59, pi0), 60, extra_after=0, n_grid=91)
    return args, r["n_first"], r["n_interp"]


def check_k2_sensitivity():
    hdr("K2 sensitivity table (sheet 14 band): pi_i0 x chi x v_c0, both crossing criteria")
    from multiprocessing import Pool
    cases = [(pi0, chi, vc, rho, rb) for pi0 in (-60.4, -80.0, -100.0) for chi in (-3.0, -3.5, -4.0)
             for vc in (1.80, 1.81, 1.82) for rho, rb in ((0.7, 0.8), (1.0, 1.0))]
    with Pool(min(len(cases), 24)) as pool:
        res = {a: (f, i) for a, f, i in pool.map(_k2_case, cases)}
    rows = []
    for pi0 in (-60.4, -80.0, -100.0):
        for chi in (-3.0, -3.5, -4.0):
            for vc in (1.80, 1.81, 1.82):
                f7, i7 = res[(pi0, chi, vc, 0.7, 0.8)]
                f1, i1 = res[(pi0, chi, vc, 1.0, 1.0)]
                ok = (None not in (f7, f1, i7, i1) and f7 < f1 and 2 <= f1 - f7 <= 6
                      and i7 < i1 and 2 <= i1 - i7 <= 6)
                rows.append(ok)
                gi = f"{i1 - i7:.2f}" if None not in (i7, i1) else "-"
                print(f"  pi0 {pi0:6.1f} chi {chi:4.1f} vc0 {vc:.2f}: first-step {f7}/{f1}, interp "
                      f"{i7 if i7 is None else round(i7, 2)}/{i1 if i1 is None else round(i1, 2)} (gap {gi})"
                      f"  {'in band' if ok else 'OUT OF BAND'}")
    print(f"  {sum(rows)}/{len(rows)} combinations satisfy ordering and gap in [2, 6] under both criteria")


def check_k2():
    hdr("K2 nominal (pi_i0 = -60.4, chi = -3.5, v_c0 = 1.81, v0 = 1.59), O1 continuum tangent")
    for rho, rb in ((0.7, 0.8), (1.0, 1.0)):
        P = Params(rho=rho, rho_bar=rb)
        s0 = initial_state(P, SIG0, 1.59, -60.4)
        t0 = time.time()
        r = k2_path(P, s0, 40)
        print(f"  rho {rho} / rho_bar {rb}: first-step criterion n = {r['n_first']}, interpolated "
              f"n = {r['n_interp']:.3f} (nearest step {r['n_interp_step']}), status {r['status']}, "
              f"{time.time() - t0:.0f} s")
        k = r["n_first"] - 1
        print(f"      normalised min det at n-2..n: {np.round(r['mindet_norm'][k - 2:k + 1], 4).tolist()},"
              f" normal at n: {np.round(r['normals'][k], 4).tolist()}")


GROUPS = dict(
    basic=[check_zeta, check_elastic_tangent, check_iso, check_loop, check_refusal],
    tangents=[check_finite_tangent, check_cont_tangent_fd],
    dissipation=[check_dissipation_and_peak],
    rtol=[check_rtol],
    caps=[check_caps],
    undrained_cs=[check_undrained_cs],
    drained_cs=[check_drained_cs],
    k2=[check_k2],
    k2_sensitivity=[check_k2_sensitivity],
)

if __name__ == "__main__":
    # python -m o1_rate.selfcheck [quick] [group ...]   (no group: all)
    t0 = time.time()
    names = [a for a in sys.argv[1:] if a != "quick"] or list(GROUPS)
    for g in names:
        for fn in GROUPS[g]:
            fn()
    print(f"\ntotal {time.time() - t0:.0f} s")
