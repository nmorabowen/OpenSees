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

from . import Params, initial_state, pi_of_eta, run_path, tangent, triaxial, k2_path, acoustic_min_det
from .integrator import integrate_increment
from .localization import finite_spatial_tangent
from .model import (ISYM, SQ6, SQ23, SQ32, c2m, continuum_tangent, elastic_strain_from_stress,
                    elastic_tangent_spectral, energy, energy_psi, floor_energy, floor_phi, floor_project, floor_target, outer,
                    plastic, zeta_fun, t2m, m2t)

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


# =======================================================================================
# Round 3 / 3b: HAR energy (sheet 2.3-2.4), p' floor (9.7, S.55), unified pi_i0 rule (S.53)
# =======================================================================================
HAR_K, HAR_G = 1889.48104361, 807.80387674        # TIMs set (sheet 2.3), n = 1/2
PA_SET = (101.0,)               # TIMs p_a = 101 kPa (campaign "Patm 101"; amendment A4: 101.325 is WITHDRAWN)


def har_params(p_a, **kw):
    base = dict(energy="HAR", k=HAR_K, g=HAR_G, n_e=0.5, p_a=p_a)
    base.update(kw)
    return Params(**base)


def tims_params(p_a, **kw):
    """K1.14b set: TIMs elastic (HAR) + M 1.3309, rho = rho_bar = 0.71 WW, fork CSL, K2 plastic constants."""
    base = dict(M=1.3309, rho=0.71, rho_bar=0.71, csl_mode="fork", e0=0.83, lambda_c=0.027, xi=0.45,
                N=0.4, N_bar=0.2, chi=-3.5, h=280.0, p_min="default")
    base.update(kw)
    return har_params(p_a, **base)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def sig_at(p, q, theta):
    """Principal stress with mean p, deviator q, Lode angle theta (0 = TXE corner, pi/3 = TXC), axes x,y,z."""
    nh = SQ23 * np.array([math.cos(theta), math.cos(theta - 2 * math.pi / 3), math.cos(theta + 2 * math.pi / 3)])
    return np.diag(p + SQ23 * q * nh), np.diag(nh)


def v_for_psi(P, pi, psi):
    """Specific volume giving psi_i at pi_i (sheet 6)."""
    if P.csl_mode == "paper":
        return psi + P.v_c0 - P.lambda_tilde * math.log(-pi)
    return 1.0 + psi + P.e0 - P.lambda_c * (-pi / P.p_a) ** P.xi


def check_har_k1():
    hdr("HAR K1.1h / K1.11 closed forms (TIMs set k 1889.48104361, g 807.80387674, n 1/2): energy and O1 rate paths")
    for pa in PA_SET:
        P = har_params(pa, rho=1.0, rho_bar=1.0)
        kn = P.k * (1 - P.n_e)
        print(f"  p_a = {pa}: domain edge 1/(k(1-n)) = {1 / kn:.9e}")
        s0 = initial_state(P, -pa * np.eye(3), 1.6, -1.0e5)
        print(f"    eps^e at sigma = -p_a 1: |eps^e| = {np.linalg.norm(s0.eps_e):.1e} (origin shift: p = -p_a at eps^e = 0)")
        for ev in (-1e-3, 5e-4):
            pex = -pa * (1 - kn * ev) ** (1 / (1 - P.n_e))
            el = energy(np.eye(3) * ev / 3, P)
            Kex = P.k * pa * (abs(pex) / pa) ** P.n_e
            st = run_path(P, s0, [np.eye(3) * ev / 3])[-1]
            prate = np.trace(st.sigma) / 3
            print(f"    13.1h eps_v = {ev:+.0e}: p = {el.p:.6f} (closed {pex:.6f}, rel {rel(el.p, pex):.1e}), K = D11 = "
                  f"{el.D11:.3f} (closed {Kex:.3f}, rel {rel(el.D11, Kex):.1e}); O1 rate path p = {prate:.9f} "
                  f"(rel {rel(prate, pex):.1e}, status {st.flags['status']}, plastic {st.flags['plastic']})")
        for es in (8.66546968e-4, 1e-3):
            fac = (1 + 3 * P.g * kn * es * es) ** (P.n_e / (2 * (1 - P.n_e)))
            pex, qex = -pa * fac, 3 * P.g * pa * es * fac
            d = np.diag([-es, 0.5 * es, 0.5 * es])        # isochoric TXC, eps_s = es
            st = run_path(P, s0, [d])[-1]
            p, q = np.trace(st.sigma) / 3, SQ32 * np.linalg.norm(st.sigma - np.trace(st.sigma) / 3 * np.eye(3))
            print(f"    13.11 eps_s = {es:.8e}: O1 p = {p:.6f} (closed {pex:.6f}, rel {rel(p, pex):.1e}), q = {q:.6f} "
                  f"(closed {qex:.6f}, rel {rel(q, qex):.1e}), eta = q/|p| = {q / -p:.8f} vs 3 g eps_s = "
                  f"{3 * P.g * es:.8f} (rel {rel(q / -p, 3 * P.g * es):.1e}), status {st.flags['status']}")
        # the same path under BA06 (K2 set): p constant (the discriminating check)
        Pb = Params()
        sb = initial_state(Pb, SIG0, 1.6, -1.0e5)
        stb = run_path(Pb, sb, [np.diag([-1e-3, 5e-4, 5e-4])])[-1]
        print(f"    BA06 on the same isochoric shear: p = {np.trace(stb.sigma) / 3:.12f} (constant -100: "
              f"{abs(np.trace(stb.sigma) / 3 + 100):.1e})")
        # inverse map spot value p = -3.5, eta = 2.1 (13.11) and round trip
        sig, _ = sig_at(-3.5, 2.1 * 3.5, math.pi / 3)
        ee = elastic_strain_from_stress(sig, P)
        el = energy(ee, P, tangent=False)
        ev, es = np.trace(ee), SQ23 * np.linalg.norm(ee - np.trace(ee) / 3 * np.eye(3))
        vp = math.sqrt(3.5 ** 2 + kn * (2.1 * 3.5) ** 2 / (3 * P.g))
        print(f"    (S.5h'') p = -3.5, eta = 2.1: varpi = {vp:.7f}, eps_v = {ev:.7e}, eps_s = {es:.7e}; round trip "
              f"|sigma(eps) - sigma|/|sigma| = {np.linalg.norm(el.sig - sig) / np.linalg.norm(sig):.1e}")
    print("  (sheet 13.1h / 13.11 values are quoted at p_a = 101: p = -381.983585 / -28.117707, K = 371129.585;"
          " eps_s 8.66546968e-4: p -166.548671, q 349.752210; 1e-3: p -183.183351, q 443.928664, eta 2.42341163;"
          " varpi 5.7714886, eps_v 9.0504735e-4, eps_s 1.2561906e-4)")


def check_har_tangent_loop():
    hdr("HAR a^e (t2 and t3/t4 terms of (S.3) live): tensor form vs spectral (S.3)+(S.33) vs FD; Psi gradient; det D; elastic loop W = 0")
    P = har_params(101.0, rho=0.7, rho_bar=0.8)
    rng = np.random.default_rng(11)
    for _ in range(3):
        sig = np.diag([-100.0, -60.0, -150.0]) + sym(rng.normal(size=(3, 3))) * 10
        ee = elastic_strain_from_stress(sig, P)
        el = energy(ee, P)
        C = el.a4
        Cs = elastic_tangent_spectral(ee, P)
        h = 1e-8
        Cfd = np.zeros_like(C)
        for k in range(3):
            for l in range(3):
                E = np.zeros((3, 3))
                E[k, l] += 0.5 * h
                E[l, k] += 0.5 * h
                Cfd[:, :, k, l] = (energy(ee + E, P, False).sig - energy(ee - E, P, False).sig) / (2 * h)
        Cref = np.einsum("ijmn,mnkl->ijkl", C, ISYM)
        gfd = np.zeros((3, 3))
        for k in range(3):
            for l in range(3):
                E = np.zeros((3, 3))
                E[k, l] = 1e-7
                gfd[k, l] = (energy_psi(ee + E, P) - energy_psi(ee - E, P)) / 2e-7
        Z = 1 + P.k * (1 - P.n_e) * (el.q / el.p) ** 2 / (3 * P.g)
        vp = abs(el.p) * math.sqrt(Z)
        detD = el.D11 * el.D22 - el.D12 ** 2
        detex = 3 * P.k * P.g * P.p_a ** 2 * (vp / P.p_a) ** (2 * P.n_e)
        print(f"  eta {el.q / -el.p:.3f}: |C-C_spec|/|C| = {np.linalg.norm(C - Cs) / np.linalg.norm(C):.1e}, |C-C_FD|/|C| = "
              f"{np.linalg.norm(Cref - Cfd) / np.linalg.norm(C):.1e}; |dPsi/deps - sigma|/|sigma| = "
              f"{np.linalg.norm(gfd - el.sig) / np.linalg.norm(el.sig):.1e}; D12 = {el.D12:.1f}, D22/(q/eps_s) = "
              f"{el.D22 / el.q_over_es:.6f} vs (1-n/Z)/(1-n) = {(1 - P.n_e / Z) / (1 - P.n_e):.6f}; det D rel err "
              f"{rel(detD, detex):.1e}")
    # BA06-swapped tangent (the t4 term dropped, AB06 eq 64 form) is visibly wrong under HAR
    ee = elastic_strain_from_stress(np.diag([-100.0, -60.0, -150.0]), P)
    el = energy(ee, P)
    from .model import _a4_from_D
    e = ee - np.trace(ee) / 3 * np.eye(3)
    C64 = _a4_from_D(e, np.linalg.norm(e), el.D11, el.D12, el.q_over_es, el.q_over_es)
    print(f"  AB06-eq-64 form (D22 -> q/eps_s) under HAR: |C64 - C|/|C| = {np.linalg.norm(C64 - el.a4) / np.linalg.norm(el.a4):.2e}")
    s0 = initial_state(P, np.diag([-120.0, -95.0, -100.0]), 1.6, -1.0e5)
    a = np.array([[-1e-4, 6e-5, 0.0], [6e-5, 2e-5, 0.0], [0.0, 0.0, 3e-5]])
    b = np.array([[2e-5, 0.0, -3e-5], [0.0, -8e-5, 5e-5], [-3e-5, 5e-5, 4e-5]])
    out = run_path(P, s0, np.array([a, b, -a, -b]))
    s = out[-1]
    W, Wa = s.flags["W"], s.flags["W_abs"]
    dPsi = energy_psi(s.eps_e, P) - energy_psi(s0.eps_e, P)
    print(f"  HAR closed non-coaxial loop: all elastic {not any(o.flags['plastic'] for o in out)}, W = {W:.3e}, "
          f"int|dW| = {Wa:.3e}, W/int|dW| = {abs(W) / Wa:.1e}; |sigma_end - sigma_0|/|sigma_0| = "
          f"{np.linalg.norm(s.sigma - s0.sigma) / np.linalg.norm(s0.sigma):.1e}; Psi_end - Psi_0 = {dPsi:.1e}")
    o1 = run_path(P, s0, np.array([a, b]))
    print(f"  half loop: W = {o1[-1].flags['W']:.10e} vs Psi(end) - Psi(0) = "
          f"{energy_psi(o1[-1].eps_e, P) - energy_psi(s0.eps_e, P):.10e}")


def check_har_parser():
    hdr("Parser refusals (sheet 2.4, 9.7, S.56 with amendments A2/A3)")
    cases = [("HAR + p0", dict(energy="HAR", k=HAR_K, g=HAR_G, n_e=0.5, p0=-100.0)),
             ("HAR + alpha0", dict(energy="HAR", k=HAR_K, g=HAR_G, n_e=0.5, alpha0=0.0)),
             ("HAR without g", dict(energy="HAR", k=HAR_K, n_e=0.5)),
             ("HAR n = 1", dict(energy="HAR", k=HAR_K, g=HAR_G, n_e=1.0)),
             ("HAR n < 0", dict(energy="HAR", k=HAR_K, g=HAR_G, n_e=-0.1)),
             ("HAR k <= 0", dict(energy="HAR", k=0.0, g=HAR_G, n_e=0.5)),
             ("HAR p_a <= 0", dict(energy="HAR", k=HAR_K, g=HAR_G, n_e=0.5, p_a=0.0)),
             ("BA06 + k", dict(k=HAR_K)),
             ("energy XYZ", dict(energy="XYZ")),
             ("p_min < 0", dict(p_min=-1.0)),
             ("smooth c1 .05 c2 .06", dict(cap="smooth", c1=0.05, c2=0.06))]
    for lab, kw in cases:
        try:
            Params(**kw)
            print(f"  {lab}: NOT refused  <-- FAIL")
        except ValueError as e:
            print(f"  {lab}: refused ({str(e)[:70]})")
    for lab, kw in (("smooth c1 .05 c2 .07", dict(cap="smooth", c1=0.05, c2=0.07)),
                    ("smooth defaults .05/.15", dict(cap="smooth")),
                    ("planar c1 = c2 = .10", dict(cap="planar", c1=0.1, c2=0.1)), ("no cap", dict(cap="none")),
                    ("HAR + fork p_a 101", dict(energy="HAR", k=HAR_K, g=HAR_G, n_e=0.5, p_a=101.0, csl_mode="fork"))):
        P = Params(**kw)
        print(f"  {lab}: accepted, W_ramp = {P.W_ramp:.4f}, p_ref = {P.p_ref}, p_min('default') = "
              f"{Params(**dict(kw, p_min='default')).pmin:.6f}")


def _floor_report(lab, out, P, s0, extra=""):
    s = out[-1]
    segs = [(sg["mode"], sg["event"]) for o in out for sg in o.flags["segments"]]
    pq = plastic(s.sigma, s.pi_i, s.v, P)
    dPsi = energy_psi(s.eps_e, P) - energy_psi(s0.eps_e, P)
    W = s.flags["W"] - s0.flags.get("W", 0.0)
    print(f"  {lab}: status {s.flags['status']}, segments {segs}")
    print(f"      p_end = {np.trace(s.sigma) / 3:.10f}, eps^f_v = {s.flags['eps_f_v']:.10e}, W_f = {s.flags['W_f']:.10e}, "
          f"F/(M|p|) = {pq.F / (P.M * abs(pq.p)):+.1e}, max floor drift {max(o.flags['max_floor_drift'] for o in out):.1e}; "
          f"energy balance Psi_end - Psi_0 - W - (W_f - W_f0) = {dPsi - W - (s.flags['W_f'] - s0.flags.get('W_f', 0.0)):.1e}"
          + extra)
    return s


def check_floor_k1_12():
    hdr("K1.12 BA06 floor (K2 set, p_min 'default' = 5e-3|p0| = 0.5 kPa): Pi_f closed form (S.49) and the rate form (S.55)")
    P = Params(p_min="default")
    ft = floor_target(0.0, P)
    print(f"  eps_v,f(0) = {ft.ev_f:.10f} (sheet 0.0529831737), p_min = {P.pmin}")
    ev_tr = P.eps_v0 - P.kappa_hat * math.log(0.25 / 100.0)
    ee_tr = np.eye(3) * ev_tr / 3
    eef, act, dvf, _ = floor_project(ee_tr, P)
    Ef, Wf1 = floor_energy(ee_tr, P)
    print(f"  split Pi_f of the trial p^tr = -p_min/2 (eps_v,tr = {ev_tr:.10f}): active {act}, Delta eps^f_v = {dvf:.12e} "
          f"(kappa ln2 = {P.kappa_hat * math.log(2):.12e}), W_f = {P.pmin * dvf:.12e}, E_f = Psi(f) - Psi(tr) = {Ef:.12e} "
          f"(closed kappa(p_min - |p_tr|) = {P.kappa_hat * 0.25:.4e}), p_f = {energy(eef, P, False).p:.15f}")
    s0 = initial_state(P, SIG0, 1.59, -300.0)
    out = run_path(P, s0, [np.eye(3) * ev_tr / 3])
    s = _floor_report("rate, one isotropic increment to eps_v,tr", out, P, s0)
    print(f"      Delta eps^f_v rel err vs kappa ln2: {rel(s.flags['eps_f_v'], P.kappa_hat * math.log(2)):.1e}; pi_i unchanged "
          f"{s.pi_i == s0.pi_i}; v = {s.v:.12f} vs v0 exp(tr eps) = {1.59 * math.exp(ev_tr):.12f}; rate-form stored-energy "
          f"gain of the floor = W_f (the (S.52) bound attained; the split's E_f is {Ef:.4e})")
    C = tangent(P, s)
    blk = np.array([[C[a, a, b, b] for b in range(3)] for a in range(3)])
    print(f"      floored tangent principal block / (2 mu0) = {np.round(blk / (2 * P.mu0), 12).tolist()}; "
          f"|1:C|/|C| = {np.linalg.norm(np.einsum('iikl->kl', C)) / np.linalg.norm(C):.1e}")
    # with deviatoric strain at the floor (q kept under BA06 alpha0 = 0)
    s1 = run_path(P, s, [np.diag([2e-5, -1e-5, -1e-5]) + np.eye(3) * 1e-3])[-1]
    q0 = SQ32 * np.linalg.norm(s.sigma - np.trace(s.sigma) / 3 * np.eye(3))
    q1 = SQ32 * np.linalg.norm(s1.sigma - np.trace(s1.sigma) / 3 * np.eye(3))
    print(f"      + expansion 3e-3 with deviator: p = {np.trace(s1.sigma) / 3:.12f}, q = {q1:.6f} (3 mu0 eps_s = "
          f"{3 * P.mu0 * SQ23 * np.linalg.norm(s1.eps_e - np.trace(s1.eps_e) / 3 * np.eye(3)):.6f}), "
          f"eps^f_v += {s1.flags['eps_f_v'] - s.flags['eps_f_v']:.6e} (= tr deps - 0 since p pinned: 3e-3), status "
          f"{s1.flags['status']}, mode {s1.flags['mode']}")
    # leaving the floor: compression gives back full stiffness (one-sided)
    s2 = run_path(P, s1, [np.eye(3) * -1e-3])[-1]
    print(f"      then compression 3e-3: mode {s2.flags['mode']}, p = {np.trace(s2.sigma) / 3:.10f} (closed "
          f"{-P.pmin * math.exp(3e-3 / P.kappa_hat):.10f}), eps^f_v unchanged {s2.flags['eps_f_v'] == s1.flags['eps_f_v']}")
    # M-F4: BA06 at p = -1 with +1e-4 / +1.1e-4: ordinary elastic steps
    sb = initial_state(P, -1.0 * np.eye(3), 1.6, -300.0)
    for dv in (1e-4, 1.1e-4):
        st = run_path(P, sb, [np.eye(3) * dv / 3])[-1]
        print(f"  M-F4 (BA06, p = -1, +{dv:.1e}): p = {np.trace(st.sigma) / 3:.10f} (closed {-math.exp(-dv / P.kappa_hat):.10f}),"
              f" floor segments {st.flags['floor_segments']}, eps^f_v = {st.flags['eps_f_v']:.1e}")


def check_floor_k1_13_14():
    hdr("K1.13 / K1.14 HAR floor (TIMs set, p_min = 5e-3 p_a): (S.50) closed forms, out-of-domain trial, rate form")
    for pa in PA_SET:
        P = har_params(pa, p_min="default", rho=1.0, rho_bar=1.0)
        kn = P.k * (1 - P.n_e)
        ft = floor_target(0.0, P)
        G = P.g * pa * (P.pmin / pa) ** 0.5
        print(f"  p_a = {pa}: p_min = {P.pmin}, eps_v,f(0) = {ft.ev_f:.8e}, eps*_f(0) = {1 / kn - ft.ev_f:.4e} "
              f"(margin (p_min/p_a)^(1-n)/(k(1-n)) = {(P.pmin / pa) ** 0.5 / kn:.4e}), G(p_min) = {G:.3f}, 2G = {2 * G:.2f}")
        ev1 = (1 - (1.0 / pa) ** (1 - P.n_e)) / kn
        s0 = initial_state(P, -1.0 * np.eye(3), 1.6, -1.0e3)
        print(f"    isotropic p = -1: eps_v = {np.trace(s0.eps_e):.8e} (closed {ev1:.8e})")
        for dv in (1e-4, 1.1e-4):
            eet = s0.eps_e + np.eye(3) * dv / 3
            elt = energy(eet, P, False)
            eef, act, dvf, _ = floor_project(eet, P)
            out = run_path(P, s0, [np.eye(3) * dv / 3])
            s = out[-1]
            Ef, Wf1 = floor_energy(eet, P)
            print(f"    +{dv:.1e}: trial in domain {elt.in_domain}, p^tr = {elt.p:.4e}; split Pi_f: Delta eps^f_v = {dvf:.8e}, "
                  f"E_f = {Ef if Ef is None else f'{Ef:.6e}'} (<= W_f {Wf1:.6e}); "
                  f"rate: eps^f_v = {s.flags['eps_f_v']:.8e} (rel {rel(s.flags['eps_f_v'], dvf):.1e}), W_f = {s.flags['W_f']:.8e}, "
                  f"p_end = {np.trace(s.sigma) / 3:.12f}, status {s.flags['status']}, modes "
                  f"{[sg['mode'] for sg in s.flags['segments']]}")
        C = tangent(P, s)
        blk = np.array([[C[a, a, b, b] for b in range(3)] for a in range(3)])
        print(f"    floored isotropic tangent block / (2G(p_min)) = {np.round(blk / (2 * G), 10).tolist()}, "
              f"|1:C|/|C| = {np.linalg.norm(np.einsum('iikl->kl', C)) / np.linalg.norm(C):.1e}")
        # K1.14: floored state under shear, eps_s = 2e-4
        es = 2e-4
        ft = floor_target(es, P)
        nh = SQ23 * np.array([1.0, -0.5, -0.5])
        ee = np.diag(np.full(3, 1.05e-3 / 3) + SQ32 * es * nh)
        eef, act, dvf, ft2 = floor_project(ee, P)
        el = energy(eef, P)
        es_f = SQ23 * np.linalg.norm(eef - np.trace(eef) / 3 * np.eye(3))
        print(f"    K1.14 eps_s = 2e-4: x = {ft.x:.11f}, eps_v,f = {ft.ev_f:.8e}, q_f = {ft.q_f:.7f}, eta_f = {ft.q_f / P.pmin:.2f}, "
              f"eps'_f = {ft.dev_f:.8f} (-D12/D11 at the floor {-el.D12 / el.D11:.8f}, rel {rel(ft.dev_f, -el.D12 / el.D11):.1e})")
        print(f"      Pi_f: p_f + p_min = {el.p + P.pmin:.1e}, q(eps_f) - q_f = {el.q - ft.q_f:.1e}, eps_s kept "
              f"{abs(es_f - es):.1e} (M-F6); x-equation residual {ft.x ** 2 - 3 * kn * P.g * es * es * ft.x - (P.pmin / pa) ** 2:.1e}")
        Phi = floor_phi(eef, P)
        Cf = np.einsum("ijmn,mnkl->ijkl", el.a4, Phi)
        Cc, _ = continuum_tangent(eef, -1.0e3, 1.6, P, plastic_branch=False, floor_branch=True)
        Phi0 = ISYM - outer(np.eye(3), np.eye(3)) / 3.0                      # M-F3b: eps'_f dropped
        C0 = np.einsum("ijmn,mnkl->ijkl", el.a4, Phi0)
        # FD of sigma(Pi_f(eps)) (S.51a) at a non-degenerate state
        ee2 = np.diag([0.40e-3, 0.30e-3, 0.38e-3]) + 0.0
        ee2[0, 1] = ee2[1, 0] = 2e-5
        e2f = floor_project(ee2, P)[0]
        A2 = np.einsum("ijmn,mnkl->ijkl", energy(e2f, P).a4, floor_phi(e2f, P))
        Afd = np.zeros_like(A2)
        h = 1e-9
        for k in range(3):
            for l in range(3):
                E = np.zeros((3, 3))
                E[k, l] += 0.5 * h
                E[l, k] += 0.5 * h
                Afd[:, :, k, l] = (energy(floor_project(ee2 + E, P)[0], P, False).sig
                                   - energy(floor_project(ee2 - E, P)[0], P, False).sig) / (2 * h)
        A2r = np.einsum("ijmn,mnkl->ijkl", A2, ISYM)
        blk0 = np.array([[C0[a, a, b, b] for b in range(3)] for a in range(3)])        # principal block (S.51)
        mf3b_blk = np.abs(blk0.sum(axis=0)).max() / np.abs(blk0).max()
        print(f"      tangent: |a^e Phi - C_continuum(S.55)|/|C| = {np.linalg.norm(Cf - Cc) / np.linalg.norm(Cc):.1e}; "
              f"|1:C|/max|C| = {np.abs(np.einsum('iikl->kl', Cf)).max() / np.abs(Cf).max():.1e}; M-F3b (eps'_f dropped): "
              f"|1:C|/max|C| = {np.abs(np.einsum('iikl->kl', C0)).max() / np.abs(C0).max():.2e} (principal block: delta:a/max|a| "
              f"= {mf3b_blk:.2e}); a^e Phi vs FD of "
              f"sigma(Pi_f(eps)) (h 1e-9) {np.linalg.norm(A2r - Afd) / np.linalg.norm(A2r):.1e}")
        # M-F3b with the sheet's own measure (floor_sympy.py (d)): principal strains (-3.1e-4, 4.7e-4, 2.2e-4) floored,
        # C = [d sigma_a / d eps_b] a^e Phi; the eps'-dropped mutant's max column sum / max|C| (correct C) = 1.45e-2
        ed = np.diag([-3.1e-4, 4.7e-4, 2.2e-4])          # p = -41 kPa: floor_sympy (d) applies the (S.48) map
        es_d = SQ23 * np.linalg.norm(ed - np.trace(ed) / 3 * np.eye(3))    # unconditionally; so do we here
        edf = ed - (np.trace(ed) - floor_target(es_d, P).ev_f) / 3.0 * np.eye(3)
        a4d = energy(edf, P).a4
        ab = np.array([[a4d[a, a, b, b] for b in range(3)] for a in range(3)])
        Phd = floor_phi(edf, P)
        Phb = np.array([[Phd[a, a, b, b] for b in range(3)] for a in range(3)])
        Ccl = ab @ Phb
        Cm = ab @ (np.eye(3) - np.ones((3, 3)) / 3.0)
        print(f"      sheet-(d) state (-3.1e-4, 4.7e-4, 2.2e-4): delta:C_f/max|C| = {np.abs(Ccl.sum(axis=0)).max() / np.abs(Ccl).max():.1e}; "
              f"M-F3b mutant delta:C/max|C| = {np.abs(Cm.sum(axis=0)).max() / np.abs(Ccl).max():.3e} (sheet 1.45e-2)")
        # general n: the bracketed scalar solve of (S.50) (n = 0.3, 0.7): p(Pi_f) = -p_min, idempotent
        for nn in (0.3, 0.7):
            Pn = har_params(pa, n_e=nn, p_min="default", rho=1.0, rho_bar=1.0)
            worst = 0.0
            for es_ in (0.0, 1e-5, 2e-4, 3e-3):
                eg = np.diag(np.full(3, 1.2 / (Pn.k * (1 - nn)) / 3) + SQ32 * es_ * nh)    # out of the domain
                gf, act, _, ftn = floor_project(eg, Pn)
                gf2, act2, _, _ = floor_project(gf, Pn)
                elg = energy(gf, Pn)
                worst = max(worst, abs(elg.p + Pn.pmin) / Pn.pmin, float(np.abs(gf2 - gf).max()) / 1e-3,
                            abs(ftn.dev_f + elg.D12 / elg.D11) / max(abs(ftn.dev_f), 1e-30) if es_ > 0 else 0.0)
            print(f"      general n = {nn}: out-of-domain trials floored, max over eps_s in (0, 1e-5, 2e-4, 3e-3) of "
                  f"|p_f + p_min|/p_min, idempotence, eps'_f vs -D12/D11: {worst:.1e} (re-projection active: {act2})")
    print("  (sheet values at p_a = 101: eps_v,f(0) 9.83645033e-4, G(p_min) 5769.156, eps_v(p=-1) 9.53167838e-4, Delta eps^f_v "
          "6.9522805e-5 / 7.9522805e-5, W_f 3.5109017e-5; K1.14 x 0.09185198375, eps_v,f 1.04102893e-3, q_f 14.8362051, "
          "eps'_f 0.08679793, M-F3b 1.45e-2)")


def _surface_state(P, p, eta, theta, psi, pmin_override=None):
    sig, nh = sig_at(p, 0.0, theta)
    z, _ = zeta_fun(theta, math.cos(3 * theta) / SQ6, P.rho, P.zeta)
    q = eta * abs(p) / z
    sig, nh = sig_at(p, q, theta)
    pi = pi_of_eta(P, p, eta)
    v = v_for_psi(P, pi, psi)
    return sig, nh, pi, v


def check_floor_fpf():
    hdr("K1.14b / A1 (HAR, dry side): the FPf case in rate form (S.55); chained-pattern analogues; tangent vs FD")
    for pa in PA_SET:
        P = tims_params(pa, p_min=0.505)
        sig, nh, pi, v = _surface_state(P, -0.6, 1.2 * P.M, 0.271, -0.10)
        s0 = initial_state(P, sig, v, pi)
        el = energy(s0.eps_e, P)
        pq = plastic(s0.sigma, s0.pi_i, s0.v, P)
        print(f"  p_a = {pa}: start p = {pq.p:.4f}, q = {pq.q:.4f}, eta = {-pq.zeta * pq.q / pq.p:.4f}, theta = {pq.theta:.3f}, "
              f"pi_i = {s0.pi_i:.5f}, v = {s0.v:.5f}, psi_i = {pq.psi:+.4f}, F/(M|p|) = {s0.flags['F_rel']:.0e}; eps_v = "
              f"{np.trace(s0.eps_e):.4e}, eps_s = {el.eps_s:.4e}, D11 {el.D11:.0f} D12 {el.D12:.0f} D22 {el.D22:.0f}")
        for lab, tr, sh, n_inc in (("FPf  (tr 2e-5, shear 2e-5)", 2e-5, 2e-5, 1),
                                   ("FPf,FPf (twice the increment)", 2e-5, 2e-5, 2),
                                   ("-P-,-Pf (tr 1.4e-5, shear 2e-5)", 1.4e-5, 2e-5, 1)):
            d = np.eye(3) * tr / 3 + sh * nh
            out = run_path(P, s0, np.array([d] * n_inc))
            s = out[-1]
            pe = plastic(s.sigma, s.pi_i, s.v, P)
            print(f"   {lab}: status {s.flags['status']}, segments "
                  f"{[(sg['mode'], sg['event'], round(sg['t1'], 4)) for o in out for sg in o.flags['segments']]}")
            print(f"      p_end = {pe.p:.10f}, q_end = {pe.q:.5f}, eta_c = {-pe.zeta * pe.q / pe.p:.4f}, pi_i = {s.pi_i:.5f}, "
                  f"eps^p_v = {s.eps_p_v:+.3e}, eps^p_s = {s.eps_p_s:.3e}, eps^f_v = {s.flags['eps_f_v']:.3e}, W_f = "
                  f"{s.flags['W_f']:.3e}, F/(M|p|) = {pe.F / (P.M * abs(pe.p)):+.1e}, min lam_f'/|deps| = "
                  f"{min(o.flags['min_lam_f_rel'] for o in out):.3e}, drift {max(o.flags['max_floor_drift'] for o in out):.1e}")
        # continuum tangent (S.55) at the FPf end state vs the integrated response (one-sided, h -> 0)
        d = np.eye(3) * 2e-5 / 3 + 2e-5 * nh
        s = run_path(P, s0, [d])[-1]
        C = tangent(P, s)
        u = d / np.linalg.norm(d)
        errs = []
        for h in (1e-7, 1e-8, 1e-9):
            s1 = integrate_increment(P, s, h * u, rtol=1e-12)
            errs.append(np.linalg.norm((s1.sigma - s.sigma) / h - np.einsum("ijkl,kl->ij", C, u))
                        / np.linalg.norm(np.einsum("ijkl,kl->ij", C, u)))
        Cnf = tangent(P, s, floor_branch=False)
        mut = np.linalg.norm(np.einsum("ijkl,kl->ij", Cnf - C, u)) / np.linalg.norm(np.einsum("ijkl,kl->ij", C, u))
        print(f"   tangent at the FPf end state (mode {s.flags['mode']}): rel err vs FD at h = 1e-7, 1e-8, 1e-9: "
              + ", ".join(f"{e:.1e}" for e in errs) + f"; |1:C|/|C| = "
              f"{np.linalg.norm(np.einsum('iikl->kl', C)) / np.linalg.norm(C):.1e}; floor mechanism dropped (plain a^ep): {mut:.2e}")
    print("  (split/BE record, sheet K1.14b at p_a 101: p^tr -0.4607 -> floored; return dlam 1.157e-5, eta_c 1.8200, pi_i -0.74335,"
          " eps^p_v +7.07e-6, eps^p_s 1.57e-5, p_c -0.4877 -> post floor -0.505, q 0.6619 -> 0.6705, F(sigma_f)/p_min -0.0041)")


def check_floor_wet():
    hdr("A5 (HAR, wet side): from the floor on the surface at eta = 0.5M, theta = 0.454, pure shear Dgamma along n^ (rate form)")
    P = tims_params(101.0, p_min=0.505)
    P0 = tims_params(101.0, p_min=0.0)
    kn = P.k * (1 - P.n_e)
    for psi in (0.0, 0.04, 0.08, 0.12):
        sig, nh, pi, v = _surface_state(P, -0.505, 0.5 * P.M, 0.454, psi)
        s0 = initial_state(P, sig, v, pi)
        s00 = initial_state(P0, sig, v, pi)
        for dg in (1e-4, 1e-3, 1e-2):
            d = dg / math.sqrt(2.0) * nh            # Dgamma = sqrt2 |dev deps|
            s = run_path(P, s0, [d])[-1]
            pe = plastic(s.sigma, s.pi_i, s.v, P)
            estar = 1 / kn - np.trace(s.eps_e)
            # unconstrained rate path (floor off) + one post Pi_f: the size of F(sigma_f) a post floor produces
            su = run_path(P0, s00, [d])[-1]
            eef, act, dvf, _ = floor_project(su.eps_e, P)
            Ff = plastic(energy(eef, P, False).sig, su.pi_i, su.v, P).F / P.pmin
            pu = np.trace(su.sigma) / 3
            print(f"  psi_i {psi:+.2f} Dgamma {dg:.0e}: rate status {s.flags['status']} modes "
                  f"{sorted(set(sg['mode'] for sg in s.flags['segments']))}, p_end = {pe.p:.8f}, F/p_min = {pe.F / P.pmin:+.1e}, "
                  f"eps^p_v = {s.eps_p_v:+.2e}, eps^f_v = {s.flags['eps_f_v']:.2e}, eps*_end = {estar:.3e}; unconstrained "
                  f"p_c = {pu:.4f} ({su.flags['status']}), +post Pi_f: F(sigma_f)/p_min = {Ff:+.3f}, eps^p_v = {su.eps_p_v:+.2e}")
    print("  (sheet round 3b, BE split: F(sigma_f)/p_min +0.27 ... +0.82, |Delta eps^p_v| <= 2.7e-5 < eps*_f 7.5e-5; Adversary"
          " +0.5 ... +1.35, 3.5e-5)")


def check_floor_rate_record():
    hdr("floor_rate.py record (S.55, BA06 K2 rho 0.7/0.8, p_min 50, eta 0.6M, wet side, both mechanisms)")
    ref = {0.0: -41.296567, 5.0: -41.384056}
    for a0 in (0.0, 5.0):
        P = Params(rho=0.7, rho_bar=0.8, alpha0=a0, p_min=50.0)
        xi = np.array([-1.2, -0.1, 1.3])
        nhv = (xi - xi.mean()) / np.linalg.norm(xi - xi.mean())
        th = math.acos(max(-1.0, min(1.0, SQ6 * float(np.sum(nhv ** 3))))) / 3
        z, _ = zeta_fun(th, math.cos(3 * th) / SQ6, P.rho, P.zeta)
        q = 0.6 * P.M * 50.0 / z
        sig = np.diag(-50.0 + SQ23 * q * nhv)
        pi = pi_of_eta(P, -50.0, 0.6 * P.M)
        s0 = initial_state(P, sig, 1.70, pi)
        d = np.diag(nhv + 0.02) * 2.0e-3
        s = run_path(P, s0, [d])[-1]
        pe = plastic(s.sigma, s.pi_i, s.v, P)
        print(f"  alpha0 {a0}: status {s.flags['status']}, modes {[sg['mode'] for sg in s.flags['segments']]}, p_end = "
              f"{pe.p:.8f}, pi_end = {s.pi_i:.6f} (floor_rate.py {ref[a0]:.6f}, diff {s.pi_i - ref[a0]:+.1e}), min lam_f'/|deps| = "
              f"{s.flags['min_lam_f_rel']:.3e}, F/(M|p|) = {pe.F / (P.M * abs(pe.p)):+.1e}, n_f_init {s0.flags['n_f_init']}")
        if a0 == 0.0:
            el = energy(s0.eps_e, P)
            from .model import PFLOOR
            PA = np.einsum("ij,ijkl->kl", PFLOOR, el.a4)
            print(f"      floor-only lam_f for an isotropic rate: {np.einsum('kl,kl->', PA, np.eye(3) / 3) / np.einsum('kl,kl->', PA, PFLOOR):.15f}"
                  f" (= tr = 1), P:a^e:P = {np.einsum('kl,kl->', PA, PFLOOR):.6f} vs K = {-el.p / P.kappa_hat:.6f}")


def check_floor_inactive():
    hdr("Floor on but inactive (p_min 'default' on K2 paths): same answer as floor off to the ODE tolerance")
    P0 = Params(rho=0.7, rho_bar=0.8)
    P1 = Params(rho=0.7, rho_bar=0.8, p_min="default")
    for kind, ax in (("drained", -0.2), ("undrained", -0.2), ("drained", 0.05)):
        a = triaxial(P0, initial_state(P0, SIG0, 1.59, -60.4), kind, ax, 40)[-1]
        b = triaxial(P1, initial_state(P1, SIG0, 1.59, -60.4), kind, ax, 40)[-1]
        print(f"  {kind} {ax:+.2f}: |sigma_on - sigma_off|/|sigma| = {np.linalg.norm(a.sigma - b.sigma) / np.linalg.norm(a.sigma):.1e},"
              f" pi_i {rel(b.pi_i, a.pi_i):.1e}; floor segments {b.flags['floor_segments']}, eps^f_v {b.flags['eps_f_v']}")


def check_pi0_rule():
    hdr("K1.15 unified pi_i0 rule (S.53), K2 set p_init = -100")
    Ps = Params(rho=1.0, rho_bar=1.0, cap="smooth", c1=0.05, c2=0.15)
    Pn = Params(rho=1.0, rho_bar=1.0)
    Pp = Params(rho=1.0, rho_bar=1.0, cap="planar", c1=0.1, c2=0.1)
    s_txc = lambda q: np.diag([-100.0 - 2 * q / 3, -100.0 + q / 3, -100.0 + q / 3])
    rows = [("smooth, isotropic (eta* = c2 M = 0.18)", Ps, SIG0, -50.995881),
            ("none, isotropic (eta* = 0, apex)", Pn, SIG0, -46.475800),
            ("smooth, eta_init 0.75", Ps, s_txc(75.0), -71.554175),
            ("smooth, eta_init = M", Ps, s_txc(120.0), -100.0),
            ("planar chi_cap 0.10, isotropic (eta* = 0.12)", Pp, SIG0, None)]
    for lab, P, sig, ex in rows:
        st = initial_state(P, sig, 1.59, None, pi_rule="S53")
        sr = initial_state(P, sig, 1.59, None)
        exs = f"sheet {ex:.6f}, rel {rel(st.pi_i, ex):.1e}" if ex is not None else \
            f"closed {pi_of_eta(P, -100.0, 0.12):.6f}"
        print(f"  {lab}: pi_i0 = {st.pi_i:.6f} ({exs}); F/(M|p|) = {st.flags['F_rel']:+.2e}; old 'surface' rule {sr.pi_i:.6f}")
    try:
        initial_state(Params(N=0.4, cap="smooth", c1=0.05, c2=0.15), s_txc(310.0), 1.59, None, pi_rule="S53")
        print("  eta* >= M/N: NOT refused  <-- FAIL")
    except ValueError as e:
        print(f"  eta* >= M/N: refused ({e})")
    # first yield: drained TXC from the isotropic state, and a constant-p TXC path
    P = Ps
    s0 = initial_state(P, SIG0, 1.59, None, pi_rule="S53")
    for lab, kind in (("drained TXC", "drained"), ("constant-p TXC (isochoric, BA06 alpha0 = 0)", "undrained")):
        da = -2e-3
        if kind == "drained":
            st = integrate_increment(P, s0, np.diag([da, 0.0, 0.0]), smask=[False, True, True, False, False, False])
        else:
            st = integrate_increment(P, s0, np.diag([da, -da / 2, -da / 2]))
        sg = st.flags["segments"][0]
        assert sg["event"] == "yield", sg
        f = sg["t1"]
        if kind == "drained":
            sy = integrate_increment(P, s0, np.diag([f * da, 0.0, 0.0]), smask=[False, True, True, False, False, False])
        else:
            sy = integrate_increment(P, s0, np.diag([f * da, -f * da / 2, -f * da / 2]))
        pq = plastic(sy.sigma, sy.pi_i, sy.v, P)
        print(f"  first yield on {lab}: q = {pq.q:.4f}, p = {pq.p:.4f}, eta = q/|p| = {pq.q / -pq.p:.4f}, w = {pq.w:.3f}, "
              f"F/(M|p|) = {pq.F / (P.M * abs(pq.p)):+.0e}")
    print("  (sheet: drained TXC q = 11.3469, p = -103.7823, eta = 0.1093, w = 0.337; constant p: eta = 0.18 exactly)")
    # initialState with the floor: Pi_f first (counted), then (S.53) at the floored p_init (HAR: q scaled by (varpi_f/varpi)^n)
    for pa in PA_SET:
        P = tims_params(pa, cap="smooth", c1=0.05, c2=0.15)
        sig, _ = sig_at(-0.2, 0.1, math.pi / 3)
        st = initial_state(P, sig, 1.75, None, pi_rule="S53")
        pq = plastic(st.sigma, st.pi_i, st.v, P)
        kn = P.k * (1 - P.n_e)
        vp0 = math.sqrt(0.04 + kn * 0.01 / (3 * P.g))
        vpf = math.sqrt(pq.p ** 2 + kn * pq.q ** 2 / (3 * P.g))
        print(f"  p_a {pa}: floored init from p = -0.2, q = 0.1: n_f_init {st.flags['n_f_init']}, p = {pq.p:.12f} (-p_min "
              f"{-P.pmin}), q = {pq.q:.8f} (0.1 (varpi_f/varpi)^n = {0.1 * (vpf / vp0) ** 0.5:.8f}), eps^f_v = "
              f"{st.flags['eps_f_v']:.6e}, W_f = {st.flags['W_f']:.4e}, eta_init = zeta q/|p| = {pq.zeta * pq.q / -pq.p:.4f} vs "
              f"c2 M = {0.15 * P.M:.4f} -> pi_i0 = {st.pi_i:.6f} (closed (S.53) at the floored p: "
              f"{pi_of_eta(P, -P.pmin, max(pq.zeta * pq.q / -pq.p, 0.15 * P.M)):.6f}); F/(M|p|) = {st.flags['F_rel']:+.1e}")
        sig, _ = sig_at(-0.2, 0.02, math.pi / 3)
        st = initial_state(P, sig, 1.75, None, pi_rule="S53")
        pq = plastic(st.sigma, st.pi_i, st.v, P)
        print(f"  p_a {pa}: floored init from p = -0.2, q = 0.02: eta_init = {pq.zeta * pq.q / -pq.p:.4f} < c2 M -> pi_i0 = "
              f"{st.pi_i:.6f} (closed {pi_of_eta(P, -P.pmin, 0.15 * P.M):.6f}), F/(M|p|) = {st.flags['F_rel']:+.2e} (inside)")


GROUPS = dict(
    basic=[check_zeta, check_elastic_tangent, check_iso, check_loop, check_refusal],
    tangents=[check_finite_tangent, check_cont_tangent_fd],
    dissipation=[check_dissipation_and_peak],
    rtol=[check_rtol],
    caps=[check_caps],
    har=[check_har_k1, check_har_tangent_loop, check_har_parser],
    floor=[check_floor_k1_12, check_floor_k1_13_14, check_floor_fpf, check_floor_wet, check_floor_rate_record,
           check_floor_inactive],
    pi0=[check_pi0_rule],
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
