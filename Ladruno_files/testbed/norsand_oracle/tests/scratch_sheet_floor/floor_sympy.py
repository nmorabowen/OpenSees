"""Sheet round 2026-10-03 (p-floor, pi_i0 rule): symbolic / 30-digit checks for the new sheet items.

  (a) HAR floor closed form (sheet S.50): at fixed eps_s, p(eps_v,f, eps_s) = -p_min exactly for n = 1/2 (closed
      form) and for general n (the scalar equation x^2 - a x^(2n) - b = 0 is equivalent to p = -p_min; unique root).
  (b) eps'_v,f = d eps_v,f / d eps_s = -D12/D11 at the floor (IFT) and its closed forms (HAR, BA06 with alpha0).
  (c) BA06 floor closed form with alpha0 (sheet S.49): p = -p_min exactly.
  (d) The floored tangent C_f = a^e(eps_f) : DPi_f (sheet S.51): equals the 30-digit central FD of sigma(Pi_f(eps))
      over the three principal strains (non-degenerate state), and delta : C_f = 0 (dp = 0), HAR and BA06 (alpha0 != 0).
  (e) Energy injected by the floor, 0 <= Psi(eps_f) - Psi(eps) <= p_min * deps^f_v (sheet S.52): symbolic for BA06
      alpha0 = 0 (1 - t <= -ln t), 30-digit grid for HAR (TIMs constants) and for BA06 with alpha0 = 5 and 0.
  (f) Unified pi_i0 rule (sheet S.53): eta(p, pi_i0) = eta* for both N branches; the N -> 0 limit of the power branch
      is the exp branch; the eta* = M and eta* = 0 limits; monotonicity in eta*.
  (g) K1 floor numbers (sheet 13.12-13.14) at 15 digits: TIMs HAR set, K2 BA06 set.
Run:  python floor_sympy.py   (local py312 venv; sympy + mpmath; a few minutes)
"""
import random

import mpmath as mp
import sympy as sp

mp.mp.dps = 30
R = sp.Rational
ok_all = True


def report(name, val, tol=sp.Float("1e-28")):
    global ok_all
    if isinstance(val, bool):
        good = val
        shown = val
    else:
        v = abs(sp.N(val, 30))
        good = bool(v <= tol)
        shown = sp.N(v, 3)
    ok_all &= bool(good)
    print(f"  [{'OK' if good else 'FAIL'}] {name}: {shown}")


# ---------------------------------------------------------------------------------------------------------------
# HAR energy (S.4h)-(S.5h') in sympy, sheet signs (p < 0)
# ---------------------------------------------------------------------------------------------------------------
k, g, n, pa, ev, es, pmin = sp.symbols("k g n p_a eps_v eps_s p_min", positive=True)


def har(ev_, es_, k_, g_, n_, pa_):
    est = 1 / (k_ * (1 - n_)) - ev_
    u = sp.sqrt(est ** 2 + 3 * g_ * es_ ** 2 / (k_ * (1 - n_)))
    Psi = pa_ / (k_ * (2 - n_)) * (k_ * (1 - n_) * u) ** ((2 - n_) / (1 - n_))
    w = (k_ * (1 - n_) * u) ** (n_ / (1 - n_))
    p = -pa_ * k_ * (1 - n_) * est * w
    q = 3 * g_ * pa_ * es_ * w
    return Psi, p, q


Psi_h, p_h, q_h = har(ev, es, k, g, n, pa)
D11_h = sp.diff(p_h, ev)
D12_h = sp.diff(p_h, es)

print("=== (a) HAR floor closed form at fixed eps_s ===")
x = sp.symbols("x", positive=True)
a_ = 3 * k * (1 - n) * g * es ** 2
b_ = (pmin / pa) ** 2
q_f = 3 * g * pa * es * x ** n
varpi2 = pmin ** 2 + k * (1 - n) * q_f ** 2 / (3 * g)
report("general n: varpi^2 - p_a^2 x^2 == p_a^2 (b + a x^(2n) - x^2)",
       sp.simplify(varpi2 - pa ** 2 * x ** 2 - pa ** 2 * (b_ + a_ * x ** (2 * n) - x ** 2)))
half = R(1, 2)
a_h = a_.subs(n, half)
x_half = (a_h + sp.sqrt(a_h ** 2 + 4 * b_)) / 2
report("n = 1/2: x^2 - a x - b == 0 at the closed form", sp.simplify(x_half ** 2 - a_h * x_half - b_))
est_f = (pmin / pa) / (k * (1 - n) * x ** n)
ev_f = 1 / (k * (1 - n)) - est_f
random.seed(3)
worst = sp.Float(0, 30)
for trial in range(6):
    kk, gg, ppa, ppm, ees = [sp.Float(v, 30) for v in (random.uniform(300, 3000), random.uniform(100, 1500),
                                                          random.uniform(50, 200), random.uniform(0.05, 2.0),
                                                          random.uniform(1e-6, 3e-3))]
    for nn in (half, R(3, 10), R(7, 10)):
        subs = {k: kk, g: gg, pa: ppa, pmin: ppm, es: ees, n: nn}
        if nn == half:
            xv = sp.N(x_half.subs(subs), 30)
        else:
            f = (x ** 2 - a_ * x ** (2 * nn) - b_).subs(subs)
            x0 = sp.N(sp.sqrt(b_.subs(subs)) + a_.subs(subs) + 1, 30)
            xv = sp.nsolve(f, x, x0, prec=30)
        evf = sp.N(ev_f.subs(subs).subs(x, xv), 30)
        pv = sp.N(p_h.subs(subs).subs(ev, evf), 30)
        worst = max(worst, abs(pv + ppm) / ppm)
report("p(eps_v,f, eps_s) + p_min = 0, 6 random states x n in {1/2, 0.3, 0.7} (rel; 30-digit nsolve floor)", worst, sp.Float("1e-26"))
# uniqueness of the root of f(x) = x^2 - a x^(2n) - b: f(0) = -b < 0, f' = 0 only at x_s^(2-2n) = n a, where
# f(x_s) = -(1-n) a x_s^(2n) - b < 0; f is decreasing on [0, x_s] and increasing after, so exactly one root, > x_s.
xs = (n * a_) ** (1 / (2 - 2 * n))
f_xs = xs ** 2 - a_ * xs ** (2 * n) - b_
worst = sp.Float(0, 30)
for nn in (half, R(3, 10), R(7, 10), R(9, 10)):
    s2 = {k: sp.Float(1889.48, 30), g: sp.Float(807.8, 30), pa: sp.Float(101, 30), pmin: sp.Float("0.505", 30),
          es: sp.Float("1e-3", 30), n: nn}
    worst = max(worst, abs(sp.N((f_xs - (-(1 - n) * a_ * xs ** (2 * n) - b_)).subs(s2), 30)))
report("f(x_s) == -(1-n) a x_s^(2n) - b (< 0) at n = 1/2, 0.3, 0.7, 0.9 (30 digits)", worst, sp.Float("1e-26"))

print("=== (b) eps'_v,f = -D12/D11 at the floor ===")
varpi2_f = pa ** 2 * x ** 2
epsp_h_closed = n * pmin * q_f / ((1 - n) * varpi2_f + n * pmin ** 2)
D11_sf = k * pa * x ** n * (1 - n + n * pmin ** 2 / varpi2_f)
D12_sf = n * k * (-pmin) * q_f * pa * x ** n / varpi2_f
report("HAR: closed form == -D12/D11 (stress form, p = -p_min)", sp.simplify(epsp_h_closed - (-D12_sf / D11_sf)))
subs = {k: sp.Float("1889.48104361", 30), g: sp.Float("807.80387674", 30), pa: sp.Float(101, 30),
        pmin: sp.Float("0.505", 30), n: half}


def evf_of(es_val):
    s2 = dict(subs)
    s2[es] = es_val
    return sp.N(ev_f.subs(s2).subs(x, x_half.subs(s2)), 30)


e0 = sp.Float("2e-4", 30)
h = sp.Float("1e-11", 30)      # eps*_f ~ 1/eps_s near the floor: FD truncation ~ h^2/eps_s^2, so h must be tiny (30 digits allow it)
fd = (evf_of(e0 + h) - evf_of(e0 - h)) / (2 * h)
s3 = dict(subs)
s3[es] = e0
s3[ev] = evf_of(e0)
ift = sp.N((-D12_h / D11_h).subs(s3), 30)
report("HAR TIMs: FD d eps_v,f/d eps_s vs -D12/D11 (strain-form Hessian), rel", (fd - ift) / ift, sp.Float("1e-14"))
p0, kap, ev0, mu0, al0 = sp.symbols("p_0 kappa eps_v0 mu_0 alpha_0", real=True)
om = -(ev - ev0) / kap
Pt = -p0 * kap * sp.exp(om)
Psi_b = Pt + R(3, 2) * (mu0 + al0 / kap * Pt) * es ** 2
p_b = sp.diff(Psi_b, ev)
q_b = sp.diff(Psi_b, es)
D11_b = sp.diff(p_b, ev)
D12_b = sp.diff(p_b, es)
print("=== (c) BA06 floor closed form with alpha0 ===")
ev_fb = ev0 - kap * sp.log(pmin / (-p0 * (1 + 3 * al0 * es ** 2 / (2 * kap))))
report("BA06: p(eps_v,f, eps_s) + p_min == 0 (symbolic)", sp.simplify(p_b.subs(ev, ev_fb) + pmin))
epsp_b_closed = 3 * al0 * es / (1 + 3 * al0 * es ** 2 / (2 * kap))
report("BA06: eps' closed form == -D12/D11 at the floor (symbolic)",
       sp.simplify(epsp_b_closed - (-D12_b / D11_b).subs(ev, ev_fb)))
report("BA06 alpha0 = 0: eps' = 0", sp.simplify(epsp_b_closed.subs(al0, 0)))

print("=== (d) floored tangent C_f = a^e(eps_f) : DPi_f vs 30-digit FD; delta : C_f = 0 ===")
e1, e2, e3 = sp.symbols("e1 e2 e3", real=True)


def principal_floor_check(sig_of_eps, ev_f_of_es, epsp_of_es, pmin_val, label):
    E = sp.Matrix([e1, e2, e3])
    evv = e1 + e2 + e3
    dev = E - evv / 3 * sp.ones(3, 1)
    ne = sp.sqrt(dev.dot(dev))
    ess = sp.sqrt(R(2, 3)) * ne
    Ef = E - (evv - ev_f_of_es(ess)) / 3 * sp.ones(3, 1)
    sig_f = sp.Matrix(sig_of_eps(*Ef))
    pt = {e1: sp.Float("-3.1e-4", 30), e2: sp.Float("4.7e-4", 30), e3: sp.Float("2.2e-4", 30)}
    hh = sp.Float("1e-12", 30)   # same truncation argument as in (b)
    C_fd = sp.zeros(3, 3)
    for b, eb in enumerate((e1, e2, e3)):
        pp = dict(pt)
        pp[eb] = pt[eb] + hh
        pm = dict(pt)
        pm[eb] = pt[eb] - hh
        C_fd[:, b] = (sp.N(sig_f.subs(pp), 30) - sp.N(sig_f.subs(pm), 30)) / (2 * hh)
    Ef_num = sp.N(Ef.subs(pt), 30)
    sig_plain = sp.Matrix(sig_of_eps(e1, e2, e3))
    ae = sp.N(sig_plain.jacobian([e1, e2, e3]).subs({e1: Ef_num[0], e2: Ef_num[1], e3: Ef_num[2]}), 30)
    es_num = sp.N(ess.subs(pt), 30)
    nh = sp.N((dev / ne).subs(pt), 30)
    epsp = sp.N(epsp_of_es(es_num), 30)
    Phi = sp.zeros(3, 3)
    for a in range(3):
        for b in range(3):
            Phi[a, b] = (1 if a == b else 0) - R(1, 3) + R(1, 3) * epsp * sp.sqrt(R(2, 3)) * nh[b]
    C_cl = ae * Phi
    cmax = max(abs(C_cl[i, j]) for i in range(3) for j in range(3))
    err = max(abs(C_fd[i, j] - C_cl[i, j]) for i in range(3) for j in range(3)) / cmax
    report(f"{label}: C_f closed form vs FD (rel)", err, sp.Float("1e-15"))
    vol = max(abs(sum(C_cl[i, j] for i in range(3))) for j in range(3)) / cmax
    report(f"{label}: delta : C_f = 0 (rel)", vol, sp.Float("1e-26"))
    pf = sp.N(sum(sig_f.subs(pt)) / 3, 30)
    report(f"{label}: p(Pi_f(eps)) + p_min = 0", pf + pmin_val, sp.Float("1e-25"))
    names = (e1, e2, e3)
    sp_ = max(abs((Ef_num[a] - Ef_num[b]) / (pt[names[a]] - pt[names[b]]) - 1) for a in range(3) for b in range(3) if a != b)
    report(f"{label}: unit spin of Pi_f", sp_, sp.Float("1e-25"))
    # the mutant with eps' dropped: delta : C != 0 when D12 != 0
    Phi_m = sp.zeros(3, 3)
    for a in range(3):
        for b in range(3):
            Phi_m[a, b] = (1 if a == b else 0) - R(1, 3)
    C_m = ae * Phi_m
    print(f"      mutant eps' dropped: delta:C_mutant/max|C| = {sp.N(max(abs(sum(C_m[i, j] for i in range(3))) for j in range(3)) / cmax, 3)}"
          f" (0 only when D12 = 0)")


def har_sig(a1, a2, a3):
    evv = a1 + a2 + a3
    dev = [a1 - evv / 3, a2 - evv / 3, a3 - evv / 3]
    ne = sp.sqrt(dev[0] ** 2 + dev[1] ** 2 + dev[2] ** 2)
    ess = sp.sqrt(R(2, 3)) * ne
    _, pp, qq = har(evv, ess, subs[k], subs[g], half, subs[pa])
    return [pp + sp.sqrt(R(2, 3)) * qq * dev[i] / ne for i in range(3)]


def har_evf(ess):
    s2 = dict(subs)
    s2[es] = ess
    return ev_f.subs(s2).subs(x, x_half.subs(s2))


def har_epsp(ess):
    s2 = dict(subs)
    s2[es] = ess
    return epsp_h_closed.subs(s2).subs(x, x_half.subs(s2))


principal_floor_check(har_sig, har_evf, har_epsp, subs[pmin], "HAR n=1/2 TIMs")

sb = {p0: sp.Float(-100, 30), kap: sp.Float("0.01", 30), ev0: 0, mu0: sp.Float(5400, 30), al0: sp.Float(5, 30),
      pmin: sp.Float("0.5", 30)}


def ba_sig(a1, a2, a3):
    evv = a1 + a2 + a3
    dev = [a1 - evv / 3, a2 - evv / 3, a3 - evv / 3]
    ne = sp.sqrt(dev[0] ** 2 + dev[1] ** 2 + dev[2] ** 2)
    ess = sp.sqrt(R(2, 3)) * ne
    pp = p_b.subs(sb).subs({ev: evv, es: ess})
    qq = q_b.subs(sb).subs({ev: evv, es: ess})
    return [pp + sp.sqrt(R(2, 3)) * qq * dev[i] / ne for i in range(3)]


principal_floor_check(ba_sig, lambda ess: ev_fb.subs(sb).subs(es, ess), lambda ess: epsp_b_closed.subs(sb).subs(es, ess),
                      sb[pmin], "BA06 alpha0=5 K2")

print("=== (e) floor energy injection: 0 <= Psi(eps_f) - Psi(eps) <= p_min deps^f_v ===")
t = sp.symbols("t", positive=True)
Ef_b0 = Psi_b.subs(al0, 0).subs(ev, ev_fb.subs(al0, 0)) - Psi_b.subs(al0, 0)
Ef_b0_t = sp.simplify(Ef_b0.subs(ev, ev0 - kap * sp.log(pmin * t / (-p0))))
report("BA06 a0=0: Psi(eps_f) - Psi(eps) == kappa p_min (1 - t), t = |p|/p_min", sp.simplify(Ef_b0_t - kap * pmin * (1 - t)))
dfv = sp.simplify((ev0 - kap * sp.log(pmin * t / (-p0))) - ev_fb.subs(al0, 0))
report("BA06 a0=0: deps^f_v == kappa ln(1/t)", sp.simplify(dfv - kap * sp.log(1 / t)))
report("1 - t <= -ln t on (0,1): d/dt[-ln t - (1-t)] == 1 - 1/t (< 0), value 0 at t = 1",
       sp.simplify(sp.diff(-sp.log(t) - (1 - t), t) - (1 - 1 / t)))


def grid_energy(Psi_expr, p_expr, ev_f_expr, subsd, label, es_vals, pfracs, har_mode):
    worst_lo, worst_hi = sp.Float(1, 30), sp.Float(1, 30)
    for esv in es_vals:
        for pv in pfracs:
            s2 = dict(subsd)
            s2[es] = esv
            target = -pv * s2[pmin]
            if har_mode:
                evv = sp.nsolve((p_expr - target).subs(s2), ev, sp.Float("1e-4", 30), prec=30)
            else:   # BA06 closed form at fixed eps_s: the same log form as (S.49) with p_min -> |target|
                evv = sp.N((ev0 - kap * sp.log(-target / (-p0 * (1 + 3 * al0 * es ** 2 / (2 * kap))))).subs(s2), 30)
            evf = sp.N(ev_f_expr.subs(s2).subs(x, x_half.subs(s2)), 30) if har_mode else sp.N(ev_f_expr.subs(s2), 30)
            dE = sp.N(Psi_expr.subs(s2).subs(ev, evf) - Psi_expr.subs(s2).subs(ev, evv), 30)
            bound = sp.N(s2[pmin] * (evv - evf), 30)
            worst_lo = min(worst_lo, dE / bound)
            worst_hi = min(worst_hi, (bound - dE) / bound)
    print(f"  {label}: min E_f/(p_min deps^f_v) = {sp.N(worst_lo, 6)}, min (1 - E_f/bound) = {sp.N(worst_hi, 6)}")
    report(f"{label}: 0 <= E_f", bool(worst_lo > 0))
    report(f"{label}: E_f <= p_min deps^f_v", bool(worst_hi > 0))


ess_vals = [sp.Float(v, 30) for v in ("0", "1e-4", "5e-4", "2e-3")]
pfrac = [sp.Float(v, 30) for v in ("0.01", "0.1", "0.5", "0.9", "0.999")]
grid_energy(Psi_h.subs(n, half), p_h.subs(n, half), ev_f.subs(n, half), subs, "HAR TIMs n=1/2 p_min=0.505", ess_vals, pfrac, True)
grid_energy(Psi_b, p_b, ev_fb, sb, "BA06 alpha0=5 K2 p_min=0.5", ess_vals, pfrac, False)
sb0 = dict(sb)
sb0[al0] = 0
grid_energy(Psi_b, p_b, ev_fb, sb0, "BA06 alpha0=0 K2 p_min=0.5", ess_vals, pfrac, False)

print("=== (f) unified pi_i0 rule (S.53) ===")
M, N, pp, eta_s = sp.symbols("M N p eta_star", positive=True)
pi_pow = -pp * ((1 - N) / (1 - eta_s * N / M)) ** ((1 - N) / N)
pi_exp = -pp * sp.exp(eta_s / M - 1)
eta_of = lambda p_, pi_: (M / N) * (1 - (1 - N) * (p_ / pi_) ** (N / (1 - N)))
eta_of0 = lambda p_, pi_: M * (1 + sp.log(pi_ / p_))
# nested fractional powers do not simplify symbolically without sign knowledge of 1 - eta* N/M; 30-digit identities
# at random admissible data (eta* < M/N) are the check, plus the exact special values.
random.seed(11)
w1 = w2 = w3 = w4 = sp.Float(0, 30)
for trial in range(10):
    Mv = sp.Float(random.uniform(0.8, 1.6), 30)
    Nv = sp.Float(random.uniform(0.05, 0.6), 30)
    pv_ = sp.Float(random.uniform(0.3, 300), 30)
    ev_ = sp.Float(random.uniform(0.0, 0.95), 30) * Mv / Nv          # eta* in [0, 0.95 M/N)
    sub = {M: Mv, N: Nv, pp: pv_, eta_s: ev_}
    w1 = max(w1, abs(sp.N((eta_of(-pp, pi_pow) - eta_s).subs(sub), 30)) / ev_ if ev_ != 0 else 0)
    w2 = max(w2, abs(sp.N((eta_of0(-pp, pi_exp) - eta_s).subs(sub), 30)) / ev_ if ev_ != 0 else 0)
    # N -> 0: power branch at N = 1e-9 vs exp branch (difference O(N))
    sub0 = dict(sub)
    sub0[N] = sp.Float("1e-9", 30)
    sub0[eta_s] = sp.Float(random.uniform(0.0, 2.0), 30) * Mv
    w3 = max(w3, abs(sp.N((sp.log(pi_pow / (-pp)) - (eta_s / M - 1)).subs(sub0), 30)))
    # monotonicity: d ln|pi_i0|/d eta* == (1-N)/(M - eta* N) (central FD at 30 digits)
    hh = sp.Float("1e-10", 30)
    sp_, sm_ = dict(sub), dict(sub)
    sp_[eta_s] = ev_ + hh
    sm_[eta_s] = ev_ - hh
    fdv = (sp.N(sp.log(-pi_pow).subs(sp_), 30) - sp.N(sp.log(-pi_pow).subs(sm_), 30)) / (2 * hh)
    w4 = max(w4, abs(fdv - sp.N(((1 - N) / (M - eta_s * N)).subs(sub), 30)) / abs(fdv))
report("N > 0: eta(p, pi_i0) == eta*, 10 random sets (rel)", w1, sp.Float("1e-27"))
report("N = 0: eta(p, pi_i0) == eta*, 10 random sets (rel)", w2, sp.Float("1e-27"))
report("N -> 0: power branch at N = 1e-9 vs exp branch (log form, O(N))", w3, sp.Float("1e-7"))
report("eta* = M: pi_i0 = p", sp.simplify(pi_pow.subs(eta_s, M) + pp))
report("eta* = 0: pi_i0 = p (1-N)^((1-N)/N)", sp.simplify(pi_pow.subs(eta_s, 0) + pp * (1 - N) ** ((1 - N) / N)))
report("d ln|pi_i0|/d eta* == (1-N)/(M - eta* N) (FD, rel)", w4, sp.Float("1e-14"))

print("=== (g) K1 floor numbers ===")
mp.mp.dps = 25
kk, gg, ppa, nn, ppm = (mp.mpf("1889.48104361"), mp.mpf("807.80387674"), mp.mpf(101), mp.mpf("0.5"), mp.mpf("0.505"))
evf0 = (1 - (ppm / ppa) ** (1 - nn)) / (kk * (1 - nn))
G_f = gg * ppa * (ppm / ppa) ** nn
K_f = kk * ppa * (ppm / ppa) ** nn
print(f"  HAR TIMs: eps_v,f(eps_s = 0) = {mp.nstr(evf0, 15)}  (domain edge 1/(k(1-n)) = {mp.nstr(1/(kk*(1-nn)), 15)})")
print(f"  HAR TIMs: G(p_min) = {mp.nstr(G_f, 15)} kPa, K(p_min) = {mp.nstr(K_f, 15)} kPa; floored shear stiffness 2G = {mp.nstr(2*G_f, 15)}")
evA = (1 - (mp.mpf(1) / ppa) ** (1 - nn)) / (kk * (1 - nn))
evT = evA + mp.mpf("1e-4")
print(f"  HAR TIMs: isotropic p = -1 kPa: eps_v = {mp.nstr(evA, 15)}; trial eps_v = eps_v + 1e-4 = {mp.nstr(evT, 15)} "
      f"({'OUT of domain' if evT >= 1/(kk*(1-nn)) else 'in domain'}); floored eps_v,f = {mp.nstr(evf0, 15)}; "
      f"deps^f_v = {mp.nstr(evT - evf0, 15)}; W_f = p_min deps^f_v = {mp.nstr(ppm*(evT-evf0), 15)} kPa")
ees = mp.mpf("2e-4")
aa = 3 * kk * (1 - nn) * gg * ees ** 2
bb = (ppm / ppa) ** 2
xx = (aa + mp.sqrt(aa ** 2 + 4 * bb)) / 2
estf = (ppm / ppa) / (kk * (1 - nn) * xx ** nn)
evf_s = 1 / (kk * (1 - nn)) - estf
qf = 3 * gg * ppa * ees * xx ** nn
epsp = nn * ppm * qf / ((1 - nn) * (ppa * xx) ** 2 + nn * ppm ** 2)
print(f"  HAR TIMs: eps_s = 2e-4 at the floor: x = varpi_f/p_a = {mp.nstr(xx, 15)}, eps_v,f = {mp.nstr(evf_s, 15)}, "
      f"q_f = {mp.nstr(qf, 15)} kPa, eta_f = q_f/p_min = {mp.nstr(qf/ppm, 15)}, eps'_v,f = {mp.nstr(epsp, 15)}")
kap_, p0_, mu0_, pm_ = mp.mpf("0.01"), mp.mpf(100), mp.mpf(5400), mp.mpf("0.5")
evf_b = -kap_ * mp.log(pm_ / p0_)
evt_b = -kap_ * mp.log(pm_ / 2 / p0_)
print(f"  BA06 K2: eps_v,f = {mp.nstr(evf_b, 15)}; trial at p = -p_min/2: eps_v,tr = {mp.nstr(evt_b, 15)}; "
      f"deps^f_v = kappa ln 2 = {mp.nstr(evt_b - evf_b, 15)}; W_f = {mp.nstr(pm_*(evt_b-evf_b), 15)} kPa; "
      f"E_f = kappa (p_min - |p_tr|) = {mp.nstr(kap_*pm_/2, 15)} kPa; floored shear stiffness 2 mu0 = {mp.nstr(2*mu0_, 10)} kPa")
print("\nALL OK" if ok_all else "\nSOME CHECK FAILED")
