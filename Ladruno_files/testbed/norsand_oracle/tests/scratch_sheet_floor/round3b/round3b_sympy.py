"""Round 3b (2026-10-03) symbolic checks for the five Adversary items.

  (a) the HAR stress forms used by har_patch.py ((S.5h') D11, D12, D22, q/eps_s) against direct differentiation of
      Psi_HAR (S.4h), symbolic and at the TIMs constants (30 digits) -- establishes the patch before it is trusted;
      the (S.5h'') inverse map round trip; the K1.13 / K1.14 numbers reproduced from the patch's closed forms.
  (b) (S.56): W_ramp as written in words (1 - pi(eta2)/pi(eta1)) is NEGATIVE; 1 - pi(eta1)/pi(eta2) equals the
      displayed formula [(1 - c2 N)/(1 - c1 N)]^((1-N)/N) and the N = 0 limit 1 - exp(-(c2 - c1)); K2 values.
  (c) (S.56) at c1 = c2 (planar) and cap = none: W_ramp = 0 identically, so an ungated 'PI_SCAN_REL <= W_ramp/10'
      refuses every planar/none model.
  (d) the HAR floor margin eps*_f(eps_s) = (p_min/p_a)/(k(1-n) x^n): the eps_s = 0 value sqrt(p_min/p_a)/(k(1-n))
      (= 7.49e-5 at TIMs) and its decrease with eps_s (K1.14 state: 1.75e-5).
Run:  python round3b_sympy.py
"""
import mpmath as mp
import sympy as sp

mp.mp.dps = 30
R = sp.Rational
ok_all = True


def report(name, val, tol=sp.Float("1e-28")):
    global ok_all
    if isinstance(val, bool):
        good, shown = val, val
    else:
        v = abs(sp.N(val, 30))
        good = bool(v <= tol)
        shown = sp.N(v, 3)
    ok_all &= bool(good)
    print(f"  [{'OK' if good else 'FAIL'}] {name}: {shown}")


print("=== (a) HAR stress forms of har_patch.py vs direct differentiation of Psi_HAR ===")
k, g, n, pa, ev, es = sp.symbols("k g n p_a eps_v eps_s", positive=True)
est = 1 / (k * (1 - n)) - ev
u = sp.sqrt(est ** 2 + 3 * g * es ** 2 / (k * (1 - n)))
Psi = pa / (k * (2 - n)) * (k * (1 - n) * u) ** ((2 - n) / (1 - n))
p_d = sp.diff(Psi, ev)
q_d = sp.diff(Psi, es)
w = (k * (1 - n) * u) ** (n / (1 - n))
p_s = -pa * k * (1 - n) * est * w
q_s = 3 * g * pa * es * w
report("p = dPsi/deps_v == (S.5h)", sp.simplify(p_d - p_s))
report("q = dPsi/deps_s == (S.5h)", sp.simplify(q_d - q_s))
# stress forms as coded in har_patch.har_D
p, q = sp.symbols("p q", real=True)
varpi2 = p ** 2 + k * (1 - n) * q ** 2 / (3 * g)
varpi = sp.sqrt(varpi2)
fac = pa * (varpi / pa) ** n
Z = varpi2 / p ** 2
D11_code = k * fac * (1 - n + n / Z)
D22_code = (3 * g / (1 - n)) * fac * (1 - n / Z)
D12_code = n * k * p * q * fac / varpi2
ratio_code = 3 * g * fac
subs_T = {k: sp.Float("1889.48104361", 30), g: sp.Float("807.80387674", 30), n: R(1, 2), pa: sp.Float(101, 30)}
worst = sp.Float(0, 30)
for (evv, ess) in ((sp.Float("3.1e-4", 30), sp.Float("2.0e-4", 30)), (sp.Float("9.9e-4", 30), sp.Float("5.0e-5", 30)),
                   (sp.Float("-2.0e-3", 30), sp.Float("1.0e-3", 30))):
    sd = dict(subs_T)
    sd[ev] = evv
    sd[es] = ess
    pv = sp.N(p_s.subs(sd), 30)
    qv = sp.N(q_s.subs(sd), 30)
    sd2 = dict(subs_T)
    sd2[p] = pv
    sd2[q] = qv
    for name_, code, direct in (("D11", D11_code, sp.diff(p_s, ev)), ("D12", D12_code, sp.diff(p_s, es)),
                                ("D21", D12_code, sp.diff(q_s, ev)), ("D22", D22_code, sp.diff(q_s, es)),
                                ("q/eps_s", ratio_code, q_s / es)):
        rel = abs(sp.N(code.subs(sd2), 30) - sp.N(direct.subs(sd), 30)) / abs(sp.N(direct.subs(sd), 30))
        worst = max(worst, rel)
report("D11, D12 = D21, D22, q/eps_s stress forms vs strain-form derivatives, 3 states, TIMs (rel)", worst, sp.Float("1e-27"))
# inverse map (S.5h'')
sd = dict(subs_T)
sd[ev] = sp.Float("3.1e-4", 30)
sd[es] = sp.Float("2.0e-4", 30)
pv, qv = sp.N(p_s.subs(sd), 30), sp.N(q_s.subs(sd), 30)
kn = subs_T[k] * (1 - subs_T[n])
vp = sp.sqrt(pv ** 2 + kn * qv ** 2 / (3 * subs_T[g]))
ev_back = (1 / kn) * (1 - (abs(pv) / subs_T[pa]) ** (1 - subs_T[n]) * (abs(pv) / vp) ** subs_T[n])
es_back = qv / (3 * subs_T[g] * subs_T[pa] * (vp / subs_T[pa]) ** subs_T[n])
report("(S.5h'') inverse map round trip eps_v", ev_back - sd[ev], sp.Float("1e-27"))
report("(S.5h'') inverse map round trip eps_s", es_back - sd[es], sp.Float("1e-27"))
# K1.13 / K1.14 numbers from the (S.50) closed form
pmin = sp.Float("0.505", 30)
kk, gg, ppa = subs_T[k], subs_T[g], subs_T[pa]
edge = 1 / kn
for ess, label, ref_evf in ((sp.Float(0, 30), "K1.13 eps_s = 0", "0.000983645033141868"),
                            (sp.Float("2e-4", 30), "K1.14 eps_s = 2e-4", "0.00104102892667979")):
    a_ = 3 * kn * gg * ess ** 2
    b_ = (pmin / ppa) ** 2
    x = (a_ + sp.sqrt(a_ ** 2 + 4 * b_)) / 2
    est_f = (pmin / ppa) / (kn * sp.sqrt(x))
    evf = edge - est_f
    report(f"{label}: eps_v,f vs sheet", evf - sp.Float(ref_evf, 30), sp.Float("1e-17"))
    print(f"      {label}: eps*_f = edge - eps_v,f = {sp.N(est_f, 6)}  (x = {sp.N(x, 8)})")

print("=== (b) (S.56) W_ramp: words vs formula ===")
N, M, c1, c2, pp = sp.symbols("N M c_1 c_2 p", positive=True)
pi_of_eta = lambda eta: pp * ((1 - N) / (1 - eta * N / M)) ** ((1 - N) / N)     # noqa: E731  (S.12) inverse, N > 0
words = 1 - pi_of_eta(c2 * M) / pi_of_eta(c1 * M)          # as the round-3 sheet wrote it
fixed = 1 - pi_of_eta(c1 * M) / pi_of_eta(c2 * M)          # the intended quantity
formula = 1 - ((1 - c2 * N) / (1 - c1 * N)) ** ((1 - N) / N)
import random
random.seed(7)
worst_b = sp.Float(0, 30)
for _ in range(8):
    pt = {N: sp.Float(random.uniform(0.05, 0.9), 30), c1: sp.Float(random.uniform(0.0, 0.3), 30),
          M: sp.Float(random.uniform(0.8, 1.6), 30), pp: sp.Float(-random.uniform(1, 500), 30)}
    pt[c2] = pt[c1] + sp.Float(random.uniform(0.01, 0.4), 30)
    worst_b = max(worst_b, abs(sp.N((fixed - formula).subs(pt), 30)))
report("1 - pi(eta1)/pi(eta2) == displayed formula [(1-c2N)/(1-c1N)]^((1-N)/N), 8 random sets (30 digits)", worst_b)
K2 = {N: sp.Float("0.4", 30), c1: sp.Float("0.05", 30), c2: sp.Float("0.15", 30), M: sp.Float("1.2", 30), pp: -100}
print(f"      K2 (c1 0.05, c2 0.15, N 0.4): words 1 - pi(eta2)/pi(eta1) = {sp.N(words.subs(K2), 6)}  (NEGATIVE);"
      f"  1 - pi(eta1)/pi(eta2) = {sp.N(fixed.subs(K2), 6)}  (the 0.0606 of the sheet)")
report("words form is negative on K2", bool(sp.N(words.subs(K2)) < 0))
report("fixed form = 0.0606 on K2 (to 4 digits)", sp.N(fixed.subs(K2), 30) - sp.Float("0.0606", 30), sp.Float("5e-5"))
# N -> 0 limit of the formula is 1 - exp(-(c2 - c1))
lim = sp.limit(formula.subs({c1: R(1, 20), c2: R(3, 20)}), N, 0)
report("N -> 0 limit == 1 - exp(-(c2 - c1)) at c1 0.05, c2 0.15", sp.simplify(lim - (1 - sp.exp(-R(1, 10)))))
for c2v, expect in (("0.07", "0.0122"), ("0.06", "0.0061")):
    val = sp.N(fixed.subs({N: sp.Float("0.4", 30), c1: sp.Float("0.05", 30), c2: sp.Float(c2v, 30), M: 1, pp: -1}), 6)
    print(f"      c1 0.05, c2 {c2v}: W_ramp = {val} (sheet {expect})")
# the ramp width measured by floor_fold.py (C): |pi(eta1) - pi(eta2)| / |pi(eta2)| at p = -100
w_meas = abs(pi_of_eta(c1 * M) - pi_of_eta(c2 * M)) / abs(pi_of_eta(c2 * M))
report("floor_fold (C) ratio |pi1 - pi2|/|pi2| == 1 - pi1/pi2 (pi1, pi2 same sign, |pi1| < |pi2|)",
       sp.simplify(w_meas.subs(K2) - fixed.subs(K2)), sp.Float("1e-28"))

print("=== (c) (S.56) at c1 = c2 and at cap = none ===")
report("c1 = c2: W_ramp == 0 identically", sp.simplify(formula.subs(c2, c1)))
print("      cap = none: w == 1, no ramp, W_ramp has no meaning (c2 := 0 in (S.53) only); an ungated refusal "
      "'PI_SCAN_REL <= W_ramp/10' with W_ramp = 0 fails for every planar/none model since PI_SCAN_REL = 1e-3 > 0.")

print("=== (d) HAR floor margin eps*_f(eps_s) ===")
es_, x_ = sp.symbols("eps_s x", positive=True)
print("      eps*_f = (p_min/p_a)/(k(1-n) x^n), x >= sqrt(b) = p_min/p_a with equality at eps_s = 0:")
pm = sp.Symbol("p_min", positive=True)
est_f_sym = (pm / pa) / (k * (1 - n) * x_ ** n)
est0 = est_f_sym.subs(x_, pm / pa)
tv = {pm: pmin, pa: ppa, k: kk, n: R(1, 2)}
report("eps_s = 0: eps*_f == (p_min/p_a)^(1-n)/(k(1-n)) (TIMs, 30 digits)",
       sp.N(est0.subs(tv), 30) - sp.N(((pm / pa) ** (1 - n) / (k * (1 - n))).subs(tv), 30))
val0 = sp.N(((pm / pa) ** (1 - n) / (k * (1 - n))).subs(tv), 8)
print(f"      TIMs, p_min = 0.505: eps*_f(0) = {val0}; a plastic contraction |deps^p_v| larger than this from the "
      f"floor leaves dom Psi (evaluation failure -> backtrack)")
print("\nALL OK" if ok_all else "\nSOME CHECK FAILED")
