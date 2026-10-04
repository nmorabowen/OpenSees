"""Fast checks of the round-3b harness wiring (seconds; run before the long smokes).

  1. HAR plug: with the plug's (n_e, g, k, p_a) the ORACLE's own elastic law (O2 kernel.elastic) has, on the axis, the
     DM04 target G(p, e) and K(p, e) at every p (n = 1/2 closed), for 'per_test' and 'global' policies; the axis Poisson
     ratio is nu; and the params carry none of the five BA06 values (refused by both oracles under HAR, sheet §2.4).
  2. BA06 plug: unchanged numbers (mu0 = G(p_rep, e), kappa_hat = p_rep/K(p_rep, e), p0 = -p_rep).
  3. The unified pi_i0 rule (S.53) computed by each ORACLE's helper equals the closed form pi0_for (the 'ramp_end'
     surface through (p_init, c2 M)) for isotropic starts, both energies, caps smooth / planar / none and N = 0 / > 0.
  4. A short PS driver run from the unified start completes under BA06 and HAR on O2 and O1 and the two agree
     (first-order level, 40 increments vs O1 truth).
  5. ElasticPolicy.bvp() ('global'): one set of energy parameters for every test of a BVP body.
  6. tcl_flags for the deck (global policy) are printed.
Run (norsand_calib/):  python -m harness.check_plugs
"""
from __future__ import annotations

import sys

import numpy as np

from . import energy as EN
from .drivers import simulate
from .model import Setup, initial, make_params, pi0_for, oracle_module
from .sand import TOYOURA_PLACEHOLDER as TOY, toyoura_dm04, TIMS_P_A

THETA = dict(chi=-3.0, h=150.0, N=0.3, N_bar=0.2, rho=0.712, rho_bar=0.75)


def ok(name, err, tol):
    good = bool(err <= tol)
    print(f"  [{'OK' if good else 'FAIL'}] {name}: {err:.3e} (tol {tol:.0e})")
    return good


def check_har():
    from o2_algo import kernel as K
    from o2_algo import api as A
    good = True
    sand = toyoura_dm04()
    assert sand.p_a == TIMS_P_A == 101.0
    for pol_name, pol in (("per_test", EN.ElasticPolicy()), ("global", EN.ElasticPolicy.bvp(e_rep=0.635))):
        for e_init in (0.635, 0.716):
            e_ref = pol.e_ref(e_init)
            setup = Setup(sand, energy="HAR", policy=pol, pi0_rule="unified")
            P = make_params("O2", setup, THETA, 49.0, e_init)
            assert P.energy == "HAR" and P.p_a == sand.p_a
            for nm in ("p0", "kappa_hat", "eps_v0", "mu0", "alpha0"):
                assert getattr(P, nm) is None, nm
            assert not P._ba06_given, P._ba06_given
            worst_G = worst_K = 0.0
            for p in (1.0, 4.9, 49.0, 101.0, 300.0):
                # the strain at isotropic p from the inverse map, then the oracle's tangent moduli on the axis
                eps = K.invert_elastic(P, np.array([-p, -p, -p]))
                el = K.elastic(P, eps)
                K_o, G_o = el.D11, el.D22 / 3.0
                worst_K = max(worst_K, abs(K_o / sand.elastic.K(p, e_ref) - 1.0))
                worst_G = max(worst_G, abs(G_o / sand.elastic.G(p, e_ref) - 1.0))
            good &= ok(f"HAR {pol_name} e_init {e_init}: oracle G(p) vs DM04 target, p 1..300", worst_G, 1e-12)
            good &= ok(f"HAR {pol_name} e_init {e_init}: oracle K(p) vs DM04 target, p 1..300", worst_K, 1e-12)
            d = EN.get("HAR").describe(sand, 49.0, e_init, pol)
            good &= ok(f"HAR {pol_name} e_init {e_init}: axis nu vs target", abs(d["nu_axis"] - d["nu_target"]), 1e-14)
    # O1 accepts the same HAR params
    setup = Setup(sand, energy="HAR")
    P1 = make_params("O1", setup, THETA, 49.0, 0.716)
    good &= ok("HAR on O1: Params constructed and energy == 'HAR'", 0.0 if P1.energy == "HAR" else 1.0, 0.0)
    # the BA06 values are refused if they sneak in
    try:
        oracle_module("O2").Params(energy="HAR", k=1.0, g=1.0, n_e=0.5, p0=-100.0).validate()
        good &= ok("HAR refuses an explicit p0", 1.0, 0.0)
    except ValueError:
        good &= ok("HAR refuses an explicit p0", 0.0, 0.0)
    return good


def check_ba06():
    good = True
    sand = toyoura_dm04()
    pol = EN.ElasticPolicy()
    d = EN.get("BA06").params(sand, 49.0, 0.716, pol)
    good &= ok("BA06 mu0 = G(p_init, e)", abs(d["mu0"] / sand.elastic.G(49.0, 0.716) - 1.0), 1e-15)
    good &= ok("BA06 kappa_hat = p/K", abs(d["kappa_hat"] / (49.0 / sand.elastic.K(49.0, 0.716)) - 1.0), 1e-15)
    good &= ok("BA06 p0 = -p_rep", abs(d["p0"] + 49.0), 0.0)
    return good


def check_unified():
    good = True
    sand = toyoura_dm04()
    worst = 0.0
    n = 0
    for oracle in ("O2", "O1"):
        for energy in ("BA06", "HAR"):
            for cap, c1, c2 in (("smooth", 0.05, 0.15), ("planar", 0.1, 0.1), ("none", 0.05, 0.15)):
                for N in (0.0, 0.3):
                    th = dict(THETA, N=N, N_bar=min(0.2, N))
                    th["rho_bar"] = 0.75 if N > 0 else th["rho"]          # N = 0: rho/rho_bar >= 1 (S.39)
                    for p in (4.9, 49.0, 400.0):
                        setup = Setup(sand, energy=energy, cap=cap, c1=c1, c2=c2, pi0_rule="unified")
                        _, st = initial(oracle, setup, th, p, 0.716)
                        ref = pi0_for(setup, th, p)
                        worst = max(worst, abs(st.pi_i / ref - 1.0))
                        n += 1
    good &= ok(f"unified pi_i0 (oracle helper) vs closed form, {n} cases (O1+O2, BA06+HAR, 3 caps, N 0 and 0.3)", worst, 1e-12)
    # the rule is the c2 M surface, not the apex, for the smooth cap
    setup = Setup(sand)
    _, st = initial("O2", setup, THETA, 49.0, 0.716)
    apex = Setup(sand, pi0_rule="on_surface")
    good &= ok("unified != apex for the smooth cap (pi_i0 ratio apex/unified differs from 1 by > 5 %)",
               0.0 if abs(pi0_for(apex, THETA, 49.0) / st.pi_i - 1.0) > 0.05 else 1.0, 0.0)
    return good


def check_short_ps():
    """O2 at n = 60 and 120 increments against the O1 truth (3 % axial, PS, unified start): both complete and the
    endpoint errors halve (first order); eps_v compared in absolute % points (it is near zero at 3 %)."""
    good = True
    sand = toyoura_dm04()
    for energy in ("BA06", "HAR"):
        setup = Setup(sand, energy=energy)
        t1 = simulate(setup, THETA, "PS", 49.0, 0.716, 0.03, 12, oracle="O1")
        cs = [simulate(setup, THETA, "PS", 49.0, 0.716, 0.03, n, oracle="O2") for n in (60, 120)]
        good &= ok(f"{energy} PS 3 %: O1 and both O2 runs complete", 0.0 if (t1.complete and all(c.complete for c in cs)) else 1.0, 0.0)
        if t1.complete and all(c.complete for c in cs):
            es = [np.linalg.norm(c.sig[-1] - t1.sig[-1]) / np.linalg.norm(t1.sig[-1]) for c in cs]
            ev = [abs(c.eps_v_pct[-1] - t1.eps_v_pct[-1]) for c in cs]
            good &= ok(f"{energy} PS 3 % endpoint sigma error (n = 60, 120: {es[0]:.2e}, {es[1]:.2e}) halves: ratio - 0.5",
                       abs(es[1] / es[0] - 0.5), 0.15)
            good &= ok(f"{energy} PS 3 % endpoint |eps_v error| in % points at n = 120 ({ev[0]:.2e} at n = 60)", ev[1], 0.05)
    return good


def check_global_policy():
    """ElasticPolicy.bvp(): one material for every test. BA06 'global': the same (p0, kappa_hat, mu0) at 4.9 and 49 kPa
    and at two initial void ratios (e_rep fixed); 'per_test' gives different ones. HAR 'global': the same (k, g, n) for
    both tests (and they do not depend on p_init at all, in either policy)."""
    good = True
    sand = toyoura_dm04()
    pol = EN.ElasticPolicy.bvp(p_rep=49.0, e_rep=0.635)
    for nm in ("BA06", "HAR"):
        a = EN.get(nm).params(sand, 4.9, 0.714, pol)
        b = EN.get(nm).params(sand, 49.0, 0.716, pol)
        good &= ok(f"{nm} 'global': identical parameters for (4.9 kPa, e 0.714) and (49 kPa, e 0.716)",
                   0.0 if a == b else 1.0, 0.0)
    a = EN.get("BA06").params(sand, 4.9, 0.714, EN.ElasticPolicy())
    b = EN.get("BA06").params(sand, 49.0, 0.716, EN.ElasticPolicy())
    good &= ok("BA06 'per_test': different parameters for the two tests", 0.0 if a != b else 1.0, 0.0)
    h1 = EN.get("HAR").params(sand, 4.9, 0.714, EN.ElasticPolicy())
    h2 = EN.get("HAR").params(sand, 49.0, 0.714, EN.ElasticPolicy())
    good &= ok("HAR 'per_test': independent of p_init (same e)", 0.0 if h1 == h2 else 1.0, 0.0)
    try:
        EN.ElasticPolicy.bvp(e_rep=0.635).rep(49.0, 0.7)
        good &= ok("BA06 'global' without p_rep is refused", 1.0, 0.0)
    except ValueError:
        good &= ok("BA06 'global' without p_rep is refused", 0.0, 0.0)
    return good


def show_flags():
    sand = toyoura_dm04()
    pol = EN.ElasticPolicy.bvp(p_rep=49.0, e_rep=0.635)
    for nm in ("BA06", "HAR"):
        print(f"  tcl_flags {nm} (global, e_rep 0.635, p_rep 49): {EN.get(nm).tcl_flags(sand, pol)}")


def main():
    print("HAR plug"); g = check_har()
    print("BA06 plug"); g &= check_ba06()
    print("unified pi_i0"); g &= check_unified()
    print("global (BVP) policy"); g &= check_global_policy()
    print("short PS drivers"); g &= check_short_ps()
    print("deck flag lists"); show_flags()
    print("CHECK_PLUGS", "PASS" if g else "FAIL")
    return g


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
