"""ADR 90 leg A0 -- regenerate every table in `Ladruno_implementation/_adr90_a0_results.md`.

Run:  python3.12 tests/_testbed/run_a0_sweep.py            (full sweep, ~2 min)
      python3.12 tests/_testbed/run_a0_sweep.py --quick    (the Zone-A subset, seconds)
      python3.12 tests/_testbed/run_a0_sweep.py --p0b      (only the P0b legs, ~1 min)

Numpy only -- no OpenSees, no build.  Everything printed here names the script, the git hash and
the date, per the fork's provenance rule.
"""
from __future__ import annotations

import datetime
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import duvaut_lions_ref as ref                                          # noqa: E402

NSTEPS = 2000                       # uniform step count for every viscous leg
U_MAX = 2.0                         # mm -- enough that even N=20 fully softens at tau=0
DE_LADDER = (1.0e-4, 3.0e-4, 5.0e-4, 1.0e-3, 3.0e-3, 1.0e-2)
DE_BRIEF = (0.01, 0.1, 1.0)         # the brief's requested set
NS = (20, 40, 80, 160)
CONV = ("one_element", "fixed_graded", "fixed_flat")


def stamp():
    try:
        h = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True,
                           text=True, check=True).stdout.strip()
    except Exception:                                                   # pragma: no cover
        h = "unknown"
    print("=" * 100)
    print("ADR 90 leg A0 -- Duvaut-Lions 1-D softening bar")
    print(f"script  tests/_testbed/run_a0_sweep.py   oracle  tests/_testbed/duvaut_lions_ref.py")
    print(f"git     {h}     date  {datetime.date.today().isoformat()}     numpy {np.__version__}")
    print(f"model   L=100 mm  A=1 mm^2  E=20000 MPa  sigY=20 MPa  H=-50 MPa (|H|/E = 1/400)")
    print(f"        residual 0.02*sigY with H_res=E/2e5; imperfection 10 % strength reduction")
    print(f"        end-displacement ramp u_max={U_MAX} mm over T with uniform dt=T/nsteps")
    print(f"        snap-back needs l_b < L|H|/E = 0.25 mm; h_min = 0.625 mm at N=160 => none")
    print("=" * 100)
    return h


def _fmt(x, w=9, p=4):
    return ("{:" + str(w) + "." + str(p) + "f}").format(x) if np.isfinite(x) else " " * (w - 3) + "n/a"


def table_h1(quick=False):
    print("\n### H1 (negative control) -- tau = 0: width proportional to h, energy -> 0 with h")
    print("w1 threshold count, w2 threshold-free second moment, w3 FWHM; W_50 = dissipated work")
    print("at the instant the load first falls to 50 % of peak (a mesh-FAIR energy measure).")
    print(f"{'convention':<14}{'N':>5}{'h':>9}{'nsteps':>8}{'w1':>10}{'w2':>10}{'w3':>10}"
          f"{'w2/h':>8}{'W_50':>10}{'W50 ratio':>11}{'P_peak':>9}")
    out = {}
    for imp in (CONV if not quick else ("one_element",)):
        prev = None
        for N in NS:
            n = NSTEPS
            r = ref.a0_bar(variant="tt", N=N, tau=0.0, T=1.0, nsteps=n, u_max=U_MAX,
                           imperfection=imp)
            while not r["converged"] and n < 64000:
                n *= 2
                r = ref.a0_bar(variant="tt", N=N, tau=0.0, T=1.0, nsteps=n, u_max=U_MAX,
                               imperfection=imp)
            ratio = (prev / r["W_50"]) if (prev and np.isfinite(r["W_50"])) else float("nan")
            print(f"{imp:<14}{N:>5}{r['h']:>9.4f}{n:>8}{_fmt(r['w1'],10)}{_fmt(r['w2'],10)}"
                  f"{_fmt(r['w3'],10)}{_fmt(r['w2']/r['h'],8,2)}{_fmt(r['W_50'],10)}"
                  f"{_fmt(ratio,11,3)}{_fmt(r['P_peak'],9)}")
            prev = r["W_50"] if np.isfinite(r["W_50"]) else None
            out[(imp, N)] = r
    return out


def table_h2(quick=False):
    print("\n### H2 -- at FIXED De, does width / curve converge under h-refinement?")
    print(f"{'convention':<14}{'variant':<8}{'De':>8}{'N':>5}{'w2':>10}{'w2 ratio':>10}"
          f"{'w3':>9}{'P_peak':>9}{'P_end':>9}{'W_end':>10}{'curveL2 vs N/2':>16}{'eps_p max':>11}")
    des = DE_LADDER if not quick else (3.0e-4, 1.0e-3)
    convs = CONV if not quick else ("fixed_graded",)
    vs = ("tt", "tdl") if not quick else ("tt",)
    out = {}
    for imp in convs:
        for var in vs:
            for De in des:
                prev_w, prev_r = None, None
                for N in NS:
                    r = ref.a0_bar(variant=var, N=N, tau=De, T=1.0, nsteps=NSTEPS,
                                   u_max=U_MAX, imperfection=imp)
                    out[(imp, var, De, N)] = r
                    if not r["converged"]:
                        print(f"{imp:<14}{var:<8}{De:>8.1e}{N:>5}   NOT CONVERGED")
                        continue
                    rat = (r["w2"] / prev_w) if prev_w else float("nan")
                    cd = ref.curve_distance(r, prev_r) if prev_r else float("nan")
                    print(f"{imp:<14}{var:<8}{De:>8.1e}{N:>5}{_fmt(r['w2'],10)}{_fmt(rat,10,3)}"
                          f"{_fmt(r['w3'],9)}{_fmt(r['P_peak'],9)}{_fmt(r['P_end'],9)}"
                          f"{_fmt(r['W_end'],10)}{_fmt(cd,16,3 if np.isfinite(cd) else 3)}"
                          f"{_fmt(r['epsp_max'],11)}")
                    prev_w, prev_r = r["w2"], r
            print()
    return out


def table_h3():
    print("\n### H3 -- De collapse: matched (tau, T) pairs at EQUAL De")
    print("(a) SAME step count.  beta = dt/(tau+dt) = 1/(1 + nsteps*De) and the strain increment")
    print("    u_max/nsteps carry NO other dependence on tau or T, so equal (De, nsteps) runs are")
    print("    BIT-IDENTICAL by construction.  Reported to full precision -- this gate cannot fail")
    print("    and must not be cited as evidence.")
    print(f"{'tau':>10}{'T':>10}{'De':>11}{'beta':>20}{'w2':>18}{'P_peak':>18}{'W_end':>18}")
    for tau, T in ((1.0e-3, 1.0), (1.0e-2, 10.0), (1.0e-4, 0.1)):
        r = ref.a0_bar(variant="tt", N=80, tau=tau, T=T, nsteps=NSTEPS, u_max=U_MAX,
                       imperfection="fixed_graded")
        print(f"{tau:>10.1e}{T:>10.3g}{r['De']:>11.1e}{r['beta']:>20.10e}"
              f"{r['w2']:>18.12f}{r['P_peak']:>18.12f}{r['W_end']:>18.12f}")
    print("\n(b) The SUBSTANTIVE version: matched De, DIFFERENT step counts (so beta differs).")
    print(f"{'tau':>10}{'T':>10}{'nsteps':>8}{'beta':>12}{'w2':>10}{'P_peak':>10}{'W_end':>10}")
    for tau, T, n in ((1.0e-3, 1.0, 1000), (1.0e-2, 10.0, 2000), (1.0e-4, 0.1, 4000)):
        r = ref.a0_bar(variant="tt", N=80, tau=tau, T=T, nsteps=n, u_max=U_MAX,
                       imperfection="fixed_graded")
        print(f"{tau:>10.1e}{T:>10.3g}{n:>8}{r['beta']:>12.4e}{_fmt(r['w2'],10)}"
              f"{_fmt(r['P_peak'],10)}{_fmt(r['W_end'],10)}")


def table_h4():
    print("\n### H4 -- TT (generic two-track wrapper) vs TDL (true Duvaut-Lions) ON THE BAR")
    print(f"{'convention':<14}{'De':>9}{'N':>5}{'w2 TT':>10}{'w2 TDL':>10}{'d w2':>10}"
          f"{'d P_peak':>10}{'d W':>10}{'curve L2':>10}")
    for imp in ("one_element", "fixed_graded"):
        for De in (0.0, 3.0e-4, 1.0e-3, 3.0e-3):
            for N in (40, 80, 160):
                a = ref.a0_bar(variant="tt", N=N, tau=De, T=1.0, nsteps=NSTEPS, u_max=U_MAX,
                               imperfection=imp)
                b = ref.a0_bar(variant="tdl", N=N, tau=De, T=1.0, nsteps=NSTEPS, u_max=U_MAX,
                               imperfection=imp)
                if not (a["converged"] and b["converged"]):
                    print(f"{imp:<14}{De:>9.1e}{N:>5}   NOT CONVERGED "
                          f"(tt={a['converged']}, tdl={b['converged']})")
                    continue
                dw = abs(a["w2"] - b["w2"]) / max(b["w2"], 1e-30)
                dp = abs(a["P_peak"] - b["P_peak"]) / b["P_peak"]
                dW = abs(a["W_end"] - b["W_end"]) / max(abs(b["W_end"]), 1e-30)
                print(f"{imp:<14}{De:>9.1e}{N:>5}{_fmt(a['w2'],10)}{_fmt(b['w2'],10)}"
                      f"{dw:>10.2e}{dp:>10.2e}{dW:>10.2e}{ref.curve_distance(a, b):>10.2e}")
            print()


def table_dt():
    print("\n### dt sensitivity at FIXED De (the transient is NOT dt-independent, even though")
    print("### PV3's steady overstress is).  N = 80, De = 1e-3, fixed_graded.")
    print(f"{'nsteps':>8}{'dt':>12}{'beta':>12}{'w2':>10}{'w3':>9}{'P_peak':>10}{'P_end':>10}{'W_end':>10}")
    prev = None
    for n in (250, 500, 1000, 2000, 4000, 8000):
        r = ref.a0_bar(variant="tt", N=80, tau=1.0e-3, T=1.0, nsteps=n, u_max=U_MAX,
                       imperfection="fixed_graded")
        print(f"{n:>8}{r['dt']:>12.3e}{r['beta']:>12.4e}{_fmt(r['w2'],10)}{_fmt(r['w3'],9)}"
              f"{_fmt(r['P_peak'],10)}{_fmt(r['P_end'],10)}{_fmt(r['W_end'],10)}")
        prev = r


def table_brief_de():
    print("\n### The brief's requested De in {0.01, 0.1, 1}: a NON-TAUTOLOGY check")
    print("`inelastic` = dissipated work / the elastic energy at peak load.  A run whose response")
    print("is essentially elastic converges trivially and proves nothing.")
    print(f"{'De':>7}{'N':>5}{'P_peak':>10}{'eps_p max':>11}{'inelastic':>11}{'w2':>10}{'w3':>9}")
    for De in DE_BRIEF:
        for N in (20, 160):
            r = ref.a0_bar(variant="tt", N=N, tau=De, T=1.0, nsteps=NSTEPS, u_max=U_MAX,
                           imperfection="fixed_graded")
            print(f"{De:>7.2f}{N:>5}{_fmt(r['P_peak'],10)}{_fmt(r['epsp_max'],11,5)}"
                  f"{_fmt(r['inelastic_ratio'],11,3)}{_fmt(r['w2'],10)}{_fmt(r['w3'],9)}")


def table_steps():
    print("\n### Steps-to-converge -- the cost of NOT regularizing (fixed_graded, TT)")
    print(f"{'N':>5}{'tau=0':>10}{'De=3e-4':>10}{'De=1e-3':>10}{'De=3e-3':>10}")
    for N in NS:
        row = []
        for De in (0.0, 3.0e-4, 1.0e-3, 3.0e-3):
            n = ref.min_steps_to_converge(variant="tt", N=N, tau=De, T=1.0, u_max=U_MAX,
                                          imperfection="fixed_graded")
            row.append("fail" if n is None else str(n))
        print(f"{N:>5}" + "".join(f"{v:>10}" for v in row))


def table_point():
    print("\n### Point-model TT vs TDL (the D2 evidence)")
    r = ref.run_tt_vs_tdl_point(verbose=False)
    print(f"  perfect plasticity, max |sig_TT - sig_TDL|          = {r['perfect_plastic_max_diff']:.2e}")
    print(f"  1-D LINEAR H/E in [-0.10, +0.10], monotonic         = {r['linear_max_d_path_rel']:.2e} (relative)")
    print(f"  1-D EXPONENTIAL law (alpha_end/alpha_f ~ 4.75)      = {r['exp_max_d_path_rel']:.2e} (relative)")
    print(f"  J2, PROPORTIONAL path                               = {r['j2_prop_max_d_rel']:.2e} (relative)")
    print(f"  J2, NON-PROPORTIONAL path (axial then shear)        = {r['j2_nonprop_max_d_rel']:.2e} (relative)")
    print(f"  1-D linear, load/UNLOAD/reload                      = {r['cyclic_max_d_path_rel']:.2e} (relative)")
    print(f"\n  {'H/E':>6}{'De':>8}{'sig_inv':>10}{'sig_TT':>10}{'sig_TDL':>10}{'d_path rel':>12}")
    for x in r["rows"]:
        print(f"  {x['H_over_E']:>6.2f}{x['De']:>8.3f}{x['sig_inv']:>10.4f}{x['sig_tt']:>10.4f}"
              f"{x['sig_tdl']:>10.4f}{x['d_path_rel']:>12.2e}")
    print(f"\n  {'H/E':>6}{'path':>10}{'De':>7}{'|s_inv|':>10}{'|s_TT|':>10}{'|s_TDL|':>10}{'rel':>11}")
    for x in r["j2_rows"]:
        print(f"  {x['H_over_E']:>6.2f}{x['path']:>10}{x['De']:>7.2f}{x['sig_inv_norm']:>10.3f}"
              f"{x['sig_tt_norm']:>10.3f}{x['sig_tdl_norm']:>10.3f}{x['d_rel']:>11.2e}")
    print(f"\n  {'H/E':>6}{'De':>8}{'sig_TT end':>12}{'sig_TDL end':>13}{'max rel diff':>14}  (cyclic)")
    for x in r["cyclic_rows"]:
        print(f"  {x['H_over_E']:>6.2f}{x['De']:>8.3f}{x['sig_tt_end']:>12.4f}"
              f"{x['sig_tdl_end']:>13.4f}{x['d_path_rel']:>14.2e}")
    return r


def table_p0b():
    """P0b -- the four measurements the 3-lens adversarial pass showed were missing."""
    print("\n\n########## P0b (a) -- STATE-DEPENDENT elastic operator ##########")
    a = {}
    for be in (0.3, 0.6, 1.2):
        print(f"\n--- E(sigma) = E0 (1 + {be} |sigma|/sigY) ---")
        a[be] = ref.run_pressure_dependent_leg(E_beta=be, verbose=True)
    print("\n  E_beta sensitivity (max relative TT-vs-TDL over the whole De ladder / over the "
          "working window 1e-4..1e-2):")
    print(f"  {'E_beta':>8}{'all De':>12}{'window':>12}{'window, non-buildable repair':>32}")
    for be, r in a.items():
        print(f"  {be:>8.1f}{r['var_max']:>12.2e}{r['window_max']:>12.2e}"
              f"{r['window_max_ttvp']:>32.2e}")

    print("\n\n########## P0b (b) -- DISSIPATION gate ##########")
    b = ref.run_dissipation_gate(verbose=True)

    print("\n\n########## P0b (c) -- is the width TAU-set or IMPERFECTION-set? ##########")
    c = ref.run_imperfection_study(verbose=True)

    print("\n\n########## P0b (d) -- non-associated DP + the BLENDED ACOUSTIC TENSOR ##########")
    d = ref.run_dp_leg(verbose=True)

    print("\n\n########## P0b VERDICTS ##########")
    print(f"  V1 generic TT over a pressure-dependent inner: max rel error "
          f"{a[0.6]['window_max']:.1e} over the working window De 1e-4..1e-2 "
          f"(E_beta=0.6), rising to {a[0.6]['var_max']:.1e} at De=0.3; the identity is EXACT "
          f"({a[0.6]['const_max']:.1e}) only for a constant elastic operator.")
    print(f"  V2 dissipation: significant negatives on {b['paths_with_significant_negatives']}; "
          f"worst step {b['worst_step']:.2e} ({b['worst_step_rel']:.1e} of that run's "
          f"cumulative); worst cumulative-negative fraction {b['neg_sum_rel_worst']:.2e}.")
    print(f"  V3 width: at FIXED De=3e-4 the converged w2 moves x{c['w2_span_over_amp']:.2f} with "
          f"imperfection AMPLITUDE and x{c['w2_span_over_len']:.2f} with ZONE LENGTH.")
    print(f"  V4 blended acoustic tensor: inviscid min det = {d['min_det_inviscid']:.4f} "
          f"(LOST); elliptic for beta < {d['beta_crit']:.4f}, i.e. "
          f"De > {(1.0/d['beta_crit']-1.0):.3e}/nsteps.")
    return a, b, c, d


def main():
    quick = "--quick" in sys.argv
    stamp()
    if "--p0b" in sys.argv:
        table_p0b()
        return
    print("\n\n########## PV1-PV6 falsification battery ##########")
    ref.run_pv_gate(verbose=True)
    print("\n\n########## Point models: TT vs TDL ##########")
    table_point()
    print("\n\n########## A0 bar ##########")
    table_h1(quick)
    table_h2(quick)
    if not quick:
        table_h3()
        table_h4()
        table_dt()
        table_brief_de()
        table_steps()
        table_p0b()


if __name__ == "__main__":
    main()
