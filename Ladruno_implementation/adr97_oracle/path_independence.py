"""ADR-97 P0 oracle 5/6 -- the gate that separates the CLOSEST-POINT map from the
shipped CUTTING-PLANE (Ortiz-Simo) ``Backward_Euler``.

WHAT IS BEING SHOWN
-------------------
1. **Both maps are consistent**: halving the step size drives both toward the
   same limit, at O(h) (the ratio of successive step-halving errors tends to 2).
2. **The CPPM per-step state is UNIQUE.**  Re-running the same step from the
   same committed start with a different Newton start guess reproduces the
   committed ``(sigma, alpha, k)`` to <= 1e-9 -- the answer is defined by the
   residual, not by the iteration.
3. **The cutting-plane per-step state is NOT.**  The shipped loop
   (``ASDPlasticMaterial3D.h:2458-2510``) accepts on ``|Phi| < tol_yf`` alone and
   accumulates the internal variables as ``sum_k deltaLambda_k * h(sigma_k)``
   along its OWN iterates.  Two iterate sequences that both satisfy the same
   acceptance test therefore commit DIFFERENT internal variables whenever ``h``
   is not constant -- i.e. for Armstrong-Frederick.  With linear hardening
   (``h`` constant along the iterates) the two agree to round-off, which is why
   ADR-94 could not see this on von Mises + linear hardening.

THE CUTTING PLANE MIRRORED HERE, LINE BY LINE
---------------------------------------------
``ASDPlasticMaterial3D.h`` Backward_Euler, per iterate::

    n   = yf.df_dsigma_ij(TrialStress, iv)            (2459)
    m   = pf(depsilon, TrialStress, iv)               (2460)
    H   = yf.hardening(depsilon, m, TrialStress, iv)  (2461)
    Phi = yf(TrialStress, iv)                         (2462)
    if |Phi| < tol_yf: break                          (2463)
    dPhi_dLambda = H - n.dot(Eelastic*m)              (2471)
    deltaLambda  = -Phi / dPhi_dLambda                (2482)
    TrialStress         -= deltaLambda * (Eelastic*m) (2494)
    TrialPlastic_Strain += deltaLambda * m            (2495)
    iv.trial_value      += deltaLambda * h(depsilon, m, TrialStress)   (2496-2502)

Note that ``h`` is evaluated with the OLD ``m`` and the IV's OWN old value, and
the stress it sees is already the UPDATED one -- reproduced exactly below.

Run::

    python3.12 Ladruno_implementation/adr97_oracle/path_independence.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from asd_common import (SQRT_2_over_3, banner, section, dev, dot_stress)
from cppm_vm import VM
from cppm_dp import DP

C23 = SQRT_2_over_3
TOL_YF = 1e-10          # |Phi| acceptance, the C++ yf_tolerance() stand-in


# ---------------------------------------------------------------------------
# cutting-plane (Ortiz-Simo) maps -- mirrors of the shipped Backward_Euler
# ---------------------------------------------------------------------------
def cp_vm(mat, state, eps_new, relax=1.0, tol=TOL_YF, max_iter=200):
    eps_n, sig_n, alpha_n, k_n = state
    sig = sig_n + mat.Ee @ (eps_new - eps_n)
    alpha, k, dl = alpha_n.copy(), float(k_n), 0.0
    if mat.f(sig, alpha, k) <= 0.0:
        return (eps_new, sig, alpha, k), 0
    for it in range(max_iter):
        Phi = mat.f(sig, alpha, k)
        if abs(Phi) < tol:
            return (eps_new, sig, alpha, k), it
        n = m = mat.m_of(sig, alpha)              # associated
        hk = mat.h_k(m)
        halpha = mat.h_alpha(m, alpha)
        H = -C23 * hk + (-n) @ halpha             # VonMises_YF::hardening
        ddl = -Phi / (H - n @ (mat.Ee @ m))
        if it == 0:
            ddl *= relax                          # a DIFFERENT, still-valid path
        if dl + ddl < 0.0:
            return (eps_new, sig_n + mat.Ee @ (eps_new - eps_n),
                    alpha_n.copy(), float(k_n)), it
        dl += ddl
        sig = sig - ddl * (mat.Ee @ m)
        alpha = alpha + ddl * halpha
        k = k + ddl * hk
    return (eps_new, sig, alpha, k), max_iter


def cp_dp(mat, state, eps_new, relax=1.0, tol=TOL_YF, max_iter=200):
    eps_n, sig_n, k_n = state
    sig = sig_n + mat.Ee @ (eps_new - eps_n)
    k, dl = float(k_n), 0.0
    if mat.f(sig, k) <= 0.0:
        return (eps_new, sig, k), 0
    for it in range(max_iter):
        Phi = mat.f(sig, k)
        if abs(Phi) < tol:
            return (eps_new, sig, k), it
        d = dev(sig)
        q = np.sqrt(0.5 * dot_stress(d, d))
        n = (np.array([1., 1., 1., 2., 2., 2.]) * d / (2.0 * q)
             + (mat.eta / 3.0) * np.array([1., 1., 1., 0., 0., 0.]))
        m = mat.m_of(sig)
        hk = mat.hk
        H = -1.0 * hk                              # df/dk = -1
        ddl = -Phi / (H - n @ (mat.Ee @ m))
        if it == 0:
            ddl *= relax
        if dl + ddl < 0.0:
            return (eps_new, sig_n + mat.Ee @ (eps_new - eps_n), float(k_n)), it
        dl += ddl
        sig = sig - ddl * (mat.Ee @ m)
        k = k + ddl * hk
    return (eps_new, sig, k), max_iter


# ---------------------------------------------------------------------------
def legs_to_steps(legs, nstep, start=None):
    out, prev = [], np.zeros(6) if start is None else start
    for leg in legs:
        for i in range(nstep):
            out.append(prev + (leg - prev) * (i + 1) / nstep)
        prev = leg
    return out


def run_vm(mat, legs, nstep, mapper, **kw):
    state = (np.zeros(6), np.zeros(6), np.zeros(6), mat.k0)
    for e in legs_to_steps(legs, nstep):
        if mapper == "cppm":
            state, _ = mat.step(state, e, verbose=False, want_tangent=False)
        else:
            state, _ = cp_vm(mat, state, e, **kw)
    return state


def run_dp(mat, legs, nstep, mapper, **kw):
    state = (np.zeros(6), np.zeros(6), 0.0)
    for e in legs_to_steps(legs, nstep):
        if mapper == "cppm":
            state, _ = mat.step(state, e, verbose=False, want_tangent=False)
        else:
            state, _ = cp_dp(mat, state, e, **kw)
    return state


VM_LEGS = [np.array([0., 0., 0., 3.0e-3, 0., 0.]),
           np.array([0., 0., 0., 3.0e-3, 3.0e-3, 0.]),
           np.array([-6.0e-4, -6.0e-4, 2.0e-3, 1.0e-3, 3.0e-3, 0.])]
DP_LEGS = [np.array([3.0e-4, 3.0e-4, -2.5e-3, 0., 0., 0.]),
           np.array([3.0e-4, 3.0e-4, -2.5e-3, 2.0e-3, 0., 0.]),
           np.array([0., 0., -1.0e-3, 2.0e-3, 1.0e-3, 0.])]


def step_halving(label, runner, mat, legs, nlist, ref_n):
    print(f"\n    {label}")
    ref = runner(mat, legs, ref_n, None)[1] if False else None
    res = {n: runner(mat, legs, n) for n in nlist + [ref_n]}
    ref = res[ref_n]
    print(f"      steps/leg   |sigma(N) - sigma(ref N={ref_n})|_inf     ratio")
    prev = None
    errs = {}
    for n in nlist:
        e = float(np.max(np.abs(res[n][1] - ref[1])))
        errs[n] = e
        r = "" if prev is None else f"{prev / e:8.3f}"
        print(f"      {n:9d}   {e:.6e}                    {r}")
        prev = e
    return res, errs


def main():
    banner("ADR-97 oracle 5/6 -- path independence: CPPM vs cutting plane")

    # =====================================================================
    section("A. von Mises + ARMSTRONG-FREDERICK, 3-leg path, step halving")
    vm = VM(H=0.0, ha=15000.0, cr=300.0)
    print("    material: E=70000 nu=0.3 k0=30 ha=15000 cr=300 (header AF form)")
    print("    legs: shear_xy -> add shear_yz -> add a triaxial + partial unload")
    res_c, err_c = step_halving(
        "CPPM (closest point)", lambda m, l, n: run_vm(m, l, n, "cppm"),
        vm, VM_LEGS, [5, 10, 20, 40], 640)
    res_p, err_p = step_halving(
        "CUTTING PLANE (mirror of Backward_Euler)",
        lambda m, l, n: run_vm(m, l, n, "cp"), vm, VM_LEGS, [5, 10, 20, 40], 640)
    d = float(np.max(np.abs(res_c[640][1] - res_p[640][1])))
    print(f"\n      the two maps' N=640 limits agree to |dsigma|_inf = {d:.6e}"
          "   (both are consistent)")
    for n in (5, 10, 20, 40):
        print(f"      N={n:3d}: |sigma_CPPM - sigma_CP|_inf = "
              f"{float(np.max(np.abs(res_c[n][1] - res_p[n][1]))):.6e}"
              f"   |alpha_CPPM - alpha_CP|_inf = "
              f"{float(np.max(np.abs(res_c[n][2] - res_p[n][2]))):.6e}")
    assert err_c[5] / err_c[40] > 4.0 and err_p[5] / err_p[40] > 4.0

    # =====================================================================
    section("B. Drucker-Prager + linear hardening, 3-leg path, step halving")
    dp = DP(etabar=0.2, H=500.0)
    res_c2, err_c2 = step_halving(
        "CPPM (closest point)", lambda m, l, n: run_dp(m, l, n, "cppm"),
        dp, DP_LEGS, [5, 10, 20, 40], 640)
    res_p2, err_p2 = step_halving(
        "CUTTING PLANE (mirror of Backward_Euler)",
        lambda m, l, n: run_dp(m, l, n, "cp"), dp, DP_LEGS, [5, 10, 20, 40], 640)
    d2 = float(np.max(np.abs(res_c2[640][1] - res_p2[640][1])))
    print(f"\n      the two maps' N=640 limits agree to |dsigma|_inf = {d2:.6e}")
    for n in (5, 10, 20, 40):
        print(f"      N={n:3d}: |sigma_CPPM - sigma_CP|_inf = "
              f"{float(np.max(np.abs(res_c2[n][1] - res_p2[n][1]))):.6e}"
              f"   |k_CPPM - k_CP| = "
              f"{abs(res_c2[n][2] - res_p2[n][2]):.6e}")

    # =====================================================================
    section("C. SAME step, SAME start, DIFFERENT iterate path")
    print("    Start state: 6 CPPM steps along leg 1+2, then ONE more step,")
    print("    solved several ways.  Every variant below satisfies the SAME")
    print("    acceptance test the C++ uses (|Phi| < tol_yf).")

    # ---- C1: von Mises + AF -------------------------------------------
    print("\n    C1  von Mises + Armstrong-Frederick (h NOT constant)")
    st = (np.zeros(6), np.zeros(6), np.zeros(6), vm.k0)
    warm = legs_to_steps(VM_LEGS[:2], 3)
    for e in warm:
        st, _ = vm.step(st, e, verbose=False, want_tangent=False)
    eps_next = warm[-1] + np.array([0., 0., 0., 0., 1.2e-3, 0.])

    base, _ = vm.step(st, eps_next, verbose=False, want_tangent=False)
    print(f"      CPPM, natural start guess:")
    print(f"        sigma = {np.array2string(base[1], precision=12)}")
    print(f"        alpha = {np.array2string(base[2], precision=12)}")
    # perturbed Newton start guesses
    from asd_common import newton
    sig_tr = st[1] + vm.Ee @ (eps_next - st[0])
    R = vm.residual(sig_tr, st[2], st[3])
    scale = float(np.max(np.abs(sig_tr)))
    worst = 0.0
    for tag, x0 in (("+30% stress, dl = 2e-4",
                     np.concatenate([sig_tr * 1.3, st[2] * 0.5, [st[3] * 1.1],
                                     [2e-4]])),
                    ("-15% stress, dl = 1e-4",
                     np.concatenate([sig_tr * 0.85, st[2] * 1.3, [st[3] * 0.95],
                                     [1e-4]]))):
        x, _ = newton(R, x0, f"CPPM alt start ({tag})", scale=scale,
                      verbose=False, iter_gate=8)
        e = max(float(np.max(np.abs(x[0:6] - base[1]))),
                float(np.max(np.abs(x[6:12] - base[2]))),
                abs(x[12] - base[3]))
        worst = max(worst, e)
        print(f"      CPPM, start = {tag:24s} -> max|state diff| = {e:.3e}")
    print(f"      => CPPM per-step state invariance: {worst:.3e}  (gate <= 1e-9)")
    assert worst <= 1e-9

    cp_states = {}
    for relax in (1.0, 0.5, 0.25, 1.5):
        s, its = cp_vm(vm, st, eps_next, relax=relax)
        cp_states[relax] = s
        print(f"      CUTTING PLANE, first corrector x {relax:<5}"
              f" ({its:2d} iters, |Phi| = {abs(vm.f(s[1], s[2], s[3])):.3e})")
        print(f"        sigma = {np.array2string(s[1], precision=12)}")
        print(f"        alpha = {np.array2string(s[2], precision=12)}")
    spread_sig = max(float(np.max(np.abs(cp_states[r][1] - cp_states[1.0][1])))
                     for r in cp_states)
    spread_alp = max(float(np.max(np.abs(cp_states[r][2] - cp_states[1.0][2])))
                     for r in cp_states)
    print(f"      => CUTTING PLANE spread over valid iterate paths: "
          f"sigma {spread_sig:.6e}, alpha {spread_alp:.6e}")
    print(f"         (all four commit |Phi| < {TOL_YF:g}: the acceptance test "
          f"cannot tell them apart)")
    assert spread_alp > 1e-6, "AF cutting plane should be path dependent"

    # ---- C2: von Mises + LINEAR hardening (the ADR-94 blind spot) -------
    print("\n    C2  von Mises + LINEAR isotropic hardening (h IS constant)")
    vml = VM(H=7000.0)
    stl = (np.zeros(6), np.zeros(6), np.zeros(6), vml.k0)
    for e in warm:
        stl, _ = vml.step(stl, e, verbose=False, want_tangent=False)
    sl = {r: cp_vm(vml, stl, eps_next, relax=r)[0] for r in (1.0, 0.5, 0.25)}
    sp = max(float(np.max(np.abs(sl[r][1] - sl[1.0][1]))) for r in sl)
    spk = max(abs(sl[r][3] - sl[1.0][3]) for r in sl)
    print(f"      CUTTING PLANE spread: sigma {sp:.3e}, k {spk:.3e}")
    print("      => with h constant the cutting plane IS path independent, which"
          " is\n         exactly why ADR-94 H6 could not see the defect on "
          "von Mises.")
    assert sp < 1e-9 and spk < 1e-9

    # ---- C3: Drucker-Prager + linear hardening --------------------------
    print("\n    C3  Drucker-Prager + linear hardening (h constant here too)")
    std = (np.zeros(6), np.zeros(6), 0.0)
    warmd = legs_to_steps(DP_LEGS[:2], 3)
    for e in warmd:
        std, _ = dp.step(std, e, verbose=False, want_tangent=False)
    eps_d = warmd[-1] + np.array([0., 0., -3e-4, 5e-4, 0., 0.])
    sd = {r: cp_dp(dp, std, eps_d, relax=r)[0] for r in (1.0, 0.5, 0.25)}
    spd = max(float(np.max(np.abs(sd[r][1] - sd[1.0][1]))) for r in sd)
    bd, _ = dp.step(std, eps_d, verbose=False, want_tangent=False)
    print(f"      CUTTING PLANE spread: sigma {spd:.3e}")
    print(f"      CPPM vs cutting plane on this single step: "
          f"|dsigma|_inf = {float(np.max(np.abs(bd[1] - sd[1.0][1]))):.6e}, "
          f"|dk| = {abs(bd[2] - sd[1.0][2]):.6e}")

    banner("PATH-INDEPENDENCE REFERENCE BLOCK (values C++ tests pin)")
    print(f"  A. VM+AF step halving, |sigma(N) - sigma(640)|_inf:")
    print(f"     CPPM         : " + "  ".join(f"N={n}:{err_c[n]:.4e}"
                                              for n in (5, 10, 20, 40)))
    print(f"     cutting plane: " + "  ".join(f"N={n}:{err_p[n]:.4e}"
                                              for n in (5, 10, 20, 40)))
    print(f"     limits agree to {d:.4e}")
    print(f"  B. DP+linear step halving, |sigma(N) - sigma(640)|_inf:")
    print(f"     CPPM         : " + "  ".join(f"N={n}:{err_c2[n]:.4e}"
                                              for n in (5, 10, 20, 40)))
    print(f"     cutting plane: " + "  ".join(f"N={n}:{err_p2[n]:.4e}"
                                              for n in (5, 10, 20, 40)))
    print(f"     limits agree to {d2:.4e}")
    print(f"  C1. VM+AF single step: CPPM invariance {worst:.3e} (<= 1e-9);")
    print(f"      cutting-plane spread sigma {spread_sig:.6e}, "
          f"alpha {spread_alp:.6e}")
    print(f"  C2. VM+linear: cutting-plane spread sigma {sp:.3e} "
          f"(path INdependent -- the ADR-94 blind spot)")
    print(f"  C3. DP+linear: cutting-plane spread sigma {spd:.3e}")


if __name__ == "__main__":
    main()
