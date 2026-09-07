"""ADR-97 P0 oracle 2/5 -- Drucker-Prager CONE + APEX closest-point return maps
with consistent tangents (associated and NON-associated).

MIRRORS
-------
``SRC/material/nD/ASDPlasticMaterial3D/YieldFunctions/DruckerPrager_YF.h``::

    f = sqrt(J2(r)) + eta * p - xi_c ,      r = dev(sigma) - alpha
    p = sigma.meanStress() = trace/3          (TENSION POSITIVE)
    sqrt(J2) = sqrt(0.5 * tensor_dot_stress_like(r, r))

and its VOIGT gradient (header lines 63-91), which is NOT flat across slots::

    n = W * r / (2 sqrt(J2)) + (eta/3) * [1,1,1,0,0,0],  W = [1,1,1,2,2,2]

``PlasticFlowDirections/DruckerPrager_PF.h`` is the same expression with
``DP_etabar`` in place of ``DP_eta``: NON-ASSOCIATED whenever
``etabar != eta``, and the consistent tangent is then UNSYMMETRIC.

Apex (header lines 138-186): ``p_apex = xi_c / eta``, ``sigma_apex = p_apex *
delta``, and the header's region test is the EUCLIDEAN one
``p - p_apex >= eta * q``, with an explicit caveat in the header that the exact
condition in the elastic metric is ``p - p_apex >= (K*etabar/G) * q``.  This
oracle derives that exact test (see ``region_of``) and PRINTS BOTH, because the
two disagree on a measurable band of trial states.

CONVENTION SURPRISE CARRIED FROM THE HEADER
-------------------------------------------
The shipped ``DruckerPrager_YF::yf`` uses the *parameter* ``DP_xi_c`` and has
its cohesion internal variable COMMENTED OUT (line 26,
``// auto eta = GET_TRIAL_INTERNAL_VARIABLE(CohesionHardeningType);``), yet
``yf_hardening`` still contributes ``df/dk = -1`` times that IV's rate.  So a
hardening Drucker-Prager in the shipped code has a hardening term in ``H``
that NO term of ``f`` matches.  This oracle pins the SELF-CONSISTENT model

    f = sqrt(J2(r)) + eta * p - (xi_c + k)

(``df/dk = -1``, exactly what ``yf_hardening`` assumes) and flags that P1 must
either restore ``k`` in ``f`` or delete the hardening term.  With ``H = 0`` the
two coincide, so the perfect-plasticity block below pins the shipped code
exactly as it stands.

Run::

    python3.12 Ladruno_implementation/adr97_oracle/cppm_dp.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from asd_common import (W_STRESS, P_DEV, banner, section, elastic_voigt, dev,
                        dot_stress, dot_strain, mean_stress, newton,
                        jac_complex_step, check_complex_step_jacobian,
                        tangent_from_jacobian, tangent_general,
                        assert_fd_tangent, rel_fro, print_matrix, print_vec)

DELTA = np.array([1., 1., 1., 0., 0., 0.])


class DP:
    """Drucker-Prager, scalar linear cohesion hardening, no back stress."""

    def __init__(self, E=30000.0, nu=0.25, xi_c=20.0, eta=0.4, etabar=0.4,
                 H=0.0):
        self.E, self.nu, self.xi_c, self.eta, self.etabar, self.H = \
            E, nu, xi_c, eta, etabar, H
        self.Ee = elastic_voigt(E, nu)
        self.Ci = np.linalg.inv(self.Ee)
        self.G = E / (2.0 * (1.0 + nu))
        self.K = E / (3.0 * (1.0 - 2.0 * nu))
        # dot_strain(m, m) is CONSTANT for this flow rule (see module doc of the
        # derivation): 0.5 from the deviatoric part + etabar^2/3 volumetric.
        self.m_dot_m = 0.5 + self.etabar ** 2 / 3.0
        self.hk = H * np.sqrt((2.0 / 3.0) * self.m_dot_m)
        # apex flow direction: pf() zeroes the deviatoric part when sqrt(J2)=0
        self.m_apex = (self.etabar / 3.0) * DELTA
        self.hk_apex = H * np.sqrt((2.0 / 3.0) * dot_strain(self.m_apex,
                                                            self.m_apex))

    # -- YF / PF -----------------------------------------------------------
    def q_of(self, sigma):
        d = dev(sigma)
        return np.sqrt(0.5 * dot_stress(d, d))

    def f(self, sigma, k):
        return self.q_of(sigma) + self.eta * mean_stress(sigma) - (self.xi_c + k)

    def m_of(self, sigma):
        d = dev(sigma)
        q = np.sqrt(0.5 * dot_stress(d, d))
        return W_STRESS * d / (2.0 * q) + (self.etabar / 3.0) * DELTA

    def p_apex(self, k):
        return (self.xi_c + k) / self.eta

    # -- region test -------------------------------------------------------
    def region_of(self, sig_tr, k_n, verbose=True):
        """Which region the ELASTIC PREDICTOR belongs to.

        EXACT test (this oracle): the cone return is admissible iff it leaves
        ``sqrt(J2) >= 0``.  Closed form for the cone (derived in ``cone_closed``)
        gives ``q_{n+1} = q_tr - G*dLambda`` with
        ``dLambda = f_tr / (G + eta*K*etabar + hk)``, so

            APEX  <=>  q_tr - G*f_tr/(G + eta*K*etabar + hk) < 0
                  <=>  p_tr - p_apex > ((K*etabar + hk/eta)/G) * q_tr   [eta>0]

        HEADER test (``DruckerPrager_YF::check_apex_region``, line 171):
        ``p_tr - p_apex >= eta * q_tr``.  The two coincide only when
        ``K*etabar/G == eta`` (and ``H == 0``).
        """
        q = self.q_of(sig_tr)
        p = mean_stress(sig_tr)
        f_tr = q + self.eta * p - (self.xi_c + k_n)
        den = self.G + self.eta * self.K * self.etabar + self.hk
        exact_apex = (q - self.G * f_tr / den) < 0.0
        header_apex = (p - self.p_apex(k_n)) >= self.eta * q
        if verbose:
            print(f"    region test: p_tr = {p:.6f}, q_tr = {q:.6f}, "
                  f"p_apex = {self.p_apex(k_n):.6f}")
            print(f"      EXACT (elastic metric, slope "
                  f"{(self.K * self.etabar + self.hk / self.eta) / self.G:.6f}): "
                  f"{'APEX' if exact_apex else 'CONE'}")
            print(f"      HEADER (Euclidean, slope {self.eta:.6f}): "
                  f"{'APEX' if header_apex else 'CONE'}"
                  + ("   <-- DISAGREES with the exact test"
                     if header_apex != exact_apex else ""))
        return ("apex" if exact_apex else "cone"), header_apex, exact_apex

    # -- cone return -------------------------------------------------------
    def cone_residual(self, sig_tr, k_n):
        Ee = self.Ee

        def R(x):
            sigma, k, dl = x[0:6], x[6], x[7]
            m = self.m_of(sigma)
            out = np.empty(8, dtype=x.dtype)
            out[0:6] = sigma - sig_tr + dl * (Ee @ m)
            out[6] = k - k_n - dl * self.hk
            out[7] = self.f(sigma, k)
            return out
        return R

    def cone_closed(self, sig_tr, k_n):
        """Closed-form cone return + closed-form consistent tangent.

        ``E @ m`` splits exactly: deviatoric ``(G/q) r``, volumetric
        ``K*etabar*delta``.  Hence ``q_{n+1} = q_tr - G dl``,
        ``p_{n+1} = p_tr - K etabar dl`` and the consistency condition is LINEAR.
        """
        r_tr = dev(sig_tr)
        q_tr = np.sqrt(0.5 * dot_stress(r_tr, r_tr))
        p_tr = mean_stress(sig_tr)
        f_tr = q_tr + self.eta * p_tr - (self.xi_c + k_n)
        den = self.G + self.eta * self.K * self.etabar + self.hk
        dl = f_tr / den
        s = 1.0 - self.G * dl / q_tr
        sigma = (p_tr - self.K * self.etabar * dl) * DELTA + s * r_tr
        k = k_n + dl * self.hk
        PE = P_DEV @ self.Ee
        dq_deps = (W_STRESS * r_tr / (2.0 * q_tr)) @ PE
        dp_deps = self.K * DELTA
        ddl = (dq_deps + self.eta * dp_deps) / den
        ds = -self.G * (ddl * q_tr - dl * dq_deps) / q_tr ** 2
        C = (np.outer(DELTA, dp_deps - self.K * self.etabar * ddl)
             + s * PE + np.outer(r_tr, ds))
        return sigma, k, dl, C

    # -- apex return -------------------------------------------------------
    def apex_residual(self, sig_tr, k_n):
        def R(x):
            sigma, k, dl = x[0:6], x[6], x[7]
            dep = self.Ci @ (sig_tr - sigma)
            out = np.empty(8, dtype=x.dtype)
            out[0:6] = sigma - self.p_apex(k) * DELTA
            out[6] = k - k_n - dl * self.hk_apex
            out[7] = dl * self.etabar - (dep[0] + dep[1] + dep[2])
            return out
        return R

    def apex_dR_deps(self):
        dR = np.zeros((8, 6))
        dR[7, :] = -(DELTA @ (self.Ci @ self.Ee))       # == -delta^T
        return dR

    def apex_closed(self, sig_tr, k_n):
        """Closed form: the apex map is affine, and its tangent is rank <= 1
        and purely VOLUMETRIC (the 'bulk-only projector'); exactly ZERO for
        perfect plasticity because the apex point does not move."""
        p_tr = mean_stress(sig_tr)
        a = self.hk_apex / (self.K * self.etabar)
        k = (k_n + a * (p_tr - self.xi_c / self.eta)) / (1.0 + a / self.eta)
        sigma = self.p_apex(k) * DELTA
        dl = (p_tr - self.p_apex(k)) / (self.K * self.etabar)
        dk_dptr = a / (1.0 + a / self.eta)
        C = np.outer(DELTA / self.eta, dk_dptr * self.K * DELTA)
        return sigma, k, dl, C

    # -- one step ----------------------------------------------------------
    def step(self, state, eps_new, name="", verbose=True, want_tangent=True,
             force_region=None):
        eps_n, sig_n, k_n = state
        sig_tr = sig_n + self.Ee @ (eps_new - eps_n)
        if self.f(sig_tr, k_n) <= 0.0:
            return (eps_new, sig_tr, float(k_n)), dict(
                region="elastic", dl=0.0, C=self.Ee.copy(),
                f=self.f(sig_tr, k_n))
        region = force_region or self.region_of(sig_tr, k_n, verbose=verbose)[0]
        if region == "cone":
            R = self.cone_residual(sig_tr, k_n)
            x0 = np.concatenate([sig_tr, [k_n], [0.0]])
        else:
            R = self.apex_residual(sig_tr, k_n)
            x0 = np.concatenate([self.p_apex(k_n) * DELTA, [k_n], [0.0]])
        scale = max(1.0, float(np.max(np.abs(sig_tr))))
        x, hist = newton(R, x0, name or f"DP-{region}", scale=scale,
                         verbose=verbose)
        sigma, k, dl = x[0:6], x[6], x[7]
        C = None
        if want_tangent:
            J = jac_complex_step(R, x)
            C = (tangent_from_jacobian(J, self.Ee) if region == "cone"
                 else tangent_general(J, self.apex_dR_deps()))
        return (eps_new, sigma, float(k)), dict(
            region=region, dl=float(dl), C=C, f=float(self.f(sigma, k)),
            hist=hist)


# ---------------------------------------------------------------------------
def run(mat, legs, nstep=10, label="", verbose_last_only=True):
    state = (np.zeros(6), np.zeros(6), 0.0)
    steps, eps_prev = [], np.zeros(6)
    for leg in legs:
        for i in range(nstep):
            steps.append(eps_prev + (leg - eps_prev) * (i + 1) / nstep)
        eps_prev = leg
    info = None
    for i, e in enumerate(steps):
        last = (i == len(steps) - 1)
        state, info = mat.step(state, e, name=f"{label} step {i+1}",
                               verbose=(last or not verbose_last_only))
    return state, info, steps


def report(mat, case, path_name, legs, nstep=10):
    section(f"DruckerPrager  |  {case}  |  path = {path_name}")
    state, info, steps = run(mat, legs, nstep=nstep, label=f"{case}/{path_name}")
    eps, sigma, k = state
    print_vec(eps, "total strain eps")
    print_vec(sigma, "final stress sigma")
    print(f"    p = {mean_stress(sigma):.10f}   sqrt(J2) = {mat.q_of(sigma):.10f}"
          f"   k = {k:.10f}")
    print(f"    region = {info['region']}   f = {info['f']:.6e}   "
          f"dLambda_last = {info['dl']:.10e}")
    C = info["C"]
    print_matrix(C, "consistent tangent dsigma/deps")
    print(f"      rank(C) = {np.linalg.matrix_rank(C, tol=1e-8 * max(1.0, np.linalg.norm(C)))}")

    st = (np.zeros(6), np.zeros(6), 0.0)
    for e in steps[:-1]:
        st, _ = mat.step(st, e, verbose=False, want_tangent=False)

    def map_fn(ep):
        s, _ = mat.step(st, ep, verbose=False, want_tangent=False)
        return s[1]

    err, _ = assert_fd_tangent(C, map_fn, steps[-1], f"{case}/{path_name}",
                               h=1e-8, ref=float(np.linalg.norm(mat.Ee)))
    return dict(sigma=sigma, k=k, f=info["f"], dl=info["dl"], C=C,
                region=info["region"], fd_err=err)


# strain paths (engineering Voigt totals)
COMPRESS = [np.array([4.0e-4, 4.0e-4, -3.0e-3, 0., 0., 0.])]
COMPRESS_SHEAR = [np.array([4.0e-4, 4.0e-4, -3.0e-3, 1.5e-3, 0., 0.])]
HYDRO_TENSION = [np.array([1.2e-3, 1.2e-3, 1.2e-3, 0., 0., 0.])]


def main():
    banner("ADR-97 oracle 2/5 -- DRUCKER-PRAGER cone + apex closest-point map")
    print("Material: E = 30000, nu = 0.25, xi_c = 20, eta = 0.4 "
          "(apex at p = +50, TENSION positive).")
    print("K = 20000, G = 12000, so the exact apex-boundary slope K*etabar/G is")
    print("0.667 (associated) or 0.333 (etabar = 0.2) against the header's 0.4.")

    section("complex-step Jacobian validation (one-time warrant)")
    m0 = DP(etabar=0.2, H=500.0)
    sig_tr = m0.Ee @ np.array([4e-4, 2e-4, -3e-3, 1.5e-3, -5e-4, 2e-4])
    Rw = m0.cone_residual(sig_tr, 0.0)
    xw, _ = newton(Rw, np.concatenate([sig_tr, [0.0], [0.0]]), "warrant",
                   scale=float(np.max(np.abs(sig_tr))))
    check_complex_step_jacobian(Rw, xw, label=" (DP cone, converged point)")
    print("    NOTE: from the natural start (sigma_tr, k_n, 0) the DP cone Newton")
    print("    converges in ONE step -- the return is exactly radial, so the")
    print("    residual is affine along the solution ray.  Started away from that")
    print("    ray the quadratic descent is visible:")
    bad = np.concatenate([sig_tr * 0.75 + 8.0, [1.0], [8e-4]])
    newton(Rw, bad, "DP cone, deliberately bad start",
           scale=float(np.max(np.abs(sig_tr))), iter_gate=8)

    cases = {
        "cone, associated (etabar = eta = 0.4), perfect": (DP(), COMPRESS),
        "cone, NON-associated (etabar = 0.2), perfect":
            (DP(etabar=0.2), COMPRESS_SHEAR),
        "cone, NON-associated (etabar = 0.2), linear hardening H = 500":
            (DP(etabar=0.2, H=500.0), COMPRESS_SHEAR),
        "apex, associated, perfect": (DP(), HYDRO_TENSION),
        "apex, NON-associated (etabar = 0.2), linear hardening H = 500":
            (DP(etabar=0.2, H=500.0), HYDRO_TENSION),
    }
    results = {}
    for case, (m, legs) in cases.items():
        results[case] = report(m, case, "compress" if legs is not HYDRO_TENSION
                               else "hydro-tension", legs)

    section("closed form vs Newton-Jacobian (cone and apex)")
    for case, (m, legs) in cases.items():
        st = (np.zeros(6), np.zeros(6), 0.0)
        _, _, steps = run(m, legs, label="x", verbose_last_only=True)
        for e in steps[:-1]:
            st, _ = m.step(st, e, verbose=False, want_tangent=False)
        sN, iN = m.step(st, steps[-1], verbose=False)
        sig_tr = st[1] + m.Ee @ (steps[-1] - st[0])
        if iN["region"] == "cone":
            sC, kC, dlC, Ccf = m.cone_closed(sig_tr, st[2])
        else:
            sC, kC, dlC, Ccf = m.apex_closed(sig_tr, st[2])
        e_s = float(np.max(np.abs(sN[1] - sC)))
        e_C = rel_fro(Ccf, iN["C"]) if np.linalg.norm(iN["C"]) > 0 else \
            float(np.linalg.norm(Ccf))
        print(f"    {case[:56]:56s} max|dsigma| = {e_s:.3e}   "
              f"rel_fro(C_cf, C_newton) = {e_C:.3e}")
        assert e_s < 1e-8 and e_C < 1e-8

    section("HEADER vs EXACT apex-region test -- states they classify "
            "DIFFERENTLY")
    print("    The header's Euclidean slope is eta; the exact elastic-metric")
    print("    slope is (K*etabar + hk/eta)/G.  Any trial state whose")
    print("    (p_tr - p_apex)/q_tr ratio falls BETWEEN the two is misclassified.")
    for m, ratios in ((DP(etabar=0.2), (0.36, 0.37)),      # exact 0.333 < 0.4
                      (DP(), (0.45, 0.60))):                # exact 0.667 > 0.4
        exact_slope = (m.K * m.etabar + m.hk / m.eta) / m.G
        print(f"\n    etabar = {m.etabar}: exact slope {exact_slope:.6f}, "
              f"header slope {m.eta:.6f}")
        for ratio in ratios:
            qt = 50.0
            pt = m.p_apex(0.0) + ratio * qt
            r = np.array([1., -1., 0., 0., 0., 0.])
            r = r * qt / np.sqrt(0.5 * dot_stress(r, r))
            sig_tr = pt * DELTA + r
            print(f"    (p_tr - p_apex)/q_tr = {ratio:.3f}")
            m.region_of(sig_tr, 0.0, verbose=True)
            _, _, dl_cone, _ = m.cone_closed(sig_tr, 0.0)
            q_after = qt - m.G * dl_cone
            print(f"      what the CONE return would give: dLambda = "
                  f"{dl_cone:.6e}, sqrt(J2)_n+1 = {q_after:.6f}"
                  + ("   <-- INADMISSIBLE (negative sqrt(J2))"
                     if q_after < 0 else ""))

    banner("DRUCKER-PRAGER REFERENCE BLOCK (values C++ tests pin)")
    for case, r in results.items():
        print(f"  {case}")
        print(f"    region = {r['region']}")
        print(f"    sigma  = {np.array2string(r['sigma'], precision=10)}")
        print(f"    k = {r['k']:.10f}   f = {r['f']:.3e}   "
              f"dLambda_last = {r['dl']:.10e}")
        print(f"    trace(C) = {np.trace(r['C']):.6f}   |C|_F = "
              f"{np.linalg.norm(r['C']):.6f}   "
              f"|C-C^T|_F/|C|_F = "
              f"{np.linalg.norm(r['C'] - r['C'].T) / max(np.linalg.norm(r['C']), 1e-30):.3e}"
              f"   FD rel err = {r['fd_err']:.3e}")


if __name__ == "__main__":
    main()
