"""ADR-97 P0 oracle 1/6 -- von Mises CLOSEST-POINT return map + consistent tangent.

MIRRORS
-------
* ``SRC/material/nD/ASDPlasticMaterial3D/YieldFunctions/VonMises_YF.h``

      f = sqrt(r_ij r_ij) - SQRT_2_over_3 * k,     r = dev(sigma) - alpha

  with ``r_ij r_ij == tensor_dot_stress_like(r, r)`` (shear weight 2) and
  ``SQRT_2_over_3`` the TRUNCATED literal 0.816496580928.  NOTE the scaling:
  this is ``sqrt(2 J2)``, i.e. ``sqrt(2/3)`` times the textbook ``sqrt(3 J2)``,
  so ``k`` is the uniaxial yield stress and ``dLambda`` is scaled accordingly.
* ``PlasticFlowDirections/VonMises_PF.h`` -- ASSOCIATED: ``m == n``, the VOIGT
  gradient (shear slots x2).
* ``AllASDHardeningFunctions.h``

      LinearHardeningForScalar : h_k     = H * sqrt(2/3 * dot_strain(m, m))
      ArmstrongFrederick       : h_alpha = ha * dev(m)
                                         - cr * sqrt(2/3*dot_strain(dev m,dev m))
                                           * dev(alpha)
                                 (h_alpha := 0 once
                                  sqrt(2/3*dot_stress(dev alpha,dev alpha))
                                  >= ha/cr)

  Two convention notes carried verbatim from the header (see ADR-97 report):
  (1) there is NO 2/3 on ``ha`` -- the ``(2./3.)*ha`` variants are commented
      out at ``AllASDHardeningFunctions.h:158-160``;
  (2) ``ha * mdev`` adds an ENGINEERING-shear (doubled) quantity to the
      STRESS-like back stress ``alpha``, which is consumed by
      ``dev(sigma) - alpha`` and by ``tensor_dot_stress_like``.  On any sheared
      path the shear slots of ``alpha`` therefore grow 2x relative to the
      normal slots compared with a tensor-consistent AF law.  The oracle
      mirrors the header (``af_voigt_h=True``, the pinned reference) and also
      prints the tensor-consistent variant for contrast.

WHAT IS NEW HERE VS ``Backward_Euler``
--------------------------------------
The map below is FULLY IMPLICIT: ``n``, ``m`` and ``h`` are evaluated at the
CONVERGED ``(sigma_{n+1}, alpha_{n+1}, k_{n+1})`` and the internal variables are
updated as ``alpha_{n+1} = alpha_n + dLambda * h(sigma_{n+1}, alpha_{n+1})``.
``Backward_Euler`` instead accumulates ``sum_k deltaLambda_k * h(sigma_k)`` over
its cutting-plane iterates (``ASDPlasticMaterial3D.h:2496-2503``), which is
exact only for constant ``h`` -- see ``path_independence.py``.

Run::

    python3.12 Ladruno_implementation/adr97_oracle/cppm_vm.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from asd_common import (SQRT_2_over_3, W_STRESS, P_DEV, banner, section,
                        elastic_voigt, dev, dot_stress, dot_strain, newton,
                        jac_complex_step, check_complex_step_jacobian,
                        tangent_from_jacobian, assert_fd_tangent, rel_fro,
                        print_matrix, print_vec)

C23 = SQRT_2_over_3


class VM:
    """von Mises with linear isotropic (k) and Armstrong-Frederick (alpha) IVs."""

    def __init__(self, E=70000.0, nu=0.3, k0=30.0, H=0.0, ha=0.0, cr=0.0,
                 af_voigt_h=True):
        self.E, self.nu, self.k0, self.H = E, nu, k0, H
        self.ha, self.cr, self.af_voigt_h = ha, cr, af_voigt_h
        self.Ee = elastic_voigt(E, nu)
        self.G = E / (2.0 * (1.0 + nu))

    # -- YF / PF -----------------------------------------------------------
    def r_of(self, sigma, alpha):
        return dev(sigma) - alpha

    def norm_r(self, r):
        t = dot_stress(r, r)
        return np.sqrt(t)

    def f(self, sigma, alpha, k):
        return self.norm_r(self.r_of(sigma, alpha)) - C23 * k

    def m_of(self, sigma, alpha):
        """VOIGT gradient == flow direction (associated)."""
        r = self.r_of(sigma, alpha)
        return W_STRESS * (r / self.norm_r(r))

    # -- hardening rates ---------------------------------------------------
    def h_k(self, m):
        return self.H * np.sqrt((2.0 / 3.0) * dot_strain(m, m))

    def h_alpha(self, m, alpha):
        if self.ha == 0.0 and self.cr == 0.0:
            return np.zeros_like(alpha)
        adev = dev(alpha)
        if self.cr != 0.0:
            anorm = np.sqrt((2.0 / 3.0) * dot_stress(np.real(adev), np.real(adev)))
            if anorm >= self.ha / self.cr:            # header saturation branch
                return np.zeros_like(alpha)
        mdev = dev(m)
        mdev_eq = np.sqrt((2.0 / 3.0) * dot_strain(mdev, mdev))
        drive = mdev if self.af_voigt_h else mdev / W_STRESS
        return self.ha * drive - self.cr * mdev_eq * adev

    # -- the closest-point map --------------------------------------------
    def residual(self, sig_tr, alpha_n, k_n):
        Ee = self.Ee

        def R(x):
            sigma, alpha, k, dl = x[0:6], x[6:12], x[12], x[13]
            m = self.m_of(sigma, alpha)
            out = np.empty(14, dtype=x.dtype)
            out[0:6] = sigma - sig_tr + dl * (Ee @ m)
            out[6:12] = alpha - alpha_n - dl * self.h_alpha(m, alpha)
            out[12] = k - k_n - dl * self.h_k(m)
            out[13] = self.f(sigma, alpha, k)
            return out
        return R

    def step(self, state, eps_new, name="", verbose=True, want_tangent=True):
        """One closest-point step.  ``state = (eps, sigma, alpha, k)``."""
        eps_n, sig_n, alpha_n, k_n = state
        sig_tr = sig_n + self.Ee @ (eps_new - eps_n)
        if self.f(sig_tr, alpha_n, k_n) <= 0.0:
            return (eps_new, sig_tr, alpha_n.copy(), float(k_n)), dict(
                plastic=False, dl=0.0, C=self.Ee.copy(), f=self.f(sig_tr, alpha_n, k_n))
        R = self.residual(sig_tr, alpha_n, k_n)
        x0 = np.concatenate([sig_tr, alpha_n, [k_n], [0.0]])
        scale = max(1.0, float(np.max(np.abs(sig_tr))))
        x, hist = newton(R, x0, name or "VM", scale=scale, verbose=verbose)
        sigma, alpha, k, dl = x[0:6], x[6:12], x[12], x[13]
        C = None
        if want_tangent:
            C = tangent_from_jacobian(jac_complex_step(R, x), self.Ee)
        return (eps_new, sigma, alpha, float(k)), dict(
            plastic=True, dl=float(dl), C=C, f=float(self.f(sigma, alpha, k)),
            hist=hist)

    # -- closed-form radial return (cases (a) and (b): alpha frozen) -------
    def closed_form(self, state, eps_new):
        eps_n, sig_n, alpha_n, k_n = state
        assert self.ha == 0.0 and self.cr == 0.0, "closed form assumes no AF"
        sig_tr = sig_n + self.Ee @ (eps_new - eps_n)
        r_tr = self.r_of(sig_tr, alpha_n)
        nrm = self.norm_r(r_tr)
        f_tr = nrm - C23 * k_n
        if f_tr <= 0.0:
            return (eps_new, sig_tr, alpha_n.copy(), float(k_n)), self.Ee.copy(), 0.0
        G2 = 2.0 * self.G
        dl = f_tr / (G2 + C23 * C23 * self.H)
        n_t = r_tr / nrm                       # unit TENSOR direction
        n_v = W_STRESS * n_t                   # VOIGT gradient
        sigma = sig_tr - G2 * dl * n_t
        k = k_n + dl * self.H * C23
        PE = P_DEV @ self.Ee
        ddl_deps = (n_v @ PE) / (G2 + C23 * C23 * self.H)
        dnt_deps = (np.eye(6) - np.outer(n_t, n_v)) @ PE / nrm
        C = self.Ee - G2 * (np.outer(n_t, ddl_deps) + dl * dnt_deps)
        return (eps_new, sigma, alpha_n.copy(), float(k)), C, float(dl)


# ---------------------------------------------------------------------------
# strain paths (total ENGINEERING Voigt strain at the end of each leg)
# ---------------------------------------------------------------------------
PATHS = {
    "triaxial": [np.array([-9.0e-4, -9.0e-4, 3.0e-3, 0., 0., 0.])],
    "simple-shear": [np.array([0., 0., 0., 4.0e-3, 0., 0.])],
    "rotating-normal": [np.array([0., 0., 0., 3.0e-3, 0., 0.]),
                        np.array([0., 0., 0., 3.0e-3, 3.0e-3, 0.])],
}
NSTEP = 10


def run_path(mat, legs, nstep=NSTEP, verbose_last_only=True, label=""):
    state = (np.zeros(6), np.zeros(6), np.zeros(6), mat.k0)
    info = None
    eps_prev = np.zeros(6)
    steps = []
    for leg in legs:
        for i in range(nstep):
            steps.append(eps_prev + (leg - eps_prev) * (i + 1) / nstep)
        eps_prev = leg
    for i, eps in enumerate(steps):
        last = (i == len(steps) - 1)
        state, info = mat.step(state, eps, name=f"{label} step {i+1}",
                               verbose=(last or not verbose_last_only))
    return state, info, steps


def report(mat, case, path_name, legs):
    section(f"VonMises  |  {case}  |  path = {path_name}")
    state, info, steps = run_path(mat, legs, label=f"{case}/{path_name}")
    eps, sigma, alpha, k = state
    print_vec(eps, "total strain eps (engineering Voigt)")
    print_vec(sigma, "final stress sigma")
    print_vec(alpha, "final back stress alpha")
    print(f"    final yield stress k = {k:.12e}")
    print(f"    f(sigma, alpha, k)   = {info['f']:.6e}   (yield residual)")
    print(f"    dLambda (last step)  = {info['dl']:.12e}")
    if mat.cr != 0.0:
        anorm = np.sqrt((2.0 / 3.0) * dot_stress(dev(alpha), dev(alpha)))
        print(f"    AF saturation check: |alpha|_eq = {anorm:.6f} "
              f"vs alpha_limit = ha/cr = {mat.ha / mat.cr:.6f}  "
              f"({'UNSATURATED' if anorm < mat.ha / mat.cr else 'SATURATED'})")
    C = info["C"]
    print_matrix(C, "consistent tangent dsigma/deps")

    # FD of the oracle's OWN map: re-run the LAST step from its own start state
    state_m1 = (np.zeros(6), np.zeros(6), np.zeros(6), mat.k0)
    for e in steps[:-1]:
        state_m1, _ = mat.step(state_m1, e, verbose=False, want_tangent=False)

    def map_fn(eps_pert):
        s, _ = mat.step(state_m1, eps_pert, verbose=False, want_tangent=False)
        return s[1]

    err, _ = assert_fd_tangent(C, map_fn, steps[-1],
                               f"{case}/{path_name}", h=1e-8,
                               ref=float(np.linalg.norm(mat.Ee)))
    return dict(sigma=sigma, alpha=alpha, k=k, f=info["f"], dl=info["dl"],
                C=C, fd_err=err)


def main():
    banner("ADR-97 oracle 1/6 -- VON MISES closest-point return map")
    print("Conventions: tension-positive, Voigt [11 22 33 12 23 13], engineering")
    print("shear, VOIGT (doubled-shear) YF/PF derivatives, SQRT_2_over_3 =")
    print(f"{SQRT_2_over_3!r} (truncated literal, ASDPlasticMaterial3DGlobals.h:41)")
    print("Material: E = 70000, nu = 0.3, k0 = 30 (kPa-like), 10 steps per leg.")

    cases = {
        "(a) perfect plasticity": VM(H=0.0),
        "(b) linear isotropic hardening H = 7000": VM(H=7000.0),
        "(c) Armstrong-Frederick ha = 15000, cr = 300 (header convention)":
            VM(H=0.0, ha=15000.0, cr=300.0, af_voigt_h=True),
    }

    # one-time warrant that the complex-step Jacobian IS the analytic one
    section("complex-step Jacobian validation (one-time warrant)")
    mat = cases["(c) Armstrong-Frederick ha = 15000, cr = 300 (header convention)"]
    st = (np.zeros(6), np.zeros(6), np.zeros(6), mat.k0)
    eps1 = np.array([1e-4, -3e-5, 2e-3, 1.5e-3, -6e-4, 3e-4])
    sig_tr = mat.Ee @ eps1
    Rw = mat.residual(sig_tr, np.zeros(6), mat.k0)
    xw, _ = newton(Rw, np.concatenate([sig_tr, np.zeros(6), [mat.k0], [0.0]]),
                   "warrant", scale=float(np.max(np.abs(sig_tr))))
    check_complex_step_jacobian(Rw, xw, label=" (VM+AF, converged point)")

    results = {}
    for case, m in cases.items():
        for pname, legs in PATHS.items():
            results[(case, pname)] = report(m, case, pname, legs)

    # closed form vs Jacobian tangent for (a) and (b)
    section("closed-form consistent tangent vs Newton-Jacobian tangent  (a)/(b)")
    for case in ["(a) perfect plasticity",
                 "(b) linear isotropic hardening H = 7000"]:
        m = cases[case]
        for pname, legs in PATHS.items():
            state = (np.zeros(6), np.zeros(6), np.zeros(6), m.k0)
            _, _, steps = run_path(m, legs, label="x")
            st = (np.zeros(6), np.zeros(6), np.zeros(6), m.k0)
            for e in steps[:-1]:
                st, _ = m.step(st, e, verbose=False, want_tangent=False)
            sN, iN = m.step(st, steps[-1], verbose=False)
            sC, Ccf, dlc = m.closed_form(st, steps[-1])
            e_s = float(np.max(np.abs(sN[1] - sC[1])))
            e_C = rel_fro(Ccf, iN["C"])
            print(f"    {case:48s} {pname:16s}  max|sigma_cf - sigma_newton| = "
                  f"{e_s:.3e}   rel_fro(C_cf, C_newton) = {e_C:.3e}")
            assert e_s < 1e-9 and e_C < 1e-9

    # AF convention contrast
    section("AF hardening convention contrast (header vs tensor-consistent)")
    for flag, tag in ((True, "header  (ha*mdev, engineering/doubled shear)"),
                      (False, "tensor  (ha*mdev/W, tensor shear)")):
        m = VM(H=0.0, ha=15000.0, cr=300.0, af_voigt_h=flag)
        st, inf, _ = run_path(m, PATHS["rotating-normal"], verbose_last_only=True,
                              label="contrast")
        print(f"    {tag}")
        print_vec(st[1], "      sigma")
        print_vec(st[2], "      alpha")

    banner("VON MISES REFERENCE BLOCK (values C++ tests pin)")
    for (case, pname), r in results.items():
        print(f"  {case} | {pname}")
        print(f"    sigma  = {np.array2string(r['sigma'], precision=10)}")
        print(f"    alpha  = {np.array2string(r['alpha'], precision=10)}")
        print(f"    k      = {r['k']:.10f}    f = {r['f']:.3e}    "
              f"dLambda_last = {r['dl']:.10e}")
        print(f"    trace(C) = {np.trace(r['C']):.6f}   |C|_F = "
              f"{np.linalg.norm(r['C']):.6f}   FD rel err = {r['fd_err']:.3e}")
    print("\n  ALL von Mises FD tangent checks <= 1e-6 and all Newtons "
          "quadratic in <= 5 iterations.")


if __name__ == "__main__":
    main()
