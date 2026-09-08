"""ADR-97 P0 oracle 3/6 -- Mohr-Coulomb closest-point return in PRINCIPAL STRESS
space (Clausen, Damkilde & Andersen 2006/2007 style): return to the face, to
either edge, or to the apex, with the Koiter consistent tangent per region
pushed back to the 6D Voigt frame through the principal-direction ("T") term.

MIRRORS THE HEADER, IN A DIFFERENT BUT PROVEN-EQUIVALENT ALGEBRA
----------------------------------------------------------------
``YieldFunctions/MohrCoulomb_YF.h`` writes f in INVARIANTS::

    f = (cos(theta) - sin(theta) sin(phi)/sqrt(3)) sqrt(J2) + I1 sin(phi)/3
        - c cos(phi)

with ``theta = asin(-3 sqrt(3) J3 / (2 J2^1.5))/3`` (``typedefs.h:381``).  This
oracle verifies NUMERICALLY (``check_invariant_form``, printed below) that the
header's f is EXACTLY, to 1e-14 relative over random states, the classical
principal-stress Mohr-Coulomb function with TENSION-POSITIVE stress and
``s1 >= s2 >= s3``::

    f = 0.5 [ (s1 - s3) + (s1 + s3) sin(phi) ] - c cos(phi)

NOTE THE SCALING: this is the ``/2`` form (half the ``(s1-s3)+...`` form used
in much of the literature), so dLambda is twice the "textbook" one.  The apex is
``s1=s2=s3 = c cos(phi)/sin(phi) = c/tan(phi)`` -- exactly
``MohrCoulomb_YF::apex_stress`` (header line 230).

FLOW DIRECTION -- the header's, not the textbook's
--------------------------------------------------
``PlasticFlowDirections/MohrCoulomb_PF.h:193`` builds

    m = deviator( dg/dsigma evaluated with PHI ) + sin(psi)/3 * delta

i.e. the DEVIATORIC shape of the phi-surface plus a psi-controlled VOLUMETRIC
part -- NOT the textbook non-associated gradient (which would use psi in the
deviatoric shape too).  In principal space, with
``a_ij = 0.5[(1+sin phi) e_i - (1-sin phi) e_j]``::

    m_ij = a_ij - (sin(phi)/3) [1,1,1] + (sin(psi)/3) [1,1,1]

which reduces to ``a_ij`` when psi == phi (associated), as it must.

TWO FURTHER HEADER WARTS RECORDED HERE
--------------------------------------
* ``MohrCoulomb_PF.h:83`` reads ``double c = GET_PARAMETER_VALUE(MC_c)*M_PI/180``
  -- the COHESION is converted to radians.  Harmless in the shipped code only
  because ``c`` enters ``g`` as the additive constant ``-c cos(phi)``, which
  differentiation kills; it would be a live bug the moment ``g`` is used for
  anything but its gradient.
* ``MohrCoulomb_YF``/``_PF`` have NO edge or apex algebra: for
  ``|theta| >= 29 deg`` they silently swap in a Drucker-Prager gradient
  (``c1 = 3*(2 sin phi/(sqrt(3)(3 - sin phi))), c2 = 1, c3 = 0``), and the
  default path is a NUMERICAL central difference of f (``MC_ds > 0``).  This
  file is the reference for what the ADR-97 ``Closest_Point`` map must do
  instead.

NO NEWTON IS NEEDED
-------------------
Every Mohr-Coulomb surface is LINEAR in principal stress, so each region's
return is a closed-form linear solve.  There is therefore no residual history
to print; instead every returned state is verified against all three surface
functions and against the header's invariant form.

Run::

    python3.12 Ladruno_implementation/adr97_oracle/cppm_mc.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from asd_common import (banner, section, elastic_voigt, dev, dot_stress,
                        voigt_to_matrix, matrix_to_voigt, assert_fd_tangent,
                        print_matrix, print_vec, rel_fro)

PAIRS = {"13": (0, 2), "23": (1, 2), "12": (0, 1)}
LINE1 = ("13", "23")        # edge s1 == s2
LINE2 = ("12", "13")        # edge s2 == s3


class MC:
    def __init__(self, E=30000.0, nu=0.25, phi_deg=30.0, c=10.0, psi_deg=10.0):
        self.E, self.nu, self.c = E, nu, c
        self.phi = np.deg2rad(phi_deg)
        self.psi = np.deg2rad(psi_deg)
        self.Ee = elastic_voigt(E, nu)
        lam = E * nu / ((1 + nu) * (1 - 2 * nu))
        mu = E / (2 * (1 + nu))
        self.D3 = lam * np.ones((3, 3)) + 2 * mu * np.eye(3)
        self.k = c * np.cos(self.phi)
        self.apex = np.full(3, c * np.cos(self.phi) / np.sin(self.phi))
        sp, sq = np.sin(self.phi), np.sin(self.psi)
        self.a, self.m = {}, {}
        for key, (i, j) in PAIRS.items():
            a = np.zeros(3)
            a[i] = 0.5 * (1 + sp)
            a[j] = -0.5 * (1 - sp)
            self.a[key] = a
            self.m[key] = a - (sp / 3.0) * np.ones(3) + (sq / 3.0) * np.ones(3)
        self._setup_boundaries()

    # -- yield functions ---------------------------------------------------
    def f_pair(self, s, key):
        i, j = PAIRS[key]
        return self.a[key] @ s - self.k

    def f_principal(self, s):
        return self.f_pair(np.sort(s)[::-1], "13")

    def f_invariant(self, sig6):
        """The header's own invariant expression, for cross-validation."""
        d = dev(sig6)
        J2 = 0.5 * dot_stress(d, d)
        J3 = (d[0] * d[1] * d[2] - d[0] * d[4] ** 2 - d[1] * d[5] ** 2
              - d[2] * d[3] ** 2 + 2 * d[3] * d[4] * d[5])
        if J2 <= np.finfo(float).eps:
            th = 0.0
        else:
            x = (-3.0 * np.sqrt(3.0) * J3) / (2.0 * J2 * np.sqrt(J2))
            th = np.arcsin(min(1.0, max(-1.0, x))) / 3.0
        I1 = sig6[0] + sig6[1] + sig6[2]
        sp = np.sin(self.phi)
        return ((np.cos(th) - np.sin(th) * sp / np.sqrt(3.0)) * np.sqrt(max(J2, 0.0))
                + I1 * sp / 3.0 - self.c * np.cos(self.phi))

    # -- Clausen boundary planes ------------------------------------------
    def _edge_dir(self, line):
        A = np.array([self.a[k] for k in line])
        _, _, Vt = np.linalg.svd(A)
        ell = Vt[-1]                                # null direction of A
        if ell[0] < ell[2]:                         # orient into s1 >= s3
            ell = -ell
        return ell

    def _setup_boundaries(self):
        rp = self.D3 @ self.m["13"]                 # face return direction
        self.rp_face = rp
        self.ell = {"line1": self._edge_dir(LINE1), "line2": self._edge_dir(LINE2)}
        # a point strictly inside the face, and a trial that returns to it
        s1 = self.apex[0] - 20.0
        s3 = (s1 * (1 + np.sin(self.phi)) - 2 * self.k) / (1 - np.sin(self.phi))
        s_ref = np.array([s1, 0.5 * (s1 + s3), s3])
        assert s_ref[0] > s_ref[1] > s_ref[2]
        tr_ref = s_ref + 0.5 * rp
        self.n = {}
        self.sgn = {}
        for name, ell in self.ell.items():
            nb = np.cross(ell, rp)
            self.n[name] = nb
            self.sgn[name] = np.sign(nb @ (tr_ref - self.apex))
            assert self.sgn[name] != 0

    def region_of(self, s_tr, verbose=False):
        d = s_tr - self.apex
        p1 = self.sgn["line1"] * (self.n["line1"] @ d)
        p2 = self.sgn["line2"] * (self.n["line2"] @ d)
        if verbose:
            print(f"      boundary-plane tests (signed, face side positive): "
                  f"p_line1 = {p1:+.6f}, p_line2 = {p2:+.6f}")
        if p1 >= 0 and p2 >= 0:
            return "face", (p1, p2)
        line = "line1" if p1 < p2 else "line2"
        s_line, dl = self._return_line(s_tr, LINE1 if line == "line1" else LINE2)
        t = (s_line - self.apex) @ self.ell[line]
        if verbose:
            print(f"      line parameter t = {t:+.6f} "
                  f"({'beyond the apex -> APEX' if t < 0 else 'on the edge'})")
        return ("apex" if t < 0 else line), (p1, p2)

    # -- returns -----------------------------------------------------------
    def _return_face(self, s_tr):
        a, m = self.a["13"], self.m["13"]
        Dm = self.D3 @ m
        dl = (a @ s_tr - self.k) / (a @ Dm)
        return s_tr - dl * Dm, np.array([dl])

    def _return_line(self, s_tr, line):
        A = np.array([self.a[k] for k in line])            # 2x3
        M = np.column_stack([self.m[k] for k in line])     # 3x2
        DM = self.D3 @ M
        dl = np.linalg.solve(A @ DM, A @ s_tr - self.k)
        return s_tr - DM @ dl, dl

    def _dydx(self, region):
        if region == "apex":
            return np.zeros((3, 3))
        if region == "face":
            a, m = self.a["13"], self.m["13"]
            Dm = self.D3 @ m
            return np.eye(3) - np.outer(Dm, a) / (a @ Dm)
        line = LINE1 if region == "line1" else LINE2
        A = np.array([self.a[k] for k in line])
        M = np.column_stack([self.m[k] for k in line])
        DM = self.D3 @ M
        return np.eye(3) - DM @ np.linalg.solve(A @ DM, A)

    # -- the map -----------------------------------------------------------
    def map6(self, sig_tr6, want_tangent=False, verbose=False):
        w, Q = np.linalg.eigh(voigt_to_matrix(sig_tr6))
        idx = np.argsort(w)[::-1]                      # s1 >= s2 >= s3
        x, Q = w[idx], Q[:, idx]
        if self.f_pair(x, "13") <= 0.0:
            return sig_tr6.copy(), dict(region="elastic", dl=np.zeros(1),
                                        C=self.Ee.copy(), x=x, y=x)
        region, planes = self.region_of(x, verbose=verbose)
        if region == "face":
            y, dl = self._return_face(x)
        elif region == "apex":
            # three active surfaces: the Koiter multipliers solve
            #   x - apex = sum_k dl_k D3 m_k     (all three must be >= 0)
            M = np.column_stack([self.m[k] for k in ("13", "23", "12")])
            y = self.apex.copy()
            dl = np.linalg.solve(self.D3 @ M, x - self.apex)
        else:
            y, dl = self._return_line(x, LINE1 if region == "line1" else LINE2)
        sig6 = matrix_to_voigt(Q @ np.diag(y) @ Q.T)
        out = dict(region=region, dl=dl, planes=planes, x=x, y=y, Q=Q)
        if want_tangent:
            out["C"] = self._tangent(x, y, Q, region)
        return sig6, out

    def _tangent(self, x, y, Q, region):
        """Koiter tangent in principal space + the principal-rotation ("T") term,
        rotated to the global Voigt frame.

        ``dsigma/dsigma_tr`` is an isotropic tensor function derivative: the
        3x3 block ``dy_i/dx_j`` plus, on each shear slot (ij), the coefficient
        ``(y_i - y_j)/(x_i - x_j)`` -- with the degenerate limit
        ``dy_i/dx_i - dy_i/dx_j`` when ``x_i == x_j``.
        """
        dydx = self._dydx(region)
        T = np.zeros((6, 6))
        T[:3, :3] = dydx
        scale = max(1.0, float(np.max(np.abs(x))))
        for slot, (i, j) in ((3, (0, 1)), (4, (1, 2)), (5, (0, 2))):
            if abs(x[i] - x[j]) > 1e-9 * scale:
                T[slot, slot] = (y[i] - y[j]) / (x[i] - x[j])
            else:                                     # degenerate limit
                T[slot, slot] = dydx[i, i] - dydx[i, j]
        Rs = np.zeros((6, 6))
        for k in range(6):
            e = np.zeros(6); e[k] = 1.0
            Rs[:, k] = matrix_to_voigt(Q @ voigt_to_matrix(e) @ Q.T)
        return (Rs @ T @ np.linalg.inv(Rs)) @ self.Ee

    def step(self, state, eps_new, want_tangent=True, verbose=False):
        eps_n, sig_n = state
        sig_tr = sig_n + self.Ee @ (eps_new - eps_n)
        sig, info = self.map6(sig_tr, want_tangent=want_tangent, verbose=verbose)
        return (eps_new, sig), info


def check_invariant_form(mat, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    worst = 0.0
    for _ in range(n):
        v = rng.normal(scale=30.0, size=6)
        a = mat.f_invariant(v)
        b = mat.f_principal(np.linalg.eigvalsh(voigt_to_matrix(v)))
        worst = max(worst, abs(a - b) / max(1.0, abs(b)))
    print(f"    header invariant f vs principal-stress f over {n} random states:"
          f"  max rel diff = {worst:.3e}")
    assert worst < 1e-12
    return worst


TRIALS = {
    "face (no shear)": np.array([-10., -40., -100., 0., 0., 0.]),
    "face (sheared, exercises the rotation term)":
        np.array([-10., -40., -100., 15., -8., 5.]),
    "edge s1 == s2": np.array([-25., -25., -140., 0., 0., 0.]),
    "edge s2 == s3": np.array([10., -95., -95., 0., 0., 0.]),
    "apex (hydrostatic tension)": np.array([45., 45., 45., 0., 0., 0.]),
    "apex (slightly deviatoric)": np.array([44., 40., 36., 0., 0., 0.]),
}


def report(mat, name, sig_tr):
    section(f"MohrCoulomb  |  trial: {name}")
    print_vec(sig_tr, "trial stress sigma_tr")
    sig, info = mat.map6(sig_tr, want_tangent=True, verbose=True)
    print(f"    principal trial x = [{info['x'][0]:.8f} {info['x'][1]:.8f} "
          f"{info['x'][2]:.8f}]")
    print(f"    REGION = {info['region']}")
    print(f"    principal return y = [{info['y'][0]:.8f} {info['y'][1]:.8f} "
          f"{info['y'][2]:.8f}]")
    print_vec(sig, "returned stress sigma")
    print("    dLambda = [" + " ".join(f"{v:.10e}" for v in info["dl"]) + "]")
    for key in ("13", "23", "12"):
        print(f"      f_{key}(y) = {mat.f_pair(np.sort(info['y'])[::-1], key):+.6e}")
    print(f"      f_invariant(sigma) = {mat.f_invariant(sig):+.6e}   "
          f"(header's own expression)")
    if info["region"] == "apex":
        print("      NOTE: at the apex the three multipliers above are the LEAST"
              "-SQUARES\n      decomposition of sigma_tr - sigma_apex onto "
              "{D3 m_13, D3 m_23, D3 m_12}.\n      For NON-associated flow "
              "(psi < phi) they need NOT all be positive: the\n      cone of "
              "return directions no longer contains the hydrostatic direction,"
              "\n      which is exactly why the apex region must be defined by "
              "the BOUNDARY\n      PLANES (Clausen) and not by an active-set "
              "search on dLambda >= 0.")
    else:
        assert np.all(info["dl"] >= -1e-12), "Koiter multiplier negative"
    assert mat.f_principal(info["y"]) <= 1e-9, "returned state inadmissible"
    C = info["C"]
    print_matrix(C, "consistent tangent dsigma/deps")
    print(f"      rank(C) = "
          f"{np.linalg.matrix_rank(C, tol=1e-8 * max(1.0, np.linalg.norm(C)))}")

    eps_tr = np.linalg.solve(mat.Ee, sig_tr)
    st = (np.zeros(6), np.zeros(6))

    def map_fn(ep):
        s, _ = mat.step(st, ep, want_tangent=False)
        return s[1]

    h = 1e-8
    err, _ = assert_fd_tangent(C, map_fn, eps_tr, name, h=h,
                               ref=float(np.linalg.norm(mat.Ee)))
    return dict(sigma=sig, y=info["y"], region=info["region"], dl=info["dl"],
                C=C, fd_err=err)


def main():
    banner("ADR-97 oracle 3/6 -- MOHR-COULOMB principal-space return (Clausen)")
    mat = MC()
    print(f"Material: E = 30000, nu = 0.25, phi = 30 deg, psi = 10 deg "
          f"(NON-associated), c = 10.")
    print(f"Apex (tension-positive) = c/tan(phi) = {mat.apex[0]:.10f} on the "
          f"hydrostatic axis.")
    section("warrant: the header's invariant f IS the principal-stress f")
    check_invariant_form(mat)
    print("    edge directions (from the null space of the two active gradients):")
    print(f"      line1 (s1 == s2): {np.array2string(mat.ell['line1'], precision=8)}")
    print(f"      line2 (s2 == s3): {np.array2string(mat.ell['line2'], precision=8)}")

    results = {}
    for name, s in TRIALS.items():
        results[name] = report(mat, name, s)

    section("region coverage")
    seen = sorted({r["region"] for r in results.values()})
    print(f"    regions exercised: {seen}")
    for want in ("face", "line1", "line2", "apex"):
        assert want in seen, f"no trial state reached region {want}"
    print("    all four regions (face, both edges, apex) covered.")

    section("associated check: psi == phi must give a SYMMETRIC face tangent")
    ma = MC(psi_deg=30.0)
    _, ia = ma.map6(TRIALS["face (no shear)"], want_tangent=True)
    Ca = ia["C"]
    asym = np.linalg.norm(Ca - Ca.T) / np.linalg.norm(Ca)
    print(f"    associated face tangent asymmetry |C-C^T|_F/|C|_F = {asym:.3e}")
    assert asym < 1e-12
    _, inas = mat.map6(TRIALS["face (no shear)"], want_tangent=True)
    Cn = inas["C"]
    print(f"    NON-associated (psi = 10 deg) asymmetry            = "
          f"{np.linalg.norm(Cn - Cn.T) / np.linalg.norm(Cn):.3e}")

    banner("MOHR-COULOMB REFERENCE BLOCK (values C++ tests pin)")
    for name, r in results.items():
        print(f"  {name}")
        print(f"    region = {r['region']}   dLambda = "
              f"{np.array2string(r['dl'], precision=10)}")
        print(f"    y      = {np.array2string(r['y'], precision=10)}")
        print(f"    sigma  = {np.array2string(r['sigma'], precision=10)}")
        print(f"    trace(C) = {np.trace(r['C']):.6f}   |C|_F = "
              f"{np.linalg.norm(r['C']):.6f}   rank = "
              f"{np.linalg.matrix_rank(r['C'], tol=1e-8 * max(1.0, np.linalg.norm(r['C'])))}"
              f"   FD rel err = {r['fd_err']:.3e}")
    print("\n  FD stencil: central difference, h = 1e-8 (relative), applied to the")
    print("  TOTAL strain of the single step; every trial sits well inside its")
    print("  region (see the printed boundary-plane margins), so the stencil does")
    print("  not cross a region boundary.")


if __name__ == "__main__":
    main()
