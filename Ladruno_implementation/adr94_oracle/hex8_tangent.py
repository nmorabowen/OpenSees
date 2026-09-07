"""ADR-94 R1-A oracle: 8-node hex stiffness from a 6x6 material tangent.

The unit-cube ``stdBrick`` with an OpenSees node order and 2x2x2 Gauss.  Voigt
order is OpenSees': ``(xx, yy, zz, gxy, gyz, gzx)`` with ENGINEERING shear
(``Brick.cpp`` strain assembly, lines 1050-1065).

Run directly to print the reference numbers::

    python3.12 Ladruno_implementation/adr94_oracle/hex8_tangent.py
"""
import numpy as np

# OpenSees stdBrick node order (matches tests/test_adr94_hlist_numerics.py)
NODES = np.array([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                  (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)], dtype=float)
_G = 1.0 / np.sqrt(3.0)
GAUSS = [(a, b, c) for a in (-_G, _G) for b in (-_G, _G) for c in (-_G, _G)]


def _shape_grad(xi, eta, zeta, coords):
    """dN/dx (8x3) and det(J) at one natural point.

    N_i = 1/8 (1+s_i0 xi)(1+s_i1 eta)(1+s_i2 zeta) with s the corner signs of
    the OpenSees stdBrick node order.
    """
    s = np.array([(-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
                  (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)], dtype=float)
    n = np.array([xi, eta, zeta])
    dN = np.empty((8, 3))
    for i in range(8):
        a = 1.0 + s[i] * n                       # (1+s0*xi, 1+s1*eta, 1+s2*zeta)
        dN[i, 0] = 0.125 * s[i, 0] * a[1] * a[2]
        dN[i, 1] = 0.125 * s[i, 1] * a[0] * a[2]
        dN[i, 2] = 0.125 * s[i, 2] * a[0] * a[1]
    J = coords.T @ dN                            # 3x3
    return dN @ np.linalg.inv(J).T, float(np.linalg.det(J))


def B_matrix(dNdx):
    """6x24 strain-displacement, engineering shear, OpenSees Voigt order."""
    B = np.zeros((6, 24))
    for i in range(8):
        bx, by, bz = dNdx[i]
        c = 3 * i
        B[0, c + 0] = bx
        B[1, c + 1] = by
        B[2, c + 2] = bz
        B[3, c + 0] = by; B[3, c + 1] = bx     # gxy
        B[4, c + 1] = bz; B[4, c + 2] = by     # gyz
        B[5, c + 0] = bz; B[5, c + 2] = bx     # gzx
    return B


def hex8_K(C, coords=NODES):
    """24x24 stiffness for a constant 6x6 tangent C."""
    K = np.zeros((24, 24))
    for gp in GAUSS:
        dNdx, detJ = _shape_grad(*gp, coords)
        B = B_matrix(dNdx)
        K += B.T @ C @ B * detJ
    return K


# --------------------------------------------------------------------------
# material tangents
# --------------------------------------------------------------------------
def C_elastic(E, nu):
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    G = E / (2 * (1 + nu))
    C = np.zeros((6, 6))
    C[:3, :3] = lam
    C[0, 0] = C[1, 1] = C[2, 2] = lam + 2 * G
    C[3, 3] = C[4, 4] = C[5, 5] = G
    return C


_W = np.array([1., 1., 1., 2., 2., 2.])          # gamma -> tensor doubling


def _dev(s):
    p = (s[0] + s[1] + s[2]) / 3.0
    d = s.copy()
    d[:3] -= p
    return d


def _norm_t(s):
    """Frobenius norm of the TENSOR whose Voigt vector is s (shear once)."""
    return float(np.sqrt(s[0]**2 + s[1]**2 + s[2]**2 + 2 * (s[3]**2 + s[4]**2 + s[5]**2)))


def vm_radial_return(sig_n, deps, E, nu, sy, Hiso):
    """Exact backward-Euler radial return for ``f = ||dev s|| - sqrt(2/3)*sy``
    (the ``VonMises_YF`` arithmetic with a zero backstress).  ``deps`` is a
    strain Voigt vector with ENGINEERING shear.  Returns (sigma, dlam, sy_new).
    """
    Ce = C_elastic(E, nu)
    tr = sig_n + Ce @ deps
    s_tr = _dev(tr)
    q = _norm_t(s_tr)
    k = np.sqrt(2.0 / 3.0)
    G = E / (2 * (1 + nu))
    f = q - k * sy
    if f <= 0.0:
        return tr, 0.0, sy
    # sigma = tr - 2G*dlam*n ;  n = s_tr/q ; iv: sy += dlam * k*Hiso
    dlam = f / (2 * G + (2.0 / 3.0) * Hiso)
    n = s_tr / q
    sig = tr - 2 * G * dlam * n
    return sig, dlam, sy + k * Hiso * dlam


def vm_consistent_tangent(sig_n, deps, E, nu, sy, Hiso):
    """dsigma/depsilon of ``vm_radial_return`` (central FD, Richardson-free:
    the map is piecewise analytic and smooth in the plastic branch)."""
    C = np.zeros((6, 6))
    h = 1e-7
    for j in range(6):
        dp, dm = deps.copy(), deps.copy()
        dp[j] += h
        dm[j] -= h
        sp = vm_radial_return(sig_n, dp, E, nu, sy, Hiso)[0]
        sm = vm_radial_return(sig_n, dm, E, nu, sy, Hiso)[0]
        C[:, j] = (sp - sm) / (2 * h)
    return C


def vm_continuum_tangent(sig, E, nu, Hiso):
    """The operator ``ComputeTangentStiffness()`` builds for ``Continuum``
    (ASDPlasticMaterial3D.h:453-465): ``E - E m (n^T E)/(n^T E m - H)``.
    ``n``, ``m`` are the Voigt GRADIENT vectors (shear entries NOT doubled),
    contracted with the doubling weights as the code's Voigt algebra does."""
    Ce = C_elastic(E, nu)
    s = _dev(sig)
    q = _norm_t(s)
    n = s / q
    nv = n * _W                     # gradient-side Voigt vector
    m = nv
    H = -np.sqrt(2.0 / 3.0) * Hiso * np.sqrt(2.0 / 3.0)   # dF/dsy * dsy/dlam
    den = float(nv @ Ce @ m) - H
    return Ce - np.outer(Ce @ m, nv @ Ce) / den


def dp_consistent_tangent(sig_n, deps, E, nu, k_alpha, k_c):
    """Two-state DP reference: f = ||dev s|| + alpha*I1 - c, associated.
    FD of the exact BE return (same structure as VM)."""
    Ce = C_elastic(E, nu)

    def ret(d):
        tr = sig_n + Ce @ d
        s = _dev(tr)
        q = _norm_t(s)
        I1 = tr[0] + tr[1] + tr[2]
        f = q + k_alpha * I1 - k_c
        if f <= 0:
            return tr
        n = s / q
        nv = n * _W + k_alpha * np.array([1., 1., 1., 0., 0., 0.])
        En = Ce @ nv
        dl = f / float(nv @ En)
        for _ in range(50):          # closest point (n re-evaluated)
            sg = tr - dl * En
            s2 = _dev(sg)
            q2 = _norm_t(s2)
            I2 = sg[0] + sg[1] + sg[2]
            f2 = q2 + k_alpha * I2 - k_c
            if abs(f2) < 1e-13 * max(1.0, k_c):
                break
            n2 = s2 / q2
            nv2 = n2 * _W + k_alpha * np.array([1., 1., 1., 0., 0., 0.])
            dl += f2 / float(nv2 @ Ce @ nv)
        return tr - dl * En
    C = np.zeros((6, 6))
    h = 1e-7
    for j in range(6):
        dp_, dm_ = deps.copy(), deps.copy()
        dp_[j] += h
        dm_[j] -= h
        C[:, j] = (ret(dp_) - ret(dm_)) / (2 * h)
    return C


if __name__ == "__main__":
    np.set_printoptions(precision=6, suppress=True, linewidth=160)
    E, nu, sy, H = 70000.0, 0.3, 30.0, 7000.0
    print("--- elastic hex K sanity (E=%g nu=%g) ---" % (E, nu))
    Ke = hex8_K(C_elastic(E, nu))
    print("trace(K_elastic) =", Ke.trace(), " sum(K) (rigid-body) =", abs(Ke.sum()))
    print("volume check sum(detJ*w) =", sum(_shape_grad(*g, NODES)[1] for g in GAUSS))

    sig0 = np.zeros(6)
    deps = np.array([0.0, 0.0, -2.0e-3, 0., 0., 0.])
    sig, dlam, syn = vm_radial_return(sig0, deps, E, nu, sy, H)
    print("\n--- VM one-step radial return, deps_zz = %g ---" % deps[2])
    print("sigma =", sig, "\n dlam =", dlam, " sy ->", syn)
    Calg = vm_consistent_tangent(sig0, deps, E, nu, sy, H)
    Ccon = vm_continuum_tangent(sig, E, nu, H)
    rel = np.max(np.abs(Ccon - Calg)) / np.max(np.abs(Calg))
    print(" max|C_continuum - C_consistent| / max|C_consistent| = %.6f" % rel)
    Kalg, Kcon = hex8_K(Calg), hex8_K(Ccon)
    print(" hex K: max|K_con - K_alg|/max|K_alg| = %.6f"
          % (np.max(np.abs(Kcon - Kalg)) / np.max(np.abs(Kalg))))

    print("\n--- H1 reference: TWO different states, one specialization ---")
    e_pl = np.array([0., 0., -3.6e-3, 0., 0., 0.])      # plastic cube
    e_el = np.array([0., 0., -4.0e-4, 0., 0., 0.])      # elastic cube
    C_pl = vm_consistent_tangent(np.zeros(6), e_pl, E, nu, sy, H)
    C_el = vm_consistent_tangent(np.zeros(6), e_el, E, nu, sy, H)
    Kp, Ke2 = hex8_K(C_pl), hex8_K(C_el)
    print(" the CORRECT assembly is blockdiag(K_plastic, K_elastic);")
    print(" they differ by max|Kp-Ke|/max|Ke| = %.4f"
          % (np.max(np.abs(Kp - Ke2)) / np.max(np.abs(Ke2))))
    print(" H1 makes OpenSees assemble blockdiag(K_x, K_x) with a single K_x.")

    print("\n--- DP two-state consistent tangents (associated, alpha=0.2) ---")
    for lbl, d in (("state A", np.array([0., 0., -6.0e-3, 0., 0., 0.])),
                   ("state B", np.array([0., 0., -6.0e-3, 4.0e-3, 0., 0.]))):
        Cdp = dp_consistent_tangent(np.zeros(6), d, E, nu, 0.2, 20.0)
        print(" %s: C[0,0]=%12.4f  C[2,2]=%12.4f  C[3,3]=%12.4f"
              % (lbl, Cdp[0, 0], Cdp[2, 2], Cdp[3, 3]))
