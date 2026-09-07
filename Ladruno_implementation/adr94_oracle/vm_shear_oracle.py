"""ADR-94 wp/94c -- closed-form von Mises radial return, in the framework's own
Voigt convention.  Reference for ``tests/test_adr94c_numerics.py``.

WHY THIS EXISTS
---------------
ADR-94 B5: ``ASDPlasticMaterial3D`` stores symmetric tensors as 6-vectors with
the shear slot held ONCE (``v11 v22 v33 v12 v23 v13``), and the review found the
catalogue split on what a yield-function gradient in that storage means.
``VonMises_YF``/``_PF`` returned the bare TENSOR derivative ``df/dsigma_12``;
``MohrCoulomb``/``HoekBrown``/``MCTC``/``StiffSoil`` returned the VOIGT one
``df/dv_12 = 2 df/dsigma_12``; and one consumption site doubled the shear terms a
second time.  wp/94c makes the convention VOIGT everywhere, which is the one the
framework actually needs, because both

    TrialPlastic_Strain += dLambda * m        (m must carry ENGINEERING shear)
    Eelastic * m                              (shear row of E is G, not 2G)

are plain Voigt operations.

WHAT THAT MAKES THE ALGEBRA
---------------------------
With ``r = dev(sigma) - alpha`` and ``||r|| = sqrt(r_ij r_ij)``,

    f  = ||r|| - c * sigma_y ,           c = SQRT_2_over_3
    n  = m = Voigt gradient of f = r/||r|| with the three shear slots DOUBLED

so ``E * m`` is exactly ``2G * (r/||r||)`` in raw Voigt storage (normal slots
``2G n_ii``, shear slots ``G * 2 n_ij``), the stress correction stays radial, and
the consistency condition is LINEAR in dLambda:

    ||r_trial|| - 2G*dLambda - c*(sigma_y0 + dLambda*H*c) = 0
    dLambda = (||r_trial|| - c*sigma_y0) / (2G + c*c*H)

The hardening rate is ``h = H * sqrt(2/3 * m_ij m_ij)`` evaluated with the
ENGINEERING contraction (wp/94c fixed that too), and for this m the tensor norm
``m_ij m_ij`` is exactly 1, so ``h = H*c`` on every stress path -- including the
sheared ones, where the pre-94c ``m.dot(m)`` gave the wrong number.

So a plastic SIMPLE-SHEAR step is a closed-form check on the whole convention:
gradient, flow direction, plastic modulus and hardening rate all appear, and any
stray factor of 2 in any of them shows up immediately.  On the pre-94c build the
same comparison is off by O(1).

CAVEAT ON THE CONSTANT
----------------------
``ASDPlasticMaterial3DGlobals.h`` defines ``SQRT_2_over_3 = 0.816496580928``, a
TRUNCATED literal, not ``sqrt(2.0/3.0)``.  The difference is 3e-13 relative -- it
would dominate a 1e-10 comparison -- so this oracle uses the same literal.
"""
import numpy as np

# Verbatim from SRC/material/nD/ASDPlasticMaterial3D/ASDPlasticMaterial3DGlobals.h:41
SQRT_2_over_3 = 0.816496580928


def elastic_voigt(E, nu):
    """Isotropic 6x6 in the framework's Voigt storage: strain slots 3..5 are
    ENGINEERING shear, so those rows/columns carry G (not 2G)."""
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    G = E / (2.0 * (1.0 + nu))
    D = np.zeros((6, 6))
    D[:3, :3] = lam
    D[0, 0] = D[1, 1] = D[2, 2] = lam + 2.0 * G
    D[3, 3] = D[4, 4] = D[5, 5] = G
    return D


def deviator(v):
    p = (v[0] + v[1] + v[2]) / 3.0
    s = v.astype(float).copy()
    s[0] -= p
    s[1] -= p
    s[2] -= p
    return s


def tensor_norm_stress(v):
    """sqrt(a_ij a_ij) for a STRESS-like Voigt vector (shear stored once)."""
    return float(np.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2
                         + 2.0 * (v[3] ** 2 + v[4] ** 2 + v[5] ** 2)))


def radial_return(eps_hist, E, nu, sigma_y0, H_iso):
    """Drive the closed-form von Mises map with linear isotropic hardening over
    a history of TOTAL engineering-Voigt strains (the rows OpenSees reports as
    ``eleResponse(ele, "strains")``), starting from a zero state.

    Returns the committed stress after each row, as an (n, 6) array in the same
    Voigt storage.
    """
    D = elastic_voigt(E, nu)
    G = E / (2.0 * (1.0 + nu))
    c = SQRT_2_over_3

    sig = np.zeros(6)
    sy = float(sigma_y0)
    eps_prev = np.zeros(6)
    out = []
    for eps in np.asarray(eps_hist, dtype=float):
        sig_tr = sig + D @ (eps - eps_prev)
        s_tr = deviator(sig_tr)
        nrm = tensor_norm_stress(s_tr)
        f = nrm - c * sy
        if f <= 0.0 or nrm == 0.0:
            sig = sig_tr
        else:
            dlam = f / (2.0 * G + c * c * H_iso)
            n_t = s_tr / nrm                 # unit TENSOR direction
            sig = sig_tr - 2.0 * G * dlam * n_t
            sy = sy + dlam * H_iso * c
        eps_prev = eps
        out.append(sig.copy())
    return np.array(out)
