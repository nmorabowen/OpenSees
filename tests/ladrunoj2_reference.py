"""Reference (oracle) for the LadrunoJ2 combined-hardening von Mises return map.

Pure-numpy mirror of SRC/material/nD/LadrunoJ2Kernel.h — the implicit
backward-Euler J2 return map with Voce+linear isotropic hardening and Chaboche
superposed Armstrong-Frederick kinematic hardening, reduced to a single scalar
Newton on dGamma (Kobayashi-Ohno 2002).

It exists to:
  * cross-check the extracted C++ kernel WITHOUT building OpenSees
    (tests/ladrunoj2_kernel_check.cpp ↔ tests/test_ladrunoJ2_kernel_cpp.py), and
  * provide the small-strain inner for the finite-strain composition proof
    (tests/test_ladrunoJ2_finite.py) — LogStrain-over-LadrunoJ2.

This implementation is deliberately written in 3×3-tensor / einsum form (the C++
kernel uses the 6-component {00,11,22,01,12,02} tensor layout with hand-rolled
factor-2 contractions), so an indexing bug in one is unlikely to be mirrored in
the other. The consistent tangent is checked independently by finite difference
in the test, so this module does NOT transcribe the analytic tangent formula.

No OpenSees / openseespy dependency — numpy only.
"""
from __future__ import annotations

import numpy as np

_I = np.eye(3)
_ROOT23 = np.sqrt(2.0 / 3.0)

# 6-component symmetric-tensor layout used by the kernel I/O.
_REP_I = (0, 1, 2, 0, 1, 2)
_REP_J = (0, 1, 2, 1, 2, 0)


def tensor_from_voigt(v6):
    """Build a symmetric 3×3 tensor from 6 TENSOR components {00,11,22,01,12,02}
    (off-diagonal entries are true tensor components, NOT engineering)."""
    T = np.zeros((3, 3))
    for a in range(6):
        i, j = _REP_I[a], _REP_J[a]
        T[i, j] = v6[a]
        T[j, i] = v6[a]
    return T


def voigt_from_tensor(T):
    """6 tensor components {00,11,22,01,12,02} of a symmetric 3×3 tensor."""
    return np.array([T[_REP_I[a], _REP_J[a]] for a in range(6)])


# --------------------------------------------------------------------------- #
#  Hardening                                                                   #
# --------------------------------------------------------------------------- #
class Params:
    def __init__(self, K, G, sig0=0.0, Qinf=0.0, bIso=0.0, Hiso=0.0,
                 C=None, gam=None):
        self.K, self.G = K, G
        self.sig0, self.Qinf, self.bIso, self.Hiso = sig0, Qinf, bIso, Hiso
        self.C = list(C) if C is not None else []
        self.gam = list(gam) if gam is not None else []
        assert len(self.C) == len(self.gam)

    @property
    def nBack(self):
        return len(self.C)

    def yield_stress(self, p):
        return self.sig0 + self.Qinf * (1.0 - np.exp(-self.bIso * p)) + self.Hiso * p

    def yield_slope(self, p):
        return self.Qinf * self.bIso * np.exp(-self.bIso * p) + self.Hiso


# --------------------------------------------------------------------------- #
#  The return map (state in 3×3-tensor form)                                   #
# --------------------------------------------------------------------------- #
def return_map(p: Params, eps, epsP_n, ebarP_n, alpha_n):
    """Combined-hardening J2 backward-Euler return map.

    Inputs (3×3 symmetric tensors; eps = TOTAL small strain, true-tensor shear):
      epsP_n   committed plastic strain (deviatoric)
      ebarP_n  committed accumulated equivalent plastic strain
      alpha_n  list of committed backstress tensors (len == nBack)

    Returns dict: stress (3×3 Cauchy), epsP, ebarP, alpha (list), dGamma, plastic.
    """
    G, K = p.G, p.K
    alpha_n = [np.asarray(a, float) for a in alpha_n]

    tr = np.trace(eps)
    edev = eps - (tr / 3.0) * _I
    s_tr = 2.0 * G * (edev - epsP_n)

    M = s_tr - sum(alpha_n) if p.nBack else s_tr.copy()
    normM = np.sqrt(np.tensordot(M, M))
    sy0 = p.yield_stress(ebarP_n)
    f_tr = normM - _ROOT23 * sy0

    stress_scale = max(normM, _ROOT23 * sy0)
    tolY = 1.0e-12 * max(stress_scale, 1.0e-300)
    norm_floor = 1.0e-10 * max(stress_scale, 1.0e-300)

    def elastic():
        s = K * tr * _I + s_tr
        return dict(stress=s, epsP=epsP_n.copy(), ebarP=ebarP_n,
                    alpha=[a.copy() for a in alpha_n], dGamma=0.0, plastic=False)

    if f_tr <= tolY:
        return elastic()

    # scalar Newton on dG
    dG = 0.0
    for _ in range(50):
        Dk = [1.0 + _ROOT23 * p.gam[k] * dG for k in range(p.nBack)]
        M = s_tr.copy()
        Mp = np.zeros((3, 3))
        for k in range(p.nBack):
            inv = 1.0 / Dk[k]
            M -= alpha_n[k] * inv
            Mp += alpha_n[k] * (_ROOT23 * p.gam[k] * inv * inv)
        normM = np.sqrt(np.tensordot(M, M))

        theta = 2.0 * G * dG
        dtheta = 2.0 * G
        for k in range(p.nBack):
            inv = 1.0 / Dk[k]
            theta += (2.0 / 3.0) * p.C[k] * dG * inv
            dtheta += (2.0 / 3.0) * p.C[k] * inv * inv

        pbar = ebarP_n + _ROOT23 * dG
        R = (normM - theta) - _ROOT23 * p.yield_stress(pbar)
        dnormM = np.tensordot(M, Mp) / normM if normM > norm_floor else 0.0
        dR = dnormM - dtheta - (2.0 / 3.0) * p.yield_slope(pbar)
        if abs(R) <= tolY:
            break
        if dR == 0.0:
            break
        dG -= R / dR
        if dG < 0.0:
            dG = 0.0

    Dk = [1.0 + _ROOT23 * p.gam[k] * dG for k in range(p.nBack)]
    M = s_tr.copy()
    for k in range(p.nBack):
        M -= alpha_n[k] / Dk[k]
    normM = np.sqrt(np.tensordot(M, M))
    if normM <= norm_floor:
        return elastic()

    n = M / normM
    pbar = ebarP_n + _ROOT23 * dG
    epsP = epsP_n + dG * n
    alpha = [(alpha_n[k] + (2.0 / 3.0) * p.C[k] * dG * n) / Dk[k]
             for k in range(p.nBack)]
    s = s_tr - 2.0 * G * dG * n
    stress = K * tr * _I + s
    return dict(stress=stress, epsP=epsP, ebarP=pbar, alpha=alpha,
                dGamma=dG, plastic=True)


# --------------------------------------------------------------------------- #
#  A stateful small-strain wrapper (mimics an OpenSees NDMaterial), used by    #
#  the finite-strain composition proof.                                        #
# --------------------------------------------------------------------------- #
class StatefulLadrunoJ2:
    """setTrialStrain(total ε) → subtracts COMMITTED εᵖ; getStress() → Cauchy.

    Mirrors how LadrunoJ2 behaves as the inner of LogStrainNDMaterial: the
    backstress/εᵖ history lives inside this object and persists across commits.
    """
    def __init__(self, p: Params):
        self.p = p
        self.epsP = np.zeros((3, 3))
        self.ebarP = 0.0
        self.alpha = [np.zeros((3, 3)) for _ in range(p.nBack)]
        self._tr = None

    def setTrialStrain(self, eps_total):
        self._tr = return_map(self.p, eps_total, self.epsP, self.ebarP, self.alpha)

    def getStress(self):
        return self._tr["stress"]

    def getInitialTangent_compliance(self, tau):
        """Cᵉ : τ for the isotropic linear-elastic inner (used by the adaptor's
        εᵉ recovery). Closed form, independent of the return map."""
        K, G = self.p.K, self.p.G
        tr = np.trace(tau)
        return (tau - (tr / 3.0) * _I) / (2.0 * G) + (tr / (9.0 * K)) * _I

    def commitState(self):
        self.epsP = self._tr["epsP"].copy()
        self.ebarP = self._tr["ebarP"]
        self.alpha = [a.copy() for a in self._tr["alpha"]]
