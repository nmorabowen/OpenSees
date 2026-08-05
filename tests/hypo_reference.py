"""Numpy oracle for the `-geom hypo` kernel (ADR 79 P0).

Independent reference implementation of LadrunoHypoKernel.h — the
Hughes-Winget midpoint objective strain increment integrated in the unrotated
(Green-Naghdi) material frame — plus the Voigt push-forwards. The polar
decomposition here is SVD-based (numpy), deliberately a DIFFERENT algorithm
from the kernel's Higham iteration, so agreement is evidence of correctness,
not of shared code.

Running this file as a script re-emits the frozen cross-check cases file
(tests/hypo_cases.txt) consumed by tests/hypo_kernel_check.cpp:

    py -3.12 tests/hypo_reference.py

Voigt order everywhere: {00, 11, 22, 01, 12, 20}; strain carries ENGINEERING
shear, stress plain tensor components.
"""
from __future__ import annotations

import os

import numpy as np

# Frozen H8 reference geometry used by the case generator: corner-node natural
# coordinates of the [-1,1]^3 cube (identity isoparametric map), gradients at
# the centroid: dN_a/dX = (xi_a, eta_a, zeta_a)/8. Satisfies linear
# completeness sum_a X_a (x) dN_a = I.
CORNERS = np.array([
    [-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
    [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1],
], dtype=float)
DNDX_CENTROID = CORNERS / 8.0


def polar_rotation(F: np.ndarray) -> np.ndarray:
    """R = polar(F) via SVD (F = U S Vt -> R = U Vt). Requires det F > 0."""
    if np.linalg.det(F) <= 0.0:
        raise ValueError("polar_rotation: det(F) <= 0")
    U, _, Vt = np.linalg.svd(F)
    R = U @ Vt
    if np.linalg.det(R) < 0.0:  # SVD sign ambiguity cannot occur for det F > 0,
        raise ValueError("polar_rotation: reflection")  # guard anyway
    return R


def disp_gradient(dNdX: np.ndarray, u: np.ndarray) -> np.ndarray:
    """H = sum_a u_a (x) dN_a/dX. dNdX, u: (nNodes, 3)."""
    return u.T @ dNdX


def hypo_increment(dNdX: np.ndarray, un: np.ndarray, u1: np.ndarray):
    """The kernel's hypoIncrement: returns (deps_mat6, R1, J1).

    deps_mat6 is the de-rotated (unrotated-frame) objective strain increment,
    Voigt engineering shear.
    """
    I = np.eye(3)
    Fn = I + disp_gradient(dNdX, un)
    F1 = I + disp_gradient(dNdX, u1)
    Fmid = 0.5 * (Fn + F1)
    if np.linalg.det(Fmid) <= 0.0:
        raise ValueError("hypo_increment: det(F_mid) <= 0")

    # midpoint spatial gradients and incremental displacement gradient
    dNdxbar = dNdX @ np.linalg.inv(Fmid)          # (n,3): dN_a/dxbar_j
    G = (u1 - un).T @ dNdxbar                      # G_ij
    eps_sp = 0.5 * (G + G.T)

    Rmid = polar_rotation(Fmid)
    eps_mat = Rmid.T @ eps_sp @ Rmid

    J1 = np.linalg.det(F1)
    if J1 <= 0.0:
        raise ValueError("hypo_increment: det(F1) <= 0")
    R1 = polar_rotation(F1)

    deps6 = np.array([
        eps_mat[0, 0], eps_mat[1, 1], eps_mat[2, 2],
        2.0 * eps_mat[0, 1], 2.0 * eps_mat[1, 2], 2.0 * eps_mat[2, 0],
    ])
    return deps6, R1, J1


def voigt_to_tensor_stress(s6: np.ndarray) -> np.ndarray:
    return np.array([
        [s6[0], s6[3], s6[5]],
        [s6[3], s6[1], s6[4]],
        [s6[5], s6[4], s6[2]],
    ])


def tensor_to_voigt_stress(S: np.ndarray) -> np.ndarray:
    return np.array([S[0, 0], S[1, 1], S[2, 2], S[0, 1], S[1, 2], S[2, 0]])


def push_stress6(R: np.ndarray, s6: np.ndarray) -> np.ndarray:
    """sigma' = R sigma R^T, in Voigt."""
    return tensor_to_voigt_stress(R @ voigt_to_tensor_stress(s6) @ R.T)


def bond6(R: np.ndarray) -> np.ndarray:
    """6x6 stress bond matrix M: push_stress6(R, s) == M @ s."""
    vi = [0, 1, 2, 0, 1, 2]
    vj = [0, 1, 2, 1, 2, 0]
    M = np.zeros((6, 6))
    for r in range(6):
        i, j = vi[r], vj[r]
        for c in range(6):
            p, q = vi[c], vj[c]
            m = R[i, p] * R[j, q]
            if p != q:
                m += R[i, q] * R[j, p]
            M[r, c] = m
    return M


def push_tangent6(R: np.ndarray, D: np.ndarray) -> np.ndarray:
    """D' = M D M^T (engineering-shear Voigt convention)."""
    M = bond6(R)
    return M @ D @ M.T


# --------------------------------------------------------------------------- #
#  Case emission (frozen, seeded)                                             #
# --------------------------------------------------------------------------- #

N_CASES = 16
SEED = 79001  # ADR 79 P0


def random_rotation(rng: np.random.Generator, max_angle: float = np.pi * 0.9):
    axis = rng.normal(size=3)
    axis /= np.linalg.norm(axis)
    th = rng.uniform(-max_angle, max_angle)
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)


def gen_case(rng: np.random.Generator):
    """One random case on the frozen centroid-gradient H8: a smooth committed
    config (moderate stretch+rotation) stepped by a further increment, plus a
    random symmetric tangent and stress for the push-forwards."""
    X = CORNERS
    # committed config: affine A_n = R_n (I + small sym) — guaranteed det > 0
    Sn = 0.15 * _sym(rng.normal(size=(3, 3)))
    An = random_rotation(rng, 0.8) @ (np.eye(3) + Sn)
    un = (X @ An.T) - X
    # trial config: a further rotation+stretch increment
    S1 = 0.08 * _sym(rng.normal(size=(3, 3)))
    A1 = random_rotation(rng, 0.4) @ (np.eye(3) + S1) @ An
    u1 = (X @ A1.T) - X

    D = _sym(rng.normal(size=(6, 6))) + 6.0 * np.eye(6)   # SPD-ish
    s6 = rng.normal(size=6)

    deps6, R1, J1 = hypo_increment(DNDX_CENTROID, un, u1)
    sout = push_stress6(R1, s6)
    Dout = push_tangent6(R1, D)
    return dict(dNdX=DNDX_CENTROID, un=un, u1=u1, D=D, s=s6,
                deps=deps6, R1=R1, J1=J1, sout=sout, Dout=Dout)


def _sym(A):
    return 0.5 * (A + A.T)


def _fmt(arr):
    return " ".join(f"{v:.17e}" for v in np.asarray(arr).ravel())


def emit_cases(path: str):
    rng = np.random.default_rng(SEED)
    lines = ["# hypo kernel cross-check cases, emitted by tests/hypo_reference.py",
             f"# seed {SEED}, {N_CASES} cases; all arrays row-major"]
    for _ in range(N_CASES):
        c = gen_case(rng)
        lines.append("case")
        for key in ("dNdX", "un", "u1", "D", "s", "deps", "R1", "J1", "sout", "Dout"):
            lines.append(f"{key} {_fmt(c[key])}")
        lines.append("endcase")
    with open(path, "w", encoding="ascii", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    return len(lines)


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(here, "hypo_cases.txt")
    emit_cases(out)
    print(f"emitted {N_CASES} cases -> {out}")
