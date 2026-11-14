"""Numpy oracle for the incremental hypo-kernel entry (WP-232b).

Independent reference for the entry frozen in tests/mpm_interface.md section 8:

    hypoIncrementFromF(Fn, F1)            -> (status, deps[6], R1[3x3], J1)
    hypoIncrementFromFPlaneStrain(Fn, F1) -> (status, deps[3], R1[2x2], J1)

written from the interface and the literature (Hughes & Winget 1980), never from SRC/. The polar
factor is taken from an SVD (numpy), so agreement with the C++ kernel's own polar algorithm is
evidence, not shared code.

Running this file as a script fills every `pending-oracle` value in
tests/hypo_incremental_cases.txt (values already present are left alone):

    python3.12 tests/hypo_incremental_reference.py

Voigt order {00, 11, 22, 01, 12, 20}; strain with ENGINEERING shear, stress with tensor shear.
Plane strain: [xx, yy, xy], same shear conventions.
"""
from __future__ import annotations

import os

import numpy as np

PENDING = "pending-oracle"
HERE = os.path.dirname(os.path.abspath(__file__))
CASES_PATH = os.path.join(HERE, "hypo_incremental_cases.txt")


# --------------------------------------------------------------------------- #
#  Kernel (batched: leading axes are particles)                               #
# --------------------------------------------------------------------------- #

def polar_rotation(F: np.ndarray) -> np.ndarray:
    """R of F = R U via SVD (R = U Vt). Valid for det F > 0; batched over leading axes."""
    U, _, Vt = np.linalg.svd(F)
    return U @ Vt


def tensor_to_voigt_strain(E: np.ndarray) -> np.ndarray:
    """3x3 symmetric strain -> [xx yy zz 2xy 2yz 2zx] (engineering shear)."""
    return np.stack([E[..., 0, 0], E[..., 1, 1], E[..., 2, 2],
                     2.0 * E[..., 0, 1], 2.0 * E[..., 1, 2], 2.0 * E[..., 2, 0]], axis=-1)


def voigt_to_tensor_stress(s: np.ndarray) -> np.ndarray:
    """[xx yy zz xy yz zx] (tensor shear) -> 3x3."""
    return np.stack([
        np.stack([s[..., 0], s[..., 3], s[..., 5]], -1),
        np.stack([s[..., 3], s[..., 1], s[..., 4]], -1),
        np.stack([s[..., 5], s[..., 4], s[..., 2]], -1),
    ], axis=-2)


def tensor_to_voigt_stress(S: np.ndarray) -> np.ndarray:
    return np.stack([S[..., 0, 0], S[..., 1, 1], S[..., 2, 2],
                     S[..., 0, 1], S[..., 1, 2], S[..., 2, 0]], axis=-1)


def voigt3_to_tensor(s: np.ndarray) -> np.ndarray:
    """Plane strain [xx yy xy] (tensor shear) -> 2x2."""
    return np.stack([np.stack([s[..., 0], s[..., 2]], -1),
                     np.stack([s[..., 2], s[..., 1]], -1)], axis=-2)


def tensor_to_voigt3(S: np.ndarray) -> np.ndarray:
    return np.stack([S[..., 0, 0], S[..., 1, 1], S[..., 0, 1]], axis=-1)


def hypo_increment_from_F(Fn: np.ndarray, F1: np.ndarray):
    """3D entry, batched. Returns (status[...], deps[..., 6], R1[..., 3, 3], J1[...]).

    status 0 = ok, -1 = det(F_mid) <= 0 or det(F1) <= 0 (then deps, R1, J1 are NaN for that
    entry: the C++ entry leaves its outputs untouched).
    """
    Fn = np.asarray(Fn, dtype=float)
    F1 = np.asarray(F1, dtype=float)
    Fmid = 0.5 * (Fn + F1)
    det_mid = np.linalg.det(Fmid)
    J1 = np.linalg.det(F1)
    ok = (det_mid > 0.0) & (J1 > 0.0)
    status = np.where(ok, 0, -1)
    eye = np.broadcast_to(np.eye(3), Fmid.shape)
    Fmid_safe = np.where(ok[..., None, None], Fmid, eye)
    F1_safe = np.where(ok[..., None, None], F1, eye)

    G = (F1_safe - np.where(ok[..., None, None], Fn, eye)) @ np.linalg.inv(Fmid_safe)
    eps_sp = 0.5 * (G + np.swapaxes(G, -1, -2))
    Rmid = polar_rotation(Fmid_safe)
    eps_mat = np.swapaxes(Rmid, -1, -2) @ eps_sp @ Rmid
    deps = tensor_to_voigt_strain(eps_mat)
    R1 = polar_rotation(F1_safe)

    deps = np.where(ok[..., None], deps, np.nan)
    R1 = np.where(ok[..., None, None], R1, np.nan)
    J1 = np.where(ok, J1, np.nan)
    return status, deps, R1, J1


def embed_plane_strain(F2: np.ndarray) -> np.ndarray:
    """2x2 -> [[F00, F01, 0], [F10, F11, 0], [0, 0, 1]] (interface section 8.3)."""
    F2 = np.asarray(F2, dtype=float)
    F3 = np.zeros(F2.shape[:-2] + (3, 3))
    F3[..., :2, :2] = F2
    F3[..., 2, 2] = 1.0
    return F3


def hypo_increment_from_F_plane_strain(Fn: np.ndarray, F1: np.ndarray):
    """Plane-strain entry, batched. Returns (status, deps[..., 3], R1[..., 2, 2], J1)."""
    status, d6, R3, J1 = hypo_increment_from_F(embed_plane_strain(Fn), embed_plane_strain(F1))
    deps = np.stack([d6[..., 0], d6[..., 1], d6[..., 3]], axis=-1)
    return status, deps, R3[..., :2, :2], J1


def push_stress(R: np.ndarray, s: np.ndarray) -> np.ndarray:
    """sigma = R shat R^T in Voigt; 3D (6 comps) or plane strain (3 comps), by R's size."""
    if R.shape[-1] == 3:
        S = voigt_to_tensor_stress(s)
        return tensor_to_voigt_stress(R @ S @ np.swapaxes(R, -1, -2))
    S = voigt3_to_tensor(s)
    return tensor_to_voigt3(R @ S @ np.swapaxes(R, -1, -2))


# --------------------------------------------------------------------------- #
#  Case file                                                                  #
# --------------------------------------------------------------------------- #

def parse_cases(path: str = CASES_PATH):
    """Return a list of dicts {key: list-of-tokens}; numeric conversion is left to the user."""
    cases, cur = [], None
    with open(path, encoding="ascii") as fh:
        for line in fh:
            tok = line.split()
            if not tok or tok[0].startswith("#"):
                continue
            if tok[0] == "case":
                cur = {}
            elif tok[0] == "endcase":
                cases.append(cur)
                cur = None
            else:
                cur[tok[0]] = tok[1:]
    return cases


def _arr(tokens):
    return np.array([float(t) for t in tokens])


def evaluate_case(case: dict) -> dict:
    """Run the oracle on one parsed case; returns {key: np.ndarray} for status/deps/R1/J1/scauchy."""
    dim = int(case["dim"][0])
    Fn = _arr(case["Fn"]).reshape(dim, dim)
    F1 = _arr(case["F1"]).reshape(dim, dim)
    if dim == 3:
        status, deps, R1, J1 = hypo_increment_from_F(Fn, F1)
    else:
        status, deps, R1, J1 = hypo_increment_from_F_plane_strain(Fn, F1)
    out = {"status": np.array([int(status)]), "deps": deps, "R1": R1.ravel(), "J1": np.array([J1])}
    if "shat" in case and int(status) == 0:
        out["scauchy"] = push_stress(R1, _arr(case["shat"]))
    return out


def _fmt(a) -> str:
    return " ".join(f"{v:.17e}" for v in np.asarray(a, dtype=float).ravel())


def fill_pending(path: str = CASES_PATH) -> int:
    """Replace every `pending-oracle` value in the cases file by the oracle's. Returns the count."""
    with open(path, encoding="ascii") as fh:
        lines = fh.read().split("\n")
    out, block, n = [], [], 0
    for line in lines:
        tok = line.split()
        if tok and tok[0] == "case":
            block = [line]
            continue
        if block:
            block.append(line)
            if tok and tok[0] == "endcase":
                case = {}
                for b in block[1:-1]:
                    t = b.split()
                    if t:
                        case[t[0]] = t[1:]
                res = None
                for b in block:
                    t = b.split()
                    if len(t) == 2 and t[1] == PENDING:
                        if res is None:
                            res = evaluate_case(case)
                        out.append(f"{t[0]} {_fmt(res[t[0]])}")
                        n += 1
                    else:
                        out.append(b)
                block = []
            continue
        out.append(line)
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("\n".join(out))
    return n


if __name__ == "__main__":
    print(f"filled {fill_pending()} pending values in {CASES_PATH}")
