"""Independent numpy oracle for the LadrunoComplexEigen reduced-pencil kernel
(ADR 46 P0).

Deliberately uses a DIFFERENT linearization from the C++ side: the C++ kernel
solves the symmetric-block generalized pencil (-B) w = lambda A w with dggev;
this oracle folds the quadratic problem into the standard companion form

    d/dt {u; v} = [[0, I], [-M^-1 K, -M^-1 C]] {u; v}

and calls numpy.linalg.eig. Agreement between the two paths therefore checks
the pencil construction AND the LAPACK call, not just a re-run of the same
recipe. scipy is not assumed (the fork's test env ships numpy only).
"""
import numpy as np


def raw_roots(M, C, K):
    """All 2p eigenvalues of (lambda^2 M + lambda C + K) z = 0, unsorted."""
    M = np.asarray(M, dtype=float)
    C = np.asarray(C, dtype=float)
    K = np.asarray(K, dtype=float)
    p = M.shape[0]
    Minv = np.linalg.inv(M)
    comp = np.zeros((2 * p, 2 * p))
    comp[:p, p:] = np.eye(p)
    comp[p:, :p] = -Minv @ K
    comp[p:, p:] = -Minv @ C
    return np.linalg.eigvals(comp)


def sorted_roots(M, C, K):
    """Raw roots sorted by (|lambda|, Re, Im) for multiset comparison."""
    lam = raw_roots(M, C, K)
    order = np.lexsort((lam.imag, lam.real, np.abs(lam)))
    return lam[order]


STRIDE = 7  # omega0, omegaD, zeta, Re, Im, kind, resid


def expand_cpp_modes(flat):
    """Expand the complexEigen -qz return (STRIDE doubles per reported mode)
    back into the raw-root multiset: an underdamped entry (kind 0)
    contributes lambda and conj(lambda); overdamped/rigid entries (kind 1/2)
    contribute their single real root."""
    flat = list(flat)
    assert len(flat) % STRIDE == 0, \
        f"expected {STRIDE} doubles per mode, got {len(flat)}"
    roots = []
    for i in range(0, len(flat), STRIDE):
        re, im, kind = flat[i + 3], flat[i + 4], int(round(flat[i + 5]))
        if kind == 0:
            roots.append(complex(re, im))
            roots.append(complex(re, -im))
        else:
            roots.append(complex(re, 0.0))
    lam = np.array(roots)
    order = np.lexsort((lam.imag, lam.real, np.abs(lam)))
    return lam[order]


def modes_table(flat):
    """The -qz return as a list of dicts, in the C++ reporting order."""
    kinds = {0: "underdamped", 1: "overdamped", 2: "rigid"}
    out = []
    for i in range(0, len(flat), STRIDE):
        out.append({
            "omega0": flat[i],
            "omegaD": flat[i + 1],
            "zeta": flat[i + 2],
            "lambda": complex(flat[i + 3], flat[i + 4]),
            "kind": kinds[int(round(flat[i + 5]))],
            "resid": flat[i + 6],
        })
    return out
