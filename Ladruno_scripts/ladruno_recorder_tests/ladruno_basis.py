"""Basis-function reconstruction for MPCO_Ladruno (schema v1, §3.3).

The whole point of the apeGmsh-native schema is that a reader can map a Gauss point
from parametric -> global coordinates for *any* element using only the self-describing
BASIS + QUADRATURE descriptor, with one basis implementation per FAMILY -- no per-class
shape-function table.

    x_global(xi) = sum_i  R_i(xi) * X_i

where X_i are the element's control points (nodes) and R_i is the declared basis
(rational if NURBS). This module is the reader-side implementation; the L3 reconstruction
test (see test_ladruno_harness.py) checks it against independently computed coordinates.

v1 implements lagrange (the legacy elements, for the parity gate) and bernstein (ready
for the planned Bezier elements). serendipity / nurbs-rational hooks are stubbed with the
shape the contract expects so they slot in when those elements land.
"""

from __future__ import annotations

import numpy as np

# ---------------------------------------------------------------------------
# Parametric-domain helpers
# ---------------------------------------------------------------------------


def _to_unit(t: float, domain: str) -> float:
    """Map a parametric coord in the declared domain to [0, 1] (Bernstein's home)."""
    if domain == "[0,1]":
        return t
    if domain == "[-1,1]":
        return 0.5 * (t + 1.0)
    raise ValueError(f"reconstruction not defined for PARAM_DOMAIN={domain!r}")


# ---------------------------------------------------------------------------
# 1-D bases
# ---------------------------------------------------------------------------


def lagrange_1d(order: int, t: float, domain: str) -> np.ndarray:
    """Equispaced Lagrange basis of given order, evaluated at t, in the declared domain.

    Returns (order+1,) values whose nodes run left->right across the domain.
    """
    lo, hi = (-1.0, 1.0) if domain == "[-1,1]" else (0.0, 1.0)
    nodes = np.linspace(lo, hi, order + 1)
    vals = np.ones(order + 1)
    for i in range(order + 1):
        for j in range(order + 1):
            if i != j:
                vals[i] *= (t - nodes[j]) / (nodes[i] - nodes[j])
    return vals


def bernstein_1d(order: int, t: float, domain: str) -> np.ndarray:
    """Bernstein basis B_{i,order}(u), u in [0,1], evaluated at t (mapped from domain)."""
    from math import comb

    u = _to_unit(t, domain)
    return np.array(
        [comb(order, i) * (u ** i) * ((1.0 - u) ** (order - i)) for i in range(order + 1)]
    )


_BASIS_1D = {"lagrange": lagrange_1d, "bernstein": bernstein_1d}


# ---------------------------------------------------------------------------
# Tensor-product assembly + control-point ordering
# ---------------------------------------------------------------------------

# Local control-point ordering per topology. For the linear cases this matches the
# OpenSees/Gmsh corner convention; higher-order ordering is pinned per the element
# contract's open item (GP_PARAM / control flattening). Listed as parametric coords in
# the [-1,1] convention; remapped per declared domain at evaluation time.
_CORNER_PARAM = {
    "line": np.array([[-1.0], [1.0]]),
    "quad": np.array([[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]]),
    "hex": np.array(
        [
            [-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
            [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1],
        ],
        dtype=float,
    ),
}


def evaluate_basis(
    family: str,
    topology: str,
    order: tuple[int, ...],
    param_domain: str,
    xi: np.ndarray,
) -> np.ndarray:
    """Return the (NUM_CTRL,) vector of basis values R_i at parametric point xi.

    For the linear corner elements this is the exact tensor-product basis in the
    canonical corner ordering. The result is *not* yet rational-normalized -- the caller
    folds in control-point weights (schema §3.3) for RATIONAL elements.
    """
    if family not in _BASIS_1D:
        raise NotImplementedError(
            f"basis family {family!r} not implemented yet "
            f"(lagrange/bernstein available; serendipity/nurbs pending their elements)"
        )
    f1d = _BASIS_1D[family]
    ndir = len(order)
    xi = np.atleast_1d(np.asarray(xi, dtype=float))

    if topology == "line":
        return f1d(order[0], xi[0], param_domain)

    if topology in ("quad", "hex") and all(o == 1 for o in order):
        # linear tensor product in the canonical corner ordering
        corners = _CORNER_PARAM[topology]
        if param_domain == "[0,1]":
            corners = 0.5 * (corners + 1.0)
        out = np.ones(len(corners))
        per_dir = [f1d(1, xi[d], param_domain) for d in range(ndir)]  # each (2,)
        node_idx = (corners > (0.5 if param_domain == "[0,1]" else 0.0)).astype(int)
        for n in range(len(corners)):
            for d in range(ndir):
                out[n] *= per_dir[d][node_idx[n, d]]
        return out

    raise NotImplementedError(
        f"tensor assembly for topology={topology!r} order={order} not implemented "
        f"in v1 (corner-linear + line covered; higher-order pinned per element contract)"
    )


def reconstruct_global(
    control_coords: np.ndarray,
    family: str,
    topology: str,
    order: tuple[int, ...],
    param_domain: str,
    xi: np.ndarray,
    control_weights: np.ndarray | None = None,
) -> np.ndarray:
    """Map a parametric point xi to global coordinates: x = sum_i R_i(xi) X_i.

    control_coords : (NUM_CTRL, ndim)
    control_weights: (NUM_CTRL,) for RATIONAL elements, else None.
    """
    R = evaluate_basis(family, topology, order, param_domain, xi)
    if control_weights is not None:
        w = np.asarray(control_weights, dtype=float)
        R = (w * R) / float(np.dot(w, R))  # rational: R_i = w_i B_i / sum_j w_j B_j
    return R @ np.asarray(control_coords, dtype=float)
