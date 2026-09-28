"""WP-150 R2: mesh-orientation variants for the WP-138 footing deck (footing_ab.py).

An ill-posed (non-elliptic) problem lets its bands follow mesh lines. The B/8 and B/16 bands run straight down the
element columns at x = +-B/2 (memo section 2.1). These variants move the INTERIOR nodes of the fine zone and nothing
else. Unchanged: the surface row (so the footprint, the footing coupling and the surcharge tributaries), the fine-zone
outline (so the graded zone), and all boundaries.

  jitter:A[:SEED[:KEEP]]  every interior fine-zone node BELOW the top KEEP element rows moves by U(-A h, A h) in x
                   and in y (default A 0.2, seed 150, KEEP 1). This breaks every straight mesh line statistically, and
                   quads stay convex for A < 0.25. KEEP >= 1 leaves the surface row of elements rectangular, so the
                   initial geostatic stress at the ring keeps the regular mesh's accuracy. Distorted bilinear quads do
                   not reproduce the 1-D self-weight field element by element: KEEP 0 gives a 26 % local error in the
                   top row.
  shear:DEG        a smooth, odd, column-bending field
                   dx = tan(DEG) * Y_FINE/pi * sin(pi x / X_FINE) * sin(pi (-y) / Y_FINE).
                   Columns under the footing edges lean outward by up to ~DEG at mid-depth. It is zero at x = 0, at
                   |x| = X_FINE, at the surface and at the fine-zone bottom, and y is unchanged.

`gauss_points` gives the four 2x2 Gauss-point coordinates and the area of a (possibly distorted) bilinear quad, in the
deck's GP_XI order, so the gravity check and the census read the true geometry.
"""
import math

import numpy as np

G = 1.0 / math.sqrt(3.0)
GP_XI = ((-1, -1), (1, -1), (1, 1), (-1, 1))


def node_coords(x, y, x_fine, y_fine, h0, spec):
    """Return X, Y arrays of shape (len(x), len(y)): the perturbed coordinates of grid node (i, j)."""
    X, Y = np.meshgrid(np.asarray(x, float), np.asarray(y, float), indexing="ij")
    if not spec:
        return X, Y
    eps = 1e-9 * h0
    interior = (np.abs(X) < x_fine - eps) & (Y < -eps) & (Y > -y_fine + eps)
    kind, *rest = spec.split(":")
    if kind not in ("jitter", "shear"):
        raise ValueError(f"unknown mesh perturbation {spec!r} (jitter:A[:SEED[:KEEP]] | shear:DEG)")
    if kind == "jitter":
        amp = float(rest[0]) if rest else 0.2
        seed = int(rest[1]) if len(rest) > 1 else 150
        keep = int(rest[2]) if len(rest) > 2 else 1
        interior = interior & (Y < -keep * h0 - eps)
        if not 0.0 < amp < 0.25:
            raise ValueError("jitter amplitude must be in (0, 0.25) h to keep the quads convex")
        rng = np.random.default_rng(seed)
        dX = rng.uniform(-amp * h0, amp * h0, X.shape)
        dY = rng.uniform(-amp * h0, amp * h0, Y.shape)
        X = np.where(interior, X + dX, X)
        Y = np.where(interior, Y + dY, Y)
    elif kind == "shear":
        deg = float(rest[0]) if rest else 15.0
        dx = (math.tan(math.radians(deg)) * y_fine / math.pi
              * np.sin(math.pi * X / x_fine) * np.sin(math.pi * (-Y) / y_fine))
        X = np.where(interior, X + dx, X)
    else:
        raise ValueError(f"unknown mesh perturbation {spec!r} (jitter:A[:SEED] | shear:DEG)")
    return X, Y


def gauss_points(xs, ys):
    """xs, ys: the 4 corner coordinates in counter-clockwise order (i,j), (i+1,j), (i+1,j+1), (i,j+1)."""
    xs = np.asarray(xs, float); ys = np.asarray(ys, float)
    gps, area = [], 0.0
    for xi, et in GP_XI:
        s, t = xi * G, et * G
        N = 0.25 * np.array([(1 - s) * (1 - t), (1 + s) * (1 - t), (1 + s) * (1 + t), (1 - s) * (1 + t)])
        dNs = 0.25 * np.array([-(1 - t), (1 - t), (1 + t), -(1 + t)])
        dNt = 0.25 * np.array([-(1 - s), -(1 + s), (1 + s), (1 - s)])
        J = np.array([[dNs @ xs, dNs @ ys], [dNt @ xs, dNt @ ys]])
        detJ = float(np.linalg.det(J))
        if not detJ > 0.0:
            raise ValueError(f"non-positive Jacobian {detJ:.3e}: quad inverted by the perturbation")
        gps.append((float(N @ xs), float(N @ ys)))
        area += detJ          # 2x2 Gauss weights are 1
    return gps, area


def min_jacobian_ratio(X, Y):
    """Worst min(detJ)/max(detJ) over all elements: a distortion measure to log (1 = parallelogram)."""
    worst = 1.0
    for i in range(X.shape[0] - 1):
        for j in range(X.shape[1] - 1):
            xs = [X[i, j], X[i + 1, j], X[i + 1, j + 1], X[i, j + 1]]
            ys = [Y[i, j], Y[i + 1, j], Y[i + 1, j + 1], Y[i, j + 1]]
            d = []
            for xi, et in GP_XI:
                s, t = xi, et            # corners give the extreme Jacobians of a bilinear map
                dNs = 0.25 * np.array([-(1 - t), (1 - t), (1 + t), -(1 + t)])
                dNt = 0.25 * np.array([-(1 - s), -(1 + s), (1 + s), (1 - s)])
                d.append(float(np.linalg.det(np.array([[dNs @ xs, dNs @ ys], [dNt @ xs, dNt @ ys]]))))
            worst = min(worst, min(d) / max(d))
    return worst
