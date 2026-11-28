"""WP-238: a correlated, truncated Gaussian random field of the initial void ratio, as a
FUNCTION OF POSITION, for LadrunoSANISAND / ManzariDafalias models.

Why a function of position.  White noise drawn per element (or per Gauss point) is a
different random field on every mesh: its variance over any fixed volume falls as the mesh
is refined, so a B/8 and a B/16 run do not see the same soil.  Here the field is

    e(x) = mean + std * Z(x),      Z(x) = sqrt(2/M) * sum_j cos(k_j . x + phi_j)

with the wave vectors k_j drawn once, from the seed, from the spectral density of the
chosen correlation function, and the phases phi_j uniform in [0, 2 pi) (random Fourier
modes; Shinozuka & Deodatis 1991, Rahimi & Recht 2007).  Z has zero mean, unit variance
and the chosen correlation, and it is defined at EVERY point: two meshes sample the same
realisation, and refining a mesh only resolves it better.

Correlation functions, both written with the scale of fluctuation theta (Vanmarcke 1983:
the integral of rho over the line), per axis:

    "gaussian"   rho(r) = exp(-pi * sum_i (r_i / theta_i)^2)      smooth
    "markov"     rho(r) = exp(-2 * sqrt(sum_i (r_i / theta_i)^2)) rough (exponential)

Truncation.  A void ratio outside [e_min, e_max] is not a state of the sand; values are
CLIPPED to that range and the clipped fraction is reported.  Keep it small (std well
inside the range): clipping lowers the variance and adds point masses at the bounds.

What this does NOT do: it does not regularise localisation.  A shear band in a local
continuum is still one element wide; the field decides WHERE it starts.  Use several
seeds and report the mean and the scatter.

Use (numpy only; any Python):

    import sanisand_void_ratio_field as vf
    e = vf.field(centroids, mean=0.635, std=0.01, theta=(0.20, 0.05), seed=1,
                 e_min=0.60, e_max=0.98)                 # centroids: (n, 2) or (n, 3)
    vf.apply(ops, ele_tags, e.values)                    # one setParameter per element

`apply` must be called BEFORE any analysis step (LadrunoSANISAND refuses it later).

CLI:  python sanisand_void_ratio_field.py centroids.csv out.csv --mean 0.635 --std 0.01 \\
          --theta 0.2 0.05 --seed 1 --emin 0.60 --emax 0.98
      (centroids.csv: tag,x,y[,z] with or without a header; out.csv: tag,e)
"""
from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass

import numpy as np

KERNELS = ("gaussian", "markov")


@dataclass
class Field:
    values: np.ndarray        # the void ratio at each point, after truncation
    z: np.ndarray             # the standard field Z at each point, before scaling
    clipped_low: int          # points clipped at e_min
    clipped_high: int         # points clipped at e_max
    mean: float
    std: float
    theta: tuple
    seed: int
    kernel: str
    n_modes: int

    @property
    def clipped_fraction(self):
        return (self.clipped_low + self.clipped_high) / max(len(self.values), 1)


def _theta(theta, ndim):
    t = np.atleast_1d(np.asarray(theta, dtype=float))
    if t.size == 1:
        t = np.repeat(t, ndim)
    if t.size != ndim or not np.all(np.isfinite(t)) or np.any(t <= 0.0):
        raise ValueError(f"theta wants one positive value, or one per axis ({ndim}); got {theta!r}")
    return t


def modes(seed, theta, kernel="gaussian", n_modes=1024):
    """(k, phi): the wave vectors (n_modes, ndim) and phases of the realisation `seed`.
    Depends on the seed, theta, the kernel and n_modes ONLY -- never on the points."""
    if kernel not in KERNELS:
        raise ValueError(f"kernel wants one of {KERNELS}, got {kernel!r}")
    t = np.atleast_1d(np.asarray(theta, dtype=float))
    ndim = t.size
    rng = np.random.default_rng(int(seed))
    g = rng.standard_normal((int(n_modes), ndim))
    if kernel == "gaussian":
        # rho = exp(-pi |s|^2), s = r/theta  <->  k_s ~ N(0, 2 pi I)
        ks = math.sqrt(2.0 * math.pi) * g
    else:
        # rho = exp(-2 |s|)  <->  k_s multivariate Cauchy of scale 2: 2 g / |n|, n ~ N(0, 1)
        n = rng.standard_normal(int(n_modes))
        ks = 2.0 * g / np.abs(n)[:, None]
    phi = rng.uniform(0.0, 2.0 * math.pi, int(n_modes))
    return ks / t[None, :], phi


def standard_field(points, theta, seed, kernel="gaussian", n_modes=1024, chunk=4096):
    """Z at `points` (n, ndim): zero mean, unit variance, the chosen correlation."""
    x = np.asarray(points, dtype=float)
    if x.ndim != 2 or x.shape[1] not in (1, 2, 3):
        raise ValueError(f"points wants shape (n, 1|2|3), got {x.shape}")
    t = _theta(theta, x.shape[1])
    k, phi = modes(seed, t, kernel, n_modes)
    z = np.empty(x.shape[0])
    amp = math.sqrt(2.0 / k.shape[0])
    for i in range(0, x.shape[0], chunk):          # (chunk, n_modes) at a time
        z[i:i + chunk] = amp * np.cos(x[i:i + chunk] @ k.T + phi[None, :]).sum(axis=1)
    return z


def field(points, mean, std, theta, seed, e_min=None, e_max=None, kernel="gaussian", n_modes=1024):
    """The void-ratio field at `points`.  Returns a Field."""
    if not (math.isfinite(mean) and math.isfinite(std) and std >= 0.0):
        raise ValueError(f"mean and std want finite values, std >= 0 (got {mean}, {std})")
    if e_min is not None and e_max is not None and not (e_min < e_max):
        raise ValueError(f"e_min < e_max wanted (got {e_min}, {e_max})")
    if (e_min is not None and mean < e_min) or (e_max is not None and mean > e_max):
        raise ValueError(f"mean {mean} is outside [e_min, e_max] = [{e_min}, {e_max}]")
    x = np.asarray(points, dtype=float)
    z = standard_field(x, theta, seed, kernel, n_modes)
    e = mean + std * z
    lo = int(np.count_nonzero(e < e_min)) if e_min is not None else 0
    hi = int(np.count_nonzero(e > e_max)) if e_max is not None else 0
    if e_min is not None or e_max is not None:
        e = np.clip(e, e_min, e_max)
    return Field(values=e, z=z, clipped_low=lo, clipped_high=hi, mean=float(mean), std=float(std),
                 theta=tuple(_theta(theta, x.shape[1]).tolist()), seed=int(seed), kernel=kernel,
                 n_modes=int(n_modes))


def correlation(kernel, lag, theta):
    """The target correlation at a lag vector (or array of lag vectors)."""
    r = np.atleast_2d(np.asarray(lag, dtype=float))
    s2 = ((r / _theta(theta, r.shape[1])[None, :]) ** 2).sum(axis=1)
    out = np.exp(-math.pi * s2) if kernel == "gaussian" else np.exp(-2.0 * np.sqrt(s2))
    return out if out.size > 1 else float(out[0])


def apply(ops, ele_tags, values, name="voidRatio"):
    """One `setParameter -val e -ele tag voidRatio` per element: every Gauss point of the
    element takes that initial void ratio.  Call it before the first analysis step.
    Returns the number of elements set."""
    tags = [int(t) for t in ele_tags]
    vals = [float(v) for v in values]
    if len(tags) != len(vals):
        raise ValueError(f"{len(tags)} element tags but {len(vals)} values")
    for t, v in zip(tags, vals):
        ops.setParameter("-val", v, "-ele", t, name)
    return len(tags)


def _read_centroids(path):
    tags, pts = [], []
    with open(path, newline="") as fh:
        for n, row in enumerate(csv.reader(fh)):
            if not row:
                continue
            try:
                vals = [float(c) for c in row[1:]]
                tag = int(float(row[0]))
            except ValueError:
                if n == 0:
                    continue                  # a header line, and only the first line may be one
                raise ValueError(f"{path}: line {n + 1} is not `tag,x,y[,z]`: {row!r}") from None
            tags.append(tag)
            pts.append(vals)
    if not tags:
        raise ValueError(f"{path}: no points read")
    return tags, np.asarray(pts, dtype=float)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("centroids")
    ap.add_argument("out")
    ap.add_argument("--mean", type=float, required=True)
    ap.add_argument("--std", type=float, required=True)
    ap.add_argument("--theta", type=float, nargs="+", required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--emin", type=float)
    ap.add_argument("--emax", type=float)
    ap.add_argument("--kernel", choices=KERNELS, default="gaussian")
    ap.add_argument("--modes", type=int, default=1024)
    a = ap.parse_args(argv)
    tags, pts = _read_centroids(a.centroids)
    theta = a.theta[0] if len(a.theta) == 1 else a.theta
    f = field(pts, a.mean, a.std, theta, a.seed, a.emin, a.emax, a.kernel, a.modes)
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["tag", "e"])
        for t, v in zip(tags, f.values):
            w.writerow([t, repr(float(v))])
    print(f"{len(tags)} points: mean {f.values.mean():.5f}, std {f.values.std():.5f}, "
          f"min {f.values.min():.5f}, max {f.values.max():.5f}, clipped {f.clipped_low} low / "
          f"{f.clipped_high} high ({100 * f.clipped_fraction:.2f} %)")


if __name__ == "__main__":
    main()
