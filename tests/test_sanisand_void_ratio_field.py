"""WP-238: the void-ratio random-field generator (Ladruno_scripts/sanisand_void_ratio_field.py).
No OpenSees binary needed (numpy only).

What is pinned:
  (a)  the field is a function of POSITION and the seed: the same point gets the same value
       whatever other points are asked for, in whatever order -- so two meshes sample one
       realisation (the property white noise per element does not have);
  (b)  refinement resolves the same field: the mean over a fixed block of a fine grid
       agrees with the coarse cell value it refines, far better than independent draws
       would, and a per-element white-noise field fails that same check (negative control);
  (c)  the statistics: zero mean and unit variance over seeds at a point, and the target
       correlation at several lags, for both kernels and for an anisotropic theta;
  (d)  different seeds are different fields; one seed is reproducible bit for bit;
  (e)  truncation: values inside [e_min, e_max], the clipped counts reported;
  (f)  bad input is refused;
  (g)  the CLI round trip.
"""
import math
import os
import sys

import pytest

np = pytest.importorskip("numpy")

_HERE = os.path.dirname(os.path.abspath(__file__))
_SCRIPTS = os.path.join(os.path.dirname(_HERE), "Ladruno_scripts")
if _SCRIPTS not in sys.path:
    sys.path.insert(0, _SCRIPTS)

import sanisand_void_ratio_field as vf  # noqa: E402


def _grid(n, L=1.0):
    c = (np.arange(n) + 0.5) * L / n
    X, Y = np.meshgrid(c, c, indexing="ij")
    return np.column_stack([X.ravel(), Y.ravel()])


# ----------------------------------------------------------------------- (a)
@pytest.mark.parametrize("kernel", vf.KERNELS)
def test_value_at_a_point_does_not_depend_on_the_other_points(kernel):
    pts = _grid(12)
    z = vf.standard_field(pts, 0.3, seed=7, kernel=kernel)
    perm = np.random.default_rng(0).permutation(len(pts))
    z_perm = vf.standard_field(pts[perm], 0.3, seed=7, kernel=kernel)
    assert np.allclose(z_perm, z[perm], rtol=0.0, atol=1e-12)
    one = vf.standard_field(pts[37:38], 0.3, seed=7, kernel=kernel)
    assert abs(one[0] - z[37]) < 1e-12
    sub = vf.standard_field(pts[::5], 0.3, seed=7, kernel=kernel)
    assert np.allclose(sub, z[::5], rtol=0.0, atol=1e-12)


# ----------------------------------------------------------------------- (b)
def test_a_refined_mesh_resolves_the_same_field_and_white_noise_does_not():
    theta, nc, r = 0.4, 8, 4                      # coarse 8x8, fine 32x32 (each coarse cell = 4x4 fine)
    zc = vf.standard_field(_grid(nc), theta, seed=3).reshape(nc, nc)
    zf = vf.standard_field(_grid(nc * r), theta, seed=3).reshape(nc, r, nc, r).mean(axis=(1, 3))
    err = np.sqrt(np.mean((zf - zc) ** 2))
    # the block mean IS the coarse value up to the cell's curvature, (h^2/24) lap Z with h = theta/3.2
    assert err < 0.1, err
    assert zf.std() > 0.9 * zc.std()              # and the variance survives the refinement
    # negative control: independent draws per element, "the same seed" on both meshes
    rng = np.random.default_rng(3)
    wc = rng.standard_normal((nc, nc))
    wf = np.random.default_rng(3).standard_normal((nc * r, nc * r)).reshape(nc, r, nc, r).mean(axis=(1, 3))
    assert np.sqrt(np.mean((wf - wc) ** 2)) > 0.5
    assert wf.std() < 0.5 * wc.std()              # and its variance falls with refinement


# ----------------------------------------------------------------------- (c)
@pytest.mark.parametrize("kernel", vf.KERNELS)
@pytest.mark.parametrize("theta", [0.5, (0.8, 0.2)])
def test_mean_variance_and_correlation(kernel, theta):
    lags = np.array([[0.0, 0.0], [0.1, 0.0], [0.25, 0.0], [0.0, 0.1], [0.0, 0.25], [0.2, 0.1], [0.6, 0.0]])
    base = np.array([0.37, -1.21])
    pts = base[None, :] + lags
    nseed = 4000
    Z = np.array([vf.standard_field(pts, theta, seed=s, kernel=kernel, n_modes=256) for s in range(nseed)])
    assert abs(Z[:, 0].mean()) < 4.0 / math.sqrt(nseed)
    assert abs(Z[:, 0].var() - 1.0) < 0.08
    got = (Z * Z[:, :1]).mean(axis=0)
    want = vf.correlation(kernel, lags, theta)
    assert np.max(np.abs(got - want)) < 0.06, (got, want)
    # the scale of fluctuation: the line integral of rho is theta (per axis)
    t = vf._theta(theta, 2)
    r = np.linspace(0.0, 40.0 * t[0], 400001)
    rho = vf.correlation(kernel, np.column_stack([r, np.zeros_like(r)]), theta)
    assert 2.0 * np.trapezoid(rho, r) == pytest.approx(t[0], rel=1e-3)


def test_the_field_is_close_to_gaussian():
    Z = np.array([vf.standard_field(np.array([[0.2, 0.7]]), 0.3, seed=s)[0] for s in range(4000)])
    assert abs(np.mean(Z ** 3)) < 0.15                       # skewness 0
    assert abs(np.mean(Z ** 4) - 3.0) < 0.3                  # kurtosis 3
    assert abs(np.mean(np.abs(Z) < 1.0) - 0.6827) < 0.03


# ----------------------------------------------------------------------- (d)
def test_seeds():
    pts = _grid(10)
    a = vf.field(pts, 0.7, 0.02, 0.3, seed=1)
    b = vf.field(pts, 0.7, 0.02, 0.3, seed=1)
    c = vf.field(pts, 0.7, 0.02, 0.3, seed=2)
    assert [float(x).hex() for x in a.values] == [float(x).hex() for x in b.values]
    assert abs(np.corrcoef(a.z, c.z)[0, 1]) < 0.9 and not np.allclose(a.values, c.values)
    assert np.allclose(a.values, 0.7 + 0.02 * a.z)


# ----------------------------------------------------------------------- (e)
def test_truncation_is_reported():
    pts = _grid(40, L=20.0)
    f = vf.field(pts, 0.70, 0.05, 0.5, seed=5, e_min=0.65, e_max=0.80)
    assert f.values.min() >= 0.65 and f.values.max() <= 0.80
    raw = 0.70 + 0.05 * f.z
    assert f.clipped_low == int((raw < 0.65).sum()) > 0
    assert f.clipped_high == int((raw > 0.80).sum()) > 0
    assert f.clipped_fraction == pytest.approx((f.clipped_low + f.clipped_high) / len(pts))
    inside = (raw >= 0.65) & (raw <= 0.80)
    assert np.array_equal(f.values[inside], raw[inside])
    g = vf.field(pts, 0.70, 0.05, 0.5, seed=5)               # no bounds: nothing clipped
    assert g.clipped_low == g.clipped_high == 0 and np.array_equal(g.values, raw)


def test_zero_std_is_the_homogeneous_field():
    f = vf.field(_grid(5), 0.6944, 0.0, 0.3, seed=9)
    assert np.all(f.values == 0.6944)


# ----------------------------------------------------------------------- (f)
@pytest.mark.parametrize("kw", [
    dict(theta=0.0), dict(theta=(0.2, 0.1, 0.3)), dict(theta=-1.0), dict(std=-0.01),
    dict(e_min=0.8, e_max=0.7), dict(e_min=0.75), dict(e_max=0.6), dict(kernel="matern"),
    dict(mean=float("nan")),
])
def test_bad_input_is_refused(kw):
    args = dict(mean=0.7, std=0.02, theta=0.3, seed=1)
    args.update(kw)
    with pytest.raises(ValueError):
        vf.field(_grid(4), **args)


def test_points_shape_is_checked():
    with pytest.raises(ValueError):
        vf.standard_field(np.zeros((5, 4)), 0.3, seed=1)
    with pytest.raises(ValueError):
        vf.apply(None, [1, 2, 3], [0.7, 0.7])


def test_three_dimensions_and_apply_calls():
    pts = np.random.default_rng(1).uniform(0, 1, (50, 3))
    f = vf.field(pts, 0.7, 0.02, (0.4, 0.4, 0.1), seed=4, e_min=0.6, e_max=0.8)
    assert f.values.shape == (50,) and f.theta == (0.4, 0.4, 0.1)

    class Rec:
        def __init__(self):
            self.calls = []

        def setParameter(self, *a):
            self.calls.append(a)
    rec = Rec()
    assert vf.apply(rec, range(101, 151), f.values) == 50
    assert rec.calls[0] == ("-val", float(f.values[0]), "-ele", 101, "voidRatio")
    assert len(rec.calls) == 50


# ----------------------------------------------------------------------- (g)
def test_cli_round_trip(tmp_path, capsys):
    pts = _grid(6)
    cen = tmp_path / "c.csv"
    out = tmp_path / "e.csv"
    with open(cen, "w") as fh:
        fh.write("tag,x,y\n")
        for i, p in enumerate(pts):
            fh.write(f"{i + 1},{float(p[0])!r},{float(p[1])!r}\n")
    vf.main([str(cen), str(out), "--mean", "0.7", "--std", "0.02", "--theta", "0.3", "0.1", "--seed", "11",
             "--emin", "0.6", "--emax", "0.8"])
    assert "36 points" in capsys.readouterr().out
    rows = [r.split(",") for r in open(out).read().split()[1:]]
    want = vf.field(pts, 0.7, 0.02, (0.3, 0.1), 11, 0.6, 0.8).values
    assert [int(r[0]) for r in rows] == list(range(1, 37))
    assert np.array_equal(np.array([float(r[1]) for r in rows]), want)
    bad = tmp_path / "bad.csv"
    bad.write_text("\n".join(["tag,x,y", "1,abc,def", ""]))
    with pytest.raises(ValueError):                          # a file with no readable point is refused
        vf.main([str(bad), str(out), "--mean", "0.7", "--std", "0.02", "--theta", "0.3", "--seed", "1"])
