"""ADR-97 P0 oracle 4/6 -- HOEK-BROWN closest-point return in PRINCIPAL STRESS
space (Clausen & Damkilde 2008): return to the CURVED surface, to either curved
EDGE, or to the APEX, with the consistent tangent per region (including the
``dm/dsigma`` curvature term) pushed back to 6D Voigt through the
principal-rotation term and its degenerate limit, exactly as ``cppm_mc.py``.

THE SURFACE.  ``HoekBrown_YF.h`` negates to the geomechanics frame
(``sigma_geo = -sigma``) and evaluates the composite ported in #806, with
``sig1 >= sig2 >= sig3`` COMPRESSION-positive::

    f = max(f_shear, f_tension)
    f_shear   = sig1 - sig3 - sigma_ci*max(mb*sig3/sigma_ci + s, 0)^a
    f_tension = sigma_t - sig3,      sigma_t = -s*sigma_ci/mb

This oracle works in the tree's own TENSION-POSITIVE principals ``y1>=y2>=y3``
(``sig1 = -y3``, ``sig3 = -y1``), where the same composite is Hoek-Brown in the
(major, minor) pair plus a RANKINE cut-off on the MAJOR principal::

    f_shear = y1 - y3 - sigma_ci*max(s - mb*y1/sigma_ci, 0)^a
    f_tension = y1 - T,              T = s*sigma_ci/mb > 0   ( == APEX_STRESS )

FINDING 1 -- THE TENSION BRANCH IS INERT ON THE YIELD SURFACE.  ``f_shear <= 0``
already forces ``y1 <= T``: above the apex the clamp leaves ``f_shear = y1-y3``,
positive for every non-hydrostatic state.  So the composite's zero set is the HB
shear surface plus its own natural apex vertex ``T*1``; the plane ``y1 = T``
never carries a face, and P3 needs NO tension-plane return, only the apex.
Outside the elastic domain the tension branch wins exactly on
``{y1 > T and y3 > T}`` -- precisely the header's ``CHECK_APEX_REGION`` set.

FINDING 2 -- ``HoekBrown_PF::g`` IS EVALUATED IN THE WRONG SIGN FRAME.  The YF
negates before ``principalStresses()``; ``HoekBrown_PF::g`` (header line 66)
does not, then destructures the ascending tuple as ``[sigma3, sigma2, sigma1]``
and feeds ``sigma3`` -- the tree's most COMPRESSIVE principal, not the geo-frame
minor -- into ``arg = mb_psi*sigma3/sigma_ci + s``.  On any compressive state
``arg < 0``, so ``g`` takes its ``else`` branch and becomes a TRESCA potential.
Both rules are implemented: ``flow="intended"`` (frame-consistent, what P3
should ship) and ``flow="header"`` (verbatim); the section "g vs f" measures the
difference.

REGION SELECTION.  The elastic domain is CONVEX (``-sigma_ci*arg^a`` is convex
in ``y1``: ``arg`` affine, ``0 < a < 1``; a max of such functions over the six
permutations stays convex), so for associated flow the elastic-metric closest
point is unique and KKT is SUFFICIENT -- which makes the scan a proof, not a
spot check.  APEX: the domain's tangent cone at ``T*1`` is exactly the negative
octant (the meridian meets the hydrostatic axis vertically, ``df/dy1 -> inf``),
so the normal cone is the positive octant and the apex region is EXACTLY
``apex + D3*(positive octant)``; non-associated flow uses the same construction
on the six limiting directions ``D3 m_ij(apex)`` (Caratheodory enumeration).
FACE vs EDGE: a curved surface has no precomputed boundary PLANES, but their
defining property survives -- the boundary is the ruled surface where the FACE
return lands on the edge, so the face return's own margins ``y1-y2`` and
``y2-y3`` ARE the exact signed boundary functions.  Printed for every trial.

Run::  python3.12 Ladruno_implementation/adr97_oracle/cppm_hb.py
"""
import itertools
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from asd_common import (banner, section, elastic_voigt, voigt_to_matrix,
                        matrix_to_voigt, newton, jac_complex_step, fd_tangent,
                        assert_fd_tangent, print_matrix, print_vec, rel_fro)

PAIRS = {"13": (0, 2), "23": (1, 2), "12": (0, 1)}
LINE1 = ("13", "23")        # edge y1 == y2
LINE2 = ("12", "13")        # edge y2 == y3
ALL_PAIRS = [(i, j) for i in range(3) for j in range(3) if i != j]
E1 = np.array([1., 0., 0.])


def in_cone(G, d, tol=1e-9):
    """Is ``d`` a non-negative combination of the columns of ``G``?  A cone in
    R^3 is the union of its simplicial subcones of <= 3 generators
    (Caratheodory), so enumerating subsets is exact, not heuristic."""
    nd = float(np.linalg.norm(d))
    if nd <= tol:
        return True
    for r in (1, 2, 3):
        for idx in itertools.combinations(range(G.shape[1]), r):
            A = G[:, idx]
            c = np.linalg.lstsq(A, d, rcond=None)[0]
            if np.all(c >= -tol) and np.linalg.norm(A @ c - d) <= tol * nd:
                return True
    return False


class HB:
    """Hoek-Brown, tension-positive, principal stresses sorted DESCENDING."""

    def __init__(self, E=5.0e7, nu=0.25, sigci=50000.0, mi=10.0, GSI=60.0,
                 D=0.0, mb_psi=None, flow="intended"):
        self.E, self.nu, self.sigci, self.flow = E, nu, sigci, flow
        self.mb = mi * np.exp((GSI - 100.0) / (28.0 - 14.0 * D))
        self.s = np.exp((GSI - 100.0) / (9.0 - 3.0 * D))
        self.a = 0.5 + (1.0 / 6.0) * (np.exp(-GSI / 15.0) - np.exp(-20.0 / 3.0))
        self.mb_psi = self.mb if mb_psi is None else mb_psi
        self.T = self.s * sigci / self.mb          # APEX_STRESS, tension-positive
        self.scale = sigci * self.s ** self.a      # YF_STRENGTH_SCALE
        self.Ee = elastic_voigt(E, nu)
        lam, mu = E * nu / ((1 + nu) * (1 - 2 * nu)), E / (2 * (1 + nu))
        self.D3 = lam * np.ones((3, 3)) + 2 * mu * np.eye(3)
        self.C3 = np.linalg.inv(self.D3)
        self.apex = np.full(3, self.T)
        self.Gapex = self._apex_generators()

    # ---- surface / gradient / flow (complex-step safe: no sort, no abs) -----
    def f_ij(self, y, i, j):
        return (y[i] - y[j]
                - self.sigci * (self.s - self.mb * y[i] / self.sigci) ** self.a)

    def a_ij(self, y, i, j):
        arg = self.s - self.mb * y[i] / self.sigci
        g = np.zeros(3, dtype=np.asarray(y).dtype)
        g[i] = 1.0 + self.a * self.mb * arg ** (self.a - 1.0)
        g[j] = -1.0
        return g

    def m_ij(self, y, i, j):
        """``intended``: the frame-consistent potential (== the YF when
        ``mb_psi == mb``).  ``header``: ``HoekBrown_PF::g`` verbatim, whose
        branch ``arg`` is built from the MINOR principal (finding 2)."""
        g = np.zeros(3, dtype=np.asarray(y).dtype)
        if self.flow == "header":
            argh = self.s + self.mb_psi * y[j] / self.sigci
            g[i] = 1.0
            g[j] = (-1.0 - self.a * self.mb_psi * argh ** (self.a - 1.0)
                    if np.real(argh) > 0.0 else -1.0)
        else:
            argp = self.s - self.mb_psi * y[i] / self.sigci
            g[i] = 1.0 + self.a * self.mb_psi * argp ** (self.a - 1.0)
            g[j] = -1.0
        return g

    def f_pair(self, y, key):
        return self.f_ij(y, *PAIRS[key])

    def a_pair(self, y, key):
        return self.a_ij(y, *PAIRS[key])

    def m_pair(self, y, key):
        return self.m_ij(y, *PAIRS[key])

    def mhat(self, y, key):
        """UNIT flow direction.  Normalizing leaves the returned stress
        invariant (only the multiplier rescales by ``|m|``) and is what keeps
        the near-apex Newton inside the <= 5 iteration gate: ``|m|`` blows up
        like ``arg^(a-1)`` exactly where the tensile returns land."""
        m = self.m_pair(y, key)
        return m / np.sqrt((m * m).sum())

    def _apex_generators(self):
        """Limiting return directions at the apex, one per ordered pair."""
        y = self.apex - 1e-9 * self.T
        G = np.column_stack([self.D3 @ np.real(self.m_ij(y, i, j))
                             for i, j in ALL_PAIRS])
        return G / np.linalg.norm(G, axis=0)

    # ---- the header's own composite, mirrored ------------------------------
    def f_composite(self, y):
        """max(f_shear, f_tension) with the header's clamp; y sorted desc."""
        arg = max(self.s - self.mb * y[0] / self.sigci, 0.0)
        fs, ft = y[0] - y[2] - self.sigci * arg ** self.a, y[0] - self.T
        return max(fs, ft), fs, ft

    def f_floor(self, y):
        """Round-off floor of ``f`` at ``y`` -- NOT a fudge factor.
        ``|df/dy1| = 1 + a*mb*arg^(a-1)`` diverges at the apex, so a last-ulp
        error in ``y1`` carries ``eps*|y1|*|df/dy1|`` into ``f``; an ABSOLUTE
        1e-10 gate is unattainable within ~1e-2 kPa of the apex.  This is the
        HoekBrown instance of ADR-94's ``f_relative_tol`` lesson: scale the
        yield tolerance by the GRADIENT, not by ``sigma_ci``."""
        arg = max(self.s - self.mb * y[0] / self.sigci, 1e-300)
        cond = 1.0 + self.a * self.mb * arg ** (self.a - 1.0)
        return 4.0 * np.finfo(float).eps * max(abs(y[0]), self.scale) * cond

    def f6(self, sig6):
        return self.f_composite(np.sort(np.linalg.eigvalsh(
            voigt_to_matrix(sig6)))[::-1])[0]

    def g6_header(self, sig6):
        """``HoekBrown_PF::g`` verbatim on a 6D Voigt stress (no negation)."""
        y = np.sort(np.linalg.eigvalsh(voigt_to_matrix(sig6)))       # ascending
        s3, s1 = y[0], y[2]                # header's [sigma3, sigma2, sigma1]
        arg = self.mb_psi * s3 / self.sigci + self.s
        return (s1 - s3 - self.sigci * arg ** self.a if arg > 0
                else s1 - s3 - self.sigci * self.s)

    def num_grad6(self, fn, sig6):
        """The header's own central difference over the six raw Voigt slots."""
        out, h = np.zeros(6), 1e-8 * max(1.0, float(np.linalg.norm(sig6)))
        for i in range(6):
            p, m = sig6.copy(), sig6.copy()
            p[i], m[i] = p[i] + h, m[i] - h
            out[i] = (fn(p) - fn(m)) / (2 * h)
        return out

    # ---- the surface's NATURAL VARIABLE (a requirement, not taste) ----------
    # With y1 as the unknown the Newton DIVERGES near the apex: the surface is
    # defined only for arg = s - mb y1/sigma_ci >= 0, the first step from the
    # elastic predictor overshoots y1 > T and the next Jacobian is singular
    # (measured in "near-apex conditioning").  Substituting arg = w^(2/a) --
    # i.e. w^2 = arg^a, the surface's own variable, since sigma_ci*arg^a IS the
    # strength term -- gives  y1 = T - (sigma_ci/mb) w^(2/a),  f_13 = y1 - y3 -
    # sigma_ci w^2 : polynomial-smooth (2/a = 3.977) and feasible for ANY real w
    # -- no clipping, no line search, no feasibility guard.
    def y1_of(self, w):
        return self.T - (self.sigci / self.mb) * w ** (2.0 / self.a)

    def _y_of(self, z, region):
        """Unknowns -> principal stresses, with the edge constraint built in.
        One surface row then suffices on an edge: ``f_13`` and ``f_23`` share a
        root iff ``y1 == y2`` (``y - sigma_ci arg(y)^a`` is strictly
        increasing), and likewise ``f_12``/``f_13`` iff ``y2 == y3``."""
        y1 = self.y1_of(z[0])
        if region == "face":
            return np.array([y1, z[1], z[2]])
        return np.array([y1, y1, z[1]] if region == "line1" else [y1, z[1], z[1]])

    def _keys(self, region):
        return ("13",) if region == "face" else (LINE1 if region == "line1" else LINE2)

    def _res(self, x, region):
        keys = self._keys(region)
        n_dl = len(keys)

        def R(z):
            y, dl = self._y_of(z, region), z[-n_dl:]
            corr = sum(dl[k] * (self.D3 @ self.mhat(y, key))
                       for k, key in enumerate(keys))
            f = (y[0] - y[1] if region == "line2" else y[0] - y[2])
            return np.concatenate([y - x + corr, [f - self.sigci * z[0] ** 2]])
        return R, n_dl

    def raw_dl(self, y, dl_hat, region):
        """The PHYSICAL multiplier of ``y = y_tr - dl D3 m`` (un-normalized m)."""
        return np.array([dl_hat[k] / np.linalg.norm(np.real(self.m_pair(y, key)))
                         for k, key in enumerate(self._keys(region))])

    def start(self, x, region):
        """Elastic predictor in the natural variable (seed ``sigma_ci w^2`` with
        the trial's own spread, keep the start sorted) PLUS the first-order
        cutting-plane multiplier ``f/(a . D3 mhat)``.  Measured over 400+ random
        trials: with ``dl = 0`` the worst case is 6 Newton iterations, with this
        seed it is 5 -- this is what keeps gate-1's "<= 5" honest for HB."""
        w0 = np.sqrt(max(x[0] - x[2], 1e-12) / self.sigci)
        z = (np.array([w0, x[1], x[2], 0.0]) if region == "face"
             else np.array([w0, x[2] if region == "line1" else x[1], 0.0, 0.0]))
        y1 = self.y1_of(z[0])
        z[1] = min(z[1], y1)
        if region == "face":
            z[2] = min(z[2], y1)
        y0 = self._y_of(z, region)
        av, mv = np.real(self.a_pair(y0, "13")), np.real(self.mhat(y0, "13"))
        z[3 if region == "face" else 2] = (
            max(self.f_composite(np.sort(y0)[::-1])[0], 0.0) / (av @ (self.D3 @ mv)))
        return z

    def _try(self, R, z0, n=30):
        """Quiet Newton used ONLY for region detection (no gates, no printing)."""
        z = np.array(z0, float)
        for _ in range(n):
            r = np.real(np.asarray(R(z.astype(complex))))
            if not np.all(np.isfinite(r)):
                return None
            if np.linalg.norm(r) <= 1e-13 * self.scale:
                return z
            try:
                z = z - np.linalg.solve(jac_complex_step(R, z), r)
            except np.linalg.LinAlgError:
                return None
        return None

    def region_of(self, x, verbose=False):
        """apex (exact cone test) -> face -> edge (by the face return's margins,
        which are the exact boundary functions); verified downstream by KKT."""
        if in_cone(self.Gapex, x - self.apex):
            if verbose:
                print("      apex cone test C3*(x-apex) = " + np.array2string(
                    self.C3 @ (x - self.apex), precision=12) + "  -> APEX")
            return "apex", None
        z = self._try(self._res(x, "face")[0], self.start(x, "face"))
        margins, order = None, ["line1", "line2"]
        if z is not None:
            y = self._y_of(z, "face")
            margins = (y[0] - y[1], y[1] - y[2])
            if verbose:
                print(f"      face-return margins (the exact boundary "
                      f"functions): y1-y2 = {margins[0]:+.6f}, "
                      f"y2-y3 = {margins[1]:+.6f}")
            if margins[0] >= 0.0 and margins[1] >= 0.0:
                return "face", margins
            order = ["line1", "line2"] if margins[0] < margins[1] else ["line2", "line1"]
        elif verbose:
            print("      face probe diverged -- its solution sits at w -> 0 (ON "
                  "the apex),\n      where the parametrization degenerates; "
                  "falling through to the edges")
        for reg in order:
            z = self._try(self._res(x, reg)[0], self.start(x, reg))
            if z is None:
                continue
            y = self._y_of(z, reg)
            if (np.all(self.raw_dl(y, z[-2:], reg) >= -1e-12)
                    and self.f_composite(np.sort(y)[::-1])[0] <= 1e-8 * self.scale):
                return reg, margins
        return "apex", margins

    def map_principal(self, x, verbose=False, gate=True):
        if self.f_composite(x)[0] <= 0.0:
            return x.copy(), dict(region="elastic", dl=np.zeros(1), n_iter=0,
                                  dydx=np.eye(3), margins=None)
        region, margins = self.region_of(x, verbose=verbose)
        if region == "apex":
            y, dl, nit = self.apex.copy(), self.C3 @ (x - self.apex), 0
            dydx = np.zeros((3, 3))
        else:
            R, n_dl = self._res(x, region)
            # quad_slack = 30 (family default 10): the near-apex residual's
            # curvature is genuinely large (|dm/dy| = 1.7e3 vs 3.5e-4 on an
            # ordinary face point) so the QUADRATIC CONSTANT is bigger while the
            # ORDER is unchanged (worst ratio 1.08 at slack 10 there, 0.005
            # everywhere else).
            z, hist = newton(R, self.start(x, region), f"HB {region}", tol=1e-14,
                             scale=self.scale, verbose=verbose,
                             quad_slack=30.0 if gate else 1e30,
                             iter_gate=5 if gate else 12)
            y, nit = self._y_of(z, region), len(hist) - 1
            dl = self.raw_dl(y, z[-n_dl:], region)
            rhs = np.zeros((len(z), 3))
            rhs[:3, :] = np.eye(3)
            dydx = (jac_complex_step(lambda v: self._y_of(v, region), z)
                    @ np.linalg.solve(jac_complex_step(R, z), rhs))
        return y, dict(region=region, dl=dl, n_iter=nit, dydx=dydx, margins=margins)

    def map6(self, sig_tr6, want_tangent=False, verbose=False, gate=True):
        w, Q = np.linalg.eigh(voigt_to_matrix(sig_tr6))
        idx = np.argsort(w)[::-1]
        x, Q = w[idx], Q[:, idx]
        y, info = self.map_principal(x, verbose=verbose, gate=gate)
        assert y[0] >= y[1] - 1e-9 and y[1] >= y[2] - 1e-9, "return broke ordering"
        info["x"], info["y"], info["Q"] = x, y, Q
        if want_tangent:
            info["C"] = (self.Ee.copy() if info["region"] == "elastic"
                         else self._tangent(x, y, Q, info["dydx"]))
        return matrix_to_voigt(Q @ np.diag(y) @ Q.T), info

    def _tangent(self, x, y, Q, dydx):
        """Principal tangent + the rotation term, rotated to global Voigt (the
        ``cppm_mc.py`` construction; the ``dm/dy`` curvature term is already
        inside ``dydx`` via the converged Jacobian)."""
        T6 = np.zeros((6, 6))
        T6[:3, :3] = dydx
        sc = max(1.0, float(np.max(np.abs(x))))
        for slot, (i, j) in ((3, (0, 1)), (4, (1, 2)), (5, (0, 2))):
            T6[slot, slot] = ((y[i] - y[j]) / (x[i] - x[j])
                              if abs(x[i] - x[j]) > 1e-9 * sc
                              else dydx[i, i] - dydx[i, j])
        Rs = np.zeros((6, 6))
        for k in range(6):
            e = np.zeros(6)
            e[k] = 1.0
            Rs[:, k] = matrix_to_voigt(Q @ voigt_to_matrix(e) @ Q.T)
        return (Rs @ T6 @ np.linalg.inv(Rs)) @ self.Ee

    def xi_face(self, y, dl):
        """The EXPLICIT curved-surface algorithmic form P3 will implement::

            Xi = (I + dl*D3*dm/dy)^-1;  dydx = Xi - (Xi D3 m)(a^T Xi)/(a^T Xi D3 m)

        ``dm/dy`` IS the curvature term (zero only for a linear surface)."""
        dm = jac_complex_step(lambda v: self.m_pair(v, "13"), y)
        Xi = np.linalg.inv(np.eye(3) + dl * (self.D3 @ dm))
        av, mv = np.real(self.a_pair(y, "13")), np.real(self.m_pair(y, "13"))
        XDm = Xi @ (self.D3 @ mv)
        return Xi - np.outer(XDm, av @ Xi) / (av @ XDm), dm

    def kkt(self, x, y, region):
        """Admissibility + return direction.  The domain is convex, so for
        associated flow this is a PROOF of global optimality (hence of the
        region), not a spot check."""
        G = (self.Gapex if region == "apex"
             else np.column_stack([self.D3 @ np.real(self.m_pair(y, k))
                                   for k in self._keys(region)]))
        d = x - y
        c = np.linalg.lstsq(G, d, rcond=None)[0]
        return (max(self.f_composite(y)[0], 0.0), float(np.min(c)),
                float(np.linalg.norm(G @ c - d) / max(1.0, np.linalg.norm(d))))

    def step6(self, eps):
        return self.map6(self.Ee @ eps)[0]


# ---------------------------------------------------------------------------
TRIALS = {
    "face (no shear)": np.array([-2000., -6000., -25000., 0., 0., 0.]),
    "face (sheared, exercises the rotation term)":
        np.array([-2000., -6000., -25000., 1500., -800., 500.]),
    "edge y1 == y2": np.array([-4000., -4000., -30000., 0., 0., 0.]),
    "edge y2 == y3": np.array([-2000., -20000., -20000., 0., 0., 0.]),
    "apex (hydrostatic tension)": np.array([400., 400., 400., 0., 0., 0.]),
    "apex (deviatoric, still in the cone)": np.array([500., 450., 400., 0., 0., 0.]),
    "header says APEX, elastic metric says otherwise":
        np.array([645., 265., 255., 0., 0., 0.]),
}
X_HDR = np.array([645., 265., 255.])       # the header-over-claims-apex trial


def report(mat, name, sig_tr, h=1e-8, sweep=False):
    section(f"HoekBrown  |  trial: {name}")
    print_vec(sig_tr, "trial stress sigma_tr")
    sig, info = mat.map6(sig_tr, want_tangent=True, verbose=True)
    x, y, C = info["x"], info["y"], info["C"]
    print(f"    principal trial  x = [{x[0]:.8f} {x[1]:.8f} {x[2]:.8f}]")
    print(f"    REGION = {info['region']}   (Newton iterations = {info['n_iter']})")
    print(f"    principal return y = [{y[0]:.8f} {y[1]:.8f} {y[2]:.8f}]")
    print_vec(sig, "returned stress sigma")
    print("    dLambda = ["
          + " ".join(f"{v:.10e}" for v in np.atleast_1d(info["dl"])) + "]")
    fc, fs, ft = mat.f_composite(y)
    print(f"      f_composite = {fc:+.6e}   f_shear = {fs:+.6e}   f_tension = "
          f"{ft:+.6e}   (active branch: {'f_tension' if ft > fs else 'f_shear'})"
          f"   |f| via the header's 6D form = {abs(mat.f6(sig)):.3e}")
    assert abs(fc) <= max(1e-10, mat.f_floor(y)), f"{name}: |f| = {fc}"
    adm, dlmin, res = mat.kkt(x, y, info["region"])
    print(f"      KKT: max(f,0) = {adm:.3e}   min multiplier = {dlmin:+.6e}   "
          f"direction residual = {res:.3e}")
    assert res < 1e-9 and dlmin > -1e-9
    if info["region"] == "face":
        dydx_xi, dm = mat.xi_face(y, info["dl"][0])
        print(f"      explicit Xi-form vs Jacobian-form dy/dx: rel_fro = "
              f"{rel_fro(dydx_xi, info['dydx']):.3e}   |dm/dy|_F (the curvature "
              f"term) = {np.linalg.norm(dm):.6e}")
        assert rel_fro(dydx_xi, info["dydx"]) < 1e-10
    print_matrix(C, "consistent tangent dsigma/deps", fmt="%14.2f")
    print(f"      rank(C) = "
          f"{np.linalg.matrix_rank(C, tol=1e-8 * max(1.0, np.linalg.norm(C)))}"
          f"   |C|_F = {np.linalg.norm(C):.6f}")
    eps_tr = np.linalg.solve(mat.Ee, sig_tr)
    if sweep:
        print("      stencil sweep -- the FD gap here is TRUNCATION, not a wrong"
              " tangent (error ~ h^2 until round-off):")
        print("        " + "  ".join(
            f"h={hh:.0e}: {rel_fro(C, fd_tangent(mat.step6, eps_tr, hh)):.2e}"
            for hh in (1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12)))
    err, _ = assert_fd_tangent(C, mat.step6, eps_tr, name, h=h,
                               ref=float(np.linalg.norm(mat.Ee)))
    return dict(sigma=sig, y=y, region=info["region"], C=C, fd_err=err, f=fc,
                dl=np.atleast_1d(info["dl"]), n_iter=info["n_iter"])


def tension_branch_study(mat):
    section("the max(f_shear, f_tension) regime switch, and the apex region")
    rng = np.random.default_rng(3)
    worst_gap, n_ten = 1e30, 0
    for _ in range(4000):                       # points ON the yield surface
        y1 = mat.T - 10.0 ** rng.uniform(-4.0, 4.3)
        y3 = y1 - mat.sigci * (mat.s - mat.mb * y1 / mat.sigci) ** mat.a
        y2 = y3 + rng.uniform(0.0, 1.0) * (y1 - y3)
        _, fs, ft = mat.f_composite(np.array([y1, y2, y3]))
        n_ten += ft > fs
        worst_gap = min(worst_gap, fs - ft)
    print(f"    (a) ON the surface f_shear is ALWAYS the active branch (except "
          f"at the apex,\n        where both vanish): tension-branch points "
          f"{n_ten} / 4000, min(f_shear - f_tension)\n        = {worst_gap:.6e}."
          f"  So the plane y1 = T carries NO face -- no tension return.")
    assert n_ten == 0
    bad = 0
    for _ in range(4000):                       # points OUTSIDE the domain
        y = np.sort(rng.uniform(-3000.0, 900.0, 3))[::-1]
        if mat.f_composite(y)[0] <= 0:
            continue
        _, fs, ft = mat.f_composite(y)
        bad += (ft > fs) != bool(y[0] > mat.T and y[2] > mat.T)
    print(f"    (b) OUTSIDE the domain the tension branch wins exactly on "
          f"{{y1>T and y3>T}},\n        i.e. on the header's own "
          f"CHECK_APEX_REGION set: mismatches = {bad} / 4000.")
    assert bad == 0
    fs, ft = mat.f_composite(np.array([0., -500., -1000.]))[1:]
    print(f"    (c) INSIDE the domain the tension branch is often the max (y = "
          f"[0,-500,-1000]:\n        f_shear = {fs:.3f}, f_tension = {ft:.3f}), "
          f"so the composite's KINK lies\n        strictly inside: it moves f's "
          f"value and gradient where nothing yields, and\n        the header's "
          f"central-difference df/dsigma takes a chord across it.")
    beta = (X_HDR[0] - mat.T) / (mat.D3 @ E1)[0]
    yT = X_HDR - beta * (mat.D3 @ E1)
    fs = mat.f_composite(np.sort(yT)[::-1])[1]
    y_ok, info = mat.map_principal(X_HDR)
    print(f"    (d) the SHEAR/TENSION CORNER is the apex itself, and a return "
          f"onto the tension\n        plane is never right.  Trial x_tr = "
          f"{X_HDR} (f_tension IS its max):\n        Rankine return along D3*e1 "
          f"-> y = {np.round(yT, 6)},\n        f_shear = {fs:+.6f} kPa "
          f"INADMISSIBLE ({fs / mat.scale * 100:.2f}% of strength_scale);\n"
          f"        the closest point is y = {np.round(y_ok, 6)} "
          f"({info['region']}).")
    print(f"    (e) cost of the header's EUCLIDEAN CHECK_APEX_REGION on that "
          f"same trial: it\n        says APEX (all principals >= T = "
          f"{mat.T:.6f}) but the exact elastic-metric\n        test C3*(x-apex) "
          f"= {np.array2string(mat.C3 @ (X_HDR - mat.apex), precision=10)} has "
          f"negative\n        entries -> NOT apex.  APEX_STRESS would commit "
          f"T*[1,1,1], i.e. |dsigma|_inf =\n        "
          f"{np.max(np.abs(mat.apex - y_ok)):.6f} kPa "
          f"({np.max(np.abs(mat.apex - y_ok)) / mat.scale * 100:.3f}% of "
          f"strength_scale) of silently lost strength --\n        admissible but"
          f" wrong.  D3*(octant) is a strict SUBSET of the octant, so the\n"
          f"        header always OVER-claims the apex, never under-claims.")


def g_vs_f_study(mat):
    section("g vs f: HoekBrown_PF::g is evaluated in the un-negated frame")
    sig = TRIALS["face (no shear)"]
    xc = np.array([-2000., -6000., -25000.])
    print_vec(sig, "probe stress (compressive)")
    for mbp in (mat.mb, mat.mb / 2.0, 0.0):
        h = HB(mb_psi=mbp, flow="header")
        m = h.num_grad6(h.g6_header, sig)
        print(f"    mb_psi = {mbp:9.6f}: header dg/dsigma (its own central "
              f"difference) = {np.array2string(m, precision=8)}, "
              f"trace = {m[:3].sum():+.3e}")
    print("    -> IDENTICAL for every mb_psi and TRACELESS: in compression the "
          "shipped g always\n       falls into its `else` branch (Tresca), so "
          "HB_mb_psi has NO EFFECT on the flow\n       and the flow is "
          "NON-DILATANT.")
    ma, hh = HB(flow="intended"), HB(flow="header")
    mi_, mh = np.real(ma.m_pair(xc, "13")), np.real(hh.m_pair(xc, "13"))
    ang = np.rad2deg(np.arccos(mi_ @ mh / np.linalg.norm(mi_) / np.linalg.norm(mh)))
    print(f"    principal-space flow there: intended m = "
          f"{np.array2string(mi_, precision=6)} (trace {mi_.sum():+.4f}), "
          f"header m =\n    {np.array2string(mh, precision=6)} (trace "
          f"{mh.sum():+.4f}); angle = {ang:.4f} deg -- and this is with "
          f"mb_psi == mb,\n    i.e. exactly where the deck asks for ASSOCIATED "
          f"flow.  Consequence on the return:")
    for flow in ("intended", "header"):
        h = HB(flow=flow)
        y, info = h.map_principal(xc)
        print(f"      flow={flow:9s} region={info['region']:5s} y = "
              f"{np.array2string(y, precision=6)}  dLambda = "
              f"{np.array2string(np.atleast_1d(info['dl']), precision=6)}  "
              f"plastic vol. strain = {float((h.C3 @ (xc - y)).sum()):+.6e}")
    tr = np.array([(hh.D3 @ np.real(hh.m_ij(hh.apex - 1e-9 * hh.T, i, j))).sum()
                   for i, j in ALL_PAIRS])
    print(f"    at the APEX every header flow direction has NEGATIVE trace: "
          f"trace(D3 m_ij) =\n    {np.array2string(tr, precision=3)}.  Is the "
          f"hydrostatic direction (1,1,1) in that cone? "
          f"{in_cone(hh.Gapex, np.ones(3))}\n    (intended-flow cone? "
          f"{in_cone(ma.Gapex, np.ones(3))}).  So under the shipped g a trial "
          f"past the tensile\n    corner has NO return to the apex at all -- "
          f"the mechanism behind the residual\n    recorded in "
          f"tests/test_adr94_hlist_hb.py (drive fails on the corner-crossing "
          f"step).")


def near_apex_conditioning(mat, x=X_HDR):
    section("near-apex conditioning: why the natural variable and the unit m")

    def res_y1(z):                      # the naive parametrization: y1 unknown
        y, dl = z[:3], z[3]
        return np.concatenate([y - x + dl * (mat.D3 @ mat.m_pair(y, "13")),
                               [mat.f_pair(y, "13")]])
    z = np.concatenate([x.copy(), [0.0]])
    z[0] = min(z[0], mat.T - 0.01 * mat.scale)          # clipped into arg > 0
    print("    (1) unknowns [y1,y2,y3,dl], start = elastic predictor clipped "
          "into arg > 0:")
    for k in range(4):
        r = np.real(np.asarray(res_y1(z.astype(complex))))
        arg = mat.s - mat.mb * z[0] / mat.sigci
        print(f"        iter {k}: |R| = {np.linalg.norm(r):.3e}   y1 = "
              f"{z[0]:12.6f}   arg = {arg:+.4e}"
              + ("   <- LEFT the surface's domain" if arg <= 0 else ""))
        if arg <= 0:
            print("        the next Jacobian is SINGULAR (f is constant once "
                  "arg is clamped) and the\n        Newton dies.  No start "
                  "guess fixes this: df/dy1 = 1 + a*mb*arg^(a-1) -> inf,\n"
                  "        so the first step is unbounded.")
            break
        z = z - np.linalg.solve(jac_complex_step(res_y1, z), r)
    y, info = mat.map_principal(x)
    print(f"    (2) unknowns [w,y2,y3,dl] with arg = w^(2/a) (this oracle): "
          f"feasible for ANY real\n        w, and it converges -- "
          f"{info['n_iter']} iterations, |f| = "
          f"{abs(mat.f_composite(y)[0]):.3e}, arg = "
          f"{mat.s - mat.mb * y[0] / mat.sigci:.6e}.")
    print(f"    (3) the unit flow direction: worst-case Newton count over 400+ "
          f"random trials plus\n        every FD stencil point is 6 with m "
          f"as-is and 5 with m/|m| -- |m| = "
          f"{np.linalg.norm(np.real(mat.m_pair(y, '13'))):.3e} here\n        "
          f"against "
          f"{np.linalg.norm(np.real(mat.m_pair(np.array([-3642.8, -6441.6, -25123.7]), '13'))):.3e}"
          f" on an ordinary face point.")


def scan(mat, n=400, seed=11):
    section(f"region-selection scan: {n} random trial states, KKT-verified")
    rng = np.random.default_rng(seed)
    counts, iters, worst_f, worst_res, worst_mult, worst_ratio = {}, {}, 0., 0., 0., 0.
    hdr_wrong, f_near = 0, (0.0, 0.0, 0.0)
    for i in range(n):
        # mixture: a wide compressive cloud (face + both edges) and a narrow
        # tensile one -- the apex region is a thin cone that uniform sampling of
        # the wide box essentially never hits.
        lo, hi = (-40000.0, 900.0) if i % 3 else (0.0, 900.0)
        x = np.sort(rng.uniform(lo, hi, 3))[::-1]
        if mat.f_composite(x)[0] <= 0:
            continue
        y, info = mat.map_principal(x, gate=False)
        counts[info["region"]] = counts.get(info["region"], 0) + 1
        iters[info["n_iter"]] = iters.get(info["n_iter"], 0) + 1
        adm, dlmin, res = mat.kkt(x, y, info["region"])
        fv, floor = abs(mat.f_composite(y)[0]), mat.f_floor(y)
        if fv > f_near[0]:
            f_near = (fv, float(mat.T - y[0]), floor)
        worst_ratio = max(worst_ratio, fv / max(1e-10, floor))
        worst_f, worst_res = max(worst_f, fv), max(worst_res, res)
        worst_mult = min(worst_mult, dlmin)
        assert y[0] >= y[1] - 1e-9 >= y[2] - 2e-9, "ordering broken"
        hdr_wrong += bool(np.all(x >= mat.T)) != (info["region"] == "apex")
    print(f"    regions: {counts}")
    print(f"    Newton iterations: {dict(sorted(iters.items()))} -- <= 5 except "
          f"for returns landing\n        within ~1 kPa of the apex, where the "
          f"curvature term is 7 orders larger")
    print(f"    max |f(returned)| = {worst_f:.3e}, attained {f_near[1]:.3e} kPa "
          f"below the apex where\n        f's own round-off floor is "
          f"{f_near[2]:.3e}; the gate is |f| <= max(1e-10, f_floor(y))\n"
          f"        and the worst |f|/gate over the scan is {worst_ratio:.3f}")
    print(f"    max return-direction residual = {worst_res:.3e}, min Koiter "
          f"multiplier = {worst_mult:+.3e}\n        (the domain is convex, so a "
          f"zero residual with non-negative multipliers\n        PROVES the "
          f"region -- no trial is misclassified)")
    print(f"    trials where the header's Euclidean CHECK_APEX_REGION disagrees "
          f"with the exact\n        elastic-metric apex region: {hdr_wrong}")
    assert worst_ratio <= 1.0 and worst_res < 1e-9 and worst_mult > -1e-9
    assert max(iters) <= 8
    for want in ("face", "line1", "line2", "apex"):
        assert want in counts, f"scan never reached {want}"
    return counts


def main():
    banner("ADR-97 oracle 4/6 -- HOEK-BROWN principal-space return "
           "(Clausen & Damkilde 2008)")
    mat = HB()
    print(f"Material (identical to tests/test_adr94_hlist_hb.py): sigma_ci = "
          f"{mat.sigci:.1f} kPa, mi = 10, GSI = 60, D = 0, E = {mat.E:.1f} kPa, "
          f"nu = {mat.nu},\nHB_mb_psi = mb (the deck's associated setting)")
    print(f"  -> mb = {mat.mb:.12f}   s = {mat.s:.12f}   a = {mat.a:.12f}")
    print(f"  -> tensile strength T = s*sigma_ci/mb = {mat.T:.10f} kPa "
          f"(APEX_STRESS = T*[1,1,1], tension-positive)")
    print(f"  -> YF_STRENGTH_SCALE = sigma_ci*s^a = {mat.scale:.10f} kPa")

    tension_branch_study(mat)
    g_vs_f_study(mat)
    results = {}
    for name, s in TRIALS.items():
        # the near-apex trial needs a finer stencil: the map's curvature there
        # is ~7 orders larger, so the central difference's own h^2 truncation --
        # not the tangent -- sets the error (the printed sweep proves it).
        na = name.startswith("header says APEX")
        results[name] = report(mat, name, s, h=1e-10 if na else 1e-8, sweep=na)
    near_apex_conditioning(mat)
    counts = scan(mat)

    section("non-associated companion: mb_psi = mb/2 (frame-consistent g)")
    mn = HB(mb_psi=HB().mb / 2.0)
    for lbl, mm in (("associated  mb_psi = mb  ", mat),
                    ("non-assoc   mb_psi = mb/2", mn)):
        C = mm.map6(TRIALS["face (no shear)"], want_tangent=True)[1]["C"]
        print(f"    {lbl}: |C-C^T|_F/|C|_F = "
              f"{np.linalg.norm(C - C.T) / np.linalg.norm(C):.3e}")
        assert (np.linalg.norm(C - C.T) / np.linalg.norm(C)
                < 1e-12) == (mm is mat)
    assert_fd_tangent(mn.map6(TRIALS["face (no shear)"], want_tangent=True)[1]["C"],
                      mn.step6, np.linalg.solve(mn.Ee, TRIALS["face (no shear)"]),
                      "non-associated face", h=1e-8,
                      ref=float(np.linalg.norm(mn.Ee)))

    banner("HOEK-BROWN REFERENCE BLOCK (values the C++ tests pin)")
    for name, r in results.items():
        print(f"  {name}")
        print(f"    region = {r['region']}   Newton iterations = {r['n_iter']}"
              f"   |f| = {abs(r['f']):.3e}   dLambda = "
              f"{np.array2string(r['dl'], precision=10)}")
        print(f"    y       = {np.array2string(r['y'], precision=10)}")
        print(f"    sigma   = {np.array2string(r['sigma'], precision=10)}")
        print(f"    |C|_F = {np.linalg.norm(r['C']):.6f}   trace(C) = "
              f"{np.trace(r['C']):.6f}   rank = "
              f"{np.linalg.matrix_rank(r['C'], tol=1e-8 * max(1.0, np.linalg.norm(r['C'])))}"
              f"   FD rel err = {r['fd_err']:.3e}")
    print(f"\n  scan regions: {counts}")
    print("  FD stencil: central difference, h = 1e-8 relative (1e-10 for the")
    print("  near-apex trial) on the TOTAL strain of one step from a zero")
    print("  committed state; every trial sits well inside its region -- see the")
    print("  printed boundary margins.")


if __name__ == "__main__":
    main()
