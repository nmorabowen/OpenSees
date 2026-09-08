"""ADR-97 P0 -- shared conventions and Newton/tangent machinery for the
closest-point (CPPM) oracles.

CONVENTIONS (all taken from the shipped headers, not from a textbook)
--------------------------------------------------------------------
* **Tension positive.**  ``VoigtVector::meanStress() == trace()/3``
  (``OTHER/eigenAPI/typedefs.h:258``), so ``p`` is the TENSION-positive mean
  stress everywhere.  The Drucker-Prager / Mohr-Coulomb apex therefore sits at
  positive ``p``.
* **Voigt order** ``[11, 22, 33, 12, 23, 13]`` -- shear stored ONCE.  Strain-like
  vectors store ENGINEERING shear ``gamma_ij = 2 eps_ij``; stress-like vectors
  store the tensor component.
* **Two contractions**, verbatim from ``typedefs.h:418`` / ``:461``::

      tensor_dot_stress_like(a,b)             = a.b with shear weight 2
      tensor_dot_engineering_strain_like(a,b) = a.b with shear weight 1/2

* **Every YF/PF derivative is a VOIGT derivative** (ADR-94 wp/94c B5): the
  gradient w.r.t. the STORED slot, so the three shear slots carry the extra
  factor 2 relative to the tensor derivative.  This is what makes
  ``TrialPlastic_Strain += dLambda*m`` and ``Eelastic*m`` plain Voigt ops.
* **Elasticity** ``LinearIsotropic3D_EL``
  (``ElasticityModels/LinearIsotropic3D_EL.h:54-61``): ``lambda = nu E/((1+nu)
  (1-2nu))``, ``mu = E/(2(1+nu))``, ``EE(0,0..2,2) = 2mu+lambda``, off-diagonal
  normal block ``lambda``, ``EE(3,3) = EE(4,4) = EE(5,5) = mu`` -- i.e. the
  ENGINEERING-shear form (G on the shear diagonal, not 2G).
* ``SQRT_2_over_3`` is a TRUNCATED literal in
  ``ASDPlasticMaterial3DGlobals.h:41`` (0.816496580928, 3e-13 relative off
  ``sqrt(2/3)``).  Mirrored here so 1e-10 comparisons stay honest.

NEWTON / TANGENT MACHINERY
--------------------------
The CPPM residuals are solved by a full Newton whose Jacobian comes from a
COMPLEX-STEP derivative (``Im f(x + i h e_j)/h`` with ``h = 1e-200``).  That is
exact to machine precision -- no subtractive cancellation -- so the Newton is a
true Newton (quadratic) and the consistent tangent obtained from the converged
Jacobian is the analytic one.  ``check_complex_step_jacobian`` cross-validates
it against a central difference once per family so the choice is warranted, not
asserted.

The consistent tangent is then ``d sigma / d eps`` from

    J @ [d sigma; d q; d dlambda] = [E @ d eps; 0; 0]

(the only eps-dependence of the residual is the elastic predictor
``sigma_tr = sigma_n + E (eps - eps_n)``, so ``dR/deps = [-E; 0; 0]``).
"""
import numpy as np

# ---------------------------------------------------------------------------
# constants and Voigt algebra (verbatim from the headers)
# ---------------------------------------------------------------------------
SQRT_2_over_3 = 0.816496580928          # ASDPlasticMaterial3DGlobals.h:41
W_STRESS = np.array([1., 1., 1., 2., 2., 2.])       # typedefs.h:418
W_STRAIN = np.array([1., 1., 1., .5, .5, .5])       # typedefs.h:461
VOIGT_LABELS = ("11", "22", "33", "12", "23", "13")


def elastic_voigt(E, nu):
    """LinearIsotropic3D_EL, engineering-shear Voigt storage."""
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    mu = E / (2.0 * (1.0 + nu))
    D = np.zeros((6, 6))
    D[:3, :3] = lam
    D[0, 0] = D[1, 1] = D[2, 2] = lam + 2.0 * mu
    D[3, 3] = D[4, 4] = D[5, 5] = mu
    return D


def dev(v):
    """VoigtVector::deviator() -- subtract trace/3 from the three normal slots."""
    p = (v[0] + v[1] + v[2]) / 3.0
    out = v.copy()
    out[0] = v[0] - p
    out[1] = v[1] - p
    out[2] = v[2] - p
    return out


def mean_stress(v):
    """VoigtVector::meanStress() -- trace/3, TENSION POSITIVE."""
    return (v[0] + v[1] + v[2]) / 3.0


def dot_stress(a, b):
    """tensor_dot_stress_like (typedefs.h:418)."""
    return (a * W_STRESS * b).sum()


def dot_strain(a, b):
    """tensor_dot_engineering_strain_like (typedefs.h:461)."""
    return (a * W_STRAIN * b).sum()


def getJ2(v):
    d = dev(v)
    return 0.5 * dot_stress(d, d)


P_DEV = np.eye(6)
P_DEV[:3, :3] -= 1.0 / 3.0          # deviatoric projector on stress-like Voigt


def voigt_to_matrix(v):
    return np.array([[v[0], v[3], v[5]],
                     [v[3], v[1], v[4]],
                     [v[5], v[4], v[2]]], dtype=float)


def matrix_to_voigt(M):
    return np.array([M[0, 0], M[1, 1], M[2, 2], M[0, 1], M[1, 2], M[0, 2]])


# ---------------------------------------------------------------------------
# complex-step Jacobian + Newton
# ---------------------------------------------------------------------------
CS_H = 1e-200


def jac_complex_step(R, x):
    """dR/dx by complex step.  R must be written with complex-safe numpy ops."""
    n = len(x)
    r0 = np.asarray(R(x.astype(complex)))
    J = np.zeros((len(r0), n))
    for j in range(n):
        xp = x.astype(complex)
        xp[j] = xp[j] + 1j * CS_H
        J[:, j] = np.imag(np.asarray(R(xp))) / CS_H
    return J


def check_complex_step_jacobian(R, x, h=1e-6, label=""):
    """Cross-validate the complex-step Jacobian against a central difference."""
    Jc = jac_complex_step(R, x)
    n = len(x)
    Jf = np.zeros_like(Jc)
    for j in range(n):
        s = h * max(1.0, abs(x[j]))
        xp = x.astype(complex).copy(); xp[j] += s
        xm = x.astype(complex).copy(); xm[j] -= s
        Jf[:, j] = (np.real(np.asarray(R(xp))) - np.real(np.asarray(R(xm)))) / (2 * s)
    err = rel_fro(Jc, Jf)
    print(f"    complex-step Jacobian vs central difference{label}: rel_fro = {err:.3e}")
    assert err < 1e-6, f"complex-step Jacobian disagrees with FD: {err}"
    return Jc


def newton(R, x0, name, tol=1e-12, max_iter=12, scale=1.0, verbose=True,
           order_floor=1e-13, quad_slack=10.0, iter_gate=5):
    """Full Newton with a complex-step Jacobian.

    Prints the residual history and asserts

    * convergence in <= ``iter_gate`` iterations (5 for the production start
      guess; relaxed only where a deliberately BAD start guess is being used to
      exhibit the quadratic descent), and
    * SECOND-ORDER convergence, as ``rho_{k+1} <= quad_slack * rho_k**2`` on the
      dimensionless residual ``rho = |R| / scale``, for every transition whose
      TARGET is still above the round-off floor ``order_floor``.

    Why not a three-point order estimate: the last residual of a converged
    Newton sits at the round-off floor of the residual evaluation itself, and
    the first is pre-asymptotic, so a log-ratio over any three consecutive
    iterates of a 3-4 iteration history is dominated by whichever end is
    contaminated (measured: 1.35 to 2.05 for histories that are visibly
    quadratic).  The squaring test above is the same statement without that
    fragility.  The three-point estimate is still PRINTED, for information.
    """
    x = np.array(x0, dtype=float)
    hist = []
    for it in range(max_iter):
        r = np.real(np.asarray(R(x.astype(complex))))
        nr = float(np.linalg.norm(r))
        hist.append(nr)
        if nr <= tol * scale:
            break
        J = jac_complex_step(R, x)
        x = x - np.linalg.solve(J, r)
    else:
        r = np.real(np.asarray(R(x.astype(complex))))
        hist.append(float(np.linalg.norm(r)))
    n_iter = len(hist) - 1
    if verbose:
        print(f"    Newton [{name}] residual history:",
              " ".join(f"{v:.3e}" for v in hist))
    assert hist[-1] <= tol * scale, f"{name}: Newton did not converge: {hist}"
    assert n_iter <= iter_gate, (f"{name}: {n_iter} iterations "
                                 f"(> {iter_gate}): {hist}")
    rho = [v / scale for v in hist]
    tested = []
    for k in range(len(rho) - 1):
        if rho[k + 1] <= order_floor:
            continue                       # target is at the round-off floor
        tested.append((k, rho[k + 1], quad_slack * rho[k] ** 2))
    if tested:
        worst = max(a / b for _, a, b in tested)
        orders = []
        for k in range(len(rho) - 2):
            if rho[k + 2] > order_floor and rho[k + 1] > order_floor:
                orders.append(np.log(rho[k + 2] / rho[k + 1])
                              / np.log(rho[k + 1] / rho[k]))
        if verbose:
            best = f"{max(orders):.2f}" if orders else "n/a"
            print(f"    Newton [{name}] 2nd-order test: max "
                  f"rho_k+1/({quad_slack:g} rho_k^2) = {worst:.3f} over "
                  f"{len(tested)} transition(s) ({n_iter} iterations); "
                  f"best 3-point order estimate = {best}")
        assert worst <= 1.0, (f"{name}: not second order -- "
                              f"rho_k+1 > {quad_slack} rho_k^2 ({hist})")
    elif verbose:
        print(f"    Newton [{name}] converged in {n_iter} iteration(s); every "
              f"target residual is already at the {order_floor:g} round-off "
              f"floor, so no order test is formed")
    return x, hist


# ---------------------------------------------------------------------------
# tangents
# ---------------------------------------------------------------------------
def tangent_from_jacobian(J, E):
    """Solve J @ [dsigma; dq; ddl] = [E; 0; 0] and return the 6x6 dsigma/deps."""
    n = J.shape[0]
    rhs = np.zeros((n, 6))
    rhs[:6, :] = E
    Z = np.linalg.solve(J, rhs)
    return Z[:6, :]


def tangent_general(J, dR_deps):
    """Consistent tangent when dR/deps is not the standard ``[-E; 0; 0]``.

    Solves ``J @ dx/deps = -dR/deps`` and returns the first six rows.
    """
    Z = np.linalg.solve(J, -dR_deps)
    return Z[:6, :]


def rel_fro(A, B):
    """Frobenius relative error of A against B."""
    den = np.linalg.norm(B)
    return float(np.linalg.norm(A - B) / (den if den > 0 else 1.0))


def fd_tangent(map_fn, eps0, h=1e-8):
    """Central-difference d sigma / d eps of the oracle's OWN map.

    ``map_fn(eps) -> sigma`` must re-run the return map from the SAME committed
    state for every call (no state carried between calls).
    """
    C = np.zeros((6, 6))
    for j in range(6):
        s = h * max(1.0, abs(eps0[j]))
        ep = eps0.copy(); ep[j] += s
        em = eps0.copy(); em[j] -= s
        C[:, j] = (map_fn(ep) - map_fn(em)) / (2 * s)
    return C


def assert_fd_tangent(C, map_fn, eps0, label, h=1e-8, tol=1e-6, ref=None):
    """Central-difference check of a consistent tangent.

    ``ref`` (typically ``|E_elastic|_F``) lets a RANK-DEFICIENT tangent pass on
    the absolute criterion: when the exact tangent is the zero matrix, the
    finite difference still returns the round-off of the map itself
    (~1e-16 * |sigma| / h), which has no meaningful relative measure.  The gate
    is then ``|C - C_cd|_F <= tol * ref``.
    """
    Cfd = fd_tangent(map_fn, eps0, h)
    err = rel_fro(C, Cfd)
    abs_err = float(np.linalg.norm(C - Cfd))
    ok = err <= tol
    if not ok and ref is not None:
        ok = abs_err <= tol * ref
    extra = "" if ref is None else f"   |C - C_cd|_F/|E|_F = {abs_err / ref:.3e}"
    print(f"    FD check [{label}]: rel_fro(C_consistent, C_cd) = {err:.3e}"
          f"   (stencil h = {h:g} relative){extra}")
    assert ok, f"{label}: consistent tangent FD error {err} (abs {abs_err}) > {tol}"
    return err, Cfd


def print_matrix(C, label, fmt="%12.4f"):
    print(f"    {label} (6x6):")
    for i in range(6):
        print("      " + " ".join(fmt % C[i, j] for j in range(6)))
    asym = np.linalg.norm(C - C.T) / max(np.linalg.norm(C), 1e-30)
    print(f"      |C - C^T|_F / |C|_F = {asym:.3e}")


def print_vec(v, label, fmt="%15.8e"):
    print(f"    {label} = [" + " ".join(fmt % x for x in np.asarray(v).ravel()) + "]")


def banner(title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def section(title):
    print()
    print("-" * 78)
    print(title)
    print("-" * 78)
