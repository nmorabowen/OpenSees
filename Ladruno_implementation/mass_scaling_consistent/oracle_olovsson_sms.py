"""Consistent / Olovsson selective mass-scaling ORACLE (Ladruno ADR 38, Phase 0).

De-risks the v2 ("consistent") increment of the mass-scaling stack before any C++.
Companion to the lumped v1 already shipped (CentralDifferenceSMS, ADR 36). This is a
self-contained numpy reference -- no OpenSees -- of the centroidal Olovsson scaling
mass, and it proves the three claims the C++ must reproduce:

  C1  FREQUENCY PRESERVATION. At the SAME per-element timestep target, conventional
      (lumped, diagonal) scaling shifts the global fundamental f1; Olovsson (consistent)
      scaling leaves f1 essentially unchanged -- because its scaling mass annihilates
      rigid-body translation (M_bar . 1 = 0) and loads only the relative/deformation
      modes. This is the entire v2 selling point.

  C2  STABLE-STEP GAIN with LESS added mass. Both schemes raise the governing element's
      stable step dt_e -> dtTarget with the SAME per-element factor beta = (dtTarget/dt_e)^2-1,
      but Olovsson injects less TOTAL mass (trace M_bar_e = beta*M_e*(1-sum m_a^2/M_e^2)
      < beta*M_e) AND places it where it does not move f1.

  C3  CHEAP SOLVE. M_tilde = M_lump + sum_e M_bar_e is SPD but non-diagonal, so the
      explicit leap-frog needs a = M_tilde^{-1} r each step. Because M_bar is a small
      perturbation of the dominant lumped diagonal, preconditioned CG (lumped-diagonal
      preconditioner) converges in a handful of iterations -- the matrix-free solve the
      C++ integrator will run inside the step (keeps `system Diagonal` as the precond).

  --------------------------------------------------------------------------------
  Centroidal Olovsson scaling mass (Olovsson, Simonsson & Unosson, CNME 2005):

      M_bar_e = beta_e * [ diag(m_a) - m m^T / M_e ]            (per spatial direction)

  with m = (m_1..m_n) the element's lumped nodal masses, M_e = sum m_a. Row sums are
  zero => M_bar_e . t_rigid = 0 (rigid translation gets no added inertia). For a 2-node
  equal-mass bar this is beta*(m/4)*[[1,-1],[-1,1]]: the [1,1] (rigid) mode is untouched,
  the [1,-1] (relative) mode mass scales by (1+beta) -> omega_e -> omega_e/sqrt(1+beta)
  -> dt_e -> dt_e*sqrt(1+beta), exactly matching the lumped factor but WITHOUT shifting
  the low modes.

Run:  pythoncore-3.12-64\\python.exe Ladruno_implementation\\mass_scaling_consistent\\oracle_olovsson_sms.py
(numpy only). Exits non-zero if any oracle assertion fails.

Written: N. Mora-Bowen (Ladruno), 2026 -- Phase 0 oracle for the consistent-SMS C++.
"""
from __future__ import annotations

import numpy as np

# ---------------------------------------------------------------------------
# 1-D bar finite-element model: a chain of 2-node axial (truss) elements.
#   element e: nodes (e, e+1), length L_e, area A=1, density rho, modulus E.
#   consistent mass M_e = rho*A*L_e/6 * [[2,1],[1,2]]
#   lumped (row-sum) mass m_a = rho*A*L_e/2 each node
#   stiffness   K_e = E*A/L_e * [[1,-1],[-1,1]]
# Node 0 is fixed (cantilever bar) so the assembled system has no rigid mode and a
# well-defined f1 -- the quantity we care about preserving.
# ---------------------------------------------------------------------------


class Bar:
    def __init__(self, lengths, E=1.0e4, rho=1.0, A=1.0):
        self.L = np.asarray(lengths, float)
        self.E, self.rho, self.A = E, rho, A
        self.nnode = len(self.L) + 1
        self.nelem = len(self.L)

    def _ke(self, e):
        k = self.E * self.A / self.L[e]
        return k * np.array([[1.0, -1.0], [-1.0, 1.0]])

    def _me_consistent(self, e):
        m = self.rho * self.A * self.L[e] / 6.0
        return m * np.array([[2.0, 1.0], [1.0, 2.0]])

    def _me_lumped(self, e):
        m = self.rho * self.A * self.L[e] / 2.0
        return np.array([m, m])  # diagonal entries for the 2 nodes

    # ---- assembly (free DOFs = nodes 1..nnode-1; node 0 clamped) ----------
    def assemble_K(self):
        n = self.nnode
        K = np.zeros((n, n))
        for e in range(self.nelem):
            ke = self._ke(e)
            for a in range(2):
                for b in range(2):
                    K[e + a, e + b] += ke[a, b]
        return K[1:, 1:]  # drop clamped node 0

    def assemble_M_lumped(self):
        n = self.nnode
        m = np.zeros(n)
        for e in range(self.nelem):
            me = self._me_lumped(e)
            m[e] += me[0]
            m[e + 1] += me[1]
        return m[1:]  # diagonal vector, clamped node dropped

    def elem_dofs(self, e):
        """Global free-DOF indices of element e's nodes; -1 for the clamped node."""
        return [e - 1, e]  # node e -> dof e-1 (node 0 -> -1, clamped)


# ---------------------------------------------------------------------------
# per-element stable step from the SAME pencil the C++ kernel uses:
#   dt_e = 2 / sqrt(lambda_max(K_e v = lambda M_e^lump v))
# (M_e lumped diagonal; the largest generalized eigenvalue.)
# ---------------------------------------------------------------------------

def elem_dt(bar, e):
    ke = bar._ke(e)
    md = bar._me_lumped(e)
    # generalized eig with diagonal M: scale rows/cols by 1/sqrt(m)
    s = 1.0 / np.sqrt(md)
    A = (ke * s[:, None]) * s[None, :]
    w = np.linalg.eigvalsh(A)
    lam = w.max()
    return 2.0 / np.sqrt(lam) if lam > 0 else np.inf


# ---------------------------------------------------------------------------
# scaling-mass builders. beta_e = (dtTarget/dt_e)^2 - 1 for the sub-target elements.
# ---------------------------------------------------------------------------

def conventional_added(bar, dtTarget):
    """Lumped (diagonal) scaling: add (s-1)*m_a to each node of a sub-target element.
    Returns the diagonal added-mass vector over free DOFs."""
    nfree = bar.nnode - 1
    dm = np.zeros(nfree)
    added = 0.0
    for e in range(bar.nelem):
        dt_e = elem_dt(bar, e)
        if dt_e >= dtTarget:
            continue
        beta = (dtTarget / dt_e) ** 2 - 1.0
        md = bar._me_lumped(e)
        for a, dof in enumerate(bar.elem_dofs(e)):
            if dof >= 0:
                dm[dof] += beta * md[a]
                added += beta * md[a]
    return dm, added


def olovsson_added(bar, dtTarget):
    """Centroidal Olovsson scaling: M_bar_e = beta_e[diag(m_a) - m m^T/M_e].
    Returns the dense added-mass matrix over free DOFs and the injected trace."""
    nfree = bar.nnode - 1
    dM = np.zeros((nfree, nfree))
    added = 0.0
    for e in range(bar.nelem):
        dt_e = elem_dt(bar, e)
        if dt_e >= dtTarget:
            continue
        beta = (dtTarget / dt_e) ** 2 - 1.0
        md = bar._me_lumped(e)
        Me = md.sum()
        Mbar = beta * (np.diag(md) - np.outer(md, md) / Me)  # 2x2, row sums = 0
        dofs = bar.elem_dofs(e)
        for a in range(2):
            for b in range(2):
                if dofs[a] >= 0 and dofs[b] >= 0:
                    dM[dofs[a], dofs[b]] += Mbar[a, b]
        added += np.trace(Mbar)  # full element trace incl. clamped-node share
    return dM, added


# ---------------------------------------------------------------------------
# diagnostics
# ---------------------------------------------------------------------------

def f1(K, M):
    """Fundamental frequency (Hz) of K phi = omega^2 M phi. M may be a diagonal
    vector or a dense SPD matrix."""
    if M.ndim == 1:
        s = 1.0 / np.sqrt(M)
        A = (K * s[:, None]) * s[None, :]
    else:
        L = np.linalg.cholesky(M)
        Linv = np.linalg.inv(L)
        A = Linv @ K @ Linv.T
    w = np.linalg.eigvalsh(A)
    w = w[w > 1e-12]
    return np.sqrt(w.min()) / (2.0 * np.pi)


def omega_max(K, M):
    if M.ndim == 1:
        s = 1.0 / np.sqrt(M)
        A = (K * s[:, None]) * s[None, :]
    else:
        L = np.linalg.cholesky(M)
        Linv = np.linalg.inv(L)
        A = Linv @ K @ Linv.T
    return np.sqrt(np.linalg.eigvalsh(A).max())


def pcg(matvec, b, Minv_diag, tol=1e-10, maxit=1000):
    """Jacobi-preconditioned CG. Minv_diag is the lumped-diagonal preconditioner
    (vector). Returns (x, iters)."""
    x = np.zeros_like(b)
    r = b - matvec(x)
    z = Minv_diag * r
    p = z.copy()
    rz = r @ z
    bnorm = np.linalg.norm(b)
    for k in range(1, maxit + 1):
        Ap = matvec(p)
        alpha = rz / (p @ Ap)
        x += alpha * p
        r -= alpha * Ap
        if np.linalg.norm(r) <= tol * bnorm:
            return x, k
        z = Minv_diag * r
        rz_new = r @ z
        p = z + (rz_new / rz) * p
        rz = rz_new
    return x, maxit


# ---------------------------------------------------------------------------
# the oracle
# ---------------------------------------------------------------------------

def run_case(name, lengths, dtTarget):
    bar = Bar(lengths)
    K = bar.assemble_K()
    Mlump = bar.assemble_M_lumped()

    dt_e = np.array([elem_dt(bar, e) for e in range(bar.nelem)])
    nsub = int((dt_e < dtTarget).sum())

    # baselines
    f1_0 = f1(K, Mlump)
    wmax_0 = omega_max(K, Mlump)
    dt_global_0 = 2.0 / wmax_0

    # conventional (lumped) scaling
    dm_conv, added_conv = conventional_added(bar, dtTarget)
    Mconv = Mlump + dm_conv
    f1_conv = f1(K, Mconv)
    dt_global_conv = 2.0 / omega_max(K, Mconv)

    # olovsson (consistent) scaling
    dM_olo, added_olo = olovsson_added(bar, dtTarget)
    Molo = np.diag(Mlump) + dM_olo
    f1_olo = f1(K, Molo)
    dt_global_olo = 2.0 / omega_max(K, Molo)

    Mtot = Mlump.sum()

    # PCG solve M_tilde a = r with the lumped-diagonal preconditioner (matrix-free)
    rng = np.random.default_rng(0)
    b = rng.standard_normal(len(Mlump))
    Minv = 1.0 / Mlump

    def matvec(x):
        return Mlump * x + dM_olo @ x

    x, iters = pcg(matvec, b, Minv)
    res = np.linalg.norm(b - matvec(x)) / np.linalg.norm(b)
    # how a plain (un-preconditioned-ish) lumped-only solve would mis-solve M_tilde:
    a_lumped_wrong = Minv * b
    relerr_lumped = np.linalg.norm(x - a_lumped_wrong) / np.linalg.norm(x)

    print(f"\n=== {name} ===")
    print(f"  elems={bar.nelem}  sub-target={nsub}  dtTarget={dtTarget:.5g}")
    print(f"  baseline   f1={f1_0:.6f} Hz   global dt_cr={dt_global_0:.5g}")
    print(f"  conv lump  f1={f1_conv:.6f} Hz ({100*(f1_conv-f1_0)/f1_0:+.3f}%)"
          f"   added {100*added_conv/Mtot:.2f}%   dt_cr={dt_global_conv:.5g}")
    print(f"  olovsson   f1={f1_olo:.6f} Hz ({100*(f1_olo-f1_0)/f1_0:+.3f}%)"
          f"   added {100*added_olo/Mtot:.2f}%   dt_cr={dt_global_olo:.5g}")
    print(f"  PCG (lumped precond): {iters} iters, rel.resid={res:.2e}; "
          f"a lumped-only solve mis-solves M_tilde by {100*relerr_lumped:.1f}%")

    return dict(
        f1_0=f1_0, f1_conv=f1_conv, f1_olo=f1_olo,
        added_conv=added_conv, added_olo=added_olo, Mtot=Mtot,
        dt_global_0=dt_global_0, dt_global_olo=dt_global_olo,
        dtTarget=dtTarget, iters=iters, res=res, nsub=nsub,
    )


def main():
    failures = []

    def check(cond, msg):
        status = "OK  " if cond else "FAIL"
        print(f"    [{status}] {msg}")
        if not cond:
            failures.append(msg)

    # Case A: one tiny interior element among a uniform bulk (the motivating refinement).
    lengths_A = [1.0] * 5 + [0.05] + [1.0] * 5   # one element 20x smaller
    bulk_dt = 2.0 / omega_max(Bar([1.0, 1.0]).assemble_K(),
                              Bar([1.0, 1.0]).assemble_M_lumped())
    rA = run_case("A: single tiny interior element", lengths_A, dtTarget=0.012)

    # Case B: a refined zone (several small elements) -- the SSI/contact-zone regime.
    lengths_B = [1.0] * 4 + [0.04] * 6 + [1.0] * 4
    rB = run_case("B: refined zone (6 small elems)", lengths_B, dtTarget=0.012)

    print("\n--- ORACLE ASSERTIONS ---")
    for r in (rA, rB):
        # C1: Olovsson preserves f1 far better than conventional lumped.
        d_olo = abs(r["f1_olo"] - r["f1_0"]) / r["f1_0"]
        d_conv = abs(r["f1_conv"] - r["f1_0"]) / r["f1_0"]
        check(d_olo < 0.005, f"C1 Olovsson f1 drift {100*d_olo:.3f}% < 0.5%")
        check(d_olo < 0.2 * d_conv + 1e-9,
              f"C1 Olovsson f1 drift {100*d_olo:.3f}% << conv {100*d_conv:.3f}% (>5x better)")
        # C2: Olovsson injects less total mass than conventional for the same target.
        check(r["added_olo"] < r["added_conv"] + 1e-12,
              f"C2 Olovsson added {r['added_olo']:.4g} < conv {r['added_conv']:.4g}")
        # C2b: both actually deliver the bigger global step.
        check(r["dt_global_olo"] > 0.9 * r["dtTarget"],
              f"C2 Olovsson global dt_cr {r['dt_global_olo']:.4g} ~>= dtTarget {r['dtTarget']:.4g}")
        # C3: PCG converges quickly and to tolerance.
        check(r["iters"] <= 30 and r["res"] < 1e-9,
              f"C3 PCG converged in {r['iters']} iters (resid {r['res']:.1e})")

    print()
    if failures:
        print(f"ORACLE FAILED: {len(failures)} assertion(s) failed")
        raise SystemExit(1)
    print("ORACLE PASSED: all consistent-SMS claims verified")


if __name__ == "__main__":
    main()
