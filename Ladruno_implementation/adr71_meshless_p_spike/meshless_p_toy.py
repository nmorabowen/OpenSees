"""ADR-71 side-spike: FEM-u + meshless-p (MLS) toy assessment.

Question (NMB, 2026-07-13): does replacing the FE pressure interpolation of the
Biot u-p element with a MESHLESS (MLS) pressure field -- particles at the FE
nodes, at element centroids, or at the Gauss points -- make sense, and does it
work?  This is the pre-code numpy-oracle idiom of ADR-70/71 P0: same block
algebra as the shipped LadrunoUP kernel (Q/H/S/H-tilde, honest-p contract,
quasi-static), only the pressure basis is swapped.

Pressure spaces compared
  FEM-eq    : Q4 nodal equal-order (shipped baseline, no stab)
  FEM-stab  : Q4 nodal + McGann alpha-Laplacian (shipped `-stab auto`)
  MLS-node  : MLS (linear basis, cubic-spline weight) over the FE NODE cloud
  MLS-cent  : MLS over ELEMENT CENTROIDS ("smoothed P0")
  MLS-GP    : MLS over the 2x2 GAUSS-POINT cloud  <-- the idea under test

Experiments
  E1  inf-sup smoke (ADR-71 P0 pinned recipe): eigenvalues of the pressure
      Schur complement S_p = Q^T K^-1 Q on a boundary-fixed patch; count
      near-zero (spurious) pressure modes.  Anchors: FEM-eq >= 2, FEM-stab == 1.
  E2  Terzaghi 1D consolidation column vs the analytic series (works-where-
      stable check + drained-BC quality for non-interpolatory bases).
  E3  undrained footing snapshot (mini-B4): checkerboard index CB +
      settlement (locking) vs the stabilized-FEM reference.
  E4  element removal mid-consolidation: fluid measure decoupled from solid
      element life-cycle (the payoff niche) -- monolithic removal vs
      fluid-persistent FEM reference vs meshless node-cloud p.
  E5  performance sweep (10/20/40 meshes): coupling density, sparse-LU fill,
      factor + per-step solve cost of MLS-node-p vs FEM-p.  NOTE: the toy's
      assembly timings are dominated by dense-matrix scatter (not
      representative); nnz/fill/factor/solve are the meaningful columns.
  E6  the staggered seam, tried end-to-end: FEM solid + persistent FEM
      pressure overlay, fixed-stress operator split (oedometric modulus).
      Accuracy vs monolithic, drained-split divergence, the E4 crack
      capability through the seam, and 40x40 cost (SPD sub-solves).

Pure numpy/matplotlib (E5 additionally uses scipy sparse LU).
Run:  python meshless_p_toy.py [all|e1|e2|e3|e4|e5]

Results 2026-07-13 are recorded in RESULTS.md next to this file.
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = __import__("os").path.dirname(__import__("os").path.abspath(__file__))

# ===================================================================== FE: Q4
G1 = 1.0 / np.sqrt(3.0)
GPTS = np.array([[-G1, -G1], [G1, -G1], [G1, G1], [-G1, G1]])
MVEC = np.array([1.0, 1.0, 0.0])  # plane-strain Voigt [exx eyy gxy]


def mesh_rect(nx, ny, Lx, Ly):
    xs = np.linspace(0.0, Lx, nx + 1)
    ys = np.linspace(0.0, Ly, ny + 1)
    nodes = np.array([[x, y] for y in ys for x in xs])
    elems = []
    for j in range(ny):
        for i in range(nx):
            n0 = j * (nx + 1) + i
            elems.append([n0, n0 + 1, n0 + nx + 2, n0 + nx + 1])
    return nodes, np.array(elems, dtype=int)


def q4(xi, eta):
    N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                         (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
    dN = 0.25 * np.array([[-(1 - eta), -(1 - xi)],
                          [(1 - eta), -(1 + xi)],
                          [(1 + eta), (1 + xi)],
                          [-(1 + eta), (1 - xi)]])
    return N, dN


def gp_data(nodes, elems):
    """Per-GP dicts: x, dv, B(3x8), udofs(8), N(4), dNx(4x2), conn, elem."""
    out = []
    for e, conn in enumerate(elems):
        xe = nodes[conn]
        udofs = np.array([[2 * n, 2 * n + 1] for n in conn]).ravel()
        for xi, eta in GPTS:
            N, dN = q4(xi, eta)
            J = xe.T @ dN
            detJ = np.linalg.det(J)
            assert detJ > 0.0
            dNx = dN @ np.linalg.inv(J)
            B = np.zeros((3, 8))
            B[0, 0::2] = dNx[:, 0]
            B[1, 1::2] = dNx[:, 1]
            B[2, 0::2] = dNx[:, 1]
            B[2, 1::2] = dNx[:, 0]
            out.append(dict(x=xe.T @ N, dv=detJ, B=B, udofs=udofs, N=N,
                            dNx=dNx, conn=conn, elem=e))
    return out


def dmat_plane_strain(E, nu):
    c = E / ((1 + nu) * (1 - 2 * nu))
    return c * np.array([[1 - nu, nu, 0.0],
                         [nu, 1 - nu, 0.0],
                         [0.0, 0.0, (1 - 2 * nu) / 2.0]])


def assemble_K(gps, nU, D, mask=None):
    K = np.zeros((nU, nU))
    for g, d in enumerate(gps):
        if mask is not None and not mask[g]:
            continue
        ke = d["dv"] * (d["B"].T @ D @ d["B"])
        K[np.ix_(d["udofs"], d["udofs"])] += ke
    return K


# ================================================================ MLS basis
def spline_w(q):
    """Cubic B-spline weight w(q) and dw/dq, support q in [0,1)."""
    w = np.zeros_like(q)
    dw = np.zeros_like(q)
    m = q <= 0.5
    w[m] = 2.0 / 3.0 - 4.0 * q[m] ** 2 + 4.0 * q[m] ** 3
    dw[m] = -8.0 * q[m] + 12.0 * q[m] ** 2
    m = (q > 0.5) & (q < 1.0)
    w[m] = 4.0 / 3.0 - 4.0 * q[m] + 4.0 * q[m] ** 2 - (4.0 / 3.0) * q[m] ** 3
    dw[m] = -4.0 + 8.0 * q[m] - 4.0 * q[m] ** 2
    return w, dw


class MLSP:
    """MLS pressure space: linear basis, cubic-spline weight, radius rho."""

    def __init__(self, cloud, rho, label):
        self.pts = np.asarray(cloud, dtype=float)
        self.rho = float(rho)
        self.nP = len(self.pts)
        self.label = label
        self.max_condA = 0.0

    def shape(self, x, grad=True):
        d = x - self.pts                       # (nP,2), x - x_j
        r = np.hypot(d[:, 0], d[:, 1])
        q = r / self.rho
        w, dwdq = spline_w(q)
        idx = np.nonzero(w > 0.0)[0]
        if len(idx) < 3:
            raise RuntimeError(
                f"{self.label}: only {len(idx)} neighbors at {x}; increase s")
        P = np.column_stack([np.ones(len(idx)),
                             self.pts[idx, 0], self.pts[idx, 1]])
        A = (P * w[idx, None]).T @ P
        condA = np.linalg.cond(A)
        self.max_condA = max(self.max_condA, condA)
        if condA > 1e12:
            raise RuntimeError(
                f"{self.label}: MLS moment matrix cond={condA:.1e} at {x}")
        px = np.array([1.0, x[0], x[1]])
        c = np.linalg.solve(A, px)
        phi = w[idx] * (P @ c)
        if not grad:
            return idx, phi, None
        rr = np.where(r[idx] > 1e-14 * self.rho, r[idx], 1.0)
        dwdx = dwdq[idx, None] * (d[idx] / (rr[:, None] * self.rho))
        dwdx[r[idx] <= 1e-14 * self.rho] = 0.0
        dphi = np.zeros((len(idx), 2))
        for k in range(2):
            Ak = (P * dwdx[:, k][:, None]).T @ P
            pk = np.zeros(3)
            pk[k + 1] = 1.0
            ck = np.linalg.solve(A, pk - Ak @ c)
            dphi[:, k] = dwdx[:, k] * (P @ c) + w[idx] * (P @ ck)
        return idx, phi, dphi

    def eval_at(self, points, grad=False):
        Phi = np.zeros((len(points), self.nP))
        dPhi = np.zeros((len(points), self.nP, 2)) if grad else None
        for i, x in enumerate(np.atleast_2d(points)):
            idx, phi, dphi = self.shape(np.asarray(x, float), grad=grad)
            Phi[i, idx] = phi
            if grad:
                dPhi[i, idx, :] = dphi
        return (Phi, dPhi) if grad else Phi

    def self_test(self, lo, hi, n=6, seed=4):
        """Partition of unity + linear consistency at random interior points."""
        rng = np.random.default_rng(seed)
        pts = lo + (hi - lo) * rng.random((n, 2)) * 0.8 + 0.1 * (hi - lo)
        for x in pts:
            idx, phi, dphi = self.shape(x, grad=True)
            assert abs(phi.sum() - 1.0) < 1e-9, "MLS PU failed"
            err = phi @ self.pts[idx] - x
            assert np.abs(err).max() < 1e-9, "MLS linear consistency failed"
            assert np.abs(dphi.sum(axis=0)).max() < 1e-7, "MLS dPU failed"
            gg = dphi.T @ self.pts[idx]      # should be identity
            assert np.abs(gg - np.eye(2)).max() < 1e-6, "MLS grad consistency"


class FEMP:
    """Q4 nodal pressure space (interpolatory)."""

    def __init__(self, nodes, label):
        self.nodes = nodes
        self.nP = len(nodes)
        self.label = label

    def eval_at(self, points):
        """Only for points that coincide with FE nodes (all the toy needs)."""
        Phi = np.zeros((len(points), self.nP))
        for i, x in enumerate(np.atleast_2d(points)):
            j = int(np.argmin(np.linalg.norm(self.nodes - x, axis=1)))
            assert np.linalg.norm(self.nodes[j] - x) < 1e-9, \
                "FEM eval_at supports node-coincident points only"
            Phi[i, j] = 1.0
        return Phi


# ======================================================== scheme construction
def build_scheme(kind, nodes, elems, gps, s=2.0):
    """Return dict: label, nP, Phi(nG,nP), dPhi(nG,nP,2), space, is_fem."""
    h = np.max(np.abs(nodes[elems[0][1]] - nodes[elems[0][0]]))  # elem size
    if kind in ("FEM-eq", "FEM-stab"):
        sp = FEMP(nodes, kind)
        loc = [(d["conn"], d["N"], d["dNx"]) for d in gps]
        return dict(label=kind, nP=sp.nP, loc=loc, space=sp,
                    is_fem=True, stab=(kind == "FEM-stab"))
    if kind == "MLS-node":
        cloud, spacing = nodes, h
    elif kind == "MLS-cent":
        cloud = np.array([nodes[c].mean(axis=0) for c in elems])
        spacing = h
    elif kind == "MLS-GP":
        cloud = np.array([d["x"] for d in gps])
        spacing = h / 2.0
    else:
        raise ValueError(kind)
    label = f"{kind}(s={s})"
    sp = MLSP(cloud, s * spacing, label)
    sp.self_test(nodes.min(axis=0), nodes.max(axis=0))
    loc = [sp.shape(d["x"], grad=True) for d in gps]
    return dict(label=label, nP=sp.nP, loc=loc, space=sp,
                is_fem=False, stab=False)


def assemble_blocks(gps, nU, sch, kbar, invQbar, qmask=None, fmask=None):
    """Q (nU x nP), H, S, Htilde(unit-alpha) -- the LadrunoUP kernel blocks.

    qmask selects the GPs carrying SKELETON coupling (Q needs a solid);
    fmask selects the GPs carrying FLUID measure (H/S/Ht -- the fluid can
    outlive the solid element: that is the E4 architecture under test).
    """
    nP = sch["nP"]
    Q = np.zeros((nU, nP))
    H = np.zeros((nP, nP))
    S = np.zeros((nP, nP))
    Ht = np.zeros((nP, nP))
    kd = np.diag(kbar)
    for g, d in enumerate(gps):
        dv = d["dv"]
        idx, phi, dphi = sch["loc"][g]
        if qmask is None or qmask[g]:
            Q[np.ix_(d["udofs"], idx)] += dv * np.outer(d["B"].T @ MVEC, phi)
        if fmask is None or fmask[g]:
            ii = np.ix_(idx, idx)
            H[ii] += dv * (dphi @ kd @ dphi.T)
            S[ii] += dv * invQbar * np.outer(phi, phi)
            Ht[ii] += dv * (dphi @ dphi.T)
    return Q, H, S, Ht


def mcgann_alpha(E, nu, h, alpha0=0.25):
    Ks = E / (3.0 * (1.0 - 2.0 * nu))
    Gs = E / (2.0 * (1.0 + nu))
    return alpha0 * h * h / (Ks + 4.0 * Gs / 3.0)


def make_solver(A):
    """Max-abs row+column equilibration, then exact LU (via inv -- toy sizes).

    The coupled u-p-lambda matrix mixes scales by ~15 orders (K ~ E, dt*H ~
    dt*kbar); a naive cond() check would misread that as singularity and a
    pinv fallback would silently Tikhonov-regularize -- i.e. FAKE-stabilize
    the very spurious modes this toy exists to expose.  Equilibrate first,
    report the equilibrated condition number, and always solve exactly.
    """
    r = np.abs(A).max(axis=1)
    r[r == 0.0] = 1.0
    A1 = A / r[:, None]
    c = np.abs(A1).max(axis=0)
    c[c == 0.0] = 1.0
    A2 = A1 / c[None, :]
    cond2 = np.linalg.cond(A2)
    Ai = np.linalg.inv(A2)

    def solve(b):
        return (Ai @ (b / r)) / c

    return solve, cond2


# ============================================================ E1: inf-sup
def experiment_infsup():
    print("=" * 78)
    print("E1  inf-sup smoke (ADR-71 P0 recipe): 8x8 patch, u fixed on the")
    print("    whole boundary, p free; eigen-count of S_p = Q^T K^-1 Q")
    print("=" * 78)
    nx = ny = 8
    nodes, elems = mesh_rect(nx, ny, 1.0, 1.0)
    gps = gp_data(nodes, elems)
    E, nu = 1.0e6, 0.3
    h = 1.0 / nx
    K = assemble_K(gps, 2 * len(nodes), dmat_plane_strain(E, nu))
    bnd = np.nonzero((np.abs(nodes[:, 0]) < 1e-12) |
                     (np.abs(nodes[:, 0] - 1) < 1e-12) |
                     (np.abs(nodes[:, 1]) < 1e-12) |
                     (np.abs(nodes[:, 1] - 1) < 1e-12))[0]
    fixed = np.array([[2 * n, 2 * n + 1] for n in bnd]).ravel()
    free = np.setdiff1d(np.arange(2 * len(nodes)), fixed)
    Kff = K[np.ix_(free, free)]
    alpha = mcgann_alpha(E, nu, h)

    kinds = [("FEM-eq", None), ("FEM-stab", None),
             ("MLS-node", 1.6), ("MLS-node", 2.4),
             ("MLS-cent", 1.6), ("MLS-cent", 2.4),
             ("MLS-GP", 1.6), ("MLS-GP", 2.4), ("MLS-GP", 3.2)]
    rows, spectra = [], []
    for kind, s in kinds:
        sch = build_scheme(kind, nodes, elems, gps,
                           s=(s if s is not None else 2.0))
        Q, Hm, S, Ht = assemble_blocks(gps, 2 * len(nodes), sch,
                                       kbar=[1.0, 1.0], invQbar=0.0)
        Qf = Q[free, :]
        Sp = Qf.T @ np.linalg.solve(Kff, Qf)
        Sp = 0.5 * (Sp + Sp.T)
        if sch["stab"]:
            Sp = Sp + alpha * Ht
        ev = np.linalg.eigvalsh(Sp)
        nzero = int((ev < 1e-10 * ev.max()).sum())
        rankQ = np.linalg.matrix_rank(Qf, tol=1e-10 * np.linalg.norm(Qf, 2))
        floor = max(0, sch["nP"] - len(free))   # structural lower bound
        rows.append((sch["label"], sch["nP"], len(free), rankQ, floor, nzero))
        spectra.append((sch["label"], ev / ev.max()))
        print(f"  {sch['label']:<18} nP={sch['nP']:>3}  freeU={len(free):>3}  "
              f"rank(Q)={rankQ:>3}  structural-floor={floor:>3}  "
              f"ZERO-MODES={nzero:>3}")
    print("  gate anchors: FEM-eq >= 2 (const + checkerboard), FEM-stab == 1")

    fig, ax = plt.subplots(figsize=(8, 5))
    for label, ev in spectra:
        ax.semilogy(np.maximum(ev, 1e-20), marker=".", ms=3, lw=1, label=label)
    ax.axhline(1e-10, color="k", ls="--", lw=0.8, label="zero-mode threshold")
    ax.set_xlabel("mode #")
    ax.set_ylabel("eig(S_p)/max")
    ax.set_title("E1: pressure Schur-complement spectra (boundary-fixed patch)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(f"{HERE}/infsup_spectra.png", dpi=150)
    plt.close(fig)
    return rows


# ============================================================ E2: Terzaghi
def terzaghi_series(Z, Tv, terms=300):
    p = np.zeros_like(Z)
    for m in range(terms):
        M = 0.5 * np.pi * (2 * m + 1)
        p += (2.0 / M) * np.sin(M * Z) * np.exp(-M * M * Tv)
    return p


def transient_run(nodes, elems, gps, sch, *, E, nu, kbar, invQbar, fvec,
                  ufix, cpoints, dt, nsteps, alpha_stab=0.0, ramp_t=None,
                  sample_steps=()):
    """Backward-Euler on  K u - Q p = f ;  Q^T du + (S+aHt) dp + dt H p = 0,
    drained BC via Lagrange multipliers at cpoints.  Returns samples {step: p}."""
    nUall = 2 * len(nodes)
    K = assemble_K(gps, nUall, dmat_plane_strain(E, nu))
    Q, H, S, Ht = assemble_blocks(gps, nUall, sch, kbar, invQbar)
    Sstar = S + alpha_stab * Ht
    if sch["is_fem"]:
        C = sch["space"].eval_at(cpoints)
    else:
        C = sch["space"].eval_at(cpoints)
    ufree = np.setdiff1d(np.arange(nUall), ufix)
    nu_, np_, nc = len(ufree), sch["nP"], len(cpoints)
    A = np.zeros((nu_ + np_ + nc, nu_ + np_ + nc))
    A[:nu_, :nu_] = K[np.ix_(ufree, ufree)]
    A[:nu_, nu_:nu_ + np_] = -Q[ufree, :]
    A[nu_:nu_ + np_, :nu_] = Q[ufree, :].T
    A[nu_:nu_ + np_, nu_:nu_ + np_] = Sstar + dt * H
    A[nu_:nu_ + np_, nu_ + np_:] = dt * C.T
    A[nu_ + np_:, nu_:nu_ + np_] = C
    solve, condA = make_solver(A)
    u = np.zeros(nu_)
    p = np.zeros(np_)
    samples = {}
    for step in range(1, nsteps + 1):
        t = step * dt
        lam = 1.0 if ramp_t is None else min(t / ramp_t, 1.0)
        rhs = np.concatenate([lam * fvec[ufree],
                              Q[ufree, :].T @ u + Sstar @ p,
                              np.zeros(nc)])
        z = solve(rhs)
        u, p = z[:nu_], z[nu_:nu_ + np_]
        if step in sample_steps:
            samples[step] = (u.copy(), p.copy())
    return samples, ufree, condA


def experiment_terzaghi():
    print("=" * 78)
    print("E2  Terzaghi column (2x20 Q4, 1x10 m), q=100 kPa step, drained top")
    print("=" * 78)
    nx, ny, Lx, Ly = 2, 20, 1.0, 10.0
    nodes, elems = mesh_rect(nx, ny, Lx, Ly)
    gps = gp_data(nodes, elems)
    E, nu = 1.0e7, 0.25
    kh, gw = 1.0e-7, 1.0e4
    kbar = [kh / gw] * 2
    Kf, n = 2.2e9, 0.4
    invQbar = n / Kf
    q = 1.0e5
    Eoed = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu))
    cv = (kh / gw) * Eoed
    t_half = 0.5 * Ly * Ly / cv
    nsteps = 400
    dt = t_half / nsteps
    tv_at = {80: 0.1, 400: 0.5}           # step -> Tv

    top = np.nonzero(np.abs(nodes[:, 1] - Ly) < 1e-9)[0]
    bot = np.nonzero(np.abs(nodes[:, 1]) < 1e-9)[0]
    side = np.nonzero((np.abs(nodes[:, 0]) < 1e-9) |
                      (np.abs(nodes[:, 0] - Lx) < 1e-9))[0]
    ufix = np.unique(np.concatenate(
        [[2 * n_, 2 * n_ + 1] for n_ in bot] +
        [[2 * n_] for n_ in side]))
    f = np.zeros(2 * len(nodes))
    dx = Lx / nx
    for n_ in top:
        w = dx if 1e-9 < nodes[n_, 0] < Lx - 1e-9 else dx / 2
        f[2 * n_ + 1] = -q * w
    cpoints = nodes[top]                   # drained top
    line = np.nonzero(np.abs(nodes[:, 0] - Lx / 2) < 1e-9)[0]
    line = line[np.argsort(nodes[line, 1])]
    Z = (Ly - nodes[line, 1]) / Ly         # depth from drained face / H
    bc_probe = np.column_stack([np.linspace(0, Lx, 21), np.full(21, Ly)])

    alpha = mcgann_alpha(E, nu, Ly / ny)
    results = {}
    for kind, s in [("FEM-eq", None), ("FEM-stab", None),
                    ("MLS-node", 2.0), ("MLS-cent", 2.0), ("MLS-GP", 2.0)]:
        sch = build_scheme(kind, nodes, elems, gps, s=(s or 2.0))
        a = alpha if sch["stab"] else 0.0
        samples, ufree, condA = transient_run(
            nodes, elems, gps, sch, E=E, nu=nu, kbar=kbar, invQbar=invQbar,
            fvec=f, ufix=ufix, cpoints=cpoints, dt=dt, nsteps=nsteps,
            alpha_stab=a, sample_steps=set(tv_at))
        profs, errs = {}, {}
        for step, tv in tv_at.items():
            _, p = samples[step]
            if sch["is_fem"]:
                pl = p[line]
            else:
                pl = sch["space"].eval_at(nodes[line]) @ p
            pa = q * terzaghi_series(Z, tv)
            profs[tv] = pl
            errs[tv] = np.linalg.norm(pl - pa) / np.linalg.norm(pa)
        if sch["is_fem"]:
            bcq = 0.0                      # interpolatory: exactly 0 at nodes
        else:
            _, p = samples[80]
            bcq = np.abs(sch["space"].eval_at(bc_probe) @ p).max() / q
        results[sch["label"]] = dict(profs=profs, errs=errs, bcq=bcq)
        print(f"  {sch['label']:<18} relL2(Tv=0.1)={errs[0.1]:.3e}  "
              f"relL2(Tv=0.5)={errs[0.5]:.3e}  "
              f"drained-face |p|/q={bcq:.2e}  cond(A)={condA:.1e}")

    fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
    for ax, tv in zip(axes, [0.1, 0.5]):
        Zf = np.linspace(0, 1, 200)
        ax.plot(terzaghi_series(Zf, tv), Zf, "k-", lw=2, label="analytic")
        for label, r in results.items():
            ax.plot(r["profs"][tv] / q, Z, marker="o", ms=3, lw=1, label=label)
        ax.set_title(f"Tv = {tv}")
        ax.set_xlabel("p / q")
        ax.invert_yaxis()
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("depth from drained face / H")
    axes[0].legend(fontsize=7)
    fig.suptitle("E2: Terzaghi consolidation, all pressure spaces")
    fig.tight_layout()
    fig.savefig(f"{HERE}/terzaghi_profiles.png", dpi=150)
    plt.close(fig)
    return results


# ============================================================ E3: footing
def experiment_footing():
    print("=" * 78)
    print("E3  undrained footing snapshot (mini-B4): 10x10 m, strip load 3 m,")
    print("    k'=1e-7 m/s, Kf=2.2e15 Pa, t=1.0 s  ->  CB index + settlement")
    print("=" * 78)
    nx = ny = 10
    Lx = Ly = 10.0
    nodes, elems = mesh_rect(nx, ny, Lx, Ly)
    gps = gp_data(nodes, elems)
    E, nu = 2.5e7, 0.3
    kh, gw = 1.0e-7, 1.0e4
    kbar = [kh / gw] * 2
    Kf, n = 2.2e15, 0.4
    invQbar = n / Kf
    q = 100.0                                # 0.1 kPa
    dt, nsteps = 0.05, 20                    # t_end = 1.0 s, ramp 0.1 s

    top = np.nonzero(np.abs(nodes[:, 1] - Ly) < 1e-9)[0]
    bot = np.nonzero(np.abs(nodes[:, 1]) < 1e-9)[0]
    lft = np.nonzero(np.abs(nodes[:, 0]) < 1e-9)[0]
    rgt = np.nonzero(np.abs(nodes[:, 0] - Lx) < 1e-9)[0]
    ufix = np.unique(np.concatenate(
        [[2 * n_, 2 * n_ + 1] for n_ in bot] +
        [[2 * n_] for n_ in np.concatenate([lft, rgt])]))
    f = np.zeros(2 * len(nodes))
    dx = Lx / nx
    for n_ in top:
        x = nodes[n_, 0]
        if x < 3.0 - 1e-9:
            w = dx if x > 1e-9 else dx / 2
        elif abs(x - 3.0) < 1e-9:
            w = dx / 2
        else:
            continue
        f[2 * n_ + 1] = -q * w
    cpoints = nodes[top]                     # drained top surface only
    corner = top[np.argmin(nodes[top, 0])]   # (0, Ly): under load center

    interior = [(i, j, j * (nx + 1) + i)
                for j in range(1, ny) for i in range(1, nx)]
    chi = np.array([(-1.0) ** (i + j) for i, j, _ in interior])
    int_ids = np.array([nid for _, _, nid in interior])

    alpha = mcgann_alpha(E, nu, dx)
    results = {}
    for kind, s in [("FEM-eq", None), ("FEM-stab", None),
                    ("MLS-node", 2.0), ("MLS-cent", 2.0),
                    ("MLS-GP", 2.0), ("MLS-GP", 2.4)]:
        sch = build_scheme(kind, nodes, elems, gps, s=(s or 2.0))
        a = alpha if sch["stab"] else 0.0
        samples, ufree, condA = transient_run(
            nodes, elems, gps, sch, E=E, nu=nu, kbar=kbar, invQbar=invQbar,
            fvec=f, ufix=ufix, cpoints=cpoints, dt=dt, nsteps=nsteps,
            alpha_stab=a, ramp_t=0.1, sample_steps={nsteps})
        u, p = samples[nsteps]
        if sch["is_fem"]:
            pnod = p.copy()
        else:
            pnod = sch["space"].eval_at(nodes) @ p
        pint = pnod[int_ids]
        cb = abs(pint @ chi) / (np.linalg.norm(pint) * np.linalg.norm(chi))
        iu = np.nonzero(ufree == 2 * corner + 1)[0][0]
        sett = u[iu]
        results[sch["label"]] = dict(pnod=pnod, cb=cb, sett=sett, condA=condA)
        print(f"  {sch['label']:<18} CB={cb:.4f}  settlement={sett:+.4e} m  "
              f"cond(A)={condA:.2e}")
    ref = results["FEM-stab"]["sett"]
    print(f"  settlement ratios vs FEM-stab reference ({ref:+.4e} m):")
    for label, r in results.items():
        print(f"    {label:<18} ratio={r['sett'] / ref:6.3f}   CB={r['cb']:.4f}")

    labels = list(results)
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    vmax = max(np.abs(results[l]["pnod"]).max() for l in labels[:4])
    for ax, label in zip(axes.ravel(), labels):
        r = results[label]
        img = r["pnod"].reshape(ny + 1, nx + 1)
        vm = min(np.abs(r["pnod"]).max(), 3 * vmax)
        im = ax.imshow(img, origin="lower", cmap="RdBu_r", vmin=-vm, vmax=vm,
                       extent=[0, Lx, 0, Ly])
        ax.set_title(f"{label}\nCB={r['cb']:.3f}  "
                     f"sett/ref={r['sett'] / ref:.2f}", fontsize=9)
        fig.colorbar(im, ax=ax, shrink=0.8)
    fig.suptitle("E3: undrained footing pore pressure at t=1 s "
                 "(nodal-lattice evaluation)")
    fig.tight_layout()
    fig.savefig(f"{HERE}/footing_maps.png", dpi=150)
    plt.close(fig)
    return results


# ==================================================== E4: element removal
def build_coupled_A(K, Q, Sstar, H, C, ufree, dt):
    nu_, np_, nc = len(ufree), Q.shape[1], C.shape[0]
    A = np.zeros((nu_ + np_ + nc, nu_ + np_ + nc))
    A[:nu_, :nu_] = K[np.ix_(ufree, ufree)]
    A[:nu_, nu_:nu_ + np_] = -Q[ufree, :]
    A[nu_:nu_ + np_, :nu_] = Q[ufree, :].T
    A[nu_:nu_ + np_, nu_:nu_ + np_] = Sstar + dt * H
    A[nu_:nu_ + np_, nu_ + np_:] = dt * C.T
    A[nu_ + np_:, nu_:nu_ + np_] = C
    return A


def experiment_removal():
    """E4 -- the payoff niche: a pressure field whose life-cycle is decoupled
    from the solid elements.  10x10 m column, uniform load, drained top only;
    at Tv=0.05 a 'crack' of 8 elements (80 % width, mid-height) is removed.

      FEM-mono   : upstream `remove element` semantics -- the element carried
                   BOTH the skeleton and the fluid, so the flow path dies too.
      FEM-fluid+ : same FEM p basis, but the fluid blocks keep integrating
                   over the crack cells (needs an architecture upstream/ADR-71
                   monolithic elements do not have -- shown as the reference).
      MLS-node   : meshless node-cloud p, fluid measure persists (the meshless
                   architecture: particle connectivity does not know elements).

    Toy simplification in the crack: fluid keeps H and S (water-filled gap,
    same k), skeleton K and coupling Q drop.  Gap-volume coupling neglected.
    """
    print("=" * 78)
    print("E4  element removal mid-consolidation: does the flow path survive?")
    print("=" * 78)
    nx = ny = 10
    Lx = Ly = 10.0
    nodes, elems = mesh_rect(nx, ny, Lx, Ly)
    gps = gp_data(nodes, elems)
    nUall = 2 * len(nodes)
    E, nu = 2.5e7, 0.3
    kh, gw = 1.0e-7, 1.0e4
    kbar = [kh / gw] * 2
    Kf, n = 2.2e9, 0.4
    invQbar = n / Kf
    q = 1.0e4
    Eoed = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu))
    cv = (kh / gw) * Eoed
    tv_end = 1.5
    t_end = tv_end * Ly * Ly / cv
    nsteps = 450
    dt = t_end / nsteps
    step_rm = int(round(0.05 / tv_end * nsteps))    # removal at Tv = 0.05
    step_map = int(round(0.20 / tv_end * nsteps))   # field snapshot Tv = 0.20

    # crack: element row j=4 (y in [4,5]), i = 0..7  (80 % width, from left)
    crack_elems = {4 * nx + i for i in range(8)}
    surv = np.array([d["elem"] not in crack_elems for d in gps])

    top = np.nonzero(np.abs(nodes[:, 1] - Ly) < 1e-9)[0]
    bot = np.nonzero(np.abs(nodes[:, 1]) < 1e-9)[0]
    sides = np.nonzero((np.abs(nodes[:, 0]) < 1e-9) |
                       (np.abs(nodes[:, 0] - Lx) < 1e-9))[0]
    ufix = np.unique(np.concatenate(
        [[2 * m, 2 * m + 1] for m in bot] + [[2 * m] for m in sides]))
    ufree = np.setdiff1d(np.arange(nUall), ufix)
    f = np.zeros(nUall)
    dx = Lx / nx
    for m in top:
        w = dx if 1e-9 < nodes[m, 0] < Lx - 1e-9 else dx / 2
        f[2 * m + 1] = -q * w
    cpoints = nodes[top]
    probe = np.array([2.0, 2.0])          # below the crack, left side
    alpha = mcgann_alpha(E, nu, dx)

    variants = [
        ("FEM-mono",   "FEM-stab", surv, surv),   # fluid dies with the solid
        ("FEM-fluid+", "FEM-stab", surv, None),   # fluid persists (reference)
        ("MLS-node",   "MLS-node", surv, None),   # meshless architecture
        ("no-removal", "FEM-stab", None, None),   # undamaged baseline
    ]
    results = {}
    for label, kind, qm, fm in variants:
        sch = build_scheme(kind, nodes, elems, gps, s=2.0)
        a = alpha if sch["stab"] else 0.0
        C = sch["space"].eval_at(cpoints)
        K0 = assemble_K(gps, nUall, dmat_plane_strain(E, nu))
        Q0, H0, S0, Ht0 = assemble_blocks(gps, nUall, sch, kbar, invQbar)
        solve_pre, _ = make_solver(
            build_coupled_A(K0, Q0, S0 + a * Ht0, H0, C, ufree, dt))
        if label == "no-removal":
            solve_post, Qpost, Spost = solve_pre, Q0, S0 + a * Ht0
        else:
            qmask = qm if qm is not None else None
            fmask = fm
            K1 = assemble_K(gps, nUall, dmat_plane_strain(E, nu), mask=surv)
            Q1, H1, S1, Ht1 = assemble_blocks(gps, nUall, sch, kbar, invQbar,
                                              qmask=surv, fmask=fmask)
            Spost = S1 + a * Ht1
            solve_post, cnd = make_solver(
                build_coupled_A(K1, Q1, Spost, H1, C, ufree, dt))
            Qpost = Q1
        nu_, np_ = len(ufree), sch["nP"]
        prow = (sch["space"].eval_at(probe[None, :])[0] if not sch["is_fem"]
                else sch["space"].eval_at(probe[None, :])[0])
        u = np.zeros(nu_)
        p = np.zeros(np_)
        curve = []
        snap = None
        Spre = S0 + a * Ht0
        for step in range(1, nsteps + 1):
            post = step > step_rm and label != "no-removal"
            Quse = Qpost if post else Q0
            Suse = Spost if post else Spre
            solve = solve_post if post else solve_pre
            rhs = np.concatenate([f[ufree],
                                  Quse[ufree, :].T @ u + Suse @ p,
                                  np.zeros(len(cpoints))])
            z = solve(rhs)
            u, p = z[:nu_], z[nu_:nu_ + np_]
            curve.append(prow @ p)
            if step == step_map:
                snap = (sch["space"].eval_at(nodes) @ p if not sch["is_fem"]
                        else p.copy())
        curve = np.array(curve)
        i90 = np.nonzero(curve <= 0.1 * q)[0]
        tv90 = tv_end * (i90[0] + 1) / nsteps if len(i90) else np.inf
        results[label] = dict(curve=curve, snap=snap, tv90=tv90)
        print(f"  {label:<12} p(probe, Tv=0.2)/q = {curve[step_map-1]/q:5.3f}   "
              f"p(end)/q = {curve[-1]/q:5.3f}   Tv @ 90% dissipated = {tv90:5.3f}")

    ref = results["FEM-fluid+"]
    mls = results["MLS-node"]
    dev = (np.linalg.norm(mls["curve"] - ref["curve"])
           / np.linalg.norm(ref["curve"]))
    print(f"  MLS-node vs FEM-fluid+ probe-curve deviation: {dev:.2%}"
          f"   (same architecture, different basis)")

    tvs = tv_end * np.arange(1, nsteps + 1) / nsteps
    fig = plt.figure(figsize=(13, 8))
    ax = fig.add_subplot(2, 1, 1)
    for label in results:
        ax.plot(tvs, results[label]["curve"] / q, lw=2, label=label)
    ax.axvline(0.05, color="k", ls="--", lw=0.8)
    ax.text(0.052, 0.55, "crack opens\n(8 elements removed)", fontsize=8)
    ax.set_xlabel("Tv")
    ax.set_ylabel("p(probe) / q")
    ax.set_title("E4: pore pressure below the crack -- does the flow path "
                 "survive element removal?")
    ax.legend()
    ax.grid(alpha=0.3)
    for k, label in enumerate(["FEM-mono", "FEM-fluid+", "MLS-node"]):
        axm = fig.add_subplot(2, 3, 4 + k)
        img = results[label]["snap"].reshape(ny + 1, nx + 1)
        im = axm.imshow(img / q, origin="lower", cmap="Reds",
                        vmin=0, vmax=1.0, extent=[0, Lx, 0, Ly])
        axm.add_patch(plt.Rectangle((0, 4), 8, 1, fill=False,
                                    edgecolor="blue", lw=1.5, ls="--"))
        axm.set_title(f"{label}  (p/q at Tv=0.2)", fontsize=9)
        fig.colorbar(im, ax=axm, shrink=0.8)
    fig.tight_layout()
    fig.savefig(f"{HERE}/removal_flowpath.png", dpi=150)
    plt.close(fig)
    return results


# ======================================================== E5: performance
def experiment_perf():
    """E5 -- price the meshless pressure field: assembly cost, coupling
    density (nnz), sparse-LU factor + per-step solve cost vs FEM, over a
    mesh sweep.  Uses scipy sparse LU (the honest solver model for the
    unsymmetric coupled system)."""
    import time
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla

    print("=" * 78)
    print("E5  performance: FEM-p vs MLS-node-p (s=2.0), coupled-system cost")
    print("=" * 78)
    E, nu = 2.5e7, 0.3
    kh, gw = 1.0e-7, 1.0e4
    kbar = [kh / gw] * 2
    invQbar = 0.4 / 2.2e9
    rows = []
    for nx in (10, 20, 40):
        nodes, elems = mesh_rect(nx, nx, 10.0, 10.0)
        gps = gp_data(nodes, elems)
        nUall = 2 * len(nodes)
        top = np.nonzero(np.abs(nodes[:, 1] - 10.0) < 1e-9)[0]
        bot = np.nonzero(np.abs(nodes[:, 1]) < 1e-9)[0]
        ufix = np.unique(np.concatenate(
            [[2 * m, 2 * m + 1] for m in bot]))
        ufree = np.setdiff1d(np.arange(nUall), ufix)
        dt = 100.0
        K = assemble_K(gps, nUall, dmat_plane_strain(E, nu))
        for kind in ("FEM-eq", "MLS-node"):
            t0 = time.perf_counter()
            sch = build_scheme(kind, nodes, elems, gps, s=2.0)
            t_shape = time.perf_counter() - t0
            t0 = time.perf_counter()
            Q, H, S, Ht = assemble_blocks(gps, nUall, sch, kbar, invQbar)
            t_asm = time.perf_counter() - t0
            C = sch["space"].eval_at(nodes[top])
            A = build_coupled_A(K, Q, S, H, C, ufree, dt)
            # equilibrate (as the transient runner does), then sparse LU
            r = np.abs(A).max(axis=1)
            r[r == 0] = 1.0
            c = np.abs(A / r[:, None]).max(axis=0)
            c[c == 0] = 1.0
            A2 = sp.csc_matrix(np.where(
                np.abs(A / r[:, None] / c[None, :]) > 1e-14,
                A / r[:, None] / c[None, :], 0.0))
            hnnz = int((np.abs(H) > 1e-12 * np.abs(H).max()).sum())
            t0 = time.perf_counter()
            lu = spla.splu(A2)
            t_fac = time.perf_counter() - t0
            b = np.ones(A2.shape[0])
            t0 = time.perf_counter()
            for _ in range(10):
                lu.solve(b)
            t_sol = (time.perf_counter() - t0) / 10
            fill = lu.L.nnz + lu.U.nnz
            rows.append((nx, sch["label"], sch["nP"], hnnz / sch["nP"],
                         A2.nnz, fill, t_shape + t_asm, t_fac, t_sol))
            print(f"  {nx:>2}x{nx:<2} {sch['label']:<16} nP={sch['nP']:>5}  "
                  f"nnz(H)/row={hnnz / sch['nP']:6.1f}  nnz(A)={A2.nnz:>8}  "
                  f"LUfill={fill:>9}  asm={t_shape + t_asm:6.2f}s  "
                  f"factor={t_fac * 1e3:7.1f}ms  step-solve={t_sol * 1e3:6.2f}ms")
    print("  ratios (MLS-node / FEM) per mesh:")
    for i in range(0, len(rows), 2):
        f_, m_ = rows[i], rows[i + 1]
        print(f"    {f_[0]:>2}x{f_[0]:<2}  nnz(A) x{m_[4] / f_[4]:4.1f}   "
              f"LU-fill x{m_[5] / f_[5]:4.1f}   assembly x{m_[6] / f_[6]:5.1f}   "
              f"factor x{m_[7] / f_[7]:4.1f}   step-solve x{m_[8] / f_[8]:4.1f}")

    ns = sorted(set(r[0] for r in rows))
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    for j, (title, idx, unit) in enumerate(
            [("LU fill (nnz of L+U)", 5, "nnz"),
             ("factorization time", 7, "s"),
             ("per-step solve time", 8, "s")]):
        ax = axes[j]
        for kind in ("FEM-eq", "MLS-node(s=2.0)"):
            ys = [r[idx] for r in rows if r[1] == kind]
            ax.loglog(ns, ys, marker="o", label=kind)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("mesh n (n x n)")
        ax.set_ylabel(unit)
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=8)
    fig.suptitle("E5: cost of the meshless pressure field (coupled system)")
    fig.tight_layout()
    fig.savefig(f"{HERE}/perf_scaling.png", dpi=150)
    plt.close(fig)
    return rows


# ==================================================== E6: staggered seam
def experiment_staggered():
    """E6 -- try the recommended architecture: FEM solid + persistent FEM
    pressure overlay coupled through the effective-stress seam by operator
    splitting (staggered), exactly the ADR-71 §6 / P7 shape.

    What crosses the seam per step (and nothing else):
      solid -> fluid : Q^T (u^{k} - u^n)   (volumetric strain increment)
      fluid -> solid : Q p^{k}             (GP pressures as skeleton load)

    Splits tested:
      drained : solid (p frozen) -> fluid.  The naive split -- known to be
                unstable when coupling is strong (near-undrained).
      fs1     : fixed-stress, single pass.  Fluid solve augmented with the
                pseudo-compressibility L = (alpha^2/K_dr) M_p (Kim-Tchelepi-
                Juanes); unconditionally stable, O(dt) splitting error.
      fs      : fixed-stress iterated to tol -- fixed point IS the
                monolithic backward-Euler solution.

    Parts: A Terzaghi accuracy, A2 near-undrained stability stress test,
    B the E4 crack capability through the seam, C 40x40 cost comparison
    (each staggered sub-solve is SPD -- the honest-p contract's unsymmetric
    coupled solver is not needed on this path)."""
    import time
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla

    print("=" * 78)
    print("E6  staggered seam: FEM solid + persistent pressure overlay "
          "(fixed-stress split)")
    print("=" * 78)

    def ops_build(nodes, elems, gps, sch, E, nu, kbar, invQbar, a_stab,
                  ufix, ptop, qmask=None, fmask=None):
        nUall = 2 * len(nodes)
        K = assemble_K(gps, nUall, dmat_plane_strain(E, nu), mask=qmask)
        Qb, Hb, Mp, Htb = assemble_blocks(gps, nUall, sch, kbar, 1.0,
                                          qmask=qmask, fmask=fmask)
        Sst = invQbar * Mp + a_stab * Htb
        ufree = np.setdiff1d(np.arange(nUall), ufix)
        pfree = np.setdiff1d(np.arange(sch["nP"]), ptop)
        ix = np.ix_
        return dict(K=K[ix(ufree, ufree)], Q=Qb[ix(ufree, pfree)],
                    H=Hb[ix(pfree, pfree)], S=Sst[ix(pfree, pfree)],
                    Mp=Mp[ix(pfree, pfree)], ufree=ufree, pfree=pfree)

    def mono_march(opsA, opsB, step_rm, ff, dt, nsteps):
        def make(o):
            n_u = o["K"].shape[0]
            n_p = o["S"].shape[0]
            A = np.zeros((n_u + n_p, n_u + n_p))
            A[:n_u, :n_u] = o["K"]
            A[:n_u, n_u:] = -o["Q"]
            A[n_u:, :n_u] = o["Q"].T
            A[n_u:, n_u:] = o["S"] + dt * o["H"]
            return make_solver(A)[0]
        sA = make(opsA)
        sB = make(opsB) if opsB is not None else sA
        n_u = opsA["K"].shape[0]
        u = np.zeros(n_u)
        p = np.zeros(opsA["S"].shape[0])
        Ph = np.zeros((nsteps, len(p)))
        for step in range(1, nsteps + 1):
            o, s_ = ((opsB, sB) if (opsB is not None and step > step_rm)
                     else (opsA, sA))
            z = s_(np.concatenate([ff, o["Q"].T @ u + o["S"] @ p]))
            u, p = z[:n_u], z[n_u:]
            Ph[step - 1] = p
        return Ph

    def stag_march(opsA, opsB, step_rm, ff, dt, nsteps, mode, lfac,
                   pscale, tol=1e-6, maxit=200):
        def make(o):
            L = lfac * o["Mp"]
            Ks = make_solver(o["K"])[0]
            Fs = make_solver(o["S"] + dt * o["H"]
                             + (L if mode != "drained" else 0.0))[0]
            return Ks, Fs, L
        mA = make(opsA)
        mB = make(opsB) if opsB is not None else mA
        u = np.zeros(opsA["K"].shape[0])
        p = np.zeros(opsA["S"].shape[0])
        Ph = np.full((nsteps, len(p)), np.nan)
        iters = []
        diverged = False
        for step in range(1, nsteps + 1):
            o, (Ks, Fs, L) = ((opsB, mB)
                              if (opsB is not None and step > step_rm)
                              else (opsA, mA))
            if mode == "drained":
                u1 = Ks(ff + o["Q"] @ p)
                p1 = Fs(o["S"] @ p - o["Q"].T @ (u1 - u))
                nit = 1
            else:
                # solid-first fixed-stress sweep (fluid-first with a stale u
                # predictor self-"converges" to p=0 on the first iterate)
                pk = p
                for nit in range(1, maxit + 1):
                    u1 = Ks(ff + o["Q"] @ pk)
                    p1 = Fs(o["S"] @ p + L @ pk - o["Q"].T @ (u1 - u))
                    dcon = (np.linalg.norm(p1 - pk)
                            / (np.linalg.norm(p1) + pscale))
                    pk = p1
                    if mode == "fs1" or dcon < tol:
                        break
                u1 = Ks(ff + o["Q"] @ p1)   # final momentum resolve
            u, p = u1, p1
            iters.append(nit)
            Ph[step - 1] = p
            if (not np.all(np.isfinite(p))
                    or np.linalg.norm(p) > 1e9 * pscale):
                diverged = True
                break
        return Ph, np.array(iters), diverged

    # ---------------- Part A: Terzaghi accuracy (compressible fluid) -------
    nx, ny, Lx, Ly = 2, 20, 1.0, 10.0
    nodes, elems = mesh_rect(nx, ny, Lx, Ly)
    gps = gp_data(nodes, elems)
    E, nu = 1.0e7, 0.25
    kh, gw = 1.0e-7, 1.0e4
    kbar = [kh / gw] * 2
    Kf, poro = 2.2e9, 0.4
    invQbar = poro / Kf
    q = 1.0e5
    Kdr = E / (3.0 * (1.0 - 2.0 * nu))
    G = E / (2.0 * (1.0 + nu))
    # fixed-stress modulus: oedometric variant Kdr+4G/3 -- measured 3.2 iters
    # vs 11.0 for classic a^2/Kdr on this problem; 0.5x oedometric diverges
    # (near-optimal, consistent with the split literature)
    lfac = 1.0 / (Kdr + 4.0 * G / 3.0)
    Eoed = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu))
    cv = (kh / gw) * Eoed
    nsteps = 400
    dt = (0.5 * Ly * Ly / cv) / nsteps
    top = np.nonzero(np.abs(nodes[:, 1] - Ly) < 1e-9)[0]
    bot = np.nonzero(np.abs(nodes[:, 1]) < 1e-9)[0]
    side = np.nonzero((np.abs(nodes[:, 0]) < 1e-9) |
                      (np.abs(nodes[:, 0] - Lx) < 1e-9))[0]
    ufix = np.unique(np.concatenate(
        [[2 * m, 2 * m + 1] for m in bot] + [[2 * m] for m in side]))
    f = np.zeros(2 * len(nodes))
    dx = Lx / nx
    for m in top:
        w = dx if 1e-9 < nodes[m, 0] < Lx - 1e-9 else dx / 2
        f[2 * m + 1] = -q * w
    sch = build_scheme("FEM-stab", nodes, elems, gps)
    a_stab = mcgann_alpha(E, nu, Ly / ny)
    ops = ops_build(nodes, elems, gps, sch, E, nu, kbar, invQbar, a_stab,
                    ufix, top)
    ff = f[ops["ufree"]]
    line = np.nonzero(np.abs(nodes[:, 0] - Lx / 2) < 1e-9)[0]
    line = line[np.argsort(nodes[line, 1])]
    Z = (Ly - nodes[line, 1]) / Ly

    def full_p(pfree_vec):
        pf = np.zeros(sch["nP"])
        pf[ops["pfree"]] = pfree_vec
        return pf

    P_m = mono_march(ops, None, 0, ff, dt, nsteps)
    runs = {"mono": (P_m, None, False)}
    for mode in ("fs", "fs1", "drained"):
        Ph, its, div = stag_march(ops, None, 0, ff, dt, nsteps, mode, lfac, q)
        runs[mode] = (Ph, its, div)
    print("  A  Terzaghi (2x20, Kf=2.2e9):")
    pa = q * terzaghi_series(Z, 0.1)
    for name, (Ph, its, div) in runs.items():
        if div:
            print(f"     {name:<8} DIVERGED at step {len(its)}")
            continue
        e_an = (np.linalg.norm(full_p(Ph[79])[line] - pa)
                / np.linalg.norm(pa))
        e_mo = np.linalg.norm(Ph - P_m) / np.linalg.norm(P_m)
        it_s = (f"iters mean={its.mean():.1f} max={its.max()}"
                if its is not None else "")
        print(f"     {name:<8} relL2 vs analytic(Tv=.1)={e_an:.3e}  "
              f"vs monolithic={e_mo:.3e}  {it_s}")
    profA = {n: full_p(r[0][79])[line] / q for n, r in runs.items()
             if not r[2]}

    # ------------- Part A2: near-undrained early-time stability ------------
    invQ2 = poro / 2.2e12
    ops2 = ops_build(nodes, elems, gps, sch, E, nu, kbar, invQ2, a_stab,
                     ufix, top)
    n2 = 60
    dt2 = (0.01 * Ly * Ly / cv) / n2
    P2m = mono_march(ops2, None, 0, ff, dt2, n2)
    print("  A2 near-undrained stress test (Kf=2.2e12, 60 steps to Tv=0.01):")
    traces = {"mono": np.abs(P2m).max(axis=1) / q}
    for mode in ("fs", "fs1", "drained"):
        Ph, its, div = stag_march(ops2, None, 0, ff, dt2, n2, mode, lfac, q)
        traces[mode] = np.nanmax(np.abs(Ph[:len(its)]), axis=1) / q
        e_mo = (np.linalg.norm(Ph - P2m) / np.linalg.norm(P2m)
                if not div else np.inf)
        print(f"     {mode:<8} {'DIVERGED step ' + str(len(its)) if div else ''}"
              f"  vs monolithic={e_mo:.3e}  iters mean={its.mean():.1f} "
              f"max={its.max()}")

    # ------------- Part B: the E4 crack capability through the seam --------
    nxB = nyB = 10
    LB = 10.0
    nodesB, elemsB = mesh_rect(nxB, nyB, LB, LB)
    gpsB = gp_data(nodesB, elemsB)
    EB, nuB = 2.5e7, 0.3
    qB = 1.0e4
    KdrB = EB / (3.0 * (1.0 - 2.0 * nuB))
    GB = EB / (2.0 * (1.0 + nuB))
    lfacB = 1.0 / (KdrB + 4.0 * GB / 3.0)   # oedometric (see Part A note)
    EoedB = EB * (1 - nuB) / ((1 + nuB) * (1 - 2 * nuB))
    cvB = (kh / gw) * EoedB
    tv_end = 1.5
    nB = 450
    dtB = (tv_end * LB * LB / cvB) / nB
    step_rm = int(round(0.05 / tv_end * nB))
    crack = {4 * nxB + i for i in range(8)}
    surv = np.array([d["elem"] not in crack for d in gpsB])
    topB = np.nonzero(np.abs(nodesB[:, 1] - LB) < 1e-9)[0]
    botB = np.nonzero(np.abs(nodesB[:, 1]) < 1e-9)[0]
    sideB = np.nonzero((np.abs(nodesB[:, 0]) < 1e-9) |
                       (np.abs(nodesB[:, 0] - LB) < 1e-9))[0]
    ufixB = np.unique(np.concatenate(
        [[2 * m, 2 * m + 1] for m in botB] + [[2 * m] for m in sideB]))
    fB = np.zeros(2 * len(nodesB))
    for m in topB:
        w = 1.0 if 1e-9 < nodesB[m, 0] < LB - 1e-9 else 0.5
        fB[2 * m + 1] = -qB * w
    schB = build_scheme("FEM-stab", nodesB, elemsB, gpsB)
    a_stabB = mcgann_alpha(EB, nuB, 1.0)
    argsB = (nodesB, elemsB, gpsB, schB, EB, nuB, kbar, invQbar, a_stabB,
             ufixB, topB)
    opsB0 = ops_build(*argsB)
    opsDead = ops_build(*argsB, qmask=surv, fmask=surv)   # fluid dies
    opsPers = ops_build(*argsB, qmask=surv, fmask=None)   # fluid persists
    ffB = fB[opsB0["ufree"]]
    probe_id = 2 * (nxB + 1) + 2                          # node at (2, 2)
    ppos = int(np.nonzero(opsB0["pfree"] == probe_id)[0][0])
    P_dead = mono_march(opsB0, opsDead, step_rm, ffB, dtB, nB)
    P_ref = mono_march(opsB0, opsPers, step_rm, ffB, dtB, nB)
    P_stg, itsB, divB = stag_march(opsB0, opsPers, step_rm, ffB, dtB, nB,
                                   "fs", lfacB, qB)
    curves = {"FEM-mono (fluid dies)": P_dead[:, ppos],
              "monolithic fluid+ (ref)": P_ref[:, ppos],
              "STAGGERED fluid+ (seam)": P_stg[:, ppos]}
    print("  B  crack rerun through the seam (fluid-persistent overlay):")
    for name, cu in curves.items():
        i90 = np.nonzero(cu <= 0.1 * qB)[0]
        tv90 = tv_end * (i90[0] + 1) / nB if len(i90) else np.inf
        print(f"     {name:<26} p(end)/q={cu[-1] / qB:5.3f}   Tv90={tv90:5.3f}")
    devB = (np.linalg.norm(P_stg[:, ppos] - P_ref[:, ppos])
            / np.linalg.norm(P_ref[:, ppos]))
    print(f"     staggered vs monolithic fluid+ probe deviation: {devB:.2%}"
          f"   (iters mean={itsB.mean():.1f} max={itsB.max()})")

    # ------------- Part C: cost, 40x40, sparse solvers ---------------------
    print("  C  cost @ 40x40 (sparse LU; staggered sub-solves are SPD):")
    nxC = nyC = 40
    nodesC, elemsC = mesh_rect(nxC, nyC, 10.0, 10.0)
    gpsC = gp_data(nodesC, elemsC)
    schC = build_scheme("FEM-eq", nodesC, elemsC, gpsC)
    topC = np.nonzero(np.abs(nodesC[:, 1] - 10.0) < 1e-9)[0]
    botC = np.nonzero(np.abs(nodesC[:, 1]) < 1e-9)[0]
    sideC = np.nonzero((np.abs(nodesC[:, 0]) < 1e-9) |
                       (np.abs(nodesC[:, 0] - 10.0) < 1e-9))[0]
    ufixC = np.unique(np.concatenate(
        [[2 * m, 2 * m + 1] for m in botC] + [[2 * m] for m in sideC]))
    opsC = ops_build(nodesC, elemsC, gpsC, schC, EB, nuB, kbar, invQbar, 0.0,
                     ufixC, topC)
    fC = np.zeros(2 * len(nodesC))
    for m in topC:
        w = 0.25 if 1e-9 < nodesC[m, 0] < 10.0 - 1e-9 else 0.125
        fC[2 * m + 1] = -q * w
    ffC = fC[opsC["ufree"]]
    dtC = (0.5 * 100.0 / cvB) / 400
    n_u = opsC["K"].shape[0]

    def sparsify(M, tol=1e-14):
        m = np.abs(M).max()
        return sp.csc_matrix(np.where(np.abs(M) > tol * m, M, 0.0))

    A = np.zeros((n_u + opsC["S"].shape[0],) * 2)
    A[:n_u, :n_u] = opsC["K"]
    A[:n_u, n_u:] = -opsC["Q"]
    A[n_u:, :n_u] = opsC["Q"].T
    A[n_u:, n_u:] = opsC["S"] + dtC * opsC["H"]
    r = np.abs(A).max(axis=1)
    r[r == 0] = 1.0
    c = np.abs(A / r[:, None]).max(axis=0)
    c[c == 0] = 1.0
    A_s = sparsify(A / r[:, None] / c[None, :])
    K_s = sparsify(opsC["K"])
    F_s = sparsify(opsC["S"] + dtC * opsC["H"] + lfacB * opsC["Mp"])
    t0 = time.perf_counter(); luA = spla.splu(A_s)
    tfA = time.perf_counter() - t0
    t0 = time.perf_counter(); luK = spla.splu(K_s)
    tfK = time.perf_counter() - t0
    t0 = time.perf_counter(); luF = spla.splu(F_s)
    tfF = time.perf_counter() - t0
    Qs = sp.csc_matrix(opsC["Q"])
    Ss = sparsify(opsC["S"])
    # monolithic: 20 steps
    u = np.zeros(n_u); p = np.zeros(opsC["S"].shape[0])
    t0 = time.perf_counter()
    for step in range(20):
        rhs = np.concatenate([ffC, Qs.T @ u + Ss @ p]) / r
        z = luA.solve(rhs) / c
        u, p = z[:n_u], z[n_u:]
    t_mono = (time.perf_counter() - t0) / 20
    # staggered fixed-stress: 20 steps, measured iterations included
    Ls = lfacB * sp.csc_matrix(opsC["Mp"])
    u = np.zeros(n_u); p = np.zeros(opsC["S"].shape[0])
    nits = []
    t0 = time.perf_counter()
    for step in range(20):
        pk = p
        for nit in range(1, 201):
            u1 = luK.solve(ffC + Qs @ pk)
            p1 = luF.solve(Ss @ p + Ls @ pk - Qs.T @ (u1 - u))
            dcon = np.linalg.norm(p1 - pk) / (np.linalg.norm(p1) + q)
            pk = p1
            if dcon < 1e-6:
                break
        u1 = luK.solve(ffC + Qs @ p1)
        u, p = u1, p1
        nits.append(nit)
    t_stag = (time.perf_counter() - t0) / 20
    nbar = np.mean(nits)
    print(f"     factor:  coupled-unsym {tfA * 1e3:7.1f}ms   "
          f"K-solid {tfK * 1e3:6.1f}ms + F-fluid {tfF * 1e3:5.1f}ms "
          f"= {(tfK + tfF) * 1e3:6.1f}ms  (x{tfA / (tfK + tfF):.2f} cheaper)")
    print(f"     per step: monolithic {t_mono * 1e3:6.2f}ms   staggered "
          f"{t_stag * 1e3:6.2f}ms at iters mean={nbar:.1f} "
          f"(x{t_mono / t_stag:.2f} {'faster' if t_stag < t_mono else 'slower'})")
    print("     note: splu does not exploit SPD; Cholesky on K/F would gain "
          "~2x more, and OpenSees gets ProfileSPD/MUMPS SYM=2 back.")

    # -------------------------------- plot ---------------------------------
    fig = plt.figure(figsize=(13, 9))
    ax = fig.add_subplot(2, 2, 1)
    Zf = np.linspace(0, 1, 200)
    ax.plot(terzaghi_series(Zf, 0.1), Zf, "k-", lw=2, label="analytic")
    for name, pr in profA.items():
        ax.plot(pr, Z, marker="o", ms=3, lw=1, label=name)
    ax.invert_yaxis()
    ax.set_xlabel("p/q")
    ax.set_ylabel("depth/H")
    ax.set_title("A: Terzaghi Tv=0.1 -- splits vs monolithic", fontsize=10)
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
    ax = fig.add_subplot(2, 2, 2)
    for name, tr in traces.items():
        ax.semilogy(np.arange(1, len(tr) + 1), np.clip(tr, 1e-6, 1e30),
                    lw=1.5, label=name)
    ax.set_xlabel("step")
    ax.set_ylabel("max |p| / q")
    ax.set_title("A2: near-undrained stability (Kf=2.2e12)", fontsize=10)
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
    ax = fig.add_subplot(2, 1, 2)
    tvs = tv_end * np.arange(1, nB + 1) / nB
    for name, cu in curves.items():
        ax.plot(tvs, cu / qB, lw=2, label=name)
    ax.axvline(tv_end * step_rm / nB, color="k", ls="--", lw=0.8)
    ax.set_xlabel("Tv")
    ax.set_ylabel("p(probe)/q")
    ax.set_title("B: crack capability through the staggered seam", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{HERE}/staggered_seam.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    import sys
    np.set_printoptions(precision=3, suppress=True)
    only = sys.argv[1] if len(sys.argv) > 1 else "all"
    if only in ("all", "e1"):
        experiment_infsup()
    if only in ("all", "e2"):
        experiment_terzaghi()
    if only in ("all", "e3"):
        experiment_footing()
    if only in ("all", "e4"):
        experiment_removal()
    if only in ("all", "e5"):
        experiment_perf()
    if only in ("all", "e6"):
        experiment_staggered()
    print("=" * 78)
    print("done -- plots: infsup_spectra.png, terzaghi_profiles.png, "
          "footing_maps.png, removal_flowpath.png, perf_scaling.png, "
          "staggered_seam.png")
