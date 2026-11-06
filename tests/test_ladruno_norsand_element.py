"""WP-144 gate G2 (Zone A): `nDMaterial LadrunoNorSand` THROUGH REAL ELEMENTS.

Scope.  The material-level shell contract is `tests/test_ladruno_norsand.py`; the maths
(return map, tangent vs the oracles, K1/K2) is Zone B.  This file is what neither of those
can see: the material inside a real element, with FREE equations, assembled by the real
solver, driven by a real Newton iteration.  Nothing here imports an oracle, scipy or sympy
(numpy only).  Four gates:

  (1) FD TANGENT THROUGH AN ELEMENT.  Assembled element stiffness (`printA`, the integrator's
      own `formTangent`, FullGeneral) against the central finite difference of the element
      resisting force (the SOE `b` vector after `ladrunoTrialResidualNorm`, i.e. the integrator's
      own `formUnbalance`) with respect to every free nodal displacement, h in {1e-7, 1e-8} x
      size, best h, per-column relative error <= 1e-6 (plan sec.5.1 bound: truncation scale
      (h/kappa_hat)^2/6 with kappa_hat = 0.01, not an observed number).  One stdBrick (3D) and
      one `quad ... PlaneStrain` (the PlaneStrain wrapper), 6 + 5 rigid-body-removed free
      equations each, at plastic states of drained TXC and of a NON-COAXIAL shear state, each
      both on a single-substep increment and on a SUBSTEPPED one (response `substeps`).
  (2) QUADRATIC NEWTON on a 2x2x2 stdBrick drained triaxial (jittered interior nodes, UmfPack,
      Newton), steps that substep included: order of the asymptotic iterations >= 1.8.
  (3) PLANE-STRAIN WRAPPER: sigma_33 is REPORTED (`stressesPlaneStrain`, 4th component per
      Gauss point) and agrees with the 3D material on the same in-plane strain path.
  (4) MASS DENSITY: `-density` enters the element mass and `-rho` (the ELLIPTICITY) does not.

WHERE EVERY EXPECTED VALUE COMES FROM (none is harvested from the shell's own output).
  * (1) The expected tangent IS the finite difference of the element's own force: the contract
    "K = dF/du of the one-step map from the committed state" is the definition of a consistent
    tangent; the tolerance is the plan sec.5.1 / G1 bound (1e-6, the oracle's own FD gate).  The
    premise of every FD point (plastic at all 8/4 Gauss points, not refused, the SAME substep
    count as the base increment: the ladder is a decision, not differentiated, sheet sec.9.6)
    is asserted.  Two CONTROLS run the same gate on deliberately wrong matrices built from the
    same data (the symmetrised tangent; the tangent of the committed, elastic state) and must
    FAIL it: they prove the gate has teeth at these states.
  * (2) The order estimate p_k = ln(e_{k+1}/e_k) / ln(e_k/e_{k-1}) of a Newton iteration with an
    exact tangent is 2 (the textbook theorem); 1.8 is the G1 convergence-test bound.  A CONTROL
    with `ModifiedNewton` (a stale tangent) must come out at ~1.
  * (3) Elastic regime: sheet (S.5) with alpha0 = 0, p = p0 exp(-(eps_v - eps_v0)/kappa_hat),
    s = 2 mu0 e, hence for eps_33 = 0:  sigma_33 = p - (2/3) mu0 eps_v,
    sigma_ii = p + 2 mu0 (eps_ii - eps_v/3), sigma_12 = mu0 gamma_12 (the sheet's own closed
    form, evaluated here in plain python).  Plastic regime: the PlaneStrain lane IS the 3D
    material restricted to eps_33 = gamma_13 = gamma_23 = 0 (class doc of LadrunoNorSand), so the
    3D brick on the same prescribed in-plane strain path is the reference, to round-off.
  * (4) Brick/quad mass: e_x' (K + c3 M) e_x = c3 * rho_mass * Volume (the rigid translation
    e_x is annihilated by the stiffness; c3 = 1/(beta dt^2) is Newmark's), a closed form.

SIGN/LAYOUT FACTS the harness relies on (each was measured, then written down):
  * `printA('-ret')` returns the matrix COLUMN-major: reshape(N, N).T is K.
  * `printA` re-forms the tangent itself (`formTangent`) from the element's CURRENT trial
    state; stdBrick/quad recompute from the nodal TRIAL displacements, so the nodal trial state
    must be injected (`ladrunoSetNodeTrial`, full vector) and the elements updated
    (`ladrunoTrialResidualNorm`) BEFORE `printA`, or K is the tangent of the previous FD point.
  * `setNodeDisp` per dof restarts from the COMMITTED vector (does not accumulate): unusable.

WALL TIME (measured on the dev box, Windows, CPython 3.12, 2026-10-01, whole file `--runslow`:
82.7 s for 18 tests).  Default tier (no `--runslow`) ~31 s.  Three cases are @pytest.mark.slow, each
because the PHYSICS is expensive (a substepped increment costs a 2^k ladder with nested pi_i scans at
every Gauss point, and the FD evaluates the force 2 x 18 x 2 times): brick `txc_4` 21 s, the perturbed
non-uniform brick 14 s, the 4-fold-substep Newton step 17 s.  The default tier keeps one m = 4 FD
case per element family (brick `shear_4` 13.6 s, quad `txc_4` 2.3 s).

Every test names, in a `KILLS:` line, the mutant that would make it fail.
"""
import math
import os

import numpy as np
import pytest

# ---- engine binding (same idiom as tests/test_ladruno_norsand.py) -----------
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
_DIST_BIN = os.path.join(_ROOT, "dist", "bin")
if os.path.isfile(os.path.join(_DIST_BIN, "opensees.pyd")):
    from _engine import bind_worktree_engine
    ops = bind_worktree_engine(_DIST_BIN)
else:
    from _testbed import ops

pytestmark = [pytest.mark.zone_a]


# ===========================================================================
#  K2 parameter set (AB06 sec.6.1, published) -- same base deck as the shell file
# ===========================================================================
_P0, _KH, _MU = -100.0, 0.01, 5400.0
_M, _N, _NB, _RHO, _RHOB, _CHI, _H = 1.2, 0.4, 0.2, 0.7, 0.8, -3.5, 280.0
_LT, _VC0, _V0, _PI0 = 0.0135, 1.81, 1.59, -60.4


def _mat_args(**over):
    d = {"p0": _P0, "kappa_hat": _KH, "mu0": _MU, "M": _M, "N": _N, "N_bar": _NB,
         "rho": _RHO, "rho_bar": _RHOB, "chi": _CHI, "h": _H, "csl": "paper",
         "lambda_tilde": _LT, "v_c0": _VC0, "v0": _V0, "pi0": _PI0}
    d.update(over)
    out = []
    for k, v in d.items():
        if v is None:
            continue
        out += ["-" + k, v]
    return out


def _make_mat(tag=1, **over):
    ops.nDMaterial("LadrunoNorSand", tag, *_mat_args(**over))


# ===========================================================================
#  Decks.  Drained triaxial by LOAD control: tractions -100 kPa on all faces
#  (lumped to the nodes, constant), plus an axial ramp Q (and a shear ramp T)
#  on the z faces.  Rigid-body motion is removed with the MINIMAL 6 constraints
#  (node 1 xyz, node 2 yz, node 4 z), so every other nodal DOF is a FREE equation
#  and every Gauss point sees the uniform-strain solution.
# ===========================================================================
_BRICK_XYZ = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
_XF1, _XF0 = [2, 3, 6, 7], [1, 4, 5, 8]
_YF1, _YF0 = [3, 4, 7, 8], [1, 2, 5, 6]
_ZF1, _ZF0 = [5, 6, 7, 8], [1, 2, 3, 4]


def _solver_stack(test_tol=1e-12, test_iter=40, system="FullGeneral"):
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system(system)                       # the tangent is NON-symmetric: never a symmetric SOE
    ops.test("NormDispIncr", test_tol, test_iter, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.1)
    ops.analysis("Static")


def _deck_brick(Q, T):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _make_mat(1)
    for i, c in enumerate(_BRICK_XYZ, 1):
        ops.node(i, float(c[0]), float(c[1]), float(c[2]))
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 0, 1, 1)
    ops.fix(4, 0, 0, 1)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Constant", 2)
    ops.pattern("Plain", 2, 2)
    for n in _XF1:
        ops.load(n, -25.0, 0.0, 0.0)
    for n in _XF0:
        ops.load(n, +25.0, 0.0, 0.0)
    for n in _YF1:
        ops.load(n, 0.0, -25.0, 0.0)
    for n in _YF0:
        ops.load(n, 0.0, +25.0, 0.0)
    for n in _ZF1:
        ops.load(n, 0.0, 0.0, -25.0)
    for n in _ZF0:
        ops.load(n, 0.0, 0.0, +25.0)
    ops.timeSeries("Linear", 3)
    ops.pattern("Plain", 3, 3)
    for n in _ZF1:
        ops.load(n, 0.0, 0.0, -Q / 4.0)
    for n in _ZF0:
        ops.load(n, 0.0, 0.0, +Q / 4.0)
    if T:          # tau_xz couple: x-forces on the z faces, z-forces on the x faces (moment-balanced)
        for n in _ZF1:
            ops.load(n, T / 4.0, 0.0, 0.0)
        for n in _ZF0:
            ops.load(n, -T / 4.0, 0.0, 0.0)
        for n in _XF1:
            ops.load(n, 0.0, 0.0, T / 4.0)
        for n in _XF0:
            ops.load(n, 0.0, 0.0, -T / 4.0)
    _solver_stack()
    return {"nodes": list(range(1, 9)), "ndf": 3, "eles": [1], "ngp": 8}


def _deck_quad(Q, T):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    _make_mat(1)
    for i, c in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], 1):
        ops.node(i, float(c[0]), float(c[1]))
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStrain", 1)
    ops.timeSeries("Constant", 2)
    ops.pattern("Plain", 2, 2)
    for n in (2, 3):
        ops.load(n, -50.0, 0.0)
    for n in (1, 4):
        ops.load(n, +50.0, 0.0)
    for n in (3, 4):
        ops.load(n, 0.0, -50.0)
    for n in (1, 2):
        ops.load(n, 0.0, +50.0)
    ops.timeSeries("Linear", 3)
    ops.pattern("Plain", 3, 3)
    for n in (3, 4):
        ops.load(n, 0.0, -Q / 2.0)
    for n in (1, 2):
        ops.load(n, 0.0, +Q / 2.0)
    if T:
        for n in (3, 4):
            ops.load(n, T / 2.0, 0.0)
        for n in (1, 2):
            ops.load(n, -T / 2.0, 0.0)
        for n in (2, 3):
            ops.load(n, 0.0, T / 2.0)
        for n in (1, 4):
            ops.load(n, 0.0, -T / 2.0)
    _solver_stack()
    return {"nodes": [1, 2, 3, 4], "ndf": 2, "eles": [1], "ngp": 4}


def _start(deck):
    """Zero-factor first step: builds the numbering (the state is in equilibrium at lambda 0)."""
    ops.integrator("LoadControl", 0.0)
    assert ops.analyze(1) == 0


def _steps(deck, steps):
    for s in steps:
        ops.integrator("LoadControl", float(s))
        assert ops.analyze(1) == 0, f"load step {s} did not converge"


def _nodal_disp(deck):
    return np.array([ops.nodeDisp(n) for n in deck["nodes"]], dtype=float)


def _free_dofs(deck):
    out = []
    for n in deck["nodes"]:
        for d, eq in enumerate(ops.nodeDOFs(n)):
            if eq >= 0:
                out.append((eq, n, d))
    return sorted(out)


def _inject(deck, U):
    z = [0.0] * deck["ndf"]
    for n, row in zip(deck["nodes"], U):
        ops.ladrunoSetNodeTrial(n, *[float(x) for x in row], *z, *z)


def _force(deck, U):
    """Free-DOF element resisting force at the injected TRIAL state (b = P - F, P constant)."""
    _inject(deck, U)
    ops.ladrunoTrialResidualNorm()
    return -np.array(ops.printB("-ret"), dtype=float)


def _gp_info(deck):
    """[refusal, plastic, vertex, cap, local_iters, pi_iters, substeps, finest, finest_sub] per Gauss point."""
    return np.array([ops.eleResponse(e, "material", g, "stepInfo")
                     for e in deck["eles"] for g in range(1, deck["ngp"] + 1)], dtype=float)


def _col_rel_err(K, Kfd):
    cn = np.linalg.norm(Kfd, axis=0)
    return float(np.max(np.linalg.norm(K - Kfd, axis=0) / cn))


_FD_H = (1e-7, 1e-8)          # x size (unit element edge)
_FD_TOL = 1e-6                # plan sec.5.1 / G1: best-h per-column relative bound


def _fd_case(make_deck, Q, T, hist, big, m_expected, pert=0.0, seed=7):
    """Run the load history to the committed state S_{n-1}, then differentiate the ONE-STEP map from
    S_{n-1} at the trial displacement of a `big` load step (obtained on an identical first run).

    Returns a dict with the per-h errors, the controls, the premise flags and K."""
    # ---- run A: the trial displacement of the big step (Newton-converged, equilibrium) ----
    deck = make_deck(Q, T)
    _start(deck)
    _steps(deck, hist)
    _steps(deck, [big])
    U = _nodal_disp(deck)
    # ---- run B: back at S_{n-1}, committed ----
    deck = make_deck(Q, T)
    _start(deck)
    _steps(deck, hist)
    fd = _free_dofs(deck)
    N = len(fd)
    free_mask = np.zeros_like(U, dtype=bool)
    for _, n, d in fd:
        free_mask[deck["nodes"].index(n), d] = True
    if pert:
        rng = np.random.default_rng(seed)
        U = U + pert * rng.standard_normal(U.shape) * free_mask    # non-uniform GP strains, free DOFs only

    def vec_of(Um):
        return np.array([Um[deck["nodes"].index(n), d] for _, n, d in fd])

    def with_vec(v):
        Um = U.copy()
        for x, (_, n, d) in zip(v, fd):
            Um[deck["nodes"].index(n), d] = x
        return Um

    u0 = vec_of(U)
    # ---- the base point: premise + K from the real integrator ----
    _force(deck, U)                                   # updates the elements to the trial state
    base = _gp_info(deck)
    assert np.all(base[:, 0] == 0), f"base increment refused at a Gauss point: {base[:, 0]}"
    assert np.all(base[:, 1] == 1), f"base increment not plastic at every Gauss point: {base[:, 1]}"
    assert np.all(base[:, 6] == m_expected), (
        f"premise: substeps at the base increment {base[:, 6]} != {m_expected} at every Gauss point")
    K = np.array(ops.printA("-ret"), dtype=float).reshape(N, N).T
    asym = float(np.max(np.abs(K - K.T)) / np.max(np.abs(K)))
    # ---- the tangent of the COMMITTED state S_{n-1} (a zero increment: elastic) as a control ----
    Uc = _nodal_disp(deck)           # run B's committed nodal displacement, U_{n-1}
    _force(deck, Uc)
    Kc = np.array(ops.printA("-ret"), dtype=float).reshape(N, N).T
    # ---- central FD of the force, per h; an h counts only if every FD point keeps the premise ----
    errs, valid = {}, {}
    for h in _FD_H:
        Kfd = np.zeros((N, N))
        ok = True
        for j in range(N):
            e = np.zeros(N)
            e[j] = h
            fp = _force(deck, with_vec(u0 + e))
            ip = _gp_info(deck)
            fm = _force(deck, with_vec(u0 - e))
            im = _gp_info(deck)
            for info in (ip, im):
                ok = ok and bool(np.all(info[:, 0] == 0) and np.all(info[:, 1] == 1)
                                 and np.all(info[:, 6] == m_expected))
            Kfd[:, j] = (fp - fm) / (2.0 * h)
        errs[h] = _col_rel_err(K, Kfd)
        valid[h] = ok
        if h == _FD_H[-1]:
            Kfd_last = Kfd
    return {"errs": errs, "valid": valid, "asym": asym, "N": N, "K": K, "Kc": Kc, "Kfd": Kfd_last,
            "err_sym": _col_rel_err(0.5 * (K + K.T), Kfd_last),
            "err_committed": _col_rel_err(Kc, Kfd_last)}


def _assert_fd(res, label):
    good = [h for h in _FD_H if res["valid"][h]]
    assert good, f"{label}: no FD step size kept the premise (plastic, same substeps) at every point: {res['valid']}"
    best = min(res["errs"][h] for h in good)
    assert best <= _FD_TOL, (
        f"{label}: assembled element stiffness vs central FD of the resisting force: best-h per-column "
        f"relative error {best:.3e} > {_FD_TOL} (per h: {res['errs']}); the tangent is not the "
        f"derivative of the one-step map")
    # CONTROLS: the same gate must reject wrong matrices built from the same data.
    assert res["asym"] > 1e-2, (
        f"{label}: premise: the tangent is expected to be clearly non-symmetric here, asym {res['asym']:.3e}")
    assert res["err_sym"] > 1e-3, (
        f"{label}: the symmetrised tangent must FAIL the FD gate (err {res['err_sym']:.3e}), else the gate "
        f"cannot see a symmetrisation mutant at this state")
    assert res["err_committed"] > 1e-2, (
        f"{label}: the committed-state (elastic) tangent must FAIL the FD gate (err {res['err_committed']:.3e}), "
        f"else the state is not plastic enough to discriminate")
    return best


# ---- (1a) stdBrick -----------------------------------------------------------------------------
# (hist, big, Q, T, substeps): load factors; Q = axial ramp (kPa at lambda = 1), T = shear ramp.
#  txc_1   drained TXC to q = 160 kPa in one 100 kPa step from q = 60 (plastic since q ~ 40): 1 substep.
#          The state is on the Willam-Warnke compression corner (sheet sec.4.3: FD is O(h) there:
#          measured 5.2e-6 / 5.2e-7 / 5.2e-8 at h = 1e-6 / 1e-7 / 1e-8), hence "best h".
#  txc_4   drained TXC to q = 300 kPa from the elastic start in ONE step: the whole increment is
#          refused and the kernel ladder takes 4 sub-increments (measured, asserted as the premise).
#  shear_1 TXC + tau_xz = 30 kPa: off the corner, non-coaxial (principal axes rotate off x,y,z).
#  shear_4 same family, one step from the elastic start, 4 sub-increments.
_BRICK_CASES = {
    "txc_1":   dict(Q=100.0, T=0.0,  hist=[0.1] * 6, big=1.0, m=1),
    "txc_4":   dict(Q=100.0, T=0.0,  hist=[],         big=3.0, m=4),
    "shear_1": dict(Q=100.0, T=30.0, hist=[0.1] * 6, big=1.0, m=1),
    "shear_4": dict(Q=100.0, T=20.0, hist=[],         big=2.5, m=4),
}


@pytest.mark.parametrize("case", [
    "txc_1",
    pytest.param("txc_4", marks=pytest.mark.slow),     # 22.6 s measured
    "shear_1",
    "shear_4",                                         # 15 s measured; the default-tier m = 4 gate
])
def test_brick_assembled_tangent_vs_fd(case):
    """Assembled stdBrick stiffness (free equations, FullGeneral) == central FD of the resisting force.

    Measured at write time (best h, per-column relative error / asymmetry of K / error of the
    symmetrised K): txc_1 6.5e-8 (the WW compression corner: O(h), 6.6e-7 at h = 1e-7) / 4.9 % / 5e-2;
    txc_4 7.7e-9 / 23.5 % / 0.26; shear_1 9.1e-11 / 11 % / 0.10; shear_4 2.2e-10 / 22 % / 0.28.  The
    committed-state (elastic) tangent is off by 0.3-1.3 at these states.  Wall times: txc_1 0.8 s,
    shear_1 0.9 s, shear_4 13.6 s, txc_4 21 s (slow tier).

    KILLS: (a) "returns the last-substep tangent" (sheet sec.9.6: 0.51 at m = 2, 0.81-0.90 at m = 4
    relative error vs this very FD): the *_4 cases; (b) a symmetrised / symmetric-assuming tangent
    (control assertion + the FD itself); (c) a wrong Voigt shear-column factor or a dropped tangent term
    in the shell (`getTangent` halves the engineering-shear columns): the shear_* cases have shear
    columns that are O(1) of the stiffness; (d) a tangent frozen at the committed/elastic state
    (control assertion); (e) vfac = v0 or v_n instead of the converged v_{n+1} in the s_k term of (S.31) (G2 owner
    decision 2026-10-01, exponential v-update, dv_{n+1}/d eps~ = v_{n+1}; sheet sec.1.2 FD-measured
    1.8e-5 (v0) and 1.7e-6 (v_n) on a~ at the K2 state, 1.1e-4 / 3.8e-7 at v = 1.45: resolved by this
    FD at the states where v has moved away from v0; the Mutate step measures it through the element)."""
    c = _BRICK_CASES[case]
    res = _fd_case(_deck_brick, c["Q"], c["T"], c["hist"], c["big"], c["m"])
    assert res["N"] == 18, "premise: 24 nodal DOFs minus the 6 rigid-body constraints = 18 free equations"
    _assert_fd(res, f"stdBrick/{case}")


@pytest.mark.slow
def test_brick_assembled_tangent_vs_fd_perturbed_nonuniform():
    """The same gate with a deterministic non-uniform perturbation (2e-4 rms) of the free nodal
    displacements on top of the substepped non-coaxial solution: the eight Gauss points now carry
    DIFFERENT strains (no shared principal axes, none of them on a corner), so the per-Gauss-point
    tangents differ and are summed by the real quadrature.  Measured: best-h error 2.9e-10, K asymmetry
    22 %.  SLOW TIER, 13.6 s.

    KILLS: a tangent that is right only for a uniform strain field or only at a shared Gauss-point
    state (e.g. a static scratch buffer in the shell shared between Gauss-point clones)."""
    res = _fd_case(_deck_brick, 100.0, 20.0, [], 2.5, 4, pert=2e-4)
    _assert_fd(res, "stdBrick/shear_4+perturbation")


# ---- (1b) quad + PlaneStrain wrapper ------------------------------------------------------------
_QUAD_CASES = {
    "txc_1":   dict(Q=100.0, T=0.0,  hist=[0.1] * 6, big=1.0, m=1),
    "txc_4":   dict(Q=100.0, T=0.0,  hist=[],         big=3.0, m=4),
    "shear_1": dict(Q=100.0, T=20.0, hist=[0.1] * 6, big=1.0, m=1),
    "shear_2": dict(Q=100.0, T=20.0, hist=[],         big=2.0, m=2),
}


@pytest.mark.parametrize("case", list(_QUAD_CASES))
def test_planestrain_quad_assembled_tangent_vs_fd(case):
    """Assembled `quad ... PlaneStrain` stiffness (5 free equations) == central FD of the resisting
    force, through the PlaneStrain wrapper (`getCopy("PlaneStrain")`).  Measured best-h error 4e-11 ..
    3e-10 in all four cases (K asymmetry 2.3-12 %, symmetrised-K error 1.7e-2 .. 0.10).  0.1-2.3 s each.

    KILLS: the same list as the brick gate, applied to the wrapper's own 3x3 reduction of the
    6x6 tangent (a wrong row/column pick of the reduced Voigt order, a shear factor taken from the
    wrong block, sigma_33 coupling dropped from the reduction)."""
    c = _QUAD_CASES[case]
    res = _fd_case(_deck_quad, c["Q"], c["T"], c["hist"], c["big"], c["m"])
    assert res["N"] == 5, "premise: 8 nodal DOFs minus 3 rigid-body constraints = 5 free equations"
    _assert_fd(res, f"quad-PlaneStrain/{case}")


# ===========================================================================
#  (2) Quadratic Newton convergence on a 2x2x2 stdBrick drained triaxial
# ===========================================================================
def _deck_mesh222(test_kind, tol, alg, jitter=0.06):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _make_mat(1)
    g = (0.0, 0.5, 1.0)
    tag, n = {}, 0
    for k in range(3):
        for j in range(3):
            for i in range(3):
                n += 1
                tag[(i, j, k)] = n
                x = [g[i], g[j], g[k]]
                rng = np.random.default_rng(100 + n)
                d = rng.uniform(-jitter, jitter, 3)
                for a, ix in enumerate((i, j, k)):
                    if ix == 1:                       # only nodes interior in that direction move
                        x[a] += d[a]
                ops.node(n, *x)
    for (i, j, k), nt in tag.items():                 # rollers on the three coordinate planes
        ops.fix(nt, int(i == 0), int(j == 0), int(k == 0))
    e = 0
    for k in range(2):
        for j in range(2):
            for i in range(2):
                e += 1
                ops.element("stdBrick", e, tag[(i, j, k)], tag[(i + 1, j, k)], tag[(i + 1, j + 1, k)],
                            tag[(i, j + 1, k)], tag[(i, j, k + 1)], tag[(i + 1, j, k + 1)],
                            tag[(i + 1, j + 1, k + 1)], tag[(i, j + 1, k + 1)], 1)

    def w(a, b):          # tributary share of a face node of a 2x2 face of unit area
        return 0.25 if (a == 1 and b == 1) else (0.125 if (a == 1 or b == 1) else 0.0625)

    ops.timeSeries("Constant", 2)
    ops.pattern("Plain", 2, 2)
    for (i, j, k), nt in tag.items():                 # confinement 100 kPa on x = 1, y = 1 and z = 1
        if i == 2:
            ops.load(nt, -100.0 * w(j, k), 0.0, 0.0)
        if j == 2:
            ops.load(nt, 0.0, -100.0 * w(i, k), 0.0)
        if k == 2:
            ops.load(nt, 0.0, 0.0, -100.0 * w(i, j))
    ops.timeSeries("Linear", 3)
    ops.pattern("Plain", 3, 3)
    for (i, j, k), nt in tag.items():                 # axial ramp, 100 kPa per unit load factor
        if k == 2:
            ops.load(nt, 0.0, 0.0, -100.0 * w(i, j))
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("UmfPack")                             # unsymmetric sparse LU
    ops.test(test_kind, tol, 60, 0)
    ops.algorithm(alg)
    ops.integrator("LoadControl", 0.1)
    ops.analysis("Static")


def _mesh_substeps():
    s = [ops.eleResponse(e, "material", g, "substeps") for e in range(1, 9) for g in range(1, 9)]
    return max(x[0] for x in s), max(x[1] for x in s)


def _asymptotic_orders(norms, floor, rel_top=1e-2):
    """Order estimates p_k = ln(e_{k+1}/e_k)/ln(e_k/e_{k-1}) over the ASYMPTOTIC triples of one
    step: e_{k-1} <= rel_top * e_0 (the pre-asymptotic iterations are excluded), e_{k+1} >= floor
    (above the round-off floor of the norm), strictly decreasing."""
    out = []
    for k in range(1, len(norms) - 1):
        a, b, c = norms[k - 1], norms[k], norms[k + 1]
        if a <= rel_top * norms[0] and c >= floor and a > b > c > 0.0:
            out.append(math.log(c / b) / math.log(b / a))
    return out


def _newton_run(steps, test_kind="NormDispIncr", alg="Newton"):
    tol, floor = (1e-12, 1e-14) if test_kind == "NormDispIncr" else (1e-10, 1e-12)
    _deck_mesh222(test_kind, tol, alg)
    ops.integrator("LoadControl", 0.0)
    assert ops.analyze(1) == 0
    out = []
    for dl in steps:
        ops.integrator("LoadControl", float(dl))
        rc = ops.analyze(1)
        norms = list(ops.testNorm()[: ops.testIter()])
        mx, nsub = _mesh_substeps()
        out.append({"dl": dl, "rc": rc, "norms": norms, "iters": len(norms), "max_substeps": mx,
                    "n_substepped_commits": nsub, "orders": _asymptotic_orders(norms, floor)})
        if rc != 0:
            break
    return out


_ORDER_MIN = 1.8


@pytest.mark.parametrize("test_kind", ["NormDispIncr", "NormUnbalance"])
def test_newton_is_quadratic_through_substepping(test_kind):
    """2x2x2 stdBrick, jittered interior nodes (so the stress field is genuinely non-uniform),
    UmfPack + Newton, load factors 0.5 (onset of yield), 1.5 and 1.0 (cumulative q = 300 kPa).  The
    1.5 step makes the material SUBSTEP (the kernel ladder, response `substeps` > 1 at some Gauss
    point); the asymptotic iterations of EVERY step must have order >= 1.8, with the measure
    NormDispIncr 1e-12 (increment norms) and, independently, NormUnbalance 1e-10 (residual force
    norms).  At least one substepped step and >= 2 asymptotic triples must be observed (the
    non-vacuity guards).  Measured: 2 asymptotic triples per measure, orders 1.98 and 2.00 (NormDispIncr),
    1.96 and 2.00 (NormUnbalance); ~6 s per measure (the substepped step dominates).

    KILLS: "last-substep tangent" (error 0.5-0.9 at m = 4 per sheet sec.9.6: linear convergence at
    best), "symmetrised tangent" (asymmetry 5-22 % measured above: linear, rate ~0.1), a tangent
    taken at the committed instead of the trial state, a dropped tangent term (all of these are
    linear-rate or divergent in Newton)."""
    res = _newton_run([0.5, 1.5, 1.0], test_kind)
    assert all(r["rc"] == 0 for r in res), f"Newton did not converge: {[(r['dl'], r['rc']) for r in res]}"
    assert max(r["max_substeps"] for r in res) >= 2, "premise: the material must substep in this run"
    assert res[-1]["n_substepped_commits"] >= 1
    orders = [p for r in res for p in r["orders"]]
    assert len(orders) >= 2, f"non-vacuity: only {len(orders)} asymptotic triples ({[r['norms'] for r in res]})"
    assert min(orders) >= _ORDER_MIN, (
        f"Newton is not quadratic ({test_kind}): orders per step "
        f"{[(r['dl'], [round(p, 2) for p in r['orders']]) for r in res]}, norms {[r['norms'] for r in res]}")


@pytest.mark.slow
def test_newton_is_quadratic_with_four_substeps_in_one_step():
    """Slow tier: one load factor 3.0 step (q = 300 kPa from the elastic start) on the 2x2x2 mesh; the
    converged increment is integrated with 4 sub-increments (response `substeps`, asserted) and Newton
    still converges quadratically (NormDispIncr; measured order 2.00).  MEASURED WALL TIME 17 s (11 Newton
    iterations; every refused whole-increment trial costs a 2^k substep ladder with the nested pi_i
    scans).  The load factor must stay <= ~3.0: the jittered mesh reaches its limit load near 3.1-3.2
    and a step past it burns minutes in refused ladders.

    KILLS: the m = 4 form of "last-substep tangent" (0.81-0.90 error, sheet sec.9.6), which the m = 2
    step in the default-tier test above can miss (a ladder whose first half is elastic has T_1 = I and
    the last-substep tangent coincides with the chain there, sheet sec.9.6 (B))."""
    res = _newton_run([3.0])
    r = res[0]
    assert r["rc"] == 0, r
    assert r["max_substeps"] >= 4, f"premise: expected a 4-fold substepped increment, got {r['max_substeps']}"
    assert r["orders"], f"non-vacuity: no asymptotic triple in {r['norms']}"
    assert min(r["orders"]) >= _ORDER_MIN, (r["orders"], r["norms"])


def test_order_estimator_sees_a_stale_tangent():
    """CONTROL for the gate above: the same mesh with `ModifiedNewton` (the tangent of the first
    iteration of the step is reused: a stale tangent, exactly what the tangent mutants amount to)
    must come out at LINEAR order, estimator ~1 (measured 0.99-1.00 over 10 triples; gate 0.8-1.5).
    Without this the >= 1.8 gate could be passing for a reason other than the tangent.  0.2 s.

    KILLS: an order estimator or a convergence harness that cannot tell quadratic from linear."""
    res = _newton_run([0.5], "NormDispIncr", "ModifiedNewton")
    r = res[0]
    assert r["rc"] == 0
    assert len(r["orders"]) >= 5, r["norms"]
    assert 0.8 < float(np.median(r["orders"])) < 1.5, r["orders"]


# ===========================================================================
#  (3) PlaneStrain wrapper: sigma_33 is reported and consistent with the 3D material
# ===========================================================================
def _prescribed_path(a_final, nsteps, e11_fac, g12_fac):
    """Strain path (e11, e22, g12) = a * (e11_fac, -1, g12_fac), a ramped 0 -> a_final in nsteps."""
    return [(a_final * (k / nsteps) * e11_fac, -a_final * (k / nsteps), a_final * (k / nsteps) * g12_fac)
            for k in range(1, nsteps + 1)]


def _deck_ps_driver(kind, e_final):
    """Zero-free-DOF unit cell, every nodal DOF prescribed: u = H x with H = [[e11, g12/2], [g12/2, e22]]
    (symmetric: no rigid rotation), in-plane only.  'quad': PlaneStrain quad; 'brick': stdBrick with the z
    displacement fixed everywhere (eps_33 = gamma_13 = gamma_23 = 0)."""
    e11, e22, g12 = e_final
    H = np.array([[e11, 0.5 * g12], [0.5 * g12, e22]])
    ops.wipe()
    if kind == "quad":
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        _make_mat(1)
        xy = [(0, 0), (1, 0), (1, 1), (0, 1)]
        for i, c in enumerate(xy, 1):
            ops.node(i, float(c[0]), float(c[1]))
        ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStrain", 1)
        pts = [(c[0], c[1], 0.0) for c in xy]
    else:
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        _make_mat(1)
        pts = [(float(c[0]), float(c[1]), float(c[2])) for c in _BRICK_XYZ]
        for i, c in enumerate(pts, 1):
            ops.node(i, *c)
            ops.fix(i, 0, 0, 1)
        ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for i, (x, y, _z) in enumerate(pts, 1):
        u = H @ np.array([x, y])
        ops.sp(i, 1, float(u[0]))
        ops.sp(i, 2, float(u[1]))
    ops.constraints("Transformation")     # `Plain` silently DROPS non-homogeneous sp ("homogeneous constraint assumed")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-12, 10, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)        # reset below
    ops.analysis("Static")


def _drive(kind, path_final, nsteps):
    """Returns per step [s11, s22, s12, s33] at Gauss point 1 and the plastic flag."""
    _deck_ps_driver(kind, path_final)
    ops.integrator("LoadControl", 1.0 / nsteps)
    rows = []
    for _ in range(nsteps):
        assert ops.analyze(1) == 0
        if kind == "quad":
            s = ops.eleResponse(1, "stressesPlaneStrain")        # [s11 s22 s12 s33] per Gauss point
            row = [s[0], s[1], s[2], s[3]]
        else:
            s = ops.eleResponse(1, "material", 1, "stress")      # [s11 s22 s33 s12 s23 s13]
            row = [s[0], s[1], s[3], s[2]]
        info = ops.eleResponse(1, "material", 1, "stepInfo")
        rows.append((row, info[1], info[0]))
    eps = ops.eleResponse(1, "material", 1, "strain")           # engineering shear, Voigt order of the lane
    got = [eps[0], eps[1], eps[2] if kind == "quad" else eps[3]]
    assert all(abs(g - w) <= 1e-12 for g, w in zip(got, path_final)), (
        f"premise: the prescribed strain did not reach the element: {got} vs {path_final}")
    return rows


def test_planestrain_sigma33_matches_3d_on_the_same_inplane_path():
    """The PlaneStrain quad and a 3D stdBrick (u_z = 0) are driven through the SAME in-plane strain path
    (compression in y, lateral expansion in x, a shear: eps_22 = -2.5 %, eps_11 = +0.5 %, gamma_12 =
    1.5 %, 50 steps, plastic from step 4 on) with every nodal DOF prescribed (no Newton tolerance can
    enter).  sigma_11, sigma_22, sigma_12 AND sigma_33 must agree at every step to 1e-11 of the largest
    stress (identical kernel, so round-off: measured 5.7e-13 absolute on stresses up to 715 kPa), sigma_33 must be finite (the wrapper REPORTS it: the NaN
    of the base `getStressZZ` would mean not reported) and non-zero, and both lanes must be plastic at
    the end (premise: the path crosses the yield surface).  0.1 s.

    KILLS: a PlaneStrain wrapper that does not override `getStressZZ` (NaN), one that returns
    sigma_11/sigma_22 as sigma_33, one that feeds eps_33 != 0 or a nonzero out-of-plane shear to the
    kernel, a wrong in-plane Voigt map (s12 <-> s33 swapped), a plane-strain lane that skips the
    return map or runs it on a differently-initialised state (v0, pi_i)."""
    path_final = (0.005, -0.025, 0.015)
    nsteps = 50
    quad = _drive("quad", path_final, nsteps)
    brick = _drive("brick", path_final, nsteps)
    scale = max(abs(v) for r, _, _ in brick for v in r)
    assert scale > 100.0, "premise: a non-trivial stress path"
    for k, ((rq, pq, fq), (rb, pb, fb)) in enumerate(zip(quad, brick), 1):
        assert fq == 0 and fb == 0, f"step {k}: refused trial (quad {fq}, brick {fb})"
        assert all(math.isfinite(v) for v in rq), f"step {k}: non-finite component in {rq}"
        assert math.isfinite(rq[3]), f"step {k}: sigma_33 not reported by the PlaneStrain lane"
        for c, name in enumerate(("s11", "s22", "s12", "s33")):
            assert abs(rq[c] - rb[c]) <= 1e-11 * scale, (
                f"step {k}: {name} quad {rq[c]!r} vs 3D brick {rb[c]!r} (tol {1e-11 * scale:.2e})")
    assert quad[-1][1] == 1 and brick[-1][1] == 1, "premise: both lanes plastic at the end of the path"
    s33_final = quad[-1][0][3]
    assert abs(s33_final) > 10.0, f"sigma_33 should carry a confining reaction, got {s33_final}"
    # sigma_33 is not just the mean of the in-plane normals (the wrapper is not a plane-stress shortcut)
    assert abs(s33_final - 0.5 * (quad[-1][0][0] + quad[-1][0][1])) > 1.0


def test_planestrain_sigma33_elastic_closed_form():
    """In the ELASTIC regime the BA06 energy with alpha0 = 0 has the closed form of sheet (S.5):
    p = p0 exp(-eps_v/kappa_hat), s = 2 mu0 e; with eps_33 = 0 and eps_v = eps_11 + eps_22:
        sigma_33 = p - (2/3) mu0 eps_v,   sigma_ii = p + 2 mu0 (eps_ii - eps_v/3),   sigma_12 = mu0 gamma_12.
    The quad (sigma_33 included) must reproduce it to 1e-9 relative at every step of a small path.
    The elastic premise is a THEOREM of the closed form, not read from the shell: the surface through
    this start has F = zeta q + p eta <= (1/rho) q + p eta (zeta <= 1/rho for WW) with eta from (S.12),
    and F < 0 along the whole path is asserted from the closed-form stress alone.  ~0.5 s.

    KILLS: a wrong sigma_33 sign/scale/reduction in `getStressZZ`, a plane-strain lane that applies
    the 3D kernel to eps_33 != 0, a hyperelastic energy with a wrong coupling in the PlaneStrain lane
    (the 3D material is tested elsewhere), an in-plane Voigt map error."""
    path_final = (0.0001, -0.0005, 0.0006)
    nsteps = 5
    rows = _drive("quad", path_final, nsteps)
    for k, (row, plastic, refused) in enumerate(rows, 1):
        a = k / nsteps
        e11, e22, g12 = path_final[0] * a, path_final[1] * a, path_final[2] * a
        ev = e11 + e22
        p = _P0 * math.exp(-ev / _KH)
        exp = [p + 2 * _MU * (e11 - ev / 3.0), p + 2 * _MU * (e22 - ev / 3.0), _MU * g12, p - 2.0 / 3.0 * _MU * ev]
        for c, name in enumerate(("s11", "s22", "s12", "s33")):
            assert abs(row[c] - exp[c]) <= 1e-9 * max(abs(exp[c]), 1.0), (
                f"step {k} {name}: shell {row[c]!r} vs closed form {exp[c]!r}")
        # theorem-level elastic premise (closed-form stress only)
        mean = (exp[0] + exp[1] + exp[3]) / 3.0
        dev2 = sum((s - mean) ** 2 for s in (exp[0], exp[1], exp[3])) + 2.0 * exp[2] ** 2
        q = math.sqrt(1.5 * dev2)
        eta = (_M / _N) * (1.0 - (1.0 - _N) * (mean / _PI0) ** (_N / (1.0 - _N)))
        assert q / _RHO + mean * eta < 0.0, f"step {k}: not provably elastic from the closed form"
        assert plastic == 0, f"step {k}: shell says plastic on a provably elastic step"
    assert abs(rows[-1][0][3] - rows[-1][0][0]) > 1.0     # sigma_33 distinct from sigma_11 (not aliased)


# ===========================================================================
#  (4) Mass density is honoured and does NOT alias the ellipticity -rho
# ===========================================================================
def _newmark_axis_mass(deck_kind, extra, dt=0.1):
    """e_x' A e_x / c3 - (nodal mass), A = the Newmark effective tangent K + c3 M assembled by the real
    integrator (printA), over ALL free DOFs of an unconstrained element with a unit nodal mass on every
    DOF (so the matrix is non-singular even at zero element density).  K annihilates the rigid
    translation e_x, hence the result is the ELEMENT mass summed over the x DOFs = rho_mass * volume,
    whatever the element's mass matrix (consistent here)."""
    ops.wipe()
    if deck_kind == "brick":
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        ops.nDMaterial("LadrunoNorSand", 1, *_mat_args(**extra))
        lx, ly, lz = 1.2, 0.8, 0.5
        for i, c in enumerate(_BRICK_XYZ, 1):
            ops.node(i, c[0] * lx, c[1] * ly, c[2] * lz)
            ops.mass(i, 1.0, 1.0, 1.0)
        ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
        vol, nn, ndf = lx * ly * lz, 8, 3
    else:
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        ops.nDMaterial("LadrunoNorSand", 1, *_mat_args(**extra))
        lx, ly, th = 1.2, 0.8, 0.25
        for i, c in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], 1):
            ops.node(i, c[0] * lx, c[1] * ly)
            ops.mass(i, 1.0, 1.0)
        ops.element("quad", 1, 1, 2, 3, 4, th, "PlaneStrain", 1)
        vol, nn, ndf = lx * ly * th, 4, 2
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-8, 5, 0)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    assert ops.analyze(1, dt) == 0
    c3 = 1.0 / (0.25 * dt * dt)
    N = nn * ndf
    A = np.array(ops.printA("-ret"), dtype=float).reshape(N, N).T
    ex = np.zeros(N)
    ex[0::ndf] = 1.0
    return float(ex @ A @ ex) / c3 - nn * 1.0, vol


@pytest.mark.parametrize("deck_kind", ["brick", "quad"])
def test_density_enters_the_mass_and_rho_does_not(deck_kind):
    """The element mass sums to density x volume (x thickness for the quad): the closed form
    e_x' (K + c3 M) e_x = c3 (nodal masses + density * volume), c3 = 1/(beta dt^2) = 400.  Checked at
    density 2.0 and 3.5 (linear in -density), with the ellipticity -rho/-rho_bar changed at the same
    time (0.7/0.8 -> 0.9/0.95: the mass must not move), and at density 0 with -rho 0.9 (the mass must
    be exactly zero: -rho is NOT a mass density).  Through `getCopy("ThreeDimensional")` (brick) and
    `getCopy("PlaneStrain")` (quad): the density must survive the Gauss-point clones.  ~1 s.

    KILLS: `getRho()` returning the ellipticity `rho` (every case), returning 0 / not honouring
    -density, a `copyFrom` that drops `density` (the clone the element actually uses), a parser that
    stores `-density` into the wrong slot, a mass that is off by the thickness/volume factor."""
    cases = [
        (dict(rho=0.7, rho_bar=0.8, density=2.0), 2.0),
        (dict(rho=0.9, rho_bar=0.95, density=2.0), 2.0),
        (dict(rho=0.9, rho_bar=0.95, density=3.5), 3.5),
        (dict(rho=0.9, rho_bar=0.95), 0.0),
        (dict(rho=0.7, rho_bar=0.8), 0.0),
    ]
    for extra, dens in cases:
        m, vol = _newmark_axis_mass(deck_kind, extra)
        assert abs(m - dens * vol) <= 1e-9 * max(dens * vol, 1.0), (
            f"{deck_kind} {extra}: element mass {m!r} != density x volume {dens * vol!r}")


def test_density_is_in_the_modal_problem_and_rho_is_not():
    """Modal view of the same fact: a unit cube fixed at the bottom, free above (no nodal masses),
    eigenvalues omega^2 of K phi = omega^2 M phi.  Doubling -density halves EVERY eigenvalue (the
    stiffness is untouched), and changing the ellipticity (-rho 0.7 -> 0.9, -rho_bar 0.8 -> 0.95)
    leaves every eigenvalue unchanged to round-off (the elastic start is independent of F).  ~0.3 s.

    KILLS: mass aliased onto -rho (eigenvalues would move with -rho and not with -density), density
    not reaching the Gauss-point clones, a mass matrix built from a stale density."""
    def eigs(extra):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        ops.nDMaterial("LadrunoNorSand", 1, *_mat_args(**extra))
        for i, c in enumerate(_BRICK_XYZ, 1):
            ops.node(i, float(c[0]), float(c[1]), float(c[2]))
        for i in (1, 2, 3, 4):
            ops.fix(i, 1, 1, 1)
        ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("FullGeneral")
        ops.test("NormDispIncr", 1e-8, 5, 0)
        ops.algorithm("Linear")
        ops.integrator("LoadControl", 0.0)
        ops.analysis("Static")
        return np.array(ops.eigen("-fullGenLapack", 12), dtype=float)

    e_a = eigs(dict(rho=0.7, rho_bar=0.8, density=2.0))
    e_b = eigs(dict(rho=0.7, rho_bar=0.8, density=4.0))
    e_c = eigs(dict(rho=0.9, rho_bar=0.95, density=2.0))
    assert np.all(e_a > 0.0)
    np.testing.assert_allclose(e_b, 0.5 * e_a, rtol=1e-9)       # linear in 1/density
    np.testing.assert_allclose(e_c, e_a, rtol=1e-9)             # independent of the ellipticity
