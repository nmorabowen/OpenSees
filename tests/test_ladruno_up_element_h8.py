"""LadrunoUP (ELE 33017) — ADR-71 P2 H8-LANE Zone-B battery (WP2.A OPUS-H8).

Gates the 3D trilinear-hex (H8) lane of the unified Biot u-p element under the
honest-p contract (ADR-71 §3.2): the nodal pressure DOF *is* p, blocks split by
the time-derivative they multiply —

    getTangentStiff  [ K , -Q ; 0 , H ]      getDamp  [ C_Ray , 0 ; Q^T , S+Ht ]
    getMass          [ M , 0  ; 0 , 0 ]      resisting force carries NO rate terms

The three gates (ADR §7 P2 row; plan §WP2.A OPUS-H8):

  1. brickUP two-leg equivalence on a 3D consolidation column (methodology
     ported from the P1 quadUP gate, tests/test_ladruno_up_element_equiv.py):
       leg 1 (tight, consistent mass): Newmark gamma=1/2 beta=1/4,
              -dynSeepage off, -stab off — the unique pair where the p-as-disp
              (ours) / p-as-vel (upstream rate-DOF) parameterizations produce
              the identical one-step rule; u AND p trajectories <= 1e-6 rel
              (achieved ~machine, see test header).  Convergence tests are
              unbalance-based (upstream's int-p-dt DOF grows unboundedly, so
              disp-increment tests mis-scale between the two contracts).
       leg 2 (production): gamma=0.6 beta=0.3025, Dt-halving MUTUAL convergence
              (1e-6 is mathematically impossible at gamma=0.6: the two
              parameterizations differ by a parasitic-root memory term damped
              only for gamma > 1/2, ADR §3.2).
       lumped leg: brickUP DOES expose -lumped, BUT it row-sum-lumps the S
              storage block (BrickUP.cpp:916) whereas LadrunoUP -lumped lumps
              the SOLID mass only (S stays consistent in the damp p-p slot,
              ADR §3.2 getMass "solid block only").  So lumped-vs-lumped is NOT
              operator-identical; it is gated at the honest measured level and
              cross-checked by OUR -lumped vs OUR consistent Dt-convergence
              pair (documented loudly — see the lumped test docstring).

  2. 3D Terzaghi 1D consolidation: H8 1x1x10 column vs the numpy-summed
     Terzaghi series (Tv sweep {0.1, 0.2, 0.5}), settlement U(Tv) check —
     mirrors the P1 2D Terzaghi test's series + tolerance ladder.

  3. B4-3D McGann checkerboard footing (ADR §7.1 3D analog): 30 m cube
     quarter-model, 10x10x10 H8 mesh (h = 3 m), 3x3 m square footing patch on
     the symmetry corner ramped 0->0.1 kPa over 0.1 s then held, drained top
     surface only, Rayleigh 0.05M + 0.02K, run to t = 1 s at gamma=0.6.  3D
     CB_lap index (normalized 6-neighbour interior-Laplacian roughness) for
     alpha0 in {off, 0.25}: CB(off) > CB(0.25) with a healthy ratio, plus the
     auto-vs-manual twin-run identity (1e-10) and the alpha pin
     alpha = 0.25 * h^2 / (Ks_drained + 4Gs/3) +/-5% vs E/nu algebra.

ARG MAPPING — brickUP (upstream), pinned here and verified against the source:

  brickUP (SRC/element/UP-ucsd/BrickUP.cpp, OPS_BrickUP + BrickUP ctor):
      element brickUP  tag N1..N8  matTag  bulk rhof  permX permY permZ
                       <b1 b2 b3>  <-lumped>
    * bulk = the PRE-COMBINED Qbar ~= Kf/n (ctor `kc(bulk)`, BrickUP.cpp:168;
      storage block divides by kc at BrickUP.cpp:916/930).  Pass Kf/poro for
      our (-Kf Kf -poro n) with alpha=1, Ks_grain -> inf (upstream hardwires
      alpha=1).
    * rhof = fluid mass density; drives the seepage body-force term only
      (BrickUP.cpp:1120: resid(jj+3) += dvol*rho*perm_i*b_i*dNp_i).
    * permX/Y/Z = k_hydraulic/gammaW per axis — SAME convention as our -perm
      (their H multiplies grad-p directly, BrickUP.cpp:681-683: damp -= dvol*
      perm_i*dNp*dNp; note the -H sign is upstream's rate-DOF block placement).
    * -lumped -> massType=1: diagonal solid mass (row-sum == HRZ on H8) AND a
      row-sum-lumped S block (BrickUP.cpp:916).  Default massType=0 = consistent.
    * python dispatch name is "brickUP" (OpenSeesElementCommands.cpp:689).
    * upstream physical p = nodal VELOCITY of the p-DOF (rate-DOF trick), read
      via ops.nodeVel(n, 4); LadrunoUP p = nodal DISPLACEMENT of the p-DOF
      (honest-p contract), read via ops.nodeDisp(n, 4).  NO sign flip (both
      compression-positive — verified in P1 for the 2D twin, same convention).
    * H8 local-node order (BrickUP shp3dv == LadrunoUP H8 provider, donor
      SRC/element/brick/shp3d.cpp): bottom face CCW then top face CCW, so the
      same 8 node tags in the same argument slots give byte-identical geometry.

CHECKERBOARD INDEX (gate 3): CB_lap over the FULL interior of the structured
(N+1)^3 node grid at t = 1 s —

  CB_lap = || p_ijk - (1/6) * sum(6 face neighbours) || / || p ||

the 3D generalisation of the P1 CB_lap (interior-Laplacian roughness; captures
ALL short-wave oscillation, monotone under stabilisation).  CB_lat (lattice
cosine with chi_ijk = (-1)^(i+j+k)) is reported as a secondary separation
metric only (it is fragile when heavy smoothing shrinks ||p|| faster than the
lattice projection — the P1 finding).

Run:  py -3.12 -m pytest tests/test_ladruno_up_element_h8.py -x -q
(the worktree dist/bin opensees.pyd is CPython 3.12 — plain `python` is 3.11
and will not import it).
"""
import math
import os
import sys
from pathlib import Path

import numpy as np
import pytest

# --------------------------------------------------------------------------
# bootstrap: THIS worktree's engine (evict any boot-.pth preloaded stale pyd)
# --------------------------------------------------------------------------
_DIST = str(Path(__file__).resolve().parents[1] / "dist" / "bin")
if not os.path.isfile(os.path.join(_DIST, "opensees.pyd")):
    pytest.skip(f"worktree engine not built: {_DIST}", allow_module_level=True)

os.environ["PATH"] = _DIST + os.pathsep + os.environ.get("PATH", "")
try:
    os.add_dll_directory(_DIST)
except (FileNotFoundError, OSError):
    pass
if _DIST not in sys.path:
    sys.path.insert(0, _DIST)
# the installer's ladruno_opensees.pth boot module PRELOADS another build's
# opensees.pyd at interpreter start — evict it so THIS worktree's pyd wins.
for _m in ("opensees", "openseespy", "openseespy.opensees"):
    sys.modules.pop(_m, None)
import opensees as ops  # noqa: E402

assert os.path.normcase(os.path.dirname(ops.__file__)) == os.path.normcase(_DIST), (
    f"wrong opensees.pyd imported: {ops.__file__} (want {_DIST})"
)

pytestmark = [pytest.mark.zone_b]

_GENERAL = ("UmfPack",)   # the honest-p tangent is unsymmetric (ADR §3.2)


# --------------------------------------------------------------------------
# numpy oracles (Terzaghi series — identical to the P1 2D battery)
# --------------------------------------------------------------------------
def terzaghi_series(Z, Tv, nterms=300):
    """Excess-pore-pressure ratio p/p0 of the Terzaghi series.

    Z = z/H measured from the DRAINED face; single-drainage path H."""
    m = np.arange(nterms)
    M = (2 * m + 1) * np.pi / 2.0
    return np.sum((2.0 / M) * np.sin(M * Z) * np.exp(-M * M * Tv))


def terzaghi_U(Tv, nterms=300):
    """Average degree of consolidation U(Tv) of the Terzaghi series."""
    m = np.arange(nterms)
    M = (2 * m + 1) * np.pi / 2.0
    return 1.0 - np.sum((2.0 / M**2) * np.exp(-M * M * Tv))


def _rel_inf(a, b):
    """NaN-proof max-abs relative difference against the max-abs of b."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    assert np.all(np.isfinite(a)) and np.all(np.isfinite(b))
    return float(np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-30))


# ==========================================================================
# shared 3D consolidation column (the brickUP equivalence rig + Terzaghi)
# ==========================================================================
# 1 m x 1 m x H column, nEz stacked H8 u-p elements; laterally confined
# (ux = uy = 0 everywhere -> 1D vertical strain), base fixed + sealed, drained
# top (p = 0), sealed sides; vertical top load on the 4 top nodes.
def _col3_nodes(nEz, Lz, b=1.0):
    """2x2x(nEz+1) node grid; returns {(ix,iy,iz): tag}."""
    nid = {}
    k = 1
    for iz in range(nEz + 1):
        for iy in range(2):
            for ix in range(2):
                ops.node(k, ix * b, iy * b, iz * Lz / nEz)
                nid[(ix, iy, iz)] = k
                k += 1
    return nid


def _col3_hex_nodes(nid, iz):
    """8 node tags of the iz-th hex in standard OpenSees order
    (bottom face CCW, then top face CCW)."""
    return [nid[(0, 0, iz)], nid[(1, 0, iz)], nid[(1, 1, iz)], nid[(0, 1, iz)],
            nid[(0, 0, iz + 1)], nid[(1, 0, iz + 1)],
            nid[(1, 1, iz + 1)], nid[(0, 1, iz + 1)]]


def _col3_bcs(nid, nEz, drained_top=True):
    # one fix per node: ux = uy = 0 everywhere (1D confinement); base uz fixed;
    # p = 0 on the drained top, p free elsewhere (sealed).
    for (ix, iy, iz), n in nid.items():
        uz = 1 if iz == 0 else 0
        fp = 1 if (drained_top and iz == nEz) else 0
        ops.fix(n, 1, 1, uz, fp)


# ==========================================================================
# 1. brickUP two-leg equivalence (3D consolidation column)
# ==========================================================================
COL3 = dict(E=1.0e4, nu=0.3, rho=2.0, poro=0.4, Kf=2.2e5, rhoF=1.0,
            perm=1.0e-4, Lz=4.0, nEz=4, load=-10.0)


def _build_col3(which, gamma, beta, lumped=False, dyn_seepage="off"):
    """which in {'ladruno','brickup'}.  Returns the node-id dict."""
    c = COL3
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 4)
    nid = _col3_nodes(c["nEz"], c["Lz"])
    ops.nDMaterial("ElasticIsotropic", 1, c["E"], c["nu"], c["rho"])
    for iz in range(c["nEz"]):
        hexn = _col3_hex_nodes(nid, iz)
        e = iz + 1
        if which == "ladruno":
            args = ["LadrunoUP", e, *hexn, 1,
                    "-Kf", c["Kf"], "-poro", c["poro"], "-rhoF", c["rhoF"],
                    "-perm", c["perm"], c["perm"], c["perm"],
                    "-dynSeepage", dyn_seepage, "-stab", "off"]
            if lumped:
                args.append("-lumped")
            ops.element(*args)
        elif which == "brickup":
            # ARG MAPPING: bulk = pre-combined Qbar ~= Kf/n (module header)
            args = ["brickUP", e, *hexn, 1,
                    c["Kf"] / c["poro"], c["rhoF"],
                    c["perm"], c["perm"], c["perm"], 0.0, 0.0, 0.0]
            if lumped:
                args.append("-lumped")
            ops.element(*args)
        else:
            raise ValueError(which)
    _col3_bcs(nid, c["nEz"], drained_top=True)
    # identical ramp for every dt: 0 -> 1 over [0, 0.05], then held
    ops.timeSeries("Path", 1, "-time", 0.0, 0.05, 1.0e6,
                   "-values", 0.0, 1.0, 1.0)
    ops.pattern("Plain", 1, 1)
    for iy in range(2):
        for ix in range(2):
            ops.load(nid[(ix, iy, c["nEz"])], 0.0, 0.0, c["load"] / 4.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system(*_GENERAL)                       # general solver: honest-p tangent
    ops.test("NormUnbalance", 1e-10, 30)
    ops.algorithm("Newton")
    ops.integrator("Newmark", gamma, beta)
    ops.analysis("Transient")
    return nid


def _run_col3(which, gamma, beta, nsteps, dt, **kw):
    """Return (u_hist, p_hist): per-step arrays over every node.

    u = (ux,uy,uz) of all nodes; p read per the CONTRACT CONVERSION —
    nodeDisp dof-4 for LadrunoUP, nodeVel dof-4 for the upstream rate-DOF
    element (module header)."""
    nid = _build_col3(which, gamma, beta, **kw)
    keys = sorted(nid)
    u_hist, p_hist = [], []
    for _ in range(nsteps):
        assert ops.analyze(1, dt) == 0, f"{which} transient step failed"
        u_row, p_row = [], []
        for kk in keys:
            n = nid[kk]
            d = ops.nodeDisp(n)
            u_row.extend(d[:3])
            p_row.append(d[3] if which == "ladruno" else ops.nodeVel(n, 4))
        u_hist.append(u_row)
        p_hist.append(p_row)
    return np.array(u_hist), np.array(p_hist)


@pytest.mark.t1
def test_brickup_equivalence_leg1_tight():
    """LadrunoUP (-dynSeepage off, -stab off, consistent mass) vs brickUP under
    Newmark gamma=1/2 beta=1/4: full u AND p trajectories (all nodes, 40 steps)
    <= 1e-6 relative.  The two parameterizations (p-as-disp vs the upstream
    p-as-vel rate-DOF trick) produce the machine-identical one-step rule at this
    unique (gamma, beta), pinning the arg mapping AND the contract conversion
    (their p = nodeVel dof-4, ours = nodeDisp dof-4) in one shot.

    First-passing-run: u 1.1e-15, p 2.8e-15 (machine — the arg map is exact)."""
    nsteps, dt = 40, 0.05
    uL, pL = _run_col3("ladruno", 0.5, 0.25, nsteps, dt)
    uU, pU = _run_col3("brickup", 0.5, 0.25, nsteps, dt)
    rel_u = _rel_inf(uL, uU)
    rel_p = _rel_inf(pL, pU)
    # non-trivial response guard (a dead model would pass any rel gate)
    assert np.max(np.abs(uU)) > 1e-6 and np.max(np.abs(pU)) > 1e-2
    assert rel_u <= 1e-6, f"leg-1 u trajectories diverge: rel {rel_u:.3e}"
    assert rel_p <= 1e-6, f"leg-1 p trajectories diverge: rel {rel_p:.3e}"


@pytest.mark.t1
def test_brickup_equivalence_leg2_production_convergence():
    """gamma=0.6 beta=0.3025 (the ZS84 production set): the p-row memory term
    makes the two parameterizations genuinely different one-step rules, so the
    gate is MUTUAL Dt-convergence, not coincidence: the end-state difference
    |u_L - u_U| must shrink under every Dt halving (~first order) and the
    finest-Dt relative difference must be <= 1%."""
    T = 2.0
    diffs_u, rels_u, rels_p = [], [], []
    for dt in (0.1, 0.05, 0.025, 0.0125):
        ns = int(round(T / dt))
        uL, pL = _run_col3("ladruno", 0.6, 0.3025, ns, dt)
        uU, pU = _run_col3("brickup", 0.6, 0.3025, ns, dt)
        diffs_u.append(float(np.max(np.abs(uL[-1] - uU[-1]))))
        rels_u.append(_rel_inf(uL[-1], uU[-1]))
        rels_p.append(_rel_inf(pL[-1], pU[-1]))
    for k in range(1, len(diffs_u)):
        assert diffs_u[k] < diffs_u[k - 1], (
            f"no mutual Dt-convergence: diffs {diffs_u}")
    order = math.log(diffs_u[0] / diffs_u[-1], 2) / (len(diffs_u) - 1)
    assert order >= 0.5, f"mutual convergence order {order:.2f} < 0.5"
    assert rels_u[-1] <= 1e-2, f"finest-Dt u rel {rels_u[-1]:.3e} > 1%"
    assert rels_p[-1] <= 1e-2, f"finest-Dt p rel {rels_p[-1]:.3e} > 1%"


@pytest.mark.t1
def test_brickup_lumped_leg():
    """The -lumped leg, both ways.

    brickUP DOES expose -lumped, but it row-sum-lumps the S STORAGE block
    (BrickUP.cpp:916) while LadrunoUP -lumped lumps only the SOLID mass (S
    stays consistent in the getDamp p-p slot, ADR §3.2 getMass "solid block
    only").  So lumped-vs-lumped is NOT the operator-identity that leg 1 is —
    the two schemes differ by the S-lumping quadrature.  This test therefore:

      (a) MEASURES LadrunoUP -lumped vs brickUP -lumped at gamma=1/2 beta=1/4
          and gates it at the honest measured level (the solid-mass lumping IS
          identical — row-sum == HRZ on H8 — so the only divergence is the S
          block, which on this consolidation column is a small perturbation);
      (b) cross-checks OUR -lumped vs OUR consistent as a Dt-convergence pair:
          lumped mass is a consistent O(h^2)-class perturbation, so the two
          settle onto the same solution as dt shrinks (finest-dt rel 1.2e-6
          at dt=0.025 — quasi-static column, mass barely participates).
    """
    # (a) lumped-vs-lumped at the tight (gamma=1/2) pair — MEASURED tolerance.
    nsteps, dt = 40, 0.05
    uLl, pLl = _run_col3("ladruno", 0.5, 0.25, nsteps, dt, lumped=True)
    uUl, pUl = _run_col3("brickup", 0.5, 0.25, nsteps, dt, lumped=True)
    assert np.max(np.abs(uUl)) > 1e-6 and np.max(np.abs(pUl)) > 1e-2
    rel_u = _rel_inf(uLl, uUl)
    rel_p = _rel_inf(pLl, pUl)
    # honest bound: solid mass identical, S-lumping differs (see docstring).
    # First-passing-run pins these (S-lumping perturbation, NOT a bug):
    # u 5.2e-3, p 4.0e-2 — i.e. leg-1's 1e-15 identity is lost EXACTLY at the
    # S-lumping quadrature, which is the honest, documented divergence.
    assert rel_u <= 1e-2, f"lumped-vs-lumped u rel {rel_u:.3e} (S-lump gap)"
    assert rel_p <= 6e-2, f"lumped-vs-lumped p rel {rel_p:.3e} (S-lump gap)"

    # (b) OUR -lumped vs OUR consistent: Dt-convergence pair.
    T = 2.0
    diffs, rels_u2 = [], []
    for dtc in (0.1, 0.05, 0.025):
        ns = int(round(T / dtc))
        uc, pc = _run_col3("ladruno", 0.6, 0.3025, ns, dtc, lumped=False)
        ul, pl = _run_col3("ladruno", 0.6, 0.3025, ns, dtc, lumped=True)
        diffs.append(float(np.max(np.abs(uc[-1] - ul[-1]))))
        rels_u2.append(_rel_inf(ul[-1], uc[-1]))
    for k in range(1, len(diffs)):
        assert diffs[k] < diffs[k - 1], (
            f"our lumped/consistent do not mutually converge: {diffs}")
    assert rels_u2[-1] <= 2e-2, f"finest-dt lumped/consistent u rel {rels_u2[-1]:.3e}"


# ==========================================================================
# 2. 3D Terzaghi 1D consolidation vs series
# ==========================================================================
@pytest.mark.t1
def test_terzaghi_3d_consolidation_vs_series():
    """H8 1x1x10 column, drained top / sealed sides+bottom, instantaneous load.

    p(z,t) vs the numpy Terzaghi series at Tv in {0.1, 0.2, 0.5} and the
    settlement degree of consolidation U(Tv) — mirrors the P1 2D Terzaghi test
    (tests/test_ladruno_up_element_analytic.py) in 3D.

    TOLERANCE RATIONALE (discretization-aware, mesh h = H/10): the spatial
    error of the trilinear pair is O(h^2) once the pressure front spans a few
    elements (gated 2% / 1% / 0.5% at Tv = 0.1 / 0.2 / 0.5, matching the 2D
    ladder).  gamma=0.6 beta=0.3025, dt = 0.02 s, -dynSeepage off (the
    default 'on' pollutes the p-field via Newmark-recovered accel noise,
    P1 finding)."""
    E, nu, rho = 1.0e4, 0.2, 2.0          # kPa, Mg/m3
    Kf, poro, rhoF, kbar = 4444.0, 0.4, 1.0, 1.0e-4
    L, q, nEz = 10.0, 10.0, 10
    Eoed = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu))
    ooQ = poro / Kf
    Qbar = 1.0 / ooQ
    cv = kbar / (1.0 / Eoed + ooQ)
    p0 = q * Qbar / (Qbar + Eoed)          # undrained (compressible-fluid) p
    s0 = q * L / (Eoed + Qbar)             # undrained settlement
    sinf = q * L / Eoed                    # drained settlement

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 4)
    nid = _col3_nodes(nEz, L)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)
    for iz in range(nEz):
        ops.element("LadrunoUP", iz + 1, *_col3_hex_nodes(nid, iz), 1,
                    "-Kf", Kf, "-poro", poro, "-rhoF", rhoF,
                    "-perm", kbar, kbar, kbar, "-stab", "off",
                    "-dynSeepage", "off")
    _col3_bcs(nid, nEz, drained_top=True)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for iy in range(2):
        for ix in range(2):
            ops.load(nid[(ix, iy, nEz)], 0.0, 0.0, -q / 4.0, 0.0)
    ops.system(*_GENERAL)
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")

    assert ops.analyze(1, 1.0e-3) == 0
    t = 1.0e-3
    dt = 0.02
    H = L
    settle, perr = {}, {}
    for Tv in (0.1, 0.2, 0.5):
        tt = Tv * H * H / cv
        nst = int(round((tt - t) / dt))
        assert ops.analyze(nst, dt) == 0
        t += nst * dt
        TvA = t * cv / (H * H)
        errs = [abs(ops.nodeDisp(nid[(0, 0, iz)], 4)
                    - p0 * terzaghi_series((L - L * iz / nEz) / H, TvA)) / p0
                for iz in range(1, nEz)]
        perr[Tv] = max(errs)
        settle[Tv] = -ops.nodeDisp(nid[(0, 0, nEz)], 3)

    assert perr[0.1] < 0.02, f"Tv=0.1 max rel err {perr[0.1]:.3%}"
    assert perr[0.2] < 0.01, f"Tv=0.2 max rel err {perr[0.2]:.3%}"
    assert perr[0.5] < 0.005, f"Tv=0.5 max rel err {perr[0.5]:.3%}"

    U_fem = {Tv: (settle[Tv] - s0) / (sinf - s0) for Tv in settle}
    for Tv in (0.1, 0.2, 0.5):
        assert abs(U_fem[Tv] - terzaghi_U(Tv)) < 0.02, (
            f"U(Tv={Tv}) fem {U_fem[Tv]:.4f} vs series {terzaghi_U(Tv):.4f}")
    ratio = U_fem[0.2] / U_fem[0.1]
    assert abs(ratio - np.sqrt(2.0)) < 0.07, (
        f"early-time settlement is not sqrt(t): U(0.2)/U(0.1) = {ratio:.3f}")


# ==========================================================================
# 3. B4-3D McGann checkerboard footing — CB index + auto/manual + alpha pin
# ==========================================================================
B4 = dict(E=25000.0, nu=0.3, rho=2.67, poro=0.4, Kf=2.2e12, rhoF=1.0,
          perm=1.0e-7, N=10, L=30.0, q=0.1, strip=3.0, dt=0.025, T=1.0)


def _b4_stab_alpha_3d(h, E, nu, alpha0):
    """alpha = alpha0 * h^2 / (Ks_drained + 4Gs/3), skeleton moduli from E,nu
    (McGann; ADR §3.3).  h = largest element EDGE."""
    Ks = E / (3.0 * (1.0 - 2.0 * nu))
    Gs = E / (2.0 * (1.0 + nu))
    return alpha0 * h * h / (Ks + 4.0 * Gs / 3.0)


B4_ALPHA_TARGET = _b4_stab_alpha_3d(B4["L"] / B4["N"], B4["E"], B4["nu"], 0.25)


def _b4_3d_build(stab):
    """30 m cube quarter-model, N^3 H8 mesh (h = 3 m), 3x3 m square footing
    patch on the symmetry corner ramped 0 -> 0.1 kPa over 0.1 s then held;
    drainage TOP surface only; symmetry rollers on x=0/y=0, far-field lateral
    restraint on x=L/y=L, fixed base uz; Rayleigh 0.05M + 0.02K; Newmark
    0.6/0.3025.  stab: 'off' | ('auto', a0) | ('manual', alpha)."""
    b = B4
    N = b["N"]
    h = b["L"] / N
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 4)
    nid = {}
    k = 1
    for iz in range(N + 1):
        for iy in range(N + 1):
            for ix in range(N + 1):
                ops.node(k, ix * h, iy * h, iz * h)
                nid[(ix, iy, iz)] = k
                k += 1
    ops.nDMaterial("ElasticIsotropic", 1, b["E"], b["nu"], b["rho"])
    e = 1
    for iz in range(N):
        for iy in range(N):
            for ix in range(N):
                hexn = [nid[(ix, iy, iz)], nid[(ix + 1, iy, iz)],
                        nid[(ix + 1, iy + 1, iz)], nid[(ix, iy + 1, iz)],
                        nid[(ix, iy, iz + 1)], nid[(ix + 1, iy, iz + 1)],
                        nid[(ix + 1, iy + 1, iz + 1)], nid[(ix, iy + 1, iz + 1)]]
                args = ["LadrunoUP", e, *hexn, 1, "-Kf", b["Kf"],
                        "-poro", b["poro"], "-rhoF", b["rhoF"],
                        "-perm", b["perm"], b["perm"], b["perm"],
                        "-dynSeepage", "off"]
                if stab == "off":
                    args += ["-stab", "off"]
                elif stab[0] == "auto":
                    args += ["-stab", "auto", stab[1]]
                else:
                    args += ["-stab", stab[1]]
                ops.element(*args)
                e += 1
    for (ix, iy, iz), n in nid.items():
        fx = 1 if (ix == 0 or ix == N) else 0     # x=0 symmetry, x=L far-field
        fy = 1 if (iy == 0 or iy == N) else 0     # y=0 symmetry, y=L far-field
        fz = 1 if iz == 0 else 0                  # fixed base
        fp = 1 if iz == N else 0                  # drained TOP only
        if fx or fy or fz or fp:
            ops.fix(n, fx, fy, fz, fp)
    ops.timeSeries("Path", 1, "-time", 0.0, 0.1, 1.0e6,
                   "-values", 0.0, 1.0, 1.0)
    ops.pattern("Plain", 1, 1)
    # 3x3 m square patch = the corner (ix,iy in {0,1}) top element; consistent
    # nodal loads q*area/4 on its 4 top nodes.
    fnode = -b["q"] * b["strip"] * b["strip"] / 4.0
    for iy in range(2):
        for ix in range(2):
            ops.load(nid[(ix, iy, N)], 0.0, 0.0, fnode, 0.0)
    ops.rayleigh(0.05, 0.02, 0.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system(*_GENERAL)
    ops.test("NormUnbalance", 1e-8, 50)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")
    return nid


def _b4_3d_cb_indices(nid):
    """(CB_lat, CB_lap) over the FULL interior of the (N+1)^3 node grid at the
    current state (definitions in the module header)."""
    N = B4["N"]
    P = {(ix, iy, iz): ops.nodeDisp(nid[(ix, iy, iz)], 4)
         for ix in range(N + 1) for iy in range(N + 1) for iz in range(N + 1)}
    inter = [(ix, iy, iz) for ix in range(1, N)
             for iy in range(1, N) for iz in range(1, N)]
    p = np.array([P[k] for k in inter])
    assert np.all(np.isfinite(p))
    chi = np.array([(-1.0) ** (ix + iy + iz) for (ix, iy, iz) in inter])
    cb_lat = abs(float(p @ chi)) / max(
        float(np.linalg.norm(p) * np.linalg.norm(chi)), 1e-30)
    lap = np.array([P[(ix, iy, iz)] - (1.0 / 6.0) * (
        P[(ix + 1, iy, iz)] + P[(ix - 1, iy, iz)] +
        P[(ix, iy + 1, iz)] + P[(ix, iy - 1, iz)] +
        P[(ix, iy, iz + 1)] + P[(ix, iy, iz - 1)])
        for (ix, iy, iz) in inter])
    cb_lap = float(np.linalg.norm(lap)) / max(float(np.linalg.norm(p)), 1e-30)
    return cb_lat, cb_lap


def _b4_3d_run(stab):
    nid = _b4_3d_build(stab)
    ns = int(round(B4["T"] / B4["dt"]))
    for s in range(ns):
        assert ops.analyze(1, B4["dt"]) == 0, f"B4-3D {stab} failed at step {s}"
    return _b4_3d_cb_indices(nid)


@pytest.fixture(scope="module")
def b4_3d():
    """One pass over {off, auto 0.25, manual-alpha twin of 0.25} shared by the
    B4-3D tests below (the alpha0 in {off, 0.25} pair per WP2.A)."""
    out = {"off": _b4_3d_run("off"),
           "auto": _b4_3d_run(("auto", 0.25))}
    out["manual"] = _b4_3d_run(("manual", B4_ALPHA_TARGET))
    return out


@pytest.mark.t1
def test_b4_3d_off_vs_auto(b4_3d):
    """-stab off checkerboards (CB_lap HIGH); -stab auto 0.25 suppresses it
    (CB_lap LOW), ratio > 2 (ADR expectation).  Gates pinned from the first
    passing run per the ADR §7.1 B4 protocol (3D analog): CB_lap off = 0.409,
    auto 0.25 = 0.128, ratio = 3.20 (CB_lat off 0.0391, auto 0.0085)."""
    lat_off, lap_off = b4_3d["off"]
    lat_auto, lap_auto = b4_3d["auto"]
    assert lap_off > lap_auto, (
        f"stabilization did not reduce CB_lap: off {lap_off:.3f} "
        f"auto {lap_auto:.3f}")
    assert lap_off / max(lap_auto, 1e-30) >= 2.0, (
        f"CB_lap suppression ratio too small: off {lap_off:.3f} / "
        f"auto {lap_auto:.3f} = {lap_off / max(lap_auto, 1e-30):.2f}")


@pytest.mark.t1
def test_b4_3d_auto_manual_identity(b4_3d):
    """auto-alpha twin-run identity (the P1 trick): '-stab auto 0.25' and
    '-stab <manual alpha>' with the SAME first-principles alpha produce
    identical CB indices to 1e-10 — the element's internal auto = the pinned
    formula, observable only through the twin run (no direct alpha response)."""
    lat_a, lap_a = b4_3d["auto"]
    lat_m, lap_m = b4_3d["manual"]
    assert abs(lat_a - lat_m) <= 1e-10 * max(abs(lat_m), 1e-30) + 1e-14, (
        f"auto/manual CB_lat differ: {lat_a} vs {lat_m}")
    assert abs(lap_a - lap_m) <= 1e-10 * max(abs(lap_m), 1e-30) + 1e-14, (
        f"auto/manual CB_lap differ: {lap_a} vs {lap_m}")


@pytest.mark.t1
def test_b4_3d_stab_alpha_value():
    """The pinned manual stabilization alpha equals 0.25 * h^2 / (Ks + 4Gs/3)
    with h = largest EDGE = 3 m exactly on this mesh, and sits within 5% of the
    ADR B4 target back-calculated from E=25000 kPa, nu=0.3 (Ks = E/(3(1-2nu)) =
    20833.3, Gs = E/(2(1+nu)) = 9615.4, Ks+4Gs/3 = 33653.8; alpha =
    0.25*9/33653.8 = 6.686e-5)."""
    h = B4["L"] / B4["N"]
    assert abs(h - 3.0) < 1e-12, f"mesh edge h = {h} (want 3 m exactly)"
    Ks = B4["E"] / (3.0 * (1.0 - 2.0 * B4["nu"]))
    Gs = B4["E"] / (2.0 * (1.0 + B4["nu"]))
    alpha_ref = 0.25 * 9.0 / (Ks + 4.0 * Gs / 3.0)
    assert abs(B4_ALPHA_TARGET - alpha_ref) <= 0.05 * alpha_ref, (
        f"pinned alpha {B4_ALPHA_TARGET:.4e} vs first-principles {alpha_ref:.4e}")
    # sanity: the value itself (6.686e-5) is the ADR B4 class
    assert abs(B4_ALPHA_TARGET - 6.686e-5) <= 0.02 * 6.686e-5
