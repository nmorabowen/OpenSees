"""Why do QUADRATIC elements lose the Newton path before their mechanism forms?

NOT a collapse-load measurement.  Note 81 already established -- and then had to
retract half of -- the collapse-load reading of this problem.  What survives it
is one empirical observable and one open question:

    every quadratic leg ever measured on the Prandtl-Reissner oracle (hex and
    tet, Lagrange and Bernstein, volumetrically relieved and locked) WALLS while
    still hardening at 10-27 % of its initial tangent, while the one LINEAR
    relieved element (`LadrunoBrick -formulation bbar`) plateaus and reproduces
    the exact answer -- in the same harness.

    WHY?

This script is the instrumentation note 81 section 5.5 asked for, run on the
cell that already exists rather than on the H27 that does not (see the H27
BLOCKAGE note below).  It answers four questions and refuses to answer a fifth.

WHAT IT MEASURES
----------------
1.  MECHANISM FORMATION AT TERMINATION (the one that carries the weight).  At
    the last converged step of each leg it dumps the DruckerPrager MOBILISATION
    field -- ||dev sigma|| / (sqrt(2/3) sigma_y - rho I1), 1.0 at yield, the
    same quantity the reference driver `r3_prandtl_tet10.py` writes -- per Gauss
    point, plus the INCREMENTAL displacement field (the mechanism IS a velocity
    field, so the increment is the more direct observable) and the surface
    heave profile.  Note 81 section 5.3's reframing is that the question is the
    SPAN of the admissible isochoric fields, not a COUNT of them, so:

      * plastic VOLUME fraction, not element count -- this mesh is graded 1.35x
        and a count is dominated by far-field elements 10x the size of the core
        ones.  Both are reported; the volume one is the headline.
      * CONNECTIVITY of the yielded set on the structured (i,k) element grid,
        and the size of its largest component.
      * BAND-NESS of that component: its area in the x-z plane against its
        bounding box, and an equivalent band thickness in elements.  A slip
        surface is 1-3 elements thick and crosses the domain; smeared
        plasticity fills its own bounding box.
      * HEAVE RATIO -- the decisive one, and the cheapest.  At psi = 0 the
        plastic flow is ISOCHORIC, so a footing that penetrates by ds through a
        collapse MECHANISM must push an equal volume of soil UP somewhere else.
        If the incremental surface displacement beside the footing is ~ +ds the
        mechanism spans; if it is ~ 0 the element is accommodating penetration
        by elastic COMPACTION instead, which is exactly the span failure.

2.  DOES THE WALL MOVE WHEN THE PATH CONTROLLER'S ALLOWANCE CHANGES?  `--floor`
    and `--budget` are the two allowances.  Note 81 section 4.2 reported that a
    10x change in the step ladder left the H20 `uri` wall where it was and read
    that as (weak) evidence of an element property.  THAT PAIR CANNOT HAVE
    MOVED IT: it scaled DS_BASE and DS_MIN TOGETHER by 10, so their RATIO -- and
    therefore the number of halvings the controller is allowed before it hits
    the floor -- was IDENTICAL in both runs (2e-4/2e-6 and 2e-5/2e-7 are both
    100:1 = 6.64 halvings).  It was a rescaling, not an allowance change.  The
    allowance experiment is to move the floor ALONE, which is what `--floor`
    does, and it is the direct analogue of the TIMs subdivision-budget sweep.

3.  TANGENT CONDITIONING AS THE WALL APPROACHES (`--cond`).  Once the step has
    been subdivided below COND_TRIGGER x DS_BASE the run switches to
    `FullGeneral` and samples the assembled tangent every `--cond-every`
    converged steps: sigma_min / sigma_max / cond by SVD (the tangent is
    UNSYMMETRIC at rho_bar != rho, so eigenvalues alone are the wrong tool), and
    the least eigenvalue of the SYMMETRIC PART, whose sign change is the
    loss-of-positive-definiteness that a real limit point produces.

4.  THE TERMINATION-MODE TABLE -- the deliverable.  The classification is the
    merged gate's (`tests/test_r3_prandtl_collapse_gate.py`), three clauses, not
    one: a leg is a CAPACITY only if it plateaued AND terminated in an
    admissible mode (TARGET or the NAMED allowance BUDGET) AND was still
    advancing freely (tail step >= 100x the floor) when it flattened.  FLOOR,
    WALL, STALL and TRUNCATED are SEIZURE and are never a capacity.

WHAT IT REFUSES TO ANSWER
-------------------------
A collapse load for any quadratic element.  Note 81 section 5.2 stands: every
quadratic leg ends on a path-controller limit with the curve still climbing, and
an independent leg has been shown to move +17 % in capacity under a change to
one solver integer.  Numbers printed here are labelled by their termination mode
and a number without its mode must not be quoted.

H27 BLOCKAGE (carried from the follow-up brief, and it is load-bearing)
----------------------------------------------------------------------
The follow-up was scoped as "H27 + selective".  THE FORK OWNS NO 27-NODE HEX.
`LadrunoBrick20` is 20-node SERENDIPITY; every "27" in `SRC/element/brick/` is a
3x3x3 Gauss-point count, not a node count.  Building one is element development,
not a measurement, so it is out of scope here and nothing below is a substitute
for it.  The recommendation on whether to open that lane is in the note.

INTERNAL / HOLD: TIMs campaign figures referenced above are the campaign's,
cited with credit, and are under the same hold as the unfiled UW report -- keep
them out of anything public or upstream-facing.

Run:
    py -3.12 quad_path_diag.py --elem h20uri --h0 1.0
    py -3.12 quad_path_diag.py --elem h20uri --h0 1.0 --floor 2e-9 --budget 800 \
                               --suffix _loosefloor
    py -3.12 quad_path_diag.py --elem h8bbar --h0 1.0 --ladder linear --budget 80
"""
import argparse
import csv
import json
import math
import os
import sys
import time
from collections import deque

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

# Every piece of problem definition, mesh generation and load algebra below is
# h20_prandtl.py's, imported rather than copied: the consistent Q8 surcharge and
# its three controls were bought expensively and must not fork.
import h20_prandtl as HP                                        # noqa: E402

# H20_NO_ENGINE=1 imports the geometry / localisation half WITHOUT the solver so
# the metrics can be unit-tested under any interpreter (h20_prandtl.py's idiom).
ops = HP.ops

SQ23 = math.sqrt(2.0 / 3.0)

# --- the capacity classification, verbatim from the merged R3 gate ----------
PLATEAU_FRAC = 0.02            # tail dq/ds < 2 % of the initial tangent
FREE_ADVANCE_FLOOR_FACTOR = 100.0
DROP_TRUNCATE = 0.02
_CAPACITY_MODES = ("TARGET", "BUDGET")
_SEIZURE_MODES = ("FLOOR", "WALL", "STALL", "TRUNCATED")

# --- diagnostic 3 trigger ---------------------------------------------------
COND_TRIGGER = 0.25            # switch to FullGeneral once ds < this x DS_BASE

# --- diagnostic 1 thresholds ------------------------------------------------
MOB_YIELD = 0.90               # "at yield" for the localisation census
MOB_HARD = 0.99
# the CORE: the region where the graded mesh is still at h0, so a localisation
# metric taken there is not reading the grading.  +-4 m = 2 B either side of the
# footing centreline, 4 m deep = 2 B.
CORE_X, CORE_Z = 4.0, -4.0

ELEMS = {"h8bbar": dict(order=1, form="bbar", label="LadrunoBrick -formulation bbar"),
         "h8std":  dict(order=1, form="std",  label="LadrunoBrick -formulation std"),
         "h20uri": dict(order=2, form="uri",  label="LadrunoBrick20 -formulation uri"),
         "h20std": dict(order=2, form="std",  label="LadrunoBrick20 -formulation std")}

log = HP.log


# ---------------------------------------------------------------------------
# the mobilisation field
# ---------------------------------------------------------------------------
def mobilisation(ncells, rho):
    """Per-element mobilisation of the Drucker-Prager cone, from EVERY Gauss
    point (the reference driver samples material point 1 only; this is a strict
    superset and reports the mean and the max separately).

        F = ||dev sigma|| + rho I1 - sqrt(2/3) sigma_y      (OpenSees sign
        convention, compression negative) so the available deviatoric radius at
        the current pressure is  cap = sqrt(2/3) sigma_y - rho I1  and
        mob = ||dev sigma|| / cap  is 1.0 exactly at yield.

    Verified against the harness's own initial-state formula: at the 1-D state
    sigma_zz = -q0, sigma_xx = sigma_yy = -K0 q0 this reduces algebraically to
    m0 = (1-K0)/(sqrt(3) alpha (1+2 K0)), which is the number h20_prandtl.py
    prints as `initial 1-D state at m0`.  The two agreeing at the surcharge step
    is control M0 below.
    """
    mean = np.full(ncells, np.nan)
    mx = np.full(ncells, np.nan)
    ngp = 0
    for e in range(1, ncells + 1):
        st = np.asarray(ops.eleResponse(e, "stress"), dtype=float)
        if st.size < 6:
            continue
        st = st.reshape(-1, 6)
        ngp = st.shape[0]
        I1 = st[:, 0] + st[:, 1] + st[:, 2]
        p = I1 / 3.0
        d0, d1, d2 = st[:, 0] - p, st[:, 1] - p, st[:, 2] - p
        ns = np.sqrt(d0 ** 2 + d1 ** 2 + d2 ** 2
                     + 2.0 * (st[:, 3] ** 2 + st[:, 4] ** 2 + st[:, 5] ** 2))
        cap = SQ23 * HP.SY - rho * I1
        m = np.where(cap > 1e-12, ns / np.maximum(cap, 1e-30), np.nan)
        mean[e - 1] = float(np.nanmean(m))
        mx[e - 1] = float(np.nanmax(m))
    return mean, mx, ngp


def element_geometry(nodes, cells, x, z):
    """Centroid, volume and structured (i,k) index of every element.

    `strip_mesh` emits cells in the order
        for i in range(nx-1): for j in range(ny-1): for k in range(nz-1)
    and this strip is ONE element thick in y (ny-1 == 1), so e == i*(nz-1)+k
    exactly.  Asserted against the geometry rather than assumed.
    """
    p = nodes[cells[:, :8]]                       # the corner sub-hex
    cen = p.mean(axis=1)
    lo, hi = p.min(axis=1), p.max(axis=1)
    vol = np.prod(hi - lo, axis=1)                # exact: the grid is axis-aligned
    nz1 = len(z) - 1
    ii = np.arange(len(cells)) // nz1
    kk = np.arange(len(cells)) % nz1
    xc = 0.5 * (x[:-1] + x[1:])
    zc = 0.5 * (z[:-1] + z[1:])
    assert np.allclose(cen[:, 0], xc[ii]) and np.allclose(cen[:, 2], zc[kk]), \
        "the (i,k) de-indexing of the cell array is wrong"
    return cen, vol, ii, kk


_NBR8 = ((1, 0), (-1, 0), (0, 1), (0, -1),
         (1, 1), (1, -1), (-1, 1), (-1, -1))


def _components(ii, kk, sel, nx1, nz1):
    """Connected components of the selected element set on the (i,k) grid.

    EIGHT-connectivity (face OR edge neighbours), and the choice matters: a
    Prandtl slip surface runs DIAGONALLY across a structured grid, so a
    one-element-thick diagonal band is not 4-connected at all and a 4-connected
    census shatters exactly the pattern this diagnostic exists to find.
    Measured on a synthetic 1-element-thick band on this mesh: 5 components
    under 4-connectivity, 1 under 8.  With 8-connectivity, "the yielded set is
    NOT connected" is therefore a strong negative result rather than an
    artefact of the stencil.

    Iterative flood fill -- a recursive one blows the stack on a
    domain-spanning plastic zone.
    """
    grid = -np.ones((nx1, nz1), dtype=np.int64)
    idx = np.where(sel)[0]
    grid[ii[idx], kk[idx]] = idx
    lab = -np.ones(len(ii), dtype=np.int64)
    comps = []
    for e in idx:
        if lab[e] >= 0:
            continue
        cid = len(comps)
        stack, members = [e], []
        lab[e] = cid
        while stack:
            c = stack.pop()
            members.append(c)
            i0, k0 = ii[c], kk[c]
            for di, dk in _NBR8:
                i1, k1 = i0 + di, k0 + dk
                if 0 <= i1 < nx1 and 0 <= k1 < nz1:
                    n = grid[i1, k1]
                    if n >= 0 and lab[n] < 0:
                        lab[n] = cid
                        stack.append(n)
        comps.append(np.array(members, dtype=np.int64))
    return comps, lab


def localisation_metrics(nodes, cells, cen, vol, ii, kk, mob, x, z, tag):
    """Is the yielded set a LOCALISED SLIP SURFACE or SMEARED plasticity?

    Every fraction here is VOLUME weighted.  A count-weighted fraction on this
    1.35x-graded mesh is dominated by the far field, where the elements are an
    order of magnitude larger and always elastic, so a count would read
    "localised" for a completely diffuse field.  Counts are reported too, as the
    control that shows the two do not agree.
    """
    nx1, nz1 = len(x) - 1, len(z) - 1
    ok = np.isfinite(mob)
    core = (np.abs(cen[:, 0]) <= CORE_X) & (cen[:, 2] >= CORE_Z)
    out = {}
    for name, thr in (("y90", MOB_YIELD), ("y99", MOB_HARD)):
        sel = ok & (mob >= thr)
        out[f"vfrac_{name}"] = float(vol[sel].sum() / vol[ok].sum())
        out[f"cfrac_{name}"] = float(sel.sum() / max(1, ok.sum()))
        out[f"vfrac_core_{name}"] = float(
            vol[sel & core].sum() / max(vol[core & ok].sum(), 1e-30))
        out[f"n_{name}"] = int(sel.sum())

    sel = ok & (mob >= MOB_YIELD)
    out["n_el"] = int(ok.sum())
    if sel.sum() == 0:
        out.update(ncomp=0, big_vfrac=0.0, big_n=0, band_fill=float("nan"),
                   band_thick_el=float("nan"), ncol=0, cells_per_col=float("nan"),
                   xmax_surface=float("nan"), xspan=0.0, zspan=0.0,
                   touches_surface=False, xmax_any=float("nan"),
                   zmin_any=float("nan"))
        return out, np.zeros(len(cells), dtype=bool)

    comps, _ = _components(ii, kk, sel, nx1, nz1)
    comps.sort(key=lambda c: vol[c].sum(), reverse=True)
    big = comps[0]
    bigsel = np.zeros(len(cells), dtype=bool)
    bigsel[big] = True

    dx = np.diff(x)
    dz = np.diff(z)
    area = dx[ii] * dz[kk]                        # the x-z plane footprint
    a_c = float(area[big].sum())
    i0, i1 = int(ii[big].min()), int(ii[big].max())
    k0, k1 = int(kk[big].min()), int(kk[big].max())
    a_bb = float((x[i1 + 1] - x[i0]) * (z[k1 + 1] - z[k0]))
    cols = np.unique(ii[big])
    per_col = np.array([(ii[big] == c).sum() for c in cols], dtype=float)
    # an equivalent band thickness: area / the length of the band's spine,
    # taken as the number of occupied columns times their mean width.
    spine = float(dx[cols].sum())
    out.update(ncomp=len(comps),
               big_vfrac=float(vol[big].sum() / vol[ok].sum()),
               big_n=int(len(big)),
               band_fill=a_c / a_bb if a_bb > 0 else float("nan"),
               band_thick_el=(a_c / spine) / HP_h0_of(z) if spine > 0 else float("nan"),
               ncol=int(len(cols)), cells_per_col=float(per_col.mean()),
               xspan=float(x[i1 + 1] - x[i0]), zspan=float(z[k1 + 1] - z[k0]))
    surf = big[kk[big] == nz1 - 1]                # elements touching z = 0
    out["touches_surface"] = bool(len(surf))
    out["xmax_surface"] = (float(np.abs(cen[surf, 0]).max()) if len(surf)
                           else float("nan"))
    out["xmax_any"] = float(np.abs(cen[big, 0]).max())
    out["zmin_any"] = float(cen[big, 2].min())
    return out, bigsel


def HP_h0_of(z):
    """The fine element size of this grid, read off the mesh rather than passed
    in, so the thickness metric cannot be given the wrong h0."""
    return float(np.min(np.diff(z)))


def heave_metrics(nodes, sets, w, du, ds_ref):
    """THE ISOCHORIC SPAN RATIO chi -- an exact volume balance, and the single
    most decisive number this script produces.

    At psi = 0 (rho_bar = 0) plastic flow is EXACTLY isochoric, so a footing
    that penetrates through a collapse MECHANISM has to push an equal volume of
    soil out somewhere, and on a strip footing "somewhere" is upward, beside the
    footing -- the Prandtl passive wedge.  That statement can be made exact
    rather than qualitative, because this domain's boundary conditions close the
    volume budget:

        dV = integral over the boundary of u.n dA

    and on this box the BASE is fixed (u = 0), the SIDES are rollers with
    u_x = 0 at x = +-XLIM, and u_y = 0 on EVERY node (that is what makes it
    plane strain).  Every one of those contributes u.n = 0, so the ONLY volume
    flux is through the free surface:

        dV = integral over the top of u_z dA = sum_n w_n du_z,n

    where `w` is the CONSISTENT surface weight vector -- the same
    `integral N_a dA` that the surcharge is built from, which is exactly the
    right operator for integrating an FE field over that surface (and is why the
    negative Q8 corner weights must NOT be replaced by tributary areas here
    either).  The rigid footing pushes its FOOTPRINT down by exactly ds, so
    splitting the surface integral at the footing edge gives
    dV = -A_foot ds + V_up with V_up the volume that came back up outside it, and

        chi = V_up / (A_foot ds) = 1 + dV / (A_foot ds),   A_foot = B * THICK.

    NOTE the split is done through `dV`, NOT by summing `w` over the footing
    NODE SET.  The first cut did the latter and it is wrong: the nodes on the
    footing edge collect `integral N dA` from the element OUTSIDE the footing
    too, so `sum(w[footing])` came out 1.1667 m^2 against a true footprint of
    1.0000 m^2 -- a 17 % error straight into the headline number.  A_foot here
    is the GEOMETRIC constant, and `dV` is exact.

        chi -> 1  the increment is ISOCHORIC: everything the footing displaces
                  comes back up beside it.  The admissible velocity fields SPAN
                  a mechanism.
        chi -> 0  the increment is swallowed as elastic volumetric COMPACTION.
                  The element cannot deliver an isochoric mechanism -- which is
                  precisely the span failure note 81 section 5.3 reframed the
                  question around.

    It needs no threshold, no connectivity, no mesh-size weighting, and it is
    immune to the artefact that killed the first cut of this measure: a plain
    `max(u_z)` beside the footing put its argmax at x = 30 m, the far BOUNDARY,
    where the box breathes elastically -- an integral over the surface cannot be
    fooled that way, and the far-field contribution is weighted by how little
    actually moves there.

    WHAT WOULD HAVE TO BE TRUE FOR chi TO BE ~1 WHILE THE MECHANISM IS NOT
    FORMING?  chi tests SPAN, not LOCALISATION: a completely DIFFUSE isochoric
    plastic flow also returns chi ~ 1.  So chi and the localisation metrics
    divide the work and neither substitutes for the other -- chi answers "can
    this element deliver an isochoric mechanism at all", `localisation_metrics`
    answers "is what it delivers a slip surface".  The converse combination --
    a large yielded set with chi ~ 0 -- is the diagnostic signature of LOCKING,
    and it is unambiguous.
    """
    top = sets["top"]
    xt, uz, wt = nodes[top, 0], du[top, 2], w[top]
    beside = np.abs(xt) > HP.B_FOOT / 2.0 + 1e-9
    a_foot = HP.B_FOOT * HP.THICK                  # GEOMETRIC, not from w
    out = dict(chi=float("nan"), heave_max=float("nan"), heave_x=float("nan"),
               heave_far=float("nan"), a_foot=a_foot,
               dV_rel=float("nan"), heave_uplift_frac=float("nan"))
    if ds_ref <= 0:
        return out
    out["dV_rel"] = float((wt * uz).sum() / (ds_ref * a_foot))
    out["chi"] = 1.0 + out["dV_rel"]
    j = int(np.argmax(uz[beside]))
    out["heave_max"] = float(uz[beside][j] / ds_ref)
    out["heave_x"] = float(xt[beside][j])
    far = beside & (np.abs(xt) > 20.0)             # the boundary-artefact monitor
    out["heave_far"] = float(uz[far].max() / ds_ref) if far.any() else float("nan")
    out["heave_uplift_frac"] = float((uz[beside] > 0).mean())
    return out


# ---------------------------------------------------------------------------
# tangent conditioning
# ---------------------------------------------------------------------------
def sample_tangent(setup, nsub_seen):
    """Assemble the tangent DENSELY, measure it, and put the sparse solver back.

    Why not simply run the whole leg on `FullGeneral`: measured, an H20 leg at
    h0 = 1.0 (2800 free DOF) steps at ~20 s/step under FullGeneral against a
    fraction of a second under PARDISO, because every Newton ITERATION pays a
    dense factorisation -- a 250-step leg becomes 80+ minutes and the
    conditioning trajectory costs more than every other diagnostic combined.

    So the leg SOLVES with PARDISO throughout and the dense matrix is assembled
    only at the sampling points: `LoadControl(0.0)` + `algorithm Linear` forms
    and factors K at the current state and commits a ZERO increment, which
    leaves displacements, stresses and pseudo-time untouched (it is
    `h20_prandtl.py`'s `leg_modes` idiom).  One dense factorisation per sample
    instead of one per iteration.
    """
    try:
        setup("FullGeneral")
        ops.test("NormUnbalance", 1.0e30, 1, 0)
        ops.algorithm("Linear")
        ops.integrator("LoadControl", 0.0)
        ops.analyze(1)
        th = tangent_health(nsub_seen)
    except Exception as exc:                              # pragma: no cover
        log(f"    [diag 3] tangent sample failed: {exc}")
        th = None
    finally:
        setup("sparse")
    return th


def tangent_health(nsub_seen):
    """SVD + symmetric-part spectrum of the assembled, restrained tangent.

    The non-associated tangent is UNSYMMETRIC (rho_bar != rho), so `eigvals`
    alone answers the wrong question: conditioning is a SINGULAR-value property.
    Both are taken.  `lam_min_sym` going negative is the loss of positive
    definiteness that a genuine limit point / localisation produces; sigma_min
    collapsing while lam_min_sym stays positive is ill-conditioning without a
    mechanism.
    """
    a = np.asarray(ops.printA("-ret"), dtype=float)
    n = int(round(math.sqrt(a.size)))
    if n * n != a.size:
        return None
    K = a.reshape(n, n)
    scale = float(np.trace(K)) / n
    sv = np.linalg.svd(K, compute_uv=False)
    S = 0.5 * (K + K.T)
    ev = np.linalg.eigvalsh(S)
    return dict(n=n, scale=scale,
                s_min=float(sv[-1]), s_max=float(sv[0]),
                cond=float(sv[0] / sv[-1]) if sv[-1] > 0 else float("inf"),
                s_min_rel=float(sv[-1] / scale),
                lam_min_sym_rel=float(ev[0] / scale),
                n_neg_sym=int((ev < 0).sum()), nsub=nsub_seen)


# ---------------------------------------------------------------------------
# one leg
# ---------------------------------------------------------------------------
def run(args):
    if ops is None:
        raise SystemExit("H20_NO_ENGINE is set; a leg needs the engine")
    spec = ELEMS[args.elem]
    order, form = spec["order"], spec["form"]
    base, floor, dmax = args.ds_base, args.floor, args.ds_max

    alpha = HP.alpha_from_phi_txc(HP.PHI_TXC)
    rho = math.sqrt(2.0) * alpha
    mc = HP.mc_from_cone(alpha)
    nq = HP.n_q(mc["ps"])
    q_exact = HP.Q0 * nq
    k0 = HP.NU / (1.0 - HP.NU)
    m0 = (1.0 - k0) / (math.sqrt(3.0) * alpha * (1.0 + 2.0 * k0))
    assoc = args.assoc

    stamp = ops.ladrunoBuild()
    log(f"[control 0] engine {ops.__file__}")
    log(f"[control 0] ladrunoBuild() = {stamp}")

    nodes, cells, sets, x, y, z = HP.strip_mesh(args.h0, order)
    w, nfaces = HP.consistent_surcharge(nodes, cells, order)
    HP.verify_surcharge(nodes, cells, w, nfaces, order)
    cen, vol, ii, kk = element_geometry(nodes, cells, x, z)

    tag = f"{args.elem}_h{str(args.h0).replace('0.', '')}" \
          f"{'_assoc' if assoc else ''}{args.suffix}"
    log(f"=== leg {tag}: {spec['label']}, h0 = {args.h0} m "
        f"({int(round(2.0 / args.h0))} el across B), {len(nodes)} nodes / "
        f"{3 * len(nodes)} DOF, {len(cells)} elements")
    log(f"    ORACLE q_u = q0*N_q = {q_exact:.3f} kPa (phi_ps = {mc['ps']:.3f} deg); "
        f"{'ASSOCIATED (control)' if assoc else 'non-associated (gate)'}")
    log(f"    ALLOWANCES (pinned, reported): ds base/floor/max = "
        f"{base * 1e3:.4g}/{floor * 1e3:.4g}/{dmax * 1e3:.4g} mm "
        f"(= {math.log2(base / floor):.2f} halvings), subdivision budget = "
        f"{args.budget}, s/B target = {args.sfrac}, wall cap = {args.tmax:.0f} s")
    log(f"    initial 1-D state m0 = {m0:.4f} of yield"
        + ("  (SAFE)" if m0 < 0.8 else "  *** VOID ***"))
    assert m0 < 0.8

    # --- the per-leg ZERO of the mechanism scale, on THIS mesh and element ---
    chi_el = elastic_chi(args.h0, order, form, HP.NU)["chi"]
    log(f"    [scale] purely ELASTIC chi on this mesh at nu = {HP.NU} is "
        f"{chi_el:+.5f}; a fully isochoric mechanism is chi = +1.  The "
        f"mechanism index below is M = (chi - {chi_el:+.5f}) / "
        f"({1.0 - chi_el:.5f}), so M = 0 is 'no mechanism, elastic' and M = 1 "
        f"is 'fully formed'.")

    HP.build_model(nodes, cells, sets, w, form, assoc)
    tol = 1.0e-5 * max(HP.Q0 * float(w.sum()), 1.0)
    HP.surcharge_step(nodes, sets, w, tol)

    # --- control 1: the 1-D elastic stress patch, on the PLASTIC model -------
    # The state after the surcharge step is still elastic (m0 = 0.29), so the
    # exact 1-D answer applies and the check costs nothing extra.  It is the
    # one that sees a load-DISTRIBUTION error; the resultant identity in
    # surcharge_step() is structurally blind to it (note 81 control 2's
    # negative control: 190 % into sigma_zz at +0.0000000 % on the resultant).
    e_zz = e_xx = 0.0
    for e in range(1, len(cells) + 1):
        st = np.asarray(ops.eleResponse(e, "stress"), dtype=float).reshape(-1, 6)
        e_zz = max(e_zz, float(np.abs(st[:, 2] / HP.Q0 + 1.0).max()))
        e_xx = max(e_xx, float(np.abs(st[:, 0] / (k0 * HP.Q0) + 1.0).max()))
    patch_err = max(e_zz, e_xx)
    log(f"    [control 1] 1-D elastic stress patch (every GP): "
        f"max|szz/-q0-1| = {e_zz:.2e}, max|sxx/-K0q0-1| = {e_xx:.2e}  -> "
        f"{'PASS' if patch_err < 1e-6 else '*** FAIL ***'}")
    if patch_err >= 1e-6:
        raise SystemExit("patch test failed; every number downstream is void")

    # --- control M0: the mobilisation field agrees with the closed form ------
    # The 1-D state sigma_zz = -q0, sigma_xx = sigma_yy = -K0 q0 gives
    #   ||dev sigma|| = q0 (1-K0) sqrt(6)/3,   I1 = -q0 (1+2 K0)
    # so the EXACT initial mobilisation is
    #   m0* = q0 (1-K0) sqrt(6)/3 / (sqrt(2/3) SY + rho q0 (1+2 K0)).
    # h20_prandtl.py's `m0` DROPS the sqrt(2/3) SY apex term, which is a fine
    # approximation for a "is the leg void?" screen but is 2.9 % off here -- and
    # this control caught exactly that on the first run, which is the reason it
    # is written against the exact form and not the harness's.
    mob0_mean, mob0_max, ngp = mobilisation(len(cells), rho)
    m0_exact = (HP.Q0 * (1.0 - k0) * math.sqrt(6.0) / 3.0
                / (SQ23 * HP.SY + rho * HP.Q0 * (1.0 + 2.0 * k0)))
    m0_err = float(np.nanmax(np.abs(mob0_mean - m0_exact)))
    log(f"    [control M0] mobilisation field at the surcharge step "
        f"({ngp} GP/element): max |mob - m0*| = {m0_err:.2e} against the exact "
        f"closed form m0* = {m0_exact:.9f} (the harness's SY-free screen prints "
        f"{m0:.6f})  -> {'PASS' if m0_err < 1e-9 else '*** FAIL ***'}")
    if m0_err >= 1e-9:
        raise SystemExit("the mobilisation formula disagrees with the 1-D state")
    log(f"    [control M0] elements at mob >= {MOB_YIELD} before the push: "
        f"{int((mob0_mean >= MOB_YIELD).sum())} of {len(cells)}  (must be 0)")

    # ---------------- the push ---------------------------------------------
    foot = [int(n) + 1 for n in sets["footing"]]
    r_corr = HP.Q0 * float(w[sets["footing"]].sum())
    area = HP.B_FOOT * HP.THICK
    uz0 = ops.nodeDisp(foot[0], 3)
    smax = args.sfrac * HP.B_FOOT
    ops.loadConst("-time", uz0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for t in foot:
        ops.sp(t, 3, 1.0)

    def setup(sysname):
        ops.wipeAnalysis()
        ops.constraints("Transformation")
        ops.numberer("Plain" if sysname == "FullGeneral" else "RCM")
        if sysname == "FullGeneral":
            ops.system("FullGeneral")
        else:
            try:
                ops.system("Pardiso")
            except Exception:
                ops.system("UmfPack")
        ops.analysis("Static")

    setup("sparse")
    ladder = ([("KrylovNewton", tol, 25, 0),
               ("NewtonLineSearch", tol, 40, 0),
               ("KrylovNewton", 10.0 * tol, 60, 1)]
              if args.ladder == "quad" else
              [("Newton", tol, 25, 0),
               ("NewtonLineSearch", tol, 40, 0),
               ("KrylovNewton", 10.0 * tol, 60, 1)])
    if args.ladder == "wide":
        # A genuinely DIFFERENT algorithm set at the SAME tolerances -- which is
        # what note 81 §4.3's tolerance rung was not.  Deliberately TRIMMED: the
        # first cut appended ModifiedNewton(200) + Broyden(60) + BFGS(60), which
        # made every subdivision cost 445 iterations against the reference
        # ladder's 125, and the leg spent 20 minutes of wall clock without
        # converging a single step.  An experiment nobody can afford to run is
        # not a stronger path algorithm either.
        ladder = [("KrylovNewton", tol, 25, 0),
                  ("NewtonLineSearch", tol, 40, 0),
                  ("ModifiedNewton", tol, 60, 0),
                  ("Newton", tol, 40, 0),
                  ("KrylovNewton", 10.0 * tol, 60, 1)]

    path = os.path.join(HERE, f"qpd_{tag}.csv")
    fh = open(path, "w", newline="")
    wr = csv.writer(fh)
    wr.writerow(["s_m", "s_over_B", "q_kPa", "ds_mm", "relaxed", "wall_s"])
    rows, ds, good, nfail, nrelax, nsub = [], base, 0, 0, 0, 0
    mode, verdict = "TARGET", "reached the target settlement"
    hist = deque(maxlen=400)                      # (s, disp) for the increment
    disp_push0 = np.array([[ops.nodeDisp(i, c) for c in (1, 2, 3)]
                           for i in range(1, len(nodes) + 1)])
    cond_log, cond_on = [], False
    t0 = time.time()
    while True:
        s_now = uz0 - ops.getTime()
        if s_now >= smax - 1e-12:
            break
        if time.time() - t0 > args.tmax:
            mode, verdict = "WALL", f"wall-clock cap at s/B = {s_now / HP.B_FOOT:.5f}"
            break
        ds = min(ds, smax - s_now)
        if args.cond and not cond_on and ds < args.cond_at:
            cond_on = True
            log(f"    [diag 3] ds fell to {ds * 1e3:.6g} mm (trigger "
                f"{args.cond_at * 1e3:.6g} mm) at s/B = {s_now / HP.B_FOOT:.5f}: "
                f"arming the tangent sampler (dense assembly at the sampling "
                f"points only; PARDISO still solves the leg)")
        ops.integrator("LoadControl", -ds)
        ok, relaxed = False, 0
        for algo, tl, it, rl in ladder:
            ops.test("NormUnbalance", tl, it, 0)
            ops.algorithm(algo)
            if ops.analyze(1) == 0:
                ok, relaxed = True, rl
                break
            nfail += 1
        if not ok:
            good, nsub = 0, nsub + 1
            ds *= 0.5
            if nsub > args.budget:
                mode = "BUDGET"
                verdict = (f"subdivision budget of {args.budget} spent at "
                           f"s/B = {s_now / HP.B_FOOT:.5f}")
                break
            if ds < floor:
                mode = "FLOOR"
                verdict = (f"step collapsed to the {floor * 1e3:.4g} mm floor at "
                           f"s/B = {s_now / HP.B_FOOT:.5f} (every ladder rung "
                           f"failed at ds = {ds * 1e3:.6g} mm)")
                break
            continue
        nrelax += relaxed
        good += 1
        if good >= HP.GROW_AFTER and ds < dmax:
            ds, good = min(2 * ds, dmax), 0
        ops.reactions()
        q = (-sum(ops.nodeReaction(t, 3) for t in foot) + r_corr) / area
        s = uz0 - ops.getTime()
        rows.append((s, s / HP.B_FOOT, q, ds * 1e3, relaxed, time.time() - t0))
        wr.writerow([f"{v:.9g}" for v in rows[-1]])
        fh.flush()
        hist.append((s, np.array([[ops.nodeDisp(i, c) for c in (1, 2, 3)]
                                  for i in range(1, len(nodes) + 1)])))
        if cond_on and len(rows) % args.cond_every == 0:
            th = sample_tangent(setup, nsub)
            if th:
                th["s_over_B"] = s / HP.B_FOOT
                th["q"] = q
                th["ds_mm"] = ds * 1e3
                cond_log.append(th)
                log(f"    [diag 3] s/B {th['s_over_B']:.5f}  q {q:8.2f}  "
                    f"ds {ds * 1e3:.5g} mm  |  sigma_min/scale "
                    f"{th['s_min_rel']:.3e}  cond {th['cond']:.3e}  "
                    f"lam_min(sym)/scale {th['lam_min_sym_rel']:.3e}  "
                    f"n_neg(sym) {th['n_neg_sym']}")
        if len(rows) > HP.STALL_WINDOW and \
                rows[-1][0] - rows[-1 - HP.STALL_WINDOW][0] < HP.STALL_ADVANCE * smax:
            mode = "STALL"
            verdict = (f"stalled: < {100 * HP.STALL_ADVANCE:.2f} % of the target "
                       f"in {HP.STALL_WINDOW} steps at s/B = {s / HP.B_FOOT:.5f}")
            break
    fh.close()
    wall = time.time() - t0
    assert len(rows) > 8, f"only {len(rows)} converged steps ({verdict})"

    # --- diagnostic 3, the sample that matters most: the tangent AT the wall --
    # The in-loop sampler only fires on CONVERGED steps, and a leg can pass from
    # its trigger to its termination without converging another one -- measured
    # on the first h20uri run, which triggered at s/B = 0.01121 and then walled.
    # So the last converged state is re-formed and sampled unconditionally here.
    # `LoadControl(0.0)` + `algorithm Linear` assembles and factors K at the
    # current state and commits a ZERO increment, so the state is untouched;
    # it is `leg_modes`' idiom from h20_prandtl.py.
    if args.cond:
        th = sample_tangent(setup, nsub)
        if th:
            th.update(s_over_B=float(rows[-1][1]), q=float(rows[-1][2]),
                      ds_mm=float(rows[-1][3]), at_wall=True)
            cond_log.append(th)
            log(f"    [diag 3] AT THE WALL, s/B {th['s_over_B']:.5f}: "
                f"sigma_min/scale {th['s_min_rel']:.3e}, cond {th['cond']:.3e}, "
                f"lam_min(sym)/scale {th['lam_min_sym_rel']:.3e}, negative "
                f"eigenvalues of the symmetric part: {th['n_neg_sym']} of "
                f"{th['n']}")

    a = np.array([r[:4] for r in rows])
    s, q, dsm = a[:, 0], a[:, 2], a[:, 3]
    peak = np.maximum.accumulate(q)
    bad = np.where(q < (1.0 - DROP_TRUNCATE) * peak)[0]
    if len(bad) and bad[0] > 8:
        mode = "TRUNCATED"
        verdict += f"; curve TRUNCATED ({len(q) - bad[0]} steps cut)"
        s, q, dsm = s[:bad[0]], q[:bad[0]], dsm[:bad[0]]

    msk = s >= 0.9 * s[-1]
    t_last = float(np.polyfit(s[msk], q[msk], 1)[0]) if msk.sum() > 2 else float("nan")
    n0 = max(4, len(s) // 50)
    t_init = float(np.polyfit(s[:n0], q[:n0], 1)[0])
    plateau = bool(abs(t_last) < PLATEAU_FRAC * abs(t_init))
    qmax = float(q.max())
    floor_mm = floor * 1e3
    ds_tail_min = float(dsm[msk].min()) if msk.sum() else float(dsm[-1])
    headroom = ds_tail_min / floor_mm
    free = bool(headroom >= FREE_ADVANCE_FLOOR_FACTOR)
    capacity = bool(plateau and free and mode in _CAPACITY_MODES)

    # ---------------- diagnostic 1, at the LAST CONVERGED STEP --------------
    # StaticAnalysis::analyze() reverts the domain to the last commit when the
    # algorithm fails, so the state standing here IS the last converged step --
    # no bookkeeping needed, and the failed attempt has left no trace.
    mob_mean, mob_max, ngp = mobilisation(len(cells), rho)
    loc, bigsel = localisation_metrics(nodes, cells, cen, vol, ii, kk,
                                       mob_mean, x, z, tag)
    s_end, disp_end = hist[-1]
    # A FINITE window, not the last two steps: near the wall the last step can be
    # 1e-7 m and its increment is round-off (see heave_metrics' docstring).  Take
    # the LATEST stored state that is at least `window` behind the end, so the
    # increment is as close to the wall as it can be while still spanning a
    # measurable settlement; if the whole stored history is shorter than that
    # (every step tiny), fall back to the earliest stored state, which is the
    # largest window available.  The first cut of this used `argmin |s - want|`
    # with a 1 %-of-run window and silently picked the LAST step itself whenever
    # one step was bigger than the window -- a zero-width increment and a NaN.
    window = max(0.02 * s_end, 50.0 * floor)
    cand = [i for i, (sv, _) in enumerate(hist) if s_end - sv >= window]
    j = cand[-1] if cand else 0
    s_ref, disp_ref = hist[j]
    du = disp_end - disp_ref
    hv = heave_metrics(nodes, sets, w, du, s_end - s_ref)
    hv["ds_window_m"] = float(s_end - s_ref)
    hv["ds_window_steps"] = int(len(hist) - j)
    hv["window_clipped"] = bool(not cand)
    # ...and the same measure over the WHOLE push, which answers a different
    # question: the incremental one asks "is a mechanism operating AT the wall",
    # the total one asks "did one ever form at all".
    tot = heave_metrics(nodes, sets, w, disp_end - disp_push0, s_end)
    hv["chi_total"] = tot["chi"]
    hv["heave_max_total"] = tot["heave_max"]
    hv["heave_x_total"] = tot["heave_x"]

    log("")
    log(f"    --- termination -------------------------------------------------")
    log(f"    {len(rows)} steps, {nfail} failed attempts, {nsub}/{args.budget} "
        f"subdivisions, {nrelax} relaxed, {wall:.0f} s")
    log(f"    MODE = {mode}   [{verdict}]")
    log(f"    q_max = {qmax:.2f} kPa = {qmax / q_exact:.4f} of exact "
        f"({q_exact:.2f} kPa); end s/B = {s[-1] / HP.B_FOOT:.5f} of "
        f"{args.sfrac}")
    log(f"    tail dq/ds = {t_last:.1f} kPa/m = {100 * t_last / t_init:.3f} % of "
        f"initial -> {'PLATEAU' if plateau else 'STILL HARDENING'}")
    log(f"    tail-min step {ds_tail_min:.6g} mm = {headroom:.1f}x the "
        f"{floor_mm:.6g} mm floor -> free-advance {'yes' if free else 'NO'}")
    log(f"    CAPACITY = {'YES' if capacity else 'NO'}"
        + ("" if capacity else "  -> this number is a CONTROLLER ALLOWANCE "
                               f"named {mode}, not an element capacity"))
    log("")
    log(f"    --- diagnostic 1: mechanism formation at termination -------------")
    log(f"    yielded (mob >= {MOB_YIELD}) : {loc['n_y90']}/{loc['n_el']} elements, "
        f"{100 * loc['vfrac_y90']:.2f} % BY VOLUME "
        f"({100 * loc['cfrac_y90']:.2f} % by count), "
        f"{100 * loc['vfrac_core_y90']:.2f} % of the core "
        f"(|x|<={CORE_X}, z>={CORE_Z})")
    log(f"    at yield  (mob >= {MOB_HARD}) : {loc['n_y99']}/{loc['n_el']}, "
        f"{100 * loc['vfrac_y99']:.2f} % by volume")
    log(f"    connectivity: {loc['ncomp']} component(s); the largest holds "
        f"{100 * loc['big_vfrac']:.2f} % of the volume ({loc['big_n']} elements)")
    log(f"    band-ness of the largest: fills {100 * loc['band_fill']:.1f} % of its "
        f"own bounding box ({loc['xspan']:.2f} x {abs(loc['zspan']):.2f} m), "
        f"{loc['ncol']} columns x {loc['cells_per_col']:.2f} cells/column, "
        f"equivalent thickness {loc['band_thick_el']:.2f} elements")
    x_prandtl = (HP.B_FOOT / 2.0
                 * math.exp(math.pi * math.tan(math.radians(mc["ps"])))
                 * math.tan(math.pi / 4 + math.radians(mc["ps"]) / 2))
    log(f"    extent: reaches |x| = {loc['xmax_any']:.2f} m and z = "
        f"{loc['zmin_any']:.2f} m; touches the free surface: "
        f"{'YES, out to |x| = %.2f m' % loc['xmax_surface'] if loc['touches_surface'] else 'NO'}"
        f"  (a complete Prandtl fan surfaces near |x| ~ {x_prandtl:.2f} m)")
    hv["chi_el"] = chi_el
    hv["M"] = (hv["chi"] - chi_el) / (1.0 - chi_el)
    hv["M_total"] = (hv["chi_total"] - chi_el) / (1.0 - chi_el)
    log(f"    ISOCHORIC SPAN RATIO chi = {hv['chi']:+.4f} AT THE WALL   ->   "
        f"MECHANISM INDEX M = {hv['M']:+.4f}   "
        f"(0 = elastic {chi_el:+.4f}, 1 = fully isochoric)")
    log(f"           over the WHOLE push: chi = {hv['chi_total']:+.4f}, "
        f"M = {hv['M_total']:+.4f}")
    log(f"           increment over {hv['ds_window_m']:.3e} m = "
        f"{hv['ds_window_steps']} steps"
        f"{', WINDOW CLIPPED to the stored history' if hv['window_clipped'] else ''}")
    log(f"           max incremental u_z beside the footing = "
        f"{hv['heave_max']:+.4f} x ds at x = {hv['heave_x']:.2f} m; at the far "
        f"boundary {hv['heave_far']:+.4f} x ds  (a peak AT the boundary is the "
        f"box breathing, not a mechanism -- chi is immune to it)")

    np.savez(os.path.join(HERE, f"qpd_{tag}_field.npz"),
             s=s, q=q, ds_mm=dsm, cen=cen, vol=vol, ii=ii, kk=kk,
             mob_mean=mob_mean, mob_max=mob_max, mob0=mob0_mean, big=bigsel,
             nodes=nodes, disp=disp_end, du=du, x=x, z=z,
             q_exact=q_exact, qmax=qmax, m0=m0, phi_ps=mc["ps"])

    res = dict(tag=tag, elem=args.elem, label=spec["label"], form=form,
               order=order, h0=args.h0, assoc=assoc, build=stamp,
               nel=int(len(cells)), ndof=int(3 * len(nodes)), ngp=int(ngp),
               qmax=qmax, q_exact=q_exact, ratio=qmax / q_exact,
               mode=mode, verdict=verdict, plateau=plateau, free=free,
               capacity=capacity, tail_pct=100 * t_last / t_init,
               t_init=t_init, t_last=t_last,
               s_end_over_B=float(s[-1]) / HP.B_FOOT, s_target=args.sfrac,
               ds_base_mm=base * 1e3, ds_floor_mm=floor_mm,
               ds_end_mm=float(dsm[-1]), ds_tail_min_mm=ds_tail_min,
               headroom=headroom, halvings=math.log2(base / floor),
               nsub=nsub, budget=args.budget,
               budget_used_frac=nsub / float(args.budget),
               steps=len(rows), nfail=nfail, nrelax=nrelax, wall_s=wall,
               patch_err=patch_err, m0_err=m0_err, ladder=args.ladder,
               cond=cond_log, **loc, **hv)
    with open(os.path.join(HERE, f"qpd_{tag}.json"), "w") as f:
        json.dump(res, f, indent=1, default=float)
    return res


def elastic_chi(h0, order, form, nu, nsteps=8, ds=1.0e-4):
    """chi for a LINEAR ELASTIC push of the same footing into the same mesh at
    Poisson ratio `nu`.  Two uses:

      * at nu -> 1/2 and nu = 0 it CALIBRATES chi (see `calibrate_chi`);
      * at the OPERATING nu it is the per-leg ZERO of the mechanism scale.  This
        is the reason chi is interpretable at all: at nu = 0.45 the elastic body
        is already fairly (not perfectly) incompressible, so "no mechanism" does
        NOT mean chi = 0, it means chi = chi_el, measured on this exact mesh and
        element.  Reporting chi against 0 instead of against chi_el would read a
        purely elastic increment as a partial mechanism.
    """
    nu_save = HP.NU
    try:
        HP.NU = nu
        nodes, cells, sets, x, y, z = HP.strip_mesh(h0, order)
        w, _ = HP.consistent_surcharge(nodes, cells, order)
        HP.build_model(nodes, cells, sets, w, form, False, elastic=True)
        tol = 1e-8 * HP.Q0 * float(w.sum())

        def sysup():
            ops.constraints("Transformation")
            ops.numberer("RCM")
            try:
                ops.system("Pardiso")
            except Exception:
                ops.system("UmfPack")
            ops.test("NormUnbalance", tol, 25, 0)
            ops.algorithm("Newton")

        sysup()
        ops.integrator("LoadControl", 1.0)
        ops.analysis("Static")
        ops.analyze(1)
        foot = [int(n) + 1 for n in sets["footing"]]
        uz0 = ops.nodeDisp(foot[0], 3)
        d0 = np.array([[ops.nodeDisp(i, c) for c in (1, 2, 3)]
                       for i in range(1, len(nodes) + 1)])
        ops.loadConst("-time", uz0)
        ops.timeSeries("Linear", 2)
        ops.pattern("Plain", 2, 2)
        for t in foot:
            ops.sp(t, 3, 1.0)
        ops.wipeAnalysis()
        sysup()
        ops.integrator("LoadControl", -ds)
        ops.analysis("Static")
        ops.analyze(nsteps)
        d1 = np.array([[ops.nodeDisp(i, c) for c in (1, 2, 3)]
                       for i in range(1, len(nodes) + 1)])
        return heave_metrics(nodes, sets, w, d1 - d0, nsteps * ds)
    finally:
        HP.NU = nu_save


def calibrate_chi(elem, h0, nsteps=8):
    """CONTROL CHI -- the two-sided calibration that makes chi worth quoting.

    chi is asserted above to read 1 for an isochoric increment and 0 for one
    swallowed by volumetric compression.  That is a claim about the measure, and
    it is checkable WITHOUT any plasticity: a LINEAR ELASTIC body at nu -> 1/2 is
    exactly incompressible, so pushing the footing into it must give chi -> 1,
    and at nu = 0 it is maximally compressible and must give chi well below 1.
    Both ends are run on the same mesh, the same element and the same push as
    the real legs, so a chi that fails here is a bug in the measure and not a
    property of the soil.

    WHAT WOULD HAVE TO BE TRUE FOR THIS TO PASS WHILE CHI IS WRONG?  chi could
    still be wrong by a CONSTANT factor and hit 1.0 at nu = 1/2 only by luck --
    which is why the nu = 0 end is run as well, and why the intermediate nu are
    printed: the sequence has to be monotone in nu and to bracket the operating
    nu = 0.45 sensibly, not just touch one endpoint.
    """
    spec = ELEMS[elem]
    log(f"=== control CHI: elastic calibration of the isochoric span ratio, "
        f"{spec['label']}, h0 = {h0}")
    rows = []
    for nu in (0.0, 0.25, 0.45, 0.49, 0.4999):
        hv = elastic_chi(h0, spec["order"], spec["form"], nu, nsteps)
        rows.append((nu, hv["chi"], hv["heave_max"], hv["heave_x"]))
        log(f"    nu = {nu:<7} chi = {hv['chi']:+.5f}   "
            f"(max heave {hv['heave_max']:+.4f} x ds at x = {hv['heave_x']:.2f} m)")
    chi_inc = rows[-1][1]
    chi_0 = rows[0][1]
    ok = abs(chi_inc - 1.0) < 0.02 and chi_0 < 0.5
    log(f"    -> chi(nu=0.4999) = {chi_inc:+.5f} (must be 1.00 +/- 0.02), "
        f"chi(nu=0) = {chi_0:+.5f} (must be < 0.5), monotone in nu: "
        f"{all(rows[i][1] <= rows[i + 1][1] + 1e-9 for i in range(len(rows) - 1))}"
        f"  -> {'PASS' if ok else '*** FAIL ***'}")
    if not ok:
        raise SystemExit("the isochoric span ratio is not calibrated; do not "
                         "quote it")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--elem", default="h20uri", choices=sorted(ELEMS))
    ap.add_argument("--h0", type=float, default=1.0)
    ap.add_argument("--sfrac", type=float, default=0.15)
    ap.add_argument("--assoc", action="store_true")
    ap.add_argument("--budget", type=int, default=200)
    ap.add_argument("--ds-base", type=float, default=None)
    ap.add_argument("--floor", type=float, default=None)
    ap.add_argument("--ds-max", type=float, default=None)
    ap.add_argument("--ladder", default=None, choices=["quad", "linear", "wide"])
    ap.add_argument("--tmax", type=float, default=5400.0)
    ap.add_argument("--cond", action="store_true")
    ap.add_argument("--cond-every", type=int, default=4)
    ap.add_argument("--cond-at", type=float, default=None, metavar="M",
                    help="switch to FullGeneral once the step falls below this "
                         "(metres).  Default COND_TRIGGER x DS_BASE, which sits "
                         "very close to the wall; for a conditioning TRAJECTORY "
                         "set it near DS_MAX so the approach is sampled too.")
    ap.add_argument("--suffix", default="")
    ap.add_argument("--calib", action="store_true",
                    help="run control CHI (elastic calibration) and stop")
    args = ap.parse_args()

    if args.calib:
        calibrate_chi(args.elem, args.h0)
        return

    order = ELEMS[args.elem]["order"]
    b, f, m = HP.DS_LADDER[order]
    args.ds_base = b if args.ds_base is None else args.ds_base
    args.floor = f if args.floor is None else args.floor
    args.ds_max = m if args.ds_max is None else args.ds_max
    if args.ladder is None:
        args.ladder = "quad" if order == 2 else "linear"
    if args.cond_at is None:
        args.cond_at = COND_TRIGGER * args.ds_base

    r = run(args)
    print()
    print(f"{r['tag']:>24} {r['elem']:>8} mode={r['mode']:<9} "
          f"ratio={r['ratio']:.4f} tail={r['tail_pct']:.3f}% "
          f"ds/floor={r['headroom']:.1f} sub={r['nsub']}/{r['budget']} "
          f"CAP={'yes' if r['capacity'] else 'NO'}")


if __name__ == "__main__":
    main()
