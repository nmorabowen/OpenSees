"""
ADR-95 H4: rank of the fully-plastic element tangent (Prandtl-Reissner /
quadratic-limit-load-ceiling side experiment).

Hypothesis (pre-registered in 95_prandtl_reissner_quadratic_root_cause_plan.md,
worktree note): with perfect plasticity, the consistent tangent D_ep at a
plastic Gauss point has rank 5 (rank-1 deficiency along the flow direction,
non-associated DP => tangent unsymmetric but still rank 5). An element whose
GPs are ALL plastic therefore has element-tangent rank <= 6(rigid) + nGP*5
(deformation modes), so its number of near-zero singular values is
>= nDOF - 6 - nGP*5 (when that is positive).

PROTOCOL: single element (unit cube / unit tet), nDMaterial DruckerPrager
(perfect plasticity: Kinf=Ko=delta1=delta2=H=theta=0). A homogeneous
(traceless, pure-shear-like) proportional strain field u_i = lambda*eps_ij*x_j
is prescribed via sp() at EVERY node (corners AND midside nodes), so the
patch-test-exact field makes every Gauss point see the identical stress state
regardless of element order -- the whole element yields at one lambda. Ramp
lambda with Newton/LoadControl to an "elastic" state (lambda=1.0, comfortably
below yield ~1.33) and a "fully plastic" state (lambda=3.0, ~2.2x yield,
verified via eleResponse(tag,'material',gp,'ladrunoBranch')[0]==1 at every
GP). Then read the element tangent:
  - PREFERRED: eleResponse(tag, 'stiff') -- the element's own getTangentStiff()
    dumped directly (LadrunoBrick, LadrunoBrick20, BezierTet10 all support
    this; bypasses the domain solve entirely, so it is immune to the
    ill-conditioning below).
  - FALLBACK: remove the sp pattern + the 6 rigid-mode fix()es, reform the
    analysis with 0 free-dof BCs, run one analyze(1) at LoadControl(0.0) with
    a Linear algorithm (forms K but does not need to converge), and read
    ops.printA('-ret'). This calls a real LAPACK dense solve on a matrix that
    is EXACTLY multiply rank-deficient in the plastic state; that solve can
    produce an exact-zero pivot and poison the whole in-place buffer with
    NaN (observed on TenNodeTetrahedron plastic and the H20-uri 2x2x2 mesh
    plastic case below) -- reported as REFUSED rather than worked around,
    per the pre-registered protocol.

Two OpenSees gotchas hit while building this harness (not fixed here, C++ is
off-limits for this experiment):
  - `constraints('Plain')` silently zeroes any non-homogeneous (pattern-scaled)
    sp() value ("non-homogeneos constraint ... homogeneous constraint
    assumed") -- use `constraints('Transformation')` for prescribed-displacement
    ramping.
  - `Domain::removeSP_Constraint(node, dof, patternTag)` (the 3-arg overload
    used by `ops.remove('sp', node, dof, patternTag)`) looks up the SP inside
    the pattern's own iterator but then calls the DOMAIN-only 1-arg remover on
    the tag it found, which never touches the pattern's own container -- so a
    pattern-scoped sp() can never actually be removed this way, only via
    `ops.remove('loadPattern', patternTag)` (drops the whole pattern, which
    correctly detaches its owned SP_Constraints).

Run: python3.12 p1_plastic_rank.py   (from the worktree root, or anywhere --
adjusts sys.path to <worktree>/dist/bin itself).
"""
import os
import sys

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
WORKTREE = os.path.abspath(os.path.join(THIS_DIR, "..", "..", ".."))
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")
sys.path.insert(0, os.path.join(WORKTREE, "dist", "bin"))

import opensees as ops
import numpy as np

DP_ARGS = (1e4, 5e3, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
# K, G, sigma_y, rho, rho_bar, Kinf, Ko, delta1, delta2, H, theta
# -- perfect plasticity (no hardening), non-associated (rho=0.2 != rho_bar=0.0)
EPS0 = np.array([1e-5, -0.5e-5, -0.5e-5])  # traceless (I1=0) proportional shear path
DLAM = 0.05
LAM_ELASTIC = 1.0   # ~0.75x the measured yield lambda (~1.33)
LAM_PLASTIC = 3.0   # ~2.2x yield -- comfortably on the fully-plastic branch


def build_common(coords, elem_fn, fixed_dofs, matTag=1):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for n, c in coords.items():
        ops.node(n, *c)
    ops.nDMaterial('DruckerPrager', matTag, *DP_ARGS)
    elem_fn(matTag)
    for n, dofs in fixed_dofs.items():
        ops.fix(n, *dofs)
    fixed_set = set()
    for n, dofs in fixed_dofs.items():
        for i, f in enumerate(dofs):
            if f:
                fixed_set.add((n, i + 1))
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for n, c in coords.items():
        x, y, z = c
        target = (EPS0[0] * x, EPS0[1] * y, EPS0[2] * z)
        for dof, val in zip((1, 2, 3), target):
            if (n, dof) not in fixed_set:
                ops.sp(n, dof, val)
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1e-14, 30)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', DLAM)
    ops.analysis('Static')
    return fixed_set


def ramp_to(lam_target, dlam=DLAM):
    for _ in range(round(lam_target / dlam)):
        ops.analyze(1)


def gp_branches(eleTag, nGP):
    return [ops.eleResponse(eleTag, 'material', gp, 'ladrunoBranch')[0] for gp in range(1, nGP + 1)]


def stiff_via_eleresponse(eleTag, ndof, name='stiff'):
    K = np.array(ops.eleResponse(eleTag, name))
    if K.size != ndof * ndof:
        return None
    return K.reshape(ndof, ndof)


def stiff_via_printA(fixed_set, ndof, patternTag=1):
    """Fallback for elements with no direct 'stiff' response. See module
    docstring for why this can come back poisoned with NaN in the plastic
    (multiply rank-deficient) state -- caller must check via sv_report."""
    ops.remove('loadPattern', patternTag)
    for n, dof in fixed_set:
        ops.remove('sp', n, dof)
    ops.wipeAnalysis()
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1e-14, 1)
    ops.algorithm('Linear')
    ops.integrator('LoadControl', 0.0)
    ops.analysis('Static')
    ops.analyze(1)
    A = np.array(ops.printA('-ret'))
    if A.size != ndof * ndof:
        return None
    return A.reshape(ndof, ndof)


def sv_report(K):
    if K is None or np.isnan(K).any() or np.isinf(K).any():
        return None
    sv = np.sort(np.linalg.svd(K, compute_uv=False))[::-1]
    nz = int(np.sum(sv < 1e-8 * sv.max()))
    return {'sv': sv, 'nz': nz, 'smallest25': sv[-25:], 'asym': None}


# ---------------------------------------------------------------------------
# Node-coordinate generators (all use the SAME 3-2-1 support pattern: fix the
# origin fully, the next axis corner in the two transverse dofs, the third
# non-collinear corner in one dof -- exactly cancels the imposed field's own
# zero there, so no conflict with the ramped sp()).
# ---------------------------------------------------------------------------

BRICK8_COORDS = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
                  5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
BRICK8_FIX = {1: (1, 1, 1), 2: (0, 1, 1), 4: (0, 0, 1)}

_CORNER_OFF = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
               (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
_EDGE_V = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
           (0, 4), (1, 5), (2, 6), (3, 7)]
_LOCAL20 = list(_CORNER_OFF) + [
    tuple((_CORNER_OFF[a][i] + _CORNER_OFF[b][i]) / 2.0 for i in range(3))
    for a, b in _EDGE_V
]
BRICK20_COORDS = {i + 1: _LOCAL20[i] for i in range(20)}
BRICK20_FIX = {1: (1, 1, 1), 2: (0, 1, 1), 4: (0, 0, 1)}

TET10_COORDS = {1: (0, 0, 0), 2: (1, 0, 0), 3: (0, 1, 0), 4: (0, 0, 1),
                5: (0.5, 0, 0), 6: (0.5, 0.5, 0), 7: (0, 0.5, 0),
                8: (0, 0, 0.5), 9: (0, 0.5, 0.5), 10: (0.5, 0, 0.5)}
TET10_FIX = {1: (1, 1, 1), 2: (0, 1, 1), 3: (0, 0, 1)}


def gen_h20_mesh_2x2x2():
    node_tags = {}
    next_tag = [1]

    def get_node(coord):
        key = tuple(round(c * 2) for c in coord)
        if key not in node_tags:
            node_tags[key] = next_tag[0]
            next_tag[0] += 1
        return node_tags[key]

    coords_out, elems, etag = {}, [], 1
    for i0 in (0, 1):
        for j0 in (0, 1):
            for k0 in (0, 1):
                nlist = []
                for (ox, oy, oz) in _LOCAL20:
                    c = (i0 + ox, j0 + oy, k0 + oz)
                    nt = get_node(c)
                    coords_out[nt] = c
                    nlist.append(nt)
                elems.append((etag, nlist))
                etag += 1
    return coords_out, elems


def find_node(coords, target):
    for n, c in coords.items():
        if all(abs(c[i] - target[i]) < 1e-9 for i in range(3)):
            return n
    raise KeyError(target)


# ---------------------------------------------------------------------------
# Case runner
# ---------------------------------------------------------------------------

def run_case(name, coords, fixed_dofs, elem_fn, ele_tags, ngp_per_ele, ndof,
             stiff_name='stiff'):
    print(f"\n=== {name} ===  nDOF={ndof}")
    out = {}
    for lam, label in [(LAM_ELASTIC, 'elastic'), (LAM_PLASTIC, 'plastic')]:
        fixed_set = build_common(coords, elem_fn, fixed_dofs)
        ramp_to(lam)
        branches = []
        for et in ele_tags:
            branches += gp_branches(et, ngp_per_ele)
        n_plastic = sum(1 for b in branches if b == 1)

        K = None
        if len(ele_tags) == 1 and stiff_name is not None:
            K = stiff_via_eleresponse(ele_tags[0], ndof, stiff_name)
        if K is None:
            K = stiff_via_printA(fixed_set, ndof)

        rep = sv_report(K)
        if rep is None:
            print(f"  {label}: plastic GPs {n_plastic}/{len(branches)} -> "
                  f"EXTRACTION REFUSED (NaN/Inf in tangent)")
            out[label] = {'nz': None, 'plastic_gp': n_plastic, 'total_gp': len(branches)}
        else:
            print(f"  {label}: plastic GPs {n_plastic}/{len(branches)} -> "
                  f"near-zero SV count = {rep['nz']}  "
                  f"(smallest 5: {np.array2string(rep['sv'][-5:], precision=3)})")
            out[label] = {'nz': rep['nz'], 'smallest25': rep['smallest25'],
                          'plastic_gp': n_plastic, 'total_gp': len(branches)}
    return out


def main():
    print("ladrunoBuild:", ops.ladrunoBuild())
    results = {}

    def brick8_bbar(matTag):
        ops.element('LadrunoBrick', 1, 1, 2, 3, 4, 5, 6, 7, 8, matTag, '-formulation', 'bbar')
    results['LadrunoBrick -bbar'] = run_case(
        'LadrunoBrick -formulation bbar', BRICK8_COORDS, BRICK8_FIX, brick8_bbar,
        [1], 8, 24)

    for form, nGP in [('uri', 8), ('std', 27)]:
        def elem_fn(matTag, form=form):
            ops.element('LadrunoBrick20', 1, *range(1, 21), matTag, '-formulation', form)
        results[f'LadrunoBrick20 -{form}'] = run_case(
            f'LadrunoBrick20 -formulation {form}', BRICK20_COORDS, BRICK20_FIX, elem_fn,
            [1], nGP, 60)

    def tet10_plain(matTag):
        ops.element('BezierTet10', 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, matTag)
    results['BezierTet10'] = run_case(
        'BezierTet10 (plain)', TET10_COORDS, TET10_FIX, tet10_plain, [1], 4, 30)

    def tet10_bbar(matTag):
        ops.element('BezierTet10', 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, matTag, '-bbar')
    results['BezierTet10 -bbar'] = run_case(
        'BezierTet10 -bbar', TET10_COORDS, TET10_FIX, tet10_bbar, [1], 4, 30)

    def ten_node(matTag):
        ops.element('TenNodeTetrahedron', 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, matTag)
    results['TenNodeTetrahedron'] = run_case(
        'TenNodeTetrahedron (vanilla)', TET10_COORDS, TET10_FIX, ten_node, [1], 4, 30,
        stiff_name=None)

    # Mesh-level H20-uri check: 2x2x2 block of 8 elements, same homogeneous field.
    mcoords, melems = gen_h20_mesh_2x2x2()
    n_origin = find_node(mcoords, (0, 0, 0))
    n_x2 = find_node(mcoords, (2, 0, 0))
    n_y2 = find_node(mcoords, (0, 2, 0))
    mfix = {n_origin: (1, 1, 1), n_x2: (0, 1, 1), n_y2: (0, 0, 1)}

    def mesh_elem_fn(matTag):
        for etag, nlist in melems:
            ops.element('LadrunoBrick20', etag, *nlist, matTag, '-formulation', 'uri')

    results['H20-uri MESH 2x2x2 (8 ele)'] = run_case(
        'LadrunoBrick20 -uri MESH 2x2x2 (8 elements)', mcoords, mfix, mesh_elem_fn,
        [e[0] for e in melems], 8, len(mcoords) * 3, stiff_name=None)

    print("\n\n=== SUMMARY ===")
    for name, r in results.items():
        e = r.get('elastic', {})
        p = r.get('plastic', {})
        print(f"{name:32s} elastic nz={e.get('nz')}   plastic nz={p.get('nz')}")


if __name__ == '__main__':
    main()
