"""WP-114 (TIMs F17) asks 2 + 3: rigid rough punch on a plane-strain
Drucker-Prager block, psi=0 vs psi=phi, BezierTri6 (-bbar / std) vs the
LadrunoQuad -bbar band.

Half model (symmetry at x=0), punch half-width b=1 m (B=2 m), block 8 m x 6 m,
graded toward the punch corner (x=1, y=0) and the surface. Stage 1 (t in [0,1])
ramps self-weight (gamma=18 kN/m3) and a surcharge q=20 kPa on the surface
outside the punch, punch nodes free in uy. Stage 2 places the rigid punch on the
settled surface: every CP on y=0, x<=b (ux fixed = rough) gets the SAME uy
increment from its own settled value (per-node Path series, displacement
control through LoadControl on pseudo time, 1 m per unit time).

Material: UW nDMaterial DruckerPrager, phi=33 deg, c=2 kPa, matched to plane-
strain Mohr-Coulomb (associated form): alpha=tan(phi)/sqrt(9+12tan^2 phi),
k=3c/sqrt(9+12tan^2 phi); UW DP f=||s|| + rho*I1 - sqrt(2/3)*sigma_y, so
rho=sqrt(2)*alpha and sigma_y=sqrt(3)*k. rho_bar = rho (psi=phi) or 0 (psi=0).
Kinf=Ko=delta1=delta2=H=0: no hardening, tension cutoff at the cone apex
I1 = To = sqrt(2/3)*sigma_y/rho (so "cutoff" and "apex" are the same point).

Stage 1 is ELASTIC by default (UW DP `materialState` 0 via setParameter, the
usual geostatic practice and the reporter's "elastic K0 stage"); the switch to
materialState 2 happens before the punch. `--stage1 plastic` keeps DP active in
the gravity stage (the pre-fix -bbar fails there already, see the scope doc).

Census: the ADR-95 `ladrunoBranch` response per GP (0 elastic, 1 cone,
2 tension cutoff, 3 corner = apex), taken after EVERY converged step so the
last converged step's census is exact (the branch after a failed step is a
trial). Location: distance of each census GP from the punch corner.

Surcharge nodal loads on the tri6 surface edges: `bernstein` (consistent for
the Bezier basis: qL/3 per CP) or `lagrange` (qL/6, 2qL/3, qL/6 -- the WRONG
loads for Bezier CPs, used to test the "Bernstein loading artefact" lead).

Usage (one case per process, binary chosen by --bin):
    python3.12 punch_nonassoc.py --bin <dir> --elem tri6_bbar --grid 1 --psi 0
"""

import argparse
import json
import math
import os
import sys
import time

B_HALF = 1.0
X_MAX, Y_MIN = 8.0, -6.0
GAMMA = 18.0
Q_SUR = float(os.environ.get('WP114_Q', 20.0))
K_BULK, G_SHEAR = 27777.78, 9259.26     # E = 25 MPa, nu = 0.35
PHI = math.radians(33.0)
COH = float(os.environ.get('WP114_C', 2.0))
_den = math.sqrt(9.0 + 12.0 * math.tan(PHI) ** 2)
RHO = math.sqrt(2.0) * math.tan(PHI) / _den
SIGY = math.sqrt(3.0) * 3.0 * COH / _den
TO = math.sqrt(2.0 / 3.0) * SIGY / RHO

# solution strategy, identical for every element: Newton, then two fallbacks,
# then halve the step (up to max_halvings) -> "walled"
TEST = ('NormDispIncr', float(os.environ.get('WP114_TOL', 1.0e-8)), 50, 0)
ALGS = [('Newton',), ('KrylovNewton',), ('NewtonLineSearch', 0.8)]

# grid levels: (punch intervals, outer x intervals, y intervals)
GRIDS = {1: (4, 10, 10), 2: (8, 20, 20), 3: (16, 40, 40)}


def _graded(x0, x1, n, h0):
    """n intervals from x0 to x1, first interval ~h0, geometric growth."""
    L = x1 - x0
    if abs(n * h0 - L) < 1e-12:
        return [x0 + L * i / n for i in range(n + 1)]
    lo, hi = 1.0 + 1e-9, 3.0
    for _ in range(200):
        r = 0.5 * (lo + hi)
        s = h0 * (r ** n - 1.0) / (r - 1.0)
        if s > L:
            hi = r
        else:
            lo = r
    xs, x, h = [x0], x0, h0
    for _ in range(n):
        x += h
        xs.append(x)
        h *= r
    xs[-1] = x1
    return xs


def grid(level):
    npunch, nout, ny = GRIDS[level]
    h = B_HALF / npunch
    xs = [B_HALF * i / npunch for i in range(npunch + 1)]
    xs += _graded(B_HALF, X_MAX, nout, h)[1:]
    ys = [-v for v in _graded(0.0, -Y_MIN, ny, h)]
    return xs, ys


def build(ops, elem, level, psi, load):
    xs, ys = grid(level)
    nx, ny = len(xs) - 1, len(ys) - 1
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    rho_bar = RHO if psi == 'phi' else 0.0
    ops.nDMaterial('DruckerPrager', 1, K_BULK, G_SHEAR, SIGY, RHO, rho_bar,
                   0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0)
    vid = {}
    coords = {}
    tag = 0
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            tag += 1
            vid[i, j] = tag
            coords[tag] = (x, y)
            ops.node(tag, x, y)
    mids = {}

    def mid(a, b):
        key = (min(a, b), max(a, b))
        if key not in mids:
            nonlocal tag
            tag += 1
            (xa, ya), (xb, yb) = coords[a], coords[b]
            coords[tag] = (0.5 * (xa + xb), 0.5 * (ya + yb))
            ops.node(tag, *coords[tag])
            mids[key] = tag
        return mids[key]

    eles = {}      # tag -> list of vertex/CP node tags
    etag = 0
    for j in range(ny):          # j: row from the surface down (ys decreasing)
        for i in range(nx):
            # CCW quad: bottom-left, bottom-right, top-right, top-left
            n1, n2, n3, n4 = vid[i, j + 1], vid[i + 1, j + 1], vid[i + 1, j], vid[i, j]
            if elem == 'quad':
                etag += 1
                ops.element('LadrunoQuad', etag, n1, n2, n3, n4, 1,
                            '-formulation', 'bbar', '-type', 'PlaneStrain', '-thick', 1.0)
                eles[etag] = [n1, n2, n3, n4]
                continue
            # union-jack split: alternate the diagonal
            if (i + j) % 2 == 0:
                tris = [(n1, n2, n3), (n1, n3, n4)]
            else:
                tris = [(n1, n2, n4), (n2, n3, n4)]
            for a, b, c in tris:
                etag += 1
                conn = [a, b, c, mid(a, b), mid(b, c), mid(c, a)]
                extra = ['-bbar'] if elem == 'tri6_bbar' else []
                ops.element('BezierTri6', etag, *conn, 1.0, 'PlaneStrain', 1, *extra)
                eles[etag] = conn

    # boundary conditions
    tol = 1e-9
    punch = []
    for t, (x, y) in coords.items():
        if abs(y - Y_MIN) < tol:
            ops.fix(t, 1, 1)
        elif abs(x) < tol or abs(x - X_MAX) < tol:
            ops.fix(t, 1, 0)
        if abs(y) < tol and x <= B_HALF + tol:
            punch.append(t)
            if abs(x) > tol:          # x=0 already has ux fixed
                ops.fix(t, 1, 0)      # rough punch
    punch.sort()

    # stage-1 ramp series and stage-2 punch series (pseudo time)
    ops.timeSeries('Path', 1, '-time', 0.0, 1.0, 1.0e6, '-values', 0.0, 1.0, 1.0)
    ops.pattern('Plain', 1, 1)
    ops.eleLoad('-ele', *eles.keys(), '-type', '-selfWeight', 0.0, -GAMMA)
    # surcharge on the surface edges outside the punch
    surf = [vid[i, 0] for i in range(nx + 1)]
    fx = {}
    for a, b in zip(surf[:-1], surf[1:]):
        xa, xb = coords[a][0], coords[b][0]
        if xa < B_HALF - tol:
            continue
        L = xb - xa
        if elem == 'quad':
            w = {a: 0.5, b: 0.5}
        else:
            m = mids[(min(a, b), max(a, b))]
            w = ({a: 1 / 3, m: 1 / 3, b: 1 / 3} if load == 'bernstein'
                 else {a: 1 / 6, m: 2 / 3, b: 1 / 6})
        for n, f in w.items():
            fx[n] = fx.get(n, 0.0) + f * L * Q_SUR
    for n, f in fx.items():
        if n in punch:
            continue
        ops.load(n, 0.0, -f)
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('UmfPack')
    ops.test('NormDispIncr', 1.0e-9, 40, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 0.05)
    ops.analysis('Static')
    return eles, coords, punch


def gp_coords(ops, elem, eles, coords):
    out = {}
    g = 1.0 / math.sqrt(3.0)
    pts = [(-g, -g), (g, -g), (g, g), (-g, g)]
    for e, conn in eles.items():
        if elem == 'quad':
            xy = [coords[n] for n in conn]
            lst = []
            for xi, eta in pts:
                N = [0.25 * (1 - xi) * (1 - eta), 0.25 * (1 + xi) * (1 - eta),
                     0.25 * (1 + xi) * (1 + eta), 0.25 * (1 - xi) * (1 + eta)]
                lst.append((sum(Ni * p[0] for Ni, p in zip(N, xy)),
                            sum(Ni * p[1] for Ni, p in zip(N, xy))))
            out[e] = lst
        else:
            v = ops.eleResponse(e, 'gaussPoint')
            out[e] = [(v[2 * k], v[2 * k + 1]) for k in range(3)]
    return out


def census(ops, elem, eles):
    ngp = 4 if elem == 'quad' else 3
    rows = []
    for e in eles:
        for g in range(1, ngp + 1):
            br = ops.eleResponse(e, 'material', g, 'ladrunoBranch')
            rows.append((e, g, int(br[0]), br[6]))
    return rows


def run(ops, elem, level, psi, load='bernstein', smax=0.10, dt0=0.4e-2 / 2.0,
        max_halvings=6, verbose=False, stage1='elastic'):
    """dt is in metres of punch displacement (series slope 1 m / unit time)."""
    t0 = time.time()
    eles, coords, punch = build(ops, elem, level, psi, load)
    # stage 1
    ele_tags = list(eles.keys())
    if stage1 == 'elastic':
        ops.setParameter('-val', 0, '-ele', *ele_tags, 'materialState')
    n_stage1 = 0
    while ops.getTime() < 1.0 - 1e-12:
        ok = ops.analyze(1)
        for alg in ALGS[1:]:
            if ok == 0:
                break
            ops.algorithm(*alg)
            ok = ops.analyze(1)
        ops.algorithm(*ALGS[0])
        n_stage1 += 1
        if ok != 0:
            st1 = census(ops, elem, eles)      # TRIAL branches of the failed step
            gpc = gp_coords(ops, elem, eles, coords)
            where = [(round(gpc[e][g - 1][0], 3), round(gpc[e][g - 1][1], 3), br)
                     for e, g, br, _ in st1 if br in (2, 3)]
            return {'error': 'stage-1 failed', 'build': ops.ladrunoBuild(),
                    't_last_converged': ops.getTime(), 'elem': elem, 'psi': psi,
                    'grid': level, 'n_cone_trial': sum(1 for r in st1 if r[2] == 1),
                    'cutoff_apex_trial_xy_branch': where, 'wall_s': time.time() - t0}
    if stage1 == 'elastic':
        ops.setParameter('-val', 2, '-ele', *ele_tags, 'materialState')
    # stage 2: punch placed on the settled surface; every punch CP then moves
    # by the same increment (rigid punch), starting from its own settled uy.
    t1 = ops.getTime()
    for k, n in enumerate(punch):
        u1 = ops.nodeDisp(n, 2)
        ops.timeSeries('Path', 100 + k, '-time', t1, t1 + 1.0e6,
                       '-values', u1, u1 - 1.0e6)
        ops.pattern('Plain', 100 + k, 100 + k)
        ops.sp(n, 2, 1.0)
    # fresh SOE/analysis after the domain change (UmfPack setSize on the old
    # SOE returns symbolic -8 here)
    ops.wipeAnalysis()
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('UmfPack')
    ops.test(*TEST)
    ops.algorithm(*ALGS[0])
    ops.integrator('LoadControl', dt0)
    ops.analysis('Static')
    dmax = smax * 2.0 * B_HALF
    dt, dt_min = dt0, dt0 / 2 ** max_halvings
    curve = [(0.0, 0.0)]
    n_steps = n_cuts = n_alg = 0
    last_census = census(ops, elem, eles)
    hist = []
    walled = False
    while True:
        d = ops.getTime() - t1
        if d >= dmax - 1e-12:
            break
        h = min(dt, dmax - d)
        ops.integrator('LoadControl', h)
        ok = ops.analyze(1)
        for alg in ALGS[1:]:
            if ok == 0:
                break
            n_alg += 1
            ops.algorithm(*alg)
            ok = ops.analyze(1)
        ops.algorithm(*ALGS[0])
        if ok != 0:
            ops.setTime(t1 + d)
            n_cuts += 1
            dt *= 0.5
            if dt < dt_min * 0.999:
                walled = True
                break
            continue
        n_steps += 1
        ops.reactions()
        F = -sum(ops.nodeReaction(n, 2) for n in punch)
        d = ops.getTime() - t1
        curve.append((d / (2 * B_HALF), F / B_HALF))
        last_census = census(ops, elem, eles)
        nc = sum(1 for r in last_census if r[2] in (2, 3))
        hist.append((d / (2 * B_HALF), nc))
        if verbose:
            print(f'  s/B={d / (2 * B_HALF):.4f} q={F / B_HALF:8.2f} dt={h:.2e} cutoff/apex={nc}',
                  flush=True)
        dt = min(dt * 2.0, dt0)
    gpc = gp_coords(ops, elem, eles, coords)
    cen = []
    for e, g, br, i1 in last_census:
        if br in (2, 3):
            x, y = gpc[e][g - 1]
            cen.append({'e': e, 'g': g, 'branch': br, 'I1': i1, 'x': x, 'y': y,
                        'r_corner': math.hypot(x - B_HALF, y)})
    qs = [q for _, q in curve]
    peak = max(qs)
    i_pk = qs.index(peak)
    return {
        'build': ops.ladrunoBuild(), 'elem': elem, 'grid': level, 'psi': psi, 'load': load,
        'n_ele': len(eles), 'n_nodes': len(coords), 'walled': walled,
        's_over_B_last': curve[-1][0], 'q_last': curve[-1][1],
        'q_peak': peak, 's_over_B_peak': curve[i_pk][0],
        'n_steps': n_steps, 'n_cuts': n_cuts, 'n_alg_fallbacks': n_alg,
        'n_cone': sum(1 for r in last_census if r[2] == 1),
        'n_cutoff': sum(1 for r in last_census if r[2] == 2),
        'n_apex': sum(1 for r in last_census if r[2] == 3),
        'n_gp': len(last_census), 'census': cen, 'gp_all': last_census, 'census_hist': hist,
        'curve': curve, 'wall_s': time.time() - t0,
        'stage1': stage1, 'tol': TEST[1], 'c': COH, 'q_sur': Q_SUR, 'rho': RHO, 'sigy': SIGY, 'To': TO,
    }


def _load_ops(bindir):
    os.environ.setdefault('LADRUNO_OPENSEES_QUIET', '1')
    sys.path.insert(0, bindir)
    import opensees as ops
    print('ladrunoBuild:', ops.ladrunoBuild(), flush=True)
    return ops


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--bin', required=True)
    ap.add_argument('--elem', required=True, choices=['tri6_bbar', 'tri6_std', 'quad'])
    ap.add_argument('--grid', type=int, default=1)
    ap.add_argument('--psi', default='0', choices=['0', 'phi'])
    ap.add_argument('--load', default='bernstein', choices=['bernstein', 'lagrange'])
    ap.add_argument('--smax', type=float, default=0.10)
    ap.add_argument('--out', default=None)
    ap.add_argument('--stage1', default='elastic', choices=['elastic', 'plastic'])
    ap.add_argument('-v', action='store_true')
    a = ap.parse_args()
    ops = _load_ops(a.bin)
    r = run(ops, a.elem, a.grid, a.psi, a.load, a.smax, verbose=a.v, stage1=a.stage1)
    brief = {k: v for k, v in r.items() if k not in ('census', 'curve', 'census_hist', 'gp_all')}
    if 'census' in r and r['census']:
        rc = sorted(c['r_corner'] for c in r['census'])
        brief['census_r_corner_min_med_max'] = (rc[0], rc[len(rc) // 2], rc[-1])
        brief['census_y_min'] = min(c['y'] for c in r['census'])
        brief['census_x_range'] = (min(c['x'] for c in r['census']), max(c['x'] for c in r['census']))
    print(json.dumps(brief, indent=1))
    if a.out:
        with open(a.out, 'w') as f:
            json.dump(r, f)
