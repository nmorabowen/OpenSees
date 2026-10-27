"""WP-216: synthetic wall tower for the per-mode RSA output (P1).

Same element types as the tower that motivated P1, built from nothing (no
project data): ASDShellQ4 walls meshed at `h`, ElasticTimoshenkoBeam columns
and beams, one LadrunoMembrane per floor (polygon = the four corner columns,
every other floor-level node a diaphragm node), a rigidDiaphragm per floor on the column and wall-end joints
whose master carries no mass, and the floor mass lumped in UX, UY on every
floor-level node (under Transformation part of it sits on MP slaves, WP-211).

    t = build_tower(nst=10, h=0.4)      # ~13k nodes
    t = build_tower(nst=45, h=0.5)      # ~50k nodes

`t` is a dict of tags: supports, pier groups (element tags + cut nodes per
storey), storey node groups (for inertia resultants), masters.
"""
import math

from ._ops import ops

E_C = 30.0e9
NU_C = 0.2
STOREY_H = 3.0
PLAN = 30.0
FLOOR_MASS = 1.0e3 * PLAN * PLAN       # kg per floor

# wall polylines in plan (x0, y0, x1, y1): a closed core with a cross wall and
# eight perimeter segments inset from the slab edge
CORE = [(11.0, 11.0, 19.0, 11.0), (19.0, 11.0, 19.0, 19.0),
        (19.0, 19.0, 11.0, 19.0), (11.0, 19.0, 11.0, 11.0),
        (15.0, 11.0, 15.0, 19.0)]
PERIM = [(3.0, 2.0, 9.0, 2.0), (21.0, 2.0, 27.0, 2.0),
         (3.0, 28.0, 9.0, 28.0), (21.0, 28.0, 27.0, 28.0),
         (2.0, 4.0, 2.0, 10.0), (2.0, 20.0, 2.0, 26.0),
         (28.0, 4.0, 28.0, 10.0), (28.0, 20.0, 28.0, 26.0)]
CORNERS = [(0.0, 0.0), (PLAN, 0.0), (PLAN, PLAN), (0.0, PLAN)]
INNER_COLS = [(6.0, 6.0), (24.0, 6.0), (24.0, 24.0), (6.0, 24.0),
              (15.0, 6.0), (15.0, 24.0), (6.0, 15.0), (24.0, 15.0)]
BEAMS = [((6.0, 6.0), (15.0, 6.0)), ((15.0, 6.0), (24.0, 6.0)),
         ((6.0, 24.0), (15.0, 24.0)), ((15.0, 24.0), (24.0, 24.0)),
         ((6.0, 6.0), (6.0, 15.0)), ((6.0, 15.0), (6.0, 24.0)),
         ((24.0, 6.0), (24.0, 15.0)), ((24.0, 15.0), (24.0, 24.0)),
         ((0.0, 0.0), (PLAN, 0.0)), ((PLAN, PLAN), (0.0, PLAN))]


def build_tower(nst=10, h=0.4, handler='Transformation', system='UmfPack',
                wall_t=0.3, membrane=True):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    nv = max(1, int(round(STOREY_H / h)))          # wall rows per storey
    dz = STOREY_H / nv
    nodes = {}

    def node(x, y, j):
        key = (round(x, 6), round(y, 6), j)
        tag = nodes.get(key)
        if tag is None:
            tag = len(nodes) + 1
            nodes[key] = tag
            ops.node(tag, x, y, j * dz)
        return tag

    ops.section('ElasticMembranePlateSection', 1, E_C, NU_C, wall_t, 0.0)
    ele = [0]

    def etag():
        ele[0] += 1
        return ele[0]

    piers = {}                          # storey s -> (shell tags, cut nodes)
    core_ele = set()
    for wi, (x0, y0, x1, y1) in enumerate(CORE + PERIM):
        L = math.hypot(x1 - x0, y1 - y0)
        nh = max(1, int(round(L / h)))
        for j in range(nst * nv):
            s = j // nv + 1
            for i in range(nh):
                a, b = i / nh, (i + 1) / nh
                p = (x0 + a * (x1 - x0), y0 + a * (y1 - y0))
                q = (x0 + b * (x1 - x0), y0 + b * (y1 - y0))
                n1, n2 = node(*p, j), node(*q, j)
                n3, n4 = node(*q, j + 1), node(*p, j + 1)
                t = etag()
                ops.element('ASDShellQ4', t, n1, n2, n3, n4, 1)
                if wi < len(CORE):
                    core_ele.add(t)
                    if j % nv == 0:     # bottom row of the storey: a pier cut
                        g = piers.setdefault(s, (set(), set()))
                        g[0].add(t)
                        g[1].update((n1, n2))

    # columns and beams (Timoshenko, rectangular 0.6 x 0.6 and 0.4 x 0.7)
    ops.geomTransf('Linear', 1, 1.0, 0.0, 0.0)      # vertical members
    ops.geomTransf('Linear', 2, 0.0, 0.0, 1.0)      # horizontal members
    G = E_C / (2 * (1 + NU_C))

    def rect(b, d):
        A = b * d
        Iy, Iz = b * d ** 3 / 12, d * b ** 3 / 12
        J = b * d * (b * b + d * d) / 12
        return A, J, Iy, Iz, 5 / 6 * A

    cA, cJ, cIy, cIz, cAv = rect(0.6, 0.6)
    bA, bJ, bIy, bIz, bAv = rect(0.4, 0.7)
    for s in range(1, nst + 1):
        jb, jt = (s - 1) * nv, s * nv
        for (x, y) in CORNERS + INNER_COLS:
            ops.element('ElasticTimoshenkoBeam', etag(), node(x, y, jb), node(x, y, jt),
                        E_C, G, cA, cJ, cIy, cIz, cAv, cAv, 1)
        for (p, q) in BEAMS:
            ops.element('ElasticTimoshenkoBeam', etag(), node(*p, jt), node(*q, jt),
                        E_C, G, bA, bJ, bIy, bIz, bAv, bAv, 2)

    supports = sorted(t for (x, y, j), t in nodes.items() if j == 0)
    for t in supports:
        ops.fix(t, 1, 1, 1, 1, 1, 1)

    masters, storey_nodes, mem_tags = [], {}, []
    nmax = len(nodes)
    for s in range(1, nst + 1):
        j = s * nv
        floor = sorted(t for (x, y, jj), t in nodes.items() if jj == j)
        corner = [nodes[(round(x, 6), round(y, 6), j)] for (x, y) in CORNERS]
        inner = [t for t in floor if t not in corner]
        if membrane:
            mt = etag()
            ops.element('LadrunoMembrane', mt, '-nodes', *corner,
                        '-diaphragmNodes', *inner, '-stiff', E_C, NU_C, 0.2)
            mem_tags.append(mt)
        m = FLOOR_MASS / len(floor)
        for t in floor:
            ops.mass(t, m, m, 0.0, 0.0, 0.0, 0.0)
        master = nmax + s
        ops.node(master, 14.7, 15.3, j * dz)
        ops.fix(master, 0, 0, 1, 1, 1, 0)
        # the diaphragm reaches the "floor object" joints only (columns and
        # wall ends/corners), not the wall-mesh nodes along the floor edge;
        # the membrane ties those in plane
        joints = [nodes[(round(x, 6), round(y, 6), j)] for (x, y) in CORNERS + INNER_COLS]
        for (x0, y0, x1, y1) in CORE + PERIM:
            joints += [nodes[(round(x0, 6), round(y0, 6), j)], nodes[(round(x1, 6), round(y1, 6), j)]]
        ops.rigidDiaphragm(3, master, *sorted(set(joints)))
        masters.append(master)
        # nodes above the storey s base: storey shear at s by inertia
        storey_nodes[s] = sorted(t for (x, y, jj), t in nodes.items() if jj > (s - 1) * nv)

    if handler == 'Penalty':
        ops.constraints('Penalty', 1.0e14, 1.0e14)
    else:
        ops.constraints(handler)
    ops.numberer('RCM')
    ops.system(system)
    ops.test('NormDispIncr', 1e-8, 10)
    ops.algorithm('Linear')
    ops.integrator('LoadControl', 0.0)
    ops.analysis('Static')
    return dict(nst=nst, nv=nv, n_nodes=len(nodes) + nst, n_ele=ele[0],
                supports=supports, piers={s: (sorted(e), sorted(c)) for s, (e, c) in piers.items()},
                core=sorted(core_ele), storey_nodes=storey_nodes, masters=masters,
                membranes=mem_tags)


def spectrum():
    """A smooth design-type spectrum (m/s^2): plateau 0.1-0.6 s, 1/T tail."""
    tn = [0.0] + [0.02 * k for k in range(1, 500)]
    sa = []
    for t in tn:
        if t < 0.1:
            v = 0.4 + 0.6 * t / 0.1
        elif t <= 0.6:
            v = 1.0
        else:
            v = 0.6 / t
        sa.append(9.81 * 0.8 * v)
    return tn, sa
