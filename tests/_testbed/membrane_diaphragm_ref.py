"""WP-204 reference harness: LadrunoMembrane in-plane response vs a refined
ASDShellQ4 membrane mesh of the same panel.

A rectangular panel (0..LX, 0..LY) with an optional rectangular hole; the
retained points (outer corners, hole corners, any extra point) are nodes of
both models. The reference is a structured ASDShellQ4 grid of spacing h whose
lines pass through every retained point; out-of-plane DOFs are fixed, so only
the membrane works. Loads and in-plane fixities act on retained points only, so
static condensation is exact and the two models must agree to discretisation
error.

WP-208: the slab edge between two consecutive element nodes on a loop is
straight (as in the linear-edge elements it connects to), so with
``boundary_step`` every outer and hole edge also carries a node every
``boundary_step`` — the membrane's loops then have exactly the boundary nodes
of a reference grid of that spacing, as ETABS slab and wall meshes share theirs.
"""
import numpy as np

from _testbed import ops


class Panel:
    def __init__(self, lx, ly, hole=None, extra=(), E=25.0e6, nu=0.2, t=0.2, boundary_step=None):
        self.lx, self.ly, self.hole = lx, ly, hole
        self.E, self.nu, self.t = E, nu, t
        self.outer = [(0.0, 0.0), (lx, 0.0), (lx, ly), (0.0, ly)]
        self.holep = ([] if hole is None else
                      [(hole[0], hole[1]), (hole[2], hole[1]), (hole[2], hole[3]), (hole[0], hole[3])])
        self.extra = list(extra)
        self.points = self.outer + self.holep + self.extra   # retained, ids 1..n
        self.outer_loop = [0, 1, 2, 3]
        self.hole_loop = [4, 5, 6, 7] if hole is not None else []
        if boundary_step:
            self.outer_loop = self._densify(self.outer, [0, 1, 2, 3], boundary_step)
            if hole is not None:
                self.hole_loop = self._densify(self.holep, [4, 5, 6, 7], boundary_step)

    def _index_of(self, x, y):
        for k, (px, py) in enumerate(self.points):
            if abs(px - x) < 1e-9 and abs(py - y) < 1e-9:
                return k
        self.points.append((x, y))
        return len(self.points) - 1

    def _densify(self, corners, corner_idx, step):
        loop = []
        for c in range(len(corners)):
            (x0, y0), (x1, y1) = corners[c], corners[(c + 1) % len(corners)]
            n = int(round(np.hypot(x1 - x0, y1 - y0) / step))
            loop.append(corner_idx[c])
            for k in range(1, n):
                loop.append(self._index_of(x0 + (x1 - x0) * k / n, y0 + (y1 - y0) * k / n))
        return loop

    # --------------------------------------------------------------- solve
    @staticmethod
    def _solve(tags, loads, fix):
        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1, 1)
        for i, (fx, fy) in loads.items():
            ops.load(tags[i], fx, fy, 0.0, 0.0, 0.0, 0.0)
        ops.system('UmfPack'); ops.numberer('RCM'); ops.constraints('Plain')
        ops.integrator('LoadControl', 1.0); ops.algorithm('Linear'); ops.analysis('Static')
        assert ops.analyze(1) == 0
        ops.reactions()
        u = np.array([[ops.nodeDisp(t, 1), ops.nodeDisp(t, 2)] for t in tags])
        r = np.array([[ops.nodeReaction(t, 1), ops.nodeReaction(t, 2)] for t in tags])
        return u, r

    def membrane(self, loads, fix, extra_args=(), diaphragm=None):
        """fix: {index: (fx, fy)}; extra points go in as -supportNodes unless
        `diaphragm` lists the indices to pass as -diaphragmNodes."""
        ops.wipe(); ops.model('basic', '-ndm', 3, '-ndf', 6)
        tags = list(range(1, len(self.points) + 1))
        for k, (x, y) in enumerate(self.points):
            fx, fy = fix.get(k, (0, 0))
            ops.node(tags[k], x, y, 0.0)
            ops.fix(tags[k], fx, fy, 1, 1, 1, 1)
        args = ['-nodes', *[tags[i] for i in self.outer_loop]]
        if self.hole is not None:
            args += ['-hole', *[tags[i] for i in self.hole_loop]]
        on_loops = set(self.outer_loop) | set(self.hole_loop)
        ext = [tags[i] for i in range(len(self.points)) if i not in on_loops]
        dia = [] if diaphragm is None else [tags[i] for i in diaphragm if i not in on_loops]
        sup = [t for t in ext if t not in dia]
        if sup:
            args += ['-supportNodes', *sup]
        if dia:
            args += ['-diaphragmNodes', *dia]
        ops.element('LadrunoMembrane', 100, *args, '-stiff', self.E, self.nu, self.t, *extra_args)
        return self._solve(tags, loads, fix)

    def reference(self, loads, fix, h=0.05, section=None):
        ops.wipe(); ops.model('basic', '-ndm', 3, '-ndf', 6)
        nx, ny = int(round(self.lx / h)), int(round(self.ly / h))
        tag = lambda i, j: 1 + i + j * (nx + 1)
        hole = self.hole
        inhole = (lambda x, y: hole is not None and hole[0] < x < hole[2] and hole[1] < y < hole[3])
        if section is None:
            ops.section('ElasticMembranePlateSection', 1, self.E, self.nu, self.t, 0.0)
        else:
            section()
        ptag = []
        for x, y in self.points:
            i, j = int(round(x / h)), int(round(y / h))
            assert abs(i * h - x) < 1e-9 and abs(j * h - y) < 1e-9, "point off the grid"
            ptag.append(tag(i, j))
        fixmap = {ptag[k]: v for k, v in fix.items()}
        used = set()
        e = 1
        for j in range(ny):
            for i in range(nx):
                if inhole((i + 0.5) * h, (j + 0.5) * h):
                    continue
                ns = [tag(i, j), tag(i + 1, j), tag(i + 1, j + 1), tag(i, j + 1)]
                for n in ns:
                    if n not in used:
                        used.add(n)
                        ii, jj = (n - 1) % (nx + 1), (n - 1) // (nx + 1)
                        ops.node(n, ii * h, jj * h, 0.0)
                        fx, fy = fixmap.get(n, (0, 0))
                        ops.fix(n, fx, fy, 1, 1, 1, 0)
                ops.element('ASDShellQ4', e, *ns, 1)
                e += 1
        return self._solve(ptag, loads, fix)


def rel_err(a, b):
    return float(np.abs(np.asarray(a) - np.asarray(b)).max() / max(np.abs(np.asarray(b)).max(), 1e-300))
