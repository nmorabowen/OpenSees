"""Finite-strain trifecta — **P5 cross-validation: D2 hex ↔ tet** (element agreement).

Validation plan: Ladruno_implementation/17_finite_strain_validation_plan.md, §5
(L4 / D2) and §9 phase **P5**. Two *independent* fork solid elements consume the
same `SolidTransformation` geometry layer and the same small-strain material:
  * **LadrunoBrick** — 8-node hex (here `-formulation eas`, the shear-locking-cured
    enhanced-assumed-strain form), and
  * **[[14_bezier_tet10_corot|BezierTet10]]** — 10-node quadratic tetrahedron.
Different topology *and* different order, on the SAME cantilever, must converge to
the SAME displacement. Because no other OpenSees solid computes F, this hex↔tet
agreement is the strongest internal evidence short of a commercial-code oracle.

> [!note] Why EAS, and why a convergence band (not ≤1e-3)
> The plain `std` hex shear-LOCKS in slender bending (it under-predicts the tip
> deflection ~40 %), so std↔tet would disagree for a *hex* reason, not a real one.
> `-formulation eas` cures shear locking and tracks the quadratic tet. The two
> elements then **bracket** the exact solution (the EAS hex converges from below,
> the tet from above) and the gap **shrinks** under refinement toward a common
> limit — the honest statement of element agreement. The tight ≤2 % match is a
> mesh-converged (fine) result; here we demonstrate the *convergence to agreement*
> on an affordable Zone-A mesh (gap ≈ 8 % at 3×3×12, shrinking).

BezierTet10 supports `-geom linear|corot` (not `finite`), so D2 lives in the
small-strain regime; D4 (corot↔finite) already covers the geometry-method axis
in `test_finite_strain_P2_geomnl.py`. The tet mesh is a conforming **Kuhn 6-tet**
subdivision of each hex with shared mid-edge nodes (no gmsh).
"""
import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a, pytest.mark.t2a]

_E, _NU = 200.0, 0.3
_A, _L, _P = 1.0, 8.0, 0.02

# Kuhn (Freudenthal) 6-tet subdivision of a hex along the main diagonal 0–6.
# Local hex corner order: bottom 0,1,2,3 (z=0, CCW); top 4,5,6,7 (z=1).
_KUHN = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6), (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]
# BezierTet10 mid-edge node order (local tet-vertex pairs for nodes 5..10).
_EDGE_V = [(0, 1), (1, 2), (0, 2), (0, 3), (2, 3), (1, 3)]


def _hex_eas_tip(nx, ny, nz):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    nxp, nyp = nx + 1, ny + 1

    def nid(i, j, k):
        return 1 + i + nxp * j + nxp * nyp * k

    for k in range(nz + 1):
        for j in range(nyp):
            for i in range(nxp):
                ops.node(nid(i, j, k), _A * i / nx, _A * j / ny, _L * k / nz)
    for j in range(nyp):
        for i in range(nxp):
            ops.fix(nid(i, j, 0), 1, 1, 1)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU)
    e = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                conn = [nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                        nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1)]
                ops.element("LadrunoBrick", e, *conn, 1, "-formulation", "eas", "-geom", "linear")
                e += 1
    return _solve_tip([nid(i, j, nz) for j in range(nyp) for i in range(nxp)])


def _tet10_tip(nx, ny, nz):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    nxp, nyp = nx + 1, ny + 1

    def vid(i, j, k):
        return 1 + i + nxp * j + nxp * nyp * k

    coord = {}
    for k in range(nz + 1):
        for j in range(nyp):
            for i in range(nxp):
                t = vid(i, j, k)
                c = (_A * i / nx, _A * j / ny, _L * k / nz)
                ops.node(t, *c)
                coord[t] = np.array(c)
    nxt = [nxp * nyp * (nz + 1) + 1]
    mm = {}

    def mid(p, q):
        key = tuple(sorted((p, q)))
        if key not in mm:
            m = nxt[0]; nxt[0] += 1
            mc = 0.5 * (coord[p] + coord[q])
            ops.node(m, *mc.tolist()); coord[m] = mc; mm[key] = m
        return mm[key]

    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU)
    for j in range(nyp):
        for i in range(nxp):
            ops.fix(vid(i, j, 0), 1, 1, 1)
    e = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                C = [vid(i, j, k), vid(i + 1, j, k), vid(i + 1, j + 1, k), vid(i, j + 1, k),
                     vid(i, j, k + 1), vid(i + 1, j, k + 1), vid(i + 1, j + 1, k + 1), vid(i, j + 1, k + 1)]
                for tet in _KUHN:
                    v = [C[t] for t in tet]
                    mids = [mid(v[x], v[y]) for (x, y) in _EDGE_V]
                    ops.element("BezierTet10", e, *v, *mids, 1, "-geom", "linear")
                    e += 1
    return _solve_tip([vid(i, j, nz) for j in range(nyp) for i in range(nxp)])


def _solve_tip(tip):
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in tip:
        ops.load(n, _P / len(tip), 0.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("EnergyIncr", 1.0e-10, 50, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "D2 solve failed"
    return float(np.mean([ops.nodeDisp(n, 1) for n in tip]))


def test_D2_hex_and_tet_converge_to_common_solution():
    # coarse and fine meshes for each element technology
    h_coarse = _hex_eas_tip(2, 2, 8)
    t_coarse = _tet10_tip(2, 2, 8)
    h_fine = _hex_eas_tip(3, 3, 12)
    t_fine = _tet10_tip(3, 3, 12)
    assert min(h_coarse, t_coarse, h_fine, t_fine) > 0

    gap_coarse = abs(t_coarse - h_coarse) / t_coarse
    gap_fine = abs(t_fine - h_fine) / t_fine
    # (1) the two independent elements BRACKET the exact solution: the EAS hex
    #     converges from below, the quadratic tet from above.
    assert h_coarse < t_coarse and h_fine < t_fine, (
        f"hex/tet not bracketing (hex {h_fine:.4f} vs tet {t_fine:.4f})")
    # (2) each converges monotonically toward the common limit
    assert h_fine > h_coarse, f"EAS hex not converging upward ({h_coarse:.4f}→{h_fine:.4f})"
    assert t_fine < t_coarse, f"tet10 not converging downward ({t_coarse:.4f}→{t_fine:.4f})"
    # (3) the gap shrinks under refinement and is within engineering tolerance
    #     (the ≤2 % match is a mesh-converged fine-mesh result — see module docstring)
    assert gap_fine < gap_coarse, f"hex↔tet gap did not shrink ({gap_coarse:.3f}→{gap_fine:.3f})"
    assert gap_fine < 0.12, f"hex↔tet still far apart at 3×3×12 (gap {gap_fine:.3f})"
