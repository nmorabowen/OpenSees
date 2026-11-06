"""WP-222: LadrunoMembrane transfer build, fast and bit-identical.

resolveAndBuildTransfer runs eagerly from Domain::addElement and cost ~86 s on
the YOO Guayaquil RevB tower (2795 membranes). Two changes, both required to
leave every output bit unchanged:

(a) a memo of the pure kernel stages (validation mesh, distribution, diaphragm
    mesh + condensation), keyed on the exact bit image of their inputs, so the
    same panel on every storey is built once (LadrunoMembraneMemo.h);
(b) an allocation-free, hypot-screened nearest-feature query and allocation-free
    feature sets in distribute_surface_load (LadrunoMembraneKernel.h).

Checked here: (a) a membrane equals its ``-noMemo`` twin to the bit, repeated
panels hit the memo, and wipe() empties it; (b) the kernel compiled as is and
with -DLADRUNO_MEMBRANE_KERNEL_REFERENCE (the pre-WP-222 query and containers)
prints identical hex floats on panels full of ties.
"""
import os
import shutil
import subprocess
from pathlib import Path

import pytest
from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_KERNEL_DIR = Path(__file__).resolve().parent.parent / "SRC" / "element" / "ladrunoMembrane"

RESPONSES = ["frameLoads", "supportReactions", "tributaryAreas", "area", "subSegmentLoads",
             "subSegments", "overhangTorques", "meshStats", "meshVertices", "meshTriangles",
             "meshSubSegments", "meshVertexFlags", "meshNodalField", "meshTriangleLoadFlow",
             "yieldLines", "membraneForces", "diaphragmVertices", "diaphragmTriangles",
             "meanMembraneForce", "diaphragmForces", "inertialForce"]

# 6 x 4 panel, a 1 x 1 opening, an interior column; (x, y) of every node
OUTER = [(0.0, 0.0), (3.0, 0.0), (6.0, 0.0), (6.0, 4.0), (0.0, 4.0)]
HOLE = [(1.0, 1.5), (2.0, 1.5), (2.0, 2.5), (1.0, 2.5)]
COLUMN = (4.5, 2.0)
OTHER = [(0.0, 0.0), (5.0, 0.0), (5.0, 3.5), (0.0, 3.5)]   # a different panel


def _floor(z, base, outer=OUTER, hole=HOLE, column=COLUMN):
    """nodes + edge frames of one panel at height z; returns (poly, hole, col, frames)"""
    poly = []
    for k, (x, y) in enumerate(outer):
        ops.node(base + k, x, y, z)
        poly.append(base + k)
    hl = []
    for k, (x, y) in enumerate(hole):
        ops.node(base + 20 + k, x, y, z)
        hl.append(base + 20 + k)
    col = []
    if column is not None:
        ops.node(base + 40, *column, z)
        col = [base + 40]
    frames = []
    for k in range(len(poly)):
        tag = base + 60 + k
        ops.element("elasticBeamColumn", tag, poly[k], poly[(k + 1) % len(poly)],
                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1)
        frames.append(tag)
    return poly, hl, col, frames


def _membrane(tag, poly, hl, col, frames, *extra):
    args = ["LadrunoMembrane", tag, "-nodes", *poly, "-frames", *frames, "-method", "tributary",
            "-stiff", 30e6, 0.2, 0.15, "-rho", 2.4, "-meshSize", 0.5]
    if hl:
        args += ["-hole", *hl]
    if col:
        args += ["-supportNodes", *col]
    ops.element(*args, *extra)


def _model():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)


def _outputs(tag):
    return {r: list(ops.eleResponse(tag, r) or []) for r in RESPONSES}


def _set_trial_state():
    for n in ops.getNodeTags():
        for dof in range(1, 7):
            ops.setNodeDisp(n, dof, 1e-3 * ((n * 7 + dof * 13) % 97 - 48) / 48.0, "-trial")
            ops.setNodeAccel(n, dof, ((n * 11 + dof * 5) % 89 - 44) / 44.0, "-trial")


def test_memo_is_bit_identical_to_no_memo_and_hits_repeats():
    _model()
    floors = []
    for f in range(3):                       # the same panel on three storeys
        floors.append(_floor(3.0 * (f + 1), 1000 * (f + 1)))
    other = _floor(12.0, 4000, outer=OTHER, hole=[], column=None)
    for f, fl in enumerate(floors):
        _membrane(100 + f, *fl)
        _membrane(200 + f, *fl, "-noMemo")   # twin on the same nodes, built cold
    _membrane(300, *other)
    _membrane(301, *other, "-noMemo")
    _set_trial_state()

    for a, b in [(100, 200), (101, 201), (102, 202), (300, 301)]:
        oa, ob = _outputs(a), _outputs(b)
        for r in RESPONSES:
            assert oa[r] == ob[r], f"{r}: membrane {a} (memo) != {b} (-noMemo)"
        assert sum(len(v) for v in oa.values()) > 1000

    hits = {t: list(ops.eleResponse(t, "transferMemo")) for t in (100, 101, 102, 200, 300)}
    assert hits[100][1] == 0 and hits[100][2] == 0          # first of its kind: built
    for t in (101, 102):                                     # repeats: reused
        assert hits[t][1] == 1 and hits[t][2] == 1
    assert hits[200][:3] == [0.0, 0.0, 0.0]                  # -noMemo never looks
    assert hits[300][1] == 0 and hits[300][2] == 0          # a different panel misses
    assert hits[300][3] == 6                                 # 2 panels x 3 stages


def test_wipe_empties_the_memo():
    _model()
    _membrane(100, *_floor(3.0, 1000))
    _membrane(101, *_floor(6.0, 2000))
    assert ops.eleResponse(101, "transferMemo")[3] == 3
    _model()                                                 # wipe(): every membrane gone
    _membrane(100, *_floor(3.0, 1000, outer=OTHER, hole=[], column=None))
    m = list(ops.eleResponse(100, "transferMemo"))
    assert m[1] == 0 and m[2] == 0 and m[3] == 3


_AB_DRIVER = r"""
#include "LadrunoMembraneKernel.h"
#include <cstdio>
using namespace LadrunoMembraneKernel;
static void p(double x) { std::printf("%a ", x); }
static void run(const char* method, std::vector<Point2D> outer, std::vector<Segment2D> fr,
                std::vector<std::vector<Point2D>> holes, std::vector<Point2D> sup, double h) {
    auto r = distribute_surface_load(outer, 1.0, method, fr, holes, sup, {1.0, 0.0}, 1.0, 0.2, h);
    std::printf("%s ", method);
    p(r.total_applied_load); p(r.total_applied_mx); p(r.total_applied_my);
    for (auto& f : r.frame_loads) {
        p(f.aL); p(f.bL); p(f.wa_eq); p(f.wb_eq); p(f.F_total); p(f.M_origin); p(f.overhang_torque);
        p(f.end_F_a); p(f.end_F_b);
        for (auto& s : f.sub_loads) { p(s.aL); p(s.bL); p(s.wa); p(s.wb); p(s.F); p(s.M_local); }
    }
    for (double x : r.support_node_reactions) p(x);
    for (auto& v : r.mesh.vertices) { p(v.x); p(v.y); }
    for (int x : r.vertex_flags) std::printf("%d ", x);
    for (double x : r.nodal_field) p(x);
    for (int x : r.triangle_support_id) std::printf("%d ", x);
    for (auto& q : r.triangle_flux_2d) { p(q.x); p(q.y); }
    for (auto& y : r.yield_lines) { p(y.x1); p(y.y1); p(y.x2); p(y.y2); std::printf("%d %d ", y.fidA, y.fidB); }
    std::printf("\n");
}
int main() {
    std::vector<Point2D> sq = {{0, 0}, {6, 0}, {6, 6}, {0, 6}};
    std::vector<Segment2D> edges = {{{0, 0}, {6, 0}}, {{6, 0}, {6, 6}}, {{6, 6}, {0, 6}}, {{0, 6}, {0, 0}}};
    std::vector<Point2D> ell = {{0, 0}, {8, 0}, {8, 3}, {3, 3}, {3, 7}, {0, 7}};
    std::vector<Segment2D> ellFr = {{{0, 0}, {8, 0}}, {{8, 0}, {8, 3}}, {{0, 7}, {0, 0}}, {{3, 3}, {3, 7}}};
    std::vector<std::vector<Point2D>> hole = {{{1, 1}, {2, 1}, {2, 2}, {1, 2}}};
    std::vector<Point2D> cols = {{3, 3}, {4.5, 1.5}};
    for (const char* m : {"tributary", "oneway", "fea_poisson", "fea_plate"}) {
        run(m, sq, edges, {}, {}, 0.5);               // symmetric: ties on both diagonals
        run(m, sq, edges, {}, {{3, 3}}, 0.75);        // + a column at the centre
        run(m, ell, ellFr, hole, cols, 0.5);          // re-entrant corner, hole, two columns
    }
    run("tributary", sq, {}, {}, sq, 0.5);            // corner supports only
    return 0;
}
"""


def _ab_run(tmp_path, gpp, name, defines):
    src = tmp_path / "ab.cpp"
    src.write_text(_AB_DRIVER, encoding="utf-8")
    exe = tmp_path / (name + (".exe" if os.name == "nt" else ""))
    subprocess.run([gpp, "-O2", "-std=c++17", *defines, f"-I{_KERNEL_DIR}", str(src), "-o", str(exe)],
                   check=True, capture_output=True, text=True, stdin=subprocess.DEVNULL)
    return subprocess.run([str(exe)], check=True, capture_output=True, text=True,
                          stdin=subprocess.DEVNULL).stdout


def test_kernel_matches_reference_query_bit_for_bit(tmp_path):
    gpp = shutil.which("g++")
    if gpp is None:
        pytest.skip("g++ not on PATH")
    fast = _ab_run(tmp_path, gpp, "fast", [])
    ref = _ab_run(tmp_path, gpp, "ref", ["-DLADRUNO_MEMBRANE_KERNEL_REFERENCE"])
    assert fast.count("\n") == 13 and len(fast) > 10000
    assert fast == ref
