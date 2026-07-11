"""GmshRecorder hex20 round-trip (Zone-A recorder fix).

Upstream GmshRecorder wrote ``Twenty_Node_Brick`` as MSH element type 12 —
which per the MSH 2.2 spec (and the recorder's own header comment table) is
the **27-node** hexahedron — and dumped ``getExternalNodes()`` raw, ignoring
that Gmsh's hex20 mid-edge ordering differs from the OpenSees serendipity
ordering (shp3dv.cpp "Local Node Pattern"). Gmsh could not read the mesh.

The fix maps the element to type **17** (hexahedron_20_node) and permutes the
mid-edge nodes on write:

  OpenSees mid-edge order (0-indexed nodes 8..19), edges:
      (0-1)(1-2)(2-3)(3-0)(4-5)(5-6)(6-7)(7-4)(0-4)(1-5)(2-6)(3-7)
  Gmsh hex20 mid-edge order (ref. manual "Node ordering"), edges:
      (0-1)(0-3)(0-4)(1-2)(1-5)(2-3)(2-6)(3-7)(4-5)(4-7)(5-6)(6-7)
  =>  write order = corners 0..7 then OpenSees nodes
      [8, 11, 16, 9, 17, 10, 18, 19, 12, 15, 13, 14]

This test builds a single unit-cube 20NodeBrick, records one step through
``recorder gmsh``, parses the emitted ``.msh`` and asserts

  * the element is written as type 17, with exactly 20 nodes;
  * matching **by coordinate** (never by tag), every written node sits where
    the Gmsh reference ordering demands: corners 0..7 on the bottom/top rings
    (CCW, bottom first) and mid-edge slots 8..19 on the midpoints of the Gmsh
    edge list above.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# Unit cube, OpenSees/shp3dv serendipity ordering (1-indexed model tags):
#   1-4 bottom ring CCW (z=0), 5-8 top ring CCW (z=1),
#   9-12 mid-edges (1-2)(2-3)(3-4)(4-1), 13-16 mid-edges (5-6)(6-7)(7-8)(8-5),
#   17-20 mid-edges (1-5)(2-6)(3-7)(4-8).
_CORNERS = [
    (0.0, 0.0, 0.0),  # 1
    (1.0, 0.0, 0.0),  # 2
    (1.0, 1.0, 0.0),  # 3
    (0.0, 1.0, 0.0),  # 4
    (0.0, 0.0, 1.0),  # 5
    (1.0, 0.0, 1.0),  # 6
    (1.0, 1.0, 1.0),  # 7
    (0.0, 1.0, 1.0),  # 8
]
_OPS_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0),
              (4, 5), (5, 6), (6, 7), (7, 4),
              (0, 4), (1, 5), (2, 6), (3, 7)]
# Gmsh hexahedron edge order (Gmsh reference manual, "Node ordering" section;
# = MHexahedron edge table in the Gmsh source).
_GMSH_EDGES = [(0, 1), (0, 3), (0, 4), (1, 2), (1, 5), (2, 3),
               (2, 6), (3, 7), (4, 5), (4, 7), (5, 6), (6, 7)]

GMSH_HEX20_TYPE = 17


def _midpoint(a, b):
    return tuple(0.5 * (a[i] + b[i]) for i in range(3))


def _expected_gmsh_coords():
    """Coordinates the .msh must list, position by position (gmsh order)."""
    return list(_CORNERS) + [_midpoint(_CORNERS[a], _CORNERS[b])
                             for a, b in _GMSH_EDGES]


def _model_node_coords():
    """Model tag -> coordinate, OpenSees ordering (tags 1..20)."""
    coords = list(_CORNERS) + [_midpoint(_CORNERS[a], _CORNERS[b])
                               for a, b in _OPS_EDGES]
    return {tag: xyz for tag, xyz in enumerate(coords, start=1)}


def _build_and_record(msh_base):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _model_node_coords().items():
        ops.node(tag, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.25)
    ops.element("20NodeBrick", 1, *range(1, 21), 1)

    for n in (1, 2, 3, 4, 9, 10, 11, 12):  # clamp the bottom face
        ops.fix(n, 1, 1, 1)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(7, 0.0, 0.0, -1.0)

    ops.recorder("gmsh", msh_base, "disp")

    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    ops.wipe()  # flush + close the recorder


def _parse_msh(path):
    """Minimal MSH 2.2 ASCII parser -> (nodes {tag: xyz}, elements list)."""
    lines = path.read_text().splitlines()
    nodes, elements = {}, []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line == "$Nodes":
            n = int(lines[i + 1])
            for j in range(n):
                parts = lines[i + 2 + j].split()
                nodes[int(parts[0])] = tuple(float(v) for v in parts[1:4])
            i += 2 + n
        elif line == "$Elements":
            n = int(lines[i + 1])
            for j in range(n):
                parts = [int(v) for v in lines[i + 2 + j].split()]
                tag, etype, ntags = parts[0], parts[1], parts[2]
                elements.append((tag, etype, parts[3 + ntags:]))
            i += 2 + n
        else:
            i += 1
    return nodes, elements


def test_hex20_written_as_type17_in_gmsh_order(tmp_path):
    base = tmp_path / "hex20"
    _build_and_record(str(base))

    msh = base.parent / (base.name + ".mesh.0.msh")
    assert msh.exists(), "GmshRecorder did not write %s" % msh

    nodes, elements = _parse_msh(msh)
    assert len(nodes) == 20
    assert len(elements) == 1

    _, etype, conn = elements[0]
    assert etype == GMSH_HEX20_TYPE, (
        "20NodeBrick written as MSH type %d, expected 17 (hexahedron_20_node); "
        "type 12 is the 27-node hex" % etype)
    assert len(conn) == 20

    # Coordinate-space check (robust to any node renumbering): the node id
    # written at gmsh position j must sit exactly at the coordinate the gmsh
    # hex20 reference ordering demands.
    expected = _expected_gmsh_coords()
    for j, node_id in enumerate(conn):
        assert node_id in nodes, "connectivity references unknown node %d" % node_id
        got = nodes[node_id]
        want = expected[j]
        dist = math.dist(got, want)
        assert dist < 1e-12, (
            "gmsh position %d: node %d at %s, expected %s (off by %.3e) — "
            "hex20 permutation broken" % (j, node_id, got, want, dist))


def test_hex8_write_path_untouched(tmp_path):
    """Regression guard: the shared writer loop must not permute other types."""
    base = tmp_path / "hex8"
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, xyz in enumerate(_CORNERS, start=1):
        ops.node(tag, *xyz)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.25)
    ops.element("stdBrick", 1, *range(1, 9), 1)
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(7, 0.0, 0.0, -1.0)
    ops.recorder("gmsh", str(base), "disp")
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    ops.wipe()

    msh = base.parent / (base.name + ".mesh.0.msh")
    assert msh.exists()
    nodes, elements = _parse_msh(msh)
    _, etype, conn = elements[0]
    assert etype == 5  # hexahedron_8_node
    assert conn == list(range(1, 9))  # raw getExternalNodes order preserved
