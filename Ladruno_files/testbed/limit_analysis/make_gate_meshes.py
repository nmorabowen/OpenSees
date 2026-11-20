"""Generate the per-angle Drucker-Prager boxes of the true-3D Upper gate (WP-237a, interface section 11).

Mesh data only. No conic program is solved and no bound is implied. The output is committed under
tests/data/limit/ and the contract records its SHA-256, so neither the tests nor any implementation need gmsh.

The box for friction angle phi holds the plane-strain Prandtl extent of a footing of half-width 0.5, taken as a
length scale (interface 16.4): half-width W >= x_extent(phi), depth D >= r1(phi). The footing stays the conforming
square [-0.5, 0.5]^2 on z = 0. Element size grows with distance from the footing edge, so that the four sizes of
each box have about the tetrahedron counts of the phi = 0 meshes.

Run with gmsh 4.15: python make_gate_meshes.py <output-directory>
"""
import sys
from pathlib import Path

import gmsh
import numpy as np

BOXES = {10: (2.5, 2.0), 20: (3.5, 2.0), 30: (5.0, 3.0)}
# name: (absolute size at the footing edge, size far away as a multiple of far_size(W))
LEVELS = {"a": (0.60, 0.92), "b": (0.42, 0.80), "c": (0.27, 0.62), "d": (0.092, 0.70)}
FOOTING_BAND = 0.1


def far_size(half_width):
    return 0.64 * half_width


def make(half_width, depth, size_min, size_max):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("General.NumThreads", 1)
    gmsh.option.setNumber("Mesh.RandomSeed", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.option.setNumber("Mesh.Optimize", 1)
    for option in ("MeshSizeFromPoints", "MeshSizeFromCurvature", "MeshSizeExtendFromBoundary"):
        gmsh.option.setNumber(f"Mesh.{option}", 0)
    gmsh.model.add("gate")
    occ = gmsh.model.occ
    box = occ.addBox(-half_width, -half_width, -depth, 2.0 * half_width, 2.0 * half_width, depth)
    pad = occ.addRectangle(-0.5, -0.5, 0.0, 1.0, 1.0)
    occ.fragment([(3, box)], [(2, pad)])
    occ.synchronize()
    curves = [tag for _, tag in gmsh.model.getEntitiesInBoundingBox(-0.51, -0.51, -0.01, 0.51, 0.51, 0.01, 1)]
    field = gmsh.model.mesh.field
    field.add("Distance", 1)
    field.setNumbers(1, "CurvesList", curves)
    field.setNumber(1, "Sampling", 100)
    field.add("Threshold", 2)
    field.setNumber(2, "InField", 1)
    field.setNumber(2, "SizeMin", size_min)
    field.setNumber(2, "SizeMax", size_max)
    field.setNumber(2, "DistMin", FOOTING_BAND)
    field.setNumber(2, "DistMax", 0.8 * half_width)
    field.setAsBackgroundMesh(2)
    gmsh.model.mesh.generate(3)
    tags, xyz, _ = gmsh.model.mesh.getNodes()
    points = xyz.reshape(-1, 3)
    index = {int(tag): i for i, tag in enumerate(tags)}
    _, connectivity = gmsh.model.mesh.getElementsByType(4)
    tets = np.array([index[int(tag)] for tag in connectivity], dtype=np.int64).reshape(-1, 4)
    gmsh.finalize()
    used = np.unique(tets)
    renumber = -np.ones(len(points), dtype=np.int64)
    renumber[used] = np.arange(len(used))
    points, tets = points[used], renumber[tets]
    planes = ((0, (-half_width, -0.5, 0.5, half_width)), (1, (-half_width, -0.5, 0.5, half_width)), (2, (-depth, 0.0)))
    for axis, values in planes:
        for value in values:
            points[np.abs(points[:, axis] - value) < 1e-10, axis] = value
    for i, tet in enumerate(tets):
        if np.linalg.det(points[tet[1:]] - points[tet[0]]) < 0.0:
            tets[i, [2, 3]] = tets[i, [3, 2]]
    return points, tets


def main(directory):
    directory = Path(directory)
    for phi, (half_width, depth) in BOXES.items():
        scale = far_size(half_width)
        for name, (near, far) in LEVELS.items():
            points, tets = make(half_width, depth, near, far * scale)
            path = directory / f"upper3d_phi{phi:02d}_{name}.mesh.npz"
            np.savez_compressed(path, nodes=points.astype("<f8"), tets=tets.astype("<i8"))
            print(path.name, len(points), len(tets), flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
