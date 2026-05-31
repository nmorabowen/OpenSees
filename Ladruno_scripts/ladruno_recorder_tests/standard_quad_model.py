"""Standard-rule QUADRATURE + GLOBAL_GP_COORDS gate -- model runner (Step B).

Exercises the new standard-rule path in writeModelElements: legacy Gauss-Legendre
elements (FourNodeQuad / stdBrick / Tri31) now emit NDIR + NUM_GP +
QUADRATURE/{GP_PARAM,GP_WEIGHT} derived from getStandardQuadrature, PLUS the
belt-and-suspenders GLOBAL_GP_COORDS computed write-side from node->getCrds().

Three independent single-element models (wiped between), each on an AXIS-ALIGNED
UNIT cell so the global GP coords are a trivial affine map of the natural coords
-- giving the checker an analytic oracle that needs no Python basis module:
    quad  : FourNodeQuad over [0,1]^2   (Quad_GL_2, 4 GP CCW)        -> sq_quad.ladruno
    tri   : Tri31 over the unit triangle (Tri_GL_1, 1 GP centroid)   -> sq_tri.ladruno
    brick : stdBrick over [0,1]^3       (Hex_GL_2, 8 GP lexicographic)-> sq_brick.ladruno

Both recorders run so the element channel can be value-diffed later.

Run with the BUILD python (no boot .pth):
    python standard_quad_model.py <dir_with_opensees_pyd> [out_dir]
"""

from __future__ import annotations

import os
import sys

TMP = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else TMP

os.add_dll_directory(TMP)
sys.path.insert(0, TMP)

import opensees as ops  # noqa: E402


def _fresh(*paths):
    for p in paths:
        if os.path.exists(p):
            os.remove(p)


# --------------------------------------------------------------------------- #
# 1) FourNodeQuad over the unit square [0,1]^2  (CCW node order)
# --------------------------------------------------------------------------- #
q_ref = os.path.join(OUT, "sq_quad_ref.mpco")
q_new = os.path.join(OUT, "sq_quad.ladruno")
_fresh(q_ref, q_new)

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)
ops.node(1, 0.0, 0.0)
ops.node(2, 1.0, 0.0)
ops.node(3, 1.0, 1.0)
ops.node(4, 0.0, 1.0)
ops.fix(1, 1, 1)
ops.fix(4, 1, 1)
ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.25)
ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStress", 1)
ops.recorder("mpco", q_ref, "-N", "displacement", "-E", "stress", "-T", "dt", 0.0)
ops.recorder("mpcoLadruno", q_new, "-N", "displacement", "-E", "stress", "-T", "dt", 0.0)
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
ops.load(2, 1.0e2, 0.0)
ops.load(3, 1.0e2, 0.0)
ops.system("BandSPD")
ops.numberer("RCM")
ops.constraints("Plain")
ops.integrator("LoadControl", 0.5)
ops.algorithm("Linear")
ops.analysis("Static")
ops.analyze(2)
ops.wipe()

# --------------------------------------------------------------------------- #
# 2) Tri31 over the unit triangle  (node order N0=xi, N1=eta, N2=1-xi-eta)
# --------------------------------------------------------------------------- #
t_ref = os.path.join(OUT, "sq_tri_ref.mpco")
t_new = os.path.join(OUT, "sq_tri.ladruno")
_fresh(t_ref, t_new)

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)
ops.node(1, 0.0, 0.0)
ops.node(2, 1.0, 0.0)
ops.node(3, 0.0, 1.0)
ops.fix(1, 1, 1)
ops.fix(3, 1, 1)
ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.25)
ops.element("tri31", 1, 1, 2, 3, 1.0, "PlaneStress", 1)
ops.recorder("mpco", t_ref, "-N", "displacement", "-E", "stress", "-T", "dt", 0.0)
ops.recorder("mpcoLadruno", t_new, "-N", "displacement", "-E", "stress", "-T", "dt", 0.0)
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
ops.load(2, 1.0e2, 0.0)
ops.system("BandSPD")
ops.numberer("RCM")
ops.constraints("Plain")
ops.integrator("LoadControl", 0.5)
ops.algorithm("Linear")
ops.analysis("Static")
ops.analyze(2)
ops.wipe()

# --------------------------------------------------------------------------- #
# 3) stdBrick over the unit cube [0,1]^3
#    node order: bottom CCW (---,+--,++-,-+-), top CCW (--+,+-+,+++,-++)
# --------------------------------------------------------------------------- #
b_ref = os.path.join(OUT, "sq_brick_ref.mpco")
b_new = os.path.join(OUT, "sq_brick.ladruno")
_fresh(b_ref, b_new)

ops.wipe()
ops.model("basic", "-ndm", 3, "-ndf", 3)
ops.node(1, 0.0, 0.0, 0.0)  # ---
ops.node(2, 1.0, 0.0, 0.0)  # +--
ops.node(3, 1.0, 1.0, 0.0)  # ++-
ops.node(4, 0.0, 1.0, 0.0)  # -+-
ops.node(5, 0.0, 0.0, 1.0)  # --+
ops.node(6, 1.0, 0.0, 1.0)  # +-+
ops.node(7, 1.0, 1.0, 1.0)  # +++
ops.node(8, 0.0, 1.0, 1.0)  # -++
for n in (1, 2, 3, 4):
    ops.fix(n, 1, 1, 1)
ops.nDMaterial("ElasticIsotropic", 2, 1000.0, 0.25)
ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 2)
ops.recorder("mpco", b_ref, "-N", "displacement", "-E", "stress", "-T", "dt", 0.0)
ops.recorder("mpcoLadruno", b_new, "-N", "displacement", "-E", "stress", "-T", "dt", 0.0)
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
ops.load(5, 0.0, 0.0, 1.0e2)
ops.load(6, 0.0, 0.0, 1.0e2)
ops.load(7, 0.0, 0.0, 1.0e2)
ops.load(8, 0.0, 0.0, 1.0e2)
ops.system("BandSPD")
ops.numberer("RCM")
ops.constraints("Plain")
ops.integrator("LoadControl", 0.5)
ops.algorithm("Linear")
ops.analysis("Static")
ops.analyze(2)
ops.wipe()

# --------------------------------------------------------------------------- #
# 4) NineNodeQuad over [0,1]^2  (higher-order: biquadratic Lagrange, Quad_GL_3)
#    node order: 4 CCW corners, 4 CCW edge-mids (bottom,right,top,left), center
# --------------------------------------------------------------------------- #
q9_ref = os.path.join(OUT, "sq_quad9_ref.mpco")
q9_new = os.path.join(OUT, "sq_quad9.ladruno")
_fresh(q9_ref, q9_new)

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)
ops.node(1, 0.0, 0.0)   # corner --
ops.node(2, 1.0, 0.0)   # corner +-
ops.node(3, 1.0, 1.0)   # corner ++
ops.node(4, 0.0, 1.0)   # corner -+
ops.node(5, 0.5, 0.0)   # edge-mid bottom
ops.node(6, 1.0, 0.5)   # edge-mid right
ops.node(7, 0.5, 1.0)   # edge-mid top
ops.node(8, 0.0, 0.5)   # edge-mid left
ops.node(9, 0.5, 0.5)   # center
ops.fix(1, 1, 1)
ops.fix(4, 1, 1)
ops.fix(8, 1, 1)
ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.25)
ops.element("quad9n", 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1.0, "PlaneStress", 1)
ops.recorder("mpco", q9_ref, "-N", "displacement", "-T", "dt", 0.0)
ops.recorder("mpcoLadruno", q9_new, "-N", "displacement", "-T", "dt", 0.0)
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
ops.load(2, 1.0e2, 0.0)
ops.load(3, 1.0e2, 0.0)
ops.load(6, 1.0e2, 0.0)
ops.system("BandSPD")
ops.numberer("RCM")
ops.constraints("Plain")
ops.integrator("LoadControl", 0.5)
ops.algorithm("Linear")
ops.analysis("Static")
ops.analyze(2)
ops.wipe()

# --------------------------------------------------------------------------- #
# 5) TenNodeTetrahedron (higher-order: quadratic tet, Tet_GL_2 4-pt rule)
#    node order: corners 0:ζ0 1:ζ1 2:ζ2 3:ζ4, edge-mids 4:0-1 5:1-2 6:2-0
#                7:0-3 8:2-3 9:1-3 (with the element's 8<->9 swap)
# --------------------------------------------------------------------------- #
t10_ref = os.path.join(OUT, "sq_tet10_ref.mpco")
t10_new = os.path.join(OUT, "sq_tet10.ladruno")
_fresh(t10_ref, t10_new)

ops.wipe()
ops.model("basic", "-ndm", 3, "-ndf", 3)
ops.node(1, 0.0, 0.0, 0.0)   # corner ζ0
ops.node(2, 1.0, 0.0, 0.0)   # corner ζ1
ops.node(3, 0.0, 1.0, 0.0)   # corner ζ2
ops.node(4, 0.0, 0.0, 1.0)   # corner ζ4
ops.node(5, 0.5, 0.0, 0.0)   # edge 0-1
ops.node(6, 0.5, 0.5, 0.0)   # edge 1-2
ops.node(7, 0.0, 0.5, 0.0)   # edge 2-0
ops.node(8, 0.0, 0.0, 0.5)   # edge 0-3
ops.node(9, 0.0, 0.5, 0.5)   # edge 2-3
ops.node(10, 0.5, 0.0, 0.5)  # edge 1-3
for n in (1, 2, 3, 5, 7):
    ops.fix(n, 1, 1, 1)
ops.nDMaterial("ElasticIsotropic", 2, 1000.0, 0.25)
ops.element("TenNodeTetrahedron", 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 2)
ops.recorder("mpco", t10_ref, "-N", "displacement", "-T", "dt", 0.0)
ops.recorder("mpcoLadruno", t10_new, "-N", "displacement", "-T", "dt", 0.0)
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
ops.load(4, 0.0, 0.0, 1.0e2)
ops.load(8, 0.0, 0.0, 1.0e2)
ops.load(9, 0.0, 0.0, 1.0e2)
ops.load(10, 0.0, 0.0, 1.0e2)
ops.system("BandSPD")
ops.numberer("RCM")
ops.constraints("Plain")
ops.integrator("LoadControl", 0.5)
ops.algorithm("Linear")
ops.analysis("Static")
ops.analyze(2)
ops.wipe()

for p in (q_new, t_new, b_new, q9_new, t10_new):
    print("NEW:", p, os.path.exists(p))
print("STANDARD_QUAD_MODEL OK")
