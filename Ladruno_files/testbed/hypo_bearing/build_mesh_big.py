"""A LARGER version of the benchmark domain, for the collapse study only.

The campaign mesh gives 4.5 B of side clearance and 4 B of depth, which scoping
finding 3 established was enough to stop the box acting as an oedometer for a
PDMY footing at ~31 % of its strength. A COLLAPSE mechanism is a different
question. The cone PDMY actually has is phi_ps = 53.7 deg, and a Prandtl
mechanism's log-spiral radius grows as exp(theta*tan phi) — a factor of 8.5
over a quarter turn at that angle. The general-shear zone is then tens of
metres wide, and 9 m of clearance cannot hold it: the collapse load would be
the box's, not the soil's.

So this rebuilds the same graded all-hex construction at

    x, y in [-30, +30] m   -> 29 m of clearance = 14.5 B   (was 4.5 B)
    z    in [-20,   0] m   -> 10 B deep                     (was 4 B)

keeping the 0.5 m footing-scale resolution and the 4-elements-across-B footing
patch EXACTLY as before, so the two meshes differ only in how far away the
boundary is. 24 x 24 x 14 = 8064 hexes (was 2816).

Run:  C:\\Users\\nmb\\venv\\opensees_env\\Scripts\\python.exe build_mesh_big.py
Writes bearing_mesh_big.npz.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import build_mesh as bm  # noqa: E402

bm.XLIM, bm.ZBOT = 30.0, -20.0
bm.XBLOCKS = [(-30.0, -1.0, 10, -1.0), (-1.0, 1.0, 4, None), (1.0, 30.0, 10, 1.0)]
bm.YBLOCKS = bm.XBLOCKS
bm.ZBLOCKS = [(-20.0, -3.0, 8, -3.0), (-3.0, 0.0, 6, None)]
bm.OUT = os.path.join(HERE, "bearing_mesh_big.npz")

if __name__ == "__main__":
    sys.exit(bm.main())
