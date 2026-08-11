"""Explicit lane, serial or MP depending on P0_MODE."""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pin
import model as M
import explicit as X

mode = os.environ.get("P0_MODE", "serial")
ops = pin.load("openseesmp" if mode == "mp" else "opensees")

pid = ops.getPID() if mode == "mp" else None
M.build(ops, rank=pid)

# mass on every node this rank declared (ghosts included -> see note in report:
# a ghost with mass would be double-counted by the shared-DOF sum, so ghosts
# are deliberately LEFT MASSLESS here; only natively-owned nodes get mass).
own = list(M.BOT_NODES) if pid == 0 else list(M.TOP_NODES)
if pid is None:
    own = list(M.BOT_NODES) + list(M.TOP_NODES)
if pid == 0 and os.environ.get("P0_GHOSTMASS"):
    own += list(M.INTERFACE_SLAVE)   # mutation: ghosts given mass -> double-count?
X.add_mass(ops, own)

if mode == "mp":
    ok = X.analyse_explicit(ops, "ParallelPlain", "MPIDiagonal")
else:
    ok = X.analyse_explicit(ops, "RCM", "Diagonal")

out = {"mode": mode, "pid": pid, "analyze": ok}
if pid in (None, 1):
    out["w15"] = ops.nodeDisp(15, 3)
    out["w11"] = ops.nodeDisp(11, 3)
if pid in (None, 0):
    out["w5"] = ops.nodeDisp(5, 3)
    if pid == 0:
        out["w11_ghost"] = ops.nodeDisp(11, 3)
print("P0EXP " + json.dumps(out))
