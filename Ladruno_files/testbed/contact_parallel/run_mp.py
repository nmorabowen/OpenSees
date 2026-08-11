"""P0 2-rank run. Each rank executes this whole script (openseespymp model).

rank 0 = master-owning rank: bottom block + BOTH contactSurface defs + the
         contact verb + GHOST declarations of the slave interface nodes 11..14.
rank 1 = slave rank: top block + the load. Emits NO contact.

The contact adapter on rank 0 assembles into DOFs of nodes 11..14, which rank 0
owns only as ghosts (no element, no mass). If ADR-78's premise holds, the
distributed SOE sums those contributions with rank 1's native ones.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pin
import model as M

ops = pin.load("openseesmp")

mutate = os.environ.get("P0_MUTATE") or None
pid, np_ = ops.getPID(), ops.getNP()
M.build(ops, rank=pid, mutate=mutate)
ok = M.analyse(ops, "ParallelPlain", "Mumps")

out = {"pid": pid, "np": np_, "mutate": mutate, "analyze": ok}
if pid == 0:
    out["w5"] = ops.nodeDisp(5, 3)
    out["w11_ghost"] = ops.nodeDisp(11, 3)
    ops.reactions()
    out["R_base"] = sum(ops.nodeReaction(t, 3) for t in M.BASE)
else:
    out["w15"] = ops.nodeDisp(15, 3)
    out["w11"] = ops.nodeDisp(11, 3)
print("P0MP " + json.dumps(out))
