import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pin
import model as M

ops = pin.load("opensees")

M.build(ops, rank=None)
ok = M.analyse(ops, "RCM", "UmfPack")

out = {
    "build": os.path.basename(os.path.dirname(ops.__file__)),
    "analyze": ok,
    "w15": ops.nodeDisp(15, 3),
    "w11": ops.nodeDisp(11, 3),
    "w5": ops.nodeDisp(5, 3),
}
ops.reactions()
out["R_base"] = sum(ops.nodeReaction(t, 3) for t in M.BASE)
print("P0SERIAL " + json.dumps(out))
