"""Make footing_ab_ns_k5.py from footing_ab_ns.py (K5 second orientation, 2026-11-11).
One change: --mesh-perturb shear:DEG is accepted for norsand. The shear field of mesh_perturb.py moves nodes in x only,
so every element row keeps its depth and the per-row initial stress keeps its meaning. jitter stays refused (it moves y).
Run in ~/wp144/deck."""
src = open("footing_ab_ns.py", encoding="utf-8").read()
a = ('        if MESH_PERTURB:\n'
     '            raise SystemExit("norsand deck: --mesh-perturb is not supported (per-row initial stress)")\n')
assert src.count(a) == 1
b = ('        if MESH_PERTURB and not MESH_PERTURB.startswith("shear"):\n'
     '            raise SystemExit("norsand deck: only --mesh-perturb shear:DEG is supported (rows keep their depth)")\n')
open("footing_ab_ns_k5.py", "w", encoding="utf-8").write(src.replace(a, b))
print("patched")
