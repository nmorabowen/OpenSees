import sys
p = sys.argv[1]
s = open(p).read()

def rep(a, b, n=1):
    global s
    assert s.count(a) == n, (s.count(a), a[:90])
    s = s.replace(a, b)

# 1. allow NorSand with the element switch; number of Gauss points per element
rep('''        if MESH_PERTURB or mat != "dp":
            raise SystemExit("ELEM: dp material and no --mesh-perturb only")''',
'''        if MESH_PERTURB or mat not in ("dp", "norsand"):
            raise SystemExit("ELEM: dp or norsand material and no --mesh-perturb only")''')
rep('''    LST = ELEM in ("lst", "bez", "bezb")          # 6-node triangles on the refined node grid''',
'''    global NGP
    NGP = {"qssp": 1, "cst": 1, "lst": 3, "bez": 3, "bezb": 3}.get(ELEM, 4)   # Gauss points per element
    LST = ELEM in ("lst", "bez", "bezb")          # 6-node triangles on the refined node grid''')
rep('''def ns_make_materials(elems, y, einit):''', '''NGP = 4


def ns_make_materials(elems, y, einit):''')

# 2. triangles: define now, create after the per-row materials exist
rep('''                    ops.element("LadrunoCST", tg, *[N(*v) for v in tri], MAT,
                                "-type", "PlaneStrain", "-thick", 1.0, "-body", 0.0, -GAMMA)''',
'''                    TRI_PENDING.append((tg, ("LadrunoCST", [N(*v) for v in tri], None)))''')
rep('''                    if ELEM == "lst":
                        ops.element("LadrunoLST", tg, *[N(*v) for v in tri], MAT,
                                    "-type", "PlaneStrain", "-thick", 1.0, "-body", 0.0, -GAMMA)
                    else:   # BezierTri6: corners then mid12, mid23, mid31 (control points at edge midpoints)
                        ops.element("BezierTri6", tg, *[N(*v) for v in tri], 1.0, "PlaneStrain", MAT,
                                    *(["-bbar"] if ELEM == "bezb" else []), "-bodyForce", 0.0, -GAMMA)''',
'''                    TRI_PENDING.append((tg, ("LadrunoLST" if ELEM == "lst" else "BezierTri6", [N(*v) for v in tri], None)))''')
rep('''    elems = element_table(x, y) if not TRI else []''', '''    elems = element_table(x, y) if not TRI else []
    TRI_PENDING = []''')
rep('''        ns_tag, ns_rows, matdesc = ns_make_materials(elems, y, einit)''',
'''        ns_tag, ns_rows, matdesc = ns_make_materials(elems, (y[0::2] if LST else y), einit)   # rows = coarse cells
    for tg, (etype, nodes, _) in TRI_PENDING:
        m = ns_tag.get(tg, MAT)
        if etype == "BezierTri6":   # corners then mid12, mid23, mid31 (control points at edge midpoints)
            ops.element("BezierTri6", tg, *nodes, 1.0, "PlaneStrain", m,
                        *(["-bbar"] if ELEM == "bezb" else []), "-bodyForce", 0.0, -GAMMA)
        else:
            ops.element(etype, tg, *nodes, m, "-type", "PlaneStrain", "-thick", 1.0, "-body", 0.0, -GAMMA)''')

# 3. Gauss-point loops of the NorSand path
rep('''        for gp in range(4):
            pi = st[4 * k + gp, 0]''', '''        for gp in range(NGP):
            pi = st[4 * k + gp, 0]''')
rep('''            for gp in (1, 2, 3, 4):
                v = ops.eleResponse(tag, "material", gp, "floor")''', '''            for gp in range(1, NGP + 1):
                v = ops.eleResponse(tag, "material", gp, "floor")''')
rep('''        for gp in (1, 2, 3, 4):
            r = ops.eleResponse(tag, "material", gp, "refusal")''', '''        for gp in range(1, NGP + 1):
            r = ops.eleResponse(tag, "material", gp, "refusal")''')
rep('''        sp = ops.eleResponse(tag, "stressesPlaneStrain")        # [s11 s22 s12 s33] per Gauss point
        for gp in range(4):''', '''        sp = ops.eleResponse(tag, "stressesPlaneStrain")        # [s11 s22 s12 s33] per Gauss point
        if len(sp) != 4 * NGP:
            raise SystemExit(f"stressesPlaneStrain: {len(sp)} values on element {tag}, want {4 * NGP}")
        for gp in range(4):
            if gp >= NGP:          # ELEM with fewer Gauss points: the slot stays NaN
                k += 1
                continue''')
open(p, "w").write(s)
print("patched")
