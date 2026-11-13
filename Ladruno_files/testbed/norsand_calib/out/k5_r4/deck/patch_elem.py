import sys
p = sys.argv[1]
s = open(p).read()

def rep(a, b):
    global s
    assert s.count(a) == 1, (s.count(a), a[:80])
    s = s.replace(a, b)

rep('''    LST = os.environ.get("ELEM_LST") == "1"     # 2026-11-13 (O15 control): 6-node triangles (LadrunoLST std)
    if LST:
        if MESH_PERTURB or mat != "dp":
            raise SystemExit("ELEM_LST: dp material and no --mesh-perturb only")
''', '''    # 2026-11-13 (O15 element control): ELEM selects the element technology (dp material only for non-default):
    #   qbbar (default, LadrunoQuad bbar) | qstd | qssp | qeas | cst (LadrunoCST, 2 per cell)
    #   lst (LadrunoLST std, 2 per cell) | bez | bezb (BezierTri6, std / -bbar, 2 per cell)
    ELEM = os.environ.get("ELEM", "lst" if os.environ.get("ELEM_LST") == "1" else "qbbar")
    if ELEM not in ("qbbar", "qstd", "qssp", "qeas", "cst", "lst", "bez", "bezb"):
        raise SystemExit(f"unknown ELEM {ELEM}")
    LST = ELEM in ("lst", "bez", "bezb")          # 6-node triangles on the refined node grid
    TRI = ELEM in ("cst", "lst", "bez", "bezb")   # triangles: no quad loop
    if ELEM != "qbbar":
        if MESH_PERTURB or mat != "dp":
            raise SystemExit("ELEM: dp material and no --mesh-perturb only")
        log(f"ELEM: {ELEM}")
    if LST:
''')
rep('''    elems = element_table(x, y) if not LST else []
    if LST:''', '''    elems = element_table(x, y) if not TRI else []
    if ELEM == "cst":
        tg = 0
        for ic in range(nx - 1):
            for jc in range(ny - 1):
                N = lambda di, dj: nt(ic + di, jc + dj)
                if (ic + jc) % 2 == 0:
                    tris = [((0, 0), (1, 0), (1, 1)), ((0, 0), (1, 1), (0, 1))]
                else:
                    tris = [((0, 0), (1, 0), (0, 1)), ((1, 0), (1, 1), (0, 1))]
                for tri in tris:
                    tg += 1
                    ops.element("LadrunoCST", tg, *[N(*v) for v in tri], MAT,
                                "-type", "PlaneStrain", "-thick", 1.0, "-body", 0.0, -GAMMA)
                    cx = [x[ic + v[0]] for v in tri]; cy = [y[jc + v[1]] for v in tri]
                    xc, yc = float(np.mean(cx)), float(np.mean(cy))
                    area = 0.5 * abs((cx[1] - cx[0]) * (cy[2] - cy[0]) - (cx[2] - cx[0]) * (cy[1] - cy[0]))
                    info[tg] = dict(xc=xc, yc=yc, gps=[(xc, yc)] * 4, area=area)
                    elems.append((tg, ic, jc))
    if LST:''')
rep('''                    ops.element("LadrunoLST", tg, *[N(*v) for v in tri], MAT,
                                "-type", "PlaneStrain", "-thick", 1.0, "-body", 0.0, -GAMMA)
''', '''                    if ELEM == "lst":
                        ops.element("LadrunoLST", tg, *[N(*v) for v in tri], MAT,
                                    "-type", "PlaneStrain", "-thick", 1.0, "-body", 0.0, -GAMMA)
                    else:   # BezierTri6: corners then mid12, mid23, mid31 (control points at edge midpoints)
                        ops.element("BezierTri6", tg, *[N(*v) for v in tri], 1.0, "PlaneStrain", MAT,
                                    *(["-bbar"] if ELEM == "bezb" else []), "-bodyForce", 0.0, -GAMMA)
''')
rep('''    for tag, i, j in (elems if not LST else []):
        ops.element("LadrunoQuad", tag, nt(i, j), nt(i + 1, j), nt(i + 1, j + 1),
                    nt(i, j + 1), ns_tag.get(tag, MAT), "-formulation", "bbar",''',
'''    for tag, i, j in (elems if not TRI else []):
        ops.element("LadrunoQuad", tag, nt(i, j), nt(i + 1, j), nt(i + 1, j + 1),
                    nt(i, j + 1), ns_tag.get(tag, MAT), "-formulation", ELEM[1:],''')
open(p, "w").write(s)
print("patched")
