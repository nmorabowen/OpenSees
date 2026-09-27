"""WP-133 (TIMs F23a) -- PDMY03 drained plane-strain decks, shared by the
byte-identity capture and the pytest battery.

One 1x1 ``quad`` (PlaneStrain) of ``PressureDependMultiYield03`` (nd=2):
isotropic confinement ``P0`` by nodal loads in stage 0 (elastic), then
``updateMaterialStage 1`` and a drained, displacement-controlled vertical
compression at constant lateral stress. Bottom on rollers, left on rollers,
top nodes tied in y. Everything is a single Gauss-point-uniform state.

Run as a script to capture the byte-identity baseline::

    python -S tests/wp133_pdmy03_deck.py <dist/bin> <out.json>

It imports nothing from opensees at module import; callers pass ``ops``.
"""
import json
import os
import sys

# dense-ish sand, PDMY03 signature: nd rho G B phi gammaPeak refP d PT
# mType ca cb cc cd ce da db dc   (then the optional positional tail)
BASE = [2, 2.0, 1.3e5, 2.6e5, 40.0, 0.1, 101.0, 0.5, 26.0,
        0, 0.013, 0.0, 0.3, 0.0, 0.0, 0.3, 3.0, 0.0]
TAIL = [20, 1.0, 0.0, 101.0, 0.1]          # NYS liq1 liq2 pa c (defaults except c)

P0 = 100.0        # isotropic confinement (kPa)
DY = -5.0e-4      # top displacement per step (1x1 element -> axial strain)
NSTEP = 100


def build_material(ops, tag, base=None, tail=None, extra=()):
    args = list(base if base is not None else BASE)
    if tail:
        args += list(tail)
    args += list(extra)
    ops.nDMaterial("PressureDependMultiYield03", tag, *args)


def run_quad(ops, mat_tag, node0=0, ele_tag=1, x0=0.0, nstep=NSTEP,
             wipe=True, materials=None):
    """Build (optionally after wipe) and run one drained quad.

    ``materials`` is a callable(ops) that defines the nDMaterial(s) after
    the model command. Returns dict(stress=[...], strain=[...]) per step.
    When several quads share one model (the interference test) the caller
    builds them with ``add_quad`` and drives them with ``drive``.
    """
    if wipe:
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        if materials is not None:
            materials(ops)
    add_quad(ops, mat_tag, node0, ele_tag, x0)
    return drive(ops, [(node0, ele_tag)], nstep)[ele_tag]


def add_quad(ops, mat_tag, node0, ele_tag, x0):
    n = [node0 + k for k in (1, 2, 3, 4)]
    ops.node(n[0], x0 + 0.0, 0.0)
    ops.node(n[1], x0 + 1.0, 0.0)
    ops.node(n[2], x0 + 1.0, 1.0)
    ops.node(n[3], x0 + 0.0, 1.0)
    ops.fix(n[0], 1, 1)
    ops.fix(n[1], 0, 1)
    ops.fix(n[3], 1, 0)
    ops.equalDOF(n[2], n[3], 2)
    ops.element("quad", ele_tag, *n, 1.0, "PlaneStrain", mat_tag)


def drive(ops, quads, nstep):
    # stage 0: elastic confinement
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for node0, _ in quads:
        n2, n3, n4 = node0 + 2, node0 + 3, node0 + 4
        ops.load(n2, -P0 * 0.5, 0.0)
        ops.load(n3, -P0 * 0.5, -P0 * 0.5)
        ops.load(n4, 0.0, -P0 * 0.5)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-9, 100, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("LoadControl", 0.1)
    ops.analysis("Static")
    assert ops.analyze(10) == 0
    ops.loadConst("-time", 0.0)
    ops.integrator("LoadControl", 0.0)
    for _, ele in quads:
        ops.updateMaterialStage("-material", _mat_of(ops, ele), "-stage", 1)
    assert ops.analyze(1) == 0  # settle into stage 1 at the same load
    ops.loadConst("-time", 0.0)

    # drained compression: prescribed top displacement, lateral load held
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for node0, _ in quads:
        ops.sp(node0 + 3, 2, 1.0)
    ops.integrator("LoadControl", DY)
    ops.analysis("Static")
    out = {ele: dict(stress=[], strain=[]) for _, ele in quads}
    for _ in range(nstep):
        rc = ops.analyze(1)
        if rc != 0:
            break
        for _, ele in quads:
            out[ele]["stress"].append(list(ops.eleResponse(ele, "stress")))
            out[ele]["strain"].append(list(ops.eleResponse(ele, "strain")))
    for _, ele in quads:
        out[ele]["n"] = len(out[ele]["stress"])
    return out


_ELE_MAT = {}


def _mat_of(ops, ele):
    return _ELE_MAT.get(ele, ele)


def set_ele_mat(ele, mat):
    _ELE_MAT[ele] = mat


def case_default(ops):
    set_ele_mat(1, 1)
    return run_quad(ops, 1, materials=lambda o: build_material(o, 1))


def case_default_tail(ops):
    set_ele_mat(1, 1)
    return run_quad(ops, 1, materials=lambda o: build_material(o, 1, tail=TAIL))


def case_user_surfaces(ops):
    # negative NYS -> user-defined backbone (r_i, Gs_i) pairs, then tail
    pairs = [1e-4, 0.9, 5e-4, 0.6, 1e-3, 0.4, 5e-3, 0.1, 1e-2, 0.06]
    set_ele_mat(1, 1)
    return run_quad(ops, 1, materials=lambda o: build_material(
        o, 1, tail=[-5] + pairs + [1.0, 0.0, 101.0, 0.1]))


CASES = {
    "default": case_default,
    "default_tail": case_default_tail,
    "user_surfaces": case_user_surfaces,
}


def to_hex(res):
    return dict(n=res["n"],
                stress=[[float(v).hex() for v in row] for row in res["stress"]],
                strain=[[float(v).hex() for v in row] for row in res["strain"]])


def main(dist, out):
    os.add_dll_directory(dist)
    sys.path.insert(0, dist)
    import opensees as ops
    assert os.path.normcase(os.path.dirname(ops.__file__)) == os.path.normcase(dist), ops.__file__
    data = {name: to_hex(fn(ops)) for name, fn in CASES.items()}
    with open(out, "w") as f:
        json.dump(data, f, indent=0)
    for k, v in data.items():
        print(k, v["n"], v["stress"][-1] if v["n"] else None)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
