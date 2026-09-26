#!/usr/bin/env python3
"""WP-124 fingerprint suite: the eight continuum element shells (LadrunoQuad, LadrunoCST,
LadrunoLST, LadrunoCSTPair, LadrunoBrick, LadrunoBrick20, BezierTri6, BezierTet10).

Imported by ../wp123_undamped/fingerprint.py (`--suite shells`). Every value an analysis can
observe through the Element-contract shell is recorded as an exact float repr, per formulation
variant:

  static      DisplacementControl into the plastic range (LadrunoJ2), Newton; disps + every
              setResponse token after each step (incl. stiff / stiffInitial / material ...)
  initial     ModifiedNewton -initial (getInitialStiff) load-controlled
  linear      algorithm Linear (tangent-only defects are invisible under Newton)
  eigen       -fullGenLapack eigenvalues (getMass + getTangentStiff)
  dyn/*       plastic preload, then Newmark with each Rayleigh factor ALONE (alphaM, betaK,
              betaK0, betaKc) and all four; HHT + Linear-algorithm all-four. A preloaded
              plastic state makes K, Kc and K0 differ, so a swapped factor changes bits.
  ground      UniformExcitation (the ground-motion sign / addInertiaLoadToUnbalance)
  explicit    CentralDifference, undamped and alphaM
  param/*     element forwarding ("G" to all GPs, "material 1 G", addToParameter, pressure,
              "rho"); each followed by a static run so the update is observed

A case whose builder raises records "BUILD-ERROR"; a failed analyze records "FAIL@step".
"""
import math

E, NU = 1000.0, 0.3
KB, GS = E / (3 * (1 - 2 * NU)), E / (2 * (1 + NU))
S0, HISO, RHO = 2.0, 40.0, 2.0

RAY = {
    "aM": (0.3, 0.0, 0.0, 0.0),
    "bK": (0.0, 2.0e-3, 0.0, 0.0),
    "bK0": (0.0, 0.0, 2.0e-3, 0.0),
    "bKc": (0.0, 0.0, 0.0, 2.0e-3),
    "all": (0.3, 1.0e-3, 7.0e-4, 5.0e-4),
}

QUAD = [(0.0, 0.0), (2.0, 0.0), (2.2, 1.1), (-0.1, 1.0)]
TRI = [(0.0, 0.0), (2.0, 0.0), (0.3, 1.2)]
TRI_MIDS = [(0, 1), (1, 2), (2, 0)]
HEX = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
       (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.1, 1.05, 1.1), (0.0, 1.0, 1.0)]
HEX_REG = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
           (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0)]
HEX20_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
               (0, 4), (1, 5), (2, 6), (3, 7)]
TET = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
TET_EDGES = [(0, 1), (1, 2), (0, 2), (0, 3), (2, 3), (1, 3)]

PLANE_TOKENS = ["force", "forces", "globalForce", "dampingForce", "dynamicForce",
                "inertialForce", "stress", "stresses", "strain", "strains",
                "stressPlaneStrain", "stiff", "stiffInitial", "charLength", "Jbar",
                ("material", "1", "stress"), ("material", "1", "strain"),
                ("material", "2", "stress"), ("materialState",)]
SOLID_TOKENS = ["force", "forces", "globalForce", "dampingForce", "dynamicForce",
                "inertialForce", "stress", "stresses", "strain", "strains", "stress3D6",
                "strain3D6", "stiff", "stiffInitial", "charLength", "hourglassEnergy",
                ("material", "1", "stress"), ("material", "2", "strain"), ("materialState",)]


def _mids(verts, edges):
    return [tuple(0.5 * (verts[a][d] + verts[b][d]) for d in range(len(verts[0])))
            for a, b in edges]


def _j2(ops, tag, rho=RHO):
    args = ["LadrunoJ2", tag, KB, GS, "-iso", "voce", S0, 0.0, 0.0, HISO]
    if rho:
        args += ["-rho", rho]
    ops.nDMaterial(*args)


# ---------------------------------------------------------------- builders
# Each returns meta: ndm, nodes (tags), ctrl (node, dof), load (list of (node, vec)),
# tokens, pressure (element accepts a "pressure" parameter).

def _plane(ops, name, coords, fixes, conn_opts, finite=False, ele_rho=False, pressure=True):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for i, (x, y) in enumerate(coords):
        ops.node(i + 1, x, y)
    for n, fx in fixes:
        ops.fix(n, *fx)
    _j2(ops, 1, 0.0 if ele_rho else RHO)
    mat = 1
    if finite:
        ops.nDMaterial("LogStrain2D", 2, 1)
        mat = 2
    n = len(coords)
    conn_opts(mat, list(range(1, n + 1)))
    return mat


def quad(form="std", geom="linear", ele_rho=False):
    def build(ops):
        def mk(mat, nodes):
            opts = ["-formulation", form, "-type", "PlaneStrain", "-thick", 0.8]
            if geom != "linear":
                opts += ["-geom", geom]
            if ele_rho:
                opts += ["-rho", RHO]
            ops.element("LadrunoQuad", 1, *nodes, mat, *opts)
        _plane(ops, "LadrunoQuad", QUAD, [(1, (1, 1)), (2, (0, 1))], mk,
               finite=(geom == "finite"), ele_rho=ele_rho)
        return dict(ndm=2, nodes=[1, 2, 3, 4], ctrl=(3, 1),
                    load=[(3, (1.0, 0.2)), (4, (1.0, -0.1))], tokens=PLANE_TOKENS, pressure=True)
    return build


def cst(geom="linear", ele_rho=False):
    def build(ops):
        def mk(mat, nodes):
            opts = ["-type", "PlaneStrain", "-thick", 0.8]
            if geom != "linear":
                opts += ["-geom", geom]
            if ele_rho:
                opts += ["-rho", RHO]
            ops.element("LadrunoCST", 1, *nodes, mat, *opts)
        _plane(ops, "LadrunoCST", TRI, [(1, (1, 1)), (2, (0, 1))], mk,
               finite=(geom == "finite"), ele_rho=ele_rho)
        return dict(ndm=2, nodes=[1, 2, 3], ctrl=(3, 1),
                    load=[(3, (1.0, 0.3))], tokens=PLANE_TOKENS, pressure=True)
    return build


def lst(geom="linear", ele_rho=False):
    def build(ops):
        def mk(mat, nodes):
            opts = ["-type", "PlaneStrain", "-thick", 0.8]
            if geom != "linear":
                opts += ["-geom", geom]
            if ele_rho:
                opts += ["-rho", RHO]
            ops.element("LadrunoLST", 1, *nodes, mat, *opts)
        _plane(ops, "LadrunoLST", TRI + _mids(TRI, TRI_MIDS),
               [(1, (1, 1)), (2, (0, 1)), (4, (0, 1))], mk,
               finite=(geom == "finite"), ele_rho=ele_rho)
        return dict(ndm=2, nodes=list(range(1, 7)), ctrl=(3, 1),
                    load=[(3, (1.0, 0.3)), (5, (0.5, 0.0))], tokens=PLANE_TOKENS, pressure=True)
    return build


def cstpair(ele_rho=False):
    def build(ops):
        def mk(mat, nodes):
            opts = ["-thick", 0.8]
            if ele_rho:
                opts += ["-rho", RHO]
            ops.element("LadrunoCSTPair", 1, *nodes, mat, *opts)
        _plane(ops, "LadrunoCSTPair", QUAD, [(1, (1, 1)), (2, (0, 1))], mk,
               finite=True, ele_rho=ele_rho)
        return dict(ndm=2, nodes=[1, 2, 3, 4], ctrl=(3, 1),
                    load=[(3, (1.0, 0.2)), (4, (1.0, -0.1))], tokens=PLANE_TOKENS, pressure=True)
    return build


def tri6(extra=()):
    def build(ops):
        def mk(mat, nodes):
            ops.element("BezierTri6", 1, *nodes, 0.8, "PlaneStrain", mat, "-rho", RHO, *extra)
        _plane(ops, "BezierTri6", TRI + _mids(TRI, TRI_MIDS),
               [(1, (1, 1)), (2, (0, 1)), (4, (0, 1))], mk, ele_rho=True)
        return dict(ndm=2, nodes=list(range(1, 7)), ctrl=(3, 1),
                    load=[(3, (1.0, 0.3)), (5, (0.5, 0.0))], tokens=PLANE_TOKENS, pressure=True)
    return build


def _solid(ops, coords, fixed, finite):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y, z) in enumerate(coords):
        ops.node(i + 1, x, y, z)
    for n in fixed:
        ops.fix(n, 1, 1, 1)
    _j2(ops, 1)
    if finite:
        ops.nDMaterial("LogStrain", 2, 1)
        return 2
    return 1


def brick(opts=(), finite=False):
    def build(ops):
        mat = _solid(ops, HEX, [1, 2, 3, 4], finite)
        ops.element("LadrunoBrick", 1, *range(1, 9), mat, *opts)
        meta = dict(ndm=3, nodes=list(range(1, 9)), ctrl=(7, 1),
                    load=[(n, (1.0, 0.2, 0.0)) for n in (5, 6, 7, 8)],
                    tokens=SOLID_TOKENS, pressure=False)
        if "-lumped" in opts:
            # C8: -lumped keeps a CONSISTENT residual inertia against the lumped tangent
            # mass, so Newton converges only linearly in transient runs
            meta["test"] = (1.0e-7, 300)
        return meta
    return build


def brick20(opts=()):
    def build(ops):
        mat = _solid(ops, HEX_REG + _mids(HEX_REG, HEX20_EDGES), [1, 2, 3, 4, 9, 10, 11, 12], False)
        ops.element("LadrunoBrick20", 1, *range(1, 21), mat, *opts)
        # uri on ONE H20 is singular by design (uncontrolled 2x2x2 modes, ADR 72 s2.2):
        # responses, eigen and the explicit lane only
        return dict(ndm=3, nodes=list(range(1, 21)), ctrl=(7, 1),
                    load=[(n, (1.0, 0.2, 0.0)) for n in (5, 6, 7, 8)],
                    tokens=SOLID_TOKENS, pressure=False, implicit="uri" not in opts)
    return build


def tet10(opts=(), finite=False):
    def build(ops):
        mat = _solid(ops, TET + _mids(TET, TET_EDGES), [1, 2, 3, 5, 6, 7], finite)
        ops.element("BezierTet10", 1, *range(1, 11), mat, "-rho", RHO, *opts)
        return dict(ndm=3, nodes=list(range(1, 11)), ctrl=(4, 1),
                    load=[(4, (1.0, 0.3, 0.0)), (8, (0.3, 0.0, 0.0))],
                    tokens=SOLID_TOKENS, pressure=True)
    return build


CASES = {
    "Quad/std": quad(), "Quad/bbar": quad("bbar"), "Quad/ssp": quad("ssp"),
    "Quad/eas": quad("eas"), "Quad/std-finite": quad("std", "finite"),
    "Quad/bbar-finite": quad("bbar", "finite"), "Quad/std-eleRho": quad(ele_rho=True),
    "CST/linear": cst(), "CST/finite": cst("finite"), "CST/eleRho": cst(ele_rho=True),
    "LST/linear": lst(), "LST/finite": lst("finite"), "LST/eleRho": lst(ele_rho=True),
    "CSTPair/matRho": cstpair(), "CSTPair/eleRho": cstpair(True),
    "Brick/std": brick(), "Brick/bbar": brick(("-formulation", "bbar")),
    "Brick/uri": brick(("-formulation", "uri", "-hourglass", "stiffness", 0.05)),
    "Brick/uri-phys": brick(("-formulation", "uri", "-hourglass", "physical", 0.05)),
    "Brick/ssp": brick(("-formulation", "ssp")), "Brick/eas": brick(("-formulation", "eas")),
    "Brick/finite": brick(("-geom", "finite"), finite=True),
    "Brick/corot": brick(("-geom", "corot")), "Brick/lumped": brick(("-lumped",)),
    "Brick20/std": brick20(), "Brick20/uri": brick20(("-formulation", "uri")),
    "Brick20/lumped": brick20(("-lumped",)),
    "Tri6/std": tri6(), "Tri6/bbar": tri6(("-bbar",)), "Tri6/cMass": tri6(("-cMass",)),
    "Tet10/std": tet10(), "Tet10/bbar": tet10(("-bbar",)), "Tet10/cMass": tet10(("-cMass",)),
    "Tet10/finite": tet10(("-geom", "finite"), finite=True),
    "Tet10/corot": tet10(("-geom", "corot")),
}


# ---------------------------------------------------------------- scenarios
def _r(vals):
    if vals is None:
        return ["None"]
    if isinstance(vals, (int, float)):
        vals = [vals]
    return [repr(float(x)) for x in vals]


def _state(ops, meta):
    out = []
    for n in meta["nodes"]:
        out += _r(ops.nodeDisp(n))
    return out


def _responses(ops, meta):
    out = []
    for tok in meta["tokens"]:
        args = (tok,) if isinstance(tok, str) else tok
        try:
            v = ops.eleResponse(1, *args)
        except Exception as exc:          # noqa: BLE001 -- the failure mode IS the datum
            v = None
            out.append(f"{'/'.join(args)}:EXC:{type(exc).__name__}")
            continue
        out.append(f"{'/'.join(args)}:{len(v) if v is not None else -1}")
        out += _r(v)
    return out


def _static_setup(ops, meta, algo=("Newton",), ctrl=True, incr=0.004):
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n, vec in meta["load"]:
        ops.load(n, *vec)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    tol, it = meta.get("test", (1.0e-10, 40))
    if algo[0] == "ModifiedNewton":
        it = max(it, 300)
    ops.test("NormDispIncr", tol, it, 0)
    ops.algorithm(*algo)
    if ctrl:
        ops.integrator("DisplacementControl", meta["ctrl"][0], meta["ctrl"][1], incr)
    else:
        ops.integrator("LoadControl", incr)
    ops.analysis("Static", "-noWarnings")


def _steps(ops, meta, n, dt=None, resp_every=False):
    rec = []
    for k in range(n):
        rc = ops.analyze(1) if dt is None else ops.analyze(1, dt)
        if rc != 0:
            rec.append(f"FAIL@{k}")
            break
        rec += _state(ops, meta)
        if resp_every:
            rec += _responses(ops, meta)
    return rec


def _preload(ops, meta):
    _static_setup(ops, meta)
    rec = _steps(ops, meta, 5)
    ops.loadConst("-time", 0.0)
    ops.wipeAnalysis()
    return rec


def _dynamic(ops, meta, integ, ray, algo="Newton", nsteps=25, dt=0.01):
    rec = _preload(ops, meta)
    ops.rayleigh(*ray)
    ops.timeSeries("Sine", 2, 0.0, 10.0, 0.13)
    ops.pattern("Plain", 2, 2)
    for n, vec in meta["load"]:
        ops.load(n, *[3.0 * x for x in vec])
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", *meta.get("test", (1.0e-10, 40)), 0)
    ops.algorithm(algo)
    ops.integrator(*integ)
    ops.analysis("Transient", "-noWarnings")
    for k in range(nsteps):
        if ops.analyze(1, dt) != 0:
            rec.append(f"FAIL@{k}")
            break
        rec += _state(ops, meta)
        for n in meta["nodes"]:
            rec += _r(ops.nodeVel(n))
    rec += _responses(ops, meta)
    return rec


def _ground(ops, meta):
    ops.rayleigh(0.2, 1.0e-3, 0.0, 0.0)
    ops.timeSeries("Path", 3, "-dt", 0.01, "-values", *[math.sin(0.7 * k) * 5.0 for k in range(40)])
    ops.pattern("UniformExcitation", 3, 1, "-accel", 3)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", *meta.get("test", (1.0e-10, 40)), 0)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient", "-noWarnings")
    rec = []
    for k in range(25):
        if ops.analyze(1, 0.01) != 0:
            rec.append(f"FAIL@{k}")
            break
        rec += _state(ops, meta)
        for n in meta["nodes"]:
            rec += _r(ops.nodeAccel(n))
    rec += _responses(ops, meta)
    return rec


def _explicit(ops, meta, ray):
    ops.rayleigh(*ray)
    ops.timeSeries("Linear", 4)
    ops.pattern("Plain", 4, 4)
    for n, vec in meta["load"]:
        ops.load(n, *vec)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-10, 5, 0)
    ops.algorithm("Linear")
    ops.integrator("CentralDifference")
    ops.analysis("Transient", "-noWarnings")
    return _steps(ops, meta, 30, dt=2.0e-4)


def _param(ops, meta, kind):
    rec = []
    try:
        if kind == "G-all":
            ops.parameter(1, "element", 1, "G")
            ops.updateParameter(1, 0.7 * GS)
        elif kind == "G-gp1":
            ops.parameter(1, "element", 1, "material", "1", "G")
            ops.updateParameter(1, 0.5 * GS)
        elif kind == "addTo":
            ops.parameter(1, "element", 1, "sigmaY")
            ops.addToParameter(1, "element", 1, "Hiso")
            ops.updateParameter(1, 3.0)
        elif kind == "rho":
            ops.parameter(1, "element", 1, "rho")
            ops.updateParameter(1, 2.0 * RHO)
        elif kind == "pressure":
            ops.parameter(1, "element", 1, "pressure")
            ops.updateParameter(1, 0.4)
        elif kind == "materialState":
            ops.parameter(1, "element", 1, "materialState")
            ops.updateParameter(1, 1.0)
        rec.append("param:ok")
    except Exception as exc:              # noqa: BLE001
        rec.append(f"param:EXC:{type(exc).__name__}")
    if kind in ("rho",):
        rec += _ground(ops, meta)
    else:
        _static_setup(ops, meta)
        rec += _steps(ops, meta, 4, resp_every=True)
    return rec


def _eigen(ops, meta):
    try:
        return _r(ops.eigen("-fullGenLapack", 2))
    except Exception as exc:              # noqa: BLE001
        return [f"eigen:EXC:{type(exc).__name__}"]


def run_case(ops, name, build):
    """Return {series_key: [str, ...]} for one element variant."""
    out = {}

    def fresh():
        return build(ops)

    try:
        meta = fresh()
    except Exception as exc:              # noqa: BLE001
        return {f"{name}/BUILD": [f"BUILD-ERROR:{type(exc).__name__}:{exc}"]}

    out[f"{name}/responses0"] = _responses(ops, meta)
    out[f"{name}/eigen"] = _eigen(ops, meta)
    if not meta.get("implicit", True):
        for label in ("off", "aM"):
            meta = fresh()
            ray = (0.0, 0.0, 0.0, 0.0) if label == "off" else RAY["aM"]
            out[f"{name}/explicit/{label}"] = _explicit(ops, meta, ray)
        ops.wipe()
        return out

    meta = fresh()
    _static_setup(ops, meta)
    out[f"{name}/static"] = _steps(ops, meta, 6, resp_every=True)

    meta = fresh()
    _static_setup(ops, meta, algo=("ModifiedNewton", "-initial"))
    out[f"{name}/initial"] = _steps(ops, meta, 6)

    meta = fresh()
    _static_setup(ops, meta, algo=("Linear",), ctrl=False, incr=0.05)
    out[f"{name}/linear"] = _steps(ops, meta, 4, resp_every=True)

    for label, ray in RAY.items():
        meta = fresh()
        out[f"{name}/dyn/newmark/{label}"] = _dynamic(ops, meta, ("Newmark", 0.5, 0.25), ray)
    meta = fresh()
    out[f"{name}/dyn/hht/all"] = _dynamic(ops, meta, ("HHT", 0.8), RAY["all"])
    meta = fresh()
    out[f"{name}/dyn/linear/all"] = _dynamic(ops, meta, ("Newmark", 0.5, 0.25), RAY["all"],
                                             algo="Linear")
    meta = fresh()
    out[f"{name}/ground"] = _ground(ops, meta)
    for label in ("off", "aM"):
        meta = fresh()
        ray = (0.0, 0.0, 0.0, 0.0) if label == "off" else RAY["aM"]
        out[f"{name}/explicit/{label}"] = _explicit(ops, meta, ray)
    kinds = ["G-all", "G-gp1", "addTo", "rho", "materialState"]
    if meta["pressure"]:
        kinds.append("pressure")
    for kind in kinds:
        meta = fresh()
        out[f"{name}/param/{kind}"] = _param(ops, meta, kind)
    ops.wipe()
    return out


def fingerprint(ops, only=None):
    out = {}
    for name, build in CASES.items():
        if only and not any(name.startswith(o) for o in only):
            continue
        out.update(run_case(ops, name, build))
    return out
