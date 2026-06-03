"""Ladruno conformal RC member builder (v1a of the RC-3D plan).

The pre-build gates (Ladruno_implementation/21_rc3d_validation_gates.md) showed the
v1 path for regular, grid-aligned reinforced-concrete members is *conformal* meshing:
longitudinal bars + ties as truss elements sharing the concrete solid's nodes
(perfect bond, NO embedded element / NO new C++). This module turns that into a
reusable parametric builder.

Design (pure openseespy; apeGmsh wiring is a later step):
  * rectangular prism column along +z, section bx x by, height H, mesh nx x ny x nz.
  * `ASDConcrete3D` solid (calibrated backbone; emergent confinement, see Gate 2),
    `stdBrick` (swap to LadrunoBrick by name).
  * longitudinal `corotTruss` rebar on every PERIMETER vertical node-line (shared
    nodes -> perfect bond); total area set by the steel ratio rho.
  * optional ties: horizontal perimeter `corotTruss` rings every `tie_every` layers.
  * returns a `BuiltRC` with node grid + base/top sets so the caller sets BCs/loads.
  * `pushover()` convenience runs axial-then-lateral with the blessed adaptive driver
    `ladruno_solve.adaptive_static` (the gates proved fixed-step Newton diverges).

THE CONFORMAL CONSTRAINT (measured in Gate 1): a bar/tie can only sit on a mesh
node-line, so perimeter long. bars = 2*(nx+ny) and tie rings ride the z node planes.
That couples reinforcement layout to mesh density; for non-grid cages / arbitrary hoop
positions, use the embedded element (ADR 20, v2). This builder owns the grid-aligned 80%.

Run:  pythoncore-3.12-64/python.exe ladruno_rc.py     # smoke test (build + pushover)
"""
from __future__ import annotations

import math
import os
import sys
from dataclasses import dataclass, field

try:
    from ladruno_solve import adaptive_static
except ImportError:  # allow running from another cwd
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from ladruno_solve import adaptive_static


# --------------------------------------------------------------------------- spec
@dataclass
class RCColumnSpec:
    bx: float = 0.4               # section width  (x)
    by: float = 0.4               # section depth  (y)
    H: float = 2.0                # height (z)
    nx: int = 2                   # cross-section elements in x
    ny: int = 2                   # cross-section elements in y
    nz: int = 5                   # layers up the height
    # concrete (consistent units; defaults ~ MPa, m)
    fc: float = 30.0e3            # compressive strength (kPa)
    E: float = 30.0e6             # modulus (kPa)
    nu: float = 0.2
    ft: float = 3.0e3             # tensile strength (kPa)
    eps_c0: float = -0.002        # strain at peak fc
    # steel
    Es: float = 200.0e6
    fy: float = 420.0e3
    b_hard: float = 0.01          # Steel02 strain-hardening ratio
    rho_long: float = 0.01        # longitudinal steel ratio (area / Ag)
    # ties
    ties: bool = True
    tie_every: int = 1            # a perimeter tie ring every N z-layers
    rho_tie: float = 0.003        # tie ratio (used to size individual tie bars)
    # element type
    brick: str = "stdBrick"       # or "LadrunoBrick" (+ formulation handled by caller)


@dataclass
class BuiltRC:
    spec: RCColumnSpec
    n_node: dict = field(default_factory=dict)   # (i,j,k) -> tag
    base: list = field(default_factory=list)     # k=0 node tags
    top: list = field(default_factory=list)      # k=nz node tags
    bricks: list = field(default_factory=list)
    longbars: list = field(default_factory=list)
    tiebars: list = field(default_factory=list)
    conc_mat: int = 0
    steel_mat: int = 0

    def nid(self, i, j, k):
        return self.n_node[(i, j, k)]


# ----------------------------------------------------------------- material build
def define_concrete(ops, tag, s: RCColumnSpec):
    """ASDConcrete3D with the Gate-2-calibrated backbone (first compression point =
    elastic limit so the unconfined peak recovers fc; confinement then emerges)."""
    E, fc, ft = s.E, s.fc, s.ft
    CE = [-0.4 * fc / E, s.eps_c0, 2.5 * s.eps_c0, 7.5 * s.eps_c0]
    CS = [-0.4 * fc, -fc, -0.5 * fc, -0.2 * fc]
    CD = [1.0 - abs(CS[i]) / (E * abs(CE[i])) for i in range(len(CE))]
    TE = [ft / E, 1.0e-3]
    TS = [ft, 0.3 * ft / 3.0]
    TD = [0.0, 1.0 - TS[1] / (E * TE[1])]
    ops.nDMaterial("ASDConcrete3D", tag, E, s.nu,
                   "-Te", *TE, "-Ts", *TS, "-Td", *TD,
                   "-Ce", *CE, "-Cs", *CS, "-Cd", *CD)
    return tag


def define_steel(ops, tag, s: RCColumnSpec):
    ops.uniaxialMaterial("Steel02", tag, s.fy, s.Es, s.b_hard, 18, 0.925, 0.15)
    return tag


# --------------------------------------------------------------------- the builder
def build_rc_column(ops, s: RCColumnSpec, *, conc_mat=1, steel_mat=2,
                    node0=1, ele0=1, origin=(0.0, 0.0, 0.0)) -> BuiltRC:
    """Build a conformal RC column into the *current* ops model (caller does
    ops.model('basic','-ndm',3,'-ndf',3) first). Returns a BuiltRC handle."""
    ox, oy, oz = origin
    dx, dy, dz = s.bx / s.nx, s.by / s.ny, s.H / s.nz
    b = BuiltRC(spec=s, conc_mat=conc_mat, steel_mat=steel_mat)

    # nodes
    t = node0
    for k in range(s.nz + 1):
        for j in range(s.ny + 1):
            for i in range(s.nx + 1):
                ops.node(t, ox + i * dx, oy + j * dy, oz + k * dz)
                b.n_node[(i, j, k)] = t
                if k == 0:
                    b.base.append(t)
                if k == s.nz:
                    b.top.append(t)
                t += 1

    define_concrete(ops, conc_mat, s)
    define_steel(ops, steel_mat, s)

    # concrete bricks
    e = ele0
    for k in range(s.nz):
        for j in range(s.ny):
            for i in range(s.nx):
                ops.element(s.brick, e,
                            b.nid(i, j, k), b.nid(i + 1, j, k),
                            b.nid(i + 1, j + 1, k), b.nid(i, j + 1, k),
                            b.nid(i, j, k + 1), b.nid(i + 1, j, k + 1),
                            b.nid(i + 1, j + 1, k + 1), b.nid(i, j + 1, k + 1),
                            conc_mat)
                b.bricks.append(e)
                e += 1

    # longitudinal bars on perimeter vertical node-lines (shared nodes => perfect bond)
    perim = [(i, j) for j in range(s.ny + 1) for i in range(s.nx + 1)
             if i in (0, s.nx) or j in (0, s.ny)]
    a_long = s.rho_long * s.bx * s.by / len(perim)
    for (i, j) in perim:
        for k in range(s.nz):
            ops.element("corotTruss", e, b.nid(i, j, k), b.nid(i, j, k + 1),
                        a_long, steel_mat)
            b.longbars.append(e)
            e += 1

    # ties: perimeter horizontal rings every tie_every layers
    if s.ties:
        a_tie = s.rho_tie * min(s.bx, s.by) * dz / max(1, 2 * (s.nx + s.ny)) * 4
        ring = [(i, j) for (i, j) in
                ([(i, 0) for i in range(s.nx)] +
                 [(s.nx, j) for j in range(s.ny)] +
                 [(i + 1, s.ny) for i in range(s.nx)] +
                 [(0, j + 1) for j in range(s.ny)])]
        ring_next = (lambda i, j: ((i + 1, j) if (j == 0 and i < s.nx) else
                                   (i, j + 1) if (i == s.nx and j < s.ny) else
                                   (i - 1, j) if (j == s.ny and i > 0) else
                                   (i, j - 1)))
        for k in range(0, s.nz + 1, max(1, s.tie_every)):
            for (i, j) in ring:
                ni, nj = ring_next(i, j)
                ops.element("corotTruss", e, b.nid(i, j, k), b.nid(ni, nj, k),
                            a_tie, steel_mat)
                b.tiebars.append(e)
                e += 1
    return b


# ----------------------------------------------------------------- analysis helper
def pushover(ops, b: BuiltRC, *, axial_ratio=0.1, drift_ratio=0.02, nstep=400,
            test_tol=1e-6, test_iter=50):
    """Fix the base, apply constant axial (N/fcAg=axial_ratio) on the top, then a
    lateral (+x) displacement-control pushover via the adaptive driver. Returns
    (result, V_peak, drift_reached). The control node is the first top node."""
    s = b.spec
    for tag in b.base:
        ops.fix(tag, 1, 1, 1)
    N = axial_ratio * s.fc * s.bx * s.by
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    for tag in b.top:
        ops.load(tag, 0, 0, -N / len(b.top))
    ops.constraints("Transformation"); ops.numberer("RCM"); ops.system("UmfPack")
    ops.test("NormDispIncr", test_tol, test_iter, 0); ops.algorithm("KrylovNewton")
    ops.integrator("LoadControl", 0.05); ops.analysis("Static")
    if ops.analyze(20) != 0:
        return None, 0.0, 0.0
    ops.loadConst("-time", 0.0)

    ctrl = b.top[0]
    ops.timeSeries("Linear", 2); ops.pattern("Plain", 2, 2)
    for tag in b.top:
        ops.load(tag, 1.0 / len(b.top), 0, 0)
    drift_max = drift_ratio * s.H
    du = drift_max / nstep
    Vpeak = [0.0]

    def _rec(_frac):
        ops.reactions()
        V = -sum(ops.nodeReaction(t, 1) for t in b.base)
        Vpeak[0] = max(Vpeak[0], abs(V))

    res = adaptive_static(
        ops, lambda sc: ops.integrator("DisplacementControl", ctrl, 1, du * sc),
        n_steps=nstep, on_step=_rec)
    drift = ops.nodeDisp(ctrl, 1) / s.H
    return res, Vpeak[0], drift


# --------------------------------------------------------------------- smoke test
def _bootstrap():
    try:
        import opensees as ops
        return ops
    except ModuleNotFoundError:
        dist = os.environ.get(
            "LADRUNO_DIST",
            r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\dist\bin")
        if os.path.isdir(dist):
            os.add_dll_directory(dist)
            sys.path.insert(0, dist)
        import opensees as ops
        return ops


def _smoke():
    ops = _bootstrap()
    s = RCColumnSpec()                       # 0.4x0.4x2.0, 2x2x5, rho=1%, ties on
    ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
    b = build_rc_column(ops, s)
    nbar = len(b.longbars) // s.nz
    print(f"built: {len(b.bricks)} bricks, {nbar} long bars x {s.nz}, "
          f"{len(b.tiebars)} tie segs, {len(b.n_node) * 3} dof")
    res, Vpeak, drift = pushover(ops, b)
    ok = (res is not None) and bool(res) and Vpeak > 0
    print(f"pushover: {res}")
    print(f"V_peak = {Vpeak:.1f} kN at drift {drift:.4f}")
    print("\nSMOKE:", "PASS — conformal RC column builds + pushes over" if ok
          else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(_smoke())
