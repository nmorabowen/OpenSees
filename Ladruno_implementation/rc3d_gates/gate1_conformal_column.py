"""RC-3D Pre-build GATE 1 — is CONFORMAL meshing viable for an RC column?

The embedded-reinforcement ADR (20) is only justified where conformal meshing (rebar
trusses sharing concrete nodes -> perfect bond, NO new element / NO bond element) is
infeasible. This gate proves the conformal path runs end-to-end for a regular column
AND measures its infeasibility boundary: with shared nodes, a bar can ONLY sit on a
mesh node-line, so the number of longitudinal bars is locked to the cross-section mesh
(perimeter bars = 4*n for an n x n element mesh). Dense cages + ties therefore force
impractically fine conformal meshes -> that is exactly where the embedded element earns
its keep.

What it does: square column (b x b x H), n x n x nz stdBrick concrete (ASDConcrete3D),
vertical Steel02 truss rebar along every perimeter vertical node-line (shared nodes =
perfect bond), constant axial load, then lateral displacement-control pushover. Reports
convergence, peak base shear, and the DOF/element cost — for n = 1/2/3 (4/8/12 bars).
NOTE: total steel is held CONSTANT (a_bar = RHO*B*B/nbar, so nbar*a_bar = RHO*B*B is the
same for every n); this isolates the mesh/bar-placement effect, it does NOT vary steel.

ACCEPTANCE (gate informational, not pass/fail): the conformal model converges and gives
a sane pushover. The TAKEAWAY is the cost coupling: with shared nodes a bar can only sit
on a mesh node-line, so perimeter bars = 4*n and the concrete mesh grows O(N^2 * nz) with
bar count. Pure longitudinal count is NOT the wall (12 grid-bars converge fine); conformal
breaks on (a) bars at NON-grid positions / non-uniform spacing, (b) ties/hoops forcing
combinatorial horizontal meshing, and (c) the O(N^2*nz) mesh-cost blow-up. Those are where
the embedded element (ADR 20) earns its keep; conformal is the v1 path for regular,
grid-aligned cages.

RESULT 2026-06-03 (build f6700775), total steel = RHO*B*B = 0.00160 m^2 for all n:
  n bars elems ndof  Vpeak[kN]
  1   4    25    72    128.1
  2   8    60   162    136.8
  3  12   105   288    129.2
All three converge under implicit pushover; Vpeak is ~flat (constant steel; the ~7% spread
is mesh refinement + corner->perimeter bar redistribution, NOT more steel). Explicit run
is a documented follow-up (mass + dt_cr), not exercised here.

Run:  pythoncore-3.12-64/python.exe gate1_conformal_column.py
"""
import os
import sys


def _bootstrap_opensees():
    try:
        import opensees as ops
        return ops
    except ModuleNotFoundError:
        pass
    dist = os.environ.get(
        "LADRUNO_DIST",
        r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\dist\bin",
    )
    if os.path.isdir(dist):
        os.add_dll_directory(dist)
        sys.path.insert(0, dist)
    import opensees as ops
    return ops


ops = _bootstrap_opensees()

# concrete (same calibrated backbone as gate 2) ------------------------------
E, NU, FC, FT = 30000.0e3, 0.2, 30.0e3, 3.0e3   # kPa, m units
CE = [-0.4 * FC / E, -0.002, -0.005, -0.015]
CS = [-0.4 * FC, -FC, -0.5 * FC, -0.2 * FC]
CD = [1.0 - abs(CS[i]) / (E * abs(CE[i])) for i in range(len(CE))]
TE = [FT / E, 1.0e-3]
TS = [FT, 0.3 * FT / 3.0]
TD = [0.0, 1.0 - (0.3 * FT / 3.0) / (E * 1.0e-3)]
# steel
ES, FY, BS = 200.0e6, 420.0e3, 0.01           # kPa
# geometry
B, H, NZ = 0.4, 2.0, 5
RHO = 0.01                                      # longitudinal steel ratio (approx)
AXIAL_RATIO = 0.10                              # N / (fc Ag)


def _nid(i, j, k, n):
    return (k * (n + 1) + j) * (n + 1) + i + 1


def build(n):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    dx = B / n
    dz = H / NZ
    for k in range(NZ + 1):
        for j in range(n + 1):
            for i in range(n + 1):
                ops.node(_nid(i, j, k, n), i * dx, j * dx, k * dz)
    ops.nDMaterial("ASDConcrete3D", 1, E, NU,
                   "-Te", *TE, "-Ts", *TS, "-Td", *TD,
                   "-Ce", *CE, "-Cs", *CS, "-Cd", *CD)
    ops.uniaxialMaterial("Steel02", 2, FY, ES, BS, 18, 0.925, 0.15)
    e = 1
    for k in range(NZ):
        for j in range(n):
            for i in range(n):
                ops.element("stdBrick", e,
                            _nid(i, j, k, n), _nid(i + 1, j, k, n),
                            _nid(i + 1, j + 1, k, n), _nid(i, j + 1, k, n),
                            _nid(i, j, k + 1, n), _nid(i + 1, j, k + 1, n),
                            _nid(i + 1, j + 1, k + 1, n), _nid(i, j + 1, k + 1, n), 1)
                e += 1
    # perimeter vertical node-lines = bar positions (the conformal constraint)
    perim = [(i, j) for j in range(n + 1) for i in range(n + 1)
             if i in (0, n) or j in (0, n)]
    nbar = len(perim)
    a_bar = RHO * B * B / nbar
    for (i, j) in perim:
        for k in range(NZ):
            ops.element("Truss", e, _nid(i, j, k, n), _nid(i, j, k + 1, n), a_bar, 2)
            e += 1
    # base fixed
    for j in range(n + 1):
        for i in range(n + 1):
            ops.fix(_nid(i, j, 0, n), 1, 1, 1)
    return n, nbar, e - 1, perim


def run(n):
    _, nbar, nele, perim = build(n)
    ndof = len(ops.getNodeTags()) * 3
    top = [_nid(i, j, NZ, n) for j in range(n + 1) for i in range(n + 1)]
    # phase 1: axial load (constant)
    N = AXIAL_RATIO * FC * B * B
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    for t in top:
        ops.load(t, 0, 0, -N / len(top))
    ops.constraints("Transformation"); ops.numberer("RCM"); ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-6, 200, 0); ops.algorithm("KrylovNewton")
    ops.integrator("LoadControl", 0.05); ops.analysis("Static")
    if ops.analyze(20) != 0:
        return dict(n=n, nbar=nbar, nele=nele, ndof=ndof, ok=False, where="axial")
    ops.loadConst("-time", 0.0)
    # phase 2: lateral pushover (displacement control on a top corner, +x)
    ctrl = top[0]
    ops.timeSeries("Linear", 2); ops.pattern("Plain", 2, 2)
    for t in top:
        ops.load(t, 1.0 / len(top), 0, 0)
    drift_max = 0.02 * H
    nstep = 400
    du = drift_max / nstep
    Vpeak = 0.0
    cur = 0.0
    fails = 0
    step = du
    while cur < drift_max - 1e-9 and fails < 12:
        ops.integrator("DisplacementControl", ctrl, 1, step)
        if ops.analyze(1) == 0:
            cur += step
            ops.reactions()
            V = -sum(ops.nodeReaction(_nid(i, j, 0, n), 1)
                     for j in range(n + 1) for i in range(n + 1))
            Vpeak = max(Vpeak, abs(V))
            step = du
        else:
            step *= 0.5
            fails += 1
    return dict(n=n, nbar=nbar, nele=nele, ndof=ndof, ok=(cur > 0.5 * drift_max),
                drift=cur / H, Vpeak=Vpeak)


def main():
    ts = RHO * B * B
    print(f"{'n':>2} {'bars':>4} {'elems':>6} {'ndof':>6} {'conv':>5} "
          f"{'drift':>6} {'Vpeak[kN]':>10} {'tot_steel':>10}")
    for n in (1, 2, 3):
        r = run(n)
        if not r.get("ok"):
            print(f"{r['n']:>2} {r['nbar']:>4} {r['nele']:>6} {r['ndof']:>6}  "
                  f"diverged@{r.get('where','push')}")
            continue
        print(f"{r['n']:>2} {r['nbar']:>4} {r['nele']:>6} {r['ndof']:>6} "
              f"{'yes':>5} {r['drift']:>6.4f} {r['Vpeak']:>10.1f} {ts:>10.5f}")
    print("\nTAKEAWAY: total steel is CONSTANT (RHO*B*B); Vpeak is ~flat across n. With "
          "shared nodes bars=4*n are locked to the mesh, so the concrete mesh grows "
          "O(N^2*nz) with bar count (4->25, 8->60, 12->105 elems). Pure count is NOT the "
          "wall (12 grid-bars converge); conformal breaks on non-grid bar positions, "
          "ties/hoops, and mesh-cost blow-up -> embedding (ADR 20) for those. Conformal is "
          "the v1 path for regular grid-aligned cages.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
