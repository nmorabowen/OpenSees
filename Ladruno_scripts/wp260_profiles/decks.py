"""Deterministic north-star benchmark decks for WP-260 Wave 0 (`260-w0`).

Covers all 5 D1 workload families on `origin/ladruno` @ 7e91aaf1c with ZERO SRC/ changes:
  1. d1a_brick_sanisand     : 3D LadrunoBrick (bbar) + LadrunoSANISAND (SAS-ME), Pardiso, Newmark
  2. d1b_brick_concrete3d   : 3D LadrunoBrick (bbar) + LadrunoConcrete3D, Pardiso, cyclic cracking/damage
  3. d2_cpml_wave           : 3D stdBrick interior + LadrunoCPML absorbing jacket, Pardiso, Newmark
  4. d3a_explicit_cdl_cpml  : Explicit CentralDifferenceLadruno + Diagonal + ssp LadrunoBrick + LadrunoCPML
  5. d3b_explicit_bathe_cpml: Explicit ExplicitBathe + Diagonal + ssp LadrunoBrick + LadrunoCPML
  6. d3c_explicit_sms_cons  : Explicit CentralDifferenceSMSConsistent + Pardiso + fine-slice LadrunoBrick
  7. d4_mpm_modeA           : Mode-A implicit MPM dynamic cantilever (2D GIMP + Pardiso + Newmark)
  8. d5_mcft_2cld12         : 2CLD12-representative RC column (LadrunoForceBeamColumn -monolithic +
                              LadrunoNDFiber + LadrunoMCFTStrip -cyclicDegradation)
"""
from __future__ import annotations

import math
from typing import Callable, Dict, Tuple


# --------------------------------------------------------------------------- #
#  Shared mesh & helper builders                                              #
# --------------------------------------------------------------------------- #
def _nid(ix: int, iy: int, iz: int, nx: int, ny: int) -> int:
    return 1 + ix + (nx + 1) * (iy + (ny + 1) * iz)


def _build_hex_nodes(ops, xs: list[float], ys: list[float], zs: list[float]) -> Tuple[int, int, int]:
    nx, ny, nz = len(xs) - 1, len(ys) - 1, len(zs) - 1
    for iz, z in enumerate(zs):
        for iy, y in enumerate(ys):
            for ix, x in enumerate(xs):
                ops.node(_nid(ix, iy, iz, nx, ny), float(x), float(y), float(z))
    return nx, ny, nz


def _hex_conn(ix: int, iy: int, iz: int, nx: int, ny: int) -> Tuple[int, ...]:
    n1 = _nid(ix, iy, iz, nx, ny)
    n2 = _nid(ix + 1, iy, iz, nx, ny)
    n3 = _nid(ix + 1, iy + 1, iz, nx, ny)
    n4 = _nid(ix, iy + 1, iz, nx, ny)
    n5 = _nid(ix, iy, iz + 1, nx, ny)
    n6 = _nid(ix + 1, iy, iz + 1, nx, ny)
    n7 = _nid(ix + 1, iy + 1, iz + 1, nx, ny)
    n8 = _nid(ix, iy + 1, iz + 1, nx, ny)
    return (n1, n2, n3, n4, n5, n6, n7, n8)


# --------------------------------------------------------------------------- #
#  D1a: Implicit 3D LadrunoBrick + LadrunoSANISAND (SAS-ME)                  #
# --------------------------------------------------------------------------- #
def run_d1a_brick_sanisand(ops, *, smoke: bool = False) -> Dict:
    """Soil block (kPa, m, Mg -> kN) under confining pressure then cyclic shear."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)

    nx, ny, nz = (4, 4, 8) if smoke else (6, 6, 10)
    Lx, Ly, Lz = 4.0, 4.0, 8.0
    xs = [Lx * i / nx for i in range(nx + 1)]
    ys = [Ly * j / ny for j in range(ny + 1)]
    zs = [Lz * k / nz for k in range(nz + 1)]
    _build_hex_nodes(ops, xs, ys, zs)

    # Laminar shear-box boundary condition: couple boundary & horizontal layer DOFs (1, 3) to corner column
    for iz in range(nz + 1):
        for iy in range(ny + 1):
            for ix in range(nx + 1):
                nd = _nid(ix, iy, iz, nx, ny)
                if iz == 0:
                    ops.fix(nd, 1, 1, 1)
                elif iy in (0, ny):
                    ops.fix(nd, 0, 1, 0)
    for iz in range(1, nz + 1):
        master = _nid(0, 0, iz, nx, ny)
        for iy in range(ny + 1):
            for ix in range(nx + 1):
                nd = _nid(ix, iy, iz, nx, ny)
                if nd != master:
                    ops.equalDOF(master, nd, 1, 3)

    rho = 2.0
    san_params = [264.32, 0.3129, 0.6944, 1.33090, 0.71, 0.027, 0.83, 0.45, 101.0,
                  0.005, 1.3, 0.968, 3.5, 0.05, 5.75, 12.5, 1100.0, rho]
    ops.nDMaterial("LadrunoSANISAND", 1, *san_params, 129, 0, 1, 1e-7, 1e-4)

    etag = 1
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                ops.element("LadrunoBrick", etag, *_hex_conn(ix, iy, iz, nx, ny), 1,
                            "-formulation", "bbar", "-b", 0.0, 0.0, -rho * 9.81)
                etag += 1

    # Stage 0 elastic overburden initialization
    ops.updateMaterialStage("-material", 1, "-stage", 0)
    top_area = (Lx / nx) * (Ly / ny)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for iy in range(ny):
        for ix in range(nx):
            f_node = -150.0 * top_area / 4.0
            for nd in (_nid(ix, iy, nz, nx, ny), _nid(ix + 1, iy, nz, nx, ny),
                       _nid(ix + 1, iy + 1, nz, nx, ny), _nid(ix, iy + 1, nz, nx, ny)):
                ops.load(nd, 0.0, 0.0, f_node)

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("Pardiso")
    ops.test("NormDispIncr", 1e-8, 25)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    if ops.analyze(1) != 0:
        raise RuntimeError("d1a_brick_sanisand gravity init failed")

    # Switch to Stage 1 elastoplastic and reset profiler BEFORE nonlinear transient excitation
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    ops.loadConst("-time", 0.0)
    ops.wipeAnalysis()

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("Pardiso")
    ops.test("NormDispIncr", 1e-6, 25)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")

    # Strong sinusoidal base excitation (3.5 m/s^2) on laminar soil column -> drives active SANISAND yielding cleanly
    dt = 0.01
    nsteps = 8 if smoke else 24
    ts_t = [i * dt for i in range(nsteps + 2)]
    ts_a = [3.5 * math.sin(2.0 * math.pi * 2.0 * t) for t in ts_t]
    ops.timeSeries("Path", 2, "-time", *ts_t, "-values", *ts_a)
    ops.pattern("UniformExcitation", 2, 1, "-accel", 2)

    return {
        "kind": "Transient",
        "dt": dt,
        "nsteps": nsteps,
        "nelm": nx * ny * nz,
        "family": "D1a-Brick-SANISAND",
        "setup_fn": None,
    }


# --------------------------------------------------------------------------- #
#  D1b: Implicit 3D LadrunoBrick + LadrunoConcrete3D (CDPM2 cracking/damage)   #
# --------------------------------------------------------------------------- #
def run_d1b_brick_concrete3d(ops, *, smoke: bool = False) -> Dict:
    """3D stub column under compression + cyclic lateral top shear past tensile cracking."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)

    nx, ny, nz = (4, 4, 8) if smoke else (6, 6, 10)
    B, D, H = 300.0, 300.0, 900.0  # mm
    xs = [B * i / nx for i in range(nx + 1)]
    ys = [D * j / ny for j in range(ny + 1)]
    zs = [H * k / nz for k in range(nz + 1)]
    _build_hex_nodes(ops, xs, ys, zs)

    for iy in range(ny + 1):
        for ix in range(nx + 1):
            ops.fix(_nid(ix, iy, 0, nx, ny), 1, 1, 1)

    # Rigid diaphragm cap at top z = H
    top_master = _nid(0, 0, nz, nx, ny)
    for iy in range(ny + 1):
        for ix in range(nx + 1):
            nd = _nid(ix, iy, nz, nx, ny)
            if nd != top_master:
                ops.equalDOF(top_master, nd, 1, 2, 3)

    # Units: N, mm, MPa; rho = 2.4e-9 tonne/mm^3
    ops.nDMaterial("LadrunoConcrete3D", 1, 30000.0, 0.2, 35.0, 3.0, 0.14, 35.0,
                   "-lch", H / nz, "-rho", 2.4e-9)

    etag = 1
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                ops.element("LadrunoBrick", etag, *_hex_conn(ix, iy, iz, nx, ny), 1,
                            "-formulation", "bbar")
                etag += 1

    # Axial compression 5% fc*A = 157.5 kN first (unprofiled preload)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(top_master, 0.0, 0.0, -0.05 * 35.0 * B * D)

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("Pardiso")
    ops.test("NormDispIncr", 1e-7, 25)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    if ops.analyze(1) != 0:
        raise RuntimeError("d1b_brick_concrete3d axial preload failed")

    ops.loadConst("-time", 0.0)
    ops.wipeAnalysis()

    # Cyclic lateral force -> DisplacementControl on top_master DOF 1 (0 -> +1.5mm -> -1.5mm)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    ops.load(top_master, 1.0e5, 0.0, 0.0)

    nsteps = 8 if smoke else 16
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("Pardiso")
    ops.test("NormDispIncr", 1e-6, 30)
    ops.algorithm("KrylovNewton")
    ops.integrator("DisplacementControl", top_master, 1, 0.20 if not smoke else 0.15)
    ops.analysis("Static")

    return {
        "kind": "Static",
        "dt": None,
        "nsteps": nsteps,
        "nelm": nx * ny * nz,
        "family": "D1b-Brick-Concrete3D",
        "setup_fn": None,
    }


# --------------------------------------------------------------------------- #
#  D2: Implicit 3D stdBrick + LadrunoCPML Wave Propagation                    #
# --------------------------------------------------------------------------- #
def _build_cpml_domain(ops, *, n_int: int, n_pml: int, lumped: bool, explicit_mat: bool) -> Tuple[int, int]:
    """3D domain centered at x,y in [-L,L], z in [-L,0] with free surface z=0 and PML jacket."""
    dx = 0.5
    nx = n_int + 2 * n_pml
    ny = n_int + 2 * n_pml
    nz = n_int + n_pml  # free surface at z = 0 (iz = nz)
    x0 = -0.5 * nx * dx
    y0 = -0.5 * ny * dx
    z0 = -nz * dx
    xs = [x0 + i * dx for i in range(nx + 1)]
    ys = [y0 + j * dx for j in range(ny + 1)]
    zs = [z0 + k * dx for k in range(nz + 1)]
    _build_hex_nodes(ops, xs, ys, zs)

    # Outer Dirichlet boundary (bottom z=z0, sides x=x0/x1, y=y0/y1)
    for iz in range(nz + 1):
        for iy in range(ny + 1):
            for ix in range(nx + 1):
                if iz == 0 or ix == 0 or ix == nx or iy == 0 or iy == ny:
                    ops.fix(_nid(ix, iy, iz, nx, ny), 1, 1, 1)

    E, nu, rho = 2.5e8, 0.25, 2000.0
    if explicit_mat:
        # Nonlinear J2 core near yield + SSP brick
        ops.nDMaterial("LadrunoJ2", 1, E / (3.0 * (1.0 - 2.0 * nu)), E / (2.0 * (1.0 + nu)),
                       "-iso", "voce", 5.0e4, 1.0e4, 5.0, 1.0e5, "-rho", rho)
    else:
        ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)

    xmin, xmax = -0.5 * n_int * dx, 0.5 * n_int * dx
    ymin, ymax = -0.5 * n_int * dx, 0.5 * n_int * dx
    zmin, zmax = -n_int * dx, 0.0
    Lpml = n_pml * dx
    cpml_opts = ["-E", E, "-nu", nu, "-rho", rho,
                 "-box", xmin, xmax, ymin, ymax, zmin, zmax,
                 "-L", Lpml, "-m", 2.5, "-R0", 1.5e-3, "-kappaMax", 1.0, "-alphaMax", 15.0]
    if lumped:
        cpml_opts.append("-lumped")

    etag = 1
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                in_int = (n_pml <= ix < nx - n_pml and
                          n_pml <= iy < ny - n_pml and
                          iz >= n_pml)
                conn = _hex_conn(ix, iy, iz, nx, ny)
                if in_int:
                    if explicit_mat:
                        ops.element("LadrunoBrick", etag, *conn, 1, "-formulation", "ssp", "-lumped")
                    else:
                        ops.element("stdBrick", etag, *conn, 1)
                else:
                    ops.element("LadrunoCPML", etag, *conn, *cpml_opts)
                etag += 1

    # Point Ricker-like vertical surface load at (0,0,0) -> (nx//2, ny//2, nz)
    src_node = _nid(nx // 2, ny // 2, nz, nx, ny)
    return src_node, etag - 1


def run_d2_cpml_wave(ops, *, smoke: bool = False) -> Dict:
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    src_node, nelm = _build_cpml_domain(ops, n_int=4 if smoke else 8, n_pml=2,
                                        lumped=False, explicit_mat=False)

    dt = 0.0005
    nsteps = 6 if smoke else 30
    t = [i * dt for i in range(nsteps + 2)]
    vals = [1.0e6 * math.sin(math.pi * ti / (6.0 * dt)) if ti <= 6.0 * dt else 0.0 for ti in t]
    ops.timeSeries("Path", 1, "-time", *t, "-values", *vals)
    ops.pattern("Plain", 1, 1)
    ops.load(src_node, 0.0, 0.0, -1.0)

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("Pardiso")
    ops.test("NormDispIncr", 1e-8, 15)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    return {
        "kind": "Transient",
        "dt": dt,
        "nsteps": nsteps,
        "nelm": nelm,
        "family": "D2-CPML-Implicit",
        "setup_fn": None,
    }


# --------------------------------------------------------------------------- #
#  D3a & D3b: Explicit Wave + Nonlinear Soil/SSP Brick + CPML                 #
# --------------------------------------------------------------------------- #
def _setup_explicit_cpml(ops, integrator: str, *, smoke: bool, full_steps: int = 100) -> Dict:
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    src_node, nelm = _build_cpml_domain(ops, n_int=6 if smoke else 8, n_pml=2,
                                        lumped=True, explicit_mat=True)
    dt = 2.0e-4
    nsteps = 30 if smoke else full_steps
    t = [i * dt for i in range(nsteps + 5)]
    vals = [1.0e5 * math.sin(math.pi * ti / (20.0 * dt)) if ti <= 20.0 * dt else 0.0 for ti in t]
    ops.timeSeries("Path", 1, "-time", *t, "-values", *vals)
    ops.pattern("Plain", 1, 1)
    ops.load(src_node, 0.0, 0.0, -1.0)

    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-8, 1)
    ops.algorithm("Linear")
    ops.integrator(integrator)
    ops.analysis("Transient")
    return {
        "kind": "Transient",
        "dt": dt,
        "nsteps": nsteps,
        "nelm": nelm,
        "family": f"D3-{integrator}",
        "setup_fn": None,
    }


def run_d3a_explicit_cdl_cpml(ops, *, smoke: bool = False) -> Dict:
    return _setup_explicit_cpml(ops, "CentralDifferenceLadruno", smoke=smoke, full_steps=100)


def run_d3b_explicit_bathe_cpml(ops, *, smoke: bool = False) -> Dict:
    return _setup_explicit_cpml(ops, "ExplicitBathe", smoke=smoke, full_steps=80)


# --------------------------------------------------------------------------- #
#  D3c: Explicit Consistent-Mass Scaling (CentralDifferenceSMSConsistent)     #
# --------------------------------------------------------------------------- #
def run_d3c_explicit_sms_cons(ops, *, smoke: bool = False) -> Dict:
    """Explicit SMS consistent-mass wave propagation with a thin fine slice to force M_bar."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)

    nx, ny = (4, 4) if smoke else (8, 8)
    # One thin interior layer (dz = 0.05) flanked by bulk layers (dz = 1.0)
    dzs = [1.0, 1.0, 1.0, 0.05, 1.0, 1.0, 1.0] if smoke else [1.0, 1.0, 1.0, 1.0, 1.0, 0.05, 1.0, 1.0, 1.0, 1.0, 1.0]
    nz = len(dzs)
    xs = [1.0 * i for i in range(nx + 1)]
    ys = [1.0 * j for j in range(ny + 1)]
    zs = [0.0]
    for dz in dzs:
        zs.append(zs[-1] + dz)
    _build_hex_nodes(ops, xs, ys, zs)

    for iy in range(ny + 1):
        for ix in range(nx + 1):
            ops.fix(_nid(ix, iy, 0, nx, ny), 1, 1, 1)

    E, nu, rho = 1.0e8, 0.25, 2000.0
    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)
    etag = 1
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                ops.element("LadrunoBrick", etag, *_hex_conn(ix, iy, iz, nx, ny), 1, "-formulation", "bbar")
                etag += 1

    # Seed small initial velocity on top plane so wave propagates cleanly
    top_vel_z = -0.1
    for iy in range(ny + 1):
        for ix in range(nx + 1):
            ops.setNodeVel(_nid(ix, iy, nz, nx, ny), 3, top_vel_z, "-commit")

    dt_target = 0.003   # forces consistent mass scaling on the 0.05 m slice
    dt = 0.0008
    nsteps = 15 if smoke else 250

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("Pardiso")
    ops.test("NormDispIncr", 1e-10, 1)
    ops.algorithm("Linear")
    ops.integrator("CentralDifferenceSMSConsistent", dt_target)
    ops.analysis("Transient")
    return {
        "kind": "Transient",
        "dt": dt,
        "nsteps": nsteps,
        "nelm": etag - 1,
        "family": "D3c-SMS-Consistent",
        "setup_fn": None,
    }


# --------------------------------------------------------------------------- #
#  D4: Mode-A Implicit MPM Dynamic Cantilever                                  #
# --------------------------------------------------------------------------- #
def run_d4_mpm_modeA(ops, *, smoke: bool = False) -> Dict:
    """2D Mode-A implicit MPM dynamic cantilever (232f cantilever_dyn scaling)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 1.2e8, 0.2, 1000.0)

    h = 0.25 if smoke else 0.10
    cells_x = 44 if smoke else 110
    cells_y = 22 if smoke else 55
    ops.mpm("grid", 232, "-origin", 0.0, -4.0, "-h", h, "-cells", cells_x, cells_y, "-basis", "gimp")
    ops.mpm("mode", "implicit")
    ops.mpm("body", 1, "box", 0.0, 0.0, 10.0, 1.0, "-ppc", 2, "-mat", 1, "-rho", 1000.0)
    ops.mpm("fix", "xmin", "xy")
    ops.mpm("gravity", 0.0, -10.0)

    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("Pardiso")
    ops.test("NormDispIncr", 1e-9, 25)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")

    # Take 1 warm-up step so domainChanged settles before timing steady steps (if desired), or time all
    nsteps = 3 if smoke else 12
    dt = 0.01
    return {
        "kind": "MPM_ModeA",
        "dt": dt,
        "nsteps": nsteps,
        "nelm": 0,
        "family": "D4-MPM-ModeA",
        "setup_fn": None,
    }


# --------------------------------------------------------------------------- #
#  D5: 2CLD12-Representative RC Column (LadrunoForceBeamColumn -monolithic    #
#      + LadrunoNDFiber + LadrunoMCFTStrip -cyclicDegradation)                #
# --------------------------------------------------------------------------- #
def run_d5_mcft_2cld12(ops, *, smoke: bool = False) -> Dict:
    """2CLD12-representative shear-flexure RC column (N, mm, MPa)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)

    # Column dimensions (Sezen & Moehle 2CLD12 style: 457x457 mm, H = 2946 mm, 8 bars #9 ~ 645 mm^2)
    b, h = 457.2, 457.2
    H = 2946.4
    n_ele = 4
    for i in range(n_ele + 1):
        ops.node(i + 1, 0.0, 0.0, H * i / n_ele)
    ops.fix(1, 1, 1, 1, 1, 1, 1)

    fc = 21.1
    ec0 = 0.0020
    # Conc parameters matching test_ladrunoMCFTStrip / Popovics-MCFT (positive compression magnitudes)
    fp_core, ep_core, e20_core = fc * 1.15, 0.0024, 0.018
    fp_cov, ep_cov, e20_cov = fc, ec0, 0.005
    smx, smy = 200.0, 200.0
    rho_tie, fy_tie, fu_tie = 0.0017, 476.0, 710.0

    ops.nDMaterial("LadrunoMCFTStrip", 1, fc, fp_core, ep_core, e20_core,
                   "-crackSpacing", smx, smy, "-ties", rho_tie, fy_tie, fu_tie,
                   "-cyclicDegradation", 300.0)
    ops.nDMaterial("LadrunoMCFTStrip", 2, fc, fp_cov, ep_cov, e20_cov,
                   "-crackSpacing", smx, smy, "-cyclicDegradation", 300.0)
    ops.uniaxialMaterial("LadrunoMCFTSteel", 10, 441.0, 645.0)
    ops.nDMaterial("LadrunoUniaxialFiber", 3, 10)

    n_strips = 24 if smoke else 40
    cov = 38.0
    dy = h / n_strips
    ops.section("LadrunoNDFiber", 1, "-shearProfile", "parabolic", "-shearDepth", h, b,
                "-shearResultant", "modelC", "-noCentroid")
    for s in range(n_strips):
        yc = -0.5 * h + (s + 0.5) * dy
        in_core = abs(yc) <= (0.5 * h - cov)
        mat_c = 1 if in_core else 2
        # Two fibers per strip at z = +/- b/4 so weak-axis bending My is non-singular in 3D
        for zf in (-0.25 * b, 0.25 * b):
            ops.fiber(yc, zf, 0.5 * b * dy, mat_c)
    # 8 #9 bars around perimeter (3 top/bottom + 2 mid-depth)
    for yb, zb in ((-175.0, -175.0), (-175.0, 0.0), (-175.0, 175.0),
                   (0.0, -175.0), (0.0, 175.0),
                   (175.0, -175.0), (175.0, 0.0), (175.0, 175.0)):
        ops.fiber(yb, zb, 645.0, 3)

    ops.geomTransf("Linear", 1, 0.0, 1.0, 0.0)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    for e in range(1, n_ele + 1):
        ops.element("LadrunoForceBeamColumn", e, e, e + 1, 1, 1,
                    "-iter", 30, 1e-12, "-monolithic", "-monoAccept", 10e3, 10e6, 1e-2)

    top_nd = n_ele + 1
    # Axial load 0.15 * fc * Ag = 661 kN (unprofiled static step)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(top_nd, 0.0, 0.0, -0.15 * fc * b * h, 0.0, 0.0, 0.0)

    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("Pardiso")
    ops.test("NormDispIncr", 1e-7, 25)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    if ops.analyze(1) != 0:
        raise RuntimeError("d5_mcft_2cld12 axial load failed")

    ops.loadConst("-time", 0.0)
    ops.wipeAnalysis()

    # Cyclic/monotonic lateral push in x with DisplacementControl past flexural & shear cracking
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    ops.load(top_nd, 1.0e5, 0.0, 0.0, 0.0, 0.0, 0.0)

    nsteps = 12 if smoke else 24
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("Pardiso")
    ops.test("NormDispIncr", 1e-6, 40)
    ops.algorithm("KrylovNewton")
    ops.integrator("DisplacementControl", top_nd, 1, 0.5)
    ops.analysis("Static")

    return {
        "kind": "Static",
        "dt": None,
        "nsteps": nsteps,
        "nelm": n_ele,
        "family": "D5-MCFT-2CLD12",
        "setup_fn": None,
    }


DECKS: Dict[str, Callable] = {
    "d1a_brick_sanisand": run_d1a_brick_sanisand,
    "d1b_brick_concrete3d": run_d1b_brick_concrete3d,
    "d2_cpml_wave": run_d2_cpml_wave,
    "d3a_explicit_cdl_cpml": run_d3a_explicit_cdl_cpml,
    "d3b_explicit_bathe_cpml": run_d3b_explicit_bathe_cpml,
    "d3c_explicit_sms_cons": run_d3c_explicit_sms_cons,
    "d4_mpm_modeA": run_d4_mpm_modeA,
    "d5_mcft_2cld12": run_d5_mcft_2cld12,
}
