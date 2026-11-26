"""Deterministic north-star benchmark decks for WP-260 Wave 0 (`260-w0`), matching D1.

Covers all 5 D1 workload families (`wp-briefs.md` §1) with size parameter `size in {"smoke", "S", "L"}`
and ZERO SRC/ edits:
  - d1a_brick_sanisand     : LadrunoBrick (bbar) + LadrunoSANISAND (SAS-ME), Pardiso +
                             KrylovNewton -reuseFactor 3, Newmark (S, L >= 8x S)
  - d1b_brick_concrete3d   : LadrunoBrick (bbar) + LadrunoConcrete3D, Pardiso +
                             KrylovNewton -reuseFactor 3, DisplacementControl (S, L >= 8x S)
  - d2_cpml_wave           : LadrunoBrick (std) interior + LadrunoCPML jacket + H5DRM input,
                             Pardiso + Linear -factorOnce, Newmark (S, L >= 1e5 DOF)
  - d3a_explicit_cdl_cpml  : CentralDifferenceLadruno, LadrunoBrick -formulation uri -lumped +
                             LadrunoCPML -lumped, Diagonal (S, L >= 8x S)
  - d3b_explicit_bathe_cpml: ExplicitBathe, LadrunoBrick -formulation uri -lumped +
                             LadrunoCPML -lumped, Diagonal (S, L >= 8x S)
  - d3c_explicit_sms_cons  : CentralDifferenceSMSConsistent, Pardiso, LadrunoBrick (bbar) (S)
  - d4_mpm_modeA           : Mode-A implicit MPM dynamic cantilever, Pardiso, Newmark (S)
  - d5_mcft_2cld12         : 2CLD12 RC column (LadrunoForceBeamColumn -monolithic +
                             LadrunoNDFiber + LadrunoMCFTStrip -cyclicDegradation) (S)

Note: This module NEVER imports `h5py` so child processes importing `opensees.pyd` avoid HDF5
library conflicts. Parent processes generate `.h5drm` files and pass `h5drm_path` via args.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple


SIZE_SPECS: Dict[str, Dict[str, Dict]] = {
    "d1a_brick_sanisand": {
        "smoke": {"nx": 4, "ny": 4, "nz": 8, "nsteps": 6},
        "S": {"nx": 6, "ny": 6, "nz": 10, "nsteps": 16},        # 360 hex
        "L": {"nx": 12, "ny": 12, "nz": 20, "nsteps": 10},     # 2,880 hex (8x S)
    },
    "d1b_brick_concrete3d": {
        "smoke": {"nx": 4, "ny": 4, "nz": 8, "nsteps": 6},
        "S": {"nx": 6, "ny": 6, "nz": 10, "nsteps": 16},        # 360 hex
        "L": {"nx": 12, "ny": 12, "nz": 20, "nsteps": 10},     # 2,880 hex (8x S)
    },
    "d2_cpml_wave": {
        "smoke": {"n_int": 6, "n_pml": 2, "nsteps": 6},        # 10x10x8 = 800 hex
        "S": {"n_int": 10, "n_pml": 2, "nsteps": 16},          # 14x14x12 = 2,352 hex
        "L": {"n_int": 30, "n_pml": 3, "nsteps": 8},           # 36x36x33 = 42,768 hex (~134k free DOF >= 1e5)
    },
    "d3a_explicit_cdl_cpml": {
        "smoke": {"n_int": 6, "n_pml": 2, "nsteps": 20},
        "S": {"n_int": 8, "n_pml": 2, "nsteps": 80},           # 12x12x10 = 1,440 hex
        "L": {"n_int": 20, "n_pml": 2, "nsteps": 25},          # 24x24x22 = 12,672 hex (8.8x S)
    },
    "d3b_explicit_bathe_cpml": {
        "smoke": {"n_int": 6, "n_pml": 2, "nsteps": 20},
        "S": {"n_int": 8, "n_pml": 2, "nsteps": 60},           # 1,440 hex
        "L": {"n_int": 20, "n_pml": 2, "nsteps": 20},          # 12,672 hex (8.8x S)
    },
    "d3c_explicit_sms_cons": {
        "smoke": {"nx": 4, "ny": 4, "nz_bulk": 3, "nsteps": 15},
        "S": {"nx": 8, "ny": 8, "nz_bulk": 5, "nsteps": 150},
    },
    "d4_mpm_modeA": {
        "smoke": {"h": 0.25, "cells_x": 44, "cells_y": 22, "nsteps": 4},
        "S": {"h": 0.10, "cells_x": 110, "cells_y": 55, "nsteps": 12},
    },
    "d5_mcft_2cld12": {
        "smoke": {"n_strips": 24, "nsteps": 6},
        "S": {"n_strips": 40, "nsteps": 24},
    },
}

DECK_SIZES: Dict[str, List[str]] = {
    k: [s for s in ("S", "L") if s in v] for k, v in SIZE_SPECS.items()
}


# --------------------------------------------------------------------------- #
#  Shared mesh helpers                                                        #
# --------------------------------------------------------------------------- #
def _nid(ix: int, iy: int, iz: int, nx: int, ny: int) -> int:
    return 1 + ix + (nx + 1) * (iy + (ny + 1) * iz)


def _build_hex_nodes(ops, xs: List[float], ys: List[float], zs: List[float]) -> Tuple[int, int, int]:
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
#  D1a: Implicit 3D LadrunoBrick + LadrunoSANISAND + Pardiso + FactorReuse   #
# --------------------------------------------------------------------------- #
def run_d1a_brick_sanisand(ops, *, size: str = "S", **_kw) -> Dict:
    sp = SIZE_SPECS["d1a_brick_sanisand"][size]
    nx, ny, nz, nsteps = sp["nx"], sp["ny"], sp["nz"], sp["nsteps"]

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)

    Lx, Ly, Lz = 4.0, 4.0, 8.0
    xs = [Lx * i / nx for i in range(nx + 1)]
    ys = [Ly * j / ny for j in range(ny + 1)]
    zs = [Lz * k / nz for k in range(nz + 1)]
    _build_hex_nodes(ops, xs, ys, zs)

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
            f_node = -100.0 * top_area / 4.0
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

    ops.updateMaterialStage("-material", 1, "-stage", 1)
    ops.loadConst("-time", 0.0)
    ops.wipeAnalysis()

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("Pardiso")
    ops.test("NormDispIncr", 1e-6, 45)
    # D1 (1): Pardiso + LadrunoFactorReuse
    ops.algorithm("KrylovNewton", "-maxDim", 20, "-reuseFactor", 3, "-reuseIter", 8)
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")

    dt = 0.01
    ts_t = [i * dt for i in range(nsteps + 4)]
    ts_a = [3.5 * math.sin(2.0 * math.pi * 2.0 * t) for t in ts_t]
    ops.timeSeries("Path", 2, "-time", *ts_t, "-values", *ts_a)
    ops.pattern("UniformExcitation", 2, 1, "-accel", 2)

    ndof_est = 3 * (nx + 1) * (ny + 1) * nz
    return {
        "kind": "Transient",
        "dt": dt,
        "nsteps": nsteps,
        "nelm": nx * ny * nz,
        "ndof_est": ndof_est,
        "size": size,
        "family": "D1a-Brick-SANISAND",
        "constraints": "Transformation",
    }


# --------------------------------------------------------------------------- #
#  D1b: Implicit 3D LadrunoBrick + LadrunoConcrete3D + Pardiso + FactorReuse  #
# --------------------------------------------------------------------------- #
def run_d1b_brick_concrete3d(ops, *, size: str = "S", **_kw) -> Dict:
    sp = SIZE_SPECS["d1b_brick_concrete3d"][size]
    nx, ny, nz, nsteps = sp["nx"], sp["ny"], sp["nz"], sp["nsteps"]

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)

    B, D, H = 300.0, 300.0, 900.0
    xs = [B * i / nx for i in range(nx + 1)]
    ys = [D * j / ny for j in range(ny + 1)]
    zs = [H * k / nz for k in range(nz + 1)]
    _build_hex_nodes(ops, xs, ys, zs)

    for iy in range(ny + 1):
        for ix in range(nx + 1):
            ops.fix(_nid(ix, iy, 0, nx, ny), 1, 1, 1)

    top_master = _nid(0, 0, nz, nx, ny)
    for iy in range(ny + 1):
        for ix in range(nx + 1):
            nd = _nid(ix, iy, nz, nx, ny)
            if nd != top_master:
                ops.equalDOF(top_master, nd, 1, 2, 3)

    ops.nDMaterial("LadrunoConcrete3D", 1, 30000.0, 0.2, 35.0, 3.0, 0.14, 35.0,
                   "-lch", H / nz, "-rho", 2.4e-9)

    etag = 1
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                ops.element("LadrunoBrick", etag, *_hex_conn(ix, iy, iz, nx, ny), 1,
                            "-formulation", "bbar")
                etag += 1

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

    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    ops.load(top_master, 1.0e5, 0.0, 0.0)

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("Pardiso")
    ops.test("NormDispIncr", 1e-6, 60)
    # D1 (1): Pardiso + LadrunoFactorReuse
    ops.algorithm("KrylovNewton", "-maxDim", 20, "-reuseFactor", 3, "-reuseIter", 8)
    ops.integrator("DisplacementControl", top_master, 1, 0.08)
    ops.analysis("Static")

    ndof_est = 3 * (nx + 1) * (ny + 1) * nz
    return {
        "kind": "Static",
        "dt": None,
        "nsteps": nsteps,
        "nelm": nx * ny * nz,
        "ndof_est": ndof_est,
        "size": size,
        "family": "D1b-Brick-Concrete3D",
        "constraints": "Transformation",
    }


# --------------------------------------------------------------------------- #
#  D2 & D3 shared CPML + LadrunoBrick domain builder                          #
# --------------------------------------------------------------------------- #
def d2_h5drm_params(size: str) -> Dict:
    sp = SIZE_SPECS["d2_cpml_wave"][size]
    dx = 0.5
    n_int, n_pml, nsteps = sp["n_int"], sp["n_pml"], sp["nsteps"]
    nx = n_int + 2 * n_pml
    ny = n_int + 2 * n_pml
    nz = n_int + n_pml
    dt = 0.0005
    return {
        "dx": dx,
        "n_int": n_int,
        "n_pml": n_pml,
        "nx": nx,
        "ny": ny,
        "nz": nz,
        "dt": dt,
        "nt": nsteps + 4,
        "x0": -0.5 * nx * dx,
        "y0": -0.5 * ny * dx,
        "z0": -nz * dx,
        "xmin": -0.5 * n_int * dx,
        "xmax": 0.5 * n_int * dx,
        "ymin": -0.5 * n_int * dx,
        "ymax": 0.5 * n_int * dx,
        "zmin": -n_int * dx,
        "zmax": 0.0,
    }


def _build_cpml_domain(ops, *, n_int: int, n_pml: int, lumped: bool,
                       explicit_form: Optional[str]) -> Tuple[int, int, int, int]:
    """3D domain centered in x,y, free surface z=0, with LadrunoBrick core + LadrunoCPML jacket."""
    dx = 0.5
    nx = n_int + 2 * n_pml
    ny = n_int + 2 * n_pml
    nz = n_int + n_pml
    x0 = -0.5 * nx * dx
    y0 = -0.5 * ny * dx
    z0 = -nz * dx
    xs = [x0 + i * dx for i in range(nx + 1)]
    ys = [y0 + j * dx for j in range(ny + 1)]
    zs = [z0 + k * dx for k in range(nz + 1)]
    _build_hex_nodes(ops, xs, ys, zs)

    fixed_nodes = 0
    for iz in range(nz + 1):
        for iy in range(ny + 1):
            for ix in range(nx + 1):
                if iz == 0 or ix == 0 or ix == nx or iy == 0 or iy == ny:
                    ops.fix(_nid(ix, iy, iz, nx, ny), 1, 1, 1)
                    fixed_nodes += 1

    E, nu, rho = 2.5e8, 0.25, 2000.0
    if explicit_form is not None:
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
                    if explicit_form is not None:
                        ops.element("LadrunoBrick", etag, *conn, 1, "-formulation", explicit_form, "-lumped")
                    else:
                        ops.element("LadrunoBrick", etag, *conn, 1, "-formulation", "std")
                else:
                    ops.element("LadrunoCPML", etag, *conn, *cpml_opts)
                etag += 1

    src_node = _nid(nx // 2, ny // 2, nz, nx, ny)
    free_dof = 3 * ((nx + 1) * (ny + 1) * (nz + 1) - fixed_nodes)
    return src_node, etag - 1, free_dof, nz


def run_d2_cpml_wave(ops, *, size: str = "S", h5drm_path: Optional[str] = None,
                     algo: str = "LinearFactorOnce") -> Dict:
    """D1 (2): LadrunoBrick (-std) interior + LadrunoCPML jacket + H5DRM input + Linear -factorOnce."""
    if not h5drm_path:
        raise RuntimeError("d2_cpml_wave requires h5drm_path (parent-generated .h5drm file)")
    sp = SIZE_SPECS["d2_cpml_wave"][size]
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _src, nelm, free_dof, _nz = _build_cpml_domain(
        ops, n_int=sp["n_int"], n_pml=sp["n_pml"], lumped=False, explicit_form=None
    )

    dt = 0.0005
    nsteps = sp["nsteps"]
    T = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    X0 = (0.0, 0.0, 0.0)
    # coordinates in .h5drm are already in model units (meters), crd_scale=1.0
    ops.pattern("H5DRM", 1, h5drm_path, 1.0, 1.0, 1.0e-4, 1, *T, *X0)

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("Pardiso")
    ops.test("NormDispIncr", 1e-8, 15)
    if algo == "Newton":
        ops.algorithm("Newton")
    else:
        ops.algorithm("Linear", "-factorOnce")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    return {
        "kind": "Transient",
        "dt": dt,
        "nsteps": nsteps,
        "nelm": nelm,
        "ndof_est": free_dof,
        "size": size,
        "family": "D2-CPML-Implicit",
        "constraints": "Transformation",
    }


# --------------------------------------------------------------------------- #
#  D3a & D3b: Explicit Wave + LadrunoBrick -formulation uri + CPML -lumped    #
# --------------------------------------------------------------------------- #
def _setup_explicit_cpml(ops, deck_key: str, integrator: str, *, size: str) -> Dict:
    sp = SIZE_SPECS[deck_key][size]
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    src_node, nelm, free_dof, _nz = _build_cpml_domain(
        ops, n_int=sp["n_int"], n_pml=sp["n_pml"], lumped=True, explicit_form="uri"
    )
    dt = 2.0e-4
    nsteps = sp["nsteps"]
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
        "ndof_est": free_dof,
        "size": size,
        "family": f"D3-{integrator}",
        "constraints": "Transformation",
    }


def run_d3a_explicit_cdl_cpml(ops, *, size: str = "S", **_kw) -> Dict:
    return _setup_explicit_cpml(ops, "d3a_explicit_cdl_cpml", "CentralDifferenceLadruno", size=size)


def run_d3b_explicit_bathe_cpml(ops, *, size: str = "S", **_kw) -> Dict:
    return _setup_explicit_cpml(ops, "d3b_explicit_bathe_cpml", "ExplicitBathe", size=size)


# --------------------------------------------------------------------------- #
#  D3c: Explicit Consistent-Mass Scaling (CentralDifferenceSMSConsistent)     #
# --------------------------------------------------------------------------- #
def run_d3c_explicit_sms_cons(ops, *, size: str = "S", **_kw) -> Dict:
    sp = SIZE_SPECS["d3c_explicit_sms_cons"][size]
    nx, ny, nz_bulk, nsteps = sp["nx"], sp["ny"], sp["nz_bulk"], sp["nsteps"]
    dzs = [1.0] * nz_bulk + [0.05] + [1.0] * nz_bulk
    nz = len(dzs)

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
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

    for iy in range(ny + 1):
        for ix in range(nx + 1):
            ops.setNodeVel(_nid(ix, iy, nz, nx, ny), 3, -0.1, "-commit")

    dt_target = 0.003
    dt = 0.0008
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
        "ndof_est": 3 * (nx + 1) * (ny + 1) * nz,
        "size": size,
        "family": "D3c-SMS-Consistent",
        "constraints": "Transformation",
    }


# --------------------------------------------------------------------------- #
#  D4: Mode-A Implicit MPM Dynamic Cantilever                                  #
# --------------------------------------------------------------------------- #
def run_d4_mpm_modeA(ops, *, size: str = "S", **_kw) -> Dict:
    sp = SIZE_SPECS["d4_mpm_modeA"][size]
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 1.2e8, 0.2, 1000.0)

    ops.mpm("grid", 232, "-origin", 0.0, -4.0, "-h", sp["h"],
            "-cells", sp["cells_x"], sp["cells_y"], "-basis", "gimp")
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

    return {
        "kind": "MPM_ModeA",
        "dt": 0.01,
        "nsteps": sp["nsteps"],
        "nelm": 0,
        "ndof_est": 0,
        "size": size,
        "family": "D4-MPM-ModeA",
        "constraints": "Plain",
    }


# --------------------------------------------------------------------------- #
#  D5: 2CLD12 RC Column (LadrunoForceBeamColumn -monolithic + MCFTStrip)      #
# --------------------------------------------------------------------------- #
def run_d5_mcft_2cld12(ops, *, size: str = "S", **_kw) -> Dict:
    sp = SIZE_SPECS["d5_mcft_2cld12"][size]
    n_strips, nsteps = sp["n_strips"], sp["nsteps"]

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)

    b, h, H, n_ele = 457.2, 457.2, 2946.4, 4
    for i in range(n_ele + 1):
        ops.node(i + 1, 0.0, 0.0, H * i / n_ele)
    ops.fix(1, 1, 1, 1, 1, 1, 1)

    fc, ec0 = 21.1, 0.0020
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

    cov = 38.0
    dy = h / n_strips
    ops.section("LadrunoNDFiber", 1, "-shearProfile", "parabolic", "-shearDepth", h, b,
                "-shearResultant", "modelC", "-noCentroid")
    for s in range(n_strips):
        yc = -0.5 * h + (s + 0.5) * dy
        mat_c = 1 if abs(yc) <= (0.5 * h - cov) else 2
        for zf in (-0.25 * b, 0.25 * b):
            ops.fiber(yc, zf, 0.5 * b * dy, mat_c)
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

    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    ops.load(top_nd, 1.0e5, 0.0, 0.0, 0.0, 0.0, 0.0)

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
        "ndof_est": 6 * n_ele,
        "size": size,
        "family": "D5-MCFT-2CLD12",
        "constraints": "Plain",
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
