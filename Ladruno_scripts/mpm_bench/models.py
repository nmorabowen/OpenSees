"""WP-253: the benchmark models. Every choice and its reason is in README.md.

mode F  column  M2's granular column collapse on a frictional floor (`mpm wall`), 3D quarter column
                (a x b x 2a, `fix xmin x`, `fix ymin y`) and its 2D plane-strain half-column
                analogue (a x 2a, `fix xmin x`). ppc 2 per axis, dt = 0.25 dt_crit.
mode A  block   M1's laterally confined block under self-weight, quasi-static, PARDISO,
                `emptyNodes renormalize`, `massCutoff 5e-2`, gravity ramp over 8 steps; 3D and 2D.

Parameters that M0b already fixed (E, nu, rho, g, the DP cone, the block ramp) are READ from
`mpm_m0b/models.py`, so these numbers continue M2's and M1's. The cell counts are chosen so N is
within 10 % of the nominal; the actual N is recorded by the harness (`mpm info numParticles`).

Units: SI (N, m, s, kg) for ElasticIsotropic and DruckerPrager (as M0b); the kPa / t / m / s deck
of the fork's soil tests for LadrunoSANISAND (tests/test_ladruno_mpm_soil_admission.py,
tests/mpm_soil_cases.txt: rho 2.0, g 9.81, the 18-parameter deck below).
"""
from __future__ import annotations

import math

from benchlib import M0B_MODELS as M0

HOST = 232
TS = 9001
MAT_TAG = 1

BASES = ("gimp", "bspline2")
MATS = ("elastic", "dp", "sanisand")

# --- SANISAND (kPa deck): the fork's test parameters, IntScheme 2 (BackwardEuler_CPPM) ----------
SANISAND_PARAMS = [264.32, 0.3129, 0.6944, 1.3309, 0.71, 0.027, 0.83, 0.45, 101.0, 0.005, 1.3,
                   0.968, 3.5, 0.05, 5.75, 12.5, 1100.0, 2.0]
SANISAND_SCHEME = 2
SAN_RHO, SAN_G = 2.0, 9.81

# --- LadrunoHS (WP-258; kPa deck of tests/test_ladruno_mpm_soil_admission.py) ----------------------
# Reported only, so NOT in MATS: the default matrix is WP-253's. Select it with `--mats hs`.
HS_PARAMS = [35000.0, 35000.0, 105000.0, 0.5, 10.0, 38.0, 8.0, "-nu", 0.2, "-pref", 100.0,
             "-Rf", 0.9, "-OCR", 1.5, "-stage", 0]
EXTRA_MATS = ("hs",)
STAGED = ("sanisand", "hs")      # kPa decks, stage 0 for the warm-up, then the plastic stage

# --- column (mode F) ------------------------------------------------------------------------------
A_COL = M0.A_COL
RHO_COL, G_COL = M0.RHO_COL, M0.G_COL
E_COL, NU_COL = M0.E_COL, M0.NU_COL
PHI_COL, DPC_SIGY = M0.PHI_COL, M0.DPC_SIGY
MU_FLOOR = M0.MU_FLOOR
DT_FACTOR = M0.DT_FACTOR

# --- block (mode A) -------------------------------------------------------------------------------
B = M0.B
BLOCK_RAMP_STEPS, BLOCK_RAMP_END = M0.BLOCK_RAMP_STEPS, M0.BLOCK_RAMP_END
MASS_CUTOFF_A = M0.MASS_CUTOFF_A
SOLVER_A = "Pardiso"


def column_cells(dim, N):
    """Cells (nx, ny, nz) of the column body so that N = 8 nx ny nz (3D) / 4 nx nz (2D) is
    within 10 % of N, with nz = 2 nx (the column is twice as tall as it is wide, as M2's).

    3D: nx and ny are free (the quarter column's y extent b = ny h need not equal a); the pair
    with ny closest to nx wins, then the smaller error. 2D: nx = round(sqrt(N / 8))."""
    if dim == 2:
        nx = max(2, round(math.sqrt(N / 8.0)))
        return nx, 1, 2 * nx
    best = None
    for nx in range(2, 400):
        ny = max(1, round(N / (16.0 * nx * nx)))
        n = 16 * nx * nx * ny
        err = abs(n - N) / N
        if err <= 0.10:
            cand = (abs(ny - nx), err, nx, ny)
            if best is None or cand < best:
                best = cand
    if best is None:
        raise ValueError(f"no column fits N = {N} within 10 %")
    return best[2], best[3], 2 * best[2]


def block_cells(dim, N):
    """Cells per side n of the block: N = 8 n^3 (3D) / 4 n^2 (2D), ppc 2 per axis."""
    return max(2, round((N / 8.0) ** (1.0 / 3.0))) if dim == 3 else max(2, round(math.sqrt(N / 4.0)))


def nominal_particles(mode, dim, N):
    if mode == "F":
        nx, ny, nz = column_cells(dim, N)
        return 8 * nx * ny * nz if dim == 3 else 4 * nx * nz
    n = block_cells(dim, N)
    return 8 * n ** 3 if dim == 3 else 4 * n * n


def material_command(mat, dim, block=False):
    """The nDMaterial command (tag MAT_TAG) and its parameters."""
    if mat == "sanisand":
        return (["nDMaterial", "LadrunoSANISAND", MAT_TAG, *SANISAND_PARAMS, SANISAND_SCHEME],
                {"deck": "kPa (tests/mpm_soil_cases.txt)", "IntScheme": SANISAND_SCHEME})
    if mat == "hs":
        return (["nDMaterial", "LadrunoHS", MAT_TAG, *HS_PARAMS],
                {"deck": "kPa (tests/test_ladruno_mpm_soil_admission.py)"})
    if block:
        E, nu = M0.E_BLOCK, M0.NU_BLOCK
        K, G = M0.block_K_G()
    else:
        E, nu = E_COL, NU_COL
        K = E / (3.0 * (1.0 - 2.0 * nu))
        G = E / (2.0 * (1.0 + nu))
    if mat == "elastic":
        return ["nDMaterial", "ElasticIsotropic", MAT_TAG, E, nu], {"E": E, "nu": nu}
    if mat == "dp":
        if block:
            rho, rhobar, sy = M0.DP_RHO, M0.DP_RHOBAR, M0.DP_SIGY
        else:
            rho, rhobar, sy = M0.dp_rho_from_phi(PHI_COL), 0.0, DPC_SIGY
        return (["nDMaterial", "DruckerPrager", MAT_TAG, K, G, sy, rho, rhobar,
                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                {"K": K, "G": G, "sigmaY": sy, "rho": rho, "rhoBar": rhobar})
    raise ValueError(mat)


def column_commands(dim, N, mat, basis):
    """[(command, args)] of the mode-F column and its parameters."""
    nx, ny, nz = column_cells(dim, N)
    h = A_COL / nx
    san = mat in STAGED
    rho = SAN_RHO if san else RHO_COL
    g = SAN_G if san else G_COL
    mcmd, mpar = material_command(mat, dim)
    cmds = [("wipe", []), ("model", ["basic", "-ndm", dim, "-ndf", dim]), (mcmd[0], mcmd[1:])]
    if dim == 3:
        cmds += [
            ("mpm", ["grid", HOST, "-origin", 0.0, 0.0, 0.0, "-h", h, "-cells", 4 * nx, 4 * ny,
                     nz + 2, "-basis", basis]),
            ("mpm", ["body", 1, "box", 0.0, 0.0, 0.0, A_COL, ny * h, nz * h, "-ppc", 2, "-mat",
                     MAT_TAG, "-rho", rho]),
            ("mpm", ["fix", "xmin", "x"]), ("mpm", ["fix", "ymin", "y"]),
            ("mpm", ["wall", "zmin", MU_FLOOR]),
            ("mpm", ["gravity", 0.0, 0.0, -g]),
        ]
        floor = f"mpm wall zmin {MU_FLOOR}"
    else:
        cmds += [
            ("mpm", ["grid", HOST, "-origin", 0.0, 0.0, "-h", h, "-cells", 4 * nx, nz + 2,
                     "-basis", basis]),
            ("mpm", ["body", 1, "box", 0.0, 0.0, A_COL, nz * h, "-ppc", 2, "-mat", MAT_TAG,
                     "-rho", rho]),
            ("mpm", ["fix", "xmin", "x"]),
            ("mpm", ["wall", "ymin", MU_FLOOR]),
            ("mpm", ["gravity", 0.0, -g]),
        ]
        floor = f"mpm wall ymin {MU_FLOOR}"
    return cmds, {"cells_body": [nx, ny, nz] if dim == 3 else [nx, nz], "h": h, "a": A_COL,
                  "height": nz * h, "rho": rho, "g": g, "basis": basis, "ppc": 2, "floor": floor,
                  "symmetry": "xmin x, ymin y (quarter column)" if dim == 3 else "xmin x (half column)",
                  "material": mpar, "dt_factor": DT_FACTOR,
                  "units": "kPa/t/m/s" if san else "SI"}


def block_commands(dim, N, mat, basis):
    """[(command, args)] of the mode-A block (quasi-static, PARDISO) and its parameters."""
    n = block_cells(dim, N)
    h = B / n
    san = mat in STAGED
    rho = SAN_RHO if san else M0.RHO
    g = SAN_G if san else M0.G_ACC
    mcmd, mpar = material_command(mat, dim, block=True)
    cmds = [("wipe", []), ("model", ["basic", "-ndm", dim, "-ndf", dim]), (mcmd[0], mcmd[1:])]
    if dim == 3:
        cmds += [
            ("mpm", ["grid", HOST, "-origin", 0.0, 0.0, 0.0, "-h", h, "-cells", n, n, n + 2,
                     "-basis", basis]),
            ("mpm", ["mode", "quasiStatic"]),
            ("mpm", ["body", 1, "box", 0.0, 0.0, 0.0, B, B, B, "-ppc", 2, "-mat", MAT_TAG,
                     "-rho", rho]),
            ("mpm", ["fix", "xmin", "x"]), ("mpm", ["fix", "xmax", "x"]),
            ("mpm", ["fix", "ymin", "y"]), ("mpm", ["fix", "ymax", "y"]),
            ("mpm", ["fix", "zmin", "xyz"]),
            ("mpm", ["gravity", 0.0, 0.0, -g]),
        ]
        fixes = "xmin/xmax x, ymin/ymax y, zmin xyz"
    else:
        cmds += [
            ("mpm", ["grid", HOST, "-origin", 0.0, 0.0, "-h", h, "-cells", n, n + 2, "-basis", basis]),
            ("mpm", ["mode", "quasiStatic"]),
            ("mpm", ["body", 1, "box", 0.0, 0.0, B, B, "-ppc", 2, "-mat", MAT_TAG, "-rho", rho]),
            ("mpm", ["fix", "xmin", "x"]), ("mpm", ["fix", "xmax", "x"]),
            ("mpm", ["fix", "ymin", "xy"]),
            ("mpm", ["gravity", 0.0, -g]),
        ]
        fixes = "xmin/xmax x, ymin xy"
    dt = BLOCK_RAMP_END / BLOCK_RAMP_STEPS
    cmds += [
        ("mpm", ["massCutoff", MASS_CUTOFF_A]),
        ("mpm", ["emptyNodes", "renormalize"]),
        ("timeSeries", ["Path", TS, "-time", 0.0, 1.0, 1.0e9, "-values", 0.0, 1.0, 1.0]),
        ("mpm", ["timeSeries", TS]),
        ("constraints", ["Plain"]), ("numberer", ["RCM"]), ("system", [SOLVER_A]),
        ("test", ["NormDispIncr", 1.0e-6, 50]), ("algorithm", ["Newton"]),
        ("integrator", ["LoadControl", dt]), ("analysis", ["Static"]),
    ]
    return cmds, {"cells": n, "h": h, "dt": dt, "B": B, "rho": rho, "g": g, "basis": basis,
                  "ppc": 2, "fixes": fixes, "material": mpar, "massCutoff": MASS_CUTOFF_A,
                  "emptyNodes": "renormalize", "solver": SOLVER_A,
                  "loading": f"gravity ramp 0->{BLOCK_RAMP_END} over {BLOCK_RAMP_STEPS} steps",
                  "units": "kPa/t/m/s" if san else "SI"}
