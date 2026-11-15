"""WP-232g (M0b): the measured models. Every choice and its reason is in README.md.

block   mode A (M1, M3): a laterally confined soft block under self-weight, 3D.
column  mode F (M2, M3): a 3D granular column collapse (quarter column) on a frictional floor.
point   the single-point driver (M2's SANISAND projection): one zero-free-DOF stdBrick.

Units SI (N, m, s, kg) except the SANISAND driver leg, which keeps the fork's SANISAND test
deck in kPa (tests/test_ladruno_sanisand.py, `_PARAMS`).
"""
from __future__ import annotations

import math

HOST = 232
TS = 9001

# --------------------------------------------------------------------------- #
#  block (mode A)                                                             #
# --------------------------------------------------------------------------- #
# n cells per side, ppc 2 per axis: N = 8 n^3.  n = 5, 11, 23 -> 1000, 10648, 97336.
BLOCK_N_CELLS = {1000: 5, 10000: 11, 100000: 23}
BLOCK_SMOKE_MEM_CELLS = {1000: 5, 2000: 6}         # M3 smoke: two sizes for a slope
B = 1.0                    # block side [m]
RHO = 1000.0               # [kg/m3]
G_ACC = 10.0               # [m/s2]
E_BLOCK = 5.0e4            # rho g B / E = 0.2: large enough that particles cross cells at n = 5
NU_BLOCK = 0.3
DP_SIGY = 400.0            # [Pa]: yields below ~43 % of the depth under K0; the top stays elastic
                           # (no apex particle). README: why 400 and not 200
DP_RHO = 0.21              # the 232f footing's cone (phi_txc = 20 deg); yields under K0 loading
DP_RHOBAR = 0.0            # non-associated, as the 232f footing (unsymmetric tangent)
GHOST_GAMMA = 1.0e-2       # the 232f study's ghost spring
MASS_CUTOFF_A = 5.0e-2     # mode A's default (ADR 5.3 amendment 232f), set explicitly for all three


def block_K_G(E=E_BLOCK, nu=NU_BLOCK):
    return E / (3.0 * (1.0 - 2.0 * nu)), E / (2.0 * (1.0 + nu))


def block_period():
    """Quarter-wave period of the confined column, 4 B / c_p, c_p = sqrt(M / rho)."""
    M = E_BLOCK * (1 - NU_BLOCK) / ((1 + NU_BLOCK) * (1 - 2 * NU_BLOCK))
    return 4.0 * B / math.sqrt(M / RHO)


BLOCK_RAMP_STEPS = 8       # QS: the run's eight steps take gravity to BLOCK_RAMP_END
BLOCK_RAMP_END = 0.75      # (the DP block's Newton degrades sharply at full gravity: 44 iterations at 1.0)


def block_dt(mode):
    """QS: a pseudo-time ramp reaching BLOCK_RAMP_END after BLOCK_RAMP_STEPS steps (the smoke
    takes the same increments and stops early). Transient: T/16, so eight steps cover the
    first loading half-cycle of the suddenly applied gravity."""
    if mode == "qs":
        return BLOCK_RAMP_END / BLOCK_RAMP_STEPS
    return block_period() / 16.0


def block_commands(N, mat, mode, solver, treatment, test_tol=1.0e-6, test_iter=50):
    """[(command, args)] for the block, analysis included."""
    n = BLOCK_N_CELLS.get(N) or BLOCK_SMOKE_MEM_CELLS[N]
    h = B / n
    K, G = block_K_G()
    cmds = [("wipe", []), ("model", ["basic", "-ndm", 3, "-ndf", 3])]
    if mat == "elastic":
        cmds.append(("nDMaterial", ["ElasticIsotropic", 1, E_BLOCK, NU_BLOCK]))
    elif mat == "dp":
        cmds.append(("nDMaterial", ["DruckerPrager", 1, K, G, DP_SIGY, DP_RHO, DP_RHOBAR,
                                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
    else:
        raise ValueError(mat)
    cmds += [
        ("mpm", ["grid", HOST, "-origin", 0.0, 0.0, 0.0, "-h", h, "-cells", n, n, n + 2,
                 "-basis", "gimp"]),
        ("mpm", ["mode", "quasiStatic" if mode == "qs" else "implicit"]),
        ("mpm", ["body", 1, "box", 0.0, 0.0, 0.0, B, B, B, "-ppc", 2, "-mat", 1, "-rho", RHO]),
        ("mpm", ["fix", "xmin", "x"]), ("mpm", ["fix", "xmax", "x"]),
        ("mpm", ["fix", "ymin", "y"]), ("mpm", ["fix", "ymax", "y"]),
        ("mpm", ["fix", "zmin", "xyz"]),
        ("mpm", ["gravity", 0.0, 0.0, -G_ACC]),
        ("mpm", ["massCutoff", MASS_CUTOFF_A]),
    ]
    if treatment == "renormalize":
        cmds.append(("mpm", ["emptyNodes", "renormalize"]))
    elif treatment == "exclude":
        cmds.append(("mpm", ["emptyNodes", "exclude"]))
    elif treatment == "ghost":
        cmds.append(("mpm", ["emptyNodes", "ghost", GHOST_GAMMA]))
    else:
        raise ValueError(treatment)
    dt = block_dt(mode)
    if mode == "qs":
        cmds.append(("timeSeries", ["Path", TS, "-time", 0.0, 1.0, 1.0e9, "-values", 0.0, 1.0, 1.0]))
        cmds.append(("mpm", ["timeSeries", TS]))
    cmds += [("constraints", ["Plain"]), ("numberer", ["RCM"]), ("system", [solver]),
             ("test", ["NormDispIncr", test_tol, test_iter]), ("algorithm", ["Newton"])]
    if mode == "qs":
        cmds += [("integrator", ["LoadControl", dt]), ("analysis", ["Static"])]
    else:
        cmds += [("integrator", ["Newmark", 0.5, 0.25]), ("analysis", ["Transient"])]
    return cmds, {"n_cells": n, "h": h, "dt": dt, "B": B, "rho": RHO, "g": G_ACC, "E": E_BLOCK,
                  "nu": NU_BLOCK, "K": K, "G": G,
                  "dp": {"sigmaY": DP_SIGY, "rho": DP_RHO, "rhoBar": DP_RHOBAR} if mat == "dp" else None,
                  "ghost_gamma": GHOST_GAMMA if treatment == "ghost" else None,
                  "massCutoff": MASS_CUTOFF_A, "test": ["NormDispIncr", test_tol, test_iter],
                  "grid_cells": [n, n, n + 2], "ppc": 2, "basis": "gimp",
                  "fixes": "xmin/xmax x, ymin/ymax y, zmin xyz",
                  "loading": f"gravity ramp 0->{BLOCK_RAMP_END} over {BLOCK_RAMP_STEPS} steps" if mode == "qs" else
                             "gravity suddenly applied at t = 0 (Newmark 0.5 0.25)",
                  "period_T": block_period()}


def top_layer(X, h):
    """Indices of the top particle layer (initial z >= B - h/2)."""
    return [p for p in range(len(X) // 3) if X[3 * p + 2] >= B - 0.5 * h]


# --------------------------------------------------------------------------- #
#  column (mode F)                                                            #
# --------------------------------------------------------------------------- #
# quarter column a x a x 2a, m cells across a, ppc 2: N = 16 m^3.  m = 4, 9, 18 -> 1024, 11664, 93312.
COLUMN_M_CELLS = {1000: 4, 10000: 9, 100000: 18}
COLUMN_SMOKE_MEM_CELLS = {1000: 4, 2000: 5}
A_COL = 1.0
RHO_COL = 2000.0
G_COL = 9.81
E_COL = 1.0e6
NU_COL = 0.3
PHI_COL = 30.0                       # DruckerPrager friction angle (triaxial-compression match)
DPC_SIGY = 10.0                      # [Pa] apex regulariser: practically cohesionless
J2_SIGY = 2.0e3                      # [Pa] LadrunoJ2 perfect plasticity
MU_FLOOR = 0.4                       # Coulomb wall on zmin
DT_FACTOR = 0.25                     # dt = 0.25 * criticalTimeStep (the host refuses > 0.5)


def dp_rho_from_phi(phi_deg):
    s = math.sin(math.radians(phi_deg))
    alpha = 2.0 * s / (math.sqrt(3.0) * (3.0 - s))
    return math.sqrt(2.0) * alpha


def column_commands(N, mat):
    m = COLUMN_M_CELLS.get(N) or COLUMN_SMOKE_MEM_CELLS[N]
    h = A_COL / m
    K = E_COL / (3.0 * (1.0 - 2.0 * NU_COL))
    G = E_COL / (2.0 * (1.0 + NU_COL))
    cmds = [("wipe", []), ("model", ["basic", "-ndm", 3, "-ndf", 3])]
    if mat == "elastic":
        cmds.append(("nDMaterial", ["ElasticIsotropic", 1, E_COL, NU_COL]))
    elif mat == "dp":
        cmds.append(("nDMaterial", ["DruckerPrager", 1, K, G, DPC_SIGY, dp_rho_from_phi(PHI_COL), 0.0,
                                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
    elif mat == "j2":
        cmds.append(("nDMaterial", ["LadrunoJ2", 1, K, G, "-iso", "voce", J2_SIGY, 0.0, 0.0, 0.0]))
    else:
        raise ValueError(mat)
    cmds += [
        ("mpm", ["grid", HOST, "-origin", 0.0, 0.0, 0.0, "-h", h, "-cells", 4 * m, 4 * m, 2 * m + 2,
                 "-basis", "gimp"]),
        ("mpm", ["body", 1, "box", 0.0, 0.0, 0.0, A_COL, A_COL, 2.0 * A_COL, "-ppc", 2, "-mat", 1,
                 "-rho", RHO_COL]),
        ("mpm", ["fix", "xmin", "x"]), ("mpm", ["fix", "ymin", "y"]),
        ("mpm", ["wall", "zmin", MU_FLOOR]),
        ("mpm", ["gravity", 0.0, 0.0, -G_COL]),
    ]
    return cmds, {"m_cells": m, "h": h, "a": A_COL, "height": 2.0 * A_COL, "rho": RHO_COL, "g": G_COL,
                  "E": E_COL, "nu": NU_COL, "K": K, "G": G, "grid_cells": [4 * m, 4 * m, 2 * m + 2],
                  "ppc": 2, "basis": "gimp", "floor": f"mpm wall zmin {MU_FLOOR}",
                  "symmetry": "xmin x, ymin y (quarter column)",
                  "dp": {"phi_deg": PHI_COL, "rho": dp_rho_from_phi(PHI_COL), "rhoBar": 0.0,
                         "sigmaY": DPC_SIGY} if mat == "dp" else None,
                  "j2": {"sigmaY": J2_SIGY, "hardening": 0.0} if mat == "j2" else None,
                  "dt_factor": DT_FACTOR}


# --------------------------------------------------------------------------- #
#  yield fractions (from particle Cauchy stress, rotation-invariant)          #
# --------------------------------------------------------------------------- #
def _dev_norm_I1(s):
    sxx, syy, szz, sxy, syz, szx = s
    I1 = sxx + syy + szz
    p = I1 / 3.0
    dn = math.sqrt((sxx - p) ** 2 + (syy - p) ** 2 + (szz - p) ** 2
                   + 2.0 * (sxy * sxy + syz * syz + szx * szx))
    return dn, I1


def yield_fraction(sig, kind, params, rel=1.0e-3):
    """Fraction of particles on (or within rel of) the yield surface. kind dp: UW
    DruckerPrager f = ||s|| + rho I1 - sqrt(2/3) sigmaY; j2: ||s|| - sqrt(2/3) sigmaY."""
    n = len(sig) // 6
    if n == 0:
        return None
    c = math.sqrt(2.0 / 3.0)
    hit = 0
    for p in range(n):
        dn, I1 = _dev_norm_I1(sig[6 * p:6 * p + 6])
        if kind == "dp":
            f = dn + params["rho"] * I1 - c * params["sigmaY"]
            scale = max(dn, abs(params["rho"] * I1), c * params["sigmaY"])
        else:
            f = dn - c * params["sigmaY"]
            scale = c * params["sigmaY"]
        if f >= -rel * scale:
            hit += 1
    return hit / n


# --------------------------------------------------------------------------- #
#  the single-point driver (M2's per-call costs)                              #
# --------------------------------------------------------------------------- #
# One unit stdBrick, every DOF prescribed (rollers on the negative faces, sp on every positive-face
# DOF): zero free equations, so a step's element work is its 8 material points and nothing else
# worth timing. The deck is the fork's confine-first SANISAND driver (tests/test_ladruno_sanisand.py,
# `_build_confined`): an isotropic compression leg, then an isochoric deviatoric ramp. The SAME
# deck runs every material; only the two strain amplitudes differ, chosen so each material sits
# at the column's stress level (mean stress about 10 kPa: rho g H / 2 = 19.6 kPa vertical at
# mid-height) and so the deviatoric leg yields early. Every step is profiled on its own and
# classified elastic or plastic (DP / J2: the point on its yield surface after the step;
# SANISAND: stage 0 = elastic, stage 1 = plastic; ElasticIsotropic: elastic), so the driver
# gives an elastic and a plastic per-call cost for each material.
SANISAND_PARAMS = [264.32, 0.3129, 0.6944, 1.33090, 0.71, 0.027, 0.83, 0.45, 101.0, 0.005, 1.3,
                   0.968, 3.5, 0.05, 5.75, 12.5, 1100.0, 2.0]          # kPa (the fork's test deck)
POINT_P_TARGET = 10.0e3                    # [Pa] the column's mean-stress level
# isotropic strain per direction: p / (3 K) for the column's K (833 kPa); SANISAND (kPa):
# about 2 p / (3 K(p)) with K(10 kPa) = 6e4 kPa (G ~ sqrt p), rounded
POINT_E_CONF = {"elastic": 4.0e-3, "dp": 4.0e-3, "j2": 4.0e-3, "sanisand": 1.0e-4}
# axial compression of the isochoric leg: DP needs q ~ 12 kPa at p = 10 kPa (3 G e = 58 kPa at
# 5e-2); SANISAND keeps the fork deck's 5e-3
POINT_E_AX = {"elastic": 5.0e-2, "dp": 5.0e-2, "j2": 5.0e-2, "sanisand": 5.0e-3}
POINT_LAT = 0.5
POINT_N_CONF = 50
POINT_N_DEV = 200
_XY = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]


def point_series(e_conf, e_ax, n_conf=POINT_N_CONF, n_dev=POINT_N_DEV, lat=POINT_LAT):
    s_lat = [i / n_conf for i in range(n_conf + 1)]
    s_ax = list(s_lat)
    for i in range(1, n_dev + 1):
        s_lat.append(1.0 - lat * e_ax / e_conf * i / n_dev)
        s_ax.append(1.0 + e_ax / e_conf * i / n_dev)
    s_lat.append(s_lat[-1])           # PathSeries returns 0 past its last point: hold it
    s_ax.append(s_ax[-1])
    return s_lat, s_ax


def point_material(mat):
    """The nDMaterial command of the driver: the column's parameters for ElasticIsotropic,
    DruckerPrager and LadrunoJ2; the fork's SANISAND deck (kPa) for LadrunoSANISAND."""
    K = E_COL / (3.0 * (1.0 - 2.0 * NU_COL))
    G = E_COL / (2.0 * (1.0 + NU_COL))
    if mat == "elastic":
        return ["ElasticIsotropic", 1, E_COL, NU_COL]
    if mat == "dp":
        return ["DruckerPrager", 1, K, G, DPC_SIGY, dp_rho_from_phi(PHI_COL), 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    if mat == "j2":
        return ["LadrunoJ2", 1, K, G, "-iso", "voce", J2_SIGY, 0.0, 0.0, 0.0]
    if mat == "sanisand":
        return ["LadrunoSANISAND", 1, *SANISAND_PARAMS]
    raise ValueError(mat)


def point_on_yield(mat, stress):
    """The driver's elastic / plastic flag after a step (DP and J2 from the point's stress)."""
    if mat == "dp":
        return yield_fraction(stress, "dp", {"rho": dp_rho_from_phi(PHI_COL), "sigmaY": DPC_SIGY}) == 1.0
    if mat == "j2":
        return yield_fraction(stress, "j2", {"sigmaY": J2_SIGY}) == 1.0
    return False


def point_build(ops, mat):
    e_conf, e_ax = POINT_E_CONF[mat], POINT_E_AX[mat]
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.node(4 * k + j + 1, x, y, float(k))
    ops.nDMaterial(*point_material(mat))
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.fix(4 * k + j + 1, 1 if x == 0.0 else 0, 1 if y == 0.0 else 0, 1 if k == 0 else 0)
    s_lat, s_ax = point_series(e_conf, e_ax)
    ops.timeSeries("Path", 1, "-dt", 1.0, "-values", *s_lat)
    ops.timeSeries("Path", 2, "-dt", 1.0, "-values", *s_ax)
    ops.pattern("Plain", 1, 1)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            n = 4 * k + j + 1
            if x == 1.0:
                ops.sp(n, 1, -e_conf)
            if y == 1.0:
                ops.sp(n, 2, -e_conf)
    ops.pattern("Plain", 2, 2)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            if k == 1:
                ops.sp(4 * k + j + 1, 3, -e_conf)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-13, 25, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    return {"e_conf": e_conf, "e_ax": e_ax, "lat": POINT_LAT, "n_conf": POINT_N_CONF,
            "n_dev": POINT_N_DEV, "p_target_Pa": POINT_P_TARGET,
            "deck": "tests/test_ladruno_sanisand.py _build_confined (stdBrick, 0 free DOF)"}
