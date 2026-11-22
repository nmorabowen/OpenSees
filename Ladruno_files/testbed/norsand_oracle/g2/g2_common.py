"""Shared plumbing of the WP-144 G2 gate (Zone B): the OpenSees bootstrap, the single-element
material-point driver, the O2 -> OpenSees parameter mapping.  Not a test file.

BOOTSTRAP (Ladruno_internal/BUILD_GOTCHAS.md sections 0 and 4).  The pyd is linked against CPython 3.12 and needs
`os.add_dll_directory(dist/bin)` plus `sys.path.insert(0, dist/bin)`.  A stale `tests/opensees.pyd` would shadow
the build (section 4b) but this folder never puts `tests/` on sys.path; `ladruno_build()` returns the build stamp
so a test can assert it is running the binary of the checked-out tree.  LADRUNO_OPENSEES_BIN overrides dist/bin.

DRIVER (what "single element with fully prescribed displacements" means here).  A unit stdBrick (small strain)
or LadrunoBrick -geom finite (large strain) cube, EVERY nodal DOF prescribed, so there are zero free equations
(the tests/test_ladruno_sanisand.py driver).  Node displacements u_i = (F - I) X_i with F uniform, so every one of
the eight Gauss points sees exactly the same strain / F: the strain at the material is the increment-by-increment
history handed in, to ~1e-18 absolute.  One Path time series per prescribed DOF (`-dt 1`), LoadControl(1.0): the
state after analyze(k) is the k-th prescribed value exactly.  A tail repeat of the last value is appended because
a Path series returns 0 at and beyond its last time unless `-useLast` is given.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NORSAND_ORACLE = os.path.dirname(HERE)
REPO = os.path.abspath(os.path.join(NORSAND_ORACLE, "..", "..", ".."))
DIST_BIN = os.environ.get("LADRUNO_OPENSEES_BIN") or os.path.join(REPO, "dist", "bin")

for _p in (NORSAND_ORACLE, os.path.join(NORSAND_ORACLE, "kernel_parity"), HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")


def import_ops():
    """The fresh dist/bin build, or None (the caller skips) when it is absent / unloadable."""
    if not os.path.isfile(os.path.join(DIST_BIN, "opensees.pyd")) and not os.path.isfile(os.path.join(DIST_BIN, "opensees.so")):
        return None
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(DIST_BIN)
    if DIST_BIN not in sys.path:
        sys.path.insert(0, DIST_BIN)
    try:
        import opensees as ops
    except ImportError:
        return None
    return ops


ops = import_ops()

XYZ = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
KERNEL_ORDER = ((0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2))     # {00,11,22,01,12,02} = Voigt {11,22,33,12,23,13}
STRAIN_NATURAL = 10.0       # strain-like quantities: zero floor 1e-6 x 10 = 1e-5, i.e. 1e-15 absolute at a 1e-10 gate
W_SHEAR = np.array([1.0, 1.0, 1.0, 0.5, 0.5, 0.5])                  # tangent shear-column factor (shell: halved)

# O2 Params field -> OpenSees flag
_COMMON = (("p0", "-p0"), ("kappa_hat", "-kappa_hat"), ("eps_v0", "-eps_v0"), ("mu0", "-mu0"),
           ("alpha0", "-alpha0"), ("M", "-M"), ("N", "-N"), ("N_bar", "-N_bar"), ("rho", "-rho"),
           ("rho_bar", "-rho_bar"), ("chi", "-chi"), ("h", "-h"))
_PAPER = (("lam_tilde", "-lambda_tilde"), ("v_c0", "-v_c0"))
_FORK = (("e0", "-e0"), ("lam_c", "-lambda_c"), ("xi", "-xi"), ("p_a", "-p_a"))


def ladruno_build():
    return ops.ladrunoBuild() if ops is not None else None


def t6(m):
    m = np.asarray(m, float)
    return np.array([0.5 * (m[i, j] + m[j, i]) for i, j in KERNEL_ORDER])


def m3(t):
    return np.array([[t[0], t[3], t[5]], [t[3], t[1], t[4]], [t[5], t[4], t[2]]], float)


def norsand_args(P, v0, pi0, sigma0):
    """O2 Params (+ initial state) -> the `nDMaterial LadrunoNorSand` argument list (after the tag).

    Round 3b: the energy option (`-energy HAR -k -g -n -p_a` replaces the five BA06 constants, which are then NOT passed: the
    shell refuses them, sheet 2.4), the p' floor (`-pmin`, always passed so the shell cannot fall back on its own default), and
    pi0 = None -> `-pi0_auto` (the unified rule (S.53) of the shell against O2's `initial_state`).  `-p_a` is passed whenever
    it is read (HAR, or the fork CSL)."""
    a = []
    if P.energy == "HAR":
        a += ["-energy", "HAR", "-k", float(P.k), "-g", float(P.g), "-n", float(P.n_e), "-p_a", float(P.p_a)]
        common = _COMMON[5:]                                    # M N N_bar rho rho_bar chi h
    else:
        common = _COMMON
    for f, flag in common:
        a += [flag, float(getattr(P, f))]
    a += ["-csl", P.csl_mode]
    for f, flag in (_PAPER if P.csl_mode == "paper" else _FORK):
        if P.energy == "HAR" and f == "p_a":
            continue                                            # the one -p_a flag was passed with the energy
        a += [flag, float(getattr(P, f))]
    a += ["-zeta", P.zeta, "-cap", P.cap]
    if P.cap != "none":
        a += ["-c1", float(P.c1), "-c2", float(P.c2)]
    a += ["-pmin", float(P.p_min)]
    a += ["-v0", float(v0)]
    a += ["-pi0_auto"] if pi0 is None else ["-pi0", float(pi0)]
    a += ["-sigma0", *[float(x) for x in t6(sigma0)]]
    return a


def _analysis_small():
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")          # the tangent is non-symmetric (never a symmetric solver)
    ops.test("NormDispIncr", 1.0e-13, 25, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


def prescribe(F_hist):
    """Prescribe u_i(k) = (F_k - I) X_i on every node of the unit cube; F_hist = [F_1, ..., F_N] (3x3)."""
    pid = 0
    I3 = np.eye(3)
    for i, x in enumerate(XYZ):
        X = np.array(x, float)
        for d in range(3):
            vals = [0.0] + [float(((F - I3) @ X)[d]) for F in F_hist]
            vals.append(vals[-1])                       # a Path series is 0 at its last time (see module doc)
            if all(v == 0.0 for v in vals):
                ops.fix(i + 1, *[1 if k == d else 0 for k in range(3)])
                continue
            pid += 1
            ops.timeSeries("Path", pid, "-dt", 1.0, "-values", *vals)
            ops.pattern("Plain", pid, pid)
            ops.sp(i + 1, d + 1, 1.0)


def build_small(mat_args, deps_hist, ele="stdBrick"):
    """unit stdBrick on the nDMaterial LadrunoNorSand tag 1; cumulative TENSOR strain history from deps_hist."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, x in enumerate(XYZ):
        ops.node(i + 1, *[float(c) for c in x])
    ops.nDMaterial("LadrunoNorSand", 1, *mat_args)
    if ele == "stdBrick":
        ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    else:                                  # LadrunoBrick forwards a material refusal at the TRIAL (WP-99 F7)
        ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1, "-formulation", "std", "-geom", "linear")
    E = np.zeros((3, 3))
    hist = []
    for d in deps_hist:
        E = E + np.asarray(d, float)
        hist.append(np.eye(3) + E)          # F - I = E for the small-strain brick: u = E X (symmetric E)
    prescribe(hist)
    _analysis_small()


def mat_response(name, gp=1, ele=1):
    return np.array(ops.eleResponse(ele, "material", gp, name), dtype=float)


# ----------------------------------------------------------------------------------------------
# finite strain: LadrunoBrick -geom finite over nDMaterial LogStrain (inner = LadrunoNorSand, tag 1)
# ----------------------------------------------------------------------------------------------
def build_finite(mat_args, F_hist, inner_tag=1, wrap_tag=2):
    """unit LadrunoBrick -geom finite, material = LogStrain(wrap_tag) over LadrunoNorSand(inner_tag); the uniform
    deformation gradients F_hist = [F_1, ..., F_N] are prescribed (total, from the reference configuration)."""
    build_finite_inner(lambda tag: ops.nDMaterial("LadrunoNorSand", tag, *mat_args), F_hist, inner_tag, wrap_tag)


def build_finite_inner(define_inner, F_hist, inner_tag=1, wrap_tag=2):
    """As build_finite, with ANY inner material: define_inner(tag) issues the `nDMaterial ...` command (a
    non-provider inner such as ElasticIsotropic / LadrunoJ2 exercises LogStrain's unchanged D0-inversion fallback)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, x in enumerate(XYZ):
        ops.node(i + 1, *[float(c) for c in x])
    define_inner(inner_tag)
    ops.nDMaterial("LogStrain", wrap_tag, inner_tag)
    ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, wrap_tag, "-formulation", "std", "-geom", "finite")
    prescribe(F_hist)
    _analysis_small()


def a4_from_stiffness(F):
    """The AB06 spatial tangent a^ep_ijkl (sheet 9.5; Kirchhoff-based, a^ep:E = d/dh [tau (1 + hE)^-T], carrying
    the tau(+)1 term) of the unit brick at its CURRENT uniform F, from the 24 x 24 element stiffness.

    Virtual work for the affine modes u = e_k (x) e_L X (every Gauss point sees the same F, unit volume):
    A_iJkL = u_iJ^T K u_kL = dP_iJ/dF_kL (first elasticity tensor) and a_ijkl = A_iJkL F_jJ F_lL (no J factor):
    verified against O2's (S.34) tangent to 6e-6 on the K2 path (the v / v0 effect), 1e-3 off with an extra /J."""
    K = np.array(ops.eleResponse(1, "stiff"), dtype=float).reshape(24, 24)
    X = np.array(XYZ, float)
    modes = {}
    for k in range(3):
        for L in range(3):
            u = np.zeros(24)
            for n in range(8):
                u[3 * n + k] = X[n, L]
            modes[(k, L)] = u
    A = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for J in range(3):
            for k in range(3):
                for L in range(3):
                    A[i, J, k, L] = modes[(i, J)] @ K @ modes[(k, L)]
    return np.einsum("iJkL,jJ,lL->ijkl", A, F, F)
