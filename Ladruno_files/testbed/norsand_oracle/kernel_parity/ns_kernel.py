"""ctypes binding of the LadrunoNORSAND C++ kernel (WP-144 P1a parity harness).

Builds kernel_parity/ns_shim.cpp (an extern "C" shim around
SRC/material/nD/LadrunoNorSandKernel.h) into a shared library with g++ and exposes the
kernel API to Python with O2-shaped inputs (o2_algo.Params, 3x3 tensors):

    lib = build(outdir)                      # path of the .so (None if no compiler)
    K = Kernel(lib)
    K.validate(P)                            # (rc, msg, warn)
    K.initial_state(P, sigma0, v0, pi_i0)    # (rc, state[12], msg); pi_i0 None -> on the surface
    K.step(P, state, deps)                   # dict(state, sigma6, C6, info)
    K.step_fractions(P, state, deps, fr, chain=True)   # detail::step_fractions (O2 api.step_fractions)
    K.stress(P, state), K.elastic_tangent(P, state)

Conventions: tensors are 6 TENSOR components {00,11,22,01,12,02}; the kernel tangent is
C[I][J] = d sigma_I / d eps_J with the shear slots varied symmetrically, i.e.
C[I][J] = C4_ijkk (J normal), C4_ijkl + C4_ijlk (J shear): c4_to_c6 maps O2's 3x3x3x3
tangent onto exactly that.
"""
from __future__ import annotations

import ctypes
import os
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
# NS_KERNEL_INCLUDE: build against a different copy of the header (mutation runs).
KERNEL_DIR = os.environ.get("NS_KERNEL_INCLUDE") or os.path.join(REPO, "SRC", "material", "nD")
KERNEL_H = os.path.join(KERNEL_DIR, "LadrunoNorSandKernel.h")
SHIM = os.path.join(HERE, "ns_shim.cpp")
DRIVER = os.path.join(REPO, "tests", "ladrunonorsand_kernel_check.cpp")

REFUSAL = ["OK", "LOCAL_NOCONV", "LOCAL_LINESEARCH", "PI_NOBRACKET", "PI_NOCONV", "B_NONPOS",
           "P_OR_PI_NONNEG", "NEGATIVE_DLAMBDA", "SUBSTEPS_EXHAUSTED"]
EVALERR = ["", "p_or_pi_nonneg", "pi_nonneg", "B_nonpos", "pi_fold", "pi_nobracket", "pi_noconv",
           "nonfinite", "singular_J"]

I6 = (0, 1, 2, 0, 1, 0)
J6 = (0, 1, 2, 1, 2, 2)

PARAM_FIELDS = ("p0", "kappa_hat", "eps_v0", "mu0", "alpha0", "M", "N", "N_bar", "rho", "rho_bar", "chi",
                "h", "lam_tilde", "v_c0", "e0", "lam_c", "xi", "p_a", "c1", "c2")
CSL = {"paper": 0, "fork": 1}
ZETA = {"WW": 0, "GA": 1}
CAP = {"none": 0, "planar": 1, "smooth": 2}


def compiler():
    return os.environ.get("CXX") or shutil.which("g++")


def build(outdir: str, extra_flags=()) -> str | None:
    """Compile the shim into outdir; returns the library path (None if no g++)."""
    cxx = compiler()
    if not cxx:
        return None
    ext = ".dll" if sys.platform == "win32" else ".so"
    lib = os.path.join(outdir, "libns_kernel" + ext)
    pic = [] if sys.platform == "win32" else ["-fPIC"]
    env_flags = os.environ.get("NS_KERNEL_CXXFLAGS", "").split()      # e.g. -fsanitize=address,undefined
    cmd = [cxx, "-std=c++17", "-O2", "-Wall", "-Wextra", "-Werror", *pic, "-shared", "-I", KERNEL_DIR,
           SHIM, "-o", lib, *env_flags, *extra_flags]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"shim build failed:\n{' '.join(cmd)}\n{res.stdout}\n{res.stderr}")
    return lib


def t6(m) -> np.ndarray:
    """3x3 (symmetric) -> 6 tensor comps {00,11,22,01,12,02}."""
    m = np.asarray(m, float)
    return np.array([m[0, 0], m[1, 1], m[2, 2], 0.5 * (m[0, 1] + m[1, 0]), 0.5 * (m[1, 2] + m[2, 1]),
                     0.5 * (m[0, 2] + m[2, 0])])


def m3(t) -> np.ndarray:
    return np.array([[t[0], t[3], t[5]], [t[3], t[1], t[4]], [t[5], t[4], t[2]]], float)


def c4_to_c6(C4) -> np.ndarray:
    """O2's 3x3x3x3 tangent -> the kernel's 6x6 (shear columns: C4_ijkl + C4_ijlk)."""
    C = np.empty((6, 6))
    for I in range(6):
        i, j = I6[I], J6[I]
        for J in range(6):
            k, l = I6[J], J6[J]
            C[I, J] = C4[i, j, k, k] if J < 3 else C4[i, j, k, l] + C4[i, j, l, k]
    return C


def parse_o2_reason(reason: str):
    """O2 reason '<finest> (substeps exhausted at 2^8)' -> (finest refusal name, eval-error name)."""
    r = reason.split(" (substeps")[0]
    if r.startswith("trial_"):
        return "P_OR_PI_NONNEG", r[len("trial_"):]
    if r.startswith("local_linesearch"):
        return "LOCAL_LINESEARCH", r.split(":", 1)[1] if ":" in r else ""
    if r == "local_noconv":
        return "LOCAL_NOCONV", ""
    if r == "local_singular_J":
        return "LOCAL_NOCONV", "singular_J"
    if r == "negative_dlambda":
        return "NEGATIVE_DLAMBDA", ""
    if r.startswith("local_"):
        e = r[len("local_"):]
        code = {"p_or_pi_nonneg": "P_OR_PI_NONNEG", "pi_nonneg": "P_OR_PI_NONNEG", "B_nonpos": "B_NONPOS",
                "pi_nobracket": "PI_NOBRACKET", "pi_noconv": "PI_NOCONV", "pi_fold": "LOCAL_LINESEARCH"}[e]
        return code, e
    raise ValueError(f"unknown O2 reason {reason!r}")


def _name(table, code):
    """Code -> name; out-of-range codes (e.g. -1: invalid fractions, or the shim's StepInfo/out-parameter
    mismatch marker) are returned as 'INVALID(<code>)' instead of wrapping around."""
    return table[code] if 0 <= code < len(table) else f"INVALID({int(code)})"


def _info(info):
    return dict(refusal=_name(REFUSAL, info[0]), plastic=bool(info[1]), vertex=bool(info[2]),
                cap_active=bool(info[3]), local_iters=int(info[4]), pi_iters=int(info[5]),
                substeps=int(info[6]), finest=_name(REFUSAL, info[7]), finest_sub=_name(EVALERR, info[8]))


_D = ctypes.POINTER(ctypes.c_double)
_I = ctypes.POINTER(ctypes.c_int)


def _dp(a: np.ndarray):
    return a.ctypes.data_as(_D)


class Kernel:
    def __init__(self, lib_path: str):
        self.lib = ctypes.CDLL(lib_path)
        L = self.lib
        L.ns_validate.argtypes = [_D, _I, ctypes.c_char_p, ctypes.c_int, _I]
        L.ns_validate.restype = ctypes.c_int
        L.ns_initial_state.argtypes = [_D, _I, _D, ctypes.c_double, ctypes.c_double, _D, ctypes.c_char_p,
                                       ctypes.c_int]
        L.ns_initial_state.restype = ctypes.c_int
        L.ns_step.argtypes = [_D, _I, _D, _D, _D, _D, _D, _I]
        L.ns_step.restype = ctypes.c_int
        L.ns_step_fractions.argtypes = [_D, _I, _D, _D, _D, ctypes.c_int, ctypes.c_int, _D, _D, _D, _I]
        L.ns_step_fractions.restype = ctypes.c_int
        L.ns_stress.argtypes = [_D, _I, _D, _D]
        L.ns_stress.restype = None
        L.ns_elastic_tangent.argtypes = [_D, _I, _D, _D]
        L.ns_elastic_tangent.restype = None

    @staticmethod
    def params(P):
        d = np.ascontiguousarray([float(getattr(P, f)) for f in PARAM_FIELDS], dtype=np.float64)
        i = np.ascontiguousarray([CSL[P.csl_mode], ZETA[P.zeta], CAP[P.cap]], dtype=np.int32)
        return d, i

    def validate(self, P):
        d, i = self.params(P)
        buf = ctypes.create_string_buffer(1024)
        warn = ctypes.c_int(0)
        rc = self.lib.ns_validate(_dp(d), i.ctypes.data_as(_I), buf, 1024, ctypes.byref(warn))
        return rc, buf.value.decode(), bool(warn.value)

    def initial_state(self, P, sigma0, v0, pi_i0=None):
        d, i = self.params(P)
        s0 = np.ascontiguousarray(t6(sigma0))
        st = np.zeros(12)
        buf = ctypes.create_string_buffer(1024)
        pi = float("nan") if pi_i0 is None else float(pi_i0)
        rc = self.lib.ns_initial_state(_dp(d), i.ctypes.data_as(_I), _dp(s0), float(v0), pi, _dp(st), buf, 1024)
        return rc, st, buf.value.decode()

    def step(self, P, st, deps):
        d, i = self.params(P)
        stn = np.ascontiguousarray(st, dtype=np.float64)
        de = np.ascontiguousarray(t6(deps))
        out = np.zeros(12)
        sig = np.zeros(6)
        C = np.zeros(36)
        info = np.zeros(9, dtype=np.int32)
        self.lib.ns_step(_dp(d), i.ctypes.data_as(_I), _dp(stn), _dp(de), _dp(out), _dp(sig), _dp(C),
                         info.ctypes.data_as(_I))
        return dict(state=out, sigma=sig, C=C.reshape(6, 6), info=_info(info))

    def step_fractions(self, P, st, deps, fractions, chain=True):
        d, i = self.params(P)
        stn = np.ascontiguousarray(st, dtype=np.float64)
        de = np.ascontiguousarray(t6(deps))
        fr = np.ascontiguousarray([float(a) for a in fractions], dtype=np.float64)
        out = np.zeros(12)
        sig = np.zeros(6)
        C = np.zeros(36)
        info = np.zeros(9, dtype=np.int32)
        self.lib.ns_step_fractions(_dp(d), i.ctypes.data_as(_I), _dp(stn), _dp(de), _dp(fr), len(fr), int(chain),
                                   _dp(out), _dp(sig), _dp(C), info.ctypes.data_as(_I))
        return dict(state=out, sigma=sig, C=C.reshape(6, 6), info=_info(info))

    def stress(self, P, st):
        d, i = self.params(P)
        sig = np.zeros(6)
        self.lib.ns_stress(_dp(d), i.ctypes.data_as(_I), _dp(np.ascontiguousarray(st, dtype=np.float64)), _dp(sig))
        return sig

    def elastic_tangent(self, P, st):
        d, i = self.params(P)
        C = np.zeros(36)
        self.lib.ns_elastic_tangent(_dp(d), i.ctypes.data_as(_I),
                                    _dp(np.ascontiguousarray(st, dtype=np.float64)), _dp(C))
        return C.reshape(6, 6)
