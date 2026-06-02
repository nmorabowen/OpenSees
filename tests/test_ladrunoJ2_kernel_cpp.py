"""Cross-check the extracted C++ LadrunoJ2 return-map kernel against the numpy
oracle WITHOUT building OpenSees.

Compiles tests/ladrunoj2_kernel_check.cpp (pure — includes only
SRC/material/nD/LadrunoJ2Kernel.h) with g++, runs it over four strain PATHS
(isotropic Voce+linear, single Armstrong-Frederick cyclic, three-term Chaboche
3D-with-shear, linear-kinematic cyclic), and asserts that, step by step:

  * the Cauchy stress and the full internal state (εᵖ, ε̄ᵖ, backstresses, Δγ)
    match the independent 3×3-tensor oracle (tests/ladrunoj2_reference.py), and
  * the algorithmic tangent Dtan matches a finite difference of the oracle stress
    w.r.t. the engineering strain — an independent check of the consistent
    tangent (the oracle does NOT transcribe the analytic tangent formula).

This is the proof that the kernel extraction (the return map moved out of
LadrunoJ2.cpp::integrate into the header-only, finite-strain-shareable
LadrunoJ2Kernel.h) is faithful. Skipped if no C++ compiler is available.
"""
import os
import shutil
import subprocess
import tempfile

import numpy as np
import pytest

import ladrunoj2_reference as lr

pytestmark = [pytest.mark.zone_a, pytest.mark.t1]

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC_ND = os.path.normpath(os.path.join(_HERE, "..", "SRC", "material", "nD"))

K, G = 1500.0, 700.0
MONO = [0.002, 0.005, 0.01, 0.02, 0.03, 0.05]
CYC = [0.03, 0.0, -0.03, 0.0, 0.03]

# Non-proportional path (8 x 6 tensor comps): axial then shear — must match the .cpp.
NONPROP = [
    [0.010, 0, 0, 0, 0, 0], [0.020, 0, 0, 0, 0, 0],
    [0.030, 0, 0, 0, 0, 0], [0.040, 0, 0, 0, 0, 0],
    [0.040, 0, 0, 0.010, 0, 0], [0.040, 0, 0, 0.020, 0, 0],
    [0.040, 0, 0, 0.030, 0, 0], [0.040, 0, 0, 0.040, 0, 0],
]

# (name, Params, direction[6 tensor comps], multipliers, explicit_path) — must
# match the .cpp. explicit_path (or None) is an nstep list of 6 tensor-comp
# vectors and overrides direction*mult.
SCEN = [
    ("elastic_iso",
     lr.Params(K, G, sig0=5.0, Qinf=4.0, bIso=30.0, Hiso=10.0),
     [1, 0, 0, 0, 0, 0], MONO, None),
    ("single_af_cyclic",
     lr.Params(K, G, sig0=5.0, C=[400.0], gam=[80.0]),
     [1, 0, 0, 0, 0, 0], CYC, None),
    ("chaboche3_3d",
     lr.Params(K, G, sig0=8.0, Qinf=3.0, bIso=20.0, Hiso=5.0,
               C=[500.0, 300.0, 150.0], gam=[100.0, 50.0, 10.0]),
     [0.5, -0.2, -0.3, 0.4, 0.1, 0.2], MONO, None),
    ("kin_linear_cyclic",
     lr.Params(K, G, sig0=5.0, C=[120.0], gam=[0.0]),
     [1, 0, 0, 0, 0, 0], CYC, None),
    ("linkin_cf",
     lr.Params(K, G, sig0=5.0, Hiso=30.0, C=[200.0, 100.0], gam=[0.0, 0.0]),
     [0.5, -0.2, -0.3, 0.4, 0.1, 0.2], MONO, None),
    ("chaboche_nonprop",
     lr.Params(K, G, sig0=6.0, C=[500.0, 250.0], gam=[60.0, 12.0]),
     [0, 0, 0, 0, 0, 0], None, NONPROP),
]


def _run_cpp():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    src = os.path.join(_HERE, "ladrunoj2_kernel_check.cpp")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "lj2k.exe")
        cc = subprocess.run([gpp, "-O2", "-std=c++17", "-I", _SRC_ND, src, "-o", exe],
                            capture_output=True, text=True)
        assert cc.returncode == 0, f"g++ failed:\n{cc.stderr}"
        rr = subprocess.run([exe], capture_output=True, text=True)
        assert rr.returncode == 0, f"run failed:\n{rr.stderr}"
    # parse into {scenario: {step: {"STRESS":[...], "STATE":[...], "DTAN":[...]}}}
    out, cur = {}, None
    for line in rr.stdout.splitlines():
        p = line.split()
        if not p:
            continue
        if p[0] == "SCENARIO":
            cur = out.setdefault(p[1], {})
        elif p[0] == "STEP":
            step, kind = int(p[1]), p[2]
            cur.setdefault(step, {})[kind] = [float(x) for x in p[3:]]
    return out


def _eng_to_tensor(eng6):
    """Engineering-shear Voigt 6-vector → symmetric 3×3 (tensor shear = eng/2)."""
    t6 = [eng6[0], eng6[1], eng6[2], 0.5 * eng6[3], 0.5 * eng6[4], 0.5 * eng6[5]]
    return lr.tensor_from_voigt(t6)


def _step_tensorvecs(direction, mult, path):
    """Per-step strain as 6 tensor-comp vectors (path overrides direction*mult)."""
    if path is not None:
        return [np.array(p, float) for p in path]
    direction = np.array(direction, float)
    return [m * direction for m in mult]


def _oracle_path(params, direction, mult, path=None):
    """Step the oracle along the path; at each step record stress6, state, and a
    finite-difference tangent (vs engineering strain) from the pre-step state."""
    epsP = np.zeros((3, 3)); ebarP = 0.0
    alpha = [np.zeros((3, 3)) for _ in range(params.nBack)]
    steps = []
    for t6 in _step_tensorvecs(direction, mult, path):
        eps_t = lr.tensor_from_voigt(t6)                    # tensor comps
        res = lr.return_map(params, eps_t, epsP, ebarP, alpha)
        stress6 = lr.voigt_from_tensor(res["stress"])

        # FD tangent: d stress6 / d (engineering strain), about this step's
        # strain, from the SAME committed (pre-step) state.
        eng0 = np.array([t6[0], t6[1], t6[2],
                         2.0 * t6[3], 2.0 * t6[4], 2.0 * t6[5]])
        h = 1.0e-7
        Kfd = np.zeros((6, 6))
        for a in range(6):
            ep = eng0.copy(); ep[a] += h
            em = eng0.copy(); em[a] -= h
            sp = lr.voigt_from_tensor(
                lr.return_map(params, _eng_to_tensor(ep), epsP, ebarP, alpha)["stress"])
            sm = lr.voigt_from_tensor(
                lr.return_map(params, _eng_to_tensor(em), epsP, ebarP, alpha)["stress"])
            Kfd[:, a] = (sp - sm) / (2.0 * h)

        steps.append(dict(stress=stress6, ebarP=res["ebarP"], dGamma=res["dGamma"],
                          epsP=lr.voigt_from_tensor(res["epsP"]),
                          alpha=[lr.voigt_from_tensor(a) for a in res["alpha"]],
                          Kfd=Kfd, plastic=res["plastic"]))
        epsP, ebarP, alpha = res["epsP"], res["ebarP"], res["alpha"]
    return steps


def test_cpp_kernel_matches_oracle_and_fd_tangent():
    cpp = _run_cpp()
    assert set(cpp) == {s[0] for s in SCEN}, f"scenario mismatch: {set(cpp)}"

    saw_plastic = False
    for name, params, direction, mult, path in SCEN:
        ref = _oracle_path(params, direction, mult, path)
        for s, r in enumerate(ref):
            saw_plastic = saw_plastic or r["plastic"]
            cs = np.array(cpp[name][s]["STRESS"])
            cstate = cpp[name][s]["STATE"]
            cdt = np.array(cpp[name][s]["DTAN"]).reshape(6, 6)

            sscale = max(np.abs(r["stress"]).max(), 1.0)
            np.testing.assert_allclose(
                cs, r["stress"], rtol=1e-10, atol=1e-9 * sscale,
                err_msg=f"{name} step {s}: Cauchy stress mismatch")

            # STATE = [ebarP, dGamma, epsP(6), alpha_k(6)...]
            assert cstate[0] == pytest.approx(r["ebarP"], rel=1e-10, abs=1e-12)
            assert cstate[1] == pytest.approx(r["dGamma"], rel=1e-10, abs=1e-12)
            np.testing.assert_allclose(cstate[2:8], r["epsP"], rtol=1e-9, atol=1e-12,
                                       err_msg=f"{name} step {s}: epsP mismatch")
            off = 8
            for k, ak in enumerate(r["alpha"]):
                np.testing.assert_allclose(
                    cstate[off:off + 6], ak, rtol=1e-9, atol=1e-10,
                    err_msg=f"{name} step {s}: backstress {k} mismatch")
                off += 6

            # consistent tangent vs independent FD of the oracle stress
            tscale = max(np.abs(cdt).max(), 1.0)
            err = np.abs(cdt - r["Kfd"]).max()
            assert err <= 1e-4 * tscale, (
                f"{name} step {s}: Dtan != FD tangent (max abs err {err:.3e} "
                f"vs scale {tscale:.3e})")
    assert saw_plastic, "no scenario yielded — paths are too small"


# --------------------------------------------------------------------------- #
#  Independent CLOSED-FORM leg for the kinematic branch (does NOT use the       #
#  scalar-Newton oracle). For gamma=0 (linear kinematic) + linear isotropic     #
#  the return map collapses to a single closed-form step                        #
#      dGamma = f_tr / (2G + (2/3)(Hiso + sum C_k)),                            #
#  so this validates the (2/3)C theta term and the AF backstress update for the #
#  COMBINED-hardening path that reduce-to-J2Plasticity (isotropic-only) never   #
#  touches — catching a shared algebraic error the mirrored oracle could hide.  #
# --------------------------------------------------------------------------- #
_I3 = np.eye(3)
_R23 = np.sqrt(2.0 / 3.0)


def _closed_form_linkin(params, tensorvecs):
    """Analytic gamma=0 linear-kinematic + linear-iso path. Returns per-step
    (stress6, ebarP, dGamma, epsP6, [alpha6...]). No scalar Newton."""
    assert all(g == 0.0 for g in params.gam), "closed form requires gamma=0"
    assert params.Qinf == 0.0, "closed form requires linear iso (Qinf=0)"
    G, K = params.G, params.K
    hsum = params.Hiso + sum(params.C)
    epsP = np.zeros((3, 3)); ebarP = 0.0
    alpha = [np.zeros((3, 3)) for _ in range(params.nBack)]
    out = []
    for t6 in tensorvecs:
        eps = lr.tensor_from_voigt(t6)
        tr = np.trace(eps)
        edev = eps - (tr / 3.0) * _I3
        s_tr = 2.0 * G * (edev - epsP)
        M = s_tr - sum(alpha) if params.nBack else s_tr.copy()
        normM = np.sqrt(np.tensordot(M, M))
        f_tr = normM - _R23 * (params.sig0 + params.Hiso * ebarP)
        if f_tr <= 1e-12 * max(normM, _R23 * params.sig0, 1e-300):
            stress = K * tr * _I3 + s_tr
            dG = 0.0
        else:
            dG = f_tr / (2.0 * G + (2.0 / 3.0) * hsum)        # closed form
            n = M / normM
            epsP = epsP + dG * n
            alpha = [alpha[k] + (2.0 / 3.0) * params.C[k] * dG * n
                     for k in range(params.nBack)]
            ebarP = ebarP + _R23 * dG
            stress = K * tr * _I3 + s_tr - 2.0 * G * dG * n
        out.append((lr.voigt_from_tensor(stress), ebarP, dG,
                    lr.voigt_from_tensor(epsP),
                    [lr.voigt_from_tensor(a) for a in alpha]))
    return out


def test_cpp_kernel_linear_kinematic_closed_form():
    """The C++ kernel's combined-hardening (gamma=0) output must match the
    closed-form algebra — an oracle that shares NO code with the scalar Newton."""
    cpp = _run_cpp()
    name, params, direction, mult, path = next(s for s in SCEN if s[0] == "linkin_cf")
    ref = _closed_form_linkin(params, _step_tensorvecs(direction, mult, path))

    saw_plastic = False
    for s, (stress6, ebarP, dG, epsP6, alphas) in enumerate(ref):
        saw_plastic = saw_plastic or dG > 0.0
        cs = np.array(cpp[name][s]["STRESS"])
        cstate = cpp[name][s]["STATE"]
        sscale = max(np.abs(stress6).max(), 1.0)
        np.testing.assert_allclose(cs, stress6, rtol=1e-9, atol=1e-9 * sscale,
                                   err_msg=f"linkin_cf step {s}: stress vs closed form")
        assert cstate[0] == pytest.approx(ebarP, rel=1e-9, abs=1e-12)
        assert cstate[1] == pytest.approx(dG, rel=1e-9, abs=1e-12)
        np.testing.assert_allclose(cstate[2:8], epsP6, rtol=1e-9, atol=1e-12)
        off = 8
        for ak in alphas:
            np.testing.assert_allclose(cstate[off:off + 6], ak, rtol=1e-9, atol=1e-10)
            off += 6
    assert saw_plastic, "linkin_cf never yielded"
