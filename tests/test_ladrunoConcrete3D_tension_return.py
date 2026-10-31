"""LadrunoConcrete3D: the DIRECT tension return near first cracking (r4 rescue of fork PR #877; review 3, MAJOR 1).

Review 3 of #877 (2026-10-10) found the direct backward-Euler hardening return failing at random once kappa_p passes
about 2: in 2000-point windows of 2e-6 relative width in virgin near-uniaxial tension (lateral -0.2, the wrapper
defaults below) it failed 645 / 904 / 1953 times at x = 1.5 / 3 / 6 eps_cr (re-measured on the PR head f19d55bd1 with
the kernel probe), and one window at lateral +0.3, x = 3 returned two different admissible roots (vertex kappa_p 23,
1984 times; regular kappa_p 113, 16 times). Cause: the plain Newton fails there and the globalized Newton's
unprojected line search walks out of the admissible set (dlam < 0, kappa_p < 0). The r4 port of db02ccc70 tries the
vertex first (a deterministic vertex-vs-regular selection) and rescues a failed plain Newton with a bracketed nested
solve whose root is accepted only on the FULL residual (LadrunoConcrete3DKernel.h returnMapHardening /
returnMapNested / polishHardeningRoot; kernel side pinned by concrete3d_kernel_check.cpp section D).

These cases pin it through the built material under `-subIncr refuse`, the mode that has NO sub-increment rescue, so a
single failed direct return is a refused step: every point of each window must integrate, with one kappa_p and one
stress per window. Every Gauss point of a LadrunoBrick under an affine displacement field sees the same strain, the
trial is imposed in one `algorithm Linear` step from the virgin state (Penalty, Constant series), as in the #2617
cases of test_ladrunoConcrete3D_element.py. Measured wall time: ~3 s.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
_E, _NU, _FC, _FT, _GF, _GC, _LCH = 30000.0, 0.2, 30.0, 3.0, 0.1, 5.0, 50.0
_ECR = _FT / _E
_N = 24            # points per 2e-6-wide window (the review used 2000 on the kernel; the failure rate was 30-97 %)


def _imposed(eps, mode="refuse"):
    """One LadrunoBrick, the review's wrapper defaults (bilinear, cdpm2 flow + drive, tcTemper proj, Df 0.85,
    hardening 0.3 0.01, lch 50, eps_fc = Gc/(fc lch)), under the homogeneous normal strain `eps` = (exx, eyy, ezz)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *c)
    ops.nDMaterial("LadrunoConcrete3D", 1, _E, _NU, _FC, _FT, _GF, _GC, "-lch", _LCH, "-gcLegacy", "-subIncr", mode)
    ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, c in _CUBE.items():
        for i in range(3):
            ops.sp(t, i + 1, eps[i] * c[i])
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Penalty", 1.0e12, 1.0e12)
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    return ops.analyze(1)


def _window(lat, x):
    out = []
    for i in range(_N):
        sc = x * _ECR * (1.0 + 2.0e-6 * i / _N)
        rc = _imposed((sc, lat * sc, lat * sc))
        kp = float(ops.eleResponse(1, "material", 1, "kappaP")[0]) if rc == 0 else float("nan")
        sxx = float(ops.eleResponse(1, "material", 1, "stress")[0]) if rc == 0 else float("nan")
        out.append((rc, kp, sxx))
    return out


@pytest.mark.parametrize("x", [1.5, 3.0, 6.0])
def test_direct_tension_return_never_fails_in_a_window(x):
    """Lateral -0.2 (near-uniaxial tension, past first cracking): the direct return integrates every point of the
    2e-6 window (before r4: 32 / 45 / 98 % of a 2000-point window failed at x = 1.5 / 3 / 6), and the window has one
    root (kappa_p and sigma_xx vary only with the 2e-6 strain change)."""
    w = _window(-0.2, x)
    refused = [i for i, (rc, _, _) in enumerate(w) if rc != 0]
    assert not refused, f"x = {x}: {len(refused)} of {_N} direct returns refused (indices {refused[:8]})"
    kps = [kp for _, kp, _ in w]
    sxx = [s for _, _, s in w]
    assert kps[0] > 2.0, f"precondition: the window must be past kappa_p ~ 2 (kappa_p = {kps[0]})"
    assert max(kps) - min(kps) <= 1.0e-4 * max(kps), f"x = {x}: kappa_p spread {min(kps)}..{max(kps)} in one window"
    assert max(sxx) - min(sxx) <= 1.0e-5 * _FT, f"x = {x}: sigma_xx spread {min(sxx)}..{max(sxx)} in one window"


def test_vertex_versus_regular_root_is_deterministic():
    """Lateral +0.3, x = 3: before r4 the window returned the vertex root (kappa_p 23) or the regular one (kappa_p 113)
    depending on rounding. The vertex is now tried first and taken whenever admissible, so the whole window returns the
    vertex root, and returns it every time."""
    w = _window(0.3, 3.0)
    assert all(rc == 0 for rc, _, _ in w), "a direct return was refused in the lateral +0.3 window"
    kps = [kp for _, kp, _ in w]
    assert max(kps) - min(kps) <= 1.0e-4 * max(kps), f"two roots in one window: kappa_p {min(kps)}..{max(kps)}"
    assert 20.0 < kps[0] < 26.0, f"expected the vertex root (kappa_p ~ 23.1), got {kps[0]}"


def test_refuse_and_deterministic_agree_on_a_direct_return():
    """Where the direct return integrates and the deterministic overshoot count is 1, `-subIncr refuse` and the
    default `deterministic` are the same map (both are the direct return): identical stress at a first-crack
    point with a small overshoot."""
    eps = (1.05 * _ECR, -0.2 * 1.05 * _ECR, -0.2 * 1.05 * _ECR)
    assert _imposed(eps, "refuse") == 0
    s_ref = list(ops.eleResponse(1, "material", 1, "stress"))
    assert _imposed(eps, "deterministic") == 0
    s_det = list(ops.eleResponse(1, "material", 1, "stress"))
    assert s_ref == pytest.approx(s_det, rel=1e-12, abs=1e-12)
