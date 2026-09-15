"""WP F8 (ADR-94 addendum) -- the ASD Drucker-Prager APEX REGION under
ASSOCIATED flow.

WHAT IS WRONG
-------------
wp/94f (#832) fixed the ZERO-DILATANCY misclassification by adding an
ELASTIC-METRIC apex test inside ``Backward_Euler`` and **unioning** it with the
yield function's own EUCLIDEAN one::

    bool be_in_apex = yf.check_apex_region(...);          // (p - p_apex) >= eta*q
    if (!be_in_apex)
        be_in_apex = cp_apex_region(...);                 // (p - p_apex) >= (K*etabar/G)*q

The union is correct at ``etabar = 0``, where the elastic-metric region
(``p >= p_apex``) strictly CONTAINS the Euclidean one, so the wider of the two
is the exact one.  It is exactly WRONG under **associated** flow, where
``K*etabar/G`` is the LARGER slope -- on the ADR-95 Prandtl deck
``K/G = 9.667`` and ``eta = 0.4457``, so the exact apex region is the cone of
half-slope ``4.308`` and the Euclidean test is a cone of half-slope ``0.4457``,
nearly TEN TIMES wider.  Taking the union therefore keeps the *too wide*
Euclidean answer and apex-projects every trial state in the wedge::

    eta * q  <=  p - p_apex  <  (K * etabar / G) * q          (the F8 WEDGE)

whose correct return is to the cone FLANK.  The committed stress is then pinned
at ``sigma_apex = p_apex * I`` -- essentially zero -- instead of the flank
state, and under ``tangent_type Continuum`` the Gauss point additionally reports
a ZERO tangent.  On the Prandtl footing that is a spurious stress collapse at
the footing-edge Gauss points, i.e. the associated leg's wall.

THE ORACLE
----------
Associated, perfectly plastic Drucker-Prager is elementary in closed form, and
needs no numpy: with ``f = q + eta*p - xi_c``, ``q = sqrt(J2)`` and an isotropic
elastic operator, the cone return from an elastic trial ``(p_tr, q_tr)`` is

    dgamma = f_tr / (G + K*eta*etabar)
    q_ret  = q_tr - G * dgamma
    p_ret  = p_tr - K * etabar * dgamma

and it is the correct (closest-point-in-the-elastic-metric) answer exactly while
``q_ret >= 0`` -- which is the same statement as ``p_tr - p_apex < (K*etabar/G)*q_tr``.
So the wedge state below has an exact answer to compare against, and "apex" is
not it.

TRAPS OBEYED (ADR-94 Sec. 8, as in tests/test_adr94f_asd_apex_fallback.py)
``LadrunoBrick`` (``stdBrick`` swallows material return codes); a fully
sp-prescribed rig so the Gauss-point strain is the prescribed one; every
measurement in a FRESH SUBPROCESS, because the ``INT_OPT_*`` option maps are
keyed by material tag and shared process-wide and because this .pyd's ``cout``
is invisible to ``capfd``.

WALL TIME: ~15 s (four subprocesses, one element each).
"""
import json
import math
import os
import subprocess
import sys
import textwrap

import pytest

_DIST = os.environ.get(
    "LADRUNO_DIST_BIN",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir,
                 "dist", "bin"))
if os.path.isdir(_DIST):
    sys.path.insert(0, _DIST)

from _testbed import ops  # noqa: E402

pytestmark = [pytest.mark.zone_a]


# the ADR-95 Prandtl deck's own cone (phi_txc = 20 deg, SY = 0.2), mapped in
# asd_path_diag.py:  eta = 3*alpha, xi_c = SY/sqrt(3)
E, NU = 45000.0, 0.45
ETA, XI_C = 0.445749, 0.115470
K_EL = E / (3.0 * (1.0 - 2.0 * NU))          # 150000.0
G_EL = E / (2.0 * (1.0 + NU))                # 15517.2413793...
P_APEX = XI_C / ETA                          # 0.2590...
SLOPE_ASSOC = K_EL * ETA / G_EL              # 4.3083... the EXACT apex slope

IV_DP = ("BackStress(TensorLinearHardeningFunction):"
         "DP_cohesion(ScalarLinearHardeningFunction):")


def _normal_strains(p_tr, q_tr):
    """Prescribed (exx, eyy, ezz) whose ELASTIC trial from a zero stress state
    has mean stress ``p_tr`` and ``sqrt(J2) = q_tr`` (see wp/94f's twin)."""
    s_sum = p_tr / K_EL                             # 2a + b = tr(eps)
    d = -q_tr * math.sqrt(3.0) / (2.0 * G_EL)       # a - b  (b > a: tension on z)
    a = (d + s_sum) / 3.0
    b = s_sum - 2.0 * a
    return a, a, b


def _cone_return(p_tr, q_tr, etabar):
    """The closed-form Drucker-Prager cone return (see the module docstring)."""
    f_tr = q_tr + ETA * p_tr - XI_C
    dgamma = f_tr / (G_EL + K_EL * ETA * etabar)
    return p_tr - K_EL * etabar * dgamma, q_tr - G_EL * dgamma


_CHILD = r'''
import json, sys
sys.path.insert(0, r"{dist}")
import opensees as ops

ops.wipe()
ops.model("basic", "-ndm", 3, "-ndf", 3)
cube = {{1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
        5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}}
fixs = {{1: (1, 1, 1), 2: (0, 1, 1), 3: (0, 0, 1), 4: (1, 0, 1),
        5: (1, 1, 0), 6: (0, 1, 0), 8: (1, 0, 0)}}
for t, c in cube.items():
    ops.node(t, float(c[0]), float(c[1]), float(c[2]))
for t, m in fixs.items():
    ops.fix(t, m[0], m[1], m[2])

ops.nDMaterial("ASDPlasticMaterial3D", 1,
    "DruckerPrager_YF", "DruckerPrager_PF", "LinearIsotropic3D_EL",
    "{iv}",
    "Begin_Model_Parameters",
    "YoungsModulus", {E}, "PoissonsRatio", {NU},
    "DP_xi_c", {xi_c}, "DP_eta", {eta}, "DP_etabar", {etabar},
    "TensorLinearHardeningParameter", 0.0,
    "ScalarLinearHardeningParameter", 0.0,
    "MassDensity", 0.0,
    "End_Model_Parameters",
    "Begin_Internal_Variables",
    "BackStress", 0., 0., 0., 0., 0., 0.,
    "DP_cohesion", 0.0,
    "End_Internal_Variables",
    "Begin_Integration_Options",
    "tangent_type", "{tangent}",
    "strict_convergence", {strict},
    "End_Integration_Options")

ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
for n in (2, 3, 6, 7):
    ops.sp(n, 1, {ex})
for n in (3, 4, 7, 8):
    ops.sp(n, 2, {ey})
for n in (5, 6, 7, 8):
    ops.sp(n, 3, {ez})
ops.constraints("Transformation")
ops.numberer("Plain")
ops.system("UmfPack")
ops.test("NormDispIncr", 1.0e-10, 50, 0)
ops.algorithm("Newton")
ops.integrator("LoadControl", 1.0 / {nsteps})
ops.analysis("Static")

codes, hist = [], []
for _ in range({nsteps}):
    rc = ops.analyze(1)
    codes.append(rc)
    if rc != 0:
        break
    hist.append(list(ops.eleResponse(1, "stresses"))[0:6])

out = {{"build": str(ops.ladrunoBuild()), "codes": codes, "hist": hist}}
sys.stdout.write("RESULT " + json.dumps(out) + "\n")
'''


def _run(etabar, ex, ey, ez, nsteps=1, strict=1, tangent="Continuum"):
    script = _CHILD.format(dist=os.path.abspath(_DIST), iv=IV_DP, E=E, NU=NU,
                           xi_c=XI_C, eta=ETA, etabar=etabar, strict=strict,
                           tangent=tangent, ex=repr(ex), ey=repr(ey),
                           ez=repr(ez), nsteps=nsteps)
    p = subprocess.run([sys.executable, "-c", textwrap.dedent(script)],
                       capture_output=True, text=True,
                       stdin=subprocess.DEVNULL, timeout=600)
    raw = p.stdout + p.stderr
    line = [ln for ln in raw.splitlines() if ln.startswith("RESULT ")]
    assert line, f"child produced no RESULT line.\n--- output ---\n{raw[-4000:]}"
    res = json.loads(line[-1][len("RESULT "):])
    res["raw"] = raw
    return res


def _pq(sig):
    p = (sig[0] + sig[1] + sig[2]) / 3.0
    sx, sy, sz = sig[0] - p, sig[1] - p, sig[2] - p
    j2 = 0.5 * (sx * sx + sy * sy + sz * sz) + sig[3] ** 2 + sig[4] ** 2 + sig[5] ** 2
    return p, math.sqrt(max(j2, 0.0))


def _f(sig):
    p, q = _pq(sig)
    return q + ETA * p - XI_C


@pytest.fixture(scope="module")
def dp_available():
    try:
        print("\nladrunoBuild() =", ops.ladrunoBuild())
    except Exception:                                   # pragma: no cover
        pytest.skip("this build has no ladrunoBuild() -- wrong/old binary")
    a, aa, b = _normal_strains(0.5 * P_APEX, 0.0)
    res = _run(0.0, a, aa, b, nsteps=2)
    if res["codes"][0] != 0:
        pytest.skip("ASDPlasticMaterial3D / DruckerPrager_YF / LadrunoBrick "
                    "not usable in this build")


# ===========================================================================
# 1. THE F8 WEDGE -- Euclidean says APEX, the elastic metric says CONE
# ===========================================================================
def test_f8_associated_wedge_returns_to_the_cone_not_the_apex(dp_available):
    """A trial state beyond the apex whose correct associated return is to the
    cone FLANK, at a mean stress well BELOW the apex.

    Pre-fix this commits ``sigma = p_apex * I`` (q = 0) because
    ``Backward_Euler`` unions the yield function's Euclidean answer with the
    elastic-metric one and the Euclidean cone is ~10x too wide here.  The
    closed-form answer is checked to 1e-6 relative -- there is no tolerance to
    tune: the two candidate answers differ by the whole deviator.
    """
    q_tr = 1.0
    p_tr = P_APEX + 2.0                        # squarely inside the wedge
    assert (p_tr - P_APEX) >= ETA * q_tr, "state is not Euclidean-apex"
    assert (p_tr - P_APEX) < SLOPE_ASSOC * q_tr, "state is not elastic-metric cone"

    p_ex, q_ex = _cone_return(p_tr, q_tr, ETA)
    assert q_ex > 0.0, "oracle says apex — the wedge was mis-constructed"

    ex, ey, ez = _normal_strains(p_tr, q_tr)
    res = _run(ETA, ex, ey, ez, nsteps=1, strict=1)
    assert res["codes"] == [0], (
        f"the wedge step failed: codes={res['codes']}\n{res['raw'][-3000:]}")

    last = res["hist"][-1]
    p, q = _pq(last)
    assert abs(_f(last)) <= 1e-6 * XI_C, f"committed state is off the surface: {last}"
    assert q > 0.5 * q_ex, (
        f"the associated wedge state was APEX-PROJECTED: committed q = {q!r} "
        f"(p = {p!r}) against the closed-form cone return q = {q_ex!r}, "
        f"p = {p_ex!r}.  Euclidean apex slope {ETA:.4f} vs exact "
        f"{SLOPE_ASSOC:.4f}; this state sits at (p - p_apex)/q = "
        f"{(p_tr - P_APEX) / q_tr:.4f}.")
    assert abs(q - q_ex) <= 1e-6 * q_ex, (
        f"cone return off the closed form: q={q!r} vs {q_ex!r}")
    assert abs(p - p_ex) <= 1e-6 * abs(p_ex), (
        f"cone return off the closed form: p={p!r} vs {p_ex!r}")


# ===========================================================================
# 2. the TRUE associated apex region is still returned to the apex
# ===========================================================================
def test_f8_associated_true_apex_region_still_projects(dp_available):
    """Beyond the EXACT elastic-metric slope the apex is the right answer, and
    must stay the right answer: this is the half of the classification the fix
    must not break."""
    q_tr = 1.0
    p_tr = P_APEX + 2.0 * SLOPE_ASSOC * q_tr   # far inside the true apex region
    p_ex, q_ex = _cone_return(p_tr, q_tr, ETA)
    assert q_ex < 0.0, "oracle says cone — the case was mis-constructed"

    ex, ey, ez = _normal_strains(p_tr, q_tr)
    res = _run(ETA, ex, ey, ez, nsteps=1, strict=1)
    assert res["codes"] == [0], (
        f"the true-apex step failed: codes={res['codes']}\n{res['raw'][-3000:]}")
    last = res["hist"][-1]
    p, q = _pq(last)
    assert q <= 1e-8 * XI_C, f"the apex return left a deviator: q={q!r}"
    assert abs(p - P_APEX) <= 1e-8 * P_APEX, (
        f"committed p={p!r} is not the apex pressure {P_APEX!r}")


# ===========================================================================
# 3. the zero-dilatancy classification wp/94f fixed is UNCHANGED
# ===========================================================================
def test_f8_zero_dilatancy_wedge_still_projects_to_the_apex(dp_available):
    """wp/94f's own case, restated: at ``etabar = 0`` the exact region is
    ``p >= p_apex`` and the Euclidean test is the narrow one, so a state the
    Euclidean test calls CONE must still be apex-projected.  A fix that made
    the elastic-metric test REPLACE the Euclidean one rather than widen it must
    leave this untouched."""
    p_tr = 2.0 * P_APEX
    q_tr = 3.0 * (p_tr - P_APEX) / ETA
    assert (p_tr - P_APEX) < ETA * q_tr, "test does not exercise wp/94f's wedge"
    ex, ey, ez = _normal_strains(p_tr, q_tr)
    res = _run(0.0, ex, ey, ez, nsteps=10, strict=1)

    assert res["codes"] == [0] * 10, (
        f"wp/94f's zero-dilatancy wedge regressed: codes={res['codes']}\n"
        f"{res['raw'][-3000:]}")
    assert "scalar Newton exhausted" not in res["raw"], res["raw"][-3000:]
    last = res["hist"][-1]
    p, q = _pq(last)
    assert q <= 1e-8 * XI_C, f"the apex return left a deviator: q={q!r}"
    assert abs(p - P_APEX) <= 1e-8 * P_APEX, f"committed p={p!r} != {P_APEX!r}"


# ===========================================================================
# 4. an ordinary associated flank step is untouched
# ===========================================================================
def test_f8_associated_cone_flank_far_from_the_apex(dp_available):
    """Deviatoric-dominated COMPRESSIVE associated step: nowhere near the apex
    under either classification, pinned against the same closed form."""
    p_tr = -20.0 * P_APEX
    q_tr = 3.0 * (XI_C - ETA * p_tr)
    p_ex, q_ex = _cone_return(p_tr, q_tr, ETA)
    ex, ey, ez = _normal_strains(p_tr, q_tr)
    res = _run(ETA, ex, ey, ez, nsteps=1, strict=1)
    assert res["codes"] == [0], (
        f"the associated flank control failed: codes={res['codes']}\n"
        f"{res['raw'][-3000:]}")
    last = res["hist"][-1]
    p, q = _pq(last)
    assert abs(q - q_ex) <= 1e-6 * q_ex, f"q={q!r} vs closed form {q_ex!r}"
    assert abs(p - p_ex) <= 1e-6 * abs(p_ex), f"p={p!r} vs closed form {p_ex!r}"
