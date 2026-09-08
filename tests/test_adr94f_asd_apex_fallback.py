"""ADR-94 wp/94f -- the FIX gate for the ASD Drucker-Prager APEX REGION TEST.

WHAT WAS WRONG
--------------
``DruckerPrager_YF::CHECK_APEX_REGION`` (shipped by #815 / wp/94c) is the
EUCLIDEAN normal-cone test in the (p, sqrt(J2)) half-plane::

    apex  <=>  (p - p_apex) >= eta * q

and says so in its own comment.  The exact test lives in the ELASTIC metric::

    apex  <=>  (p - p_apex) >= (K * etabar / G) * q

At ``etabar = 0`` -- non-associated, zero dilatancy, which is the ADR-95
Prandtl-Reissner footing deck -- the exact test degenerates to ``p >= p_apex``,
because a non-dilatant flank return CANNOT MOVE p AT ALL.  Every over-apex trial
state with ``q > eta*(p - p_apex)`` was therefore routed to the flank map, which
has no admissible solution (``f = q + eta*p - xi_c >= eta*p - xi_c > 0`` for
every ``q >= 0``), so its scalar Newton exhausted and, under
``strict_convergence``, the step was REFUSED.  Measured on the ADR-95 deck: 435
refusals and a FLOOR at s/B 0.01122.

THE FIX (both layers live in the integrator; the YF signature cannot see K, G
or etabar)
* **(a) elastic-metric pre-check.**  ``Backward_Euler`` unions the YF's Euclidean
  answer with ``cp_apex_region`` -- ADR-97's family-agnostic elastic-metric
  classification, reused verbatim -- for yield functions that opt in through the
  new ``yf_apex_elastic_metric`` trait (today: Drucker-Prager only).
* **(b) flank-first apex fallback.**  At exactly the sites that would otherwise
  return ``LADRUNO_MATERIAL_REFUSED`` (Newton exhaustion under
  ``strict_convergence``, singular local tangent, NaN, strict plastic
  inconsistency), a trial whose MEAN STRESS is beyond the apex takes the apex
  projection instead, under the same ``|f(sigma_apex)| <= tol_yf`` guard.

WHAT IS MEASURED HERE
---------------------
1. the misclassified wedge at ``etabar = 0`` is ACCEPTED and lands ON the apex,
   with no "scalar Newton exhausted" refusal anywhere in the child's output;
2. the same path at ``etabar = eta`` (associated) commits an ADMISSIBLE state --
   the two classifications may legitimately differ there, the answer may not;
3. the pure hydrostatic-tension case (#815's own case) is unchanged;
4. a cone-flank step far from the apex is untouched, pinned against the closed
   form of the non-dilatant Drucker-Prager return.

Byte-identity of ``Backward_Euler`` on everything that is NOT an apex/failure
state is gated separately and far more strongly by
``tests/test_adr97_p4_inertness.py`` (23 decks, baseline taken at 3622d6214,
i.e. long before this change).

TRAPS OBEYED (ADR-94 Sec. 8)
----------------------------
``LadrunoBrick`` (``stdBrick`` swallows material return codes); ``UmfPack``
(non-associated => unsymmetric tangent); fully sp-prescribed rig; every
measurement in a FRESH SUBPROCESS, because the ``INT_OPT_*`` option maps are
keyed by material tag and shared process-wide, and because this .pyd's ``cout``
is invisible to ``capfd``.
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


# ---------------------------------------------------------------------------
# the ADR-95 Prandtl deck's own Drucker-Prager, mapped in asd_path_diag.py
# (phi_txc = 20 deg, SY = 0.2):  eta = 3*alpha, xi_c = SY/sqrt(3)
# ---------------------------------------------------------------------------
E, NU = 45000.0, 0.45
ETA, XI_C = 0.445749, 0.115470
K_EL = E / (3.0 * (1.0 - 2.0 * NU))          # 150000.0
G_EL = E / (2.0 * (1.0 + NU))                # 15517.2413793...
P_APEX = XI_C / ETA                          # 0.2590...

IV_DP = ("BackStress(TensorLinearHardeningFunction):"
         "DP_cohesion(ScalarLinearHardeningFunction):")


def _normal_strains(p_tr, q_tr):
    """Prescribed (exx, eyy, ezz) whose ELASTIC trial from a zero stress state
    has mean stress ``p_tr`` and ``sqrt(J2) = q_tr``.

    With exx = eyy = a, ezz = b:  p = K*tr(eps) = K*(2a+b) and, writing
    d = a - b, the deviatoric stress is s = (2G/3)*(d, d, -2d), so
    sqrt(J2) = sqrt(0.5*(s:s)) = (2/sqrt(3)) * G * |d|.
    """
    s_sum = p_tr / K_EL                             # 2a + b = tr(eps)
    d = -q_tr * math.sqrt(3.0) / (2.0 * G_EL)       # a - b  (b > a: tension on z)
    a = (d + s_sum) / 3.0
    b = s_sum - 2.0 * a
    return a, a, b


# ---------------------------------------------------------------------------
# child process: build, drive, and report.  Prints ONE json line prefixed with
# "RESULT " so the (noisy) material diagnostics can be scanned separately.
# ---------------------------------------------------------------------------
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


def _run(etabar, ex, ey, ez, nsteps=10, strict=1):
    script = _CHILD.format(dist=os.path.abspath(_DIST), iv=IV_DP, E=E, NU=NU,
                           xi_c=XI_C, eta=ETA, etabar=etabar, strict=strict,
                           ex=repr(ex), ey=repr(ey), ez=repr(ez),
                           nsteps=nsteps)
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
    """Also the place the build hash is reported (ladrunoBuild provenance)."""
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
# 1. the misclassified wedge at etabar = 0 -- THE ADR-95 FAILURE
# ===========================================================================
def test_94f_zero_dilatancy_over_apex_wedge_is_accepted(dp_available):
    """Trial state beyond the apex with SMALL BUT NONZERO shear, chosen so the
    EUCLIDEAN test says CONE (``q > eta*(p - p_apex)``) while the exact
    elastic-metric test at ``etabar = 0`` says APEX (``p > p_apex``).

    With zero dilatancy the flank map cannot move p, so no flank return exists:
    ``f = q + eta*p - xi_c >= eta*p_tr - xi_c > 0`` for every q >= 0.  The apex
    is the ONLY admissible answer, and pre-94f the step was refused instead.
    """
    p_tr = 1.5 * P_APEX
    q_tr = 3.0 * ETA * (p_tr - P_APEX)          # squarely inside the wedge
    assert q_tr > ETA * (p_tr - P_APEX), "test does not exercise the wedge"

    ex, ey, ez = _normal_strains(p_tr, q_tr)
    res = _run(0.0, ex, ey, ez, nsteps=10, strict=1)

    assert res["codes"] == [0] * 10, (
        f"the over-apex wedge path still fails: codes={res['codes']}\n"
        f"{res['raw'][-3000:]}")
    assert "scalar Newton exhausted" not in res["raw"], (
        "the flank map is still being asked to solve an over-apex state:\n"
        + res["raw"][-3000:])
    assert "rejecting step" not in res["raw"], (
        "a refusal was issued on a path that completed:\n" + res["raw"][-3000:])

    last = res["hist"][-1]
    assert all(v == v and abs(v) < 1e30 for v in last), f"non-finite: {last}"
    p, q = _pq(last)
    assert abs(p - P_APEX) <= 1e-8 * P_APEX, (
        f"committed mean stress {p!r} is not the apex pressure {P_APEX!r}")
    assert q <= 1e-8 * XI_C, f"the apex return left deviatoric stress q={q!r}"
    assert abs(_f(last)) <= 1e-6 * XI_C, (
        f"committed state is inadmissible: f={_f(last)!r}")


# ===========================================================================
# 2. the same path, ASSOCIATED (etabar = eta) -- classifications may differ,
#    the committed state may not be inadmissible
# ===========================================================================
def test_94f_associated_dilatancy_still_commits_an_admissible_state(dp_available):
    p_tr = 1.5 * P_APEX
    q_tr = 3.0 * ETA * (p_tr - P_APEX)
    ex, ey, ez = _normal_strains(p_tr, q_tr)
    res = _run(ETA, ex, ey, ez, nsteps=10, strict=1)

    assert res["codes"] == [0] * 10, (
        f"the associated over-apex path fails: codes={res['codes']}\n"
        f"{res['raw'][-3000:]}")
    last = res["hist"][-1]
    assert all(v == v and abs(v) < 1e30 for v in last), f"non-finite: {last}"
    assert _f(last) <= 1e-6 * XI_C, (
        f"committed state is OUTSIDE the yield surface: f={_f(last)!r}")


# ===========================================================================
# 3. #815's own case -- pure hydrostatic tension -- unchanged
# ===========================================================================
def test_94f_pure_hydrostatic_tension_is_unchanged(dp_available):
    """Zero deviator the whole way: the Euclidean and the elastic-metric tests
    AGREE here (both say apex), so this path must behave exactly as #815 left
    it -- driven to 2x the apex volumetric strain, it commits sigma = p_apex*I.
    """
    ev = 2.0 * P_APEX / (3.0 * K_EL)        # p = K*tr(eps) = 2*p_apex
    res = _run(0.0, ev, ev, ev, nsteps=10, strict=1)

    assert res["codes"] == [0] * 10, (
        f"hydrostatic tension regressed: codes={res['codes']}\n"
        f"{res['raw'][-3000:]}")
    last = res["hist"][-1]
    p, q = _pq(last)
    assert q <= 1e-8 * XI_C, f"hydrostatic path grew a deviator: {last}"
    assert abs(p - P_APEX) <= 1e-8 * P_APEX, (
        f"committed p={p!r} is not xi_c/eta = {P_APEX!r}")


# ===========================================================================
# 4. a cone-flank state FAR from the apex -- the fix must not touch it
# ===========================================================================
def test_94f_cone_flank_far_from_the_apex_is_untouched(dp_available):
    """Deviatoric-dominated COMPRESSIVE path: p_tr = -20*p_apex, so the apex is
    nowhere near, `check_apex_region` is false and `cp_apex_region` is false,
    and the scalar Newton converges as it always did.

    The non-dilatant Drucker-Prager flank return is available in closed form
    (etabar = 0 => the return is a pure radial scaling of the deviator at
    FROZEN p), which is what this pins:  q_ret = xi_c - eta*p_tr, p_ret = p_tr.
    Byte-identity of the whole committed history is gated separately by
    tests/test_adr97_p4_inertness.py (baseline 3622d6214, pre-94f).
    """
    p_tr = -20.0 * P_APEX
    q_tr = 3.0 * (XI_C - ETA * p_tr)            # 3x the flank radius: plastic
    ex, ey, ez = _normal_strains(p_tr, q_tr)
    res = _run(0.0, ex, ey, ez, nsteps=10, strict=1)

    assert res["codes"] == [0] * 10, (
        f"the cone-flank control path fails: codes={res['codes']}\n"
        f"{res['raw'][-3000:]}")
    last = res["hist"][-1]
    p, q = _pq(last)

    # p is FROZEN at the elastic trial (zero dilatancy) ...
    assert abs(p - p_tr) <= 1e-8 * abs(p_tr), (
        f"a non-dilatant flank return moved the mean stress: {p!r} vs {p_tr!r}")
    # ... and the deviator sits exactly on the cone.
    q_exact = XI_C - ETA * p_tr
    assert abs(q - q_exact) <= 1e-6 * q_exact, (
        f"flank return is off the cone: q={q!r} vs closed form {q_exact!r}")
    assert abs(_f(last)) <= 1e-6 * XI_C, (
        f"committed flank state is inadmissible: f={_f(last)!r}")
