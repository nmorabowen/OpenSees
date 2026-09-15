"""WP F8 (ADR-94 addendum) -- the ASD Drucker-Prager APEX REGION under
DILATANT flow (`etabar > eta*G/K`), which starts at about psi = 2.3 deg.

WHAT IS WRONG
-------------
wp/94f (#832) fixed the ZERO-DILATANCY misclassification by adding an
ELASTIC-METRIC apex test inside ``Backward_Euler`` and **unioning** it with the
yield function's own EUCLIDEAN one::

    bool be_in_apex = yf.check_apex_region(...);          // (p - p_apex) >= eta*q
    if (!be_in_apex)
        be_in_apex = cp_apex_region(...);                 // (p - p_apex) >= (K*etabar/G)*q

A union keeps the WIDER of the two regions, so it is correct exactly while the
EXACT one is the wider, i.e. while::

    K * etabar / G  <=  eta      <=>      etabar  <=  eta * G / K

On the ADR-95 Prandtl deck ``G/K = 0.10345`` and ``eta = 0.4457``, so the
crossover is ``etabar = 0.0461`` -- **psi ~ 2.3 deg**.  The union is therefore
unsafe not "under associated flow" but under essentially ANY dilatancy;
``etabar = 0`` is the single case where it is safe, and it is the single case
wp/94f measured:

    etabar = 0        exact slope 0       (p >= p_apex)   exact region is WIDER -> safe
    etabar = eta*G/K  exact slope 0.4457                  the crossover
    etabar = eta/2    exact slope 2.1545                  Euclidean 4.8x too wide
    etabar = eta      exact slope 4.3089                  Euclidean ~10x too wide

Above the crossover the union keeps the *too wide* Euclidean answer and
apex-projects every trial state in the wedge::

    eta * q  <=  p - p_apex  <  (K * etabar / G) * q          (the F8 WEDGE)

whose correct return is to the cone FLANK.  The committed stress is then pinned
at ``sigma_apex = p_apex * I`` -- essentially zero -- instead of the flank
state, and under ``tangent_type Continuum`` the Gauss point additionally reports
a ZERO tangent.  On the Prandtl footing that is a spurious stress collapse at
the footing-edge Gauss points, i.e. the dilatant leg's wall.

AND WHAT NARROWING THE REGION EXPOSED (review round 1)
-------------------------------------------------------
Routing those states to the flank Newton instead is only right if the flank
Newton answers them correctly.  With cohesion SOFTENING it does not: its
``dPhi/dlambda`` carries the shipped yield function's ``df/dk = -1`` term for an
internal variable that ``f`` does not contain (ADR-97 P0 header finding 2,
pinned not fixed), and it CONVERGES -- ``|Phi| ~ 1e-7``, ``rc = 0``, no
exhaustion, no refusal -- onto a state whose deviator points OPPOSITE the trial
deviator.  A Drucker-Prager return is a non-negative radial scaling of the trial
deviator plus a pressure change, so that is inadmissible for ANY parameters, and
no yield-function tolerance can see it (the state is ON the surface).  Cases 6
and 7 below pin the geometric guard that catches it.

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

WALL TIME: MEASURED 0.71 s for the whole file, 12 tests (one-element
subprocesses, ~0.05 s each).  zone_a, NOT slow tier.
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
    "TensorLinearHardeningParameter", {ht},
    "ScalarLinearHardeningParameter", {hs},
    "MassDensity", 0.0,
    "End_Model_Parameters",
    "Begin_Internal_Variables",
    "BackStress", {a0}, {a1}, {a2}, {a3}, {a4}, {a5},
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


def _run(etabar, ex, ey, ez, nsteps=1, strict=1, tangent="Continuum", hs=0.0,
         ht=0.0, alpha0=(0.0,) * 6):
    a = [repr(float(v)) for v in alpha0]
    script = _CHILD.format(dist=os.path.abspath(_DIST), iv=IV_DP, E=E, NU=NU,
                           xi_c=XI_C, eta=ETA, etabar=etabar, strict=strict,
                           tangent=tangent, ex=repr(ex), ey=repr(ey),
                           ez=repr(ez), nsteps=nsteps, hs=repr(float(hs)),
                           ht=repr(float(ht)), a0=a[0], a1=a[1], a2=a[2],
                           a3=a[3], a4=a[4], a5=a[5])
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


# ===========================================================================
# 5. the condition is NOT "associated" — it is  etabar > eta*G/K
# ===========================================================================
# Review round 1, finding (1).  The union is unsafe whenever the EXACT slope
# exceeds the Euclidean one, `K*etabar/G > eta`, i.e. `etabar > eta*G/K`.  On
# this deck `G/K = 0.10345`, so the threshold is `etabar = 0.046`, which is
# `psi ~ 2.3 deg` — every ordinary dilatant sand, not only `psi = phi`.  At
# `etabar = eta/2` (`psi ~ phi/2`) the exact slope is 2.1545 against the
# Euclidean 0.4457: a 4.8x-wide wedge that was apex-projected pre-fix.
ETABAR_HALF = 0.5 * ETA
SLOPE_HALF = K_EL * ETABAR_HALF / G_EL          # 2.15446...
ETABAR_UNSAFE_THRESHOLD = ETA * G_EL / K_EL     # 0.046108... (psi ~ 2.3 deg)


def test_f8_half_dilatancy_wedge_returns_to_the_cone(dp_available):
    """`etabar = eta/2`: the Euclidean test says APEX, the elastic metric says
    CONE, and the closed form says which is right."""
    assert ETABAR_HALF > ETABAR_UNSAFE_THRESHOLD, (
        "eta/2 is below the unsafe threshold on this cone — the case is void")
    q_tr = 1.0
    p_tr = P_APEX + 1.0                         # 0.4457 < ratio 1.0 < 2.1545
    assert (p_tr - P_APEX) >= ETA * q_tr, "state is not Euclidean-apex"
    assert (p_tr - P_APEX) < SLOPE_HALF * q_tr, "state is not elastic-metric cone"

    p_ex, q_ex = _cone_return(p_tr, q_tr, ETABAR_HALF)
    assert q_ex > 0.0, "oracle says apex — the wedge was mis-constructed"

    ex, ey, ez = _normal_strains(p_tr, q_tr)
    res = _run(ETABAR_HALF, ex, ey, ez, nsteps=1, strict=1)
    assert res["codes"] == [0], (
        f"the eta/2 wedge step failed: codes={res['codes']}\n{res['raw'][-3000:]}")
    last = res["hist"][-1]
    p, q = _pq(last)
    assert q > 0.5 * q_ex, (
        f"a psi ~ phi/2 state was APEX-PROJECTED: committed q = {q!r} (p = {p!r}) "
        f"against the closed-form cone return q = {q_ex!r}, p = {p_ex!r}. "
        f"Euclidean slope {ETA:.4f}, exact slope at etabar = eta/2 "
        f"{SLOPE_HALF:.4f}; this state sits at (p - p_apex)/q = "
        f"{(p_tr - P_APEX) / q_tr:.4f}.")
    assert abs(q - q_ex) <= 1e-6 * q_ex, f"q={q!r} vs closed form {q_ex!r}"
    assert abs(p - p_ex) <= 1e-6 * abs(p_ex), f"p={p!r} vs closed form {p_ex!r}"


def test_f8_half_dilatancy_true_apex_region_still_projects(dp_available):
    """Just beyond the exact slope at `etabar = eta/2` the apex is the answer."""
    q_tr = 1.0
    p_tr = P_APEX + 2.1760 * q_tr               # the reviewer's apex row
    assert (p_tr - P_APEX) > SLOPE_HALF * q_tr, "row is not inside the apex region"
    p_ex, q_ex = _cone_return(p_tr, q_tr, ETABAR_HALF)
    assert q_ex < 0.0, "oracle says cone — the row was mis-constructed"

    ex, ey, ez = _normal_strains(p_tr, q_tr)
    res = _run(ETABAR_HALF, ex, ey, ez, nsteps=1, strict=1)
    assert res["codes"] == [0], (
        f"the eta/2 apex step failed: codes={res['codes']}\n{res['raw'][-3000:]}")
    last = res["hist"][-1]
    p, q = _pq(last)
    assert q <= 1e-8 * XI_C, f"the apex return left a deviator: q={q!r}"
    assert abs(p - P_APEX) <= 1e-8 * P_APEX, f"committed p={p!r} != {P_APEX!r}"


# ===========================================================================
# 6. COHESION SOFTENING — the narrowed classification must not hand the flank
#    Newton a state it answers with a SIGN-FLIPPED deviator
# ===========================================================================
# Review round 1, finding (3).  Narrowing the apex region routes near-boundary
# trials into the flank scalar Newton, whose `dPhi/dlambda` carries the pinned
# vanilla `df/dk = -1` term for a cohesion internal variable that `f` itself
# does not contain (ADR-97 P0 header finding 2).  With SOFTENING (HS < 0) that
# Newton can converge to `|f| ~ 1e-7` on a state whose deviator points OPPOSITE
# the trial deviator — it walked through the vertex and out the other side.  A
# Drucker-Prager return is a non-negative radial scaling of the trial deviator
# plus a pressure change, so a sign flip is inadmissible for ANY parameters, and
# no yield-function tolerance can see it.  Layer (b) cannot catch it either:
# `be_exhausted` is only computed under `strict_convergence`, which is OFF by
# default, and this state does not exhaust anyway — it CONVERGES, to the wrong
# root.
HS_SOFTENING = -20000.0
HS_HARDENING = 2000.0


def _run_softening(ratio, strict, hs=HS_SOFTENING, etabar=None):
    q_tr = 1.0
    p_tr = P_APEX + ratio * q_tr
    ex, ey, ez = _normal_strains(p_tr, q_tr)
    res = _run(ETA if etabar is None else etabar, ex, ey, ez,
               nsteps=1, strict=strict, hs=hs)
    return p_tr, q_tr, res


@pytest.mark.parametrize("ratio", [3.0, SLOPE_ASSOC])
@pytest.mark.parametrize("strict", [0, 1])
def test_f8_softening_never_commits_a_flipped_deviator(dp_available, ratio, strict):
    """The admissibility the return map owes regardless of hardening law."""
    p_tr, q_tr, res = _run_softening(ratio, strict)
    if res["codes"] != [0]:
        return                      # a refusal is an acceptable outcome here
    last = res["hist"][-1]
    p, q = _pq(last)
    # the trial deviator has sigma_zz TENSILE relative to the mean (see
    # _normal_strains), so an admissible return keeps s_zz - s_xx >= 0
    s_gap_tr = 1.0                  # by construction, positive
    s_gap = last[2] - last[0]
    assert s_gap * s_gap_tr >= -1e-12 * max(abs(q), XI_C), (
        f"the return committed a SIGN-FLIPPED deviator at ratio {ratio:.4f}, "
        f"strict={strict}, HS={HS_SOFTENING}: committed (p, q) = ({p!r}, {q!r}) "
        f"with s_zz - s_xx = {s_gap!r} against a trial deviator of the opposite "
        f"sign. A Drucker-Prager return is a NON-NEGATIVE radial scaling of the "
        f"trial deviator plus a pressure change; this state walked through the "
        f"vertex. |f| = {_f(last)!r}")
    assert _f(last) <= 1e-6 * XI_C, f"committed state is outside f: {_f(last)!r}"


# The HS = 0 and HS = +2000 rows, measured on the PRE-guard build `3324485f7`.
# These are the CORRECT answers (ratio 3.0 is a cone return at etabar = eta,
# whose exact slope is 4.3089; ratio 4.3089 is the apex) and they are pinned so
# the softening guard is proved INERT on every non-softening law.  The HS = 0
# row is pinned to 1e-15 relative (perfectly plastic, exact); the HS = +2000 row
# to 3e-6, which is the size of the `df/dk = -1` term's own effect on that row
# (3e-8 in f), not a fitted tolerance.
_NONSOFT_ROWS = {
    (0.0, 3.0): (-0.18910264022674624, 0.19976231277843354, 1e-14),
    (0.0, None): (0.25904713190607503, 0.0, 1e-14),
    (2000.0, 3.0): (-0.18910259602908464, 0.19976232303571295, 3e-6),
    (2000.0, None): (0.2590471871367433, 0.0, 3e-6),
}


@pytest.mark.parametrize("hs", [0.0, HS_HARDENING])
def test_f8_nonsoftening_rows_are_unaffected(dp_available, hs):
    """The companion control: the guard must not touch a perfectly plastic or
    HARDENING cohesion law.  Both rows are pinned to values measured BEFORE the
    guard existed."""
    for ratio, key in ((3.0, 3.0), (SLOPE_ASSOC, None)):
        p_ref, q_ref, rtol = _NONSOFT_ROWS[(hs, key)]
        p_tr, q_tr, res = _run_softening(ratio, strict=1, hs=hs)
        assert res["codes"] == [0], (
            f"HS={hs}, ratio={ratio}: step failed {res['codes']} "
            f"{res['raw'][-2000:]}")
        last = res["hist"][-1]
        p, q = _pq(last)
        assert abs(p - p_ref) <= rtol * max(abs(p_ref), XI_C), (
            f"HS={hs}, ratio={ratio}: p moved, {p!r} vs pre-guard {p_ref!r}")
        assert abs(q - q_ref) <= rtol * max(abs(q_ref), XI_C), (
            f"HS={hs}, ratio={ratio}: q moved, {q!r} vs pre-guard {q_ref!r}")
        assert last[2] - last[0] > 0.0 or q <= 1e-6 * XI_C, (
            f"HS={hs}, ratio={ratio}: deviator sign flipped: {last!r}")


# ===========================================================================
# 8. NONZERO BACK STRESS -- the apex region and the flip guard both live in the
#    RELATIVE deviator r = dev(sigma) - alpha
# ===========================================================================
# Review round 2.  ``DruckerPrager_YF::check_apex_region`` has always measured in
# ``r`` (its own line 146); ``cp_apex_region``, the elastic-metric twin the
# integrator uses, measured in the raw ``dev(sigma)``.  With ``alpha``
# antiparallel to the trial deviator the two disagree: the raw deviator can cross
# zero while ``sqrt(J2(r))`` stays positive, so a cone state is classified APEX.
# ``be_apex_project`` then refuses the step, because the vertex ``apex_stress()``
# names is not on the surface either -- ``|f(sigma_apex)| = sqrt(J2(alpha))``,
# measured 0.034641 for the alpha below.  Same variable, same mistake, in the
# deviator-flip guard.
#
# The reviewer's two trials, with the oracle computed HERE rather than
# transcribed: the elastic trial from a zero stress state built by
# ``_normal_strains`` has deviator ``(-k, -k, 2k)`` with ``k = q_tr/sqrt(3)``, so
# ``r_tr`` and the closed-form cone return are both elementary.
ALPHA0 = (0.02, 0.02, -0.04, 0.0, 0.0, 0.0)
BACKSTRESS_TRIALS = [(0.05, 0.5103), (0.02, 0.3810)]
# the reviewer's quoted answer, which is trial 2's (trial 1 lands 5.0e-6 away in
# sqrt(J2(r)) -- the two trials do NOT share a return point)
REVIEWER_QREL, REVIEWER_P = 0.017321, 0.220190


def _q_rel(sig, alpha):
    """sqrt(J2) of the RELATIVE deviator r = dev(sigma) - alpha."""
    p = (sig[0] + sig[1] + sig[2]) / 3.0
    r = [sig[0] - p - alpha[0], sig[1] - p - alpha[1], sig[2] - p - alpha[2],
         sig[3] - alpha[3], sig[4] - alpha[4], sig[5] - alpha[5]]
    j2 = 0.5 * (r[0] ** 2 + r[1] ** 2 + r[2] ** 2) + r[3] ** 2 + r[4] ** 2 + r[5] ** 2
    return math.sqrt(max(j2, 0.0))


def _relative_trial(p_tr, q_tr, alpha):
    """sqrt(J2(r)) of the elastic trial `_normal_strains(p_tr, q_tr)` produces."""
    k = q_tr / math.sqrt(3.0)
    dev = (-k, -k, 2.0 * k, 0.0, 0.0, 0.0)
    r = [dev[i] - alpha[i] for i in range(6)]
    j2 = 0.5 * (r[0] ** 2 + r[1] ** 2 + r[2] ** 2) + r[3] ** 2 + r[4] ** 2 + r[5] ** 2
    return math.sqrt(max(j2, 0.0))


@pytest.mark.parametrize("q_tr,p_tr", BACKSTRESS_TRIALS)
@pytest.mark.parametrize("strict", [0, 1])
def test_f8_backstress_admissible_return_is_not_refused(dp_available, q_tr, p_tr,
                                                        strict):
    """A legitimate cone return whose ABSOLUTE deviator has crossed zero.  It
    must be accepted, and it must land on the closed form."""
    q_rel_tr = _relative_trial(p_tr, q_tr, ALPHA0)
    p_ex, q_ex = _cone_return(p_tr, q_rel_tr, ETA)       # the same closed form,
    assert q_ex > 0.0, "oracle says apex -- the trial was mis-constructed"

    ex, ey, ez = _normal_strains(p_tr, q_tr)
    res = _run(ETA, ex, ey, ez, nsteps=1, strict=strict, ht=0.0, alpha0=ALPHA0)
    assert res["codes"] == [0], (
        f"an ADMISSIBLE return was refused at (q_tr, p_tr) = ({q_tr}, {p_tr}), "
        f"strict={strict}, alpha0={ALPHA0}: codes={res['codes']}. Its exact "
        f"return is sqrt(J2(r)) = {q_ex!r}, p = {p_ex!r}. Both the apex "
        f"classification and the flip guard must measure in the RELATIVE "
        f"deviator; in the raw one this state reads as 'impossible' because "
        f"dev(sigma) crosses zero when alpha is antiparallel to it.\n"
        f"{res['raw'][-2500:]}")
    last = res["hist"][-1]
    p = (last[0] + last[1] + last[2]) / 3.0
    q_rel = _q_rel(last, ALPHA0)
    assert abs(q_rel - q_ex) <= 1e-6 * q_ex, (
        f"sqrt(J2(r)) = {q_rel!r} vs closed form {q_ex!r}")
    assert abs(p - p_ex) <= 1e-6 * abs(p_ex), (
        f"p = {p!r} vs closed form {p_ex!r}")


def test_f8_backstress_oracle_matches_the_reviewers_quoted_answer():
    """Provenance, not physics: the closed form used above reproduces the figure
    the review quoted, so the two are the same calculation.  It is trial 2's --
    trial 1 returns to a point 5.0e-6 away in sqrt(J2(r))."""
    q_tr, p_tr = BACKSTRESS_TRIALS[1]
    p_ex, q_ex = _cone_return(p_tr, _relative_trial(p_tr, q_tr, ALPHA0), ETA)
    assert abs(q_ex - REVIEWER_QREL) <= 5e-7, (q_ex, REVIEWER_QREL)
    assert abs(p_ex - REVIEWER_P) <= 5e-7, (p_ex, REVIEWER_P)
