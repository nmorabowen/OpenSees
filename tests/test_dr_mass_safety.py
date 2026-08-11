"""LadrunoDynamicRelaxation `-massSafety` — the Gershgorin stability margin (#728).

The raw Gershgorin fictitious mass `m_i = (dt²/4)Σ_j|K_ij|` bounds `ω_max·dt ≤ 2`
with EQUALITY: it marches EXACTLY ON the central-difference stability boundary,
where the amplification matrix is defective and round-off is amplified instead of
damped. A quasi-static DR run takes orders of magnitude more steps than a physical
explicit run, so it sits there long enough for the growth to become the answer —
and nothing reports an error while it happens (note 83 §3: a zero-push hold on an
EXACT equilibrium walked to an 87 kN residual on a 300 kN problem, silently).

`-massSafety f` divides the mass prefactor by `f²` so `ω_max·dt ≤ 2f`, at a cost of
`f²` relaxation per step. The default is 0.5, NOT 1.

The instrument here is the note's own control DR-0: push a block to a converged
static equilibrium with Newton, freeze the load, then relax with DR applying NO
further push. The state IS an equilibrium, so it must STAY there — every unit of
residual growth is manufactured by the integrator.

  MS-1  hold — H8/H20 × elastic/plastic stay put under the default margin.
  MS-2  NON-VACUITY — the same hold at the old `f = 1` blows up. Without this,
        MS-1 could pass for any reason at all (e.g. if DR never moved anything).
  MS-3  `ladrunoDR stabilityMargin` reads (ω_max·dt/2)² = f² on an unchanged tangent.
  MS-4  parser: forms, clamping, and the ignored-mode warning path.

See Ladruno_implementation/83_dynamic_relaxation_quadratic_probe.md §3 and the
LEDGER_quirks row it ships.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

DR = "LadrunoDynamicRelaxation"

_E, _NU = 2.0e5, 0.3
_K, _G = _E / (3.0 * (1 - 2 * _NU)), _E / (2.0 * (1 + _NU))
_SIG0 = 250.0                 # yield stress -> yield strain 1.25e-3
_L = 1.0                      # unit cube: sigma*A == E*|u| in the elastic range
_U_ELASTIC = -2.0e-4          # eps = 2e-4  (0.16 of yield)
_U_PLASTIC = -5.0e-3          # eps = 5e-3  (4x yield -> deeply plastic)

_DEFAULT_SAFETY = 0.5         # must track LADRUNO_DR_DEFAULT_MASS_SAFETY


# --------------------------------------------------------------------------
# model: a unit cube of LadrunoBrick (H8) or LadrunoBrick20 (H20), base fixed,
# top face driven down in z by an SP. 2x2x2 elements either way.
# --------------------------------------------------------------------------
def _build(kind, plastic, n=2):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    if plastic:
        ops.nDMaterial("LadrunoJ2", 1, _K, _G, "-iso", "voce", _SIG0, 0.0, 0.0, 200.0)
    else:
        ops.nDMaterial("ElasticIsotropic", 1, _E, _NU)

    # one half-index lattice serves both orders: for H20 a site is a node iff at
    # most one of its half-indices is odd (serendipity drops face/body centres).
    nn = n + 1 if kind == "H8" else 2 * n + 1
    h = _L / (nn - 1)
    nid, tag = {}, 1
    for k in range(nn):
        for j in range(nn):
            for i in range(nn):
                if kind == "H20" and (i % 2) + (j % 2) + (k % 2) > 1:
                    continue
                nid[(i, j, k)] = tag
                ops.node(tag, i * h, j * h, k * h)
                tag += 1
    top = nn - 1
    s = 1 if kind == "H8" else 2

    eid = 1
    for k in range(0, nn - 1, s):
        for j in range(0, nn - 1, s):
            for i in range(0, nn - 1, s):
                c = [nid[(i, j, k)], nid[(i + s, j, k)],
                     nid[(i + s, j + s, k)], nid[(i, j + s, k)],
                     nid[(i, j, k + s)], nid[(i + s, j, k + s)],
                     nid[(i + s, j + s, k + s)], nid[(i, j + s, k + s)]]
                if kind == "H8":
                    ops.element("LadrunoBrick", eid, *c, 1)
                else:
                    m = [nid[(i + 1, j, k)], nid[(i + s, j + 1, k)],
                         nid[(i + 1, j + s, k)], nid[(i, j + 1, k)],
                         nid[(i + 1, j, k + s)], nid[(i + s, j + 1, k + s)],
                         nid[(i + 1, j + s, k + s)], nid[(i, j + 1, k + s)],
                         nid[(i, j, k + 1)], nid[(i + s, j, k + 1)],
                         nid[(i + s, j + s, k + 1)], nid[(i, j + s, k + 1)]]
                    ops.element("LadrunoBrick20", eid, *(c + m), 1,
                                "-formulation", "uri")
                eid += 1

    tops = [t for (i, j, k), t in nid.items() if k == top]
    for (i, j, k), t in nid.items():
        if k == 0:
            ops.fix(t, 1, 1, 1)
    # a genuinely free interior/side witness node at mid height: its displacement
    # is the state, and a stable hold must not move it at all.
    kmid = max(1, top // 2)
    witness = next(t for (i, j, k), t in nid.items() if k == kmid)
    return tops, witness


def _push(tops, u_target, nsteps=10):
    """Newton to a CONVERGED static equilibrium at the imposed top settlement."""
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in tops:
        ops.sp(t, 3, u_target)
    ops.system("UmfPack")
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1e-10, 60, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    for _ in range(nsteps):
        assert ops.analyze(1) == 0, "static push failed to converge"
    ops.reactions()
    return sum(ops.nodeReaction(t, 3) for t in tops)


def _hold(safety, nsteps, chunk=500):
    """Relax with DR applying NO further push. Returns (ok, r_first, r_last, margin)."""
    ops.loadConst("-time", ops.getTime())
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("Diagonal")
    ops.test("NormUnbalance", 1e30, 1, 0)     # DR's gate is script-side, not here
    ops.algorithm("Linear")
    args = [DR] + (["-massSafety", safety] if safety is not None else [])
    ops.integrator(*args)
    ops.analysis("Transient")

    r_first = r_last = None
    done = 0
    while done < nsteps:
        if ops.analyze(chunk, 1.0) != 0:
            return False, r_first, r_last, ops.ladrunoDR("stabilityMargin")
        done += chunk
        r_last = ops.ladrunoDR("residualNorm")
        if r_first is None:
            r_first = r_last
    return True, r_first, r_last, ops.ladrunoDR("stabilityMargin")


def _case(kind, plastic, safety, nsteps=3000):
    tops, witness = _build(kind, plastic)
    u = _U_PLASTIC if plastic else _U_ELASTIC
    R0 = _push(tops, u)
    if plastic:
        # non-vacuity for the "plastic" label: the block must actually have yielded,
        # or this is just the elastic case with a different material name.
        assert abs(R0) < 0.9 * _E * abs(u) * _L * _L, (
            f"{kind} 'plastic' case did not yield: |R|={abs(R0):.4e} vs elastic "
            f"{_E * abs(u):.4e}"
        )
    w0 = [ops.nodeDisp(witness, d) for d in (1, 2, 3)]
    ok, r_first, r_last, margin = _hold(safety, nsteps)
    w1 = [ops.nodeDisp(witness, d) for d in (1, 2, 3)] if ok else [float("nan")] * 3
    ops.reactions()
    R1 = sum(ops.nodeReaction(t, 3) for t in tops) if ok else float("nan")
    return dict(ok=ok, r_first=r_first, r_last=r_last, margin=margin,
                R0=R0, R1=R1, w0=w0, w1=w1, u=u)


_GRID = [("H8", False), ("H8", True), ("H20", False), ("H20", True)]
_IDS = [f"{k}-{'plastic' if p else 'elastic'}" for k, p in _GRID]


# --------------------------------------------------------------------------
# MS-1 : the hold STAYS a hold under the default margin
# --------------------------------------------------------------------------
@pytest.mark.parametrize("kind,plastic", _GRID, ids=_IDS)
def test_ms1_zero_push_hold_stays_put(kind, plastic):
    c = _case(kind, plastic, None)          # None => the integrator's own default
    assert c["ok"], "DR aborted during a zero-push hold on an exact equilibrium"

    # the residual must stay at the level it started at, not walk. Scale it against
    # the problem's own load: this is the note's "87 kN on a 300 kN problem" made
    # into a gate.
    scale = abs(c["R0"])
    assert c["r_last"] < 1e-9 * scale, (
        f"{kind} residual walked to {c['r_last']:.4e} on a {scale:.4e} problem "
        f"(started {c['r_first']:.4e})"
    )
    # the STATE must not move either — a residual can stay small while the model
    # creeps, and the note's symptom was a drifting reaction, not a loud one.
    drift = max(abs(a - b) for a, b in zip(c["w0"], c["w1"]))
    assert drift < 1e-9 * abs(c["u"]), (
        f"{kind} witness node drifted {drift:.4e} during a zero-push hold"
    )
    assert abs(c["R1"] - c["R0"]) < 1e-9 * scale, (
        f"{kind} reaction drifted {c['R0']:.8e} -> {c['R1']:.8e} during a hold"
    )
    # ...and the margin the integrator reports must be the safe one it promised.
    assert c["margin"] == pytest.approx(_DEFAULT_SAFETY ** 2, rel=1e-6), (
        f"default stability margin {c['margin']} != {_DEFAULT_SAFETY}^2"
    )


# --------------------------------------------------------------------------
# MS-2 : NON-VACUITY — the old marginal factor is detectably worse
# --------------------------------------------------------------------------
@pytest.mark.parametrize("kind,plastic", _GRID, ids=_IDS)
def test_ms2_boundary_factor_is_detectably_worse(kind, plastic):
    """`-massSafety 1` is the pre-#728 behaviour: the bare stability boundary.

    MS-1 is only worth anything if this fails, so assert the failure directly.
    """
    c = _case(kind, plastic, 1.0)
    blew_up = (not c["ok"]) or (c["r_last"] > 1e6 * max(c["r_first"], 1e-300)) \
        or not math.isfinite(c["r_last"])
    assert blew_up, (
        f"{kind} at -massSafety 1 did NOT diverge (r {c['r_first']:.4e} -> "
        f"{c['r_last']:.4e}); MS-1 would then be passing vacuously"
    )
    # and it must be worse in the PHYSICAL quantity too, not just the residual
    if c["ok"]:
        assert not (abs(c["R1"] - c["R0"]) < 1e-9 * abs(c["R0"])), (
            "reaction held steady at the stability boundary — re-derive the gate"
        )
    # the boundary is exactly where the integrator says it is
    assert c["margin"] >= 1.0 - 1e-9, f"margin at f=1 read {c['margin']}, expected >= 1"


# --------------------------------------------------------------------------
# MS-3 : the reported margin IS (omega_max*dt/2)^2 = f^2 on an unchanged tangent
# --------------------------------------------------------------------------
@pytest.mark.parametrize("f", [1.0, 0.75, 0.5, 0.25])
def test_ms3_margin_equals_safety_squared(f):
    # elastic, so the tangent never changes and the margin must be exactly f^2:
    # any other reading means the prefactor and the diagnostic disagree.
    tops, _ = _build("H8", False)
    _push(tops, _U_ELASTIC)
    ok, _, _, margin = _hold(f, 100, chunk=50)
    assert ok
    assert margin == pytest.approx(f * f, rel=1e-9), (
        f"stabilityMargin {margin} != f^2 = {f * f} for -massSafety {f}"
    )


# --------------------------------------------------------------------------
# MS-4 : parser
# --------------------------------------------------------------------------
@pytest.mark.parametrize("args,expect", [
    ((DR,), _DEFAULT_SAFETY ** 2),                       # default is NOT 1
    ((DR, "-massSafety", 0.25), 0.0625),
    ((DR, "-massSafety", 1.0), 1.0),                     # opt back in (warns)
    ((DR, "-massSafety", 4.0), 1.0),                     # clamped to the limit
    ((DR, "-massSafety", 0.0), _DEFAULT_SAFETY ** 2),    # non-positive -> default
    ((DR, "-massSafety", -1.0), _DEFAULT_SAFETY ** 2),
    ((DR, "-mass", "gershgorin", "-massSafety", 0.5), 0.25),
    ((DR, "-massSafety", 0.5, "-damping", "viscous"), 0.25),
])
def test_ms4_parser(args, expect):
    tops, _ = _build("H8", False)
    _push(tops, _U_ELASTIC)
    ops.loadConst("-time", ops.getTime())
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("Diagonal")
    ops.test("NormUnbalance", 1e30, 1, 0)
    ops.algorithm("Linear")
    ops.integrator(*args)
    ops.analysis("Transient")
    assert ops.analyze(20, 1.0) == 0, f"first steps failed for {args}"
    assert ops.ladrunoDR("stabilityMargin") == pytest.approx(expect, rel=1e-9)


def test_ms4b_massSafety_ignored_for_non_gershgorin_mass():
    """-massSafety only scales the gershgorin mass; unity/lumped have no such bound.

    The integrator warns and reports margin = -1 ("not available") rather than
    quoting a number it cannot stand behind.
    """
    tops, _ = _build("H8", False)
    _push(tops, _U_ELASTIC)
    ops.loadConst("-time", ops.getTime())
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("Diagonal")
    ops.test("NormUnbalance", 1e30, 1, 0)
    ops.algorithm("Linear")
    ops.integrator(DR, "-mass", "unity", "-massSafety", 0.25)
    ops.analysis("Transient")
    assert ops.analyze(5, 1.0) == 0
    assert ops.ladrunoDR("stabilityMargin") == -1.0
