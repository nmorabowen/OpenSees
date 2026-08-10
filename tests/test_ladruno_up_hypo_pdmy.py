"""LadrunoUP ``-geom hypo`` + PressureDependMultiYield — the ADR-79 P3 smoke.

The driving application of the whole hypo lane (ADR 79 §1): rate-form
pressure-dependent soil plasticity at genuine large strain on the saturated
u-p element — the combination `-geom finite` can never serve (the LogStrain
adaptor requires a linear-elastic inner; PDMY's G ~ p'^0.5 violates it), and
the reason the unrotated-frame design was chosen (PDMY's nested-surface
tensorial state is private and can never be co-rotated by an element).

A TIMs-style footing-in-miniature (3x3x3 saturated PDMY box) through the full
staging idiom — elastic gravity (stage 0) -> updateMaterialStage 1 -> plastic
re-equilibration -> DISPLACEMENT-CONTROLLED footing push to 5% penetration —
run under ``-geom hypo`` with ``-geom corot`` (the ADR-78-validated method on
exactly this regime) as the anchor twin.

Two scoping findings from the probe runs, both load-bearing for the design:

* The push must be DISPLACEMENT-controlled. A force-controlled footing push
  on stage-1 PDMY diverges from the first increment: the near-surface GPs sit
  at ~zero confinement, G ~ p'^d makes their tangent nearly singular, and the
  first Newton iterate is unbounded — under -geom hypo the huge iterate
  inverts the top-layer hexes and updateHypo (correctly) cuts the step;
  under -geom linear the iterate is representable and Newton simply grinds
  without converging. Displacement control bounds the iterates and is how a
  bearing backbone is traced anyway.

* ``-geom linear`` is NOT a usable twin on this problem (it grinds, see
  above) — corot is the anchor, and the hypo-vs-corot gap IS the measured
  geometric-nonlinearity content: probe values 0.3% at 0.25% penetration,
  1.8% at 1%, 15.4% at 5% (hypo stiffer: UL assembly on the compacted
  current configuration + n(J) storage). The gates pin exactly that shape —
  tight agreement early, bounded divergence late.

The full bearing-limit-point campaign on the real SFIM mesh is analysis
work, not CI.

Run:  py -3.12 -m pytest tests/test_ladruno_up_hypo_pdmy.py -x -q
"""
import os
import sys
from pathlib import Path

import numpy as np
import pytest

_DIST = str(Path(__file__).resolve().parents[1] / "dist" / "bin")
if not os.path.isfile(os.path.join(_DIST, "opensees.pyd")):
    pytest.skip(f"worktree engine not built: {_DIST}", allow_module_level=True)

from _engine import bind_worktree_engine  # noqa: E402
ops = bind_worktree_engine(_DIST)

pytestmark = [pytest.mark.zone_b]

_G = 9.81
_RHO_SAT = 2.0
_RHO_W = 1.0
_KF = 2.2e6
_PORO = 0.4

# medium-loose sand, UCSD modern 16-arg signature, nd=3 (the 2D column
# battery's parameter set — test_ladruno_up_element_pdmy.py — lifted to 3D;
# liquefaction OFF: monotonic push, not shaking)
PDMY = dict(rho=_RHO_SAT, G=5.5e4, B=1.5e5, phi=33.0, gammaPeak=0.1,
            refPress=80.0, d=0.5, PT=27.0, contract=0.05, dil1=0.0, dil2=0.0,
            liq1=0.0, liq2=0.0, liq3=0.0, NYS=20)

_N = 3          # 3x3x3 box of H8 (compact; ~20 s per leg measured)
_L = 3.0
_SMAX = 0.05    # footing penetration: 5% of the box height
_NSTEP = 20
_DT = 25.0


def _mk_pdmy(tag):
    p = PDMY
    ops.nDMaterial("PressureDependMultiYield", tag, 3, p["rho"], p["G"],
                   p["B"], p["phi"], p["gammaPeak"], p["refPress"], p["d"],
                   p["PT"], p["contract"], p["dil1"], p["dil2"], p["liq1"],
                   p["liq2"], p["liq3"], p["NYS"])


def _box_mesh(n, l):
    tags = {}
    t = 1
    for k in range(n + 1):
        for j in range(n + 1):
            for i in range(n + 1):
                ops.node(t, l * i / n, l * j / n, l * k / n)
                tags[(i, j, k)] = t
                t += 1
    conns = []
    for k in range(n):
        for j in range(n):
            for i in range(n):
                conns.append([tags[(i, j, k)], tags[(i + 1, j, k)],
                              tags[(i + 1, j + 1, k)], tags[(i, j + 1, k)],
                              tags[(i, j, k + 1)], tags[(i + 1, j, k + 1)],
                              tags[(i + 1, j + 1, k + 1)], tags[(i, j + 1, k + 1)]])
    return tags, conns


def _footing_run(geom):
    """Gravity(stage 0) -> flip -> settle -> displacement-controlled footing
    push. Returns (reaction backbone, hypoState of the center-column element
    or None)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 4)
    tags, conns = _box_mesh(_N, _L)
    _mk_pdmy(1)
    for e, conn in enumerate(conns, start=1):
        ops.element("LadrunoUP", e, *conn, 1,
                    "-Kf", _KF, "-poro", _PORO, "-rhoF", _RHO_W,
                    "-perm", 1e-4, 1e-4, 1e-4, "-stab", "auto", 0.25,
                    "-body", 0.0, 0.0, -_G, "-geom", geom)
    for (i, j, k), t in tags.items():
        fx = 1 if (i == 0 or i == _N) else 0
        fy = 1 if (j == 0 or j == _N) else 0
        fz = 1 if k == 0 else 0
        if k == 0:
            fx = fy = 1
        fp = 1 if k == _N else 0            # drained top
        ops.fix(t, fx, fy, fz, fp)

    # stage 0: elastic gravity (big-dt Newmark transient, the upstream idiom)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-7, 100, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")
    for _ in range(10):
        assert ops.analyze(1, 500.0) == 0, f"{geom}: elastic gravity failed"
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    for _ in range(10):
        assert ops.analyze(1, 500.0) == 0, f"{geom}: plastic settle failed"

    # displacement-controlled central footing push (2x2-node patch)
    ops.loadConst("-time", 0.0)
    c = _N // 2
    foot = [tags[(c + di, c + dj, _N)] for di in (0, 1) for dj in (0, 1)]
    times = [_DT * (s + 1) for s in range(_NSTEP)] + [_DT * (_NSTEP + 5)]
    vals = [-_SMAX * (s + 1) / _NSTEP for s in range(_NSTEP)] + [-_SMAX]
    ops.timeSeries("Path", 2, "-time", *times, "-values", *vals, "-useLast")
    ops.pattern("Plain", 2, 2)
    for t in foot:
        ops.sp(t, 3, 1.0)
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-8, 60, 0)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")
    backbone = []
    for s in range(_NSTEP):
        ok = ops.analyze(1, _DT)
        if ok != 0:                      # substep fallback (pdmy-battery policy)
            ops.test("NormDispIncr", 1e-7, 100, 0)
            ops.algorithm("KrylovNewton")
            ok = 0
            for _ in range(10):
                if ops.analyze(1, _DT / 10) != 0:
                    ok = -1
                    break
            ops.test("NormDispIncr", 1e-8, 60, 0)
            ops.algorithm("Newton")
        assert ok == 0, f"{geom}: push step {s + 1} hard-failed"
        ops.reactions()
        backbone.append(sum(ops.nodeReaction(t, 3) for t in foot))

    hypo_state = None
    if geom == "hypo":
        hs = []
        for e in range(1, len(conns) + 1):
            hs.append(np.array(ops.eleResponse(e, "hypoState"),
                               dtype=float).reshape(-1, 3))
        hypo_state = np.vstack(hs)
    for e in range(1, len(conns) + 1):
        s = np.array(ops.eleResponse(e, "stresses"), dtype=float)
        assert np.all(np.isfinite(s)), f"{geom}: non-finite stress in ele {e}"
    return np.array(backbone), hypo_state


def test_pdmy_footing_smoke_hypo_vs_corot():
    bb_hyp, hs = _footing_run("hypo")
    bb_cor, _ = _footing_run("corot")

    # backbones: finite, compressive, monotonically growing bearing reaction
    for name, bb in (("hypo", bb_hyp), ("corot", bb_cor)):
        assert np.all(np.isfinite(bb)), f"{name} backbone not finite"
        assert bb[-1] < 0.0, f"{name} reaction has the wrong sign"
        assert np.all(np.diff(np.abs(bb)) > 0.0), f"{name} backbone not monotone"

    # hypo kinematic state responded and stayed physical
    J, n = hs[:, 0], hs[:, 1]
    assert np.all(J > 0.5) and np.all(J < 1.5), (
        f"unphysical J range: {J.min()}..{J.max()}")
    assert np.all(n > 0.0) and np.all(n < 1.0)
    assert J.min() < 1.0 - 1e-4, "no compaction anywhere — the push did nothing"
    n_want = 1.0 - (1.0 - _PORO) / J
    assert np.abs(n - n_want).max() <= 1e-9, "porosity-J identity broke under PDMY"

    # the measured geometric-nonlinearity shape (module docstring): tight
    # agreement with corot at small penetration, bounded smooth divergence at
    # 5% — hypo STIFFER (UL on the compacted configuration).
    rel_early = np.abs(bb_hyp[:4] - bb_cor[:4]).max() / np.abs(bb_cor[:4]).max()
    rel_final = abs(bb_hyp[-1] - bb_cor[-1]) / abs(bb_cor[-1])
    print(f"PDMY footing: R_final hypo={bb_hyp[-1]:.4e} corot={bb_cor[-1]:.4e} "
          f"rel_early={rel_early:.3e} rel_final={rel_final:.3e}; "
          f"J range [{J.min():.4f}, {J.max():.4f}]")
    assert rel_early <= 0.05, (
        f"hypo diverges from corot already at small penetration ({rel_early:.1%})")
    assert rel_final <= 0.30, (
        f"hypo vs corot final bearing reaction differs by {rel_final:.1%} — "
        "beyond the documented geometric-content band")
    assert abs(bb_hyp[-1]) > abs(bb_cor[-1]), (
        "hypo backbone not stiffer than corot — the measured UL-compaction "
        "signature inverted")
