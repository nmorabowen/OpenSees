# 22 — Conformal RC member recipe (v1a) + adaptive driver (v1b)

**Status:** shipped — `Ladruno_scripts/ladruno_rc.py` (builder) + `Ladruno_scripts/ladruno_solve.py`
(driver), both self-testing, pure openseespy, **no `SRC/` change**. This is items 1–2 of the
[[20_ladruno_embedded_reinforcement_adr]] §8 v1 build order, sitting on the passed
[[21_rc3d_validation_gates]] (Gate 1 conformal viable, Gate 2 confinement real).

## What v1 gives you today

For **regular, grid-aligned** RC members you can build a full 3D solid model with discrete
rebar **right now, with no new C++** — the embedded element (ADR 20, v2) is only for the
non-grid / dense-hoop cases conformal can't mesh.

- **`ladruno_rc.build_rc_column(ops, RCColumnSpec(...))`** — rectangular RC column: an
  `ASDConcrete3D` `stdBrick` solid (Gate-2-calibrated backbone → emergent confinement),
  longitudinal `corotTruss` bars on every perimeter vertical node-line (**shared nodes ⇒
  perfect bond, no bond element**), and optional perimeter tie rings every `tie_every`
  layers. Returns a `BuiltRC` with the node grid + base/top sets so you own BCs/loads.
- **`ladruno_rc.pushover(ops, built, axial_ratio=…, drift_ratio=…)`** — fixes the base,
  applies constant axial, runs a lateral displacement-control pushover through the adaptive
  driver. Returns `(AdaptiveResult, V_peak, drift)`.
- **`ladruno_solve.adaptive_static(ops, set_step, n_steps, …)`** — the blessed driver:
  algorithm ladder (KrylovNewton → Newton → NewtonLineSearch → ModifiedNewton),
  halve-on-fail, grow-on-success, `on_step` history hook. Integrator-agnostic via the
  `set_step(scale)` callback. **Use this for any RC softening run** — plain Newton
  fixed-step diverges past peak (Gate evidence).

## Minimal example
```python
from ladruno_rc import RCColumnSpec, build_rc_column, pushover
ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
b = build_rc_column(ops, RCColumnSpec(bx=0.4, by=0.4, H=2.0, nx=2, ny=2, nz=5,
                                      rho_long=0.01, ties=True))
res, Vpeak, drift = pushover(ops, b, axial_ratio=0.1, drift_ratio=0.02)
```

## Smoke result (2026-06-03, build f6700775)
```
built: 20 bricks, 8 long bars x 5, 48 tie segs, 162 dof
pushover: <AdaptiveResult OK 400.000/400 substeps=433 min_scale=0.125 algos={'KrylovNewton':405}>
V_peak = 148.6 kN at drift 0.0200   -> SMOKE PASS
```
The driver took 433 substeps for 400 nominal (halving to 1/8 at the hard points) — the ties
make this stiffer than the bare Gate-1 column, so the adaptive halving is doing real work.

## Scope & next
- **Covers:** rectangular columns with grid-aligned longitudinal bars + perimeter ties,
  monotonic pushover, implicit. Confinement is emergent (validate `Kc`+backbone per
  [[LEDGER_quirks]]).
- **Not yet:** beams/walls (same builder pattern, add next), explicit conformal run,
  apeGmsh front-end (build the model from a `Part`), cyclic protocol helper, bond-slip
  (that needs the v2 embedded element + `LadrunoBondSlip`).
- **Reproduce:** `pythoncore-3.12-64/python.exe Ladruno_scripts/ladruno_rc.py` (and
  `…/ladruno_solve.py`) → SMOKE / SELF-TEST PASS.
