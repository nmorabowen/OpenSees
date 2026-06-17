---
title: Robust nonlinear-solution driver — user guide
project: Ladruno
status: shipped
tags:
  - guide
  - solver
  - robust-solve
  - integrator
---

# Robust nonlinear-solution driver — user guide

How to drive a hard nonlinear OpenSees analysis to completion with the Ladruno
fork's robust-solve driver (`robust_drive`) and the C++ solver toolkit it
orchestrates. Design rationale: [[31_ladruno_robust_solve_driver_adr]]. Session
map / cold-resume: [[robust_solve_handoff]]. Honest gotchas: [[LEDGER_quirks]].

> **One-line mental model.** `robust_drive` is an *outer loop over `analyze(1)`*
> that automatically resizes the increment and, when the solver stalls, escalates
> through a fixed **ladder of rungs** — each a heavier (and less "clean") way to
> get past the hard spot. It tells you honestly *how* it got through (the
> `verdict`), because a result that needed artificial dissipation or a dynamics
> excursion is **not** a publishable equilibrium.

---

## 1. Quick start

```python
import opensees as ops
from robust_drive import robust_drive          # Ladruno_scripts/ on sys.path

# ... build model, loads; set constraints/numberer/system/test/analysis once ...
ops.integrator("LoadControl", 0.1)             # nominal start; the driver resizes
ops.analysis("Static")

res = robust_drive(
    ops,
    done=lambda: ops.nodeDisp(2, 1) <= -0.010,  # termination predicate
    load_increment=0.1,                         # nominal load-factor step
    control=(2, 1, -5e-5),                      # (node, dof, du) -> enables rung-3
    stabilize=True,                             # enables rung-4
    dynamics=True,                              # enables rung-5
    log_path="run.jsonl")                       # auditable decision log (optional)

if res:                                         # truthy ONLY for verdict=="equilibrium"
    print("clean equilibrium:", res)
else:
    print("did not reach a clean equilibrium:", res.verdict, res)
```

`robust_drive` does **not** call `constraints`/`numberer`/`system`/`analysis` —
set those once before calling it (it re-issues only `integrator`/`algorithm`/`test`
as it works). It reaches the target iff `done()` becomes true.

---

## 2. The solver toolkit (what the driver orchestrates)

The driver is pure-Python orchestration; the hard numerics are shipped fork C++.
You can also use these directly without the driver.

| Tool | Command | Role | Class tag |
|---|---|---|---|
| **LadrunoArcLength** | `integrator LadrunoArcLength $arc $alpha <-adapt Jd lo hi> <-p e>` | adaptive arc-length radius (Ramm desired-iteration); **`-stabilize [f]`** adds Abaqus-style viscous regularization (rung-4) | 33004 |
| **LadrunoStabilizedUnbalance** | `test LadrunoStabilizedUnbalance $tol $maxIter` | true-equilibrium test for `-stabilize` (norms `‖λp−f_int‖` *without* the viscous `f_v`) | — |
| **LadrunoDynamicRelaxation** | `integrator LadrunoDynamicRelaxation <-mass …> <-dt $dt> <-damping kinetic\|viscous>` | matrix-free quasi-static dynamics (Gershgorin fictitious mass + Cundall kinetic damping); the rung-5 floor | 33005 |
| **LadrunoIndirectControl** | `integrator LadrunoIndirectControl $incr -dof $n $dof $coef …` | CMOD / indirect displacement control — snap-back path tracking | — |
| **Explicit integrators** | `integrator ExplicitBatheLNVD` / `CentralDifferenceLadruno` | true explicit dynamics (with the queryable critical `dt`) | — |
| **adaptive_static** | `ladruno_solve.adaptive_static(ops, set_step, n_steps)` | the increment-sizing wrapper (rungs 0-2) used standalone, integrator-agnostic | — |

The two **runtime query commands** that make the driver's decisions observable:

```text
ladrunoArcLength <sub>     # active LadrunoArcLength (StaticIntegrator):
   reduceStep f | increaseStep f | setArcLength v | scaleCVisc f   (actuators)
   arcLength | deltaLambdaStep | currentLambda | sign | deltaUstepNorm
   dissipationRatio | dissipatedEnergy | referenceEnergy           (rung-4 gate)

ladrunoDR <sub>            # active LadrunoDynamicRelaxation (TransientIntegrator):
   residualNorm | kineticEnergy                                    (rung-5 settling)
```

> **Python quirk (important).** In openseespy, any command that writes a WARNING
> to `opserr` raises `opensees.OpenSeesError` instead of returning `-1`. So
> `ladrunoArcLength('scaleCVisc', 0.0)` (factor must be > 0) and a bogus
> `ladrunoDR` subcommand **raise** — pre-validate args or wrap in `try/except`.
> (Tcl sees the `-1` return; Python does not.) See [[LEDGER_quirks]].

---

## 3. The rung ladder

The driver is a phase machine: `phase_implicit` (rungs 0-3) → `phase_stabilized`
(rung-4) → `phase_dynamics` (rung-5). Escalation happens at the **cutback floor**
(when the increment has been halved below `min_scale` and the algorithm ladder is
exhausted), in this order:

| Rung | Mechanism | Fires when | Verdict it can yield |
|---|---|---|---|
| 0 | increment cutback / grow (target-iteration sizing) | always (the spine) | `equilibrium` |
| 1 | line search | a step fails | `equilibrium` |
| 2 | modified / Krylov / Newton ladder | a step fails | `equilibrium` |
| 3 | **constraint switch** load → displacement control | floor reached **and** `control` given | `equilibrium` |
| 4 | **auto-stabilization** `LadrunoArcLength -stabilize` | floor reached, rung-3 unavailable/exhausted, `stabilize` on | `regularized` / `unverified` |
| 5 | **dynamics fall-through** `LadrunoDynamicRelaxation` + R-HANDOFF | rung-4 also stalls, `dynamics` on | `regularized` |

**Preference order matters.** Rung-3 (constraint switch) is the *clean* way past a
limit point and is preferred whenever a control DOF exists. Rungs 4 and 5 are
heavier last resorts for when there is no single control DOF (distributed
buckling, contact chatter) or the tangent is unusable.

---

## 4. `robust_drive(...)` parameter reference

| Parameter | Default | Meaning |
|---|---|---|
| `done` | — | `callable() -> bool`; True when the target state is reached |
| `load_increment` | — | nominal LoadControl load-factor step (resized by the spine) |
| `control` | `None` | `(node, dof, du)` enabling rung-3; `None` disables it |
| `algorithms` | KrylovNewton→Newton→NewtonLineSearch→ModifiedNewton | the rung-1/2 ladder (`ops.algorithm(*spec)` tuples) |
| `min_scale` | `1/1024` | cutback floor; below it the driver escalates |
| `grow` | `2.0` | scale multiplier after a clean step (capped at 1.0) |
| `peak_cutbacks` | `3` | consecutive no-commit cutbacks → early rung-3 switch (limit-point heuristic) |
| `max_substeps` | `20000` | global `analyze(1)` budget (R-THRASH guard) |
| `stabilize` | `False` | enable rung-4 |
| `stab_f` | `1e-3` | `-stabilize` dissipated-energy fraction. **Elevated** vs the Abaqus `2e-4` default (too weak to cross a hard limit); `-adaptStab` is deliberately NOT used (it re-weakens `c` and prevents crossing) |
| `stab_tol` | `1e-8` | `LadrunoStabilizedUnbalance` tolerance |
| `stab_hard_gate` | `0.05` | dissipation ratio above which the run is flagged `over_diffused` (informs the verdict; does **not** abort — crossing inherently dissipates) |
| `stab_rampdown_window` | `8` | clean stabilized steps before `scaleCVisc(0.5)` decays `c` (R-RAMPDOWN) |
| `stab_max_cutbacks` | `12` | stabilized `reduceStep` failures before escalating to rung-5 |
| `dynamics` | `False` | enable rung-5 |
| `dr_settle_tol` | `1e-4` | DR is quasi-static once `residualNorm < dr_settle_tol · ‖P‖` (mass-free) |
| `dr_max_steps` | `4000` | DR step budget per excursion |
| `dr_pref` | `None` | reference load norm `‖P‖` for the settle gate; `None` → proxy `max(1, |λ|)` |
| `dr_setup` | `None` | `callable(ops)` to configure the transient analysis (constraints/numberer/system) for DR; `None` → Plain/RCM/BandGen defaults |
| `log_path` | `None` | JSONL decision log path (one event per line; auditable) |
| `verbose` | `False` | print each decision |

### The `RobustResult`

`bool(res)` is **True only for `verdict == "equilibrium"`** — a stabilized or
dynamics result is falsy by construction, so a paper-grade consumer cannot mistake
it for clean.

| Field | Meaning |
|---|---|
| `completed` | did `done()` become true |
| `verdict` | `equilibrium` \| `regularized` \| `unverified` \| `incomplete` \| `aborted` |
| `mode` | last active mode (`LoadControl`/`DisplacementControl`/`Stabilized`/`Dynamics`) |
| `switches` | number of rung escalations |
| `substeps` | total `analyze(1)` calls |
| `min_scale_used` | smallest increment scale that was needed |
| `algos` | `{algorithm_name: success_count}` |
| `dissipation_ratio` / `stab_dissipated` | rung-4 `W_stab/Estrain0` and `W_stab` (the L0 energy-partition `SW` term) |
| `over_diffused` | rung-4 dissipation crossed `stab_hard_gate` |
| `peak_drift` | c-reduction diffusion bound (`None` = not computed → verdict stays `unverified`) |
| `dr_settled` / `dr_lambda` | rung-5 reached quasi-static; load factor frozen across the handoff |

### The honest verdicts

- **`equilibrium`** — reached via rungs 0-3 only. Truthy. Publishable.
- **`regularized`** — reached via rung-4 *with* a computed c-reduction drift, or a
  rung-5 DR rest state. A real solution point, but path-fidelity is qualified.
- **`unverified`** — reached via rung-4 *without* the diffusion bound (the current
  default, since that pass is deferred): "I got there, but I did not bound how
  much the artificial viscosity moved the answer." Never trust it as-is.
- **`incomplete` / `aborted`** — did not reach `done()`.

---

## 5. Worked examples

### Rung-3 — softening (the clean win)

A softening Concrete02 truss: load control dies at the strength peak (no
equilibrium above it); the constraint switch traces the descending branch.

```python
res = robust_drive(ops,
    done=lambda: ops.nodeDisp(2, 1) <= -0.010,
    load_increment=30.0/200, control=(2, 1, -5e-5))
assert res and res.verdict == "equilibrium"      # truthy; rung-3 traced softening
assert res.mode == "DisplacementControl"
```

### Rung-4 — engaging stabilization (control=None)

Same model with **no** control DOF: rung-3 is unavailable, so load control
escalates to stabilized load control. (Pure softening cannot actually be passed by
`-stabilize` — it *is* load control — so this honestly ends `incomplete`, having
engaged and refused to claim success.)

```python
res = robust_drive(ops,
    done=lambda: ops.nodeDisp(2, 1) <= -0.010,
    load_increment=30.0/200, control=None, stabilize=True, stab_f=1e-3)
assert res.mode == "Stabilized" and not bool(res)   # engaged, honest non-equilibrium
```

### Rung-5 — dynamics fall-through + R-HANDOFF

With `control=None, stabilize=False, dynamics=True`, the driver freezes the load,
relaxes the matrix-free dynamics to a quasi-static rest state (gated on the
mass-free `ladrunoDR residualNorm`), and restores the load factor **exactly**.

```python
res = robust_drive(ops,
    done=lambda: ops.nodeDisp(2, 1) <= -0.010,
    load_increment=30.0/200, control=None, stabilize=False, dynamics=True)
assert res.mode == "Dynamics" and res.dr_settled
assert ops.getTime() == res.dr_lambda            # R-HANDOFF: lambda restored exactly
assert res.verdict == "regularized" and not bool(res)
```

---

## 6. Honest limitations & gotchas

- **R-LOG-MASK — adaptive load control can *jump* a snap-through.** With step
  halving, plain load control can dynamically jump to the far stable branch of a
  snap-through (it "succeeds" while skipping the true unstable path). It therefore
  often never reaches rung-4 on snap-through models. Softening (no far branch) is
  the unambiguous load-control killer.
- **No clean-WIN fixture for rung-4/5 in the shipped battery.** Pure softening
  defeats stabilize too; snap-through gets jumped before rung-4. So rung-4/5 are
  exercised for *escalation + honest verdict + the handoff contract*, not a win.
  A distributed-buckling model with no control DOF would make rung-4 a genuine
  win; a tangent-pathology model static can't pass would make rung-5 one.
- **`-cVisc` is overwritten** by the first-commit calibration in the shipped
  ADR-20 `-stabilize` code (a latent quirk; the explicit-coefficient path does not
  take effect today). See [[LEDGER_quirks]].
- **Deferred** (tracked in the ADR + handoff): the rung-4 **c-reduction diffusion
  bound** (until built, stabilized successes stay `unverified`, never
  `regularized`-with-evidence) and the rung-5 **indirect-control polish tail**
  (re-land on the true branch with a few CMOD steps after a DR excursion).

---

## 7. Running the acceptance battery

```bat
set LADRUNO_DIST=...\dist\bin            REM the built opensees.pyd
set LADRUNO_OPENSEES_QUIET=1             REM silence the splash banner
python -m pytest Ladruno_scripts\robust_solve_tests -q     REM 15 green
```

Fixtures: `torture_softening.py` (clean load-control killer), `torture_snapthrough.py`
(limit point), `torture_stabilize.py` (the rung-4 dissipation seam),
`torture_dynamics.py` (rung-5 + R-HANDOFF). The suite skips cleanly if no `dist` is
present.

---

## 8. References

- ADR: [[31_ladruno_robust_solve_driver_adr]] — design, rejected alternatives, the
  full risk register (R-OBS, R-DIFFUSION, R-RAMPDOWN, R-HANDOFF, R-LOG-MASK …).
- Siblings: [[20_ladruno_arclength_stabilized_adr]] (the `-stabilize` seam),
  [[21_ladruno_dynamic_relaxation_adr]] (the DR floor),
  [[22_ladruno_dissipation_arclength_adr]] (the scoped-but-unbuilt snap-back keystone).
- Handoff / cold-resume: [[robust_solve_handoff]]. Gotchas: [[LEDGER_quirks]].
- Code: `Ladruno_scripts/robust_drive.py`, `Ladruno_scripts/ladruno_solve.py`,
  `Ladruno_scripts/robust_solve_tests/`.
