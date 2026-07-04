# ADR-65 Route B — P0 headroom oracle (RUN 2026-07-03)

Build-free go/no-go measurement for the **adaptive-Δt driver** (ADR-65 Route B).
Question the oracle answers: *on the softening/contact regime that motivates an
adaptive driver, how much wall-clock could a variable-Δt driver actually save by
growing Δt as stiffness degrades?*

Method: run `integrator CentralDifferenceLadruno -recompute 1` (tangent-tracked
critical step) on two 1-D truss bars where the tangent evolves, log
`criticalTimeStep()` every step, and compare adaptive step-counts against a
**constant-safety** baseline (a careful user's fixed `0.9·Δt_cr,elastic`). Holding
the safety factor equal is the whole trick — otherwise the "speedup" just measures
the safety-factor choice, not the nonlinearity.

- **Case A** — `Concrete01` softening bar, prescribed tip displacement (ramp past
  peak → unload → hold at degraded stiffness).
- **Case B** — `Steel01` (b=0.02) tip-force yield pulse (elastic → plastic flow →
  elastic unload → hold).

Run: `<pythoncore-3.12>\python.exe -S headroom_oracle.py` (needs the local
`dist/bin/opensees.pyd`; the script wires `sys.path` + `os.add_dll_directory` and
asserts it loaded the local build — see [[../../project_opensees_test_env]] for the
worktree `-S` trap).

## Result — growth headroom is NIL on reloadable nonlinearity

| case | max Δt_cr / elastic | SAFE grow (equal safety) | UNSAFE grow-to-tangent |
|---|:---:|:---:|:---:|
| A (concrete softening) | **12.51×** | **1.00×** | 1.19× (mirage) |
| B (steel b=0.02) | **7.07×** (=√(1/b), exact) | **1.00×** | 1.05× (mirage) |

Both cases: **0 invalid/nonpositive Δt_cr samples** across 2200 steps — the tangent
tracking works and stays positive even on the softening branch (it returns a
*larger* Δt_cr because |tangent| drops). That is exactly the trap:

- The **safe** adaptive step never exceeds `Δt_cr,elastic`, so at equal safety the
  growth speedup is **1.00× (nil)**. Reason: a reloadable material (plasticity,
  damage-with-reload) can stiffen back to `E0` in one step; a wave hitting the
  softened/yielded element reloads it elastic, so the only step safe against reload
  is the *undamaged* one.
- Trusting the instantaneous tangent Δt_cr over-reports the safe step by
  **√(E0/E_tan)** — measured 7× (steel, = √(1/b)) and 12× (concrete). A naive
  driver that grows to `0.9·criticalTimeStep()` here would blow up on the first
  elastic-reload wave. Case B's max ratio matching √(1/b) to 3 sig-figs is the
  internal-consistency check that the tangent tracking is real and correct.

## Verdict for ADR-65 Route B

The oracle **disproves the adaptive-growth value** for the motivating
softening/contact regime — those materials reload elastic, so growth headroom is
nil at equal safety. Route B's growth trigger (≥2× wasted steps) does **not** fire
on these cases.

Route B's *defensible* value is therefore **not** adaptive growth but:
1. the **shrink / safety** direction — protecting a fixed run when Δt_cr *drops*
   mid-run (contact engagement, geometric stiffening), which these monotone-softening
   cases don't exercise; and
2. **permanent, non-reloadable** stiffness loss (element erosion), where the softened
   step is genuinely safe because reload can't happen.

Not tested here (honest scope): a Δt_cr-*dropping* case, and non-reloadable erosion.
Those are where a positive Route B trigger would have to come from — build a
dropping-Δt_cr bundle before promoting Route B past accept-on-trigger.

Feeds the third quirk (tangent Δt_cr over-reports the safe explicit step for
reloadable materials) in `../LEDGER_quirks.md` and the Route B §P0 update in
`../65_ladruno_explicit_dt_strategies_adr.md`.
