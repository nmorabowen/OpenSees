# Damping sensitivity — FD vs DDM, and the zero-vanilla-touch route

Companion to [`12_damping_channels.md`](../12_damping_channels.md). That note maps
*where damping enters* OpenSees; this one answers *can we differentiate the
response w.r.t. damping*, and **how much vanilla code that costs**.

## TL;DR

- **The DDM (analytic) chain term for damping is already wired in vanilla.** The
  production `Newmark::formSensitivityRHS` (`SRC/analysis/integrator/Newmark.cpp:697`)
  already assembles `−(dC/dh)·v`, and the generic Rayleigh chain rule
  `dC/dh = αM·dM/dh + βK·dKt/dh + βK0·dK0/dh + βKc·dKc/dh` is implemented in
  `Element::getDampSensitivity` (`SRC/element/Element.cpp:610`). For an element
  that already supplies mass/stiffness sensitivity, **element-Rayleigh damping
  sensitivity is free** — *with a caveat*: the `βK` (current-tangent) and `βKc`
  (committed-tangent) terms call `getTangentStiffSensitivity`/
  `getCommittedStiffSensitivity`, which are **base-class stubs that silently
  return zero and warn "not implemented"** for most elements (`Element.cpp:530/550/569`).
  So "free" covers the `αM` and `βK0` (initial-stiffness) terms cleanly; the
  current/committed-tangent stiffness-proportional parts need per-element work.
- **The gaps are coverage, not framework:** nodal `αM` sensitivity is commented
  out (`Node::getDampSensitivity`, `SRC/domain/node/Node.cpp:1261` returns zero);
  the `SRC/damping/` `Damping` objects and the viscous/dashpot materials
  (`ViscousMaterial`, `DamperMaterial`, oil dampers) expose **no** `*Sensitivity`
  methods at all.
- **Autograd is the wrong tool *at whole-solver scope*.** AD-ing the *entire*
  forward response means re-typing every `double` on the forward path (elements,
  materials, `Matrix`/`Vector`, BLAS/LAPACK) — that touches essentially all of
  vanilla, more invasive than DDM (which is additive). The claim is scoped: a
  *local* dual-number/operator-overload wrapper around a **single fork-authored
  material/element** is fine and is exactly the recommended route below — that is
  AD applied narrowly, not autograd-through-the-solver.
- **The zero-vanilla-touch route that works today = finite differences.** Treat
  the forward solve as a black box, perturb the damping parameter, re-run,
  difference the response. Needs no derivative methods, no `Damping`/material
  changes, and works on **every** damping channel — including the ones with no
  DDM derivative. OpenSees even ships an in-engine version
  (`SRC/reliability/analysis/gradient/FiniteDifferenceGradient`), but a
  ~30-line Python wrapper is simpler and just as accurate.

## Decision

For damping **sensitivity studies / feasibility / a handful of parameters**: use
black-box finite differences. Zero vanilla edits, zero ledger debt, validated
exact below.

Reach for analytic derivatives only when FD's `N+1`-solves-per-gradient cost
bites (hundreds of parameters in a tight optimization/reliability loop). Two
sub-cases:
- **A new Ladruno element/material** → compute its `getMassSensitivity` /
  `getTangentStiffSensitivity` (and thus `getDampSensitivity` for free) inside
  the fork-authored class, optionally with a *local* operator-overloading AD tool.
  Plugs into the already-present vanilla DDM hook. **Still zero vanilla edits.**
- **A vanilla damping channel** (`Damping` objects, `ViscousMaterial`, nodal αM)
  → there is no no-touch analytic route; you must edit those vanilla files (and
  ledger them). FD remains the only no-touch option for these.

## Validation (all runs zero vanilla edits, pure openseespy + Python FD)

Run with the CPython 3.12 that matches `dist/bin/opensees.pyd`
(`pythoncore-3.12-64\python.exe`; see `project_opensees_test_env`). Each script
bootstraps the DLL dir itself.

All six OpenSees damping channels are now demonstrated (FD is **source-agnostic** —
it perturbs a knob and re-runs, never inspecting how the damping is implemented,
which is exactly why it covers the channels DDM cannot):

| Script / test | Channel(s) | What it proves | Result |
|---|---|---|---|
| `damping_fd_sensitivity_demo.py` | element Rayleigh `αM` | FD `dU/dξ` vs closed form `−F0/(2kξ²)` | **0.00%** at plateau; clean O(h²) |
| `viscous_damper_sensitivity_demo.py` | **`Viscous` dashpot** (no DDM derivative) | FD `dU/dC` vs `−F0/(C²ωₙ)` | **0.00%** at plateau |
| `damping_drift_2dof_demo.py` | element Rayleigh `[αM, βK]` | 2-DOF: `d(peak drift)/dξ` + gradient vector | A: plateau −0.5298; B: `[−1.68e-2, −4.18e+0]` |
| `remaining_channels_demo.py` | **modal**, **`UniformDamping` object**, **nodal `αM`** | FD gradient vs parameter-free oracle `dU/dp = −U/p` | **0.00–0.01%** at plateau, all three |
| `tests/test_fd_damping_sensitivity.py` | all of the above + `FDSensitivity` unit tests | Zone-A regression (class math + 5 channels + the UniformDamping quirk) | **11 passed** |

Channel coverage: element-Rayleigh `αM`/`βK` ✓, `Viscous` material ✓, modal ✓,
`UniformDamping` object ✓, nodal `αM` ✓. The sixth — **numerical (integrator)
damping** (HHT `α`, Newmark `γ`) — is reachable by the identical mechanism but is
algorithmic dissipation, not a physical source, so it is left documented-not-run.

The parameter-free oracle `dU/dp = −U/p` (exact whenever the resonant amplitude
`U(p) = A/p`, i.e. the damping coefficient scales linearly with the knob) validates
*all* the resonance channels at once — immune to per-channel closed-form typos and
to the fact that a `UniformDamping` object's **input `zeta` is not the realised SDOF
ξ** (band-fit operator: input 0.05 → realised ξ_eff ≈ 0.029). Two findings worth
noting: the FD caught a 10× typo in a hand-derived `dU/dαM` oracle, and the
`UniformDamping` input≠realised-ξ fact — both exactly the hand-derivative error
classes that argue *for* black-box FD over hand-coded DDM.

Representative output (Viscous damper — the no-DDM-derivative channel):

```
 dU/dC : central FD vs analytic -F0/(C^2 wn) = -1.574782e-01
   rel.step         FD dU/dC    rel.err
      1e-01    -1.590668e-01     1.01%   <- truncation
      1e-02    -1.574919e-01     0.01%   <- plateau (trust here)
      3e-04    -1.574761e-01     0.00%
```

### Bonus: FD surfaces real modeling quirks DDM would miss

The 2-DOF demo first reported `d(peak drift)/d(βK) = exactly 0.0`. Not a bug —
the FD gradient *correctly* detected that `zeroLength` defaults to
`-doRayleigh 0`, which **silently discards** stiffness-proportional Rayleigh
damping (see `12_damping_channels.md` / `LEDGER_quirks`). The response genuinely
did not depend on βK. A hand-coded DDM `dC/dh` would have happily returned a
nonzero derivative of a damping term *the analysis never applied* — a silent
correctness trap the black-box route is immune to. Fix: `-doRayleigh 1`.

## The reusable helper

`fd_sensitivity.py` is model-agnostic (no OpenSees dependency). Primary API is the
`FDSensitivity` class — it holds the `forward` callable + step config once,
**memoizes** forward evaluations (each is a full transient solve — the expensive
thing) so repeated gradients and the step-plateau sweep don't re-solve identical
models, and **counts** solves (`n_solves`) for honest cost accounting:

```python
from fd_sensitivity import FDSensitivity
fd = FDSensitivity(forward, rel_step=1e-2, scheme="central")  # forward(params)->scalar
grad    = fd.gradient(x0)          # vector; central = 2N solves, forward = N+1
plateau = fd.step_study(x0, comp=0)  # sweep the step to find the trust region
fd.n_solves                        # actual forward solves (cache hits excluded)
```

One-shot free functions `fd_gradient(forward, x0, ...)` / `step_study(...)` are kept
for quick use (they build a throwaway `FDSensitivity`). Point `forward` at any
parameter the model exposes: a damping ratio, raw `αM`/`βK`, a `Viscous` coefficient,
a modal/`UniformDamping` ζ, a Ladruno element/material property — the helper does
not care. Always confirm the step-size plateau before trusting a number.

## Limitations / when FD misleads (from adversarial review)

FD assumes the response is a **smooth** function of the parameter. Watch for:

- **Non-smooth response functionals.** A `max(|u|)`-over-time peak pick is only
  piecewise-smooth: if the timestep that holds the peak *switches* between the `+h`
  and `−h` runs, or under hysteretic/contact/yielding response, the difference
  quotient is noisy or ill-defined. The validated demos are linear + at resonance
  (smooth), so this doesn't bite here — but a nonlinear model needs a smoother
  functional (energy, RMS, response at a fixed time) or a wider step.
- **Steady-state proxy needs the transient dead.** The "max over last 10 cycles ≈
  steady-state amplitude" identity relies on `ξ·2π·(n_cycles − tail_cycles) ≫ 1`.
  At ξ=0.08–0.10, n_cycles=120 the transient is killed to ~1e-24; but for **light
  damping** (e.g. ξ=0.005) the same run leaves ~3% residual transient in the tail
  window, biasing the measured amplitude **high**. Lowering ξ or n_cycles → re-check.
- **Sampling bias is benign but real.** Discrete `max` over dt-spaced samples
  under-reads the true crest by ≈`cos(π·dt/period)` (~2e-5 at 500 samples/cycle).
  It's a near-constant multiplicative factor that **cancels in the FD ratio** —
  which is why the *gradient* error (0.00%) beats the *amplitude* error (0.01%).
- **"Works on every channel" is true at the channel level.** Modal damping, rate-
  dependent path-memory materials, etc. are all differentiable by FD — but they
  inherit the smoothness caveat above just like any other functional.

## Cost / when this stops being enough

`N+1` forward solves per gradient. Fine for a few damping parameters. If you ever
need gradients w.r.t. many parameters inside an optimizer or FORM/SORM reliability
loop, that linear-in-parameters cost is when exact DDM (for elements that have it)
or AD-inside-a-Ladruno-class earns its keep.
