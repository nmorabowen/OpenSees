---
title: Explicit Noh–Bathe integrators + energy-balance recorder
project: Ladruno
status: implemented
priority: high
owner: nmora
tags:
  - implementation
  - integrator
  - recorder
  - explicit
---

# Explicit Noh–Bathe integrators + energy-balance recorder

> **ADR** (Architecture Decision Record). Retrospective record of the explicit-dynamics
> diagnostics subsystem: the `EnergyBalanceRecorder`, the `ExplicitBathe` /
> `ExplicitBatheLNVD` integrators, and the critical-time-step (`dt_cr`) machinery.
> PR #2 (recorder + integrators) **merged to `ladruno` 2026-05-30**; PR #3 (queryable
> `dt_cr` + tangent/periodic recompute + aliasing fix) **open**. Captures the settled
> decisions and — importantly — the marked extension points for future work.

## What

A first-class, *verified* explicit-dynamics core "with discipline" (per
[[Ladruno_explicit_roadmap]]):

- **`EnergyBalanceRecorder`** (classTag 26) — per step, for the whole model and
  optionally per `MeshRegion`: kinetic energy `KE`, internal energy `IE`, damping
  work `DW`, unbalanced-load work `ULW`, the closure residual `RES = ULW−(KE+IE+DW)`,
  and a normalized `ERR%`. Plain-text sidecar output.
- **`ExplicitBathe`** (classTag 61) — the Noh–Bathe two-sub-step explicit scheme
  (2nd order, controllable high-frequency dissipation via `p`).
- **`ExplicitBatheLNVD`** (classTag 63) — same scheme + FLAC local non-viscous
  damping for dynamic relaxation / quasi-static solving.
- **`dt_cr` machinery** — per-element critical-time-step estimate, queryable from
  Python/Tcl (`criticalTimeStep()`), optionally tangent-based and periodically
  refreshed.

**Not in scope** (live in the roadmap / future): mass scaling, hourglass-energy
monitoring, contact/penalty energy, explicit shells, sub-cycling *inside* the
integrator, MPCO/STKO visualization of energy.

## Why

Explicit integration is the roadmap's path to robust collapse / contact / wave /
SSI problems (implicit convergence is fragile there). But explicit is "silently
wrong" without diagnostics — there is no Newton residual to warn you. So the
deliverables are inseparable: the integrators **and** the energy-balance
instrument that proves a run stayed physical (energy closed to ~1%), plus a sound
`dt_cr` so users can pick a stable step. `ExplicitBatheLNVD` additionally gives a
quasi-static path (dynamic relaxation) that the implicit static solvers can't
always reach (snap-through, near-singular tangents).

## Decisions (settled 2026-05-30)

| #  | Decision | Rationale | Consequence / extension point |
|----|----------|-----------|-------------------------------|
| D1 | Energy balance writes a **plain-text sidecar**; MPCO is **left frozen** | `MPCORecorder` is byte-identical to upstream (deferred even a crash fix to keep it so); MPCO has only `ON_NODES`/`ON_ELEMENTS`, no global/scalar bucket; STKO can't render a non-nodal/element result | Energy parsed with numpy/pandas, **not** STKO_to_python; not inside the `.mpco`. **Extend:** an `ON_REGIONS/ENERGY_BALANCE` group only if STKO/[[03_mpco_ladruno]] gains global-result support |
| D2 | Energy = **incremental work integrals + closure residual** (`KE` instantaneous; `IE`/`DW`/`ULW` trapezoidally integrated; `RES`, `ERR%`); per-region + global | The residual is the actual engineering diagnostic; incremental `IE` makes `KE+IE` cancel correctly for free vibration; `ERR%` normalized by summed-component magnitudes so it doesn't blow up when terms cancel | **Coverage gaps (documented; residual exposes them):** modal damping (applied in the solve, not via `getDamp`) and element/distributed `eleLoad` loads are NOT captured. **Extend:** independent external-work-input channel from load patterns; elastic-vs-dissipated `IE` split |
| D3 | Keep the tree's **refined `ExplicitBathe`**; **wire up LNVD's FLAC damping** (was commented-out dead code) | The in-tree base was more evolved than jaabell's; LNVD's entire reason to exist was disabled | LNVD applies the local damping **symmetrically to both sub-steps** via a virtual `formUnbalance()` override (the algorithm calls it for sub-step 1, `update()` for sub-step 2); reuses `ExplicitBathe`'s eigensolve |
| D4 | `dt_cr` = **per-element generalized eigensolve** (`K v = λ M v`, DGGEV); report the **conservative central-difference limit `2/ω`** AND the **Noh–Bathe `~2×`** bound | Per-element avoids a global eigensolve (theorem: global `ω_max ≤ max element ω_max`); the CD limit is provably safe for Noh–Bathe | **Caveats:** ignores constraints/MP; row-sum lumping is non-conservative for rotational DOFs (beams/shells); needs **element** mass+stiffness (pure nodal-mass models → `dt_cr ≤ 0`). **Extend:** `DSYGV` (symmetric pencil, ~2× faster, no complex-eigenvalue fragility); a cheap `ℓ_e/c_e` per-element estimate; constraint-aware bound |
| D5 | Time-history adaptivity is **driver-level sub-stepping**, NOT integrator-internal sub-cycling | Each real step must **commit** material state (path-dependent plasticity/damage evolves), and OpenSees couples `commitDomain()` = material-commit + recorder-fire + clock. Internal sub-cycling would corrupt nonlinear history (commit-once) or flood recorders (commit-each) | User does `analyze(n_sub, dt/n_sub)` with `n_sub = ceil(dt/(safety·dt_cr))` (queried). Correct for nonlinear; loads sub-sampled. **Extend:** a built-in sub-cycling integrator only **if** the commit/record coupling is decoupled framework-side |
| D6 | Expose `dt_cr` via a **virtual `getCriticalTimeStep()` on `TransientIntegrator`** + a `criticalTimeStep()` command (Py/Tcl) | Avoids downcasting in the interpreter; base default `-1.0` leaves the ~30 other integrators untouched | Requires `-cfl`/`-tangent`/`-recompute` AND ≥1 step before a valid value. **Extend:** compute `dt_cr` in `domainChanged()` so the query is valid pre-`analyze` (drops the `analyze(1,1e-9)` priming); distinct sentinels for not-applicable vs not-yet-computed |
| D7 | **Register classTags** 26/61/63 in `FEM_ObjectBrokerAllClasses` + `TclPackageClassBroker` | `openseessp` / database-restart reconstruct objects via the broker | Benign for `openseesmp` (each rank constructs locally from the command) |
| D8 | **Lesson (cross-cutting):** copy an element matrix **before** fetching another from the same element | Many elements (Truss, FourNodeQuad, DispBeamColumn, anything using base `Element::getDamp` which calls `getMass`/`getTangentStiff`) return a reference to a **shared static `theMatrix`**; the second call clobbers the first reference | Caught **two real bugs** — `computeCriticalTimestep` (`dt_cr=inf` for all distributed-mass models) and `EnergyBalanceRecorder::addElementEnergy` (KE from the damping matrix, ~0). Both fixed by value-copying the mass. **Apply this pattern anywhere two element matrices are held at once** |

## Where

- New: `SRC/recorder/EnergyBalanceRecorder.{h,cpp}`,
  `SRC/analysis/integrator/ExplicitBatheLNVD.{h,cpp}`
- Modified: `SRC/analysis/integrator/ExplicitBathe.{h,cpp}`,
  `SRC/analysis/integrator/TransientIntegrator.h` (new virtual),
  `SRC/interpreter/{OpenSeesCommands.{h,cpp},PythonWrapper.cpp,TclWrapper.cpp}`
  (the `criticalTimeStep` command), `SRC/recorder/TclRecorderCommands.cpp`,
  `SRC/interpreter/OpenSeesOutputCommands.cpp`, `SRC/tcl/commands.cpp`,
  `SRC/classTags.h`, the two brokers, plus `CMakeLists.txt`/`Makefile`.
- Reference patterns: `CentralDifference.cpp` (explicit contract), `NodeRecorder.cpp`
  (recorder/column metadata), `OPS_numIter`/`OPS_systemSize` (query command).
- Build: no new target/dep. Full installer build still **blocked by the
  [[03_mpco_ladruno]] `MPCORecorderLadruno` link error** — see [[../Ladruno_internal/01_compilation_journal]].
- Examples: `examples/explicit_bathe_example.py` (recipe + LNVD + adaptive sub-step),
  `Ladruno_scripts/_verify_explicit.py` (9-test numerical battery).

## How

**Explicit recipe (required):** lumped/element mass, `system Diagonal`,
`algorithm Linear` (exactly 2 solves/step), `dt < dt_cr`.

```python
# diagnostics-instrumented explicit run
ops.integrator('ExplicitBathe', 0.54, '-cfl')      # p=0.54; -cfl computes dt_cr
ops.recorder('EnergyBalance', '-file', 'e.txt', '-time', '-region', 1)
# time-history adaptive sub-stepping (driver-level, D5)
ops.analyze(1, 1e-9)                                # prime: triggers dt_cr compute
n = max(1, ceil(dt / (0.9 * ops.criticalTimeStep())))
ops.analyze(n, dt/n)
# quasi-static relaxation
ops.integrator('ExplicitBatheLNVD', 0.54, 0.8)     # alpha=0.8 (FLAC default)
```

Options: `-verbose` (gate per-step output), `-cflAbort` (stop if `dt`>NB limit),
`-divergence f` (abort on KE growth), `-tangent` / `-recompute N` (track stiffness
changes). Stability constant `EB_NB_STABILITY_FACTOR = 2.0` (Noh & Bathe 2013).

**Verification** — `_verify_explicit.py`, 9/9: order-of-accuracy slope **1.99**;
stability `ωΔt ≥ 3.0` vs central difference 1.95; cross-check vs Newmark/CD <2.6%;
spectral damping light/well-behaved; 1D wave speed exact `= √(E/ρ)`; LNVD→static
field match 0.0%; rigid-body momentum exact; `criticalTimeStep() = ℓ/c` exact;
recorder element-mass KE correct.

## Future work / extension points

Ordered roughly by value; many are **roadmap-gated** (need element/integrator hooks).

1. **Mass scaling (selective).** The single biggest efficiency lever for 3D/SSI.
   `dt_cr` already *identifies the critical element*; the EnergyBalance recorder is
   the added-mass monitor. Design seam: inflate the diagonal `M` for sub-target
   elements, feed the same lumped `M` the integrator reads, recompute `dt_cr` on the
   *scaled* mass, hard-warn >5% added mass. Roadmap §5.1. **Separate integrator/decorator.**
2. **`dt_cr` correctness/robustness.** ~~`DSYGV` for the symmetric pencil + relative-`β`
   threshold (A3 from review)~~ **done (2026-05-30; DGGEV fallback retained; the β fix
   corrected a silent `dt_cr→0` bug)**; ~~`-lump` for rotational DOFs~~ **partial — `-lump
   diagonal` added (diagonal-of-consistent; a true mass-conserving HRZ is still TODO)**;
   a `ℓ_e/c_e` per-element wave-speed estimate (cheaper; still **deferred**); compute in
   `domainChanged()` (drops the prime, D6) — **deferred** (would rework PR #3's merged
   state model + the public `criticalTimeStep()` contract); make `-recompute` *enforce* —
   **already satisfied**: PR #3's `-cflAbort` aborts in `newStep` against the refreshed
   `dt_cr`.
3. **Energy recorder additions.** External-work-input channel (independent check);
   elastic-vs-dissipated `IE` split (needs a standardized element/material "recoverable
   vs dissipated energy" response); hourglass / contact / added-mass energy columns
   (need element/integrator hooks — roadmap §5.3/§5.7/§5.1). Reserve the response-string
   contracts now so those subsystems emit them when built.
4. **Parallel correctness.** ~~Cross-rank `dt_cr` reduction~~ **code added (2026-05-30,
   guarded `MPI_Allreduce(MPI_MIN)`); NOT yet runtime-verified — confirm `OPS_Analysis`
   carries a parallel macro in the `openseesmp` build (the sequential build does not, so
   the reduce is inert there).** MPI reduce for a global energy balance still per-partition.
5. **Sub-cycling integrator.** Only after the commit/recorder coupling (D5) is solved
   framework-side; otherwise stays driver-level.
6. **MPCO/STKO energy visualization.** Blocked on STKO support for global/region results
   (D1); revisit with [[03_mpco_ladruno]].

## Risks / open questions

> [!question]
> Row-sum lumping for rotational DOFs (beams/shells) can be non-conservative — should
> `dt_cr` switch to an HRZ/diagonal-of-consistent lumping, or an `ℓ_e/c_e` estimate, for
> those element classes?

> [!question]
> Should `-recompute`/`-tangent` *enforce* a smaller `dt` (sub-cycle or abort) rather than
> only reporting? Today only `-cflAbort` enforces; everything else is advisory + driver-side.

- ~~The `extern computeCriticalTimestep` is shared between `ExplicitBathe.cpp` and
  `ExplicitBatheLNVD.cpp` (hand-copied) — fragile; hoist to a small header.~~ **Resolved
  2026-05-30:** hoisted to `SRC/analysis/integrator/CriticalTimeStep.{h,cpp}`.
- `-lump diagonal` (diagonal-of-consistent) avoids the row-sum rotational pathology but is
  **NOT mass-conserving** — it can be non-conservative for translational DOFs. A true
  mass-conserving **HRZ** lump is the proper follow-up.
- `dt_cr` ignores constraints/MP and pure nodal mass; safe-but-non-binding there.
- Full build / installer still blocked by the [[03_mpco_ladruno]] link error.

## Implementation log

- 2026-05-29 — energy-balance recorder + integrator port, 3-agent review, hardening.
- 2026-05-30 — runtime smoke + 7-test numerical battery (caught the OpenSeesPy arg-parsing
  bug); **PR #2 merged**. Queryable `dt_cr` + tangent/recompute; runtime test caught the
  `dt_cr` mass-aliasing bug; **PR #3 open**. Pre-merge 3-agent review caught the same
  aliasing bug in `EnergyBalanceRecorder::addElementEnergy`; fixed + polish (PR #3, 9/9).
- 2026-05-30 — **PR #3 + #4 merged to `ladruno`** (queryable `dt_cr`, advisory
  `-tangent`/`-recompute`, `criticalTimeStep()` command + base virtual, recorder
  aliasing fix, ADR docs).
- 2026-05-30 — **`dt_cr` hardening reconciled on top of PR #3** (a follow-up branch was
  started in parallel, before PR #3 landed; reconciled via reset-and-reapply, not rebase,
  after a 3-agent overlap/design/mechanics analysis). Shipped: shared
  `CriticalTimeStep.{h,cpp}` (hoist; kills the `extern`); **`DSYGV`→`DGGEV` fallback +
  relative-β threshold** — this fixes a *silent* bug in the merged code where PR #3's
  exact `beta != 0` test let a near-massless DOF drive `dt_cr → 0` and abort a stable run
  under `-cflAbort`; `-lump rowsum|diagonal` (diagonal-of-consistent; non-conservative,
  not true HRZ — see risks); guarded `MPI_Allreduce(MPI_MIN)` cross-rank reduction.
  **Deliberately not changed** (PR #3's just-merged choices kept): enforcement stays in
  `newStep` under `-cflAbort` (PR #3 already aborts there; our `commit()`/`-7` variant
  dropped); single `<=0` sentinel kept (distinct sentinels deferred); compute stays in
  `newStep` (domainChanged-compute deferred — it would rework the merged state model and
  the public command's contract). Each TU compile-verified standalone via `cl.exe`;
  `_verify_explicit.py` carries tests 8/10–12 (query, `<=0` contract, `-lump`, enforcement).
