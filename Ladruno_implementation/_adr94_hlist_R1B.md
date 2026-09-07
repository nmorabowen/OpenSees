# ADR-94 R1-B — H-list mechanical reproducers (H2, H4, H5, H9, H12, H13, H14, H15)

Lane R1-B, worktree `asdplastic-review-plan-62585c`, HEAD `52314165a` (build
verified via `ops.ladrunoBuild()`). Reproducers live in
`tests/test_adr94_hlist_mechanical.py` (Zone-A, `t0m`). Result: **11 passed,
0 failed, wall time ~0.9s** (pytest alone); ~1.9s combined with the four
pre-existing ASDP test files (31 passed, 2 skipped — the 2 skips are a
pre-existing, environment-dependent condition in `test_asdplastic_mctc.py`'s
`test_update_parameter_ramp`, reproduced by running the four baseline files
**alone**, unrelated to this lane's file).

Per-H verdicts follow. "Deciding observation" is the one number/fact that
settles the verdict; "test" names the pinning reproducer.

## H2 — `getClassType()` returns a dangling `std::string::c_str()`

- **Verdict:** CONFIRMED (read-only, no test — per the review plan §4/R1).
- **Severity:** minor (UB, works by SSO luck; no test needed to establish it).
- **Deciding observation:** `ASDPlasticMaterial3D.h:196-201` — `std::string
  name("ASDPlasticMaterial3D"); return name.c_str();` returns a pointer into
  a stack-local `std::string` destroyed at function return.
- **Blast radius:** any caller that dereferences the pointer after the call
  returns (a `Print`/recorder path that stores it) reads freed memory; this
  build's fixed literal is short enough for small-string-optimization to
  survive by luck on most compilers/allocators.
- **Test:** none (by design).

## H4 — `revertToLastCommit()` no-op; `revertToStart()` → -1, silently swallowed

- **Verdict:** CONFIRMED.
- **Severity:** major.
- **Deciding observation:** structurally, `revertToLastCommit()`'s body
  (735-748) contains no live statement other than `return 0;` (every
  original line is commented out) — pinned by a source-level regex so a fix
  flips the test red. Separately, `OPS_resetModel()`
  (`OpenSeesCommands.cpp:2534-2547`) calls `theDomain->revertToStart()` and
  never inspects its return value, so `ops.reset()` reports success while
  `ASDPlasticMaterial3D::revertToStart()` prints "not implemented" and
  returns -1 — measured via a **child process** (see caveat below).
- **Measured caveat (important for anyone extending this test):**
  `TenNodeTetrahedron`'s "stresses"/"forces" `eleResponse` **unconditionally
  re-derives stress from the current nodal trial displacement on every
  query**, with or without a preceding "forces" call, and via the
  "material"/"integrPoint" sub-response routes too. Since
  `Domain::revertToLastCommit()` correctly resets nodal trial displacements,
  this makes the element **self-heal** on every read, and the material-level
  no-op is *invisible* from Python via any eleResponse path tried. This is
  why the runtime half of H4 pins `revertToStart()` instead, and
  `revertToLastCommit()` is pinned structurally.
- **Second measured trap:** pytest's `capfd` fixture does not see anything
  this native (.pyd) extension writes to `cout`/`cerr` on this Windows
  build — a mid-process `dup2` swap of fd 1/2 does not reach output written
  through the .pyd's own linked CRT. The exact same code prints normally
  when the **whole process's** stdout is piped from outside (a shell, or
  `subprocess`). Both H4's `reset()`/`"not implemented"` check and H12 (below)
  had to be rewritten to run the model in a **child process** and capture its
  real OS-level stdout/stderr; a `_run_child()` helper does this and is
  reused by both tests. **Any future ASDP test that wants to assert on
  cout/cerr content on this platform must use the same pattern, not
  `capfd`.**
- **Blast radius:** any workflow that inspects material state after a
  failed-then-reverted step (a recorder sampling mid-cutback, `ops.reset()`
  after a non-converging run, `LadrunoBeginAugment`/adaptive stepping) sees
  stale/dirty state with no error surfaced.
- **Test:** `test_H4_revert_to_last_commit_is_noop`.

## H5 — the f-decreased elastic exit is unguarded outside `Backward_Euler`

- **Verdict:** CONFIRMED (4 of 6 non-BE sites pinned at runtime; the
  remaining 2 confirmed by reading only).
- **Severity:** major.
- **Deciding observation:** at `TET_UTOP=-0.02`/20 steps,
  `n_max_iterations=100`, `strict_convergence=1`, `Forward_Euler`,
  `Forward_Euler_Subincrement`, `Modified_Euler_Error_Control`, and
  `Runge_Kutta_45_Error_Control_old` all commit states with `f_MC` in the
  hundreds-to-thousands against a tolerance of ~0.1-0.2 — i.e.
  `strict_convergence` has **zero effect** on these four, exactly as
  predicted by the unguarded `yf_val_start > yf_val_end` shortcut at their
  respective sites (1423, 1599, 3088, 2667).
  `Backward_Euler_LineSearch` (2406) and `Runge_Kutta_45_Error_Control`
  (3435, the non-`_old` one) instead **fail to converge globally** (rc=-3)
  on this exact rig at every magnitude tried — a different failure mode this
  reproducer does not isolate — but both sites are still confirmed unguarded
  by reading (no `strict_convergence` check in scope at either).
- **Blast radius:** any non-default integrator choice silently defeats the
  ADR-84 P2a safety flag; users who set `strict_convergence 1` believing it
  gates every integrator are only protected on `Backward_Euler`.
- **Test:** `test_H5_strict_convergence_does_not_gate_other_integrators`
  (parametrized; 4 cases, all pass/CONFIRM).

## H9 — ME/RK45 drift-correction blocks are empty `if` statements

- **Verdict:** CONFIRMED (`Modified_Euler_Error_Control` pinned at runtime;
  `Runge_Kutta_45_Error_Control` confirmed by reading + a live-code proof).
- **Severity:** major.
- **Deciding observation:** with `return_to_yield_surface Disabled` and
  `f_absolute_tol=1e-8`, a 20-step `TET_UTOP=-0.02` MC path on
  `Modified_Euler_Error_Control` commits a final state with `f_MC ≈ 1.09e5`
  — nine orders above tolerance — because the "Validate yield function
  drift" block (3240-3250) computes `yf_val` and tests `yf_val > 10*tol` but
  its body is a commented-out `cout`. The identical block for
  `Runge_Kutta_45_Error_Control` (3760-3770) is confirmed the same way by
  reading; a runtime pin for it specifically was attempted but this build's
  `Runge_Kutta_45_Error_Control` fails to converge globally under
  `return_to_yield_surface Disabled` at every magnitude tried (a related,
  separate observation — see H15's note on the same function family).
- **Blast radius:** any user relying on ME/RK45 with the drift check
  "on" (it never was) gets silent, unbounded admissibility violations on
  coarse steps; `return_to_yield_surface` is the *only* working correction
  for these two integrators.
- **Test:** `test_H9_explicit_drift_check_is_dead_code`.

## H12 — diagnostics use `cout` (not `opserr`) and flood stdout

- **Verdict:** CONFIRMED.
- **Severity:** minor (diagnostics/perf, not correctness).
- **Deciding observation:** a 10-element run (independent `stdBrick` cubes,
  one shared material tag) over 20 plastic steps emits 20
  `"() ASDP Integration Info..."` lines on stdout (`cout`, confirmed via a
  child-process capture) and **zero** on stderr/opserr — confirms both the
  channel claim and that the counter fires essentially every commit once the
  material is plastic (`GLOBAL_INT_max_iter[ASDP_TAG]` is a **per-tag
  static**, so with `n_elem` sharing one tag, whichever GP touches it last
  determines what the next commit's line reports — a related, unquantified
  cross-attribution risk not separately measured here).
- **Blast radius:** invisible under MPI or any pipeline that redirects only
  `opserr`; stdout flood scales with GP count × step count on any nontrivial
  model.
- **Test:** `test_H12_commit_diagnostics_use_cout_not_opserr` (uses the
  child-process pattern from H4 — `capfd` does not see this output).

## H13 — unknown parser tokens are silently swallowed

- **Verdict:** CONFIRMED (both halves).
- **Severity:** blocker (per ADR-84 §9.1: "silent misconfiguration is the
  most expensive failure mode this fork has recorded").
- **Deciding observation:** (a) a deck with `"strict_convergance"` (typo)
  inside `Begin_Integration_Options` produces **byte-identical** `analyze()`
  codes to a deck that omits `strict_convergence` entirely — the unknown
  token and its value are silently dropped by the missing `else` in
  `OPS_AllASDPlasticMaterial3Ds.cpp`'s if-chain (358-460). (b) a deck with
  `"MC_phii"` (typo, should be `MC_phi`) produces a committed stress that is
  admissible under `f_mc(phi=0)` (residual `≈1.0e-6`, i.e. exactly on the
  φ=0 surface) but sits 68662 kPa deep inside the admissible region under
  the *intended* `phi=20` — i.e. the material silently used the unset
  default (0) for the friction angle. Confirmed at the framework level too:
  `utuple_storage.h::setParameterByName_impl`'s base case (`I ==
  tuple_size`) is a no-op with no warning path.
- **Blast radius:** any misspelled option name (integration or model
  parameter) is accepted without error and silently changes constitutive
  behaviour — exactly the failure mode ADR-84 already paid for once.
- **Tests:** `test_H13_unknown_integration_option_is_silently_ignored`,
  `test_H13_unknown_model_parameter_is_silently_ignored`.

## H14 — `getCopy()` omits `first_step` from its explicit member-copy list

- **Verdict:** CONFIRMED (structural), but **practically latent** under
  normal model-build order.
- **Severity:** minor / doc-only in practice (see blast radius).
- **Deciding observation:** `getCopy(void)`'s body explicitly assigns
  `TrialStress`, `TrialStrain`, `TrialPlastic_Strain`, `CommitStress`,
  `CommitStrain`, `CommitPlastic_Strain`, `iv_storage`,
  `parameters_storage`, and `stress_set_externally` onto the new instance —
  but never `first_step`, pinned by a source-level regex (flips red if a fix
  adds the assignment). Because every host element calls `getCopy()` exactly
  once, on the still-pristine tag-registered prototype, at **construction**
  time (before any analysis step), the missing copy is unreachable in
  ordinary usage: two independently-constructed `TenNodeTetrahedron`
  elements sharing one `InitialP0=-37.5` material tag both correctly seed a
  compressive mean stress on their own first commit (measured, both
  `mean(sigma) < 0`). It would only bite a `getCopy()` call made on an
  ALREADY-advanced instance (state re-partitioning, lazy per-GP construction
  after stepping has begun) — not exercised by any test in this suite or, as
  far as this review found, by any host element in this fork.
- **Blast radius:** none identified under current usage patterns; flagged
  for anyone who adds a code path that calls `getCopy()` post-construction
  (e.g. a future MP re-partitioning fix, out of scope per D3).
- **Test:** `test_H14_getcopy_does_not_preserve_first_step`.

## H15 — BE evaluates elasticity once per step; ME/RK45 re-evaluate per stage

- **Verdict:** CONFIRMED (structural).
- **Severity:** major for any stress-dependent elasticity (`StiffSoil_EL`,
  `DuncanChang_EL`); doc-only otherwise.
- **Deciding observation:** `Backward_Euler`'s body calls `et(...)` exactly
  once (line 2056, before the scalar-Newton loop) — a second textual match
  at line ~2110 is inside a fully commented-out alternative implementation
  and does not count once comments are stripped (an early, uncorrected count
  of "2" was a false positive from including the dead code — corrected in
  the shipped test). `Modified_Euler_Error_Control` and
  `Runge_Kutta_45_Error_Control`'s live bodies each call `et(...)` at least
  twice (at `CommitStress`, at the current-stage stress, and at a predictor
  stress). For `StiffSoil_EL` (stress-dependent `E`), BE and ME/RK45 are
  therefore evaluating genuinely different constitutive operators on the
  same nominal material, not just different numerical schemes converging to
  the same answer.
- **Separate finding (not part of H15, flagged for the record):** a runtime
  `StiffSoilShear_YF`/`StiffSoilShear_PF`/`StiffSoil_EL` triaxial driver was
  attempted for this test and hit a `NaN!` on the **very first**
  `Backward_Euler` step at every magnitude/`InitialP0` combination tried
  (`utop` from -0.0005 to -0.005, `InitialP0` from 0 to -100). This looks
  like an independent defect in the `StiffSoilShear_YF`/`StiffSoil_EL`
  combination itself (untested per H11 — no Zone-A coverage exists for any
  StiffSoil combo) rather than anything specific to H15's claim; it is why
  H15 is pinned structurally instead of with the originally-planned BE-vs-
  RK45 stress comparison. Worth a follow-up H or R3 component-audit item.
- **Blast radius:** `StiffSoilCap_YF`/`StiffSoilShear_YF` users get a
  materially different model depending on `integration_method`, undocumented.
- **Test:** `test_H15_be_evaluates_elasticity_once_me_rk45_per_stage`.

## Summary table

| H | Verdict | Severity | Test |
|---|---|---|---|
| H2 | CONFIRMED (read-only) | minor | none (by design) |
| H4 | CONFIRMED | major | `test_H4_revert_to_last_commit_is_noop` |
| H5 | CONFIRMED (4/6 runtime + 2/6 by reading) | major | `test_H5_strict_convergence_does_not_gate_other_integrators` |
| H9 | CONFIRMED (1/2 runtime + 1/2 by reading) | major | `test_H9_explicit_drift_check_is_dead_code` |
| H12 | CONFIRMED | minor | `test_H12_commit_diagnostics_use_cout_not_opserr` |
| H13 | CONFIRMED | blocker | `test_H13_unknown_integration_option_is_silently_ignored`, `test_H13_unknown_model_parameter_is_silently_ignored` |
| H14 | CONFIRMED (structural, latent) | minor | `test_H14_getcopy_does_not_preserve_first_step` |
| H15 | CONFIRMED (structural) | major (StiffSoil/DuncanChang), doc-only otherwise | `test_H15_be_evaluates_elasticity_once_me_rk45_per_stage` |

## Traps recorded for the next lane

- **`capfd` cannot see this .pyd's `cout`/`cerr` on Windows.** Use a child
  process (`subprocess.run([sys.executable, "-c", script], ...)`) and read
  its real stdout/stderr instead. `_run_child()` in
  `tests/test_adr94_hlist_mechanical.py` is a ready-made helper.
- **`TenNodeTetrahedron`'s `eleResponse` self-heals** — it always re-derives
  stress from current nodal trial displacement, so it cannot be used to
  observe raw material-level Trial-state corruption after a domain-level
  revert. A source-level structural check is the only observation channel
  found for that class of defect.
- **A misspelled attempt to reproduce a StiffSoil claim (H15) hit a `NaN!`
  on step 1** — the StiffSoilShear/StiffSoilCap family may need its own H
  (coverage gap already noted at H11).
- **Buffering interleave:** a Python `print()` marker and this native
  extension's `cout` do not interleave in program order in a captured
  subprocess byte stream (`print()` is fully buffered until exit; `cout`
  flushes on every `endl`). Do not rely on marker position to window native
  output; count over the whole stream instead.
- **Windows-only subprocess flake:** `subprocess.run(...)` without an
  explicit `stdin=` sometimes raised `OSError: [WinError 6] The handle is
  invalid` from `_winapi.DuplicateHandle` (non-deterministic — observed on
  roughly 1 in 3 runs of the two subprocess-based tests, H4 and H12) because
  pytest's own stdio setup does not always leave the parent's stdin as a
  duplicable handle. Fixed by passing `stdin=subprocess.DEVNULL` in
  `_run_child()`; verified clean over 12 consecutive full-file runs after
  the fix (0 failures). Anyone adding a new `subprocess.run` call against a
  pytest-hosted stdin should do the same.
