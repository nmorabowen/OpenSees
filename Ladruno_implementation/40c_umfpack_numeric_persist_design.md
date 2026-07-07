---
title: UMFPACK Numeric-persist — factorization reuse across solves (ADR-40 rank 8/10)
project: Ladruno
status: design → implementing
priority: high
owner: nmora
amends: 40_ladruno_performance_adr
tags:
  - performance
  - solver
  - umfpack
  - factorization-reuse
  - rank-8
  - rank-10
  - G-BYTE
---

# ADR-40c — UMFPACK Numeric-persist (factorization reuse)

> The program's **top implicit lever**, named precisely by the #512 factor/trisolve
> split ([[40b_phase0_dominance_report]], rank-8/10 addenda):
> - Lane B (3D solid J2, full/ModifiedNewton): factorization ≈ **89–90 % of `linearSolve`
>   ≈ 37–41 % of the step**; trisolve ≈ 4–5 %.
> - Lane E (DisplacementControl): **59 % of the step is UMFPACK numeric factorization,
>   two-thirds of it OUTSIDE `linearSolve`** — `DisplacementControl::update()` does
>   `setB(phat); solve()` for the reference-displacement `dUhat` *every iteration*, and the
>   current solver **refactorizes the very K the algorithm factored moments earlier**.
> - The ModifiedNewton smoke run proved it directly: `formTangent` 10 calls but `soe.factor`
>   **20** calls — tangent reuse buys **zero** factorization reuse today.
>
> **Root cause (source-verified):** `UmfpackGenLinSolver::solve()`
> (`SRC/system_of_eqn/linearSOE/umfGEN/UmfpackGenLinSolver.cpp:139`) allocates a **local**
> `void* Numeric = 0`, calls `umfpack_*_numeric` (full LU) into it, trisolves, then
> `umfpack_*_free_numeric` — **every call, unconditionally**. It is the *only* direct solver
> in OpenSees that ignores the SOE `factored` flag; every Lapack solver
> (BandGen/ProfileSPD/FullGen) has honored it for 25 years.

## Decision

Persist the UMFPACK `Numeric` object in the solver across `solve()` calls and reuse it
while the assembled matrix `A` (`Ax`) is unchanged, invalidating on exactly the events
that change `A` or its structure. **Gate = G-BYTE** (bit-identical recorder output) — this
is work-*skipping*, not algorithm replacement: the reused `Numeric` **is** the LU
factorization of the current `Ax`, and LU is a deterministic function of `Ax`, so
reuse and refactorize produce identical `X` bit-for-bit.

## The invalidation contract — reuse the 25-year-old `factored` contract, verbatim

The correctness argument is **not novel**. Every Lapack direct solver already reuses its
factorization via a `bool factored` flag on its SOE:
`BandGenLinSOE::factored` — init `false`; **reset `false` in `zeroA()` and `setSize()`**;
the solver sets it `true` after a successful factorization and, when it is already `true`,
**skips the factor and calls only the triangular solve** (`BandGenLinLapackSolver.cpp:124`,
`:155`). The load-bearing invariant is:

> **`factored == true`  ⟺  the persisted factorization is the LU of the *current* `Ax`.**

This holds because **the only writers of `Ax` are `zeroA()` and `addA()`**, and
`IncrementalIntegrator::formTangent()` *always* calls `theSOE->zeroA()` before its `addA`
assembly loop. So resetting `factored` in `zeroA()` (+ `setSize()`) captures every tangent
reassembly. UMFPACK is simply the one direct solver that never wired into this contract.

**We add nothing to the contract — we extend UMFPACK into it:**

| Event | `Ax` / structure | `factored` | Numeric action |
|---|---|---|---|
| `zeroA()` (start of every tangent reassembly) | zeroed | → `false` | (kept; freed lazily at next factor) |
| `addA()` assembly loop | filled | stays `false` | — |
| first `solve()` after assembly | valid, `factored==false` | solver sets → `true` | **factor** into persisted Numeric (free any stale first) |
| subsequent `solve()`, no reassembly (ModifiedNewton iters; DisplacementControl `dUhat`; sensitivity `setB;solve`) | unchanged | `true` | **skip factor, trisolve only** ← the win |
| `setB/addB/zeroB` then `solve()` | `A` unchanged | `true` | trisolve only |
| `setSize(Graph&)` (domainChanged / first assembly) | structure rebuilt, Symbolic re-done | → `false` | **free Numeric** (stale vs new Symbolic) |
| failed factor (singular / UMFPACK status ≠ OK) | — | stays `false` | free Numeric; return −1 |
| solver swap (`setUmfpackGenLinSolver`) | new solver, `Numeric==0` | — | belt-and-suspenders guard `Numeric!=0` forces refactor |

### Belt-and-suspenders guard
Reuse fires only if **`factored && Numeric != 0`**. A freshly constructed or swapped-in
solver has `Numeric == 0`, so even a spuriously-`true` `factored` cannot trisolve against a
null/foreign factorization — it refactorizes instead. This makes the reuse robust to any
path that sets `factored` without our Numeric being valid.

## What changes (4 files, all in `SRC/system_of_eqn/linearSOE/umfGEN/`)

1. **`UmfpackGenLinSOE.h`** — add `bool factored;` (private), `friend`-visible to the solver
   (already `friend class UmfpackGenLinSolver`). Add a public accessor is unnecessary — the
   solver is a friend, mirroring BandGen where the solver reads `theSOE->factored` directly.
2. **`UmfpackGenLinSOE.cpp`** — init `factored(false)` in both constructors; set
   `factored = false` at the end of `setSize(Graph&)` and in `zeroA()`. (Mirror BandGen
   exactly. `addA` is intentionally NOT touched — `zeroA` always precedes reassembly.)
3. **`UmfpackGenLinSolver.h`** — add `void* Numeric;` member (init 0); add `bool persist;`
   (default true) for the `-noNumericPersist` escape.
4. **`UmfpackGenLinSolver.cpp`** —
   - `solve()`: reuse-or-factor logic above; the `soe.factor` profiler scope now wraps only
     the *actual* factor call (skipped calls simply don't enter it — the measurement signal),
     `soe.trisolve` always runs. On success set `theSOE->factored = true`. On any error free
     Numeric + leave `factored == false`. If `persist == false`, free Numeric after trisolve
     and set `factored = false` (⇒ every solve refactorizes = today's behavior, bit-for-bit).
   - `setSize()`: free the persisted Numeric before/after rebuilding Symbolic (structure
     changed ⇒ Numeric stale).
   - destructor: free the persisted Numeric.
   - `OPS_UmfpackGenLinSolver()`: parse `-noNumericPersist` → `persist = false`.
   - Both index paths (`useLongIndices` DLONG + default DI) get the identical treatment; the
     persisted `Numeric` is index-mode-homogeneous (a solver instance never switches mode).

## Why it is bit-identical at EVERY solve site (the G-BYTE argument)

- **Same matrix ⇒ same factorization.** Single-threaded UMFPACK numeric factorization is a
  deterministic function of `(Symbolic, Ax, Control)`. The persisted Numeric equals what a
  refactorization of the identical `Ax` would produce — the reused-vs-refactored `X` are
  bit-for-bit equal. (BLAS-threading nondeterminism, if any, sits inside the *single* first
  factorization and is shared by both paths; it does not distinguish reuse from refactor.)
- **`dUhat` at `newStep`/`update`.** These solve against whatever `Ax`/Numeric is current
  (possibly last step's final tangent). The *current* code refactorizes that same `Ax`; the
  invariant guarantees the persisted Numeric is the factorization of exactly that `Ax`. No
  change in *which* matrix is used — only the recompute is elided.
- **The whole claim rests on one lemma:** no `Ax` mutation ever leaves `factored == true`.
  `Ax` is written only by `zeroA`/`addA`; `zeroA` resets the flag and always precedes an
  `addA` reassembly. Constraint handlers, integrators, elements all funnel matrix
  contributions through `addA`. ∎

## Serial-only — no MPI surface
`UmfpackGenLinSOE`/`Solver` is a **serial** SOE (the parallel implicit path is MUMPS,
[[project_phase0_dominance_measured]] MUMPS lane). No `sendSelf`/collective concerns like
P-NEW-1 had; `sendSelf/recvSelf` are already no-ops. `wipe`/`clearAll` destroy the
instance (destructor frees Numeric) — Numeric is a per-instance member, not static, so no
process-sticky wart.

## Memory footprint
Persisting Numeric holds the LU factors resident between solves. Peak memory is **unchanged**
(the non-persist path allocates the identical factors transiently inside every `solve()`);
steady-state resident memory rises by one factorization's worth of fill. For memory-bound
runs `-noNumericPersist` restores the free-every-solve behavior. Worth a one-line note in the
solver docstring; not a default concern (a solver that can't hold one factorization can't
factor at all).

## Validation — Zone-A gate tests (bit-identity vs `-noNumericPersist`)
`tests/test_umfpack_numeric_persist.py` (mirrors `test_cdl_mass_cache.py`), each running the
identical sequence twice — persist default-ON vs `-noNumericPersist` — comparing the full
recorded trajectory with **exact float equality**:

- **NP-1 Newton bit-identity** — full-Newton static push (A reassembled every iter ⇒ refactor
  every iter): reuse never fires, must still be identical (proves no regression on the
  no-reuse path).
- **NP-2 ModifiedNewton reuse bit-identity** — the primary reuse lane (tangent formed once,
  iters trisolve-only): identical trajectory + a convergence-count check that the run really
  iterated (reuse actually exercised).
- **NP-3 DisplacementControl bit-identity** — the `dUhat` intra-iteration reuse lane (lane E):
  identical load-displacement path incl. the `newStep`/`update` reference-displacement solves.
- **NP-4 domainChanged invalidation** — mid-run `remove element` ⇒ `setSize` ⇒ Numeric freed +
  Symbolic rebuilt: identical to `-noNumericPersist`, and the response actually changes across
  the removal (invalidation truly happened, not a stale-factor solve).
- **NP-5 failed-step revert + retry** — a step driven to non-convergence, then
  `revertToLastStep` + retry: identical, proving the error path leaves `factored==false` /
  Numeric freed so the retry refactorizes.

**Perf gate:** lane-B (ModifiedNewton) + lane-E (DisplacementControl) A/B on the built
worktree via the profiler `soe.factor` call-count drop (the direct reuse signal) + median-of-≥5
interleaved real walls (PHASE0_RECIPE machine-noise protocol). Expected: lane E ≈ −35 % (dUhat
→ trisolve), lane B ModifiedNewton material reuse; **full Newton unchanged** (no reuse to give).

## Adversarial gate (before merge — same ceremony as P-NEW-1/2)
Opus panel, distinct lenses: (a) invariant completeness — enumerate every `Ax` writer and
prove none leaves `factored==true`; (b) bit-identity — is UMFPACK reuse provably identical to
refactor, incl. index-mode and Control interplay; (c) lifecycle/edge — setSize/domainChanged
free ordering, failed-factor + revert, solver swap, `-noNumericPersist` parity, the
`saveSparseA`/`getSparseA`/`setX` read-only paths. Ship only on unanimous defensible-default-ON,
else default OFF like `-commitSolveState`.

### Gate outcome — UNANIMOUS DEFENSIBLE DEFAULT-ON (2026-07-07, 3 Opus reviewers)
- **(a) Invariant completeness → DEFENSIBLE DEFAULT-ON.** Grep-verified the ONLY writers of
  `theSOE->Ax` are `zeroA()`, `addA()`, `setSize()` (all other `Ax` tokens are unrelated
  locals; `getSparseA`/`saveSparseA` are read-only; `setX` writes `X`). Every `addA`-bearing
  reassembly is provably bracketed by a preceding `zeroA` that clears `factored`:
  `IncrementalIntegrator::formTangent` (zeroA:100 unconditional, flag-independent),
  `TransientIntegrator`, `PFEMIntegrator`, `KRAlphaExplicit(_TP)`, the **Staged\*** lonely-DOF
  identity patch (rides inside the base `formTangent`'s window), and
  **LadrunoArcLength/FictitiousMass** diagonal poke (`zeroFirst` both callers preceded by a
  zeroA). `factored` set `true` ONLY by the solver post-factor. `statusFlag` CURRENT/INITIAL/
  NO_TANGENT all route through the unconditional zeroA ⇒ no cross-tangent stale reuse.
- **(b) Bit-identity → G-BYTE HOLDS.** `umfpack_*_solve(UMFPACK_A,…)` (default IRSTEP=2
  iterative refinement) reads the SAME current `Ax` in reuse and refactor paths; `Control`
  set once in `setSize`, only read thereafter; stale `Info[]` is diagnostic-only (solve
  re-clears its slots; the singular/irstep decision reads `Numeric->rcond/nnzpiv`, not `Info`).
  The reused Numeric IS one deterministic realization of `factorize(current Ax)`. Sole caveat
  (threaded-BLAS run-to-run nondeterminism) afflicts the `-noNumericPersist` baseline equally —
  reuse adds no new nondeterminism; the persist-vs-noPersist gate test is itself the detector.
- **(c) Lifecycle → SAFE.** setSize frees Numeric BEFORE freeing/rebuilding Symbolic (UMFPACK-
  correct: Numeric references its parent Symbolic); failed factor AND failed solve both
  `freeNumeric()` + `factored=false` ⇒ retry refactorizes; `needFactor` frees before overwrite
  (no leak), `freeNumeric` nulls after free (no double-free); dtor frees Numeric before Symbolic;
  `wipe` destroys the per-instance member (not static); `-noNumericPersist` reproduces the legacy
  free-every-solve lifecycle byte-for-byte.
- **Hardening applied from the panel (non-blocking notes):** hoisted `factored=false` to the
  TOP of `UmfpackGenLinSOE::setSize` (covers the early-error returns); added an explicit
  `factored=false` in `setUmfpackGenLinSolver` (solver swap) so the invariant is locally
  self-evident, not only guard-protected.

## Measured results (2026-07-07, this worktree build, idle machine)
- **Bit-identity (Zone-A gate `test_umfpack_numeric_persist.py`): 5/5 PASS** — NP-1 Newton,
  NP-2 ModifiedNewton, NP-3 DisplacementControl, NP-4 domainChanged, NP-5 revert-retry all
  exact-float-identical vs `-noNumericPersist`.
- **Wall A/B (median-of-5 interleaved, 10³ LadrunoBrick+LadrunoJ2, UmfPack, profiler OFF):**
  | algorithm | persist | -noNumericPersist | delta |
  |---|---|---|---|
  | Newton (control) | 1.864 s | 1.842 s | **+1.2 %** (no reuse possible ⇒ no regression, no spurious reuse) |
  | **ModifiedNewton** | 1.247 s | 2.953 s | **−57.8 %** (tangent formed once/step, iters trisolve-only) |
  | **DisplacementControl** | 2.656 s | 3.532 s | **−24.8 %** (dUhat → trisolve, lane-E) |
  The `soe.factor` call-count proof was inconclusive (in-process profiler accumulation across
  configs double-books counts); the profiler-OFF wall A/B is the clean proof and is decisive:
  reuse fires exactly where the mechanism predicts (ModifiedNewton, DisplacementControl) and
  nowhere it shouldn't (full Newton unchanged).

## Ledger
- `LEDGER_vanilla_files.md`: `UmfpackGenLinSOE.{h,cpp}` + `UmfpackGenLinSolver.{h,cpp}` —
  extend UMFPACK into the standard `factored` factorization-reuse contract (rank 8/10). No new
  class tag. All edits marked `// Ladruno (ADR-40 rank 8/10)`.
</content>
</invoke>
