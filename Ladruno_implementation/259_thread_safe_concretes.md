# WP-259 — thread-safe Ladruno concretes (`LadrunoConcrete3D`, `LadrunoRCConcrete`)

Board #937 · stacked on WP-258 (#1153, itself stacked on WP-253 #1147) · follows the WP-258 pattern
(`258_thread_safe_soil_materials.md`).

**Why.** Neither Ladruno concrete answers `Material::ladrunoThreadSafeUpdate()` /
`ladrunoThreadSafeForm()` (the WP-107 / WP-192 contract, `Material.h`, `Element.h` items 1–5), so
every host that holds one (`stdBrick`, `LadrunoQuad`, `ASDShellQ4` + `LayeredShell`, MPM mode F)
falls back to its serial loop. `ASDConcrete3D` already threads (WP-192). This WP makes both
Ladruno concretes answer true wherever the audit below proves it, without changing a single
constitutive result: with 1 thread every record stays bit-identical, and the threaded result is
bit-identical to 1 thread.

**Out of scope:** `LadrunoRCFiniteStrain` (a separate `FiniteStrainNDMaterial` class, classTag
33018; it does not inherit either material and keeps the default false); new elements opting in
(`LadrunoNDFiberSection3d`, `LadrunoBrick`, `LadrunoSolidShell` do not thread today); MPI; any
constitutive change.

## 1. Pre-registered targets

None. The measurements (§5) are reported, as the brief allows. The only hard gate is determinism
(§3).

## 2. Audit — shared mutable state on the threaded call graphs

The worker-thread methods are the contract's: `setTrialStrain` / `setTrialStrainIncr` (Update) and
`getStress`, `getTangent`, `getInitialTangent`, `commitState`, `revertToLastCommit` (Form). MPM
mode F calls `setTrialStrain`, `getStress`, `commitState` and `getInitialTangent` on workers. The
parsers, `getCopy`, `sendSelf` / `recvSelf`, `Print`, `setResponse` / `getResponse`,
`revertToStart` and `NDMaterial::wipe`'s budget reset run on the main thread between loops.
Shared infrastructure already per-thread inside a region (WP-107 / WP-192 / WP-258):
`ops_TheActiveElement` (`thread_local`, `Element.cpp:61`; set per element by `Domain`'s threaded
update loop `Domain.cpp:2776`, by the threaded tangent / residual loops
`LadrunoThreadedAssembly.cpp:67`, `:239`, and zeroed per step by MPM mode F, W9,
`LadrunoMPMPattern.cpp:131`, the same value its serial path sees); `opserr` under openseespy
(`PythonStream::err_out` writes to C `stderr` off the GIL, WP-258 §2.4); `ladrunoNoteCommitRefusal`
(`LadrunoMaterialStatus.h`, atomic).

Line numbers are the pre-change ones.

### 2.1 `LadrunoConcrete3D` (`LadrunoConcrete3D.{h,cpp}`, `LadrunoConcrete3DKernel.h`)

The kernel is header-only and OpenSees-free (`#include <cmath>` only): `inline` functions on stack
arrays and a `const Params&`. Its namespace-scope `static const` scalars and the function-local
`static const double W6[6]` (`:2294`, `:2344`, `:2426`) are constants; it has no `Matrix` /
`Vector`, `opserr`, `printf`, `throw`, `exit` or `getenv`; it solves its own 3×3 / 4×4 systems
(`solve3`). Every dimension mode, every option (`-implex`, `-eta`, `-autoRegularization`,
`-tensionLaw`, `-flowPotential`, `-compressionDrive`, `-ctTemper`, `-tcTemper`, `-subIncr`,
`-deadThreshold`, `-hoop`, `-epsFc` / `-gcLegacy`, `-verbose`) is a branch on per-instance data or
on `Params` built on the stack in `integrate()` (`:448-475`).

| Item | Where | Kind | Reached from | Fix | Serial bits |
|---|---|---|---|---|---|
| `s_returnMapWarnPrinted` + two `opserr` lines | `LadrunoConcrete3D.cpp:57`, `:529-536` | file static `long`, non-atomic check-then-increment, `opserr` | `integrate()` ← `setTrialStrain` (every final return-map failure, every option, every view) | `std::atomic<long>` claim with `fetch_add` (the same 20-line budget: claim *k* prints when *k* < 20, the suppression line when *k* + 1 = 20 — exactly the serial sequence, and no overshoot under threads); `LadrunoOutputGuard`; the `wipe` reset becomes a `store(0)` | unchanged (same lines, same order) |
| `-verbose` final-failure line | `:530` | `opserr` (no budget) | same | `LadrunoOutputGuard` | unchanged |
| σ₂₂ condensation "pivot vanished" / "did not converge" | `:627`, `:633` | `opserr` (unthrottled) | `setTrialStrain`, condensed views (PlaneStress, PlateFiber) | `LadrunoOutputGuard` | unchanged |
| Gc/lch out-of-range warning | `:872` | `opserr`; per-instance once-latch `gcWarned` | `compressiveEpsFc()` ← `integrate()` | `LadrunoOutputGuard` (the latch is per instance) | unchanged |
| `-verbose` commit summary | `:780` | `opserr` | `commitState` | `LadrunoOutputGuard` | unchanged |
| `gcLch`, `gcEpsFc`, `gcWarned`, `gcTabReady`, `gcEfc[]`, `gcG[]` | `.h:161-164` | per-instance cache (written in `compressiveEpsFc` / `ensureGcTable`) | `integrate()` | none (each instance belongs to one element / particle; `getCopy` builds the table on the prototype, main thread, before copying it) | — |
| `nSub`, `nFail`, `nSubStep`, `nFailStep` | `.h:160` | per-instance diagnostics | `integrate()`, `commitState` | none | — |
| `kernelCallCounters()` | `Kernel.h:217` | `static thread_local` | every return map | none (already per thread; read only by the g++ kernel self-check) | — |
| `ops_TheActiveElement->getCharacteristicLength()` | `:463-466` | `thread_local` global read | `integrate()` with `-autoRegularization` | none: set per element by every threaded loop before the element's calls (above); MPM zeroes it per step on master and worker alike, so `lch` falls back to `-lch` in both | — |
| `ops_Dt` | `:491` | global read | `integrate()` (`-implex`, `-eta`) | none (written once per step before the loops, `Domain.cpp:2512`) | — |
| `ladrunoNoteCommitRefusal()` | `:767` | process-wide `std::atomic<int>` | `commitState` (refused trial) | none (integer, order-free) | — |
| `Vector ne(ncomp)` | `setTrialStrainIncr`, `:650` | heap temporary | `setTrialStrainIncr` | none (re-entrant) | — |
| `static Vector data(LC3D_NDATA)` | `:927`, `:968` | function static | `sendSelf` / `recvSelf` | none (main thread) | — |
| `throw` / `exit` / `Matrix::Solve` / `Invert` | — | — | not on any worker path | — | — |

**Verdict:** safe for every option once the five `opserr` rows are guarded and the budget is
atomic. Opt-in per view (§4): the views a threaded host reaches and the gate runs.

### 2.2 `LadrunoRCConcrete` (`LadrunoRCConcrete.{h,cpp}`, `LadrunoRCKernel.h`)

The kernel is header-only (`#include <math.h>` only): its function-local statics are
`static const` (`:1353`, `:1399`, `:1460`, `:1558`, `:1608`); no namespace-scope mutable
variable, `opserr`, `throw` or `Matrix`. `Params` holds the backbones as fixed arrays
(`MAXPTS = 32`), so a copy owns its own; `regularize()` mutates the instance's `P.ht` / `P.hc`.
Options: `-implex` (+ `-implexAlpha`, `-implexControl`), `-autoRegularization` / `-lch`, `-beta`,
`-interlock` / `-cyclic` / `-xcrack` / `-shearRetention`, `-tensStiff*`, `-betaC`, `-crackedNu`,
`-addonPlane`, `-crackGate`, `-ties` / `-crackCheck` / `-lateralNoRoot` / `-tieHardening`
(BeamFiber), `-secant` / `-numericalTangent`, views 3D / PlaneStress / PlateFiber / BeamFiber.

| Item | Where | Kind | Reached from | Fix | Serial bits |
|---|---|---|---|---|---|
| `-autoRegularization` loud failure | `LadrunoRCConcrete.cpp:451` | `opserr`; per-instance once-latch `regWarned` | `regularizeIfNeeded()` ← `setTrialStrain` | `LadrunoOutputGuard` | unchanged |
| BeamFiber lateral non-convergence | `:490-493` | `opserr`; per-instance budget `bfWarned` | `setTrialStrain` (BeamFiber) | `LadrunoOutputGuard` | unchanged |
| σ₃₃ condensation "did not converge" | `:528` | `opserr` (unthrottled) | `setTrialStrain` (PlaneStress, PlateFiber) | `LadrunoOutputGuard` | unchanged |
| IMPL-EX error warning | `:612` | `opserr` | `commitState` (`-implexControl`) | `LadrunoOutputGuard` | unchanged |
| `implexError`, `dtime_n`, `dtime_n_commit`, `dtime_0`, `commitDone` | `.h:122-124` | per instance | `setTrialStrain`, `commitState` | none (no process-wide ledger, unlike SANISAND / ASDConcrete) | — |
| `regularizationDone`, `regLch`, `P.ht`, `P.hc` | `.h:127-129` | per instance (latched once) | `regularizeIfNeeded()` | none | — |
| `ops_TheActiveElement->getCharacteristicLength()` | `:445` | `thread_local` read | first `setTrialStrain` with `-autoRegularization` and no `-lch` | none (as 2.1) | — |
| `ops_Dt` | `:473` | global read | `setTrialStrain` (`-implex`) | none | — |
| `static bool warned` | `:684` | function static, plain RMW, `opserr` | `getCopy("BeamFiber")` with `-interlock` | `std::atomic<bool>` `exchange` (main thread today — element / section construction — but cheaper to make race-free than to argue; lint-exempt) | unchanged |
| `static Vector data(RC_DATA)` | `:728`, `:793` | function static | `sendSelf` / `recvSelf` | none (main thread) | — |
| `throw` / `exit` / `Matrix::Solve` / `Invert` | — | — | not on any worker path | — | — |

**Verdict:** safe for every option once the four `opserr` rows are guarded.

### 2.3 Hosts (not changed here; recorded so the opt-in reaches a loop)

| Host | View requested | Loops it threads when its materials answer true |
|---|---|---|
| `stdBrick` (`Brick.cpp:2295`) | ThreeDimensional | update, tangent, residual, commit |
| `LadrunoQuad` (`LadrunoQuad.cpp:826`, not EAS, not `-geom finite`) | PlaneStrain, PlaneStress | update |
| `ASDShellQ4` + `LayeredShell` (`ASDShellQ4.cpp:2448`, `LayeredShellFiberSection.cpp:850`) | PlateFiber | update, tangent, residual, commit |
| MPM mode F (WP-253) | ThreeDimensional, PlaneStrain | the particle material calls |

No threaded host requests AxiSymmetric or BeamFiber (`LadrunoNDFiberSection3d` does not opt in).

## 3. Determinism gate

`tests/test_wp259_concrete_threads.py` (zone_a): every threaded result bit-identical to 1 thread
at 2, 4 and 8 threads (`float.hex` of every nodal displacement after every step, plus the load
factor and every step's return code), shuffled element tags, and the deck must really crack /
crush (asserted from the material's own damage response). Meshes: `stdBrick` 3D, `LadrunoQuad`
PlaneStrain and PlaneStress, `ASDShellQ4` + `LayeredShell` (PlateFiber), and MPM mode F 3D / 2D.
Options: default, `-implex`, `-eta`, `-autoRegularization` (LadrunoConcrete3D); default,
`-implex`, `-autoRegularization`, `-interlock -cyclic` (LadrunoRCConcrete). Neither material keeps
a process-wide FP accumulator, so no tolerance exception is needed: everything is gated exactly.
Mutation evidence: re-introduce a shared static on a gated path, build, show red; never committed.

## 4. Opt-in matrix

| Material | View | Option | `ladrunoThreadSafeUpdate` / `Form` | Reason |
|---|---|---|---|---|
| `LadrunoConcrete3D` | ThreeDimensional, PlaneStrain, PlaneStress, PlateFiber | every option | **true** | §2.1 + gate |
| `LadrunoConcrete3D` | AxiSymmetric, BeamFiber | any | false | no threaded host requests them; not gated |
| `LadrunoRCConcrete` | ThreeDimensional, PlaneStress, PlateFiber | every option | **true** | §2.2 + gate |
| `LadrunoRCConcrete` | BeamFiber | any | false | no threaded host requests it; not gated |
| `LadrunoRCFiniteStrain` | — | — | false (unchanged) | out of scope |

## 5. Results

(filled in after the build)
