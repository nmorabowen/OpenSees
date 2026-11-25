# WP-258 — thread-safe soil materials (DruckerPrager, ManzariDafalias / LadrunoSANISAND, LadrunoHS / HSSmall)

Program #1075 · board #937 · ADR 232 wall W14 · follows WP-253 (#1147) · absorbs WP-239 (#1129).

**Why (owner's decision, 2026-11-25).** WP-253 threads MPM mode F's transfer loops, but a
material's calls run serially unless the material answers `ladrunoThreadSafeUpdate()` and
`ladrunoThreadSafeForm()` (the WP-107 / WP-192 contract, `Material.h`, `Element.h` items 1–5).
No soil model answers it (W14). So the soil rows are capped by their material share: DruckerPrager
42–65 % of a step → S(8) 1.3–1.7 (2D DP missed WP-253's T5 1.4, owner-accepted); SANISAND
85–95 % → S(8) 1.0–1.1. This WP removes that cap, for the FE element loops (WP-107 / WP-192) and
MPM mode F (WP-253), without changing a single constitutive result.

**Out of scope:** MPI; new materials; NorSand (#896); any change to constitutive results.

## 0. WP-239 (#1129) — absorbed

WP-239 (draft #1129, stale since 2026-11-19) audited the materials that already opt in and added
the `thread-static` quirk lint. It also edits `LadrunoSANISAND.cpp` (six `thread-ok` waivers that
say "never threaded"). Two branches must not edit the same file, and this WP retires exactly those
waivers, so WP-239 is merged into this branch as-is (its three commits keep their authorship) and
#1129 is closed as superseded here. Its own scope (ASDSteel1D Newton history and IMPL-EX
accumulator, the shell-section `getType` latch, the lint, its test) ships with this PR unchanged.

## 1. Pre-registered speedup targets (committed before any code change, 2026-11-25)

Source of the bounds: WP-253 `253_mpm_performance.md` ("The #1097 target", Results, "What soil
models stay serial means"), binary of record `3d2cb6330`, nmb (i7-13700KF, 8 P-cores, affinity
`0xFFFF`), MKL 2026.1. WP-253 measured, at 1 thread, the material share of a mode-F step:
DP 49.6 % (3D `gimp` 10⁵), 61–65 % (2D); SANISAND 85–95 %. The non-material part of the step
(the transfer loops) scales by p ≈ 5.4 at 8 threads in 2D (from 2D elastic S(8) = 3.5 with its
serial 12.5 % material share) and ≈ 5.3 in 3D (3D elastic S(8) = 5.0, material threaded).

**Amdahl, once the material threads:** with material share m, material scaling p_m and transfer
scaling p_t ≈ 5.3, S(8) = 1 / ((1 − m)/p_t + m/p_m). The material update is compute-bound but
allocates (WP-239: 30–50 Vector/Matrix temporaries per substep, #933), so p_m is not assumed to
reach 8: with p_m = 4 the bounds are DP 3D 4.6, DP 2D 4.4, SANISAND 4.1; with p_m = 3, 3.9 / 3.5 /
3.2. The targets sit at the p_m ≈ 3 end:

| # | Row (mode F, nmb, WP-253 harness `Ladruno_scripts/mpm_bench/`) | Before (WP-253) | Target |
|---|---|---|---|
| U1 | DP 3D, `gimp` and `bspline2`, 10⁵, 8 threads | 1.54 / 1.71 | **S(8) ≥ 3.0** |
| U2 | DP 2D, `gimp` and `bspline2`, 10⁵, 8 threads | 1.27 / 1.31 | **S(8) ≥ 2.5** (lifts WP-253's T5 miss) |
| U3 | SANISAND (IntScheme 2) 2D and 3D, `bspline2`, 10⁵, 8 threads | 1.02–1.11 | **S(8) ≥ 2.5** |
| U4 | DP and SANISAND rows above, 4 threads | — | S(4) ≥ 2.0 |
| U5 | every DP / SANISAND row at 10⁴ | — | reported only (WP-253 found 10⁴ rows noisy: 23–45 % spread); S(n) ≥ 0.95 |
| U6 | every row at 1 thread, after vs WP-253's baseline | — | G(1) within noise of 1.0, and bit-identical results |
| U7 | `regress.py` against WP-253's `baseline_nmb.json` | — | no row regresses beyond noise |

Read on medians with WP-253 §4's noise rule (3 repeats, a 5-repeat rerun above 10 % spread; a
median under target by less than its spread is "within noise of the target", not a pass). A miss
goes to the owner. HS is not in WP-253's matrix; it gets the determinism gate (§3) and one
reported mode-F row, no target.

## 2. Audit — shared mutable state on the threaded call graphs

The worker-thread methods are the WP-107 / WP-192 contract's: `setTrialStrain` (Update) and
`getStress`, `getTangent`, `getInitialTangent`, `commitState`, `revertToLastCommit` (Form). MPM
mode F (WP-253) calls `setTrialStrain`, `getStress`, `commitState` and `getInitialTangent` on
workers. Everything else (parsers, `getCopy`, `sendSelf`/`recvSelf`, `setParameter` /
`updateParameter`, `updateMaterialStage`, the MPM stage walk, `Print`, `getResponse`) runs on the
main thread between loops. Shared infrastructure already per-thread inside a threaded region
(WP-192): `Matrix::Solve` / `Invert` work areas (`Matrix.cpp:55-77`, `:425`, `:533`, `:642`);
every loop (`Domain::update`, the integrators' tangent / residual loops, `Domain::commit`, MPM mode
F) runs through `ladruno_parallelFor`, which sets `ladruno_tlInRegion`.

### 2.1 DruckerPrager (vanilla UW, `DruckerPrager3D`, `DruckerPragerPlaneStrain`)

| Item | Where | Kind | Reached from | Fix | Serial bits |
|---|---|---|---|---|---|
| `strain`, `stress`, `tangent` | `DruckerPragerPlaneStrain.h:75-77`, `.cpp:29-31` | class-static return buffers | `getStrain`, `getStress`, `getTangent`, `getInitialTangent` (PlaneStrain) | `static thread_local` | unchanged |
| `ladNormCepMsgCount` + `opserr "NormCep"` | `DruckerPrager.cpp:798-803` | function static, non-atomic RMW, `opserr` | `plastic_integrator` (every apex point after ADR-95) | `std::atomic` claim of the same 10-line budget, `LadrunoOutputGuard` | unchanged (same lines, same order) |
| `opserr "Jact = … count = …"` | `DruckerPrager.cpp:638-641` | `opserr` | the corner-loop bailout (`count > 3`) | `LadrunoOutputGuard` | unchanged |
| `numDruckerPragerMaterials` | `:56` | file static | parser | none (main thread) | — |
| `static Vector data(45)` | `:1117`, `:1182` | function static | `sendSelf` / `recvSelf` | none (main thread) | — |
| `mElastFlag`, `mFlag` | `DruckerPrager.h:133-134` | **per instance** | read in `plastic_integrator`, `updateElasticParam` | none; `materialState` is written by `updateParameter` (main thread) | — |
| `g.Invert(g_contra)` | `:549`, `:584` | Matrix work area | the corner return | already per-thread in a region (WP-192) | — |
| locals `Vector`/`Matrix` | `plastic_integrator` | heap temporaries | every call | none (re-entrant; allocator contention only) | — |
| `ops_InitialStateAnalysis` | `:322` | global read | `revertToStart` (main thread) | none | — |

No `throw`/`exit` on the worker paths (`exit` only in the base's "subclass responsibility" stubs,
which the leaves override). **Verdict: safe for every option once the three rows above are
fixed; both leaves opt in.**

### 2.2 ManzariDafalias / LadrunoSANISAND

Options: IntScheme 1 (ModifiedEuler), 2 (BackwardEuler_CPPM), 129 (SAS-ME,
`LadrunoSANISANDSasME.cpp`), `-implex` (ADR-92); ThreeDimensional and PlaneStrain. "All" = every
one of these. Line numbers are the pre-change ones.

| Item | Where | Kind | Reached from | Fix | Serial bits |
|---|---|---|---|---|---|
| **`opserr` → `PythonStream::err_out` → `PySys_FormatStderr` on a worker** | `PythonStream.h:64-75` | a CPython call from a thread with no Python thread state (and one shared `std::string msg`) | every worker print below, under openseespy | `err_out` writes a local string to C `stderr` when `!PyGILState_Check()`; the GIL-holding thread keeps the Python path. **Root cause of the WP-107 segfault** (§2.4) | unchanged |
| ME low-p clamp print | `ManzariDafalias.cpp:1781-1792` | `opserr` (atomic budget since WP-107) | update: 1; 2 (ME fallback); `-implex` (companion); 129 when SAS-ME is not allowed | `LadrunoOutputGuard` | unchanged |
| ME substep-cap print | `:1910-1936` | `opserr` (atomic budget) | same, with `-maxSubsteps` | `LadrunoOutputGuard` | unchanged |
| RK45 scratch, `do_once`, RK45 clamp budget | `:2396-2461` | function statics, plain budget, `opserr` | IntScheme 45 only | **not admitted** (scheme 45 refuses); lint waivers | — |
| `mElastFlag` | `ManzariDafalias.h:557`, `.cpp:119` | class static (stage flag) | **read** on workers (`GetElasticModuli`, `LadrunoSANISAND.cpp:3304`, `:3398`); **written** by `updateMaterialStage` / `updateParameter`, SANISAND replay, the MPM stage walk — all main thread | none; never written from a worker (ADR 232 232k: it stays class-wide) | — |
| `mI1`, `mIIco`, `mIIcon`, `mIImix`, `mIIvol`, `mIIdevCon`, `mIIdevMix`, `mIIdevCo` | `.cpp:126-134` | class statics | constant after static initialisation | none | — |
| `mSigma_M`, `mEpsilon_M`, `rSigma`, `mTangent`, `mTangent_init` | `ManzariDafaliasPlaneStrain.cpp:27-31`; `LadrunoSANISANDPlaneStrain.cpp:53-57` | class-static return buffers | `getStress`, `getStrain`, `getTangent`, `getInitialTangent` (PlaneStrain) | `static thread_local` | unchanged |
| `mSigma_M`, `mEpsilon_M` | `ManzariDafalias3D.cpp:27-28`; `LadrunoSANISAND3D.cpp:50-51` | class-static return buffers | `getStress`, `getStrain` (3D; `getTangent` returns per-instance `mCep`) | `static thread_local` | unchanged |
| `LadrunoImplexGlobals` refusal / guard counters (`long ++`) | `LadrunoSANISAND.cpp:2770-3000`; bumped at `:3383`, `:3521`, `:3684`, `:3700` (the reversal-noise guard, **also without `-implex`**), `:4019`, `:4096`, `:4199`, `:4368`, `:4399`, `:4516`, `:4547`, `:4593`, `:4631`, `:4765`, `:4963`, `:5130` | process-wide singleton, non-atomic RMW | trial and commit | `std::atomic<long>` (integer, order-free: exact) | unchanged |
| `LadrunoImplexGlobals` error ledger (`noteCommitRound` first-committer marker, `maxError`, `sumError`, `count`) | `:2778-2824`, called at `:4841-4843` | process-wide singleton, FP sum, a round boundary that assumes the marker commits FIRST | `commitState` with `-implex` | per-thread partials keyed by the commit loop's id (`ladruno_currentLoopId()`, new in `LadrunoThreads`), folded under a mutex when the thread moves to another loop and before every read; a partial older than the latest marker loop is dropped; the globals reset once per marker round. Max exact; the average may differ in its last bits under threads (a diagnostic — the ASDConcrete WP-192 precedent). Serial paths fold first, then run the original code | unchanged (flush is a no-op when no partial was written) |
| WP-239's six `thread-ok` "never threaded" statics: post-latch warn, round-off α_in warn, D2 warn, `-implexControl` rung `lastDt` / `rungLines` / count, companion warn | `:3348`, `:3787`, `:4021`, `:4412-4422`, `:4633` | function statics (plain RMW) + `opserr` | trial (latch: all; round-off: `-flipAlphaIn vanilla`; the rest `-implex`); companion at commit | `std::atomic` (lint-exempt; the waivers are removed — their reason became false) and `LadrunoOutputGuard` | unchanged |
| per-instance commit-refusal print | `:4918` | `opserr` | `commitState` (all) | `LadrunoOutputGuard` | unchanged |
| SAS-ME refusal print (per-instance warn-once) | `LadrunoSANISANDSasME.cpp:1459` | `opserr` | update (129) | `LadrunoOutputGuard` | unchanged |
| `OPS_PROFILE_SCOPE` | `LadrunoSANISANDSasME.cpp:685-1475` | profiler | 129 | none (per-thread registry, `Profiler.h:347-383`) | — |
| `ladrunoNoteCommitRefusal` | `LadrunoMaterialStatus.h:123-159` | atomic + `thread_local` depth | commit | none | — |
| `ops_Dt` | `:3467`, `:3682`, `:4592`, `:4950` | global read | trial / commit | none (written once before the loops, `Domain.cpp:2512`) | — |
| `Matrix::Solve` / `Invert` | the live CPPM path inverts on the stack (`ladrunoInvertLocal`, `ManzariDafalias.cpp:81-111`) | — | 2 | none (and every Matrix overload has the WP-192 per-thread work area inside a region) | — |
| substep trace `push_back` | `ManzariDafalias.cpp:1684` | per instance (`getCopy` nulls it) | replay only | none | — |
| `getPStrain` `static Vector result` | `:6216` | function static | `setResponse` / `getResponse` only | lint waiver (main thread) | — |
| `numManzariDafaliasMaterials`, constructor latches | `:137`, `:272`, `:371` | file / function statics | parser, constructors | none (main thread) | — |

`getCopy` is `*clone = *this` (a deep copy); nothing on a worker path writes another instance or
the prototype. `ManzariDafaliasRO` and its wrappers were not audited: they now refuse explicitly,
because the base admits IntScheme 1/2.

**Verdict:** IntScheme 1, 2, 129 (129 only where SAS-ME is allowed, i.e. `LadrunoSANISAND`), with
or without `-implex`, 3D and PlaneStrain: **admitted**. IntScheme 45: refused (function-static
scratch). IntSchemes 0 and 3–9 share ModifiedEuler's helpers and hold no statics, but were not in
the determinism gate: refused, pending a gate.

### 2.4 The WP-107 segfault — mechanism

WP-107 found that ManzariDafalias IntScheme 1 segfaulted (`0xC0000005`) on an OpenMP worker as
soon as the plastic branch ran in volume, **even with the whole `update()` serialised under
`omp critical`**, and never on the elastic branch. So it was not a data race. Its logs
(`Ladruno_files/testbed/perf/wp107/run1/`) show what differs: the serial twin prints the
ModifiedEuler substep-cap warning (`ManzariDafalias.cpp:1912`) and the threaded run stops dead
after the `THREADED` line. The elastic branch never prints.

Under openseespy `opserr` is a `PythonStream`, whose `err_out` calls `PySys_FormatStderr`. On a
worker thread, which has no Python thread state, CPython dereferences null. Serialising with
`omp critical` does not give a worker a thread state, which is why it changed nothing.
`LadrunoOutputGuard` is only a mutex, so it does not help either. This was WP-131's open item S1
(`131_sanisand_threaded_inventory.md` §0), never tested. Taking the GIL from the worker would
deadlock, because the main thread holds it for the whole region. The fix is in `err_out`: a
thread that does not hold the GIL writes to C `stderr`. The same defect sat behind every WP-192
guarded print (ForceBeamColumn3d, LadrunoDispBeamColumn) under openseespy; Tcl (`OpenSees.exe`,
`StandardStream`) was never affected. The verification is in §5.

### 2.3 LadrunoHS / HSSmall

HSSmall is the same class with `use_small_strain` (`OPS_LadrunoHSSmall`, `LadrunoHS.cpp:326`).
The kernel `LadrunoHSKernel.h` is pure: `inline` functions on stack arrays and a
`const HSParams&`; its only statics are `static constexpr` (`:58`, `:168`, `:768-769`, `:795`,
`:1707`); no `Matrix`/`Vector`, `opserr`, `throw` or `exit`; it solves its own 2×2 systems
(`solve_2x4_neg`, `:1006`) and its eigenvalues by Jacobi (`:318`).

| Item | Where | Kind | Reached from | Fix | Serial bits |
|---|---|---|---|---|---|
| tag stage map, ref-count map, mutex | `LadrunoHS.cpp:60-73` | immortal registry (function-static `new`) | **written**: parser, constructor (`getCopy`), destructor (`wipe`), `recvSelf`, `updateParameter` (incl. the MPM stage walk) — all main thread. **Read**: `syncStageFromRegistry` from `integrate3D`, `getStress`, `getStressZZ`, `getStrain`, `getTangent`, `getInitialTangent`, `commitState`, `revertToLastCommit` | atomic **generation** counter bumped by every writer under the mutex; `syncStageFromRegistry` returns when the generation equals the instance's last-seen value, else takes the lock as before | unchanged: the lazy transition runs at exactly the same calls (invariant in the code comment) |
| the same lock | `syncStageFromRegistry` | one process-wide `std::mutex` taken 3–4× per point per iteration (up to 25× in condensed modes) | every worker method | removed from the hot path by the generation | — (a scalability defect, WP-192 trap 4, not a race: every access was locked) |
| `ladrunoNoteCommitRefusal()` | `:833` → `LadrunoMaterialStatus.h:155-159` | process-wide `std::atomic<int>` | `commitState` | none (integer, order-free); WP-253 adds the per-thread tally | — |
| `stressOut_`, `strainOut_`, `tangentOut_`, `initTangentOut_` | `LadrunoHS.h:142-145` | per-instance return buffers | getters | none | — |
| IMPL-EX state | `:635-679`, `:837-841` | per instance | `-implex` | none (no global accumulator, unlike SANISAND/ASDConcrete) | — |
| `-sig0`/`-sig03`, `-sigmat` | `sig0_comp_`, `buildInitialTangent6x6` (`:467-498`) | per instance / pure | FE (MPM refuses `-sig0`) | none | — |
| `static Vector data` | `:990`, `:1063` | function static | `sendSelf` / `recvSelf` | none (main thread) | — |
| `opserr` | parser, send/recv, `updateParameter` only | — | never from a worker method | none | — |
| `ops_Dt`, `ops_TheActiveElement`, Domain | not referenced | — | — | — | — |

**Verdict: safe for every option (HS, `-smallStrain`, `-implex`, `-sigmat`, `-sig0`/`-sig03`,
every dimension mode); opts in unconditionally.** The registry's three immortal objects carry
`thread-ok` waivers (written only on the main thread, under the mutex).

## 3. Determinism gate

Every threaded result bit-identical to 1 thread, at 2, 4 and 8 threads, for FE (WP-192 element
loops: update, tangent, residual, commit) and MPM mode F (WP-253 loops), including a
shuffled-particle-order case; every option admitted in §4. With 1 thread every existing record
stays bit-identical. Mutation evidence: each determinism test goes red when a shared static is
re-introduced on purpose.

## 4. Opt-in matrix

| Material | Option | `ladrunoThreadSafeUpdate` / `Form` | Reason |
|---|---|---|---|
| `DruckerPrager3D`, `DruckerPragerPlaneStrain` | every option (`materialState` 0/1/2) | **true** | §2.1 |
| `ManzariDafalias` (vanilla), 3D / PlaneStrain | IntScheme 1, 2 | **true** | §2.2 |
| `ManzariDafalias` | IntScheme 0, 3–9 | false | not in the determinism gate |
| `ManzariDafalias` | IntScheme 45 (RungeKutta45) | false | function-static scratch (`:2414-2421`) |
| `ManzariDafaliasRO` and its 3D / PlaneStrain wrappers | any | false | not audited |
| `LadrunoSANISAND`, 3D / PlaneStrain | IntScheme 1, 2, 129; with or without `-implex` | **true** | §2.2 |
| `LadrunoSANISAND` | other IntSchemes | false | as the base |
| `LadrunoHS` / HSSmall | every option, every dimension mode | **true** | §2.3 |

Elements still answer for themselves: `stdBrick` threads every loop, `LadrunoQuad` only the update
loop, and MPM mode F threads a particle's material calls when its material answers both.

## 5. Results

*(after the build of record)*
