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

**Verdict:** `LadrunoSANISAND` IntScheme 1, 2 and 129, with or without `-implex`, 3D and
PlaneStrain: **admitted**. Vanilla `ManzariDafalias`: IntScheme 1 only (after the critic, §6.1). IntScheme 45: refused (function-static
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
| `ManzariDafalias` (vanilla), 3D / PlaneStrain | IntScheme 1 | **true** | §2.2 |
| `ManzariDafalias` (vanilla) | IntScheme 2 | false | no substep cap: the gate cannot run it in bounded time (§6.1 #3); LadrunoSANISAND admits it |
| `ManzariDafalias` | IntScheme 0, 3–9 | false | not in the determinism gate |
| `ManzariDafalias` | IntScheme 45 (RungeKutta45) | false | function-static scratch (`:2414-2421`) |
| `ManzariDafaliasRO` and its 3D / PlaneStrain wrappers | any | false | not audited |
| `LadrunoSANISAND`, 3D / PlaneStrain | IntScheme 1, 2, 129; with or without `-implex` | **true** | §2.2 |
| `LadrunoSANISAND` | other IntSchemes | false | as the base |
| `LadrunoHS` / HSSmall | every option, every dimension mode | **true** | §2.3 |

Elements still answer for themselves: `stdBrick` threads every loop, `LadrunoQuad` only the update
loop, and MPM mode F threads a particle's material calls when its material answers both.

## 5. Results

**Binary of record:** `efc4113f1` (WP-258 on WP-253 `faae5c51c`). `ladrunoBuild` of
`OpenSees.exe` and `opensees.pyd` were checked before every run. Later commits change tests,
scripts and documents only. Machine nmb (i7-13700KF, affinity `0xFFFF`), MKL 2026.1.0.0, Windows
11. Data: `Ladruno_scripts/mpm_bench/results/after_F_efc4113f1_wp258/` (+ `_rerun`, the §4
5-repeat reruns; + `_hs`).

### 5.1 Speedup targets (§1) — every target met

| # | Row | Before (WP-253) | After S(n), median (spread) | Verdict |
|---|---|---|---|---|
| U1 | DP 3D 10⁵, 8 threads | 1.54 `gimp` / 1.71 `bspline2` | **5.29** (5 %) / **4.97** (1 %) | pass (≥ 3.0) |
| U2 | DP 2D 10⁵, 8 threads | 1.27 / 1.31 (WP-253 T5 miss) | **5.17** (9 %) / **5.20** (9 %) | pass (≥ 2.5); WP-253's T5 miss lifted |
| U3 | SANISAND (IntScheme 2, `bspline2`) 10⁵, 8 threads | 1.08 (3D) / 1.11 (2D) | **4.47** (4 %) / **2.63** (3 %) | pass (≥ 2.5) |
| U4 | the same rows, 4 threads | — | DP 3.06–3.17; SANISAND 3D 2.82, 2D 2.44 | pass (≥ 2.0) |
| U5 | 10⁴ rows (reported) | — | DP S(8) 5.25–5.93; SANISAND 4.70 (3D), 3.23 (2D) | reported; every S(n) ≥ 0.95 |
| U6 | 1 thread, after against WP-253's baseline | — | G(1) 0.99–1.05, each "faster" or "within noise"; **particle digests bit-identical** to `baseline_F_t1` (origin `a3820d5ac` + scopes) and to WP-253's `after_F_3d2cb6330`, 12/12 rows | pass |
| U7 | `regress.py` against `baseline_nmb.json` | — | no regression (3D `gimp` elastic 8 threads +20 %, other rows within ±4 %) | pass |

T8-style determinism: in the bench, every 4- and 8-thread digest equals its 1-thread digest
(24/24 rows).

**What the remaining gap is.** 2D SANISAND reaches only 2.63 at 8 threads against 2.44 at 4. The
material share is 83–89 % of a 2D SANISAND step at one thread, so the bound is set by the
material update's own scaling. At 10⁵ its 2D call is about 3× cheaper than in 3D
(2903 vs 9798 ns per particle-step), so the per-call `Vector`/`Matrix` heap traffic (30–50
temporaries per substep, #933) is a larger fraction. The likely ceiling is allocator contention,
not a lock. No measurement here isolates it; it is #933's work.

### 5.2 HS (reported, no target)

`--mats hs` (added to the harness as a reported-only material: the default matrix stays WP-253's
81 rows). The kPa deck is the one from `test_ladruno_mpm_soil_admission.py`: stage 0 for the
warm-up, then stage 1. On `bspline2` at 10⁵, S(8) is **4.52** in 3D (0.599 → 2.71 M p-steps/s)
and **4.74** in 2D (1.07 → 5.09). The material share at one thread is 34 % (3D) and 50 % (2D).
The 8-thread digests equal 1 thread.

### 5.3 Determinism gates (Windows, `efc4113f1`)

- `tests/test_wp258_soil_threads.py` (FE): bit-identical to 1 thread at 2, 4 and 8 threads.
  - Materials: DP, vanilla ManzariDafalias IS1, LadrunoSANISAND IS1/IS2/IS129/`-implex`, and HS, `-smallStrain`, `-implex`, `-sigmat`.
  - Meshes: 2D (`LadrunoQuad`, update loop) and 3D (`stdBrick`, update / tangent / residual / commit loops).
  - Element tags are shuffled.
  - The stage flips happen inside the threaded loops.
  - Every deck is plastic: the footing secant ratio over the push is 0.54–0.92.
  - Plus a worker-print test under openseespy: DP at the apex prints `NormCep` from workers.
- `tests/test_ladruno_mpm_wp258_soil_threads.py` (mode F): every particle array (`x v V F R epsHat sigHat sig`) is float64-equal at 2, 4 and 8 threads, and `mpm debug threads` reads [n, n].
  - Soils: SANISAND IS1/IS2/IS129/`-implex` and HS (plain, `-smallStrain`, `-implex`), in 2D and 3D on `bspline2`: a staged, confined column with a rigid plate pushing down after `mpm stage 1`.
  - DP: WP-253's collapsing column (`gimp`, 2D and 3D).
  - A shuffled-particle SANISAND 3D case.
- `tests/test_ladruno_mpm_wp253_threads.py` (WP-253's gate, updated): `3d_druckerprager` now threads both loops; `3d_ladrunoj2_split` keeps the split path covered. All of the file passes.
- **Mutation evidence** (working-tree mutants, built and run, never committed):
  - M1, the `PythonStream` off-GIL path disabled: the worker-print child dies with `0xC0000005`. This reproduces WP-107's segfault mechanism.
  - M2′, `LadrunoSANISAND3D::mSigma_M` back to a plain shared static: all four SANISAND 3D FE cases go red.
  - M2, the same on vanilla `ManzariDafalias3D`, stayed green because the matrix did not yet run vanilla ManzariDafalias. That is why `manzari_is1` was added.
- **Race stress** (`Ladruno_scripts/wp258_race_stress.py`, 8 threads, 10 repeats per case, each
  checked against its 1-thread digest; `results/after_F_efc4113f1_wp258/race_stress.json`):
  0 mismatches in 110 runs.
  - FE: SANISAND IS2 3D, `-implex` 3D, IS129 2D; HS `-smallStrain` 3D; DP 2D; vanilla ManzariDafalias IS1 3D.
  - MPM: SANISAND IS2 3D shuffled, `-implex` 3D, IS1 2D; HS `-implex` 3D; DP 3D.

### 5.4 Not done here, and why

- **TSan.** GCC's libgomp is not TSan-instrumented, so OpenMP's own barriers report as races
  unless the runtime is LLVM libomp built with Archer. That is a toolchain project on esmeralda,
  not a WP step. The brief's fallback was run instead: repeated 8-thread runs compared by digest
  (`Ladruno_scripts/wp258_race_stress.py`). The class of bug that the TSan item in the WP-107
  quirk was meant to catch turned out not to be a race (§2.4).
- **Mode A soil particles stay serial.** The mode-A particle element (`LadrunoMPMPoint`) does not
  opt into the WP-107 / WP-192 loops (ADR 232 §6). `test_s7_k2_threads_honoured_as_one` still
  sees its SERIAL warning. Opting the element in is a separate audit.
- **Vanilla ManzariDafalias IntScheme 2** is admitted (the same base path LadrunoSANISAND IS2
  exercises) but is not in the FE matrix: it has no substep cap, and one serial run of the deck
  exceeds the child timeout (the uncapped-substepping quirk, not a threading one).

## 6. Review of record

### 6.1 Fable critic (the other family), on `efc4113f1`: AMEND

It found no data race and no change to serial arithmetic on any admitted call graph. It traced
the ledger through multi-loop commit rounds, serial/threaded mixes, a marker in either kind of run,
stale partials, `wipe` and a refused marker commit: no data is lost or misattributed, no lock
nests, and the 1-thread path is the old arithmetic plus one atomic load. It also checked HS's
generation counter (all seven map writers bump under the mutex), the DP and ManzariDafalias call
graphs, SAS-ME and the `thread_local` conversions. Findings and what was done:

| # | Finding | Action |
|---|---|---|
| 1 | `avgImplexError` is folded per thread, so under a threaded commit loop it is schedule-dependent in its last bits, and no test read the ledger | **Owner decision (§7).** Gated: `test_wp258_soil_threads.py` now reads `implexError` (max), `implexRefusals` and `avgImplexError` after the threaded push. Max and counters must match exactly; the average must be within 1e-12 relative. |
| 2 | `thread-class-static` (and `thread-static`) scanned only files that spell the opt-in: the inheriting wrappers, DruckerPrager's base `.cpp`, `LadrunoSANISANDSasME.cpp` were unscanned | Fixed: `_thread_material_files` adds descendants (unless they refuse), bases up to the framework classes, and every file with a column-0 method definition of those classes. Four self-tests cover the cases. Reverting `LadrunoSANISAND3D::mSigma_M`'s `thread_local` now fails the lint. |
| 3 | The gates did not assert that the push ran or was plastic; vanilla ManzariDafalias IntScheme 2 was admitted but ungated | Fixed: every push step must return 0 at one thread, and the footing secant ratio must be < 0.95 in every case. Vanilla ManzariDafalias now admits IntScheme 1 only; LadrunoSANISAND admits 1, 2 and 129, which is what is gated. Not done: no FE element threads a plane-strain wrapper's `getTangent` (LadrunoQuad threads only its update loop), so those `thread_local` buffers are reached concurrently only through MPM 2D's `getInitialTangent` / `getStress`. The worker-print test cannot prove that a worker, rather than the master, printed. Its mutant M1 crashed, so in that run a worker printed. |
| 4 | A worker's `opserr` goes to fd 2; the master's goes to `sys.stderr` | Recorded as a limitation (quirk fragment). Where `sys.stderr` is not fd 2 (Jupyter on Windows), a worker's warn-once line is lost. A per-thread queue flushed on the caller after the loop would fix it; that is a follow-up. |
| 5 | Stale comments (`ManzariDafalias.h`, `LadrunoSANISAND.h`, `Material.h`) | Fixed. |
| 6 | Nits: check-then-increment warn budgets can print up to 10 + T − 1 lines under threads; the partial max was seeded from the first error, not 0 | The budgets are cosmetic and recorded here. The partial max is now seeded at 0 with `e > max`, as in the serial code. |

### 6.2 `/code-review high` (#1153, on `a6dd610f4`)

Ten findings: two fixed, eight recorded.

- **Fixed in `3d8ffc432`:**
  - the vanilla waiver comments lacked a `// Ladruno` marker, so `grep -rn Ladruno` did not find them;
  - the MPM gate re-spawned its thread-count probe for every case.
- **Recorded, not fixed:**
  - **The threaded first-round marker is a CAS winner.** If the thread count drops to 1 mid-run, serial rounds reset at a mid-loop point. Diagnostic only.
  - **`avgImplexError`:** decision 7.1.
  - **Worker prints go to fd 2:** quirk limitation; a per-thread queue is the follow-up.
  - **The worker-print test is probabilistic.** Mutant M1 did crash it.
  - **A fourth copy of the per-thread partials pattern.** A shared helper in `LadrunoThreads` is the follow-up.
  - **Warn budgets can overshoot under threads** (cosmetic).
  - **`tlLoopId` is zeroed, not restored.** Unreachable: no nested `parallelFor`.
  - **The lint cache never invalidates.** No long-lived caller.

## 7. Decisions for the owner

1. **`avgImplexError` under threads.** With LadrunoSANISAND `-implex` on an element whose commit
   loop threads (e.g. `stdBrick`; MPM mode F), the process-wide average implex error may differ
   from the 1-thread value in its last bits, and between runs. `implexError` (the max), every
   refusal and guard counter, and every displacement, stress and state stay bit-identical. ASDConcrete
   (WP-192) set the precedent. Options: (a) accept it as a diagnostic exception (this PR's
   default); (b) keep `-implex` serial on threaded commit loops, which costs `-implex`
   its commit-loop speedup; (c) an ordered replay of every commit's error, which needs a
   per-round buffer of every point's error (≈1 MB per 10⁵ points per round) and a flush hook.
