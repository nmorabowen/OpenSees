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

*(see below — filled from the read-only audit, before any SANISAND code change)*

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

*(filled from §2)*

## 5. Results

*(after the build of record)*
