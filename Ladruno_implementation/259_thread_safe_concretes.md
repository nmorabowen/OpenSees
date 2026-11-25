# WP-259 — thread-safe Ladruno concretes (`LadrunoConcrete3D`, `LadrunoRCConcrete`)

Board #937 · stacked on WP-258 (#1153, itself stacked on WP-253 #1147) · follows the WP-258 pattern
(`258_thread_safe_soil_materials.md`).

**Why.** Neither Ladruno concrete answers `Material::ladrunoThreadSafeUpdate()` /
`ladrunoThreadSafeForm()` (the WP-107 / WP-192 contract, `Material.h`, `Element.h` items 1–5), so
every host that holds one (`stdBrick`, `LadrunoQuad`, `ASDShellQ4` + `LayeredShell`) falls back
to its serial loop. `ASDConcrete3D` already threads (WP-192). This WP makes both
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
| any of the above under `recorder EnergyBalance` (WP-196, `EnergyBalanceRecorder.cpp:784`) | as the host | the element sweep (on `ladrunoThreadSafeForm`; calls the same Form methods) |
| `PlaneStressMaterial` / `PlateFiberMaterial` / `PlateFromPlaneStressMaterial` wrappers | ThreeDimensional / PlaneStress copies | whatever their host threads (they delegate; every view they request is admitted) |
| MPM mode F (WP-253) | — | none: `mpm body` refuses both concretes at admission (ADR 232 §6 admits ElasticIsotropic, LadrunoJ2, DruckerPrager, LadrunoSANISAND, LadrunoHS; found while writing the gate) |

No threaded host requests AxiSymmetric or BeamFiber (`LadrunoNDFiberSection3d` does not opt in).
**A condition of the opt-in (critic, §6).** `ops_TheActiveElement` is set per element by the
threaded update loop and by the tangent / residual loops, but the threaded commit loop
(`Domain.cpp:2884-2890`) leaves workers with a stale value. Neither `commitState` reads it today
(LadrunoConcrete3D reads `lch` only in `integrate()` from `setTrialStrain`; RC latches it once in
the first `setTrialStrain`, and its IMPL-EX commit re-integration does not re-read it). A future
change that reads the element at commit must set it there, or 1-thread identity breaks.

The opt-in therefore reaches the FE element loops only; if MPM ever admits a concrete, the §2
audit already covers the four methods mode F calls on workers, but a mode-F gate must come with
the admission.

## 3. Determinism gate

`tests/test_wp259_concrete_threads.py` (zone_a): every threaded result bit-identical to 1 thread
at 2, 4 and 8 threads (`float.hex` of every nodal displacement after every step, plus the load
factor and every step's return code), shuffled element tags, and the deck must really crack /
crush (asserted from the material's own damage response). Meshes: a `stdBrick` cantilever in
bending (`3d`), the same prism crushed axially (`3dc`, compressive damage asserted) and, for
LadrunoConcrete3D, pushed past the crush (`3dcr`: return maps fail in volume, steps are refused
and the process-wide warning budget saturates — the one shared counter, gated by its printed
lines); `LadrunoQuad` PlaneStrain (`pe`) and PlaneStress (`ps`, the condensation Newton);
`ASDShellQ4` + `LayeredShell` (`shell`, PlateFiber). Options: default, `-implex`, `-eta`,
`-autoRegularization` (LadrunoConcrete3D); default, `-implex`, `-autoRegularization`,
`-interlock -cyclic` (LadrunoRCConcrete). Twenty cases; the digest is every nodal displacement
and the load factor after every step, every step's return code, every Gauss point's damage
(plus `substeps` / `returnFailures`, or RC `implexError`), and the warning-line counts.

After the critic (§6) the gate also: asserts the tangent AND the residual announcement for every
form-loop host (one string used to be satisfied by the tangent line alone); counts every other
worker print (σ₂₂ / σ₃₃ condensation, `-verbose` commit summaries, RC `-implexControl`) against
1 thread; adds `c3d_verbose` on `3dcr` (no budget: every failure prints), `rc_implexctl` on
`shell` (the commit loop prints from workers) and the EnergyBalance sweep (`c3d_energy` `3d+eb`,
`rc_energy` `shell+eb`, the recorder file compared byte for byte). The crush decks were moved off
their cliff: `3dc` requires only its first 8 of 10 steps to converge and "cracked or crushed"
(LadrunoConcrete3D's axial crush is sharp: no compressive damage after step 9, ω_c 0.17 after
step 10, refusals from step 11 on nmb), and `3dcr` refuses from step 7 of 9 (increment 0.12 mm).
**Not provable from the output:** the threaded commit loop prints no announcement
(`Domain::ladrunoThreadedCommit`), so the gate shows bit-identity with the commit loop enabled
but cannot assert that it ran threaded; the mutation evidence (§5) targets the tangent loop. Neither material keeps
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

**Binary of record:** `7c01546a3` (`ladrunoBuild` of `dist/bin/OpenSees.exe` and
`dist/bin/opensees.pyd` checked before every run; built with `Ladruno_scripts\build.bat`, all five
targets). Later commits change tests and documents only. Machine nmb (i7-13700KF), Windows 11,
MKL 2026.1. **The machine was not idle:** other sessions held it at 28–34 % total CPU during the
bench and stress runs, so the speedups below are lower bounds of what an idle nmb gives. Data:
`Ladruno_files/testbed/perf/wp259/`.

### 5.1 Determinism gate (Windows, `7c01546a3`)

`tests/test_wp259_concrete_threads.py`: **24 / 24 cases bit-identical to 1 thread at 2, 4 and 8
threads.** The comparison covers:

- every nodal displacement and the load factor after every step;
- every return code;
- every Gauss point's diagnostics;
- the warning-line counts;
- in the EnergyBalance cases, the recorder file.

Every case asserts that the update loop announced itself on n threads; on `stdBrick` and
`ASDShellQ4` cases, the tangent and residual loops too. Every base run cracks (ω_t 0.48–0.94) or
crushes (RC ω_c ≈ 0.3; LadrunoConcrete3D κ_p 0.93). The exception is `3dcr`, which is meant to
refuse: steps 7–9 of 9 refuse, printing 20 warning lines plus the suppression line (4 371 lines
under `-verbose`).

### 5.2 Mutation evidence (working-tree mutants, built and run, never committed)

- **M1 + M2, built together.** `LadrunoConcrete3D::getTangent` and
  `LadrunoRCConcrete::getTangent` copy into, and return, a function-static `Matrix` shared by
  every instance: the classic shared return buffer.
- **Result: 15 red, 4 green**, against the gate as it stood then (19 cases, `3dcr` deselected).
  - Red: every `stdBrick` and `ASDShellQ4` case of both materials, at 2 threads. Examples:
    `c3d-3d` "nodal displacements differ at push step 0"; `rc-shell` load factor
    `0x1.9612405e55abcp+15` vs `…ab9p+15` at step 0.
  - Green: the four `LadrunoQuad` cases. This is expected: `LadrunoQuad` threads only its
    update loop, which never calls `getTangent` concurrently.
- The mutants were reverted (`git checkout`) and the binary of record rebuilt from a clean tree.
- The mutants target the tangent loop. The commit loop prints no announcement (§3), so it is gated
  by bit-identity only.

### 5.3 Speed (reported; no pre-registered target)

Source: `Ladruno_scripts/wp259_bench.py`, `bench_7c01546a3.json`.

- **Deck:** the gate's `3d` cantilever refined to 40 × 8 × 12 = 3 840 `stdBrick` (30 720 material
  points), 6 DisplacementControl steps.
- **Statistic:** 3 repeats; median wall time of the push only.
- **`mnewton`:** ModifiedNewton `-initial` + UmfPack, FixedNumIter 20. One factorization; every
  iteration is the threaded update and residual loops plus one triangular solve.
- **`newton`:** full Newton + `system Pardiso`, FixedNumIter 6. The tangent loop threads too, and
  there is one factorization per iteration.

| Material | Algorithm | 1 thread (s) | 4 threads (s) / S(4) | 8 threads (s) / S(8) | spread (max) | 1-vs-n digest |
|---|---|---|---|---|---|---|
| LadrunoConcrete3D | mnewton | 53.7 | 28.9 / **1.86** | 23.0 / **2.34** | 10 % | equal |
| LadrunoConcrete3D | newton | 418.6 | 111.6 / **3.75** | 79.8 / **5.24** | 24 % | n/a (PARDISO) |
| LadrunoRCConcrete | mnewton | 16.0 | 11.1 / **1.45** | 10.5 / **1.53** | 2 % | equal |
| LadrunoRCConcrete | newton | 11.9 | 6.5 / **1.81** | 5.2 / **2.30** | 8 % | n/a (PARDISO) |

**Reading it.**

- Before this WP, every row ran its element loops serially.
- **LadrunoConcrete3D `newton` is the material-dominated row.** Its last two steps refuse: return
  maps fail in volume, which is the expensive path. So the tangent and update loops dominate,
  and the row scales 5.2× at 8 threads.
- **The `mnewton` rows are capped** by the serial UmfPack triangular solve and the serial scatter.
- **LadrunoRCConcrete is cheap per call**, so the same serial remainder caps it at 1.5×.
- No row is slower threaded.
- In every `mnewton` row the threaded digests equal 1 thread. The PARDISO rows are not compared:
  PARDISO is not byte-reproducible run to run by design (ADR-75 P1f).

### 5.4 Race stress

`Ladruno_scripts/wp259_race_stress.py --repeats 10 --threads 8`: **0 mismatches in 80 runs.** Eight cases, each 10 runs at 8 threads against its 1-thread digest (`race_stress_7c01546a3.json`): `c3d-3d`, `c3d-3dcr` (refusals + saturated budget), `c3d_implex-shell`, `c3d-ps`, `rc-3d`, `rc_implex-shell`, `rc_interlock-shell`, `rc-ps`.

### 5.5 Not done here

- **Zone-A (GCC / Ubuntu).** It is dispatched at the ready flip; its run id is in the handoff.
- **TSan.** Not run: GCC's libgomp is not TSan-instrumented (the WP-258 §5.4 reasoning holds). The
  race stress is the fallback.
- **A commit-loop announcement**, so that a gate can assert the commit loop threaded. It needs a
  vanilla `Domain.cpp` edit, which is out of scope; recorded as a follow-up.
- **A shared warn-budget helper** in `LadrunoThreads.h`; a follow-up. The exact `fetch_add` claim
  of §2.1 is now written out in LadrunoConcrete3D. DruckerPrager and SANISAND keep
  check-then-increment budgets, which can overshoot by T − 1 lines.

## 6. Review of record

### 6.1 Fable critic (the other family), on `0724afa48`: AMEND

**Verdict on the code.** No shared mutable state was missed and the opt-in conditions are right.
Every wrapper that delegates to these concretes (`PlaneStressMaterial`, `PlateFiberMaterial`,
`PlateFromPlaneStressMaterial`) lands on an admitted view. Serially, the warning-budget rewrite
prints the same lines in the same order.

**What it verified:**

- the kernel statics;
- the per-instance buffers, and that `getCopy` shares no pointers;
- `ops_TheActiveElement` in each loop;
- every `getCharacteristicLength` the hosts reach;
- the output guards (none nest);
- MPM's class-tag admission.

Findings and what was done:

| # | Finding | Action |
|---|---|---|
| 1 | #1154 targets `ladruno`, so its diff shows WP-258's commits; WP-258's Zone-A was red at the time | Base kept at `ladruno`, the brief's rule; the PR says it is stacked and which commits to review. #1153's Zone-A was re-running when checked. This PR's Zone-A is dispatched at the flip. |
| 2 | §5 empty; no mutation or stress evidence recorded | §5 written (this section). |
| 3 | The gate's "FE_Elements computed in parallel" string was satisfied by the tangent line alone; the commit loop prints nothing | Fixed: the tangent and residual lines are asserted separately. The commit loop is recorded as not provable from the output (§3, §5.5). |
| 4 | `3dc` / `3dcr` sit at a cliff; Zone-A (glibc) may move it | Fixed. `3dc` stops at 8 steps, three short of the crush cliff, and gates on κ_p / RC ω_c. `3dcr` uses a larger increment and refuses from step 7 of 9, a margin of 2–3 steps. |
| 5 | The WP-196 EnergyBalance sweep threads on `ladrunoThreadSafeForm` and is neither audited nor gated | Fixed: a §2.3 row, plus two cases, `c3d_energy` (`3d+eb`) and `rc_energy` (`shell+eb`). They compare the recorder file byte for byte and assert the sweep announcement. |
| 6 | `-verbose`, RC `-implexControl`, the Gc/lch warning, the wrappers and an RC refusal are not exercised under threads | Added `-verbose` (`c3d_verbose` on `3dcr`) and `-implexControl` (`rc_implexctl`), with every worker print counted. Still ungated: the Gc/lch warning, explicit wrapper decks and an RC refusal. The `ASDShellQ4` PlateFiber cases run the wrappers' shared code. |
| 7 | The 20th warning and the suppression line took separate guards | Fixed: one guard holds both. |
| 8 | The 32-bit `long` counter could wrap after 2³¹ failures | Fixed: the claim saturates (no `fetch_add` once the budget is spent). |
| 9 | Stale "MPM" wording | Fixed (comments, test docstring). |
| 10 | Commit workers keep a stale `ops_TheActiveElement` | Recorded as a condition of the opt-in (§2.3). |

### 6.2 `/code-review high`, on `7c01546a3`

Five findings: three fixed in the next commit, two recorded.

- **Fixed:**
  - `3dc`'s "cracked or crushed" check could pass on tensile damage alone. It now checks κ_p
    (LadrunoConcrete3D) or ω_c (RC).
  - The EnergyBalance temp file was created for every run.
  - The substring that counts commit summaries was loose; it is now `) commit: `.
- **Recorded:**
  - the commit loop cannot be asserted threaded (§5.5);
  - the exact warn-budget claim should be a shared helper (§5.5).
