# WP-260 Dispatch Briefs

Before starting any WP, read `SKILL.md`, `references/contracts.md`, `references/fork-checklist.md`, and `Ladruno_implementation/260_ladruno_solver_acceleration_adr.md` (§2.2 open-PR overlaps, §2.3 prior closed lanes). If a brief's **Blocked by** names an open PR or an owner decision, do not start: post a handoff on the `program` issue.

**Only the packages the owner picked after `260-w0` may start** (plan §0.5: phase share ≥ 10 %, projected whole-run gain ≥ 1.10×, top 2–3 at a time). The other briefs are kept ready, not scheduled.

---

## `260-w0` — Phase-0 Profiler Baseline (Read/Bench Only)
- **Branch:** `wp/260w0-phase0-profiles` from `origin/ladruno`.
- **Blocked by:** D1.
- **Goal:** measure where the time goes on the owner-approved D1 decks (3D implicit nonlinear solid/soil, implicit `LadrunoCPML` + `H5DRM`, explicit `CDL`/`ExplicitBathe` + `CPML`, Mode-A Implicit MPM, Model C MCFT `2CLD12`) and rank the WP-260 menu by projected whole-run gain. Zero `SRC/` edits; the deck scripts and a small summarizer go under `Ladruno_scripts/wp260_profiles/`.
- **How to measure (two runs per deck, from one 5-target build whose `ladrunoBuild` SHA is recorded):**
  1. **Phase shares:** `profiler start -deep -perStep` at **1 thread**. The deep gate cannot be on while a loop is threaded, so shares come from the 1-thread run only.
  2. **Wall time at production threads:** coarse `profiler start -perStep` at the thread count the deck is normally run with (state it).
  3. Run each 3 times, interleaved across decks; report the median and the spread.
- **Read the scopes carefully:** `elem.tangent` times the whole loop **including the `soe.addA` scatter**; element cost is the `elem_by_type` rows. Split assembly into kernel (`elem_by_type`) and scatter (`soe.addA`). Report `soe.symbolic` / `soe.factor` / `soe.trisolve` separately, plus `formUnbalance`, `update`, `commit.elements`, `commit.record`, and `dc.setSize` / `domainChanged` where present.
- **Deliverable (one comment on the `program` issue + the same tables in `Ladruno_scripts/wp260_profiles/README.md`):**
  1. *Phase-share table:* per deck, the % of wall time per scope above (1-thread run) and the production-thread wall time.
  2. *Phase → package map:* `soe.factor`/`soe.symbolic`/`soe.trisolve` → 260a (PARDISO levers), 260g, 260m; `soe.addA`, `dc.setSize` → 260a (MUMPS), 260f; `elem_by_type` (by element/material) → 260c1, 260c2, 260f, 260h, 260k; `update`/`formUnbalance` in explicit decks → 260b, 260l; `commit.record` → recorders (outside WP-260; flag it if it is large).
  3. *Ranking:* for every Wave-1/2 package and each deck, `p`, the conservative `s` used, and `G = 1 / ((1 − p) + p / s)`; mark which packages pass the gate (`p ≥ 10 %` and `G ≥ 1.10×` on at least one deck).
  4. *Proposed first pick:* the top 2–3 by `G`, for the owner to confirm.
- **Done When:** the comment is posted with all four parts, the scripts reproduce it from the recorded SHA, and Claude has reviewed the ranking (review protocol §4).

---

## `260a` — Wave 1: Implicit SOE & Outer-Loop Hygiene (`G-BYTE` / `G-OPT`)
- **Branch:** `wp/260a-soe-scatter-hygiene` from `origin/ladruno`.
- **Blocked by:** —
- **Owned Files:** `SRC/system_of_eqn/linearSOE/mumps/Mumps{,Parallel}{SOE,Solver}.{h,cpp}`, `SRC/system_of_eqn/linearSOE/pardiso/PARDISOGenLin{SOE,Solver}.{h,cpp}`, `SRC/analysis/algorithm/equiSolnAlgo/accelerator/KrylovAccelerator.cpp`, `SRC/analysis/algorithm/equiSolnAlgo/AcceleratedNewton.cpp`, `SRC/system_of_eqn/linearSOE/OpenSeesSystemCommands.cpp` (if CLI flags added).
- **Targets:**
  1. Replace linear scan in `MumpsSOE::addA` (`:313-397`) with binary search (`G-BYTE`).
  2. Store 1-based CSR/COO permanently in `MumpsSOE`/`MumpsParallelSOE`; eliminate `2*nnz` `rowA++/--` passes on `JOB=3` (`G-BYTE`).
  3. Replace insertion-sort `setSize()` in `MumpsSOE`, `MumpsParallelSOE` (also iterate local vertices only, not `0..N_global-1`), and `PARDISOGenLinSOE` with copy + `std::sort` (`G-BYTE`).
  4. Delete dead $O(k^2 N)$ loop in `KrylovAccelerator.cpp:200-216` (`G-BYTE`).
  5. Fix the `PARDISOGenLinSolver.cpp:366-368` comment (quote the oneMKL 2024.2 `iparm[24]` table text in the PR) and re-prove `-deterministic` (WP-132 CNR) bit identity. Expose opt-in `iparm[23]=1|10`, `iparm[1]=3`, `iparm[24]=2`, `mtype -2` `iparm[10/12]`, `pardiso_64` (`G-OPT`); each is refused under `-deterministic` unless proved deterministic there.
  6. Serialize `theReuse` in `AcceleratedNewton::sendSelf` (`OpenSeesSP`).
- **Done When:** 5-target build green; all static gates green; CI job `MKL 2024.2 compat (WP-148)` green; `G-BYTE` MD5 bit-identity proved; `-deterministic` identity re-proved; `-perStep` speedup logged in PR; Zone-A green.

---

## `260b` — Wave 1: Explicit Wave Hygiene + `LadrunoCPML` Threading (`G-BYTE` / `G-EQUIV`)
- **Branch:** `wp/260b-explicit-wave-hygiene` from `origin/ladruno`.
- **Blocked by:** —
- **Owned Files:** `SRC/analysis/integrator/{ExplicitBathe,LadrunoConsistentRefine,LadrunoMassScaling,CriticalTimeStep}*`, `SRC/system_of_eqn/linearSOE/diagonal/MPIDiagonalSOE*`, `SRC/analysis/handler/LadrunoConstraintProjector*`, `SRC/element/LadrunoElementShell.h`, `SRC/element/ladrunoBrick/LadrunoBrick{,20}.cpp`, `SRC/element/ladrunoPlane/LadrunoQuad.cpp`, `SRC/element/ladrunoCPML/LadrunoCPML*`, `SRC/analysis/integrator/TransientIntegrator.cpp`.
- **Targets:**
  1. Port CDL P-NEW-1 (`formTangent` mass cache, opt-out `-noMassCache`) + opt-in P-NEW-2 (`-commitSolveState`) + hoisted `-lnvd` (`V_fake` fast read, preallocated `lnvdDamped`) to `ExplicitBathe`.
  2. Ship ADR-67 P1.1/P1.2′: hoist `ConsistentPCGScratch`, cache MPI multiplicity `w` per `domainChanged`, preallocate `MPI_Request`/`Status` in `MPIDiagonalSOE`, pack 2-double `Allreduce`, add rank-1 $O(n)$ centroidal matvec for consistent SMS, add massless short-circuit before `DGGEV` in `CriticalTimeStep.cpp`.
  3. Add T7 zero-accel inertia skip to `LadrunoShell::addNodalInertia`, `LadrunoBrick20`, and `LadrunoCPML`; fix `LadrunoBrick::formUri(0)` `stiff.Zero()` / `getTangent()` & cache bulk-viscosity `c_d`.
  4. Pre-factor `LadrunoConstraintProjector` Cholesky $L$; add global $O(N)$ diagonal `alphaM` fast path; opt `LadrunoCPML` into `ladrunoThreadSafeUpdate/Form/Commit`.
- **Done When:** 5-target build green; static gates green; `G-BYTE`/`G-EQUIV` port of CDL/SMS/CPML/projector tests + `removeElement` cache invalidation test green; Zone-A green.

---

## `260c1` — Wave 1: `LogStrainNDMaterial` 5-Bug Fix + Stress-Only Contract #1a
- **Branch:** `wp/260c1-logstrain-stressonly` from `origin/ladruno`.
- **Blocked by:** —
- **Owned Files:** `SRC/material/nD/{NDMaterial.h,FiniteStrainNDMaterial.h,LogStrainNDMaterial*,LogStrainKernel.h}`, `SRC/element/ladrunoUP/LadrunoUP.cpp`.
- **Targets:**
  1. Add Contract #1a (`setTrialStrainStressOnly`/`setTrialFStressOnly`/`setTrialFincrStressOnly`, `ladrunoSupportsStressOnly`) with defaults that forward; implement it on `LogStrainNDMaterial` covering **both** `getTangent()` and `getSpatialTangentTensor()` (`G-BYTE`).
  2. Fix all 5 `LogStrainNDMaterial::setTrialF()` bugs (`:165` cache `epsN6`+`FnInv`, `:172/:200/:206` remove statics, `:201-205` cache `Ce = D0^-1`, `:210` defer `be_from_hencky_voigt` to commit (cache `tauV`/`epsEnp1`), `:191/:234` one spatial tangent, produced on demand) (`G-BYTE`).
  3. Opt `LadrunoUP` into OpenMP after guarding `:599` `windingNoteShown`.
- **Done When:** 5-target build green; static gates green; `G-BYTE` identity on the finite-strain battery; a test that a tangent requested after a stress-only trial equals the full-trial tangent bit for bit; Zone-A green.

---

## `260c2` — Wave 1: Concrete3D / HS / RCConcrete Stress-Only + IMPL-EX Deferral
- **Branch:** `wp/260c2-concrete-hs-implex` from `origin/ladruno` **after #1153 and #1154 merge**.
- **Blocked by:** #1153 (WP-258), #1154 (WP-259) — both edit `LadrunoHS*`, `LadrunoConcrete3D*`, `LadrunoRCConcrete*`.
- **Owned Files:** `SRC/material/nD/{LadrunoConcrete3D*,LadrunoConcrete3DKernel.h,LadrunoHS*,LadrunoHSKernel.h,LadrunoRCConcrete*}`.
- **Targets:**
  1. Contract #1a on these materials (all tangent channels; thread-safe per WP-258/259) (`G-BYTE`).
  2. `doTangent && !mp.implex` in `LadrunoConcrete3DKernel.h:3253` + commit-deferred `-implex` in `LadrunoConcrete3D` & `LadrunoHS` (with `HSsmall` `Ge_implex` freeze option), on the existing `LadrunoImplexClock` / `getImplexError()` (Contract #1b) (`G-BYTE`/`G-OPT`).
- **Done When:** as 260c1, plus the WP-258/259 race-stress tests still green.

---

## `260d` — Wave 1: `LadrunoAdaptiveTransient` + Effective-Tangent Key + Modal Vectorization
- **Branch:** `wp/260d-adaptive-transient` from `origin/ladruno`.
- **Blocked by:** —
- **Owned Files:** `SRC/analysis/integrator/IncrementalIntegrator.{h,cpp}`, `SRC/analysis/integrator/TRBDF2.{h,cpp}`, `SRC/analysis/integrator/TRBDF3.{h,cpp}`, new `SRC/analysis/integrator/LadrunoAdaptiveTransient*`, `SRC/analysis/algorithm/equiSolnAlgo/LadrunoFactorReuse*`.
- **Targets:** Contract #2 `getEffectiveTangentKey` (base + `TRBDF2`/`TRBDF3` + `LadrunoFactorReuse` keying); `LadrunoAdaptiveTransient` (HAFTOL half-step residual LTE + Hairer–Wanner PI + quantized $2^k/\sqrt{2}^k$ rungs that keep the factor-reuse key stable); vectorized flat `addModalDampingForce` (`IncrementalIntegrator.cpp`).
- **Not in scope:** `LadrunoImplicitBathe`, Contract #3, `Domain.{h,cpp}` — see `260d-B`.

## `260d-B` — Gated: Implicit Bathe (only after D6)
- **Blocked by:** owner decision **D6**. ADR-77 G3 closed implicit Bathe with the benchmark on record.
- **Step 1 (the D6 benchmark, pre-registered):** shared factor across both substeps under `LadrunoFactorReuse` vs `TRBDF2`/`GeneralizedAlpha`, on the ADR-77 T3 deck plus a contact/wave deck, cost per simulated second, bar ≥ 1.3×. Post the result on the `program` issue. Below the bar: stop, the lane stays closed.
- **Step 2 (only if it clears):** `LadrunoImplicitBathe` ($\rho_\infty$-Bathe, $\gamma_0(\rho_\infty)$), Contract #3 (`ScopedRecorderSuppression` separate from `contactAugmenting`, the commit-side-effect table, substep 1 held as trial state), and the fork-owned substep driver. Owns `SRC/domain/domain/Domain.{h,cpp}`.

---

## `260e` — Wave 1: Windows Single OpenMP (`libiomp5md`), `mimalloc`, `/LTCG`
- **Branch:** `wp/260e-allocator-lto-omp` from `origin/ladruno`.
- **Blocked by:** D4.
- **Targets:** Eliminate `vcomp140.dll` dual-pool linking (`CMakeLists.txt:774-781`, verified via `dumpbin /dependents` on every DLL/EXE in `dist\`), wire `mimalloc` + `/GL /LTCG` with `/fp:precise`. **No `/arch:AVX2`** (ADR-77 V0).

---

## Wave 2 (`260f–260j`) & Wave 3 (`260k–260o`)
Cut fresh branches from `origin/ladruno` using the contracts in `references/contracts.md`:
- **`260f`:** Full 3D/MPM-A/MCFT de-static + `LadrunoFastScatter` (Contract #4 NVI `setSize`/`doSetSize`) + persistent `LadrunoThreadedAssembly` + sparse `TransformationFE` $T^T K T$. *Blocked by #1147 (MPM, `LadrunoThreads`) and #1153 (fiber sections).*
- **`260g`:** FP32 PARDISO (`float*` shadow) + MUMPS GMRES-IR + `LadrunoSymPrecondKrylovSOE`. oneMKL 2024.2.2 floor.
- **`260h`:** Analytic `Concrete3D`/`HS` return maps, `SAS-RK23/45`, `ASDPlastic3D` closed-form VM/DP, Eisenstat–Walker (`LinearSOESolver::setRelativeTolerance`) + `LadrunoAndersonNewton` + `LadrunoLBFGS`. *Blocked by #1153/#1154.*
- **`260i`:** `LadrunoImplexAlgo` (Contract #1b) + `LadrunoSuperelement` (`OPS_LadrunoCMS`) + `H5DRMLoadPattern` v2.
- **`260j`:** `LadrunoContact` 3D Mortar bucket sort/BVH + `LadrunoContactFE` threading.
- **`260k–260o`:** `LadrunoGPStore` SoA + 7×7 Schur CPPM; `LadrunoExplicitEngine` + Clustered LTS (CPML memory variables FP64); `LadrunoAMGSOE` (`AMGCL` + Biot Schur) + `ClusterPardiso` + weighted METIS; `LadrunoIMEX` (Hughes–Liu → Prakash–Hjelmstad); `LadrunoSpectralHex` + GPU explicit/MPM (GPU only per D3).
