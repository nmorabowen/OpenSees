# WP-260 Dispatch Briefs

Before starting any WP, read `SKILL.md`, `references/contracts.md`, `references/fork-checklist.md`, and `Ladruno_implementation/260_ladruno_solver_acceleration_adr.md` (§2.2 open-PR overlaps, §2.3 prior closed lanes). If a brief's **Blocked by** names an open PR or an owner decision, do not start: post a handoff on the `program` issue.

**Only the packages the owner picked after `260-w0` may start** (plan §0.5: phase share ≥ 10 %, projected whole-run gain ≥ 1.10×, top 2–3 at a time). The other briefs are kept ready, not scheduled.

---

## `260-w0` — Phase-0 Profiler Baseline (Read/Bench Only)
- **Branch:** `wp/260w0-phase0-profiles` from `origin/ladruno`.
- **Blocked by:** D1. The **ranking run** also waits for #1153 (thread-safe soils) and #1154 (thread-safe concretes) to merge. Before they merge, SANISAND and Concrete3D refuse the threaded loops, so no production-thread number for D1a/D1b is real. The scripts and the small decks can be fixed and smoke-tested now.
- **Goal:** measure where the time goes on the owner-approved D1 decks, at production size and with production settings, and rank the WP-260 menu by projected whole-run gain. Zero `SRC/` edits; the deck scripts and the summarizer go under `Ladruno_scripts/wp260_profiles/`.
- **Why this brief is strict:** the first sweep (`a7a139a76`, #1158) got the scope arithmetic right but ranked on decks 10–100× below production, on a linear deck that re-factored every iteration, on an "8T" column that ran serial, and on `s` values outside rule 2. The review is on #1156. Each rule below closes one of those holes.

### 1. Decks — the D1 configuration, at two sizes
Every deck runs as D1 states it and as a user would run it in production. A deviation from D1 needs an owner amendment recorded on the `program` issue **before** the run, not a note after it.

| Deck | D1 item | Required configuration | Sizes |
|---|---|---|---|
| D1a | (1) | `LadrunoBrick` + `LadrunoSANISAND`, `Pardiso` **with `LadrunoFactorReuse`**, `Newmark` | S (≈ 360 hex) and L (≥ 8× S elements) |
| D1b | (1) | `LadrunoBrick` + `LadrunoConcrete3D`, `Pardiso` **with `LadrunoFactorReuse`** | S and L (≥ 8× S) |
| D2 | (2) | `LadrunoBrick` interior + `LadrunoCPML` jacket + `H5DRM` input. A linear deck runs `algorithm Linear -factorOnce` (`Linear.cpp:67`) after a check that the tangent is constant at fixed `dt`; a nonlinear one runs `LadrunoFactorReuse` | S and **L ≥ 10⁵ DOF** (the solver-share probe) |
| D3a / D3b | (3) | `CentralDifferenceLadruno` / `ExplicitBathe`, `LadrunoBrick -formulation uri` + `LadrunoCPML -lumped`, `Diagonal` | S and L (≥ 8× S) |
| D3c | (3) variant | `CentralDifferenceSMSConsistent`, `Pardiso` | S |
| D4 | (4) | Mode-A `LadrunoMPMImplicit` | S (the MPM program, #1075, owns its scaling) |
| D5 | (5) | `2CLD12` column, `LadrunoForceBeamColumn -monolithic` + `LadrunoNDFiber` + `LadrunoMCFTStrip` | as tested (a single column: its size is the specimen) |

The size pair turns each share into a trend. Factorization grows about as n² in 3D and element work as n, so the post reports `p` at S and at L, and ranks on L.

### 2. Runs — three modes, one build, the production environment
From one 5-target build; every record stores its own `ops.ladrunoBuild()` and the summarizer refuses mixed SHAs (the `.pyd` and `OpenSees.exe` must match).

| Mode | Profiler | Threads | Used for |
|---|---|---|---|
| `1t_deep` | `-deep -perStep` | 1 | phase shares (the deep gate refuses threaded loops) |
| `1t_coarse` | `-perStep` | 1 | the serial wall time; deep overhead = `1t_deep / 1t_coarse` |
| `nt_coarse` | `-perStep` | production count (state it; 8 on the i7-13700KF P-cores) | production wall time; scaling = `1t_coarse / nt_coarse` |

- **Environment:** no `OMP_WAIT_POLICY`, `KMP_BLOCKTIME` or affinity override that the fork does not set by default. Record every `OMP_*`, `KMP_*`, `MKL_*` and `LADRUNO_*` variable, and the effective affinity mask (checked with `GetProcessAffinityMask`), in the manifest.
- **Threading engaged:** keep each run's stderr. The `ladrunoThreads` announcement says whether the loops threaded or refused; the table shows it per deck. A "production" row that ran serial is labelled so.
- **Repetitions:** interleave decks and modes; drop step 1 from every `-perStep` series (symbolic factorization, allocation, first-touch); 3 reps, raised to 5 for any deck whose spread exceeds 5 % of the median. Report median and min..max.

### 3. Reading the scopes
- `elem.tangent` times the whole loop, including the `soe.addA` scatter. Element cost is the `elem_by_type` rows, whose per-element timer wraps only `getTangent` (`IncrementalIntegrator.cpp:119-131`). So `elem_by_type + soe.addA` is a disjoint split. Count the `elem_by_type` rows under `elem.tangent`, `elem.residual` and `elem.update` (`Domain.cpp:2521`).
- Report `soe.symbolic` / `soe.factor` / `soe.trisolve` / `soe.cgs` separately, plus `formUnbalance`, `update`, `newStep`, `commit.elements`, `commit.record` and `domainChanged` / `dc.setSize`.
- **MPM:** report `mpm.a.prepare` in its own column. It is MPM host work for the MPM program (#1075), not `domainChanged`.
- **Accounting check:** print `step_wall − Σ(parts)` per deck. If the parts sum to more than 100.5 %, something is double counted: fail. Never clamp `other` at zero.
- **Unexplained time:** any remainder (`other`, or a scope's time outside its named children) above 5 % of a deck is resolved from `summarize.py --dump-trees` before that deck is ranked. The first sweep's D1b `update` remainder (28.9 %) is the example.

### 4. The ranking — `p` and `s` by rule, not by choice
- **`p` is only the phases the package changes,** named per package in the summarizer. A constant `p` (whole step, `p = 1.0`) is not allowed. If the gain cannot be written as a phase share, the package is "not phase-rankable" and needs its own whole-run benchmark: `260d` (adaptive stepping changes the step count) and `260e` (allocator/LTO).
- **`s` follows plan §0.5 rule 2:** take the low end of the plan's hypothesis. If it was never measured, halve the excess, `s = 1 + (s_low − 1) / 2`. A package with no number in the plan is "unranked" until the post states an `s` with its source (a micro-benchmark, or a cited measurement). The summarizer keeps the table: package, plan hypothesis, `s` used, source.

  | Pkg | Plan hypothesis | `s` (unmeasured) | `p` = |
  |---|---|---:|---|
  | `260a` | 1.3–2× PARDISO factor | 1.15 | `soe.factor` on `Pardiso` decks (the `G-OPT` `iparm` levers); `soe.addA` + `dc.setSize` only on a MUMPS deck |
  | `260b` | 1.4–2.2× explicit | 1.20 | `LadrunoCPML` `elem_by_type` + `commit.elements`; plus `ExplicitBathe` `formTangent` on D3b; consistent-SMS `soe.trisolve` on D3c |
  | `260c1` | 3–6× finite-strain GP | 2.00 | `LogStrain` material share (needs a finite-strain deck; none in D1 → no `p`) |
  | `260c2` | 1.5–4× | 1.25 | `LadrunoConcrete3D` / `LadrunoHS` material share |
  | `260d` | 2–5× adaptive runs | — | not phase-rankable |
  | `260e` | to be measured | — | not phase-rankable |
  | `260f`–`260j` | no number | unranked until sourced | `260f`: `soe.addA` + assembly only; element threading that #1147/#1153/#1154 deliver is theirs, not `260f`'s |
  | `260k`–`260o` | Wave 3 | — | not eligible (rule 6); listed only for context |

- **Per deck:** the post gives the full `p` / `s` / `G` matrix (package × deck), not only each package's peak deck.

### 5. Deliverable
- One comment on the `program` issue, and the same tables in `Ladruno_scripts/wp260_profiles/RESULTS.md`. `README.md` stays the how-to (dependencies, commands, where outputs go). The summarizer writes `RESULTS.md` and never overwrites `README.md`.
- The comment has five parts:
  1. *Phase-share table:* per deck and size, the share of each scope above (from `1t_deep`), the residual, and wall time for the three modes with spread. It also has a "threaded?" column.
  2. *Size trend:* `p` at S and at L for the solver and element phases.
  3. *Phase → package map:* `soe.*` → 260a, 260g, 260m; `soe.addA`, `dc.setSize` → 260a (MUMPS), 260f; `elem_by_type` (by element and material) → 260c1, 260c2, 260f, 260h, 260k; explicit `update` / `formUnbalance` / `commit.elements` → 260b, 260l; `mpm.a.prepare` → the MPM program; `commit.record` → recorders (outside WP-260; flag it if ≥ 10 %). Name phases ≥ 10 % that no package targets (for example D5's `LadrunoMCFTStrip` kernel) as gaps.
  4. *Ranking:* the §4 matrix; mark pass (`p ≥ 10 %` and `G ≥ 1.10×` on at least one deck), not eligible, unranked, and not phase-rankable.
  5. *Proposed pick:* the top 2–3 passing, eligible packages by `G` (ties: `G-BYTE` before `G-EQUIV` before `G-OPT`), with each pick's blocker.
- Commit `manifest.json`, `trees.json` (each run's parsed scope tree with its `elem_by_type` rows) and the `--dump-trees` text under `Ladruno_scripts/wp260_profiles/results/<sha9>/`. The `.h5` files are large and stay local; state their path.

### 6. Done When
- The comment is posted with all five parts.
- `python Ladruno_scripts/wp260_profiles/summarize.py --from-json Ladruno_scripts/wp260_profiles/results/<sha9>` reproduces it from the committed files alone, without the `.h5` files or a build.
- Every deck's residual is within ±0.5 %, and no unexplained remainder above 5 % is left.
- Claude's review protocol §4 verdict on the ranking is `ready`.

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
