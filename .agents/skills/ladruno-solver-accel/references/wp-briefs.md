# WP-260 Dispatch Briefs

Before starting any WP, read `SKILL.md`, `references/contracts.md`, `references/fork-checklist.md`, and `Ladruno_implementation/260_ladruno_solver_acceleration_adr.md`.

---

## `260-w0` — Phase-0 Profiler Baseline (Read/Bench Only)
- **Branch:** `wp/260w0-phase0-profiles` from `origin/ladruno`.
- **Goal:** Run `-perStep` + `-deep` profiling on the owner-approved D1 benchmark decks (3D implicit nonlinear solid/soil, implicit `LadrunoCPML` + `H5DRM`, explicit `CDL`/`ExplicitBathe` + `CPML`, Mode-A Implicit MPM, Model C MCFT `2CLD12`) and post exact phase time percentages to the `WP-260` `program` issue. Zero `SRC/` edits.

---

## `260a` — Wave 1: Implicit SOE & Outer-Loop Hygiene (`G-BYTE` / `G-OPT`)
- **Branch:** `wp/260a-soe-scatter-hygiene` from `origin/ladruno`.
- **Owned Files:** `SRC/system_of_eqn/linearSOE/mumps/Mumps{,Parallel}{SOE,Solver}.{h,cpp}`, `SRC/system_of_eqn/linearSOE/pardiso/PARDISOGenLin{SOE,Solver}.{h,cpp}`, `SRC/analysis/algorithm/equiSolnAlgo/accelerator/KrylovAccelerator.cpp`, `SRC/analysis/algorithm/equiSolnAlgo/AcceleratedNewton.cpp`, `SRC/system_of_eqn/linearSOE/OpenSeesSystemCommands.cpp` (if CLI flags added).
- **Targets:**
  1. Replace linear scan in `MumpsSOE::addA` (`:313-397`) with binary search (`G-BYTE`).
  2. Store 1-based CSR/COO permanently in `MumpsSOE`/`MumpsParallelSOE`; eliminate `2*nnz` `rowA++/--` passes on `JOB=3` (`G-BYTE`).
  3. Replace insertion-sort `setSize()` in `MumpsSOE`, `MumpsParallelSOE` (also iterate local vertices only, not `0..N_global-1`), and `PARDISOGenLinSOE` with copy + `std::sort` (`G-BYTE`).
  4. Delete dead $O(k^2 N)$ loop in `KrylovAccelerator.cpp:200-216` (`G-BYTE`).
  5. Fix wrong `:367` `iparm[24]` comment in `PARDISOGenLinSolver.cpp`; expose opt-in `iparm[23]=1|10`, `iparm[1]=3`, `iparm[24]=2`, `mtype -2` `iparm[10/12]`, `pardiso_64` (`G-OPT`).
  6. Serialize `theReuse` in `AcceleratedNewton::sendSelf` (`OpenSeesSP`).
- **Done When:** 5-target build green; all static gates green; `G-BYTE` MD5 bit-identity proved; `-perStep` speedup logged in PR; Zone-A green.

---

## `260b` — Wave 1: Explicit Wave Hygiene + `LadrunoCPML` Threading (`G-BYTE` / `G-EQUIV`)
- **Branch:** `wp/260b-explicit-wave-hygiene` from `origin/ladruno`.
- **Owned Files:** `SRC/analysis/integrator/{ExplicitBathe,LadrunoConsistentRefine,LadrunoMassScaling,CriticalTimeStep}*`, `SRC/system_of_eqn/linearSOE/diagonal/MPIDiagonalSOE*`, `SRC/analysis/handler/LadrunoConstraintProjector*`, `SRC/element/LadrunoElementShell.h`, `SRC/element/ladrunoBrick/LadrunoBrick{,20}.cpp`, `SRC/element/ladrunoPlane/LadrunoQuad.cpp`, `SRC/element/ladrunoCPML/LadrunoCPML*`, `SRC/analysis/integrator/TransientIntegrator.cpp`.
- **Targets:**
  1. Port CDL P-NEW-1 (`formTangent` mass cache, opt-out `-noMassCache`) + opt-in P-NEW-2 (`-commitSolveState`) + hoisted `-lnvd` (`V_fake` fast read, preallocated `lnvdDamped`) to `ExplicitBathe`.
  2. Ship ADR-67 P1.1/P1.2′: hoist `ConsistentPCGScratch`, cache MPI multiplicity `w` per `domainChanged`, preallocate `MPI_Request`/`Status` in `MPIDiagonalSOE`, pack 2-double `Allreduce`, add rank-1 $O(n)$ centroidal matvec for consistent SMS, add massless short-circuit before `DGGEV` in `CriticalTimeStep.cpp`.
  3. Add T7 zero-accel inertia skip to `LadrunoShell::addNodalInertia`, `LadrunoBrick20`, and `LadrunoCPML`; fix `LadrunoBrick::formUri(0)` `stiff.Zero()` / `getTangent()` & cache bulk-viscosity `c_d`.
  4. Pre-factor `LadrunoConstraintProjector` Cholesky $L$; add global $O(N)$ diagonal `alphaM` fast path; opt `LadrunoCPML` into `ladrunoThreadSafeUpdate/Form/Commit`.
- **Done When:** 5-target build green; static gates green; `G-BYTE`/`G-EQUIV` port of CDL/SMS/CPML/projector tests + `removeElement` cache invalidation test green; Zone-A green.

---

## `260c` — Wave 1: `LogStrainNDMaterial` 5-Bug Fix + Stress-Only Contract + IMPL-EX Deferral
- **Branch:** `wp/260c-logstrain-stressonly-implex` from `origin/ladruno` (rebased after `#1153`/`#1154` if touching `Concrete3D`/`HS`).
- **Owned Files:** `SRC/material/nD/{NDMaterial.h,FiniteStrainNDMaterial.h,LogStrainNDMaterial*,LogStrainKernel.h,LadrunoConcrete3D*,LadrunoConcrete3DKernel.h,LadrunoHS*,LadrunoHSKernel.h,LadrunoRCConcrete*}`, `SRC/element/ladrunoUP/LadrunoUP.cpp`.
- **Targets:**
  1. Add Contract #1 `setTrialStrainStressOnly`/`setTrialFStressOnly` + lazy `getTangent()` materialization on Ladruno materials (`G-BYTE`).
  2. Fix all 5 `LogStrainNDMaterial::setTrialF()` bugs (`:165` cache `epsN6`+`FnInv`, `:172/:200/:206` remove statics, `:201-205` cache `Ce = D0^-1`, `:210` defer `be_from_hencky_voigt` to commit, `:191/:234` lazy single spatial tangent) (`G-BYTE`).
  3. Immediate `doTangent && !mp.implex` in `LadrunoConcrete3DKernel.h:3253` + commit-deferred `-implex` in `LadrunoConcrete3D` & `LadrunoHS` (with `HSsmall` `Ge_implex` freeze option) after `#1153/#1154` merge (`G-BYTE`/`G-OPT`).
  4. Opt `LadrunoUP` into OpenMP after guarding `:599` `windingNoteShown`.

---

## `260d` — Wave 1: `LadrunoImplicitBathe` + `LadrunoAdaptiveTransient` + Modal Vectorization
- **Branch:** `wp/260d-implicit-bathe-adaptive` from `origin/ladruno`.
- **Targets:** Implement Noh–Bathe $\rho_\infty$-Bathe with $\gamma_0(\rho_\infty)$, Contract #2 `getEffectiveTangentKey`, Contract #3 `Domain::ScopedRecorderSuppression`, Substep-2 rollback; `LadrunoAdaptiveTransient` (HAFTOL half-step residual LTE + Hairer–Wanner PI + quantized $2^k/\sqrt{2}^k$ rungs); vectorized flat `addModalDampingForce`.

---

## `260e` — Wave 1: Windows Single OpenMP (`libiomp5md`), `mimalloc`, `/LTCG` & AVX2
- **Branch:** `wp/260e-allocator-lto-omp` from `origin/ladruno` (gated on owner D2/D4).
- **Targets:** Eliminate `vcomp140.dll` dual-pool linking (`CMakeLists.txt:774-776`, verified via `dumpbin`), wire `mimalloc` + `/GL /LTCG` + `/arch:AVX2` (or per-TU dispatch) with `/fp:precise`.

---

## Wave 2 (`260f–260j`) & Wave 3 (`260k–260o`)
Cut fresh branches from `origin/ladruno` using the contracts in `references/contracts.md`:
- **`260f`:** Full 3D/MPM-A/MCFT de-static + `LadrunoFastScatter` + persistent `LadrunoThreadedAssembly` + sparse `TransformationFE` $T^T K T$.
- **`260g`:** FP32 PARDISO (`float*` shadow) + MUMPS GMRES-IR + `LadrunoSymPrecondKrylovSOE`.
- **`260h`:** Analytic `Concrete3D`/`HS` return maps, `SAS-RK23/45`, `ASDPlastic3D` closed-form VM/DP, Eisenstat–Walker + `LadrunoAndersonNewton` + `LadrunoLBFGS`.
- **`260i`:** `LadrunoImplexAlgo` + `LadrunoSuperelement` (`OPS_LadrunoCMS`) + `H5DRMLoadPattern` v2.
- **`260j`:** `LadrunoContact` 3D Mortar bucket sort/BVH + `LadrunoContactFE` threading.
- **`260k–260o`:** `LadrunoGPStore` SoA + 7×7 Schur CPPM; `LadrunoExplicitEngine` + Clustered LTS; `LadrunoAMGSOE` (`AMGCL` + Biot Schur) + `ClusterPardiso` + weighted METIS; `LadrunoIMEX` (Hughes–Liu → Prakash–Hjelmstad); `LadrunoSpectralHex` + GPU explicit/MPM.
