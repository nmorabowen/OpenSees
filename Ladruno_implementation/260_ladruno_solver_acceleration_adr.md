---
title: "WP-260 — Solver & Wave Acceleration for OpenSees Ladruno: Architecture, Contracts, and Phased Roadmap"
project: Ladruno
status: draft — NO C++ code (bootstrap ADR for architectural review)
priority: high
owner: nmora
wp: "WP-260 (program #1156)"
tags:
  - adr
  - solvers
  - explicit
  - implicit
  - wave-propagation
  - threading
  - return-mapping
aliases:
  - WP-260
  - Solver & Wave Acceleration
  - FastSolve
updated: 2026-11-25
---

# WP-260 — Solver & Wave Acceleration Program (`FastSolve`)

## 1. Purpose & Scope

Accelerate **both implicit and explicit** nonlinear mechanics and wave-propagation analyses across the `OpenSees_ladruno` fork while preserving full numerical rigor (ADR-68 gate discipline: `G-BYTE`, `G-EQUIV`, `G-OPT`). Because all `Ladruno*` elements, materials, integrators, and SOE wrappers are fork-owned, internal layouts, batching APIs, return-mapping kernels, and assembly/solver pipelines can be redesigned freely without breaking upstream OpenSees compatibility.

This ADR establishes:
1. The code-verified baseline and live de-duplication against parallel fork WPs (`origin/ladruno` at `7e91aaf1c`, WP-257),
2. The complete bottleneck inventory across constitutive return mapping, global SOE/assembly, explicit wave dynamics, modern algorithms, and build/constraint/contact/partitioning gaps,
3. Six locked cross-pillar C++ interface contracts so independent WP branches compose cleanly,
4. The unified mixed-precision policy and explicitly rejected Tier-3 traps, and
5. Open owner decisions (`D1–D4`) before Wave 2/3 commits.

---

## 2. Grounded Baseline & Live De-Duplication (`origin/ladruno` @ `7e91aaf1c`)

### Already Shipped in `OpenSees_ladruno` (Do Not Redo)
- **Sparse Direct & Reuse:** `PARDISOGenLinSOE`/`Solver` phase split, `-matrixType 0|1|2`, `-krylov` (`mtype 11/2`), `-deterministic`/`-cbwr`, `LadrunoSymmetryGuard`; `MumpsParallelSOE`/`Solver` `-matrixType`/`-BLR`, `LadrunoMumpsMixed` (`-precision single -refine k`, WP-194); `LadrunoFactorReuse` (WP-188/193); `OPS_LadrunoCMS` (`eigen` only, ADR-1000).
- **Explicit & Wave:** `CentralDifferenceLadruno` (`formTangent` mass cache ADR-67 P-NEW-1 + `-commitSolveState` P-NEW-2), `ExplicitBathe`, `CentralDifferenceSMS`/`Consistent`, HRZ lumping, `CriticalTimeStep`, `LadrunoBrick::formInertiaTerms(0)` skip (ADR-68 T7), `LadrunoCPML` (WP-199; single-phase is **already strictly SPD**; only `-up` is unsymmetric), `H5DRMLoadPattern`, `LadrunoPorousOverlay` (ADR-73 staggered fixed-stress + `-fluidUpdate explicit`), `LadrunoConstraintProjector`, `LadrunoContactFE`, `LadrunoEmbeddedNode`, and `mpm` (WP-232).
- **Constitutive:** Material `-implex` in `LadrunoSANISAND`, Lemaitre damage, uniaxial hinges + `LadrunoImplexClock`, `ASDConcrete3D`; stack-array `LadrunoSANISANDSasME` (WP-231).

### Live De-Duplication Against Pending/Active Fork WPs (Drop From WP-260)
- **`WP-253` (#1147):** MPM Mode-F transfer-loop threading (`LadrunoMPMThreads.h`) — already built.
- **`WP-258` (#1153):** Thread-safe soil materials (`DruckerPrager`, `ManzariDafalias`, `LadrunoSANISAND*`, `LadrunoHS`) + CPython GIL/`opserr` worker-thread crash fix (`PythonStream.h`) — ready on GitHub.
- **`WP-259` (#1154):** Thread-safe concretes (`LadrunoConcrete3D`, `LadrunoRCConcrete`) — ready on GitHub.
- **`WP-146` (#933):** Harmonized — ADR-77 V0 already showed bare compiler vectorization flags gain ≈ 0 on AoS loops; the true single-core vector lever is fixed-size sparse brick `B^T D B` (`34,560 → 17,280` MACs, bit-identical) + SoA batching (`LadrunoGPStore`). `WP-146` P0/P1 findings (`KMP_BLOCKTIME`, CRT heap contention, physical-core clamping) fold directly into WP-260.

---

## 3. Verified Bottlenecks & Scoped Solutions

### 3.1 Constitutive Return Mapping, Wrappers & Material-Point Layout
1. **Unused 6×6 tangents computed unconditionally:** `LadrunoConcrete3D` (`:587-619`), `LadrunoRCConcrete` (`:328`), `LadrunoHSKernel.h` (`:1701`), `LadrunoSANISANDSasME` (`:1475`), `ASDPlasticMaterial3D.h` (`:4430`), and `LogStrainNDMaterial` always compute 6×6 tangents during `setTrialStrain`/`setTrialF` even in explicit dynamics, `formUnbalance()`, line-search, and factor-reuse iterations. *(Note: standalone `LadrunoJ2` `tang_flag` and closed-form plane-stress J2 were demoted by ADR-68 Phase-0 T1/T2 at 3.9% / ~2% ceilings.)*
2. **`LogStrainNDMaterial::setTrialF()` 5-bug bottleneck (`SRC/material/nD/LogStrainNDMaterial.cpp:158-236`):** (a) re-runs `hencky_voigt(Be_n, epsN6)` eigensolve + `log` on committed `Be_n` every trial call; (b) executes dense 6×6 `D0.Invert(Ce)` every trial call; (c) runs `be_from_hencky_voigt` every trial call though only read in `commitState()`; (d) computes spatial tangent twice in implicit mode and discards it in explicit; (e) holds 3 `static` scratch objects (`:172, :200, :206`) blocking threading.
3. **Broken `-implex` deferral in `LadrunoConcrete3D` & `LadrunoHS`:** `LadrunoConcrete3DKernel.h:3253` calls `effectiveReturn()` with `doTangent=true` and runs full implicit return + damage every Newton iteration; `LadrunoHS.cpp:624-678` runs full implicit return + tangent and calls `hssmall_predictor` twice (`Ge_implex` also varies at trial `epsComp` in `HSsmall` unless frozen/extrapolated).
4. **Expensive return-map internals:** `LadrunoConcrete3DKernel.h:2227` runs 36–48 Jacobi eigensolves via micro-FD in `damagedTangent` + $3^6$ tensor loops (`:1925`); `LadrunoHSKernel.h:1386-1455` uses FD Jacobians + 60×60 nested bisection in Regime 4/5 corners; `LadrunoSANISAND` 2nd-order `SAS-ME` can drop 2.5–5× substeps via stack-array `SAS-RK23/45`, and 19×19 `BackwardEuler_CPPM` reduces via Schur condensation to **7×7** `(sig, dGamma)`; `ASDPlasticMaterial3D` runs up to 26×26 `FullPivLU` on VM/DP and re-factors `J` for the algorithmic tangent.
5. **Unthreaded elements & SoA ceiling:** `LadrunoBrick*`, `BezierTri6/Tet10`, `LadrunoSolidShell`, `LadrunoCPML` (already hot-path static-free at `LadrunoCPML.cpp:725-830`!), `LadrunoUP` (only 1 `static bool` at `:599`), Mode-A Implicit MPM (`LadrunoMPMImplicit/Point/TwoPhase`), and Model C RC fiber/MCFT (`LadrunoNDFiberSection3d`, `LadrunoMCFTStrip`, `2CLD12` = 780 s) are unthreaded. Long-term SIMD/GPU scaling requires a 64-byte-aligned SoA `LadrunoGPStore` with two-stage plastic active-set compaction. *(Note: `ops_TheActiveElement` is already `thread_local` at `Domain.cpp:2611`.)*

### 3.2 Global SOE, Sparse Assembly, Linear Solvers & Nonlinear Outer Loops
1. **`MumpsSOE` & `setSize()` quadratic costs:** `MumpsSOE::addA` (`:313-397`, inherited by `MumpsParallelSOE`) still uses an $O(d^2 \cdot \text{collen})$ linear scan; `PARDISOGenLinSOE::setSize`, `MumpsSOE::setSize`, and `MumpsParallelSOE::setSize` use $O(\text{nnz}\cdot\text{deg})$ insertion sort (hitting Mode-A MPM repeatedly on `domainChanged()`); `MumpsParallelSolver` shifts `2*nnz` indices even on `JOB=3`; `MumpsParallelSOE` allocates 4 global $N$-vectors per rank with rank-0 point-to-point gather/bcast; no Element→CSR/CSC `int32_t` offset scatter map exists (`LadrunoFastScatter`).
2. **`LadrunoThreadedAssembly` & `Domain` churn:** Re-snapshots FEs/nodes, re-audits `ladrunoThreadSafe*()`, fork-joins OpenMP every 512 FEs, deep-copies `K`, and scatters serially (~200 ms floor at 8T).
3. **Untapped PARDISO controls & spec traps:** `PARDISOGenLinSolver.cpp:367` falsely says `iparm[24]=0` is sequential (per Intel oneMKL spec `0` is already parallel; `1` forces sequential!). Real levers: `iparm[23]=1|10`, `iparm[1]=3`, `iparm[24]=2`, `mtype -2` `iparm[10/12]`, `pardiso_64` (present in local oneMKL 2025.1). `iparm[27]=1` (Fortran `IPARM(28)`) requires **`float*` `a, b, x`** (needs a `std::vector<float>` shadow + outer FP64 GMRES-IR via `mkl_rci.h` `dfgmres`). Upgrade `LadrunoMumpsMixed` from classical IR ($\kappa \lesssim 10^7$) to Carson–Higham GMRES-IR ($\kappa \lesssim 10^{11}$).
4. **3D Fill-In Wall & Unsymmetric Penalties:** Factoring $K_{\text{sym}}=\frac{1}{2}(K+K^T)$ with symmetric PARDISO/MUMPS + right-preconditioned FGMRES (`LadrunoSymPrecondKrylovSOE`) recovers ~2× speed and −43% RAM on non-associated models (`SANISAND`, `HS`, `Concrete3D`, contact, `CPML -up`). Past ~200k–1M 3D DOFs, refresh vendored header-only `OTHER/AMGCL` → `LadrunoAMGSOE` (6 rigid-body near-nullspace modes + 3×3 block SA-AMG + prolongator reuse + White–Castelletto–Tchelepi fixed-stress Schur preconditioner for monolithic `LadrunoUP` + `CPML -up`) alongside `LadrunoClusterPardisoSOE` and optional GPU `cuDSS`.
5. **Outer-loop leaks:** `KrylovAccelerator::accelerate` (`:200-216`) runs a dead $O(k^2 N)$ loop + unblocked scratch `DGELS`; `BFGS.cpp` calls `theSOE->solve()` twice per iteration (`:346, :442`) and `break`s the entire loop on one `sdotr[i] < eps` (`:457, :507`); `AcceleratedNewton::sendSelf` drops `theReuse`. Address via Eisenstat–Walker Inexact Newton, safeguarded QR `LadrunoAndersonNewton`, and two-loop Powell-damped `LadrunoLBFGS`.

### 3.3 Explicit Dynamics, Matrix-Free Fast Paths & Wave Propagation
1. **Seven unbanked explicit costs:** `ExplicitBathe` lacks CDL P-NEW-1 mass cache, P-NEW-2 `-commitSolveState` (3 `updateDomain` passes/step), and allocates 2 `neq` vectors/step in `-lnvd`; `LadrunoBrick::formUri(0)` zeroes `stiff` (4.6 KB) and calls `getTangent()` for 1 scalar; `LadrunoBrick`/`LadrunoQuad` bulk viscosity queries `getInitialTangent()` every GP/step for $c_d$; `LadrunoShell::addNodalInertia`, `Brick20`, and `CPML` lack T7 zero-accel skips; `LadrunoConsistentRefine.h`/`LadrunoMassScaling.h`/`MPIDiagonalSOE` heap-allocate 6–7 vectors + run an extra `assembleSharedSum(ones)` MPI collective every solve (`ConsistentBlock` can also use an $O(n)$ rank-1 centroidal matvec); `LadrunoConstraintProjector` re-solves `LtML` every step instead of caching Cholesky $L$; `Element::getRayleighDampingForces` builds dense matrices for $\alpha_M$.
2. **Matrix-free `LadrunoExplicitEngine`:** Following the WP-232 Mode-F `LoadPattern` host (`explicitAnalyze $N $dt -outputEvery $K`), fuse 4 passes into 1 SIMD pass over SoA `LadrunoGPStore` inside 1 persistent OpenMP team (and optionally VRAM-resident GPU) with native flat adapters for `LadrunoCPML`, `ConstraintProjector`, `Contact` (`SOFT=1/2`), `EmbeddedNode`, and `PorousOverlay -fluidUpdate explicit`.
3. **CFL & wave dispersion:** Add power-of-two clustered LTS subcycling inside `LadrunoExplicitEngine` (momentum-conserving interface sync curing Daniel bands), `LadrunoSpectralHex` ($p=3\text{–}5$ GLL SEM for 3D wave basins: ~5 pts/$\lambda$, exact diagonal mass, $O(p^4)$ sum factorization), and pure-$M$ Hulbert–Chung / 2-stage Soares integrators on CDL's diagonal chassis (never stock `HHTGeneralizedExplicit`, which adds `C` to the LHS at `:238-257`).

### 3.4 Modern Algorithms Missing in OpenSees
1. **`LadrunoImplicitBathe` ($\rho_\infty$-Bathe):** Fixes `TRBDF2.cpp` (`:109`, `:128-146`, `:179`); optimal $\gamma_0(\rho_\infty)$ keeps $c_2,c_3$ identical across both substeps → **1 LU/AMG setup shared across substeps**, zero overshoot, superior wave dispersion.
2. **`LadrunoAdaptiveTransient`:** HAFTOL half-step residual LTE + Hairer–Wanner PI + quantized hysteretic $\Delta t$ rungs ($2^k/\sqrt{2}^k$) preserving factor reuse.
3. **`LadrunoImplexAlgo`:** Global 1-SPD-solve/step driver across unified Ladruno materials with cross-step factor caching and mandatory stress-gap error control ($O(\Delta t)$ viscosity disclosure).
4. **`LadrunoIMEX`:** Phase A Hughes–Liu (pruned sparse LU over `I` elements only); Phase B Prakash–Hjelmstad dual-Schur multi-rate subcycling ($\Delta t_E = \Delta t_I/m$) with exact interface energy conservation.
5. **`LadrunoSuperelement` & `H5DRM` v2:** Wrap shipped `OPS_LadrunoCMS` as a nonlinear-coupled `Element`; upgrade `H5DRMLoadPattern.h:47-59` with velocity/damping terms, Biot $u\text{–}p$, and cached matrices. Optional Go/No-Go hooks for EMC and many-query POD+ECSW.

### 3.5 Toolchain, Contact, Damping, Constraints & Partitioning Gaps
1. **Toolchain (`CMakeLists.txt`):** Unify Windows OpenMP on `libiomp5md` (eliminating `vcomp` dual-pool spin contention), evaluate `mimalloc` (solving CRT heap contention flagged in `#933`), `/GL /LTCG`, and `/arch:AVX2` (or per-TU dispatch).
2. **`LadrunoContact`:** Upgrade 3D Mortar brute-force broad phase (`LadrunoContactHandler.cpp:2171,2345-2351`) to slave-aware bucket sort/BVH, coordinate MPI `-reemit`, and thread `LadrunoContactFE`.
3. **Modal damping & constraints:** Vectorized `addModalDampingForce` + Sherman–Morrison–Woodbury low-rank tangent; sparse/block-aware `TransformationFE` $T^T K T$ (`TransformationFE.cpp:476+`); weighted `vwgts` in `MetisWrapper.cpp:215`; internal Morton/RCM cache reordering in `LadrunoGPStore`; optional static PTC globalization.

---

## 4. Six Locked Cross-Pillar C++ Interface Contracts

All WP-260 work packages compose through six frozen seams:

1. **Stress-Only & Unified IMPL-EX Hooks (`NDMaterial.h`, `FiniteStrainNDMaterial.h`):**
   - Non-breaking virtuals: `setTrialStrainStressOnly(const Vector&)`, `setTrialFStressOnly(const Matrix&)`, `ladrunoSupportsStressOnly()`. Lazy tangent materialization inside `getTangent()` if `!tangentValid_` (never silent stale tangents).
   - Unified IMPL-EX contract: `supportsImplex()`, `setImplexMode(bool enable, double alpha = 1.0)`, `getImplexStressGap(double &gapNorm2, double &refNorm2)`, `getImplexElasticDrift()`.
2. **Effective Tangent Key (`IncrementalIntegrator.h`):**
   - `virtual bool getEffectiveTangentKey(double &c1, double &c2, double &c3) const { return false; }`. `LadrunoFactorReuse` keys reuse on `(c1, c2, c3)` when true (falling back to `dt` when false).
3. **Commit-vs-Recorder Decoupling (`Domain.h/.cpp`):**
   - Generalize `contactAugmenting` into ref-counted `int suppressRecorderDepth = 0` + RAII `Domain::ScopedRecorderSuppression` for composite substeps and subcycled explicit/IMEX domains.
4. **Fast CSR Scatter & Inexact-Newton Hooks (`LinearSOE.h`):**
   - `uint64_t patternGeneration_` incremented in `setSize()`; `locateEntry(row, col)`, `addAMapped(...)`, `supportsMappedScatter()`; `setRelativeTolerance(double eta)` for Eisenstat–Walker Inexact Newton.
5. **64-Byte SoA `LadrunoGPStore` & `LADRUNO_HD` Kernels:**
   - Pure POD `*Kernel.h` headers compiled `__host__ __device__`; `LadrunoGPStore` holds 64-byte-aligned `eps[6][N]`, `sig[6][N]`, `Dtan[36][N]` (null in stress-only), optional lazy `statePrevCommit` (for $O(1)$ pointer-swap rollback on Substep-2 failure in `LadrunoImplicitBathe`), `syncToDomain()` (for explicit→implicit staging), `evalBatchStressOnly(W, eps, state, sig)`, and two-pass active-set compacted `plasticIdx[]`.
6. **Windows Build / OpenMP / ISA Contract (`CMakeLists.txt`):**
   - Single OpenMP runtime (`libiomp5md`, zero `vcomp140.dll`), `/fp:precise` / `-ffp-contract=off` for `G-BYTE`, plus AVX2 + allocator (`mimalloc`) & LTO controls.

---

## 5. Unified Mixed-Precision Policy

- **Safe & High-ROI:**
  - Implicit: FP32 MUMPS/PARDISO factors (`float*` shadow `A_f` for PARDISO `iparm[27]=1`) paired with **FP64 outer GMRES-IR** (Carson–Higham, valid up to $\kappa \lesssim 10^{11}$ with automatic FP64 fallback).
  - Explicit: FP32 element kinematics / elastic / CPML / J2 forces paired with **FP64 nodal state accumulation** (`u += dt * v_half`, 1.7–2.0× explicit speedup).
- **Strictly Forbidden (FP64 Only):** `LadrunoSANISAND` (`p'→p_min`), `LadrunoConcrete3D` softening, `LadrunoHS` corner/cap returns, and un-accumulated long-duration wave time-stepping.

---

## 6. Tier-3 Traps Explicitly Rejected

1. **Parareal / MGRIT / PFASST** for hyperbolic wave propagation or path-dependent plasticity (phase instability + history mismatch).
2. **Dynamic h-AMR / p-adaptivity** inside OpenSees `Domain`/`LinearSOE` (use static graded `apeGmsh` meshes + mortar `LadrunoTie`).
3. **High-order implicit SDIRK / ESDIRK** ($C^0$ ground motions and yield corners cap accuracy at 2nd order).
4. **Unpreconditioned JFNK** for stiff heterogeneous SSI/RC/contact ($\kappa \sim 10^{10}\text{–}10^{14}$).
5. **Rosenbrock / W-methods** over Global IMPL-EX (noisy unsymmetric tangents + multiple solves/step).
6. **GPU offload for implicit OpenSees** before SoA/assembly fixes (PCIe + <300k DOF sparse LU loses to CPU PARDISO).
7. **Stock `HHTGeneralizedExplicit`** as a diagonal explicit wave integrator (`C` added to LHS).

---

## 7. Verification & Review Discipline (ADR-68 + ADR-87 D4)

- Every C++ WP must pass:
  - Local 5-target `Ladruno_scripts\build.bat` (verified via `ops.ladrunoBuild()` hash on both `OpenSees.exe` and `dist\bin\opensees.pyd`),
  - All static CI gates (`check_classtags`, `check_manifest`, `check_ladruno_commands`, `check_ledger_fragments`, `check_quirk_patterns`, `check_duplicate_defs`, `ruff F811`, `stamp_headers --check`),
  - ADR-68 `G-BYTE` (MD5/signbit identity), `G-EQUIV` (independent oracle + fail-closed/refusal + Linux TSan), or `G-OPT` gates,
  - Before/after `-perStep` + `-deep` profiler measurements pasted into the PR, and
  - Zone-A Ubuntu CI (`workflow_dispatch`) + out-of-family **Claude review** (Gemini author + Claude critic satisfies ADR-87 D4).

---

## 8. Open Owner Decisions (`D1–D4`)

- **D1 — Phase-0 Benchmark Decks:** Confirm the 4–5 representative decks (suggested: 100k–500k DOF 3D nonlinear `LadrunoBrick` + `SANISAND`/`Concrete3D`, implicit `LadrunoCPML` + `H5DRM`, explicit `CDL`/`ExplicitBathe` + `CPML`, Mode-A Implicit MPM, and Model C MCFT `2CLD12`).
- **D2 — Windows CPU ISA Baseline:** Can Windows `dist\` require AVX2 (Haswell 2013+), or do we need runtime AVX2 dispatch?
- **D3 — GPU / CUDA Envelope:** Target workstation GPU model/VRAM & whether CUDA (`cuDSS` + GPU explicit/MPM) is in Wave 3 scope.
- **D4 — Allocator & LTO Policy:** Approve testing `mimalloc` + `/GL /LTCG` in Wave 1 (`WP-260e`) subject to `G-BYTE` + clean Windows/Linux CI build times.
