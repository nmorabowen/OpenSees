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

> **Revision 2026-11-25b (Claude review of PR #1157).** Re-scoped against prior measured
> decisions (ADR-77 closed implicit Bathe and `/arch:AVX2`; see §2.3), rewrote contracts
> #1, #3 and #4 against the code on `origin/ladruno` (existing IMPL-EX API, the side effects
> of `Domain::commit()`, `setSize()` being per-subclass), corrected the de-duplication status
> of #1147/#1153/#1154/#933 and their file overlaps, pinned the MKL floor (oneMKL 2024.2.2),
> fixed the ADR-87 D4 wording, and added owner decisions D5–D7.

## 1. Purpose & Scope

Accelerate **both implicit and explicit** nonlinear mechanics and wave-propagation analyses across the `OpenSees_ladruno` fork while preserving full numerical rigor (ADR-68 gate discipline: `G-BYTE`, `G-EQUIV`, `G-OPT`). Because all `Ladruno*` elements, materials, integrators, and SOE wrappers are fork-owned, internal layouts, batching APIs, return-mapping kernels, and assembly/solver pipelines can be redesigned freely without breaking upstream OpenSees compatibility.

This ADR establishes:
1. The code-verified baseline, the prior decisions it builds on or must overturn, and live de-duplication against parallel fork WPs (`origin/ladruno` at `7e91aaf1c`, WP-257),
2. The complete bottleneck inventory across constitutive return mapping, global SOE/assembly, explicit wave dynamics, modern algorithms, and build/constraint/contact/partitioning gaps,
3. Six cross-pillar C++ interface contracts so independent WP branches compose cleanly (frozen once the owner accepts this ADR; until then, drafts),
4. The unified mixed-precision policy and explicitly rejected Tier-3 traps, and
5. Open owner decisions (`D1–D7`) before Wave 1/2/3 commits.

**Selection principle (adopted 2026-11-25).** §3's inventory is a menu, not a schedule, and its speedup figures are hypotheses until measured. `260-w0` profiles the D1 decks first; a work package is scheduled only if the phase it targets takes ≥ 10 % of a D1 run and its projected whole-run gain `1 / ((1 − p) + p / s)` is ≥ 1.10×; the owner picks the top 2–3 at a time and the menu is re-ranked after each merge with the measured gain. Rules: `_wp260_orchestration_plan.md` §0.5.

---

## 2. Grounded Baseline & Live De-Duplication (`origin/ladruno` @ `7e91aaf1c`)

### 2.1 Already Shipped in `OpenSees_ladruno` (Do Not Redo)
- **Sparse Direct & Reuse:** `PARDISOGenLinSOE`/`Solver` phase split, `-matrixType 0|1|2`, `-krylov` (`mtype 11/2`), `-deterministic`/`-cbwr` (WP-132 CNR), `LadrunoSymmetryGuard`; `MumpsParallelSOE`/`Solver` `-matrixType`/`-BLR`, `LadrunoMumpsMixed` (`-precision single -refine k`, WP-194); `LadrunoFactorReuse` (WP-188/193; reuse is keyed on `(SOE, integrator, neq, dt)` with `dt` recovered from domain time between converged steps, `LadrunoFactorReuse.cpp:487-498`); `OPS_LadrunoCMS` (`eigen` only, ADR-1000).
- **Explicit & Wave:** `CentralDifferenceLadruno` (`formTangent` mass cache ADR-67 P-NEW-1 + `-commitSolveState` P-NEW-2), `ExplicitBathe`, `CentralDifferenceSMS`/`Consistent`, HRZ lumping, `CriticalTimeStep`, `LadrunoBrick::formInertiaTerms(0)` skip (ADR-68 T7), `LadrunoCPML` (WP-199; single-phase is **already strictly SPD**; only `-up` is unsymmetric), `H5DRMLoadPattern`, `LadrunoPorousOverlay` (ADR-73 staggered fixed-stress + `-fluidUpdate explicit`), `LadrunoConstraintProjector`, `LadrunoContactFE`, `LadrunoEmbeddedNode`, and `mpm` (WP-232).
- **Constitutive:** Material `-implex` in `LadrunoSANISAND`, Lemaitre damage, uniaxial hinges; the IMPL-EX API already in the tree is `getImplexError()` (`const Vector&` in `ASDConcrete3DMaterial`/`ASDConcrete1DMaterial`/`ASDSteel1DMaterial`; `double` in `Bilin`, `ModIMKPeakOriented`, `LadrunoCohesiveHinge(Biaxial)`, consumed by `LadrunoIMKHinge.h`) plus `LadrunoImplexClock` modes and `-implexControl` trial refusal. `ASDConcrete3D`; stack-array `LadrunoSANISANDSasME` (WP-231).
- **Threading infrastructure:** `Element/Material::ladrunoThreadSafeUpdate/Form/Commit` opt-ins, `LadrunoThreadedAssembly`, `SRC/utility/LadrunoThreads.{h,cpp}`.

### 2.2 Live De-Duplication Against Open Fork PRs (Drop From WP-260)

All three PRs below are **open and stacked** (`#1154` contains `#1153` contains `#1147`); none is merged at `7e91aaf1c`. A WP-260 branch that touches any file listed here waits for that PR to merge and is cut (or rebased) after it.

| PR | State | Scope WP-260 must not redo | Files it edits that WP-260 also targets |
|---|---|---|---|
| **`WP-253` (#1147)** | open, **draft** | MPM Mode-F transfer-loop threading (`LadrunoMPMThreads.h`) | `SRC/utility/LadrunoThreads.{h,cpp}`, `LadrunoMaterialStatus.h`, `ladrunoMPM/*` (→ 260f Mode-A MPM threading, 260l) |
| **`WP-258` (#1153)** | open | Thread-safe soil materials (`DruckerPrager`, `ManzariDafalias`, `LadrunoSANISAND*`, `LadrunoHS/HSSmall`) + CPython GIL/`opserr` worker-thread fix (`PythonStream.h`) | `SRC/material/Material.h` (→ Contract #1 neighbours), `LadrunoHS.{h,cpp}` (→ 260c), `LadrunoSANISANDSasME.cpp` (→ 260h `SAS-RK23/45`), `ASDConcrete3DMaterial.cpp` (→ 260h), `FiberSection3d.{h,cpp}`, `FiberSection2d.h`, `LayeredShellFiberSection.cpp`, `ElasticSection{2,3}d.h` (→ 260f Model C threading) |
| **`WP-259` (#1154)** | open | Thread-safe concretes (`LadrunoConcrete3D`, `LadrunoRCConcrete`) | `LadrunoConcrete3D.{h,cpp}`, `LadrunoRCConcrete.{h,cpp}` (→ 260c) |

- **`WP-146` (#933):** labelled `human-gate` ("scoped, awaiting go"). ADR-77 V0 already showed bare compiler vectorization flags gain ≈ 0 on AoS loops; the true single-core vector lever is fixed-size sparse brick `B^T D B` (`34,560 → 17,280` MACs, bit-identical) + SoA batching (`LadrunoGPStore`). Folding #933's P0/P1 findings (`KMP_BLOCKTIME`, CRT heap contention, physical-core clamping) into WP-260 **is owner decision D5**; until then #933 stays open and WP-260 does not close it.

### 2.3 Prior Measured Decisions This Program Builds On (or Must Overturn)

| Prior decision | What it closed | WP-260 position |
|---|---|---|
| **ADR-77 T3 / gate G3** (§6d of `77a_c0_t0_results_2026-07-26`) | Implicit Bathe β1/β2 lane **closed with the benchmark on record**: TRBDF2 1.18–1.25×, best damped scheme 1.00–1.07× against a ≥1.3× bar; numerical damping bought **zero** iteration reduction. | `LadrunoImplicitBathe` is **not** a Wave-1 WP. It is re-opened only through **D6**: a pre-registered benchmark of the one hypothesis ADR-77 did *not* test — a shared factor across both substeps under `LadrunoFactorReuse` (Contract #2) — on the ADR-77 T3 deck plus a contact/wave deck, against the same ≥1.3× cost-per-simulated-second bar. If it misses the bar, the lane stays closed. |
| **ADR-77 V0** (§6h) | `/arch:AVX2` probe: within noise (≤1–3 % slower), zero bit drift at 1 thread (MSVC `/fp:precise` does not contract to FMA). "Do not add the flag." | `/arch:AVX2` is **dropped** from 260e. Vector speed comes from fixed-size kernels + SoA batching (260k), not a global ISA flag. A per-TU dispatch for a measured hot kernel is allowed only with its own `G-OPT` measurement. |
| **ADR-67 P1.2 → P1.2′** | Original P1.2 dropped by cross-review; replaced by P1.2′ (massless short-circuit). | 260b ports P1.1/P1.2′ only. |
| **ADR-68 T1/T2** | Standalone `LadrunoJ2` `tang_flag` (3.9 %) and closed-form plane-stress J2 (~2 %) demoted. | Not in scope (Tier-3 trap 8). |

---

## 3. Verified Bottlenecks & Scoped Solutions

### 3.1 Constitutive Return Mapping, Wrappers & Material-Point Layout
1. **Unused 6×6 tangents computed unconditionally:** `LadrunoConcrete3D` (`:587-619`), `LadrunoRCConcrete` (`:328`), `LadrunoHSKernel.h` (`:1701`), `LadrunoSANISANDSasME` (`:1475`), `ASDPlasticMaterial3D.h` (`:4430`), and `LogStrainNDMaterial` always compute 6×6 tangents during `setTrialStrain`/`setTrialF` even in explicit dynamics, `formUnbalance()`, line-search, and factor-reuse iterations. *(Note: standalone `LadrunoJ2` `tang_flag` and closed-form plane-stress J2 were demoted by ADR-68 Phase-0 T1/T2 at 3.9% / ~2% ceilings.)*
2. **`LogStrainNDMaterial::setTrialF()` 5-bug bottleneck (`SRC/material/nD/LogStrainNDMaterial.cpp:158-236`):** (a) re-runs `hencky_voigt(Be_n, epsN6)` eigensolve + `log` on committed `Be_n` every trial call (`:165`); (b) executes dense 6×6 `D0.Invert(Ce)` every trial call (`:201-205`); (c) runs `be_from_hencky_voigt` every trial call (`:210`) though only read in `commitState()`; (d) computes the spatial tangent twice in implicit mode (`:191` and `getSpatialTangentTensor` `:234`) and discards it in explicit; (e) holds 3 `static` scratch objects (`:172, :200, :206`) blocking threading.
3. **Broken `-implex` deferral in `LadrunoConcrete3D` & `LadrunoHS`:** `LadrunoConcrete3DKernel.h:3253` calls `effectiveReturn()` with `doTangent=true` and runs full implicit return + damage every Newton iteration; `LadrunoHS.cpp:624-678` runs full implicit return + tangent and calls `hssmall_predictor` (`:662`) on every trial (`Ge_implex` also varies at trial `epsComp` in `HSsmall` unless frozen/extrapolated).
4. **Expensive return-map internals:** `LadrunoConcrete3DKernel.h:2227` runs 36–48 Jacobi eigensolves via micro-FD in `damagedTangent` + $3^6$ tensor loops (`:1925`); `LadrunoHSKernel.h:1386-1455` uses FD Jacobians + 60×60 nested bisection in Regime 4/5 corners; `LadrunoSANISAND` 2nd-order `SAS-ME` can drop 2.5–5× substeps via stack-array `SAS-RK23/45`, and 19×19 `BackwardEuler_CPPM` reduces via Schur condensation to **7×7** `(sig, dGamma)`; `ASDPlasticMaterial3D` runs up to 26×26 `FullPivLU` on VM/DP and re-factors `J` for the algorithmic tangent.
5. **Unthreaded elements & SoA ceiling:** `LadrunoBrick*`, `BezierTri6/Tet10`, `LadrunoSolidShell`, `LadrunoCPML` (already hot-path static-free at `LadrunoCPML.cpp:725-830`), `LadrunoUP` (only 1 `static bool` at `:599`), Mode-A Implicit MPM (`LadrunoMPMImplicit/Point/TwoPhase`), and Model C RC fiber/MCFT (`LadrunoNDFiberSection3d`, `LadrunoMCFTStrip`, `2CLD12` = 780 s) are unthreaded. Long-term SIMD/GPU scaling requires a 64-byte-aligned SoA `LadrunoGPStore` with two-stage plastic active-set compaction. *(Note: `ops_TheActiveElement` is already `thread_local`, see the H7 note at `Domain.cpp:2611`.)*

### 3.2 Global SOE, Sparse Assembly, Linear Solvers & Nonlinear Outer Loops
1. **`MumpsSOE` & `setSize()` quadratic costs:** `MumpsSOE::addA` (`:313-397`, inherited by `MumpsParallelSOE`) still uses an $O(d^2 \cdot \text{collen})$ linear scan; `PARDISOGenLinSOE::setSize`, `MumpsSOE::setSize`, and `MumpsParallelSOE::setSize` use $O(\text{nnz}\cdot\text{deg})$ insertion sort (hitting Mode-A MPM repeatedly on `domainChanged()`); `MumpsParallelSolver` shifts `2*nnz` indices even on `JOB=3`; `MumpsParallelSOE` allocates 4 global $N$-vectors per rank with rank-0 point-to-point gather/bcast; no Element→CSR/CSC offset scatter map exists (`LadrunoFastScatter`).
2. **`LadrunoThreadedAssembly` & `Domain` churn:** Re-snapshots FEs/nodes, re-audits `ladrunoThreadSafe*()`, fork-joins OpenMP every 512 FEs, deep-copies `K`, and scatters serially (~200 ms floor at 8T).
3. **Untapped PARDISO controls & spec traps:** the comment at `PARDISOGenLinSolver.cpp:366-368` (inside the WP-132 CNR block) says `iparm[23]/[24] = 0` are "classic factorization / sequential forward-backward solve" and "already the deterministic choices"; per the Intel oneMKL `iparm` table `iparm[24]=0` selects the parallel solve (`1` forces sequential) — **cite the oneMKL 2024.2 table text in the 260a PR** and re-prove `-deterministic` bit identity before changing anything in that block. Candidate levers (all opt-in, `G-OPT`, each refused or separately proved under `-deterministic`): `iparm[23]=1|10`, `iparm[1]=3`, `iparm[24]=2`, `mtype -2` `iparm[10/12]`, `pardiso_64`. **MKL floor:** every MKL routine, macro or `iparm` value used must exist in **oneMKL 2024.2.2** (esmeralda, the oldest fork MKL), proved by `ci/check_mkl_compat.py` — "present in the local oneMKL 2025.1" is not evidence (WP-132/WP-148, #886). `iparm[27]=1` (Fortran `IPARM(28)`) requires **`float*` `a, b, x`** (a `std::vector<float>` shadow + outer FP64 GMRES-IR via `mkl_rci.h` `dfgmres`). Upgrade `LadrunoMumpsMixed` from classical IR ($\kappa \lesssim 10^7$) to Carson–Higham GMRES-IR ($\kappa \lesssim 10^{11}$).
4. **3D Fill-In Wall & Unsymmetric Penalties:** Factoring $K_{\text{sym}}=\frac{1}{2}(K+K^T)$ with symmetric PARDISO/MUMPS + right-preconditioned FGMRES (`LadrunoSymPrecondKrylovSOE`) recovers ~2× speed and −43% RAM on non-associated models (`SANISAND`, `HS`, `Concrete3D`, contact, `CPML -up`). Past ~200k–1M 3D DOFs, refresh vendored header-only `OTHER/AMGCL` → `LadrunoAMGSOE` (6 rigid-body near-nullspace modes + 3×3 block SA-AMG + prolongator reuse + White–Castelletto–Tchelepi fixed-stress Schur preconditioner for monolithic `LadrunoUP` + `CPML -up`) alongside `LadrunoClusterPardisoSOE` and optional GPU `cuDSS`.
5. **Outer-loop leaks:** `KrylovAccelerator::accelerate` (`:200-216`) runs a dead $O(k^2 N)$ loop + unblocked scratch `DGELS`; `BFGS.cpp` calls `theSOE->solve()` twice per iteration (`:346` in `solveCurrentStep`, `:442` in `BFGSUpdate`) and `break`s the entire loop on one `sdotr[i] < eps` (`:457, :507`); `AcceleratedNewton::sendSelf` drops `theReuse`. Address via Eisenstat–Walker Inexact Newton, safeguarded QR `LadrunoAndersonNewton`, and two-loop Powell-damped `LadrunoLBFGS`.

### 3.3 Explicit Dynamics, Matrix-Free Fast Paths & Wave Propagation
1. **Seven unbanked explicit costs:** `ExplicitBathe` lacks CDL P-NEW-1 mass cache, P-NEW-2 `-commitSolveState` (3 `updateDomain` passes/step), and allocates 2 `neq` vectors/step in `-lnvd`; `LadrunoBrick::formUri(0)` zeroes `stiff` (4.6 KB) and calls `getTangent()` for 1 scalar; `LadrunoBrick`/`LadrunoQuad` bulk viscosity queries `getInitialTangent()` every GP/step for $c_d$; `LadrunoShell::addNodalInertia`, `Brick20`, and `CPML` lack T7 zero-accel skips; `LadrunoConsistentRefine.h`/`LadrunoMassScaling.h`/`MPIDiagonalSOE` heap-allocate 6–7 vectors + run an extra `assembleSharedSum(ones)` MPI collective every solve (`ConsistentBlock` can also use an $O(n)$ rank-1 centroidal matvec); `LadrunoConstraintProjector` re-solves `LtML` every step instead of caching Cholesky $L$; `Element::getRayleighDampingForces` builds dense matrices for $\alpha_M$.
2. **Matrix-free `LadrunoExplicitEngine`:** Following the WP-232 Mode-F `LoadPattern` host (`explicitAnalyze $N $dt -outputEvery $K`), fuse 4 passes into 1 SIMD pass over SoA `LadrunoGPStore` inside 1 persistent OpenMP team (and optionally VRAM-resident GPU) with native flat adapters for `LadrunoCPML`, `ConstraintProjector`, `Contact` (`SOFT=1/2`), `EmbeddedNode`, and `PorousOverlay -fluidUpdate explicit`.
3. **CFL & wave dispersion:** Add power-of-two clustered LTS subcycling inside `LadrunoExplicitEngine` (momentum-conserving interface sync curing Daniel bands), `LadrunoSpectralHex` ($p=3\text{–}5$ GLL SEM for 3D wave basins: ~5 pts/$\lambda$, exact diagonal mass, $O(p^4)$ sum factorization), and pure-$M$ Hulbert–Chung / 2-stage Soares integrators on CDL's diagonal chassis (never stock `HHTGeneralizedExplicit`, whose `formEleTangent`/`formNodTangent` add `C` to the LHS at `:238-257`).

### 3.4 Modern Algorithms Missing in OpenSees
1. **`LadrunoImplicitBathe` ($\rho_\infty$-Bathe) — gated by D6 (§2.3), not a Wave-1 WP.** `TRBDF2` alternates full TR and BDF2 steps across separate `analyze()` calls with different effective-stiffness weights (`:128-146`), restarts the TR/BDF2 pair whenever `dt` changes (`:109`) or after `revertToLastStep()` (`:179`). With $\gamma_0(\rho_\infty)$ the two Bathe substeps share $c_2,c_3$, so *if* a shared factor pays (the hypothesis ADR-77 did not test) one LU/AMG setup serves both substeps. Needs Contract #2, Contract #3 and the substep driver of §4.3.
2. **`LadrunoAdaptiveTransient`:** HAFTOL half-step residual LTE + Hairer–Wanner PI + quantized hysteretic $\Delta t$ rungs ($2^k/\sqrt{2}^k$) preserving factor reuse.
3. **`LadrunoImplexAlgo`:** Global 1-SPD-solve/step driver across unified Ladruno materials with cross-step factor caching and mandatory stress-gap error control ($O(\Delta t)$ viscosity disclosure), built on the existing `getImplexError` / `LadrunoImplexClock` API (Contract #1b).
4. **`LadrunoIMEX`:** Phase A Hughes–Liu (pruned sparse LU over `I` elements only); Phase B Prakash–Hjelmstad dual-Schur multi-rate subcycling ($\Delta t_E = \Delta t_I/m$) with exact interface energy conservation.
5. **`LadrunoSuperelement` & `H5DRM` v2:** Wrap shipped `OPS_LadrunoCMS` as a nonlinear-coupled `Element`; upgrade `H5DRMLoadPattern` (header notes `H5DRMLoadPattern.h:47-59`: 3-DOF/8-node limit, velocities unused) with velocity/damping terms, Biot $u\text{–}p$, and cached matrices. Optional Go/No-Go hooks for EMC and many-query POD+ECSW.

### 3.5 Toolchain, Contact, Damping, Constraints & Partitioning Gaps
1. **Toolchain (`CMakeLists.txt`):** Unify Windows OpenMP on `libiomp5md` (removing the `vcomp` dual pool documented at `CMakeLists.txt:774-781`), evaluate `mimalloc` (CRT heap contention flagged in `#933`) and `/GL /LTCG`. **No global `/arch:AVX2`** (ADR-77 V0, §2.3).
2. **`LadrunoContact`:** Upgrade the 3D Mortar brute-force broad phase (`LadrunoContactHandler.cpp:2345-2352`; the 2D loop at `:2169-2173` mirrors it) to a slave-aware bucket sort/BVH, coordinate MPI `-reemit`, and thread `LadrunoContactFE`.
3. **Modal damping & constraints:** Vectorized `addModalDampingForce` (`IncrementalIntegrator.cpp:377-574`) + Sherman–Morrison–Woodbury low-rank tangent; sparse/block-aware `TransformationFE` $T^T K T$ (`TransformationFE.cpp:476+`); weighted `vwgts` in `MetisWrapper.cpp:215`; internal Morton/RCM cache reordering in `LadrunoGPStore`; optional static PTC globalization.

---

## 4. Six Cross-Pillar C++ Interface Contracts

All WP-260 work packages compose through six seams. They are frozen when the owner accepts this ADR; a WP that finds a contract wrong stops and amends the ADR (stop rule), it does not patch around it.

### 4.1 Contract #1a — Stress-Only Evaluation (`NDMaterial.h`, `FiniteStrainNDMaterial.h`)
- Non-breaking virtuals: `setTrialStrainStressOnly(const Vector&)`, `setTrialFStressOnly(const Matrix&)`, `setTrialFincrStressOnly(const Matrix&)`, `ladrunoSupportsStressOnly() const` (default `false`). Defaults forward to `setTrialStrain` / `setTrialF` / `setTrialFincr`, so un-upgraded materials and wrappers are unchanged.
- **Every tangent channel is covered**, not only `getTangent()`: `getTangent()`, `getSpatialTangentTensor()` (the finite-strain element channel; `LogStrainNDMaterial` reads the inner `getTangent()` there) and any wrapper that reads its inner tangent. A material that took a stress-only trial must, on any of these calls, return a tangent consistent with **that** trial — never the previous iteration's.
- **How the tangent is produced is fixed by the contract:** the stress-only trial caches what the tangent needs (return-map result, active set, trial elastic `B^e`) so materializing the tangent is a pure function of the cached trial. It never re-runs the return map from committed state.
- **Thread safety:** materialization writes only the material instance's own members, and is allowed only for materials whose `ladrunoThreadSafeForm()` is `true` with the instance owned by one thread during form (the WP-258/259 rule). An element that calls a tangent channel from more than one thread on the same instance must not use stress-only mode.
- `tangentValid_` resets on `commitState()`/`revertToLastCommit()`/`revertToStart()`; a tangent request after commit materializes from the committed state.

### 4.2 Contract #1b — IMPL-EX Control Built on the Existing API
- Builds on what ships: `getImplexError()` and `LadrunoImplexClock`. A global driver (`LadrunoImplexAlgo`) reads the material's existing `getImplexError()`; a new base virtual is added only for materials that lack one, with the **same name and a `double` return** (`virtual double getImplexError() const`) so the uniaxial hinges already conform. `ASDConcrete3D`/`ASDConcrete1D`/`ASDSteel1D` keep their `const Vector&` accessor and gain the `double` norm through a differently-named adapter (no overload-by-return-type).
- `supportsImplex()` and `setImplexMode(bool enable, double alpha = 1.0)` are added as the driver's switch. **Precedence:** a material declared with `-implex` / `-implexControl` in the deck keeps its own `LadrunoImplexClock` mode; `setImplexMode` on such a material is a no-op that warns once. The driver never runs a second extrapolation on top of a material-level one.
- `getImplexElasticDrift()` is dropped until a WP needs it with a defined meaning.

### 4.3 Contract #2 — Effective Tangent Key (`IncrementalIntegrator.h`)
- `virtual bool getEffectiveTangentKey(double &c1, double &c2, double &c3) const { return false; }`. `LadrunoFactorReuse` keys reuse on `(SOE, integrator, neq, c1, c2, c3)` when it returns `true`, falling back to today's `(SOE, integrator, neq, dt)` key when `false`.
- Integrators whose coefficients change between steps at constant `dt` (`TRBDF2`, `TRBDF3`) should implement it, so the kept factor is not reused across a TR↔BDF2 switch.

### 4.4 Contract #3 — Substep Commit Without Recorders (`Domain.h/.cpp`)
- **Recorder suppression is separate from contact augmentation.** `contactAugmenting` stays exactly as it is (it also drives `theContactDomain->commit(contactAugmenting, dT)` augment-mode updates and gates the ADR-60 `needsResort` re-emit). A new ref-counted `int suppressRecorderDepth = 0` + RAII `Domain::ScopedRecorderSuppression` skips **only** the recorder loop and the `commitTag++`.
- **Every other commit-time side effect is specified for a substep commit:**

  | `Domain::commit()` side effect | Under a substep commit |
  |---|---|
  | node/element `commitState` | runs |
  | WP-99 commit-refusal latch | runs (a refused substep fails the step) |
  | `theContactDomain->commit(..., dT)` | runs with the substep `dT`; never in augment mode |
  | ADR-60 `needsResort` re-emit | deferred to the step-final commit |
  | `LadrunoPorousOverlay::onDomainCommit` (fluid backward-Euler advance) | **skipped**; the step-final commit advances the fluid once over the full step |
  | `committedTime` | set to the substep time |
  | recorders + `commitTag++` | skipped |

- **Rollback:** after a substep commit, `revertToLastCommit()` can only return to the substep state. A composite integrator therefore must either (a) not commit substep 1 (hold it as trial state and evaluate substep 2 from it), or (b) keep a step-start snapshot it can restore for **all** materials, not only `LadrunoGPStore` ones. `LadrunoImplicitBathe` uses (a) unless D6's benchmark proves (b) is needed.
- **Substep driver:** `DirectIntegrationAnalysis::analyze` calls `newStep` → `solveCurrentStep` → `commit` once per step. A composite implicit integrator needs a driver that runs two Newton solves per step; WP-260 adds it as a fork-owned analysis (or a `LadrunoStaggeredAnalyze`-style driver), not by having the integrator call the algorithm.

### 4.5 Contract #4 — Fast CSR Scatter & Inexact-Newton Hooks (`LinearSOE.h`, `LinearSOESolver.h`)
- **The generation bump cannot be skipped:** `LinearSOE::setSize(Graph&)` becomes a non-virtual wrapper that calls the subclass's `doSetSize(Graph&)` and then increments `uint64_t patternGeneration_` (NVI). Subclasses migrate to `doSetSize`; any subclass still overriding `setSize` directly reports `supportsMappedScatter() == false`.
- **Plan key:** a scatter plan is keyed on `(const LinearSOE*, patternGeneration_)`, so a new `system` command (a fresh SOE whose generation restarts at 0) never matches an old plan.
- `locateEntry(row, col) const`, `addAMapped(const double *mColMajor, int nr, const IndexT *map, double fact = 1.0)`, `supportsMappedScatter() const` (default `false`; default `addAMapped` falls back to `addA`). `IndexT` is the SOE's own index type (`int32_t` today, `int64_t` for a `pardiso_64` / ILP64 SOE).
- `setRelativeTolerance(double eta)` (Eisenstat–Walker) lives on **`LinearSOESolver`**, where the tolerance is used; default ignores it (direct solvers).

### 4.6 Contract #5 — 64-Byte SoA `LadrunoGPStore`
- `LadrunoGPStore` holds 64-byte-aligned `eps[6][N]`, `sig[6][N]`, `Dtan[36][N]` (null in stress-only), optional lazy `statePrevCommit` (for $O(1)$ pointer-swap rollback), `syncToDomain()` (for explicit→implicit staging), `evalBatchStressOnly(W, eps, state, sig)`, and two-pass active-set compacted `plasticIdx[]`. Zero dynamic allocations inside `evalBatch`. Always retain the virtual `NDMaterial*` fallback for non-Ladruno materials.
- `*Kernel.h` headers stay plain C++ POD kernels. A `LADRUNO_HD` (`__host__ __device__`) annotation is added **only if D3 puts CUDA/HIP in scope**; until then the contract is "no allocation, no `opserr`, no statics, no virtual calls inside a kernel" so a later annotation is mechanical.

### 4.7 Contract #6 — Windows Build / OpenMP Contract (`CMakeLists.txt`)
- Single OpenMP runtime (`libiomp5md`, zero `vcomp140.dll`), verified with `dumpbin /dependents` on every Windows DLL/EXE in `dist\`; `/fp:precise` (MSVC) and `-ffp-contract=off` (GCC/Clang) for `G-BYTE` targets; `mimalloc` + `/GL /LTCG` per D4. No global `/arch:AVX2` (§2.3).

---

## 5. Unified Mixed-Precision Policy

- **Safe & High-ROI (opt-in, `G-OPT`/`G-EQUIV`):**
  - Implicit: FP32 MUMPS/PARDISO factors (`float*` shadow `A_f`, `b_f`, `x_f` for PARDISO `iparm[27]=1`) paired with **FP64 outer GMRES-IR** (Carson–Higham, valid up to $\kappa \lesssim 10^{11}$ with automatic FP64 fallback).
  - Explicit: FP32 element kinematics / elastic / CPML / J2 **force evaluation** paired with **FP64 nodal state accumulation** (`u += dt * v_half`). CPML auxiliary memory variables (the recursive-convolution state, updated every step) are **accumulated state and stay FP64**.
- **Strictly Forbidden (FP64 Only):** `LadrunoSANISAND` (`p'→p_min`), `LadrunoConcrete3D` softening, `LadrunoHS` corner/cap returns, any state accumulated over time steps (nodal state, CPML memory variables, plastic/damage internal variables), and un-accumulated long-duration wave time-stepping.

---

## 6. Tier-3 Traps Explicitly Rejected

1. **Parareal / MGRIT / PFASST** for hyperbolic wave propagation or path-dependent plasticity (phase instability + history mismatch).
2. **Dynamic h-AMR / p-adaptivity** inside OpenSees `Domain`/`LinearSOE` (use static graded `apeGmsh` meshes + mortar `LadrunoTie`).
3. **High-order implicit SDIRK / ESDIRK** ($C^0$ ground motions and yield corners cap accuracy at 2nd order).
4. **Unpreconditioned JFNK** for stiff heterogeneous SSI/RC/contact ($\kappa \sim 10^{10}\text{–}10^{14}$).
5. **Rosenbrock / W-methods** over Global IMPL-EX (noisy unsymmetric tangents + multiple solves/step).
6. **GPU offload for implicit OpenSees** before SoA/assembly fixes (PCIe + <300k DOF sparse LU loses to CPU PARDISO).
7. **Stock `HHTGeneralizedExplicit`** as a diagonal explicit wave integrator (`C` added to LHS).
8. **Global `/arch:AVX2`** (ADR-77 V0: no speed) and **standalone `LadrunoJ2` `tang_flag`** / closed-form plane-stress J2 (ADR-68 T1/T2).

---

## 7. Verification & Review Discipline (ADR-68 + ADR-87 D4)

- Every C++ WP must pass:
  - Local 5-target `Ladruno_scripts\build.bat` (verified via `ladrunoBuild` on both `dist\bin\OpenSees.exe` and `dist\bin\opensees.pyd`),
  - The static gates exactly as `.github/workflows/ladruno.yml` runs them (list in the skill's `references/fork-checklist.md`), plus `python ci/check_mkl_compat.py` against oneMKL 2024.2.2 for any MKL/PARDISO/FEAST edit,
  - ADR-68 `G-BYTE` (MD5/signbit identity), `G-EQUIV` (independent oracle + fail-closed/refusal + Linux TSan), or `G-OPT` gates; anything inside the WP-132 `-deterministic` block re-proves CNR bit identity,
  - Before/after `profiler start -perStep -deep` measurements pasted into the PR, and
  - Zone-A Ubuntu CI (`workflow_dispatch` on drafts).
- **Review.** Gemini authors; Claude runs an adversarial review on each draft PR (a bug-finding pass, logged in the PR). **ADR-87 D4** is a separate, owner-run step for family-level PRs and for accepting this ADR: the agent writes the handoff pack `Ladruno_implementation/reviews/handoff_260<x>.md` (the ADR section, key diffs, specific questions), the owner runs it through an out-of-family model, and the verdict is committed as `Ladruno_implementation/reviews/adr260<x>_verdict.md`. A Gemini-author + Claude-critic pass does not replace that unless the owner amends ADR-87.

---

## 8. Open Owner Decisions (`D1–D7`)

- **D1 — Phase-0 Benchmark Decks:** Confirm the 4–5 representative decks (suggested: 100k–500k DOF 3D nonlinear `LadrunoBrick` + `SANISAND`/`Concrete3D`, implicit `LadrunoCPML` + `H5DRM`, explicit `CDL`/`ExplicitBathe` + `CPML`, Mode-A Implicit MPM, and Model C MCFT `2CLD12`).
- **D2 — Windows CPU baseline:** Confirm `dist\` stays at the current ISA baseline (no global `/arch:AVX2`, per ADR-77 V0), with per-TU runtime dispatch allowed only for a measured hot kernel.
- **D3 — GPU / CUDA Envelope:** Target workstation GPU model/VRAM & whether CUDA (`cuDSS` + GPU explicit/MPM) is in Wave 3 scope; decides whether Contract #5 gains `LADRUNO_HD`.
- **D4 — Allocator & LTO Policy:** Approve testing `mimalloc` + `/GL /LTCG` in Wave 1 (`WP-260e`) subject to `G-BYTE` + clean Windows/Linux CI build times.
- **D5 — #933 (WP-146):** Fold #933's P0/P1 findings into WP-260 and close #933, or keep it as its own gated WP.
- **D6 — Re-open implicit Bathe?** Authorize the pre-registered shared-factor benchmark of §2.3; `LadrunoImplicitBathe` (260d's integrator part) proceeds only if it clears ADR-77's ≥1.3× bar.
- **D7 — Review pack:** Confirm ADR-87 D4 applies to this ADR's acceptance and to which WP-260 PRs (family-level), or amend ADR-87.
