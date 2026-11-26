---
title: "WP-260 orchestration plan — work packages, waves, parallelism, and dual-agent protocol"
project: Ladruno
status: draft — awaiting owner D1-D7 after Bootstrap ADR review (revised 2026-11-25b after the PR #1157 review)
program_issue: "#1156"
companion: 260_ladruno_solver_acceleration_adr
updated: 2026-11-25
---

# WP-260 Orchestration Plan (`Solver & Wave Acceleration`)

## 0. Dual-Agent Mechanics (Gemini Author + Claude Reviewer)

- **Gemini / Antigravity** discovers workspace skills in `.agents/skills/ladruno-solver-accel/` and runs in-workspace worktrees under `.agents/worktrees/wp260*` (gitignored via `.gitignore`: `.agents/*` except `.agents/skills/`), keeping all file edits, `build.bat`, static gates and `pytest` inside its worktree.
- **Claude** does not author WP-260 code. It runs the adversarial bug-finding review on each Draft PR before ready-flip. That review is **not** the ADR-87 D4 out-of-family review: D4 is the owner-run handoff pack + committed verdict (ADR §7, decision D7).
- **GitHub holds all live state** (`Ladruno_internal/ISSUES.md`): one `program` issue (#1156), board `#937`, and one Draft PR per WP (`wp/260<x>-<slug>`) with `### Handoff` comments.
- **Local Windows builds are serialized:** never run two `Ladruno_scripts\build.bat` builds at once on the same workstation; run static gates and coding in parallel, then queue local 5-target builds one worktree at a time before firing parallel Zone-A Ubuntu CI runs.

---

## 1. Work Package Map & File Ownership

Ownership is disjoint **between WP-260 packages**. It is **not** disjoint from the open PRs #1147/#1153/#1154 (ADR §2.2): a package that owns a file those PRs edit is marked *blocked-by* and is cut from `origin/ladruno` only after that PR merges.

### Wave 0 — Bootstrap & Phase-0 Profiler Baseline
| WP | Branch | Content | Dependencies | Gate |
|---|---|---|---|---|
| **260** | `wp/260-solver-accel-program` | Program ADR, this plan, `ladruno-solver-accel` skill, `.gitignore`/`.agents/` setup, `AGENTS.md` row, quirk-pointer lint widened to `.agents/skills/`, 3 quirk fragments | None (merges first) | Static gates green; Claude review; owner D1–D7; ADR-87 D4 verdict if D7 says so |
| **260-w0** | `wp/260w0-phase0-profiles` | `profiler start -perStep -deep` sweep on the D1 decks (`origin/ladruno` + `2CLD12`); exact phase % ceilings posted on the `program` issue (no `SRC/` edits) | `260` + D1 | Baseline table posted on `program` issue |

### Wave 1 — Parallel Surgical Wins (Mostly `G-BYTE`)
| WP | Branch | Scope & File Ownership | Blocked by | Speedup Target | Gate |
|---|---|---|---|---|---|
| **260a** | `wp/260a-soe-scatter-hygiene` | `MumpsSOE::addA` binary search + skip `rowA++/--` on `JOB=3` + `std::sort` `setSize()` (incl. local-vertex `MumpsParallelSOE`); `PARDISOGenLinSOE::setSize` `std::sort`; `KrylovAccelerator.cpp:200-216` dead loop; fix the `PARDISOGenLinSolver.cpp:366-368` comment (re-prove `-deterministic`) & expose opt-in `iparm[23]/iparm[1]/iparm[24]=2/iparm[10/12]` (refused or proved under CNR; oneMKL 2024.2.2 floor via `check_mkl_compat.py`); `AcceleratedNewton::sendSelf` reuse serialization | — | 3–20× `Mumps`/`setSize`, 1.3–2× PARDISO factor (to be measured), faster `KrylovNewton` | `G-BYTE` (default) + `G-OPT` (flags) + MKL compat |
| **260b** | `wp/260b-explicit-wave-hygiene` | `ExplicitBathe` CDL mass-cache + `-commitSolveState` + hoisted `-lnvd`; ADR-67 P1.1/P1.2′ consistent-SMS & `MPIDiagonalSOE` alloc/halo hoist + rank-1 $O(n)$ centroidal matvec + `CriticalTimeStep` massless check; T7 inertia skips (`LadrunoElementShell.h`, `Brick20`, `CPML`); `LadrunoBrick::formUri`/BV $c_d$ cache; `ConstraintProjector` Cholesky $L$; global $O(N)$ $\alpha_M$ (`TransientIntegrator.cpp`); `LadrunoCPML` OpenMP opt-in | — | 1.4–2.2× explicit (`CDL`/`ExplicitBathe`), 4–8× consistent-SMS PCG | `G-BYTE` / `G-EQUIV` |
| **260c1** | `wp/260c1-logstrain-stressonly` | `LogStrainNDMaterial::setTrialF()` 5-bug fix; Contract #1a on `NDMaterial.h`/`FiniteStrainNDMaterial.h` + `LogStrainNDMaterial`/`LogStrainKernel.h` (all tangent channels incl. `getSpatialTangentTensor`); `LadrunoUP` OpenMP opt-in after guarding the `:599` static | — | 3–6× finite-strain GP | `G-BYTE` / `G-EQUIV` |
| **260c2** | `wp/260c2-concrete-hs-implex` | Contract #1a on `LadrunoConcrete3D*`/`LadrunoHS*`/`LadrunoRCConcrete*`; `doTangent && !mp.implex` (`LadrunoConcrete3DKernel.h:3253`); commit-deferred `-implex` in `Concrete3D` & `HS` (`Ge_implex` freeze) on the existing `LadrunoImplexClock` (Contract #1b) | **#1153, #1154** | 1.5–4× unused-tangent skips; 4–15× `-implex` | `G-BYTE` / `G-EQUIV` |
| **260d** | `wp/260d-adaptive-transient` | `LadrunoAdaptiveTransient` (HAFTOL half-step LTE + quantized $2^k/\sqrt{2}^k$ rungs) + Contract #2 `getEffectiveTangentKey` (incl. `TRBDF2`/`TRBDF3`) + vectorized `addModalDampingForce`. **Owned:** `SRC/analysis/integrator/IncrementalIntegrator.{h,cpp}`, `TRBDF2.{h,cpp}`, `TRBDF3.{h,cpp}`, new `LadrunoAdaptiveTransient*`, `SRC/analysis/algorithm/equiSolnAlgo/LadrunoFactorReuse*`. `LadrunoImplicitBathe` + Contract #3 are **not** in 260d: they wait for D6 (260d-B below) | — | 2–5× adaptive runs (to be measured) | `G-EQUIV` / `G-OPT` |
| **260e** | `wp/260e-allocator-lto-omp` | Unify Windows OpenMP on `libiomp5md` (remove the `vcomp` dual pool, `CMakeLists.txt:774-781`, verified via `dumpbin`); evaluate/wire `mimalloc` + `/GL /LTCG` (per D4). **No `/arch:AVX2`** (ADR-77 V0) | D4 | to be measured (#933's heap-contention finding is the hypothesis) | `G-BYTE` + `dumpbin` check |

**Gated (not scheduled until its decision):**
- **260d-B** `wp/260dB-implicit-bathe` — D6 benchmark first (shared factor across substeps under `LadrunoFactorReuse`, on the ADR-77 T3 deck + a contact/wave deck, ≥1.3× bar). Only if it clears: `LadrunoImplicitBathe` + Contract #3 (`ScopedRecorderSuppression`, commit-side-effect table, substep driver, no-commit substep 1). Owns `SRC/domain/domain/Domain.{h,cpp}` and the new integrator/driver files.

### Wave 2 — Core Kernels, Solvers & SSI Upgrades (3–8 Weeks)
- **`260f` (3D/MPM/MCFT Threading + Persistent Assembly):** De-static & opt-in `LadrunoBrick*`, `Bezier*`, `LadrunoSolidShell`, Mode-A Implicit MPM, Model C `LadrunoNDFiberSection3d`/`MCFTStrip`; `LadrunoFastScatter` (Contract #4, NVI `setSize`/`doSetSize`) + persistent `LadrunoThreadedAssembly` + colored scatter + sparse `TransformationFE` $T^T K T$ + `Domain` snapshot cache. *Blocked by #1147 (MPM, `LadrunoThreads`) and #1153 (fiber sections).*
- **`260g` (Mixed-Precision GMRES-IR & Sym-Preconditioned Krylov):** FP32 PARDISO (`iparm[27]=1` `float*` shadow) + MKL RCI `dfgmres` FP64 GMRES-IR; `LadrunoMumpsMixed` GMRES-IR upgrade; `LadrunoSymPrecondKrylovSOE` ($K_{\text{sym}}$ factor + FGMRES). oneMKL 2024.2.2 floor via `check_mkl_compat.py`.
- **`260h` (Analytic Return Maps & Outer Loop):** `Concrete3D` analytic `damagedTangent` & Voigt lift + Cardano/Smith 3×3 eigensolver; `LadrunoHS` analytic R4/R5 Jacobian & IFT Newton–Brent; `SANISAND` `SAS-RK23/45`; `ASDPlastic3D` closed-form VM/DP; Eisenstat–Walker Inexact Newton (Contract #4 `LinearSOESolver::setRelativeTolerance`) + `LadrunoAndersonNewton` + two-loop Powell-damped `LadrunoLBFGS`; optional modal Woodbury & static PTC. *Blocked by #1153/#1154 (SasME, ASDConcrete3D, Concrete3D, HS).*
- **`260i` (Global IMPL-EX + `LadrunoSuperelement` + `H5DRM` v2):** `LadrunoImplexAlgo`/`LadrunoImplexLinear` on Contract #1b; `LadrunoSuperelement` wrapping `OPS_LadrunoCMS`; velocity-consistent & Biot $u\text{–}p$ `H5DRMLoadPattern` v2.
- **`260j` (`LadrunoContact` Acceleration):** Slave-aware corner-span bucket sort / BVH for 3D Mortar (`LadrunoContactHandler.cpp:2345-2352`), MPI `-reemit`, threaded `LadrunoContactFE`, ALM active-set audit.

### Wave 3 — Transformative Engines & HPC Frontier (2–4+ Months)
- **`260k` (`LadrunoGPStore` SoA + 7×7 Schur SANISAND CPPM):** 64-byte SoA + two-pass active-set compaction + Morton/RCM internal cache ordering.
- **`260l` (`LadrunoExplicitEngine` + Clustered LTS):** WP-232 Mode-F `LoadPattern` flat explicit engine + flat `CPML`/`Contact`/`EmbeddedNode`/`PorousOverlay` + power-of-two LTS subcycling + FP32-force/FP64-state (CPML memory variables FP64).
- **`260m` (`LadrunoAMGSOE` + Monolithic Biot Schur + Cluster Direct):** Refreshed Boost-free `OTHER/AMGCL` (6 rigid-body near-nullspace modes + 3×3 block SA-AMG + monolithic `LadrunoUP` fixed-stress Schur); `LadrunoClusterPardisoSOE`; weighted METIS `vwgts`; optional `cuDSS` (per D3).
- **`260n` (`LadrunoIMEX` Spatial Partitioning):** Phase A Hughes–Liu pruned-LU IMEX → Phase B Prakash–Hjelmstad dual-Schur multi-rate subcycling atop `260l`.
- **`260o` (`LadrunoSpectralHex` + Pure-$M$ Integrators + GPU Explicit/MPM):** $p=3\text{–}5$ GLL SEM + pure-$M$ Hulbert–Chung / 2-stage Soares + `LADRUNO_HD` CUDA/HIP backend (only if D3 puts GPU in scope).

---

## 2. Landing & Review Cadence

1. Merge **`WP-260`** first after Claude review + owner D1–D7 answers (+ the D4 verdict file if D7 requires it).
2. Run **`260-w0`** baseline profiling and post numbers to the `program` issue. Then launch the unblocked Wave-1 packages (**`260a`**, **`260b`**, **`260c1`**, **`260d`**, and **`260e`** after D4) in parallel worktrees under `.agents/worktrees/`. **`260c2`** starts after #1153/#1154 merge; **`260d-B`** only after D6.
3. Serialize local `build.bat` per worktree → run the `ladruno.yml` static gates (+ `check_mkl_compat.py` for MKL edits) + local pytest + `-perStep`/`-deep` benchmarks → push Draft PRs → run Zone-A Ubuntu CI in parallel → Claude adversarial review → address findings → owner merge commit.
