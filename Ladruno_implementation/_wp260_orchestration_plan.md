---
title: "WP-260 orchestration plan — work packages, waves, parallelism, and dual-agent protocol"
project: Ladruno
status: draft — awaiting owner D1-D4 confirmation after Bootstrap ADR review
program_issue: "#1156"
companion: 260_ladruno_solver_acceleration_adr
updated: 2026-11-25
---

# WP-260 Orchestration Plan (`Solver & Wave Acceleration`)

## 0. Dual-Agent Mechanics (Gemini Author + Claude Reviewer)

- **Gemini / Antigravity** discovers workspace skills in `.agents/skills/ladruno-solver-accel/` and runs in-workspace worktrees under `.agents/worktrees/wp260*` (gitignored via `.gitignore`: `.agents/*` except `.agents/skills/`), keeping all file edits, `build.bat`, `ci/check_*.py`, and `pytest` inside the sandbox.
- **Claude** stays 100% undisturbed (`CLAUDE.md`, `.claude/*`, `.claude/worktrees/*` untouched) and serves as the out-of-family critic (ADR-87 D4) on each Draft PR before ready-flip.
- **GitHub holds all live state** (`Ladruno_internal/ISSUES.md`): one `program` issue (`Solver & Wave Acceleration (WP-260…)`), board `#937`, and one Draft PR per WP (`wp/260<x>-<slug>`) with `### Handoff` comments.
- **Local Windows builds are serialized:** Never run two `Ladruno_scripts\build.bat` builds simultaneously on the same workstation; run static checks and coding in parallel, then queue local 5-target builds one worktree at a time before firing parallel Zone-A Ubuntu CI runs.

---

## 1. Work Package Map & Disjoint File Ownership

### Wave 0 — Bootstrap & Phase-0 Profiler Baseline
| WP | Branch | Content | Dependencies | Gate |
|---|---|---|---|---|
| **260** | `wp/260-solver-accel-program` | Program ADR, this plan, `ladruno-solver-accel` skill, `.gitignore`/`.agents/` setup, `AGENTS.md` row, 3 quirk fragments | None (merges first) | Static gates green; Claude review + owner D1–D4 |
| **260-w0** | `wp/260w0-phase0-profiles` | Run `-perStep` + `-deep` profiler sweep on 5 north-star decks (`origin/ladruno` + `2CLD12`), record exact phase % ceilings on the `program` issue (no `SRC/` edits) | `260` + owner D1 | Baseline table posted on `program` issue |

### Wave 1 — Parallel Surgical Wins (1–2 Weeks, Mostly `G-BYTE`, Disjoint Files)
| WP | Branch | Scope & Disjoint File Ownership | Speedup Target | Gate |
|---|---|---|---|---|
| **260a** | `wp/260a-soe-scatter-hygiene` | `MumpsSOE::addA` binary search + skip `rowA++/--` on `JOB=3` + `std::sort` `setSize()` (incl. local-vertex `MumpsParallelSOE`); `PARDISOGenLinSOE::setSize` `std::sort`; `KrylovAccelerator.cpp:200-216` dead loop; fix `PARDISOGenLinSolver.cpp:367` comment & expose `iparm[23]/iparm[1]/iparm[24]=2/iparm[10/12]`; `AcceleratedNewton::sendSelf` reuse serialization | 3–20× `Mumps`/`setSize`, 1.3–2× PARDISO factor, faster `KrylovNewton` | `G-BYTE` (default) + `G-OPT` (flags) |
| **260b** | `wp/260b-explicit-wave-hygiene` | `ExplicitBathe` CDL mass-cache + `-commitSolveState` + hoisted `-lnvd`; ADR-67 P1.1/P1.2′ consistent-SMS & `MPIDiagonalSOE` alloc/halo hoist + rank-1 $O(n)$ centroidal matvec + `CriticalTimeStep` massless check; T7 inertia skips (`LadrunoElementShell.h`, `Brick20`, `CPML`); `LadrunoBrick::formUri`/BV $c_d$ cache; `ConstraintProjector` Cholesky $L$; global $O(N)$ $\alpha_M$; `LadrunoCPML` OpenMP opt-in | 1.4–2.2× explicit (`CDL`/`ExplicitBathe`), 4–8× consistent-SMS PCG | `G-BYTE` / `G-EQUIV` |
| **260c** | `wp/260c-logstrain-stressonly-implex` | `LogStrainNDMaterial::setTrialF()` 5-bug fix; Contract #1 `setTrialStrainStressOnly`/`setTrialFStressOnly`; immediate `doTangent && !mp.implex` (`LadrunoConcrete3DKernel.h:3253`); post-`#1153/#1154` commit-deferred `-implex` in `Concrete3D` & `HS` (`Ge_implex` freeze); `LadrunoUP` OpenMP opt-in | 3–6× finite-strain GP; 1.5–4× unused-tangent skips; 4–15× `-implex` | `G-BYTE` / `G-EQUIV` |
| **260d** | `wp/260d-implicit-bathe-adaptive` | `LadrunoImplicitBathe` ($\rho_\infty$-Bathe + `getEffectiveTangentKey` + `ScopedRecorderSuppression`) + `LadrunoAdaptiveTransient` (HAFTOL half-step LTE + quantized $2^k/\sqrt{2}^k$ rungs) + vectorized `addModalDampingForce` | 1.5–2× fewer steps + 2× fewer factors vs `TRBDF2`; 2–5× adaptive runs | `G-EQUIV` / `G-OPT` |
| **260e** | `wp/260e-allocator-lto-omp` | Unify Windows OpenMP on `libiomp5md` (remove `vcomp` dual-pool); evaluate/wire `mimalloc` + `/GL /LTCG` + `/arch:AVX2` (per D2/D4) | 10–30% global across all targets | `G-BYTE` + `dumpbin` check |

### Wave 2 — Core Kernels, Solvers & SSI Upgrades (3–8 Weeks)
- **`260f` (3D/MPM/MCFT Threading + Persistent Assembly):** De-static & opt-in `LadrunoBrick*`, `Bezier*`, `LadrunoSolidShell`, Mode-A Implicit MPM, Model C `LadrunoNDFiberSection3d`/`MCFTStrip`; `LadrunoFastScatter` + persistent `LadrunoThreadedAssembly` + colored scatter + sparse `TransformationFE` $T^T K T$ + `Domain` snapshot cache.
- **`260g` (Mixed-Precision GMRES-IR & Sym-Preconditioned Krylov):** FP32 PARDISO (`iparm[27]=1` `float*` shadow) + MKL RCI `dfgmres` FP64 GMRES-IR; `LadrunoMumpsMixed` GMRES-IR upgrade; `LadrunoSymPrecondKrylovSOE` ($K_{\text{sym}}$ factor + FGMRES).
- **`260h` (Analytic Return Maps & Outer Loop):** `Concrete3D` analytic `damagedTangent` & Voigt lift + Cardano/Smith 3×3 eigensolver; `LadrunoHS` analytic R4/R5 Jacobian & IFT Newton–Brent; `SANISAND` `SAS-RK23/45`; `ASDPlastic3D` closed-form VM/DP; Eisenstat–Walker Inexact Newton + `LadrunoAndersonNewton` + two-loop Powell-damped `LadrunoLBFGS`; optional modal Woodbury & static PTC.
- **`260i` (Global IMPL-EX + `LadrunoSuperelement` + `H5DRM` v2):** `LadrunoImplexAlgo`/`LadrunoImplexLinear`; `LadrunoSuperelement` wrapping `OPS_LadrunoCMS`; velocity-consistent & Biot $u\text{–}p$ `H5DRMLoadPattern` v2.
- **`260j` (`LadrunoContact` Acceleration):** Slave-aware corner-span bucket sort / BVH for 3D Mortar (`LadrunoContactHandler.cpp:2171,2345`), MPI `-reemit`, threaded `LadrunoContactFE`, ALM active-set audit.

### Wave 3 — Transformative Engines & HPC Frontier (2–4+ Months)
- **`260k` (`LadrunoGPStore` SoA + 7×7 Schur SANISAND CPPM):** 64-byte SoA + two-pass active-set compaction + Morton/RCM internal cache ordering.
- **`260l` (`LadrunoExplicitEngine` + Clustered LTS):** WP-232 Mode-F `LoadPattern` flat explicit engine + flat `CPML`/`Contact`/`EmbeddedNode`/`PorousOverlay` + power-of-two LTS subcycling + FP32-force/FP64-state.
- **`260m` (`LadrunoAMGSOE` + Monolithic Biot Schur + Cluster Direct):** Refreshed Boost-free `OTHER/AMGCL` (6 rigid-body near-nullspace modes + 3×3 block SA-AMG + monolithic `LadrunoUP` fixed-stress Schur); `LadrunoClusterPardisoSOE`; weighted METIS `vwgts`; optional `cuDSS`.
- **`260n` (`LadrunoIMEX` Spatial Partitioning):** Phase A Hughes–Liu pruned-LU IMEX → Phase B Prakash–Hjelmstad dual-Schur multi-rate subcycling atop `260l`.
- **`260o` (`LadrunoSpectralHex` + Pure-$M$ Integrators + GPU Explicit/MPM):** $p=3\text{–}5$ GLL SEM + pure-$M$ Hulbert–Chung / 2-stage Soares + `LADRUNO_HD` CUDA/HIP backend (per D3).

---

## 2. Landing & Review Cadence

1. Merge **`WP-260`** first after Claude review + owner D1–D4 answers.
2. Run **`260-w0`** baseline profiling, post numbers to the `program` issue, then launch **`260a`**, **`260b`**, **`260c`** (and `260d`/`260e`) in parallel worktrees under `.agents/worktrees/`.
3. Serialize local `build.bat` per worktree → run local pytest + static gates + `-perStep`/`-deep` benchmarks → push Draft PRs → run Zone-A Ubuntu CI in parallel → submit each PR for **intensive Claude review** → address findings → merge via owner merge commit.
