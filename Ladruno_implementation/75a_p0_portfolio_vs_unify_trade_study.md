---
title: ADR-75 P0 trade study — portfolio (PARDISO+MUMPS) vs unify-on-MKL
project: Ladruno
status: DECIDED — portfolio confirmed (decisive finding: MUMPS is load-bearing for 4 subsystems)
priority: medium
owner: nmora
amends: 75_ladruno_sparse_direct_strategy_adr
tags:
  - adr
  - performance
  - solver
  - pardiso
  - mumps
  - trade-study
  - p0
---

# ADR-75a — P0 trade study: portfolio vs unify-on-MKL

> The P0 decision [[75_ladruno_sparse_direct_strategy_adr]] §9 elevated from a footnote after the
> adversarial review weakened the original "MUMPS is already proven" argument. **Resolved by a
> code-verified fact that neither the ADR nor the review had checked.**

## The question

ADR-75 chose a **portfolio**: MKL PARDISO on desktop (shared-memory), MUMPS on cluster (MPI). The
review objected that this was under-argued — once PARDISO turned out to be a *2019 prototype needing
real work* (§12), the "don't replace a proven path" defence weakened, and **unify-on-MKL** looked
attractive:

- **Option A — Portfolio:** PARDISO (desktop) + MUMPS (cluster). Two solver families.
- **Option B — Unify on MKL:** PARDISO (desktop) + `cluster_sparse_solver` / MKL distributed PARDISO
  (cluster). One family, near-identical `iparm` interface, one dependency, one test surface.

Option B's whole case rests on **"one dependency, one test surface."**

## Decisive finding — that premise is false

**MUMPS cannot be removed from this fork.** It is load-bearing for four subsystems beyond
`system Mumps`, verified by `dmumps_c`/`DMUMPS_STRUC_C` call sites:

| Subsystem | Files | Coupling |
|---|---|---|
| **CMS eigensolver** (ADR-1000) | `ladrunoCMS/LadrunoCMSMumps.{cpp,h}` | **`FATAL_ERROR "LADRUNO_CMS=ON requires MUMPS_DIR"`** (`CMakeLists.txt:645-646`) — a *hard build gate* |
| **Distributed FEAST / modal** (ADR-43/45) | `eigenSOE/LadrunoDistBlockZKernel.{cpp,h}` | distributed inner solve; the `-commSplit` sub-communicator work exists *for this* |
| **PFEM** (fluid) | `sparseGEN/PFEMSolver_Mumps`, `PFEMCompressibleSolver_Mumps`, `PFEMUnifiedSolver_Hybrid` | 3 solver back-ends |
| **`system Mumps`** | `linearSOE/mumps/*` | the linear SOE (the only part this ADR is about) |

So under **Option B the build still links MUMPS** — for CMS, FEAST, and PFEM. You would **not** shed
the dependency, **not** shrink the test surface, and **not** retire the MUMPS build recipe
(`scivision/mumps`, LP64, the `intsize64=OFF` gotcha). You would simply **add a third solver family**
(`cluster_sparse_solver`) alongside the two you already maintain.

Option B's premise inverts: it *increases* dependencies rather than unifying them.

## Evaluation

| Criterion | A — Portfolio | B — Unify on MKL |
|---|---|---|
| Dependencies | MKL + MUMPS (MUMPS **mandatory anyway**) | MKL + MUMPS (**still**) + `cluster_sparse_solver` path |
| New unproven code | PARDISO desktop only | PARDISO desktop **+ a new distributed solver** |
| Cluster path maturity | **ADR-74-hardened, proven at 18.6 M hex / np240** | unproven here; would need its own scaling campaign |
| Distributed-assembly risk | none new (`MumpsParallelSOE` exists, `ICNTL18=3`) | **new distributed SOE** — the highest-risk seam (§11.5) |
| Interface consistency | two idioms (`iparm` vs `ICNTL`) | one idiom — **B's only real win** |
| Reuses ADR-43 `-commSplit` / FEAST wiring | yes | no — would need re-doing |
| Effort | PARDISO (medium) | PARDISO (medium) **+ distributed port (large)** |

**B's single genuine advantage** is one `iparm` idiom instead of two. That is a developer-ergonomics
gain, purchased with a large distributed-solver port and the abandonment of a hardened, measured
cluster path — while *still* shipping MUMPS.

## Decision — **Option A (portfolio), confirmed**

Not "leaning" any more: the objection that reopened this question does not survive the finding that
MUMPS is mandatory regardless. Concretely:

1. **MUMPS stays the cluster solver.** It is already in the build (CMS/FEAST/PFEM force it), already
   hardened, already measured at scale. Using it for `system Mumps` too is *free* by comparison.
2. **PARDISO is added for desktop only** — its scope stays exactly ADR-75 §3 Lane 1 (serial build,
   shared-memory). No distributed PARDISO.
3. **`cluster_sparse_solver` is formally OUT** — moved from "open question" to **anti-goal**, unless
   a future measured deck shows MUMPS-distributed is the binding constraint *and* MKL demonstrably
   beats it. (Even then, MUMPS stays linked for CMS/FEAST/PFEM.)

### What would reopen it
Only a *measured* result: MUMPS-distributed factorization becoming the dominant cluster cost on a
real deck, **and** a spike showing `cluster_sparse_solver` materially better on the same model. Note
even that would not remove MUMPS from the build — it would only change which solver `system` uses.

## Consequences for ADR-75

- §9 P0 open question → **closed** (this document); `cluster_sparse_solver` → §8 anti-goals.
- P1 (PARDISO desktop) is **unblocked** and its scope is *narrowed and fixed*: serial/shared-memory
  only. The medium-effort estimate from §12 stands unchanged.
- P2 (MUMPS cluster tuning) is confirmed as tuning of a solver the fork must ship anyway.
- **Bonus insight for §11.6 (the non-MKL desktop gap):** since MUMPS is mandatory in *any* build with
  CMS/FEAST/PFEM enabled, a non-MKL desktop could in principle fall back to a sequential MUMPS after
  all. This does **not** revive serial-MUMPS as P1 work (PARDISO remains the desktop tool where MKL
  exists), but it is a cheaper fallback than previously assumed, and worth noting if a non-MKL
  desktop lane ever becomes a real requirement.

## Provenance
Call sites verified by `grep -rln "dmumps_c\|DMUMPS_STRUC_C" SRC/`; the CMS hard gate read at
`CMakeLists.txt:643-650`. No source modified by this study.
