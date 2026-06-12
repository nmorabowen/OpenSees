---
title: LadrunoParallelNumberer — distributed DOF numbering (kill the P0 gather/merge serial term)
project: Ladruno
status: draft
priority: medium
owner: nmora
tags:
  - implementation
  - solver
  - mpi
  - adr
---

# LadrunoParallelNumberer — distributed DOF numbering

**Study-first ADR.** The motivating measurements are done (see *Why*);
the design below is a sketch with open questions that must be resolved
against the source before implementation. Do not start coding from
this document alone.

## What

A distributed DOF numberer for OpenSeesMP that eliminates the
rank-0-serialized gather → merge → number → scatter inside
`ParallelNumberer::numberDOF`. Each rank numbers the DOFs it primarily
owns from an `MPI_Exscan` prefix offset and resolves shared-node DOFs
by a neighbor (halo) exchange — O(N/np) + O(log np) instead of
O(N_total) serialized through P0.

Scope v1: explicit/diagonal analyses (`MPIDiagonal`,
`CentralDifferenceLadruno`, `ExplicitBathe`), where equation-ordering
quality is irrelevant. Out of scope (follow-up study): implicit +
`Mumps` (MUMPS does its own internal fill-reducing ordering, so the
numberer's ordering quality is *probably* also irrelevant there — but
that claim must be verified, not assumed). A secondary target,
separable: `MPIDiagonalSOE::setSize`'s all-pairs shared-DOF discovery.

## Why

Measured on Esmeralda (2026-06-12, jobs 143859–143869; 66,564-hex /
227k-DOF apeGmsh plane-wave model, in-deck `clock milliseconds`
segment timers; artifacts
`Downloads/results/ladruno_wave_propagation/bench_adr0061/` on the
Windows machine, `bench2.py` is the rerunnable harness):

| run | np | model parse+build | domainChange | **analyze step 1 (setup)** | steady step |
|---|---|---|---|---|---|
| ParallelRCM    | 16 | 130 ms | ~0 | **2655 ms** | 21 ms |
| ParallelPlain  | 16 | 130 ms | ~0 | **3927 ms** | 75 ms |
| ParallelRCM ×2 | 64 | ~57 ms | ~0 | **2373 / 2502 ms** | ~20 ms |
| ParallelPlain ×2 | 64 | ~54 ms | ~0 | **3792 / 3805 ms** | ~17 ms |

Three facts fall out:

1. **The per-run fixed cost is first-analyze analysis setup**, and its
   dominant piece is DOF numbering. Deck parse, `domainChange` (stamp
   bump), and chain construction are all ~free.
2. **It does not shrink with ranks** (2655 → ~2440 ms going 16 → 64):
   the merged global graph P0 processes is the same 227k DOF
   regardless of np. `SRC/analysis/numberer/ParallelNumberer.cpp`
   (`numberDOF`, comment at line ~111) is explicit about the
   mechanism: P0 collects every partition's DOF graph over the
   channels, merges them into one graph, numbers it (through the
   provided GraphNumberer — RCM for `numberer ParallelRCM`), and sends
   each rank back its dofTag → startEqn mapping. O(N_total) work and
   communication serialized through one rank — a true Amdahl term.
3. **Anti-finding: `ParallelPlain` is ~1.4–1.5 s SLOWER than
   `ParallelRCM`** (reproducible, ±3% across reps). The
   no-GraphNumberer fallback path ("number by subdomain order") is
   less efficient than handing the merged graph to RCM. So the cheap
   mitigation "drop RCM for explicit runs" is wrong; there is no
   stock-numberer escape — only a structural fix helps.

Projection: setup is roughly linear in total DOF. At production
fidelity targeted by the apeGmsh pipeline (P=8/f=8 ≈ 5.6M hexes ≈ 18M
DOF, ~80× the bench), that is **~3 minutes of serialized numbering per
(domainChange + analyze)** — i.e., **per stage**, since apeGmsh emits
an unconditional per-stage `domainChange` (required for recorder
MODEL_STAGE splitting, apeGmsh PR #633). A 10-stage staged-SSI model
would burn ~half an hour in numbering alone before any physics. This
is the successor serial term now that per-rank deck emission (apeGmsh
ADR 0061) has cleared deck text RAM/NFS off the critical path.

For today's 66k–1M-hex single-stage runs the 2.5 s amortizes and
nothing needs to change — which is why this is *study-then-implement,
demand-gated on multi-M-DOF staged production*, not urgent.

## Where

- Reference (current impl, the thing being replaced):
  `SRC/analysis/numberer/ParallelNumberer.{h,cpp}` — note
  `numberDOF()` (P0 merge path and the slower no-GraphNumberer
  fallback), `mergeSubGraph()`.
- New code: `SRC/analysis/numberer/LadrunoParallelNumberer.{h,cpp}`
  (subclass of `DOF_Numberer`, sibling of `ParallelNumberer`).
- Registration: interpreter `numberer` dispatch in
  `SRC/interpreter/OpenSeesCommands.cpp` (~line 1572, where
  `ParallelPlain` / `ParallelRCM` dispatch) **and classic Tcl** — heed
  the fork PR #233 lesson: classic-Tcl registration goes through the
  flat factory-table chain slot, not a new else-if branch (MSVC C1061
  nesting limit).
- Secondary target (separate slice): `SRC/system_of_eqn/linearSOE/
  diagonal/MPIDiagonalSOE.cpp` `setSize()` — `maxNeighbors =
  numProcesses` all-pairs sorted-ID exchange for shared-DOF discovery.
  Fine at np ≤ 256; revisit only if numbering stops dominating.

## How (sketch — to be validated in the study phase)

The contract the replacement MUST reproduce (this is what makes
OpenSeesMP work at all): **a DOF shared by several ranks gets the SAME
equation number on every owning rank.** Today that falls out of the
P0 merge keying vertices by node tag. The distributed version has to
re-derive it:

1. **Ownership rule.** Primary owner of a shared node = lowest rank
   that owns it (the same convention apeGmsh's partitioned emit
   already uses for additive mass/load dedup). Each rank classifies
   its DOF_Groups as owned vs ghost. Open question Q1 below: where
   does a rank *learn* which of its nodes are shared and with whom —
   today the numbering itself is that mechanism, so the new numberer
   needs its own discovery exchange (the sorted-ID handshake in
   `MPIDiagonalSOE::setSize` is the pattern to copy, or a
   tag-range/hash-bucketed exchange to avoid all-pairs).
2. **Count + prefix.** Each rank counts equations of its *owned*
   DOF_Groups; `MPI_Exscan(MPI_SUM)` gives the rank's starting
   equation number.
3. **Number owned.** Walk owned DOF_Groups in local (graph) order,
   assigning sequential equation numbers. No ordering quality is
   attempted — v1 is explicitly for diagonal systems.
4. **Halo exchange.** For each shared node, the primary owner sends
   (nodeTag, startEqn) to the other owners; they stamp their local
   DOF_Groups. This is the same information P0's scatter carries
   today, restricted to actual neighbors.
5. **Constraint DOFs.** `-2`/`-3` ID codes (constrained/condensed,
   handled by the Transformation handler before numbering) follow the
   same owned/ghost classification — study item Q2.

```tcl
# proposed surface — drop-in alternative
numberer LadrunoParallel
```

Testing / acceptance gates:

- **Identity gate:** same model, `ParallelRCM` vs `LadrunoParallel`,
  np {2, 8, 16}: results identical to round-off (the apeGmsh
  plane-wave gates run at ≤1.5e-15 of peak; reuse those decks).
- **Shared-DOF gate:** a deliberately interface-heavy partition
  (apeGmsh `tests/opensees` partitioned fixtures) — every shared node
  numbered identically on all owners (assert via a debug dump).
- **Scaling gate:** bench2-style segment timers; setup must drop
  ~np× at fixed N and stay flat under np at fixed N/np (weak scaling).
- **Staged gate:** multi-stage deck (unconditional per-stage
  domainChange) renumbers correctly every stage.

## Risks / open questions

> [!question]
> Q1 — Shared-node discovery: what is the cheapest correct way for a
> rank to learn which of its node tags exist on other ranks (and
> who)? Options: all-pairs sorted-ID exchange (MPIDiagonalSOE
> pattern, O(np²) messages), hash-bucketed owner lookup (O(np) with
> two alltoallv), or having the *deck* declare it (apeGmsh already
> knows the shared-node sets at emit time — a `# shared-nodes` table
> per rank would make discovery free, at the cost of coupling the
> numberer to deck metadata; probably wrong layer, but cheap to
> prototype).

> [!question]
> Q2 — Interaction with constraint handlers: Transformation condenses
> slave DOFs before numbering. Are the -2/-3 codes' semantics
> identical across ranks for a shared constrained node, and does the
> halo stamp need to carry per-DOF codes rather than just startEqn?

> [!question]
> Q3 — Implicit/Mumps follow-up: verify (don't assume) that MUMPS's
> internal ordering makes numberer ordering irrelevant for fill-in,
> in which case LadrunoParallel covers the implicit path too and
> ParallelRCM can be retired from the apeGmsh auto-emit default.

> [!question]
> Q4 — Where exactly do the remaining setup milliseconds go? The
> bench brackets numbering + handler + SOE setSize + integrator
> domainChanged together (~2.4 s RCM at 227k DOF). Before designing,
> instrument `ParallelNumberer::numberDOF` entry/exit (opserr +
> MPI_Wtime, the debug pattern already half-present in
> MPIDiagonalSOE::setSize) to confirm numbering is ≥80% of it — if
> handler or setSize is a co-dominant cost, the design should widen.

- No new dependencies; pure MPI + existing classes.
- Backwards compatibility: additive — `ParallelPlain`/`ParallelRCM`
  untouched; apeGmsh would gate emission of the new numberer on a
  fork-feature probe (tags-≥33000-style point-of-use gating).

## Implementation log

*(empty — study phase not started)*
