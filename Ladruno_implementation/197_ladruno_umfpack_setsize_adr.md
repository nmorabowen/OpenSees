---
title: "WP-197 — UmfpackGenLinSOE::setSize rebuilds the sparsity pattern from scratch (#981)"
project: Ladruno
type: bug fix (vanilla file, upstreamable)
status: "implemented, PR #TBD — not merged"
priority: high
owner: nmora
related:
  - "fork issue #981 (reported by Patricio Palacios' agent, WP-138 strip footing, apeChat #ladruno msgs 32/38/43)"
  - "ADR-40 rank 8/10 (UmfPack Numeric persistence, the `factored` flag)"
  - "ADR-75 P1g (addA binary search, the DISCARDED miss detector that made this visible)"
tags: [solver, umfpack, sparsity, upstream-bug, wp-197]
updated: 2026-10-15
---

# WP-197 — UmfPack re-analysis after a DOF-graph change

> [!summary] The short version
> `UmfpackGenLinSOE::setSize(Graph&)` called `Ap.reserve`, `Ai.reserve`, `Ax.resize` and then
> `push_back`ed the new CSC pattern **onto the old one**. A second `setSize` after the DOF graph
> changed left `Ap[0..n]` describing the old pattern. The fix clears `Ap`/`Ai` and sizes `Ax`
> with `assign(nnz, 0.0)`; the solver also drops `Symbolic` on the empty-system early return.
> An unchanged-graph re-analysis and the `wipeAnalysis` workaround are byte-identical to the
> base build. Upstream master has the same code.

## 1. The defect

`setSize` runs at the start of an analysis and again whenever the domain changes (`sp`/`mp`
added or removed, element added or removed, `domainChange`). The second call appended:

- **pattern grows, `n` unchanged** (element added): `addA` binary-searches the old column
  ranges, finds no slot for the new couplings, drops them and prints the ADR-75 P1g warning
  "an element entry has no slot in the sparsity pattern; its contribution was DISCARDED".
  Symbolic factorization accepts the old prefix, so the solve uses a wrong tangent. Under
  `NormUnbalance` Newton stalls; under a displacement or energy test it can look converged.
- **`n` shrinks** (`sp` on a free DOF, `constraints Transformation`, the reported case): the
  old prefix holds row indices `>= n`, UMFPACK's symbolic analysis rejects the matrix,
  `Symbolic` stays 0 and every later solve prints "setSize has not been called".
- **nnz shrinks, `n` unchanged** (element removed): `Ax` is resized down while `Ap` still
  spans the old, larger pattern, so `addA` writes past the end of `Ax`.

## 2. Repro (tests/test_wp197_umfpack_setsize.py)

A 2-storey, 1-bay elastic frame, classic `OpenSees.exe`, Newton + `NormUnbalance 1e-9`,
two `LoadControl 0.5` steps, a graph change, two more steps. Cases: `sp_add`, `elem_add`,
`elem_remove`; each compared with `system FullGeneral` on the same deck and with UmfPack after
`wipeAnalysis` (a fresh SOE), relative tolerance 1e-9. `none` forces `setSize` on an
unchanged graph with `domainChange` and must print the same bytes as the run without it;
with `LADRUNO_WP197_BASE_DIST` set, the unchanged-graph runs are byte-compared with a base
build too.

Measured on the local Windows build (classic `OpenSees.exe`, base `4387c7d14` vs the fix):

| Test | base (before) | WP-197 |
|---|---|---|
| `sp_add` vs FullGeneral and vs wipeAnalysis | FAIL: `symbolic analysis returns -8`, "solver failed setSize()", second `analyze` returns -1 | pass (1e-9) |
| `elem_add` | FAIL: "DISCARDED", NormUnbalance fails to converge, `analyze` -3 | pass |
| `elem_remove` (nnz shrinks) | FAIL: Newton fails to converge, `analyze` -3 (no crash seen, but `addA` indexed past `Ax`) | pass |
| `sp_add` imposes `u = 0.002 t` exactly | FAIL | pass |
| unchanged graph, forced re-setSize prints the same bytes as without | FAIL (last-bits difference, see §3) | pass |
| byte-identical to base: unchanged graph, and `wipeAnalysis` for the three cases | (reference) | pass (md5-equal text) |

Module: 9 passed in 0.9 s with `LADRUNO_WP197_BASE_DIST` set (5 run in CI, 4 skip without
a base dist). On the base build the five non-base tests fail.

## 3. The fix

`UmfpackGenLinSOE::setSize`: `Ap.clear(); Ai.clear();` before the `reserve` calls, and
`Ax.assign(nnz, 0.0)` instead of `Ax.resize(nnz, 0.0)`. The old values were not only
meaningless for a new pattern, they were read: `umfpack_di_symbolic` takes `Ax`, and after
`resize` it held the last assembled tangent, while the first `setSize` sees zeros. So on the
base build a forced re-`setSize` on an UNCHANGED graph (`domainChange`) gave an answer that
differed in the last bits from the same run without it. With `assign` both print the same
bytes. This is the one place the fix changes bits on a graph that did not change, and it
changes them towards the run that never re-sized. A re-analysis without a `setSize` (the
common case) does not go through this code and is byte-identical to the base.

`UmfpackGenLinSolver::setSize`: the non-empty path already freed the persisted `Numeric`
(ADR-40) and the old `Symbolic` and ran a new symbolic analysis on every call, so a size change
always forces one; it failed only because it was fed the corrupt pattern. The empty-system
early return (`n == 0 || nnz == 0`) now frees `Symbolic` too, so no symbolic object outlives
its pattern. Nothing else changes: `missWarned` stays once per SOE, the ascending-column check
stays.

## 4. Sibling SOEs

Grep of `setSize` under `SRC/system_of_eqn/linearSOE/` for `reserve`/`push_back`: only
`UmfpackGenLinSOE` builds its pattern by appending without clearing. `PFEMGeneralLinSOE` also
appends but clears `rowInd`/`colPtr`/`A` first; `SparsePython*` clear first. The others
(`SparseGenCol/Row`, `SymSparse`, `Mumps`, `PARDISOGen/Sym`, `Itpack`, `Petsc`, band, profile,
full, diagonal) allocate arrays and fill them by index on every call. Not this defect, but
noted while reading: `SymSparseLinSOE::setSize` allocates `colA` with `new` on every call
without freeing the previous one (a leak per re-setSize, not a wrong answer).

## 5. Upstream

OpenSees/OpenSees master `UmfpackGenLinSOE::setSize` (lines ~95-134) is the same code, so
vanilla has the bug. The two-line SOE change is upstreamable as is (ledger
`vanilla/WP-197-umfpack-setsize-upstreamable.md`).
