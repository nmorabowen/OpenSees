---
title: "WP-195 — stdBrick tangent cache; zeroLength springs and dashpots on the WP-192 thread opt-in"
project: Ladruno
type: performance work package (element lane)
status: "implemented, draft PR — not merged; cache default ON, threads default 1"
priority: medium
owner: nmora
related:
  - "[[192_ladruno_threaded_element_loops_adr]] (WP-192: the 3C profile that found the brick tangent; the per-type opt-in zeroLength joins)"
  - "[[40_ladruno_performance_adr]] (lever 5: remove redundant recompute; anti-goal: dirty-flag memoization in a material)"
  - "paper repo: review/engstruct-2026/30_campaigns/timehistory/perf/PROFILE_4D.md"
tags: [performance, element, brick, cache, zerolength, openmp, wp-195]
updated: 2026-10-14
---

# WP-195 — stdBrick tangent cache; zeroLength springs and dashpots threaded

> [!summary] The short version
> **Item 1.** WP-192's 3C profile found `stdBrick` at **54 % of element time**: it re-forms its
> 24×24 stiffness (8 Gauss points × 8 × 8 Bᵀ D B products, 63 µs per element per call) at
> every tangent call, though the Paper 1 soil is elastic (rev3 3C and 4D: every brick is
> `nDMaterial ElasticIsotropic 21 3400.0 0.262 2.4e-09`, interior and DRM layer alike).
> The element now keeps its last stiffness together with the **raw bits** of the eight
> material tangents and of the nodal coordinates that built it, and returns it whenever, at
> the next tangent call, those bits are unchanged. Elasticity is *observed* through the
> material, never assumed; the result is **bit-identical to the uncached path by
> construction** and verified byte for byte on the rev3 3C window.
>
> **Item 2.** The Tier 2 springs and dashpots (`zeroLength` + `Parallel(Elastic, Viscous)`,
> 2 487 in 2C) were un-audited in WP-192 and ran serially in place. zeroLength's
> state path wrote **class-wide** matrices shared by every instance; they are per thread now,
> and zeroLength / Parallel / Viscous join the per-type opt-in. The rev3 2C window at
> 4 threads equals 1 thread byte for byte.
>
> Branch `wp/195-brick-tangent-cache`, **based on `wp/192-threaded-element-loops`** (item 2
> needs WP-192's opt-in). Item 1 is one self-contained commit that cherry-picks onto
> `ladruno` alone.

## 1. Item 1 — the cache rule

`Brick::getTangentStiff()` (no Damping object, cache on, fewer than 16 consecutive misses):

1. **Residual** by `formResidAndTangent(0)` — the uncached code itself (same `getStress()`
   per Gauss point, same arithmetic). It also sets `xl` (`computeBasis`).
2. **Material tangents.** `getTangent()` of each of the eight points; the 36 raw doubles
   are copied out.
3. **Hit test.** `memcmp` of `xl` (24 doubles) and of each point's 36 tangent doubles
   against the bits stored with the cached stiffness. All equal ⇒ copy the cached 24×24
   into `stiff`.
4. **Miss.** Rebuild `stiff` from the fetched tangents with the **same loop nest, operations
   and order** as `formResidAndTangent(1)` (point *i*, node *j*, node *k*; same
   `addMatrixProduct` calls), then store `stiff`, `xl` and the tangent bits (36 doubles when
   all eight points carry the same tangent, else 288).

Why bit-identical: on a hit the cached matrix was computed by that same code from inputs whose
bits are equal to today's; on a miss it *is* that code. The comparison is on **bits**, not
`==`: `+0.0 == -0.0`, and both propagate signed zeros into the products.

The one ordering difference from `formResidAndTangent(1)`: there, point *i*'s `getTangent()`
follows point *i*'s `getStress()` immediately; here all eight `getStress()` calls come first.
Each point is its own material instance, so this matters only for a material whose
`getTangent()` depends on *another* instance's `getStress()` — shared mutable state, which is
exactly what WP-192's thread-safety contract already excludes.

### Why the bits and not an "is elastic" flag

No nD material in the fork exposes "is elastic" (`NDMaterial`, `Material` checked), and a
flag would not be enough even if one did:

* **A parameter update never reaches the element.** `parameter 1 element 5 E` forwards to
  the eight materials, and each *material* registers itself with the `Parameter`;
  `updateParameter 1 v` then calls the material directly. `Brick::updateParameter` is not on
  that path. A cache invalidated only by element hooks returns the old stiffness —
  converged, and wrong. (Quirk filed: `ledger/quirks/WP-195-element-tangent-cache-keyed-on-raw-bits.md`.)
* `setNodeCoord` changes the geometry without telling the element; a staged material
  changes its tangent at `updateMaterialStage`; a plastic point leaves the elastic range.

The bit compare needs no cooperation from any of them. The element's own hooks still
invalidate (belt and braces): `revertToStart` (also resets the give-up count), `setDomain`,
`setDamping`, `setParameter`, `updateParameter`, `recvSelf`.

### Scope and cost

| | |
|---|---|
| Elements | `stdBrick` (`Brick`, classTag 56). Not BbarBrick / SSPbrick / LadrunoBrick / Brick20 |
| Off | a Damping object on the element (its stiffness multiplier is not a material tangent); `LADRUNO_BRICK_TANGENT_CACHE=0` (read with `getenv` at element construction); a non-6×6 material tangent (disables for good, plain path) |
| Give-up | 16 consecutive misses of a valid cache (a tangent that keeps changing): the buffer is freed and the plain path runs until `revertToStart` |
| Memory | 636 doubles (4.9 KB) per element with one tangent for all points, 888 (6.9 KB) otherwise, allocated at the first tangent call. 3C: 27 000 bricks ≈ 137 MB over the model |
| Hit cost | 8 `getTangent()` (elastic: fills 12 entries), 312 doubles compared, 576 copied |
| Query | element response `ladrunoTangentCache` → `[enabled, state, hits, builds, consecutive misses]` |
| Threads | per-element storage, written only by the element's own call; the scratch is `thread_local` as WP-192 left it, so the cache is safe in the WP-192 loops |

### Same pattern elsewhere (checked, not changed)

* **ASDAbsorbingBoundary3D**: `getTangentStiff()` and `getDamp()` rebuild K and C every call,
  but both are pure functions of the stage, `G/v/rho`, the boundary type and the nodal
  coordinates — all element-owned and changed only through `updateParameter`/`recvSelf`/
  `setDomain`. A cache would be exact with an element-level dirty flag (no material
  involved). Its share in WP-192's 3C profile: 101 ms/step of 2.6 s (4 %); memory would be
  up to 2 × 72² doubles per element. A follow-up if the element loops become the floor.
* **ASDEmbeddedNodeElement**: forms K every tangent *and* every residual (`F = K·U`); K is a
  function of the geometry only. 22 ms/step (0.9 %). Not worth the memory today.

## 2. Item 2 — zeroLength springs and dashpots on the opt-in

Tier 2 decks (rev3 `Tier_2/2C/base/elements.tcl`, `materials.tcl`): 3 512 ASDShellQ4, 192
forceBeamColumn (both opted in by WP-192) and **2 487 `zeroLength`**:

```
element zeroLength <tag> <ground> <node> -mat P1 P2 P3 -dir 1 2 3 -doRayleigh 0 -orient 1 0 0 0 1 0
uniaxialMaterial Parallel P 9001 9002 -factors k c      ;# 7 461 of them
uniaxialMaterial Elastic  9001 1.0
uniaxialMaterial Viscous  9002 1.0 1.0
```

Audit over update / tangent / damp / mass / residual / commit (WP-107 `Element.h` items 1–5):

| Class | Shared state on the path | Fix | Opt-in |
|---|---|---|---|
| ZeroLength | `theMatrix`/`theVector` of **every** instance point at the class-wide `ZeroLengthM2/4/6/12`, `V2/4/6/12`: tangent, damp, mass and residual all write them | pools `static thread_local`; the state path reads `ladrunoM()`/`ladrunoV()` = this thread's copy, selected by `numDOF` exactly as `setDomain` pairs them; `displaySelf` statics `thread_local` | update + form + commit; refused with a Damping object; every uniaxial (incl. the `-doRayleigh 2` damping copies) must opt in |
| ParallelMaterial | none (trial strain/rate, children, factors per instance) | — | = its children's answer |
| ViscousMaterial | none (statics only in send/recvSelf) | — | true |
| ElasticMaterial | — | — | true (WP-192) |

`-doRayleigh 1` uses `Element::getDamp()`/Rayleigh forces, whose pool WP-192 made per thread;
`t1d`, `fd`, `d0/v0` and the ADR-96 passenger buffers are per element; `update()` builds its
difference vectors locally; nothing on these paths writes `opserr` or throws.

The WP-107/WP-192 tests used a `zeroLength` as their "un-audited element" probe; they now
use a `Truss` (still un-audited).

## 3. Verification

### 3.1 Item 1 — cached = uncached, byte for byte

| Configuration | Comparison | Result |
|---|---|---|
| rev3 3C serial emission (`apesees_3C.tcl`, Pardiso, `MKL_NUM_THREADS=1`; 47 524 FEs incl. 27 000 stdBrick), 4 gravity stages × 3 increments + 5 transient steps, every node dof 1–3 every step at 17 digits | WP-192 base binary vs WP-195 cache off vs WP-195 cache on | **identical** (one md5 for all three, 17 records) |
| `tests/test_wp195_brick_tangent_cache.py`: distorted 27-brick block; full-field md5 (disp + reaction of every node + element stresses) every step; each setting in a fresh child (the switch is read by `getenv`) | cache on vs off for: elastic static + Newmark/Newton; yielding J2 (rebuilds, then gives up); `updateParameter` (material route) and `setParameter -ele` (element route) mid-run; `setNodeCoord` mid-run; `reset` + replay | **identical** in all 6; counters show the cache hit (elastic: 1 build, ≥ 10 hits per element) and rebuilt where it must |
| Mutation M1: the hit test ignores the material-tangent bits | — | **3 tests red** (J2, `updateParameter`, reset); `setParameter -ele` stays green because the element hook invalidates — which is why the bit compare, not the hooks, carries correctness |
| The fork's 73 other test files that build a stdBrick | cache on | 849 passed; the 3 failures + 1 error are environmental or pre-existing: `test_ladruno_sanisand_cppm_newton` (2) fails identically with the cache off **and** on the WP-192 base binary; `test_adr94_matrix` (relative path) and `test_ladruno_sanisand_sasme` (child without pytest) are runner-environment issues |

### 3.2 Item 2 — 2C threaded = serial, byte for byte

| Configuration | Comparison | Result |
|---|---|---|
| rev3 2C serial emission (`tools.campaign_t2 emit 2C --rev3`, Pardiso, `MKL_NUM_THREADS=1`), 4 gravity stages × 3 increments + 10 transient steps (motion shifted to its strong part) | `LADRUNO_THREADS=4` vs 1 thread | **identical**: every node dof 1–3 every step (17 digits), `energy_balance` (EnergyBalance -v2, 10 regions), `th_steps`, and all 266 numeric datasets of the five `recorder ladruno` files; the loops announce `6191 of 6191 FE_Elements` threaded, update loop `THREADED` |
| `tests/test_wp195_zerolength_threads.py`: ASDShellQ4 slab on 35 zeroLength Parallel(Elastic, Viscous) springs/dashpots, static + KrylovNewton/Newmark shake, full-field md5 (disp, vel, reaction, a spring force) every step | 1 vs 2 vs 4 threads | **identical**; no serial-in-place element |
| Mutation M2: the ZeroLength pools back to plain `static` | — | the zeroLength test splits into **three** digests (1/2/4 threads all differ) |
| WP-192 + WP-107 suites | — | 24/24 |

## 4. Measured — 3C window, 2 ranks (as WP-192's profile)

WP-192's own 3C deep-profile window (`cases.case_t3 --cell 3C --rev3 --ranks 2`, MUMPS,
Transformation, KrylovNewton, Newmark), shortened to 3 gravity increments per stage and
**30 transient steps**; `profiler start -deep -perStep -memory` around the transient; local
i7-13700KF, ranks pinned to P-cores 2 and 4, 1 thread per rank; a shared box (another
agent's 4-rank OpenSeesMP job ran throughout). ms/step, rank 1 (critical):

| ms/step, rank 1 | cache off | cache on | change |
|---|---|---|---|
| stdBrick, all element loops | 1 065 | 316 | −749 (3.4× faster) |
| stdBrick tangent (16 600 bricks, 1 call/step) | 847 (51 µs/element) | 132 (8 µs/element) | 6.4× faster |
| element loops, all types (update + tangent + residual) | 1 945 | 1 095 | **−850 (−44 %)** |
| linearSolve (MUMPS) | 18 440 | 14 181 | not this WP: shared-box noise |
| **step** | **20 629** | **15 491** | see note |
| transient wall, 30 steps | 619 s | 465 s | |

Rank 0 is within 1 % of rank 1 on every row. The element saving is exact work removed
(0.85 s/step, ~4 % of the cache-off step); the rest of the measured step difference is
MUMPS time, which this WP does not touch and which moved by 4.3 s between the two runs
because another agent's MPI jobs shared the box. Same iteration count (465) in both runs;
the 2-rank answers agree to ~1e-9 relative (MUMPS under OpenSeesMP is not run-to-run
reproducible — WP-188 quirk; byte identity is shown on the serial window, §3.1).

In WP-192's profile of the same window stdBrick was 54 % of element time; with the cache it
is 29 % (316 of 1 095 ms), the residual (8 calls' worth of `getStress` + Bᵀσ) now its larger
part. On today's 3C the step is still MUMPS-bound (≥ 90 %); the saving matters once
factorization reuse / METIS / fewer ranks shrink the solve (WP-188/192 rationale).

## 5. Not done here

| Item | Note |
|---|---|
| ASDAbsorbingBoundary3D K/C cache | exact with an element dirty flag (no material); 4 % of 3C element time; memory up to 2 × 72² doubles per element |
| ASDEmbeddedNodeElement K cache | 0.9 %; K is formed twice per iteration (tangent and residual) |
| Other bricks (BbarBrick, SSPbrick, LadrunoBrick, Brick20) | same rule applies; not in the Paper 1 decks |
| A cache-hit path that skips `getTangent()` entirely | needs a material-side "tangent unchanged since" stamp; not worth the contract while an elastic `getTangent()` is ~12 stores |
| zeroLength with a Damping object | refused, as for Brick/ASDShellQ4 (Damping not audited) |

## 6. Caveats

* The cache is **on by default**. It changes no result (bit-identical by construction and by
  test), costs ≈ 5 KB per brick, and can be switched off per process with
  `LADRUNO_BRICK_TANGENT_CACHE=0` (read when the element is built; on Windows, Python's
  `os.environ` does not reach the C runtime's environment — set it before starting Python).
* The give-up after 16 consecutive misses is per element and permanent until
  `revertToStart`: a plastic brick that later unloads elastically does not re-enable it.
* Thread safety of item 2 is evidenced by byte identity at 1/2/4 threads (tests) and 1 vs 4
  (2C window) plus the M2 mutation; no ThreadSanitizer run (WP-192 §6.2 follow-up).
