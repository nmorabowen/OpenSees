---
title: "WP-196 — the commit phase's cost: EnergyBalance sweep threaded and de-duplicated, recorders measured"
project: Ladruno
type: performance work package (commit phase)
status: "implemented, draft PR #984 — not merged; numbers unchanged (byte-identical)"
priority: medium
owner: nmora
related:
  - "[[192_ladruno_threaded_element_loops_adr]] (WP-192: the loops and the measurement this WP starts from; based on its branch)"
  - "issue #959 (energy balance sweep ~10 µs per element-step)"
  - "WP-184 (the damping-outside-F check), ADR-69 (EnergyBalance -v2), WP-164 (recorder flush/compress)"
tags: [performance, recorder, energy, openmp, wp-196]
updated: 2026-10-14
---

# WP-196 — the commit phase's cost

> [!summary] The short version
> WP-192 left the non-element part of `Domain::commit` as the next floor (~850 ms of a
> 3.75 s step on the rev3 1D window, 46 % of the step at 8 threads). Profiled here, it is
> **almost all `EnergyBalance`** (~807 of 822 ms on the critical rank); the four
> `recorder ladruno` files cost ~15 ms/step together, MPI nothing. Inside the energy sweep,
> per element and step: the `getDamp()` tangent (44 %), the WP-184 "damping outside F"
> check (40 %: a second residual, a second tangent, two more masses), the residual (10 %),
> the mass copy (3 %).
>
> Three changes, **every recorded number unchanged** (byte-identical files; HDF5 datasets
> identical value for value):
>
> 1. the WP-184 check is skipped for elements that declare their damping is exactly
>    Rayleigh (`Element::ladrunoEnergyRayleighOnly`, ASDShellQ4 `-noeas` without a Damping
>    object), for which it always answers "no";
> 2. the per-element terms are formed into one slot per element and summed in the serial
>    order, so the sweep runs on the **WP-192 threads** (branch based on `wp/192`);
> 3. the element mass goes into a per-thread scratch instead of a fresh heap `Matrix`.
>
> Critical rank, 1D window, 2 ranks: **step 4 244 → 3 796 ms at 1 thread (−11 %)** and
> **2 323 → 1 654 ms at 8 threads (−29 %, commit 921 → 194 ms)**.

## 1. Profile — where the commit phase goes

Window: the WP-192 1D rev3 window (`tools.campaign_deck emit 1D --rev rev3 --ranks 2`,
12 968 ASDShellQ4 `-noeas -drillingNL` + 192 forceBeamColumn, Rayleigh αM = 0.098,
βK = 0.0028, MUMPS, KrylovNewton, Newmark, dt 0.0025), recorders **as in the campaign**:
the four `recorder ladruno` at f32 (nodes and mass every step, elements `-T dt 0.0249`,
fibers `-T dt 0.0999`) and `EnergyBalance -v2 -ownedNodes 900` with ten `-region` blocks.
Local i7-13700KF, P-cores only (ranks pinned), shared box (other agents' jobs ran).
New profiler scopes (`profiler start -deep -perStep`, 20 transient steps, rank 1 = critical,
`LADRUNO_EB_LEGACY=1` = the pre-WP-196 work):

| Scope | ms/step (legacy) | after WP-196 |
|---|---|---|
| **commit** | **1 057** | **758** |
| `commit.elements` (element commitState; WP-192 threads it) | 235 | 245 |
| `commit.nodes` | 1.3 | 1.5 |
| `rec.energyBalance` | **807** | **497** |
| — `eb.damp` (`getDamp()`: αM·M + βK·K_t, i.e. one tangent) | 358 | 364 |
| — `eb.check` (WP-184: `getResistingForceIncInertia` = residual + Rayleigh tangent + 2×mass) | 320 | 0.1 |
| — `eb.resid` (`getResistingForce()`) | 79 | 81 |
| — `eb.mass` (`getMass()` + copy) | 26 | 28 |
| — `eb.nodes` (nodal KE / DW / ULW, owned nodes) | 3.5 | 3.4 |
| — `eb.write` (one text row) | 0.1 | 0.1 |
| `rec.ladruno.*` (four files: gather + HDF5 write + deflate + flush) | 14 | 15 |
| MPI | 0 (nothing in commit communicates) | 0 |

Per element (6 538 shells on rank 1): legacy ≈ 123 µs per element-step, of which the
tangent ×2 ≈ 2 × 55 µs. The four Ladruno recorders: `lr.nodes` 0.1–0.6 ms per written
step, `lr.elements` 1.6 ms (elements file) / 3.9 ms (fibers file) per written step, and
`lr.flush` (H5Fflush: the chunk deflate at level 1 and the write) 3–4 ms per flush, every
10 s of wall (WP-164). **Per-step work when not writing: none** — a non-recording step of
the fibers recorder costs ~5 µs (the `-T` gate; nothing before it, the `-G energy` advance
is not requested in these decks). File I/O and compression together are < 1.5 % of the
commit phase, so no recorder change was made (§4 lists what the recording itself costs).

## 2. What ships

| Piece | Where | Notes |
|---|---|---|
| `ebkernel::elementEnergyTerms` | `SRC/recorder/EnergyBalanceKernel.h` | one element's (KE, IE, DW, KE_ms), bitwise what `addElementEnergy` accumulates into zeroed locals; per-thread mass/velocity scratch; `sweepDomain/sweepRegion` (the `-G energy` source) use it too |
| WP-184 check skipped for opted-in elements | same + `Element::ladrunoEnergyRayleighOnly()` (default false) | `LADRUNO_EB_LEGACY=1` forces it everywhere (A/B and verification) |
| ASDShellQ4 opt-in | `SRC/element/shell/ASDShellQ4.{h,cpp}` | no Damping object **and `-noeas`** (see §3.2) |
| Threaded sweep | `EnergyBalanceRecorder::sweepElements/sweepThreads` | snapshot of the live element order; maximal runs of `ladrunoThreadSafeForm()` elements on `ladruno_parallelFor` (Rayleigh slot prepared serially), others alone in place; the reduction loop is the old one over the slots |
| Gates | `sweepThreads` | WP-192's (`ladrunoThreads n`, OpenMP, MPI hybrid acknowledgement, plain Domain, no deep profiler) + `LADRUNO_THREADS_LOOPS` bit 16 (default mask now 31) + no consistent-SMS registry; announced once: `EnergyBalanceRecorder: element sweep threaded on 8 threads: 6538 of 6538 elements formed in parallel, summed in serial order` |
| Profiler scopes | `Domain::commit`, `EnergyBalanceRecorder`, `LadrunoRecorder`, the kernel | `commit.nodes/elements/record`, `rec.energyBalance`, `rec.ladruno.<i>`, `rec.node`, `rec.other.<i>`; `eb.elements/nodes/write` (coarse), `eb.mass/damp/resid/check` (deep); `lr.nodes/elements/flush` |

### 2.1 Why the numbers cannot move

* **Reduction order.** The old loop added each element's terms to the running sums as it
  went. Now each element's terms land in their slot and the same loop adds the slots in the
  same (live iterator) order: same operands, same order, same sums — at any thread count.
  The per-element arithmetic is the old code on zeroed locals (`0.0 + dot` reproduces
  `internal_rate += dot`, `0.0 + gRate` likewise); `KE_ms` is added only when the SMS
  registry is active, as before.
* **The skipped check.** For an element whose damping force is exactly Rayleigh,
  `g = F_inc − F − M a` equals `C v` up to round-off, and the check's tolerance is
  `1e-12 × scale` — it returns "no", and DW books `vᵀ C v`, which is what the fast path
  books. The element answers for that (Element.h contract (a)–(c)).
* **Fewer evaluations.** Skipping the check evaluates the element fewer times after the
  commit. That is invisible only when an evaluation leaves nothing behind that the next
  step reads — contract (c), and the reason EAS shells are excluded (§3.2).

## 3. Verification

### 3.1 Byte-identity

| Configuration | Comparison | Result |
|---|---|---|
| 1D rev3 serial window (`--ranks 1`, Pardiso, `MKL_NUM_THREADS=1`), 30 gravity increments + 5 transient steps, recorders f64, `recorder Node -precision 17` on every node | base (WP-192 binary) vs WP-196 at 1 thread | `disp`, `energy_balance`, `th_steps` **byte-identical**; the four `.ladruno` files: all 257 datasets identical (raw bytes) |
| same | base vs WP-196 at **8 threads** (energy sweep threaded: 13 160 of 13 160) | **identical** |
| 1D 2 ranks (MUMPS) | — | not a byte-identity oracle: the base binary is not run-to-run reproducible under MUMPS (WP-192 §6.1); timing only |
| `tests/test_wp196_energy_commit_phase.py` (new, zone_a): shell wall (`-noeas -drillingNL`) + fiber columns, Rayleigh, `EnergyBalance -v2 -region` at 17 digits | `LADRUNO_EB_LEGACY=1` vs default; 1 vs 2 / 4 threads; with a zeroLength Viscous dashpot (takes the check, books v·g; runs serial in place: "38 of 39"); an EAS shell (not opted in) | 6/6 green, file md5 equal |
| existing: `test_energyBalanceRecorder`, `test_energy_material_dashpot` (WP-184), `test_massScaling_consistent_energy`, `test_ladrunoBrick_hourglass_energy`, `test_ladruno_recorder_{hardening,stage_semantics,write_path}`, `test_wp185_ladruno_schema`, `test_relative_mass_damping`, `test_wp192_threaded_element_loops`, `test_wp107_threaded_update` | WP-196 build | 108 passed |

Compression: no setting changed, so the HDF5 payloads are the same bytes after
decompression (the dataset comparison above reads them decompressed).

### 3.2 Found on the way: ASDShellQ4 with AGQI (EAS) is not evaluation-count neutral

The first cut opted every Damping-free ASDShellQ4 in. The test's EAS arm (the element's
default) failed: legacy and new energy files differed from the first step, and so did the
**trajectory** (KE_nod, nodal displacements). Cause: every `calculateAll()` re-forms the
AGQI condensation state `m_eas->{KQU, KUQ, KQQ_inv, Q_residual}`, an LHS-only call leaves
`Q_residual` zero, and the next `update()` reads it. The legacy sweep ends each element on a
tangent (inside `getRayleighDampingForces`); without the check it ends on a residual — a
different `Q_residual` carried into the next step. (With the legacy sweep the energy
recorder is neutral: same displacements with and without it, checked.) So EAS shells keep
the check; the Paper 1 decks are `-noeas` and get the saving. Quirk fragment
`WP-196-asdshellq4-eas-last-evaluation-state`.

## 4. Speed

1D window, 2 ranks, f32 recorders as in the campaign, 30 gravity increments + 50 transient
steps, coarse `-perStep`, ranks pinned (`I_MPI_PIN_DOMAIN=[1,100]` at 1 thread,
`[FF,FF00]` at 8: cores 0–3 | 4–7 with HT), shared box. Rank 1 (critical), ms/step:

| threads/rank | binary | step | commit | of which EB | window wall |
|---|---|---|---|---|---|
| 1 | base (WP-192) | 4 244 | 1 120 | (~800) | 366 s |
| 1 | WP-196 | **3 796** | **760** | 487 | 356 s |
| 8 | base (WP-192) | 2 323 | 921 | (~750) | 207 s |
| 8 | WP-196 | **1 654** | **194** | 84 | 154 s |

At 8 threads the commit phase falls from 40 % to 12 % of the step; the step is now
`linearSolve` (535 ms) + the element loops. What is left of the sweep at 1 thread is the
tangent inside `getDamp()` (~56 µs per shell) and the residual (~12 µs): both genuinely
needed for DW and IE as defined (the committed-state tangent and force), so only threads
reduce them without changing a number.

Not done, with reasons: **caching the element mass** (`eb.mass` is 28 ms/step at 1 thread,
4 µs per shell, and `getDamp()` re-reads the mass anyway; a cache needs invalidation on
`setParameter`/staged changes for < 3 % of the commit); **caching C** (βK K_t changes every
step in these nonlinear shells; for linear elements this is WP-195's tangent cache, which
`getDamp()` would inherit through `getTangentStiff()`); **skipping zero terms** of M and C
in the quadratic forms (~20 ms/step for every loop and the reduction together; exact only
for finite velocities).

## 5. Options that need the recording to change (NOT done — the author's call)

Measured on the same window (rank 1, ms/step) and the 50-step file sizes (f32, both ranks):

| Option | Saving per step | Changes | Against the extraction |
|---|---|---|---|
| A. Drop `EnergyBalance` from production runs (keep it on a probe subset) | 487 ms at 1 thread (13 %), 84 ms at 8 (5 %) | no energy file | EXTRACTION_DESIGN's energy closure / per-component DW |
| B. Energy at a cadence `-T` (sweep every N steps, rectangle over N·dt) | ≈ (1 − 1/N) × the EB row above | IE/DW/ULW integrated over N·dt — different numbers, closure residual grows | only if the extraction reads cumulative totals, not RES/ERR |
| C. Element-side energy accumulation (issue #959 (b), Radioss design: σ:Δε·vol per GP at commit, external work as a running scalar) | most of EB (est. 450 of 487 ms at 1 thread) | IE becomes Σ σ:Δε (not ∫F·v dt), DW from the element's own Rayleigh force | a new definition: re-validate closure decks (ADR-69) |
| D. Drop the `-region` blocks | ~0 (binning is a hash lookup per element) | fewer columns | no reason to |
| E. Ladruno recorders: `results_mass` / `results_nodes` at a coarser `-T`, or drop `results_fibers` | ≤ 1 ms/step (time); size: over the 50-step window fibers 80 MB, elements 14 MB, mass 9 MB, nodes 0.2 MB | fewer samples | size, not time — the 20 GB per run is fibers and elements |
| F. `-compress 0` | ≤ 3 ms per flush (every 10 s) | none to the numbers; files larger | none (bit-identical decompressed) |

## 6. Follow-ups

* Opt in more element types to `ladrunoEnergyRayleighOnly` (stdBrick, ForceBeamColumn3d
  without Damping) after auditing contract (c) for each — the 3C/4D decks' bricks.
* The same threaded sweep for the Ladruno recorder's `-G energy` source (serial here; the
  campaign decks do not use it).
* An ASDShellQ4 that is evaluation-count neutral with EAS (keep `Q_residual` from the last
  RHS evaluation) would let EAS shells opt in — an upstream-element change, out of scope.
