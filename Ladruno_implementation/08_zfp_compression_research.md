---
title: ZFP compression for the Ladruno recorder — feasibility & gains assessment
project: Ladruno
status: decided — REJECTED (benchmark 2026-07-06, see Implementation log)
priority: medium
owner: nmora
tags:
  - research
  - recorder
  - compression
---

# ZFP compression for the Ladruno recorder — research assessment

> Research-phase decision doc. **No code yet.** Goal: decide whether error-bounded
> floating-point compression (ZFP) is worth adopting in the Ladruno recorder, what it
> would cost architecturally, and whether it is the best option. Sources cited inline.

## TL;DR verdict

ZFP is **technically sound, BSD-licensed, mature, and feasible** to wire in via the
LLNL **H5Z-ZFP** HDF5 filter (filter id `32013`). But three facts reshape the decision:

1. **We are not compressing raw f64.** The recorder *already* ships `shuffle + gzip(4)`
   on every chunked `DATA` dataset, **and** an opt-in `-precision f32` lossy mode
   (measured **0.496×**, worst relative error 4.2e-7, zero violations). ZFP must be
   judged on its **marginal** gain over `f32 + shuffle + gzip`, not over f64. The
   headline "100:1" ZFP numbers are against uncompressed f64 and do not transfer.
2. **The headline ratio is gated by our data layout, which is frozen.** ZFP's gains
   come from spatial/temporal correlation inside 4ᵈ blocks. Our `DATA[T×nIds×nComp]`
   layout interleaves components on the fast axis (the exact anti-pattern ZFP docs warn
   against) and orders `nIds` arbitrarily. Realizing ZFP's potential implies restructuring
   the schema — which collides with the **frozen `FORMAT_VERSION=1`** and the apeGmsh-native
   contract.
3. **The biggest cost is read-side, not write-side.** `gzip`/`shuffle` are HDF5 **built-in**
   filters — every `h5py` / `STKO_to_python` / apeGmsh reader decompresses transparently
   with zero install. ZFP is **not** built in: any downstream reader must have the filter
   present (`pip install hdf5plugin`). The recorder's whole reason to exist is to be read by
   that Python tooling. `-precision f32` has none of this problem (native f32 reads everywhere).

**Where we land (2026-07-06, superseding the 2026-05-31 "qualified yes"): REJECTED.**
The [benchmark](#benchmark-do-this-first) ran (harness + inputs in
`Ladruno_scripts/zfp_benchmark/`, adversarially reviewed before interpretation) and the answer
is unambiguous: **no ZFP configuration beats `f32+shuffle+gzip` at matched fidelity — on the
frozen layout OR on ZFP's preferred component-separated/time-major layout.** Best
fidelity-matched marginal: **0.65×** (current layout) / **0.97×** (schema-v2 layout) — a 3×
shortfall against the ≥2× adopt bar, consistent across three models spanning ZFP's optimistic
case (smooth 3D wave field) to its pessimistic one (fiber hysteresis, where `f32+gzip` hits
~13:1 and ZFP is 4–8× *bigger*). ZFP-reversible *expands* past raw f64 on all three models.
SZ3 also loses everywhere it was measured. Even the decision rule's fallback presumption died:
the frozen interleaved layout is *favorable* to the shipped codec, so there is **no
compression-motivated FORMAT_VERSION v2** — the question is closed, not deferred. Full numbers,
eight mandatory caveats, and reopening criteria in the [Implementation log](#implementation-log).

---

## What we gain

- **Error-bounded lossy ratio beyond f32.** Fixed-accuracy ZFP truncates to a chosen
  *absolute* error tolerance. On smooth correlated FP64 the literature reports up to ~100:1
  vs raw f64; realistic engineering tolerances on transient histories land in the **4–10×**
  range vs raw f64, i.e. plausibly **~2–4× beyond our current `f32+gzip`**. (zfp docs / LLNL.)
- **A principled error knob.** Unlike f32 truncation (uncontrolled ~7-digit *relative* error),
  fixed-accuracy gives a stated bound: `|f−g| ≤ 2^floor(log2(tol))`, typically several×
  tighter in practice. This is a better story for "publication-grade but smaller" output.
- **A lossless option too.** ZFP **reversible** mode is bit-exact (and the only mode that safely
  handles ±Inf/NaN), giving ~1.5–4× lossless — a drop-in for users who refuse any loss.
- **Speed is a non-issue.** ~1–2 GB/s/core serial/OpenMP, "6–170× faster" than gzip. It will
  not be the I/O bottleneck; gzip(4) already is the slower stage today.
- **Permissive license.** zfp = BSD-3-Clause, H5Z-ZFP = BSD, hdf5plugin wrapper = MIT.
  No copyleft/redistribution concerns for the fork.

## What we lose / risk

- **Read-side hard dependency (the big one).** Files written with ZFP are **unreadable** by any
  HDF5 reader lacking filter 32013. Python consumers must `pip install hdf5plugin; import
  hdf5plugin`. This breaks the "any h5py reads it" property that `gzip`/`f32` preserve, and the
  canonical consumer is exactly that Python tooling (apeGmsh, STKO_to_python).
- **Layout-limited gains under the frozen schema.** Our fast axis is `nComp` (interleaved,
  e.g. `ux,uy,uz,rx,ry,rz`) and `nIds` order is arbitrary. ZFP shares one exponent per 4ᵈ block,
  so mixing unrelated components/ids in a block **caps the ratio well below the headline**.
  Time-local smoothness *is* captured along the chunk's slow axis, so it is not zero — but the
  best wins need separated scalar datasets + a long time axis, i.e. a schema change.
- **NaN/Inf are UB in lossy modes.** Diverged/uninitialized steps in a structural solve can emit
  non-finite values; in fixed-rate/precision/accuracy modes that is undefined behavior and a
  single sentinel poisons its whole 4ᵈ block. Requires scrubbing or a reversible fallback.
- **Absolute-error trap.** Fixed-accuracy bounds *absolute*, not *relative*, error. A dataset
  mixing small (1e-3 displacements) and large (1e6 forces) magnitudes loses all relative
  significance on the small values. Our schema already keeps result *families* separate, which
  helps, but per-quantity tolerances would be needed.
- **Build/toolchain surface.** A new external dep (zfp + H5Z-ZFP) into the
  Windows/Intel oneAPI/Conan/MUMPS toolchain that has historically been finicky. H5Z-ZFP CMake
  build requires zfp built with CMake; there's a documented HDF5≥1.12
  `--disable-memory-alloc-sanity-check` quirk for library-mode use.

## Performance & size (concrete numbers)

The figures below separate **measured** from **estimated**; the [benchmark](#benchmark-do-this-first)
is what pins the estimates to *our* models.

### Size — the thing to beat is `f32+gzip`, not raw f64

| Config | File size vs f64 | Status | Fidelity |
|---|---|---|---|
| f64 + shuffle + gzip(4) | ~0.6–0.8× | est. (lossless, smooth FP64) | exact |
| **f32 + shuffle + gzip(4)** | **0.496×** | **measured** (commit d88d51349) | ~7 digits, worst rel err 4.2e-7 |
| ZFP fixed-accuracy, *current* layout | ~0.15–0.30× | est. — **needs benchmark** | tunable bound |
| ZFP fixed-accuracy, *separated/time-major* layout | ~0.08–0.20× | est. — needs schema v2 | tunable bound |
| ZFP reversible (lossless) | ~0.35–0.65× | est. | exact |

Marginal win that matters — **ZFP over `f32+gzip`**: ~**1.5–3×** on the current
(interleaved-component, arbitrary-node-order) layout; ~**3–6×** only with a restructured schema.
The spread is wide because ZFP's ratio is entirely data-dependent: smooth/low-frequency histories
compress hard, noisy/near-impulse content barely at all. A single number is meaningless without
real models — hence the benchmark.

### Performance — favorable, and firmer than the size numbers (set by the codec, not the data)

| Stage | Current (gzip-4) | ZFP | Net |
|---|---|---|---|
| Write / compress | ~30–80 MB/s/core | **1–2 GB/s/core** (~20–60×) | ZFP-only **speeds up** writes; ZFP+gzip ≈ neutral |
| Read / decompress | ~200–400 MB/s | **1–2 GB/s** | faster reads |
| Disk I/O | dominated by file size | smaller files | **less I/O = faster end-to-end** on big runs |
| Memory | chunk buffer | + small per-chunk scratch | negligible |
| Read-side init | none | one-time `H5Z_zfp_initialize()` / plugin load | negligible |

Key point: **gzip(4) is currently the slowest stage in the write path**, and ZFP is ~20–60×
faster per byte — so adopting ZFP isn't the usual speed/size trade; it tends to be *faster AND
smaller*. The cost is **not** CPU — it's the read-time plugin dep (accepted) + the lossy
correctness work. Caveat: the recorder appends one `[1×nIds×nComp]` slab per step; for very small
slabs (few DOFs) ZFP's per-call overhead shrinks the throughput edge, but absolute time is trivial
there anyway.

## Integration shape — an option flag, **not** a new recorder

ZFP rides the exact rails `-precision f32` already proved; no new recorder, class, or schema break.

```
-compression zfp:accuracy=1e-6   →  ProcessInfo.compression{mode,tol}
                                  →  INFO/COMPRESSION attr (additive, like STORED_PRECISION)
                                  →  StreamingSink::begin() passes it to createTimeSeries3d()
                                  →  H5Pset_filter(ZFP) instead of/alongside H5Pset_deflate
```

- **Orthogonal & composable with `-precision`:** precision picks the on-disk *type* (f32/f64),
  compression picks the *filter* (none/gzip/zfp). `-precision f32 -compression zfp` = ZFP over
  f32; `-precision f64 -compression zfp:accuracy=1e-6` = error-bounded ZFP over doubles. ZFP
  handles both single and double.
- **Scope (updated 2026-07-06 — three exclusions, all structural):**
  - Applies only to the chunked per-step `DATA` time-series (the `createTimeSeries3d` path) —
    coords/metadata untouched, same as f32.
  - **Envelope mode: force-disabled, mirroring f32.** Since this doc was written, `-precision f32`
    is force-disabled in envelope mode with a warning (`LadrunoRecorder.cpp` ~L530 — envelope
    datasets are always f64). `-compression zfp` MUST adopt the same contract: envelope mode →
    warn + stamp `COMPRESSION=none`. Structurally right anyway: `EnvelopeSink` delete-recreates
    its small leaf datasets on every periodic flush — worst-possible fit for a per-chunk codec.
  - **Modal/eigen output: excluded by construction.** The eigen fix (post-doc) made
    `recordModeChannel` bypass the sink entirely — modal `DATA` is a *group* of small contiguous
    `MODE_<k>` datasets (no dcpl, no filters, `LadrunoRecorder.cpp` ~L1823). ZFP never touches it.
  - **Monitor sink: MUST stay filter-free.** `LadrunoMonitorRecorder`/`LadrunoMonitorSink`
    (#485/#487, post-doc) is a separate SWMR live-tail file whose datasets are chunked but
    deliberately uncompressed — `profiler_monitor.py` and the `monitor_viewer` dashboard tail it
    with plain h5py, no install. Do NOT wire `-compression` through it "for consistency"; that
    would break the live-monitoring story.
- More than f32 needed: richer flag parsing (mode + tolerance) and the NaN/Inf guard below.

## Dependencies — write side compiles it, read side pip-installs it

- **Write side (OpenSees build):** link **zfp + H5Z-ZFP** (two small BSD libs) — *library mode*
  recommended (filter baked into the binary + a one-time `H5Z_zfp_initialize()`; nothing extra to
  ship, no env var; clean fit for the Inno installer). *Plugin mode* (ship a `.dll` + set
  `HDF5_PLUGIN_PATH`) avoided — more moving parts on Windows. It is the **OpenSees compilation**
  that carries the dep, because the *recorder* writes — not apeGmsh.
- **Read side (any consumer):** `pip install hdf5plugin` + `import hdf5plugin`, then h5py reads ZFP
  transparently. Notes: (a) it's **not only apeGmsh** — STKO_to_python and any h5py script need it
  too, so apeGmsh should declare `hdf5plugin` a hard requirement; (b) **no version-lock** between
  the two sides — hdf5plugin bundles its own zfp build and the ZFP bitstream is portable across
  versions; (c) the dependency is **conditional on opt-in** — files written *without*
  `-compression zfp` stay readable by plain h5py with zero plugin. ZFP's cost only "infects" files
  that used it.

## Architectural cost

The recorder's write path makes the *mechanical* insertion cheap; the *correctness/UX*
work is where the cost lives.

- **Chokepoint already exists.** Every float write funnels through
  `ladruno::h5::dataset::createTimeSeries3d()` / `appendSlab3d()` in
  [`SRC/recorder/Ladruno_Hdf5.h`](../SRC/recorder/Ladruno_Hdf5.h). Filters are set in
  `createTimeSeries3d()` (today: `H5Pset_shuffle` + `H5Pset_deflate(4)`). ZFP is one more
  `H5Pset_filter()` call here.
- **Chunking is already on** (required by ZFP) — `H5Pset_chunk` targets ~256 KiB chunks.
  No new buffering layer needed; ZFP works per-chunk transparently.
- **Opt-in precedent is in place.** `-precision f32` (parsed at
  [`LadrunoRecorder.cpp`](../SRC/recorder/LadrunoRecorder.cpp), `store_data_f32` in
  `ProcessInfo`, stamped as `INFO/STORED_PRECISION`) is the exact template: add a parallel
  `-compression zfp[:accuracy=…]` flag + an `INFO/COMPRESSION` attr.
- **Real work items:** (a) CMake wiring of zfp + H5Z-ZFP for the Windows/oneAPI build;
  (b) NaN/Inf guard or reversible fallback; (c) per-quantity tolerance plumbing; (d) extend the
  bounded-error gate (`Ladruno_scripts/ladruno_recorder_tests/precision_check.py`) to validate
  ZFP bounds; (e) **decide and document the read path** (bundle/document `hdf5plugin`).
- **Schema tension:** maximal ratio wants component-separated scalar datasets + time-major
  layout. That is a `FORMAT_VERSION` bump and an apeGmsh-contract change — out of scope for a
  v1 opt-in; note it as a v2 layout question.

## Is this the best we can do?

| Option | Lossy? | vs ZFP | Read portability |
|---|---|---|---|
| **f32 truncation** (have it) | yes (~2×) | strictly dominated by ZFP fixed-rate@32, but trivially portable | **native** (no plugin) |
| **gzip+shuffle** (have it) | no (~1.5–2×) | far lower ratio, slower | **built-in** |
| **ZFP** | yes/lossless | best-*balanced*: ratio + bound + license + python path | needs filter (`hdf5plugin`) |
| **SZ3** | yes | **higher ratio at equal error bound**, comparable speed | needs filter (`hdf5plugin`) |
| **SPERR** | yes | highest ratio on smooth data, **~3× slower**, newer | weaker HDF5 support |
| **Blosc2+bitshuffle** | no | fast lossless, modest ratio | needs filter (`hdf5plugin`) |
| **fpzip** | lossless/prec | simpler, lower lossy quality, no abs-error bound | `hdf5plugin` |

**Answer:** for a recorder valuing turnkey-ish integration + good lossy ratio + an error bound +
permissive license + a Python read path, **ZFP is the best-balanced** choice. **SZ3** edges it on
pure ratio-at-bound if we accept a heavier dep. But note the cheapest 2× (f32) is **already
banked at zero portability cost** — so the live decision is whether ZFP's marginal 2–4× justifies
a hard read-time plugin dependency.

## Feasibility

**High, mechanically.** zfp 1.0.1 (2023-12) and H5Z-ZFP 1.1.1 both build with CMake, both list
Windows in CI, both are pure-C/permissive with HDF5 as the only real dep. The recorder already has
the chunked-dataset + filter-chain + opt-in-flag scaffolding. Python read path is solved by
`hdf5plugin`. The non-trivial parts are the Windows/oneAPI/Conan build wiring and the NaN/Inf
correctness guard — both bounded, known-shape work.

## Maintainability

**Moderate cost, mostly at the seams:**

- +1 external dependency pair to track/upgrade in the (historically finicky) Windows build.
- A new lossy code path needs the bounded-error gate extended and kept green — but that machinery
  exists for f32.
- The read-time `hdf5plugin` requirement is a **documentation + support** burden forever (every
  "I can't open the .h5" ticket). This is the durable maintenance tax, not the C++.
- If we ever restructure the schema for better ratios, that's a `FORMAT_VERSION` migration with
  reader fan-out — a larger, separate commitment.

## Benchmark (do this first)

Before any integration, settle the *marginal* gain empirically on 2–3 representative models
(a nonlinear transient with many steps; a fiber-section element history; a large nodal field).
**Input files are free (2026-07-06):** the regression battery
(`Ladruno_scripts/ladruno_recorder_tests/run_regression.bat`, 19 gates) already produces real
`.ladruno` files from real models — the re-encode spike consumes those (or files freshly written
by the *installed* Ladruno build via the venv) instead of inventing fixtures. Configs:

1. `f64 + shuffle + gzip(4)`  ← current lossless baseline
2. `f32 + shuffle + gzip(4)`  ← current lossy baseline (the real thing to beat)
3. `ZFP fixed-accuracy` at 2–3 tolerances, current layout
4. `ZFP reversible` (lossless ceiling)
5. *(optional)* `ZFP fixed-accuracy` on a component-separated / time-major prototype layout
   — to quantify how much the frozen layout is costing us.

Record, per config: file size ratio, worst absolute & relative error vs f64 oracle, write
wall-time. Decision rule: adopt only if (3) beats (2) by a margin that justifies the read-time
plugin dependency for our users.

## Risks / open questions

> [!check] RESOLVED (2026-05-31)
> A hard `hdf5plugin` read dependency **is acceptable** — "any h5py opens it" is *not* a hard
> requirement. ZFP stays on the table. (Downstream Python tooling will document/require
> `pip install hdf5plugin`.) The benchmark still gates whether the marginal gain over `f32+gzip`
> justifies adoption.

> [!check] DECIDED (2026-07-06) — NaN/Inf guard is a requirement, not an open question.
> The recorder only records on commit, but a falsely-converged step can still carry non-finite
> values, and one sentinel poisons its whole 4ᵈ ZFP block (UB in lossy modes). Since ZFP runs at
> ~1–2 GB/s, a per-slab `std::isfinite` scan before compress is noise. Requirement: **scan each
> slab; on any non-finite value, warn once per dataset + fall back to reversible mode for that
> dataset.** No scrubbing (never rewrite user data).

> [!question]
> Do we accept the frozen-layout ratio ceiling, or is a v2 component-separated/time-major schema
> on the table? The benchmark's config (5) answers whether that's worth a `FORMAT_VERSION` bump.

## Sources

zfp modes/FAQ/install (zfp.readthedocs.io, github.com/LLNL/zfp), LLNL zfp & floating-point-compression
project pages, H5Z-ZFP docs/repo (h5z-zfp.readthedocs.io, github.com/LLNL/H5Z-ZFP),
hdf5plugin (silx.org/doc/hdf5plugin), error-bounded-compressor survey (arXiv:2404.02840).
Codebase refs: `SRC/recorder/Ladruno_Hdf5.h`, `Ladruno_Sinks.cpp`, `LadrunoRecorder.cpp`,
`Ladruno_implementation/ladruno_schema_v1.md`.

## Implementation log

### 2026-07-06 — benchmark ran; DECISION: REJECT `-compression zfp` (and SZ3)

**Harness:** `Ladruno_scripts/zfp_benchmark/` (`reencode_bench.py` + three model
generators + `check_inputs.py` yield gate + `results.json`). Zero recorder code — re-encodes
the f64 `DATA` arrays of real `.ladruno` files (written by the installed build) under 16
configs. Methodology was **adversarially reviewed (Opus 4.8) before any result was
interpreted**; the review caught and we fixed three blockers: (1) the f32 baseline must be
re-encoded with *recorder-true* chunking (chunk budget is in on-disk bytes —
`Ladruno_Hdf5.h::chunkSlabsPerBlock` — so real f32 chunks span 2× the time-slabs of f64;
copying f64 chunk geometry under-measures the baseline and inflates ZFP's margin on exactly
the gate metric); (2) the time-major layout prototype needs gzip controls (full 2×2
codec×layout matrix) or layout-vs-codec is unattributable; (3) the gate must be per-family
worst-case, not aggregate, so one big smooth nodal field can't mask element-history losses.

**Inputs** (all gated by `check_inputs.py`): **A** nonlinear J2 2D plate, *broadband*
random-phase 0.5–25 Hz load, 3000 steps, yielded (191.8 MB raw); **B** RC fiber column
(dispBeamColumn 8×3 IP, Concrete02+Steel02), cyclic to ductility ≈1.1, 2000 steps
(199.0 MB raw — fiber stress/strain are 756-column datasets); **C** 3D elastic brick
Ricker wave field, 2601 nodes, 600 steps (112.4 MB raw) — ZFP's optimistic case.

**Result — marginal ratio vs `f32+shuffle+gzip(4)` (>1 = smaller than baseline; adopt bar
was ≥2×), "aggregate | worst-family":**

| Config | A | B (fiber) | C |
|---|---|---|---|
| baseline `f32_gzip4` (abs) | 0.300× raw, rangerel 5.5e-8 | **0.077× raw (13:1)** | 0.187× raw |
| ZFP acc 1e-6, frozen layout | 0.63 \| 0.21 | 0.22 \| 0.14 | 0.56 \| 0.42 |
| ZFP **f32-matched** (rr 5.5e-8), frozen | 0.65 \| 0.32 | 0.18 \| 0.17 | 0.47 \| 0.46 |
| ZFP f32-matched, **time-major (v2)** | 0.97 \| 0.71 | 0.20 \| 0.17 | 0.68 \| 0.62 |
| ZFP rr 1e-6, time-major | 1.23 \| **0.85** | 0.24 \| 0.21 | 0.86 \| 0.75 |
| ZFP reversible | 0.26 (1.17× raw!) | 0.07 | 0.17 |
| SZ3 abs 1e-6 | 0.57 \| 0.30 | 0.47 \| 0.29 | 0.55 \| 0.33 |
| `f32_gzip4` on time-major layout | 0.76 | 0.20 | 0.87 |

**Grounds for REJECT:** (1) at f32-matched fidelity on the shipped layout ZFP's best is
0.65× — a 3× shortfall vs the adopt bar; (2) even granting a FORMAT_VERSION v2 layout AND
oracle-assisted per-quantity tolerances, best f32-matched is 0.97× — still a loss, so
**v2-for-compression is closed permanently** (the "<1.5× ⇒ layout is the bottleneck"
presumption is *refuted*: ZFP loses on its preferred layout too, and the interleaved frozen
layout is actively favorable to shuffle+gzip — decisively so for fiber families, where
cross-fiber byte redundancy gives gzip ~13:1 and component separation destroys it); (3) at
every reachable size point, `f32+gzip` offers equal-or-better worst-case range-relative
error in a smaller file — there is no error-size operating point where ZFP dominates; (4)
the sole >1 aggregate (1.23×, one model, 5–7× worse realized error, v2 required) fails its
own worst-family gate at 0.85×.

**Mandatory caveats (from the results-interpretation adversarial pass, all
verdict-preserving):**
1. *rr floor:* ZFP's power-of-2 tolerance flooring made rr=5.5e-8 rows realize ~1e-8
   rangerel (4–5× better than f32); interpolated true-fidelity-matched ≈ 0.72–0.75×
   (frozen) / ~1.05× (tm) — still decisive losses.
2. *Oracle-assisted rr:* per-dataset tol used the posterior whole-run max, which a
   streaming appender can't know — the benchmark bounds real per-quantity plumbing from
   ABOVE.
3. *Mixed-unit datasets:* rr matches per-dataset range-relative error only; small
   components in mixed-unit families keep worse relative fidelity under ZFP-rr than under
   per-value f32.
4. *SZ3 rows are not fidelity-matched* (global abs bound: over-tight on stress, loose on
   strain); SZ3 is rejected on the per-family bound + 10–20× slower encode + same plugin
   dependency, not those size numbers. Re-proposal requires sz3_rr + sz3_tm_rr,
   per-family-gated.
5. *B's 13:1 is the upper end* (mild ductility μ≈1.1 + dispBeamColumn smoothness inflate
   gzip's cross-fiber win); deep hysteresis would cut it (plausibly 4–6:1) but degrades
   ZFP's transform at least as much — ranking robust, magnitude not.
6. *tm-gzip magnitudes* (0.76–0.87×) used untuned (64×1024) chunks — sign is
   mechanism-backed, magnitude untuned; tm configs should be quoted by file size
   (storage-sum undercounts ~6% metadata on B's 756-dataset tm layout).
7. *Timing is Python wall-clock end-to-end* — must not be quoted for/against the
   literature "20–60× faster than gzip" claim (observed: ZFP ~2–4× faster than f64-gzip
   encode; SZ3 ~10× slower; order-of-magnitude only).
8. *Reopening criteria:* a NEW codec class (predictor-based, fidelity-matched, per-family
   gated) — not a layout change, not a tolerance retune.

**Standing decision: keep `f32 + shuffle + gzip(4)` as the only lossy mode. No new
dependency, no plugin read-path, no schema change.**
