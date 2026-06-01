---
title: ZFP compression for the Ladruno recorder — feasibility & gains assessment
project: Ladruno
status: research
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

**Where we land (2026-05-31):** a **qualified yes — but not built yet**. The case *for* is
solid: it's an **opt-in flag, not a new recorder** (sibling to `-precision f32`); performance is
*favorable* not a trade (ZFP is ~20–60× faster than the gzip stage it replaces → tends to be
faster **and** smaller); and the read-side plugin dependency — the one real downside — is
**accepted** (see [§Resolved](#risks--open-questions)). The case for *waiting*: the cheap 2× is
**already banked** by `f32+gzip` at zero cost, ZFP's *marginal* gain on the **current frozen
layout** is only **~1.5–3×**, the headline 3–6× needs a `FORMAT_VERSION` v2, and there's real
correctness (NaN/Inf) + Windows-build work. The magnitude is uncertain and that uncertainty lives
entirely in *our data + layout*.

**Move:** run the [benchmark](#benchmark-do-this-first) **before** any C++ — and it needs **zero
recorder code / zero build risk** (re-encode an existing f64 `.h5` run's `DATA` arrays in Python
via `hdf5plugin`'s ZFP at a few tolerances; measure size + error vs the f64 oracle). That
half-day offline spike answers the one number the whole decision hinges on. If ZFP beats
`f32+gzip` by ≳2× at acceptable error → implement the `-compression zfp` flag (fixed-accuracy
default, per-quantity tolerance, NaN/Inf-guarded, never default-on). If <1.5× → the frozen layout
is the bottleneck; stick with `f32+gzip` until/unless schema v2 is on the table. **Decision
deferred until the benchmark runs.**

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
- Applies only to the per-step `DATA` arrays — envelopes/coords/metadata untouched, same as f32.
- More than f32 needed: richer flag parsing (mode + tolerance) and a NaN/Inf guard / reversible
  fallback for lossy modes.

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
(a nonlinear transient with many steps; a fiber-section element history; a large nodal field):

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

> [!question]
> Can result `DATA` ever contain NaN/Inf (diverged steps)? If yes, the lossy path must scrub or
> fall back to reversible — confirm before shipping lossy-by-default anything.

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

*(empty — research phase)*
