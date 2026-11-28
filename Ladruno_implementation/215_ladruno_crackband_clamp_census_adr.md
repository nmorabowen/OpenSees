---
title: "ADR 215 — Crack-band snap-back clamps: a per-tag run-time census in every crack-band concrete (WP-218)"
project: Ladruno
type: ADR
status: "Proposed (2026-11-01) — design only, no C++"
priority: medium
owner: nmora
related:
  - "apeGmsh ADR 0129 (nmorabowen/apeGmsh#2661, Proposed): the fork reports, apeGmsh predicts, a live test reconciles. This ADR is its D1."
  - "apeGmsh ADR 0044 (ASDConcrete regularization contract), apeGmsh #2648 (DR-6 calibration), #2659 (the emit-time sweep)"
  - "WP-99 (commit-time refusal seam, `LadrunoMaterialStatus.h`), WP-104 (`wipe` is the zero point of a process-wide census), WP-107 (threaded update/commit; per-model warn-once latches), WP-152 (count transitions at commit)"
  - "WP-168 (`LadrunoCommandTable.h`, one row registers a command in all three interpreters)"
tags: [adr, concrete, crack-band, regularization, snap-back, diagnostics, wp-218]
updated: 2026-11-01
---

# ADR 215 — Crack-band snap-back clamps: a per-tag run-time census

## 1. Context

- **What a clamp is.** A crack-band concrete spreads its post-peak energy
  over the element length `lch`, rescaling the fracture energy as
  `g = G · lch_ref / lch` (or `G / lch`). Past a ceiling `l_max` the rescaled
  energy would fall below the elastic energy at peak, so the stress-strain
  curve would snap back. The material then floors the energy (or, in one
  case below, does nothing at all). The element is more brittle than the
  material and the answer depends on the mesh. Nothing fails, and nothing in
  the results shows it.
- **Who knows.** Only the binary has the exact `lch` (per integration point,
  from the element, curved geometry and beam tributary lengths included) and
  knows whether the point ever softened. apeGmsh mirrors each element's
  `getCharacteristicLength` in Python to warn at emit time (apeGmsh #2659).
  apeGmsh ADR 0129 makes the fork the authority (its D1) and adds a live test
  that compares the two (its D3). That test needs a count it can read from
  the fork, and the fork has none today.

### 1.1 Today's state, read from the source

FORK_PIN `23d937bc4` is an ancestor of `origin/ladruno` (`f37e2ff19`). The
clamp code itself did not move between them. `LadrunoConcrete3D.cpp` grew by
~160 lines, so its citations differ. Each row gives `pin` / `HEAD` lines.

| Material | Side | Where it clamps | Criterion | Reported? |
|---|---|---|---|---|
| `LadrunoConcrete3D` | compression | `compressiveEpsFc()`, `LadrunoConcrete3D.cpp` 702-720 / 862-880. The clamp itself is `epsFcFromGc()`, `LadrunoConcrete3DKernel.h:3300-3306`, status `-1` | `Gc/lch` below the table's smallest energy `g[0]` (`eps_fc = 0.01 fc/E`, `GC_TABLE_LO`, kernel `:3256`); `eps_fc` is clamped to `efc[0]` | **Yes, but per integration point, not per tag.** The `opserr` at 711-716 / 871-876 is latched by `gcWarned`, a member that every copy constructs `false` (pin 357 / HEAD 357, 387). `getCopy` (HEAD 885-890) copies the table but not the latch, so each Gauss point prints the line once. The message is printed from `integrate()`, which can run under the WP-107 threaded update. It also fires when `lch` is first seen, whether or not the point ever softens. |
| `LadrunoConcrete3D` | tension | **No clamp exists.** The bilinear law is `solveOmegaBilinear()` (kernel `:613-640`), with `wf = Gf/(0.225 ft)` and opening `w = lch·(κ1 + ω κ2)`; `omegaT()` (`:642-647`) caps ω only at `OMEGA_MAX`. The legacy law is `solveOmegaBracketed(…, Gf/(ft·lch))` (`:644`) | `lch > l_max,t = k·E·Gf/ft²`, `k = 20/21` bilinear, `k = 1` exponential (apeGmsh `tension_snapback_factor`) | **No.** Nothing caps the energy and nothing is said. What the ω solve returns past `l_max,t` is not characterized. S1's oversized test records it. |
| `ASDConcrete3D` | tension and compression (`ht`, `hc`) | `HardeningLaw::regularize()`, `ASDConcrete3DMaterial.cpp` 985-1040 (same at pin and HEAD); the floor is `gmin = 1.01·peak.y·peak.x/2`, `gnew = max(gnew, gmin)` at **1002-1003**. Called once per instance from `setTrialStrain`, 1685-1693 | `g0·lch_ref/lch < gmin` | **No.** Silent. apeGmsh ADR 0044 cites `:961-965`, which is a stale line number at the pin. |
| `ASDConcrete1D` | tension and compression | `HardeningLaw::regularize()`, `ASDConcrete1DMaterial.cpp` 640-695; floor at **657-658**. Called from `setTrialStrain` 1057-1065 | same | **No.** Silent. In a fiber section `lch` is the beam's per-IP length (`ForceBeamColumn3d.cpp:1183-1188`). |
| `LadrunoRCConcrete` | tension and compression | `ladruno_rc_kernel::regularize()`, `LadrunoRCKernel.h` 253-283 / 254-284; floor at **260-261 / 261-262**. Called by `regularizeIfNeeded()`, `LadrunoRCConcrete.cpp` 353-354 / 374-375 | same rule as ASD (a faithful clone) | **No.** Silent. `regWarned` (HEAD 263, 285) is only for the "no `lch` available" refusal. |
| `LadrunoRCFiniteStrain` | tension and compression | same kernel; `LadrunoRCFiniteStrain.cpp` 280-281 / 300-301 | same | **No.** Silent. |

Three facts follow from the table and shape the design:

1. **Every clamp is decided once per instance**, when the point first sees
   its `lch`. The decision is "the floor is in force at this point". It is
   not "this point softened". A point can be clamped and never crack. So the
   census has two levels: **armed** (the floor is in force) and **realized**
   (an armed point committed a softening state on that side).
2. **The criteria are not the ones apeGmsh mirrors.**
   - ASD and RC floor at `1.01×` the peak triangle. Their `l_max` is
     `l_ref·g0/gmin`, which for ASD tension is `2EGf/(1.01 ft²)`. apeGmsh
     uses `2EGf/ft²`, so an `lch` within 1 % below apeGmsh's ceiling is
     clamped by the fork and not flagged by apeGmsh.
   - `LadrunoConcrete3D` compression clamps at the table floor `Gc/g[0]`.
     That is a calibration-range limit, not the analytic
     `2EGc/(fc²(1+Hp))` that apeGmsh uses.
   - apeGmsh does not check the ASD compression side at all.
   - D3 would find all three. The census also reports the fork's own
     `l_max` per tag and side (D2), so the live test can check apeGmsh's
     formula as well as its `lch` rule.
3. **The update and the commit can both run threaded** (WP-107:
   `Domain::ladrunoThreadedUpdate`, `Domain.cpp:2609`, and
   `Domain::ladrunoThreadedCommit`, `:2840`, which runs `commitState()`
   under `ladruno_parallelFor`, `:2896`). `ops_TheActiveElement` is
   `thread_local` (H7, `:2603`). Note that `LadrunoMaterialStatus.h:112-116`
   still says the commit phase is serial; that comment is stale. Any shared
   write here must be locked.

### 1.2 The existing counters, and which pattern this reuses

- **`ladrunoCommitRefusalCounter()`** (`LadrunoMaterialStatus.h:117-155`) is
  per-commit. `Domain::commit()` clears it after acting on it. It has no
  read path, so it is not a census.
- **`LadrunoImplexGlobals`** (`LadrunoSANISAND.cpp`, around 2900-3015) *is* a
  census: a process-wide singleton, zeroed on `wipe` by
  `ladrunoSanisandResetImplexGlobals()` from `OPS_clearAllNDMaterial()`
  (`NDMaterial.cpp:110-119`). `ops.reset()` / `revertToStart()` leave it
  alone, so it is read as deltas inside one model (WP-104 quirk).
  It is read through any one Gauss point's material response
  (`implexRefusals`, `implexGuards`, `LadrunoSANISAND.cpp:5483`, `:5529`).
- **This ADR reuses its lifetime rule exactly.** Wipe is the zero point;
  reset/revertToStart are not; the hook is the family's `OPS_clearAll*`, not
  `Domain::clearAll()`, which also runs on an MP rank's `recvSelf`. It does
  **not** reuse the read-through-a-point channel (D5 says why).
- **It reuses the header-only shape of `LadrunoMaterialStatus.h`.** An inline
  function with a function-local static gives one instance across TUs, with
  no new library edge. That matters because the recorder library
  (`OPS_Recorder`) and the interpreters must read it too.

## 2. Decisions

### D1 — one census, keyed by (material tag, side)

A new header, `SRC/material/LadrunoCrackBandCensus.h`, holds one
process-wide registry. Entries are keyed by `(matTag, side)`, with side ∈
{tension, compression}. Each entry holds:

| Field | Meaning |
|---|---|
| `className` | `"ASDConcrete3D"`, `"ASDConcrete1D"`, `"LadrunoConcrete3D"`, `"LadrunoRCConcrete"`, `"LadrunoRCFiniteStrain"` |
| `nArmed` | integration points whose floor is in force (counted once per point, D3) |
| `nRealized` | armed points that committed a softening state on this side (once per point) |
| `worstRatio` | max `lch / l_max` over the armed points (> 1 by construction) |
| `worstEle` | the element tag at the worst point (`ops_TheActiveElement`, -1 without one) |
| `lMax` | the fork's own `l_max` for this tag and side: `lch`-independent for every material here, so one number |
| `warned` | the per-tag warning latch (D4) |

Only armed points create entries, so a clean run leaves the registry empty.

### D2 — what each material counts as armed, and its `l_max`

Each is computed where the clamp is already decided, once per instance,
with no new work per iteration.

| Material / side | Armed when | `l_max` | `lch / l_max` |
|---|---|---|---|
| ASD 3D/1D, RC, RC-FS (both sides) | `gnew_raw = g0·lch_ref/lch < gmin` in `regularize()` | `lch_ref·g0/gmin` | `gmin / gnew_raw` |
| `LadrunoConcrete3D` compression | `epsFcFromGc` status `-1` | `Gc / g[0]` | `lch·g[0] / Gc` |
| `LadrunoConcrete3D` tension | `lch > k·E·Gf/ft²` (`k` from `tensionLaw`) | `k·E·Gf/ft²` | `lch·ft² / (k·E·Gf)` |

- `regularize()` gains a return value, or an out-parameter: the ratio when
  it floors, 0 otherwise. The two kernels (`ASDConcrete*::HardeningLaw`,
  `ladruno_rc_kernel`) stay byte-identical in what they compute.
- `LadrunoConcrete3D` tension only **counts**. It does not add a clamp,
  because a clamp would change the response. Whether to add one is an open
  question (§6 Q2).
- The `epsFcFromGc` status `+1` case ("above the calibrated range") is not
  snap-back. Its existing warning moves to the same per-tag latch, but it
  is not counted here.

### D3 — per instance, latch then register; never per iteration

The cost has to stay negligible on the hot path, and the counts have to be
safe under the threaded update and commit. So the work splits in two:

1. **At the latch** (the existing once-per-instance `regularization_done` /
   `regularizationDone` / `gcLch` site). The instance stores, in its own
   members, `armedT/armedC` (bool), `ratioT/ratioC` (double) and the active
   element's tag. These are instance-local writes, safe under the threaded
   update. It then calls `ladrunoCrackBandNoteArmed(tag, side, ratio, ele,
   lMax, className)`. That call takes the registry's `std::mutex`, creates
   or updates the entry (`worstRatio`, `worstEle`, `lMax`), and emits the
   per-tag warning if it has not been emitted (D4). This happens once per
   instance in the instance's lifetime, so the lock is taken O(points)
   times per run. **It does not bump `nArmed`.**
2. **At `commitState()`**, one bool test per commit. The first commit of an
   armed instance bumps `nArmed`, under the same mutex. The first commit in
   which that side's committed damage is non-zero bumps `nRealized`:
   - ω_t / ω_c for `LadrunoConcrete3D`;
   - the tension/compression damage of the ASD state;
   - the backbone's post-peak branch for the RC kernel.
   Two instance flags (`reportedArmed`, `reportedRealized[side]`) make each
   bump happen once. Counting at commit follows WP-152: a trial counts
   calls, a commit counts states. Both bumps skip while
   `ladrunoProbeScopeDepth() > 0`, so LadrunoBrick/LadrunoQuad hourglass
   shadow copies, which commit under `LadrunoProbeCommitScope`, are never
   counted as points.

Why register at the latch and count at the commit:

- A shadow copy also latches. Counting at the latch would inflate
  `nArmed` by the shadow count.
- Registering at the latch still produces the warning and `worstRatio` when
  the run dies before its first commit. That is the case a brittle-drop user
  most needs. The entry then reads `nArmed = 0`, `worstRatio > 1`, and D5
  flags that state.

**Cost.**
- Steady state: one predictable branch per `commitState()` per point.
- Per run: one mutex acquisition per armed point at the latch, and at most
  three at commit.
- Each instance gains a few bytes of state (two doubles, an int, a few
  bools).
- No new work in `setTrialStrain`'s iteration path. The WP-200/107 profiler
  budget does not need re-measuring. S1 still checks a 10k-point brick deck
  for < 0.5 % wall-time change.

### D4 — one warning per material tag, per model

- One `opserr` per material tag, at the first armed registration of that
  tag, whichever side comes first:

  ```
  WARNING ASDConcrete3D (tag 3): crack-band snap-back clamp (tension): lch = 412.0 > l_max = 337.5
  (lch/l_max = 1.22) at element 1017. The softening energy is floored at the peak triangle: the
  response there is more brittle than Gf and mesh-dependent. Further clamps for this tag are
  counted, not printed; read them with `ladrunoCrackBandClamps` after analyze.
  ```

- It replaces `LadrunoConcrete3D`'s per-point compressive line (§1.1) and
  its "above the calibrated range" line. Those become per tag too.
- The latch lives in the registry entry. So it resets on `wipe` with
  everything else: the WP-107 rule that a warn-once latch belongs to the
  owning object, not the process. A `static bool` would mute every model
  after the first.
- Never per step, iteration or point.
- An opt-in `ladrunoCrackBandClamps -quiet` setting is not proposed. One
  line per tag is already the floor.

### D5 — the read channel: a fork query command, plus a per-point response

**Primary: `ladrunoCrackBandClamps`**, one row in `LadrunoCommandTable.h`,
so it registers in openseespy, interpreter-Tcl and classic `OpenSees.exe`
with nothing else to wire (WP-168). The registry test
`test_ladruno_command_registry.py` proves the three-interpreter presence.

```tcl
ladrunoCrackBandClamps <-tag $matTag> <-print>
```

- It returns a flat list, 7 numbers per entry, in `(matTag, side)` order:
  `matTag side nArmed nRealized worstRatio worstEle lMax`. `side` is 0 for
  tension and 1 for compression.
- An entry with `nArmed == 0` and `worstRatio > 1` means a point latched
  but never committed, usually because the run failed first.
- `-print` also writes a table to `opserr` that names each class.
- It is read-only. It returns an empty list when nothing was clamped.
- openseespy returns a Python list. apeGmsh wraps it into records.

**Secondary: a per-point material response `crackBandClamp`** →
`[armedT, armedC, ratioT, ratioC]`, on all five classes. It is reached by
`eleResponse` / `recorder Element … material $i crackBandClamp` (or the
section path for fibers and layers), and the `.ladruno` recorder can record
it like any other material response.

- This is what apeGmsh ADR 0129 D3 needs. D3 compares *which elements*
  apeGmsh flags with which the fork clamped, and a per-tag count cannot
  name elements.
- `worstEle` names only one element. The per-point response names all of
  them.

**`.ladruno` recorder.** At close (`LadrunoRecorder::~LadrunoRecorder`,
`LadrunoRecorder.cpp:439`), the recorder writes the registry as a table:
`INFO/DIAGNOSTICS/CRACK_BAND_CLAMPS`, an `N×7` dataset plus a `CLASS`
string column, beside the existing `INFO` group (`:985`). One file per
rank, as today. This is an addition to schema v1 that readers can ignore.
Slice S3 bumps the schema minor and updates `ladruno_schema_v1`.

**Alternatives rejected.**
- *Read through any point's response, the `implexRefusals` precedent.*
  - It needs an element and a point index to address a number that belongs
    to a tag.
  - It means five classes, both nD and uniaxial, and fiber or layered
    sections where the path is several tokens deep.
  - For SANISAND that was acceptable: one class, one solid element. Here it
    is not.
- *`getParameter` / `Parameter`.* That mechanism exists to update
  parameters for sensitivity, not to read state. It would need a parameter
  object per tag and a domain object to own it.
- *Recorder only.* Not every run attaches a recorder, and the
  reconciliation test should not depend on one.
- *Reuse `ladrunoCommitRefusalCounter`.* It is per-commit and cleared by
  `Domain::commit()`; its semantics are "fail this step", which a clamp
  must never do.

### D6 — lifetime: `wipe` is the zero point, as WP-104

- `ladrunoCrackBandCensusResetOnWipe()` is called from
  `OPS_clearAllNDMaterial()` (`NDMaterial.cpp:110`) **and** from
  `OPS_clearAllUniaxialMaterial()` (`UniaxialMaterial.cpp:79`), because
  ASDConcrete1D is uniaxial. The reset is idempotent, so a wipe that runs
  both hooks is fine.
- `ops.reset()` / `revertToStart()` do **not** clear the registry. Inside
  one model it is read as deltas (WP-104 rule).
- The instance flags are not reset by `revertToStart()` either. The
  regularization latch they shadow survives it too
  (`regularization_done`), so a reverted point is not counted again.
- `remove element` leaves the counts alone. They describe the run, not the
  live domain.

### D7 — copies, sections and MPI

- **Copies → tag.** Every Gauss-point copy has its parent's `getTag()`, so
  the registry key aggregates them without extra plumbing. That holds for
  ASD's per-tag `HardeningLawStorage` too.
  - Wrappers (`PlaneStress*`, `LogStrain`, `InitDefGrad`, `StagedStrain`,
    `LayeredShellFiberSection` layers) hold a copy of the inner material,
    which reports under the **inner** tag. That is the tag apeGmsh emitted
    for the concrete.
  - Fiber sections hold uniaxial copies, which aggregate the same way.
- **Threads.** D3's mutex covers the threaded update and commit.
  The instance members are written only by the thread that owns the
  instance.
- **openseesmp / OpenSeesMP.** Each rank is its own process and its own
  registry, and each element lives on one rank. So:
  - `ladrunoCrackBandClamps` on a rank reports **that rank**;
  - the global count is the sum of `nArmed` and `nRealized` across ranks,
    and `worstRatio` is the max;
  - the `.ladruno` recorder already writes one file per rank, and apeGmsh
    sums the per-rank tables.
  An MPI reduction (`-allRanks`) is §6 Q3.
- **OpenSeesSP.** The elements live in the worker ranks' subdomains, so P0's
  registry is empty. On a `PartitionedDomain` the query prints a warning
  that the counts live on the subdomain ranks, and returns what P0 holds.
  Aggregating SP is out of scope (§6 Q3).

## 3. Verification plan

Each case is one element with **one integration point**, so "exactly one" is
a statement about points, not about a Gauss rule. Each case runs to
softening on the side under test, in displacement control.

- **In-range (`lch / l_max ≈ 0.5`)**: `ladrunoCrackBandClamps` returns `[]`,
  and the per-point response is `[0,0,0,0]`.
- **Oversized (`lch / l_max ≈ 1.5`)**:
  - exactly one entry for the tag and side, `nArmed == 1`, `nRealized == 1`;
  - `worstRatio` equal to the analytic `lch / l_max` to 1e-12 (it is
    computed from the same doubles);
  - `worstEle` the element's tag, and `lMax` equal to D2's formula;
  - exactly one `WARNING … crack-band snap-back clamp` line, captured from
    a child process.

| Case | Host (1 point) | `lch` rule (apeGmsh D3 cross-check) |
|---|---|---|
| ASDConcrete3D T, C | `LadrunoBrick -formulation uri` (single point, stiffness hourglass), plus a twin with the C3a shadow-damage mode on. That twin builds 8 shadow copies (`LadrunoBrick.cpp:323-331`) and must still report exactly 1, which exercises the probe-scope guard | `cbrt(V)` |
| ASDConcrete1D T, C | `truss` (`getCharacteristicLength` = L) and a 1-IP `forceBeamColumn` with a single-fiber section | L; `wt·L` |
| `LadrunoConcrete3D` T, C | `LadrunoBrick -formulation uri` | `cbrt(V)` |
| `LadrunoRCConcrete` T, C | `LadrunoCST` plane stress (1 point) | `sqrt(2A)` |
| `LadrunoRCFiniteStrain` T, C | `LadrunoCST -geom finite` | `sqrt(2A)` |

Plus:

- **Armed, never softened.** An oversized element loaded only on the
  *other* side gives `nArmed == 1`, `nRealized == 0` on the untouched side.
- **The 1.01 band.** ASD at `lch = 2EGf/ft² · 0.995` is clamped by the
  fork. This pins §1.1 fact 2 so that apeGmsh can adopt the fork's `lMax`.
- **Aggregation.** A 2×2×2 brick mesh with an 8-point rule, all oversized,
  gives `nArmed == 64`. A mesh where half the elements are oversized gives
  half.
- **Wipe.** Model A clamps, then `wipe`, then model B with a different tag
  and in range: `[]`, and B's oversized twin warns again (the latch reset).
  `ops.reset()` between two legs of one model does not clear it. Classic-Tcl
  twin of the test, as WP-104.
- **Threads.** The aggregation deck under `LADRUNO_THREADS=4`, where both the
  threaded update and the threaded commit engage, gives the same counts as
  serial, repeated 20×.
- **Before any commit.** An oversized deck whose first step is refused
  still prints the one warning and reads `nArmed == 0`, `worstRatio > 1`.
- **Three interpreters.** The command-registry test covers the row.
  Classic `OpenSees.exe` runs the ASD oversized case.
- **Mutations.** Dropping the commit bump fails the count tests. Dropping
  the probe-scope guard fails the shadow-damage twin. Dropping the wipe
  call fails the wipe test. Making the latch a `static bool` fails the
  re-warn.

These one-element cases are exactly what apeGmsh ADR 0129 D3 replays. For
each family, apeGmsh's emit-time flag set must equal the set of elements
whose `crackBandClamp` response is armed, and apeGmsh's `l_max` must equal
the census `lMax`.

## 4. Slices

- **S1 — census core + ASD 3D/1D.** The header, the wipe hooks, the query
  row, the per-point response, the D4 warning, and the ASD cases from §3.
  This slice also records what `LadrunoConcrete3D` tension does past
  `l_max,t`, as a characterization only, and updates the stale
  "commit is serial" comment in `LadrunoMaterialStatus.h`.
- **S2 — the Ladruno concretes.** `LadrunoConcrete3D`: both sides, replacing
  the per-point warning. `LadrunoRCConcrete` / `LadrunoRCFiniteStrain`
  through the shared kernel's return value.
- **S3 — `.ladruno` recorder.** The `INFO/DIAGNOSTICS/CRACK_BAND_CLAMPS`
  table, a schema minor bump, and an apeGmsh reader stub.
- **S4 — hand-off.** A ledger fragment per slice (`implementations/`,
  `vanilla/` for `NDMaterial.cpp`, `UniaxialMaterial.cpp` and the two ASD
  files; a `quirks/` entry for the 1.01 band and the per-point
  `LadrunoConcrete3D` warning). Banner line: "crack-band snap-back census:
  `ladrunoCrackBandClamps`". apeGmsh then bumps FORK_PIN and lands its D3.

Vanilla footprint: `ASDConcrete3DMaterial.cpp`, `ASDConcrete1DMaterial.cpp`
(already fork-touched, #2044), `NDMaterial.cpp`, `UniaxialMaterial.cpp`.
Each change is a call or a return value, marked `// Ladruno WP-218`.

## 5. Consequences

- A silent clamp becomes one line per tag and a number that can be read
  after `analyze`. apeGmsh gets a fact to grade its prediction against.
- `LadrunoConcrete3D` stops printing one compressive warning per Gauss
  point. Logs from runs that hit it get shorter. Any test that counts those
  lines must change (S2 greps for them).
- The census reports armed points, not damage. A clamped point that never
  cracks counts in `nArmed`, not `nRealized`. Users read `nRealized` as "it
  mattered".
- Responses stay the same. No clamp is added or removed in this ADR; the
  only change to an existing number is the return value of the two
  `regularize()` functions.

## 6. Open questions (the owner decides)

- **Q1 — what apeGmsh D3 compares against: `nArmed` or `nRealized`?** The
  proposal is armed, because the prediction is a geometric statement ("would
  clamp if it softened"). `nRealized` is what a user is told mattered.
- **Q2 — `LadrunoConcrete3D` tension: count only, or also clamp?** Today it
  has no floor. The proposal is to count only in this ADR. A floor at
  `l_max,t` (the ASD rule) would change responses and needs its own
  measured WP.
- **Q3 — MPI aggregation.** Is per-rank, summed by the reader, enough? Or
  should the query take `-allRanks` (an MPI allreduce on openseesmp)? And is
  OpenSeesSP in scope at all? The proposal is per-rank, with SP out of
  scope.
- **Q4 — should the 1.01 factor count as part of `l_max`?** The proposal is
  yes: the census reports the fork's actual threshold, and apeGmsh adopts
  `lMax` rather than its own `2EGf/ft²`. apeGmsh would also need to mirror
  the ASD compression side and the `LadrunoConcrete3D` table floor `Gc/g[0]`,
  or read them from the binary.
- **Q5 — beams (apeGmsh ADR 0129 D4).** The per-point response gives each
  integration point's state, so the fork already separates end and interior
  points. Should the census also split `nArmed` into end and interior points
  for beams? The proposal is no; the response is enough.
