---
title: EnergyBalanceRecorder v2 — explicit energy channels + parallel-safe reduction
project: Ladruno
status: draft
priority: medium
owner: nmora
tags:
  - implementation
  - recorder
  - explicit
  - energy-balance
  - parallel
---

# EnergyBalanceRecorder v2 — explicit energy channels + parallel-safe reduction

> **ADR** (Architecture Decision Record). Successor scoping to
> [[04_explicit_dynamics_and_energy_balance]] (the v1 `EnergyBalanceRecorder`,
> classTag 33000). v1 is a **closure** instrument — four accounted channels
> `KE IE DW ULW` + residual `RES` + `ERR%`. This ADR scopes turning it into an
> **attribution** instrument: named channels for the energy sinks/sources that
> today either silently pollute `RES` (closure gaps) or hide inside an existing
> channel (attribution gaps) — absorbing boundaries, hourglass, contact, modal /
> LNVD damping, mass scaling — and fixing the **cross-rank double-counting** on
> partitioned models. No code yet; this pins the schema and the two mechanisms so
> the sign conventions and the parallel rule are settled before implementation.

## What

Extend `EnergyBalanceRecorder` from six columns to a **channel-aware** balance in
which each significant energy sink/source is either (a) accounted in its own
column or (b) *provably* single-counted inside an existing one, so `RES` means
"genuinely unexplained energy" rather than "everything we didn't wire up."

Two orthogonal problems, two mechanisms:

1. **Attribution & closure of physical channels** — modal/LNVD damping, absorbing
   radiation vs. incident injection, hourglass, contact. Resolved by a
   generalized **energy-channel registry** (mechanism A), a **region-isolation**
   convention (mechanism B), and a small **`Element::getEnergy()`** query
   (mechanism C).
2. **Parallel reduction correctness** — nodal terms at shared partition-boundary
   nodes are multiply-counted when per-rank sidecars are summed. Resolved by a
   **count-where-owned** rule plus an **element/nodal split** in the output.

**Not in scope:** MPCO integration (v1 decision D1 stands — plain-text sidecar
only, [[project_mpco_exit_crash]] frozen); a global auto-reduced single MPI file
(deferred — see mechanism (4) below); contact-energy *internals* (defer to
`opensees-contact`; this ADR only defines the channel it plugs into).

## Why

Explicit dynamics has no equilibrium iteration to fail loudly. v1 catches a run
that *fails to close*. It does **not** catch the two failure modes that matter
most in practice on this fork's workloads:

- **"Closed but wrong."** A run can close `RES` to <1% while 40% of `IE` is
  hourglass ringing, or while mass-scaling has shifted the period. The energy
  balance is blind to this because the polluting energy *is* accounted — just in
  the wrong bucket. Only a dedicated column exposes it. (v1's own header warns:
  *"a 'closed' RES computed without an existing sink's channel is falsely
  reassuring."*)
- **SSI / soil models lie today — and may lie while *closing*.** Verified from
  source: a Lysmer compliant-base boundary's incident-injection force leaks into
  `IE` through a stale member (`getResistingForceIncInertia` fills
  `internalForces = C·v_gnd`; stage-0 `getResistingForce` returns that same
  member, and recorders run at commit, after residual formation). With
  resisting-force sign the leak is ≈ `−W_inject`, so `RES` can *accidentally
  close* while `IE` goes strongly negative — the worst failure mode, because the
  tripwire stays green. ASD's injection (`addBaseActions`/`addRffToSoil`) is in
  `getResistingForce` by design — same class. Either way the injection needs a
  first-class **source** channel and an IE exclusion. Given the fork's
  SSI/DRM/absorbing-boundary direction ([[project_absorbing_pml_guide]]), this is
  the driving use case.
- **Partitioned totals are silently wrong on reduction.** v1's header names only
  `ULW` as double-counted on shared nodes; in fact nodal `KE` and nodal `DW` ride
  the identical shared-node path. Any offline sum of per-rank files is suspect.

## Where

- Modify: `SRC/recorder/EnergyBalanceRecorder.{h,cpp}` — column schema, output
  split, ownership gate, `OPS_` parser flags.
- Modify: `SRC/recorder/EnergyBalanceKernel.h` (`namespace ebkernel`) — the
  single source of truth for the per-entity math; add channel accumulation.
- Generalize: `SRC/analysis/integrator/LadrunoMassScalingEnergy.h`
  (`Ladruno::MassScalingEnergyRegistry`) → a channel-keyed
  `Ladruno::EnergyChannelRegistry`. This is the **proven template** — the
  consistent-SMS KE correction already uses exactly this publish/consume pattern
  ([[project_ladruno_mass_scaling]] ADR-38 / V4).
- Producers that publish into the registry:
  - `SRC/analysis/integrator/ExplicitBatheLNVD.cpp` — LNVD/FLAC local damping
    work (its `formUnbalance` override `r_i -= α|r_i|sign(v_i)`; fires per
    sub-step — see the cumulative-publish requirement in mechanism A).
  - `SRC/analysis/integrator/TransientIntegrator.cpp` — modal damping
    (`modalDampingFactors` lives here; this is the application site to
    instrument). **Vanilla file** → `// Ladruno` marker + ledger row.
  - `SRC/element/absorbentBoundaries/LysmerTriangle.cpp`,
    `ASDAbsorbingBoundary{2,3}D.cpp` — incident-injection work. **Vanilla
    files** → `// Ladruno` markers + ledger rows, strictly-additive hooks per
    [[../CLAUDE|the fork change-budget rule]].
- Mechanism C (element energy split — hourglass/contact):
  - **Preferred, zero-vanilla:** the existing `setResponse/getResponse` path —
    the recorder builds `Response` objects once in `buildCache()` (e.g.
    `setResponse({"energy","hourglass"})`) and queries per record. No new
    virtual, no core-header edit.
  - **Alternative, cleaner API:** a new `Element::getEnergy(int type)` virtual —
    but that touches `SRC/element/Element.h`, a **core vanilla header**
    (tree-wide recompile, ledger row, upstream-divergence cost). Note upstream
    already has an energy idiom at the material level
    (`UniaxialMaterial::getEnergy()`, consumed by `DispBeamColumn`/
    `ForceBeamColumn`) — if the virtual route is taken, align with that naming,
    don't invent a parallel one.
  - Overriders/responders: `SRC/element/.../LadrunoBrick*` (hourglass / EAS
    stabilization energy — fork-owned, free to change); contact elements
    (deferred to `opensees-contact`).
- Reference (the pattern to copy): `LadrunoMassScalingEnergy.h`
  `MassScalingEnergyRegistry::elementScalingKE` + its consumer in
  `ebkernel::addElementEnergy`.
- **Ledger impact:** every vanilla producer file above → `LEDGER_vanilla_files.md`
  row in the same PR; new fork-authored files → `LEDGER_implementations.md` +
  LADRUNO header stamp.

## How

### The governing balance and the two gap types

Physical work–energy over the run:

$$ W_{\text{ext}} + W_{\text{inj}} \;=\; E_{\text{kin}} + E_{\text{int}} + E_{\text{damp}}
   + E_{\text{hg}} + E_{\text{cont}} + E_{\text{rad}} + \ldots $$

Classify every term by *how it fails today*, because the fix differs:

- **Closure gap** — reaches no channel → pollutes `RES`. Balance truly doesn't
  close. Fix = give it a channel (mechanism A or C).
- **Attribution gap** — is in a channel but mislabeled → `RES` closes while a
  channel lies. Fix = a *separate* column so the fraction is visible; `RES` is
  unaffected.

### Accounting map (verified against current source)

| Sink / source | Lands today in | Gap | Fix mechanism |
|---|---|---|---|
| Lumped mass scaling | `KE` (`Node`/`Element::getMass`) | none | — (done) |
| Consistent SMS (Olovsson) | `KE` via `MassScalingEnergyRegistry` | none | registry (**shipped**, ADR-38) |
| MS *fidelity* (fictitious inertia) | `KE`, correctly | not energy — period/wave-speed error | advisory fraction column (P3) |
| Absorbing **radiation** (Lysmer/ASD dashpot) | `DW` (`getDamp` · v_node) | attribution (mixed w/ Rayleigh) | region-isolate (B) or channel (A) |
| Absorbing **injection** (Lysmer `gnd_velocity`) | `IE` — stale `internalForces = C·v_gnd` returned by stage-0 `getResistingForce` | attribution + **sign** (leak ≈ `∓W_inj`; RES may accidentally close w/ IE lying) | registry-publish `E_inject` + **IE exclusion** (A) |
| Absorbing **injection** (ASD `addBaseActions`/`addRffToSoil`) | `IE` (folded into `getResistingForce`) | attribution + sign | registry-publish `E_inject` + **IE exclusion** (A) |
| Hourglass / EAS stabilization | `IE` (`getResistingForce`) | attribution | `Element::getEnergy` (C) |
| Contact penalty (recoverable) | `IE` | attribution | `Element::getEnergy` (C) |
| Contact friction (dissipated) | `IE` | attribution | `Element::getEnergy` (C) |
| Contact chatter / penalty jump | **`RES`** | closure | registry (A), hard |
| Modal damping | **`RES`** (integrator solve, no `getDamp`) | closure | integrator-publish (A) |
| LNVD / FLAC local damping | **`RES`** (`formUnbalance` override) | closure | integrator-publish (A) |
| Numerical dissipation (Noh–Bathe, gen-α) | **`RES`** | closure (intentional) | integrator-publish (A), optional |

### The decisive absorbing-boundary finding

`LysmerTriangle::getDamp()` returns the real dashpot `C = diag(ρV_s, ρV_s, ρV_p)`
(transformed) — so the recorder's `DW = ∫ vᵀC v dt`, dotted with **node**
velocities, **does** capture radiated energy (the node-velocity damping force
reaches the residual via the assembled global C, so `DW` is the physically
radiated power). The injection path is subtler — traced through the stale-member
mechanics (`LysmerTriangle.cpp:479,505`):

1. Each step the integrator calls `getResistingForceIncInertia()`, which fills the
   member `internalForces = C·(0·v_node + gnd_velocity) = C·v_gnd` — the incident
   forcing only.
2. Stage-0 `getResistingForce()` zeroes `springForces` but **returns that same
   `internalForces` member** unmodified.
3. Recorders run at commit, *after* residual formation — so the kernel's
   `IE_rate = F·v` picks up `(C·v_gnd)ᵀ v_node`: the **injection work, with
   resisting-force sign**, inside `IE`.

Net:

- **Radiation OUT** → captured in `DW` ✓ (but indistinguishable from structural
  Rayleigh → attribution gap, fixed by region-isolation).
- **Incident IN** → leaks into `IE` with resisting-force sign. If the leak is
  `−W_inj` (resisting convention), then `IE_rec = IE_true − W_inj` and
  `RES = W_inj − (KE + IE_true + DW) ≈ 0`: the balance **accidentally closes
  while `IE` goes strongly negative** — a silent lie, worse than an open
  residual. If the loader's sign convention flips it, `RES ≈ −2·W_inj` and the
  balance visibly fails. **The sign must be pinned by the Ricker test, not
  assumed.** Observable signature of the leak either way: `IE` trending negative
  in an elastic compliant-base run.
- **Radiation-only Lysmer** (no `LysmerVelocityLoader`) has `v_gnd = 0` → no
  leak → correct today with no changes.

ASD's injection (`addBaseActions`, `addRffToSoil`) is inside `getResistingForce`
by design — the same attribution+sign class. Both need the injection published as
a first-class **source** term, *and the corresponding leak excluded from `IE`*
(see the double-counting requirement in mechanism A).

### Mechanism A — generalize the registry

Promote `MassScalingEnergyRegistry` (single-purpose `elementScalingKE`) to a
channel-keyed singleton any producer publishes into and the kernel consumes:

```cpp
namespace Ladruno {
enum class EnergyChannel { MassScaleKE, ModalDamp, LNVD, AbsorbInject,
                           NumericalDiss, /* ... */ };

class EnergyChannelRegistry {
public:
  static EnergyChannelRegistry& instance();
  bool active(EnergyChannel c) const;         // no producer => zero-cost skip
  // producer side (integrator/element): ACCUMULATE energy at the producer's own
  // cadence (per sub-step, per iteration — whatever the producer's natural unit
  // of work is). The registry stores a cumulative running total per channel.
  void addEnergy(EnergyChannel c, double dE);
  // consumer side (ebkernel): read the running total. NO integration in the
  // recorder for these channels.
  double energy(EnergyChannel c) const;
  void   reset();                             // on domain wipe / analysis reset
};
}
```

**Publish cumulative energy, not power.** The obvious alternative — producers
publish a per-step *power* that the recorder trapezoid-integrates — is wrong for
two structural reasons:

1. **Recorder cadence ≠ step cadence.** Recording every N steps trapezoids a
   sampled rate → aliasing error in exactly the channels we are adding to chase
   sub-1% closure. (This aliasing already exists in v1's `IE/DW/ULW` when
   recording sparsely — documented below as a v1 caveat; don't replicate it.)
2. **Sub-stepped producers have no single per-step rate.** LNVD's
   `formUnbalance` override fires in *both* Noh–Bathe sub-steps with different
   forces; the only well-defined quantity is the work increment per firing.

The SMS registry works precisely because it publishes an *instantaneous*
quantity (`M̄_e`, consumed as KE at the recorder's own sample time). Work-type
channels get the equivalent: the producer integrates at its own cadence, the
recorder reads a running total. `active(c)==false` ⇒ no column, no cost — the
lumped path and every vanilla integrator stay byte-identical, as the SMS registry
already guarantees.

**No double-counting (hard requirement).** For the absorbing-injection channel,
publishing `E_inject` is *not sufficient*: the injection force is **already
inside `IE`** via `getResistingForce` (Lysmer stale-member, ASD by design).
Lighting up the channel without removing the leak counts the energy twice and
*re-opens* the residual. P1 must pair every published channel with its
exclusion: either the kernel skips the `F·v` contribution of these element
classes (classTag test in `addElementEnergy` — cheap, but hard-codes a list), or
the element exposes the injection component separately so the kernel can subtract
it (cleaner, touches vanilla files). Decide per-element in P1; the Ricker closure
test is the gate either way.

**Sign convention (open question Q1, leaning source-side):** injection is a
*source*, entered on the `W_ext` side so `RES` stays interpretable:

$$ \text{RES} = (\text{ULW} + E_{\text{inject}}) - (\text{KE} + \text{IE} + \text{DW}
   + E_{\text{modal}} + E_{\text{lnvd}}) $$

### Mechanism B — region isolation (works today, zero code)

For any spatially-localizable sink, put the elements in a `MeshRegion` and read
the region block. **Absorbing radiation works now**: one region for all
Lysmer/ASD elements ⇒ that region's `DW` = radiated energy, cleanly separated from
structural Rayleigh. Ship as a documented modeling convention (P0). Does **not**
solve injection closure (that's A).

### Mechanism C — element energy split (hourglass / contact)

For attribution gaps *inside* `IE` where only the element knows the split
(hourglass vs. real strain energy; contact penalty vs. friction). Two routes:

- **Preferred (zero vanilla footprint):** the `setResponse/getResponse` path.
  The recorder resolves `Response*` objects once per participating element in
  `buildCache()` and queries them per record. Fork elements (`LadrunoBrick`,
  contact) add an `"energy"` response — fork-owned files only. Cost: one virtual
  dispatch + Vector copy per element per record; acceptable at recorder cadence.
- **Fallback (only if the response path proves too awkward):** a
  `virtual double getEnergy(int type)` on `Element` — core vanilla header edit,
  tree-wide recompile. If taken, align with the existing upstream
  `UniaxialMaterial::getEnergy()` idiom rather than inventing a sibling.

Either way the kernel emits `E_HG`, `E_CONT` columns and (for elements that also
fold the split component into `getResistingForce`) subtracts it from `IE` — same
no-double-counting rule as mechanism A. This is how LS-DYNA reports hourglass
energy — a first-class channel, never a residual.

### Parallel-safe reduction

Root cause (verified): `Node` carries **no owner/rank field**, so at the flat
per-rank `Domain` level a recorder cannot tell a shared boundary node from an
interior one. Therefore:

- **Element terms** (`KE/IE/DW` from `ele->…`) — elements are uniquely
  partitioned ⇒ **single-counted** ⇒ safe under naive cross-rank summation.
- **Nodal terms** (`KE/DW/ULW` from `node->…`) at shared nodes — the `Node` object
  exists on every rank touching it ⇒ **multiply-counted** on reduction. v1 names
  only `ULW`; nodal `KE` and `DW` are equally exposed.

Two-part fix:

1. **Element/nodal output split (P1, cheap).** Emit
   `KE_ele IE DW_ele | KE_nod DW_nod ULW` so the element block reduces by naive
   summation (always correct) and the nodal block is quarantined for careful
   handling. Turns a silent error into a visible, separable one.
2. **Count-where-owned gate (P2, correct form).** Skip a node's nodal terms unless
   this rank owns it. Clean in the Subdomain model (skip `getExternalNodes()`
   ghosts); in the flat MP model needs the parallel shared-node map plumbed to the
   recorder (`Node` has no owner query today). Same plumbing a true in-recorder
   `MPI_Allreduce` would need ⇒ pinned here, not bolted on.

### Proposed column schema

The v2 layout is **opt-in behind a flag** (`-v2` or `-split`); the default
remains the legacy six columns, byte-identical, so existing parsers and the
regression battery are untouched. (These two goals — an always-on element/nodal
split and byte-identical legacy output — are mutually exclusive; the flag
resolves it in favor of compatibility.)

With `-v2`, prefer a **flexible, stable header** the registry/regions populate by
active channels (avoids a wall of zeros; stable names keep parsers simple):

```
time | KE_ele KE_nod IE DW_ele DW_nod ULW
     [ E_inject E_modal E_lnvd E_hg E_cont ]   # only columns for active producers
     | RES ERR%
     [ per-region blocks ... ]
```

`RES`/`ERR%` recomputed against the *full* active set so closure improves as
channels light up. P1 channels are **global-only** (no per-region channel
attribution — that needs per-element granularity in the registry and is
deferred; regions still split the four legacy channels as today).

### Testing

- **Closure regression (must stay green):** existing `_verify_explicit.py` /
  `_smoke_energy_bathe.py` — without `-v2` and with all producers inactive the
  six legacy columns must be byte-identical (lumped path, vanilla integrators).
- **Sign-pinning FIRST (P0.5, before any channel code):** 1-D elastic soil column
  on a Lysmer compliant base, Ricker pulse, *current* recorder. Record `IE` and
  `RES`. Predicted signatures: either `IE → strongly negative` with `RES ≈ 0`
  (leak = `−W_inj`, accidental closure) or `RES ≈ −2·W_inj` (opposite loader
  sign). This single run resolves Q1's sign convention empirically and is the
  baseline the P1 fix is measured against.
- **Absorbing injection closes:** same model; after `E_inject` publish + IE
  exclusion, assert `|ERR%| < 1%` **and** `IE ≥ 0` throughout (the leak's
  signature gone). Reference: input energy = analytic pulse energy.
- **Double-count guard:** with `E_inject` active, assert closure does not
  *degrade* relative to the P0.5 baseline in the accidental-closure sign case —
  the regression that catches publishing without excluding.
- **Region radiation:** same model, all boundary elements in one region; assert
  region `DW` ≈ radiated energy (energy leaving a control volume).
- **Hourglass channel:** under-integrated `LadrunoBrick -uri` in an hourglassing
  mode; assert `E_hg` grows while total `RES` stays closed (attribution, not
  closure).
- **Parallel:** 2-partition model, sum per-rank files; assert element block sums
  exactly and nodal block flags/handles shared nodes; compare to serial.

## Risks / open questions

> [!question] Q1 — injection sign/ownership
> `E_inject` as a **source** on the `W_ext` side (recommended, keeps `RES`
> interpretable) vs. its own signed term. The *empirical* sign of the current IE
> leak is settled by the P0.5 sign-pinning run (see Testing) **before** any
> channel code; the remaining decision is only the reporting convention.

> [!question] Q2 — first-phase target
> SSI/absorbing closure (⇒ do A first for injection, defer C) vs. general explicit
> hygiene incl. hourglass (⇒ C matters equally)? Fork direction says SSI ⇒ A first.

> [!resolved] Q3 — schema: **resolved** — v2 layout opt-in behind `-v2`, legacy
> six columns stay the byte-identical default; with `-v2`, flexible active
> columns under a stable header. (An always-on split and byte-identical legacy
> output were mutually exclusive; compatibility wins.)

> [!question] Q4 — parallel scope for v2
> Ship only the element/nodal split (P1) and document the nodal caveat, or commit
> to the count-where-owned gate (P2) — which needs `Node` owner plumbing — in the
> same PR?

- **Backwards compat:** legacy six-column consumers must keep working. Achieved
  by the `-v2` flag (layout) plus the registry `active()` short-circuit
  (channels); guard with the byte-identical regression above.
- **Pre-existing v1 caveat, inherited:** `IE/DW/ULW` are trapezoid-integrated at
  *recorder* cadence — recording every N steps aliases the work integrals, and a
  spurious `RES` drift can be a cadence artifact, not physics. Chasing sub-1%
  closure requires record-every-step (document) or integrator-side accumulation
  (future). New channels avoid this by design (cumulative publish).
- **Producer discipline:** producers accumulate into the registry; the registry
  needs a single well-defined `reset()` on domain wipe / new analysis (same
  failure class as the SMS registry — reuse its reset hook). No per-step clear:
  channels are running totals by design.
- **Numerical dissipation is intentional** — if `E_numerical` is published, label
  it a scheme term, not a leak; a dissipative Noh–Bathe run *should* show it.
- `HDF5_USE_FILE_LOCKING=FALSE` still required for any HDF5-backed reading
  ([[project_mpco_ladruno]] quirk) — unchanged.

## Implementation log

*(filled in once execution begins; move to `Ladruno_internal/implemented_<name>.md`
when done)*

- **Phasing:** P0 region-isolation convention + docs (no code) →
  **P0.5 sign-pinning Ricker run** (current recorder, no code — pins Q1's
  empirical sign and baselines the leak) → P1 registry generalization +
  `E_inject`/`E_modal`/`E_lnvd` **with paired IE exclusions** + `-v2`
  element/nodal split → P2 element energy split via `getResponse` (hourglass
  first) + count-where-owned gate → P3 mass-scaling fidelity advisory column.
- **2026-07-05 self-review (pre-implementation):** corrected the Lysmer
  injection classification (leaks into `IE` via stale `internalForces`, NOT a
  pure RES gap — may accidentally close with `IE` lying); added the
  no-double-counting requirement (publish + exclude); switched registry API from
  per-step power to cumulative energy (recorder-cadence aliasing + LNVD
  sub-step cadence); made the v2 layout opt-in (`-v2`); surfaced the vanilla-file
  ledger footprint and the zero-vanilla `getResponse` route for mechanism C.
- **2026-07-05 P1 BUILT** (this worktree, pending battery + PR):
  - `SRC/analysis/integrator/LadrunoEnergyChannels.h` — cumulative-energy
    channel registry (`ABSORB_LEAK`, `LNVD_WORK`, `MODAL_WORK` reserved);
    `declare()`/`addEnergy()`; recorder baselines at first record so totals are
    monotone with no reset hook needed.
  - `EnergyBalanceKernel.h` — `EnergyAccumulatorV2` (ele/nodal split + channel
    math); V1 accumulator frozen. `EnergyBalanceRecorder` — `-v2` flag; legacy
    accumulation variables untouched (byte-identical guarantee is structural,
    not just tested); columns fixed at initialize from declared channels; ASD
    warnOnce.
  - Producers: `LysmerTriangle::commitState` publishes the leak (R_inj
    RECOMPUTED via `getDamp()·v_gnd`, not read from the stage-3-mutating
    `internalForces`; trapezoid per commit; dt≤0 commits resync only) —
    **the IE exclusion is the global `IE_display = IE_raw − L` rebucket in the
    accumulator, not a kernel classTag skip** (uniform for future ASD publish).
    `ExplicitBathe -lnvd` accumulates per-sub-step work
    (`p·dt` / `(1−p)·dt` attribution) into a scratch zeroed in `newStep`
    (failed steps never publish), published in `commit()` before
    `commitDomain()`.
  - **Scope cuts, deliberate:** ASD2D/3D publish deferred to **P1.5** — their
    injection is interleaved with four genuine mechanisms across ~2300-line
    files; a wrong split corrupts IE, worse than the documented gap (recorder
    warns). Modal publisher deferred to P2 (`TransientIntegrator.cpp`
    `modalDampingFactors` is the site; enum slot reserved).
  - Known accepted limitation: `-v2` channel columns are fixed at first
    record — an integrator constructed *after* the first record would miss its
    column (echoed layout makes it visible; not a plausible script order).
  - Known accepted wart: channel *declarations* are process-sticky (monotone,
    like the totals — no reset on `wipe`, by design: the recorder baselines
    deltas so totals never corrupt a later model). Consequence: a `-v2` run in
    a process that previously declared a channel shows that channel as a benign
    all-zero column. Zone-A tests order accordingly
    (`test_energyBalanceRecorder.py` note).
- **2026-07-06 P0.5 RAN + P1 VERIFIED (`Ladruno_implementation/energy_v2/`):**
  - **Q1 sign pinned empirically: outcome (a).** Leak = `−W_inject`
    (surrogate run: `IE_end/W = −0.9999`, `ULW = 0` exactly) — the
    accidental-closure mode. So published `L < 0` and `E_inject = −L > 0`
    reads as positive energy-in, as designed. (`p05_lysmer_sign_pin.py`,
    `RESULTS.md`.)
  - **P0.5 headline blocker: `LysmerVelocityLoader` was UNREACHABLE** — no
    `eleLoad` hook existed in any interpreter (and Tcl's `eleLoad` silently
    returns TCL_OK for unknown `-type` flags → 15 years of silent no-op).
    Fixed in-scope: additive `-lysmerVelocityLoader <dir>` branch in
    `TclModelBuilder.cpp` (ledgered). `LysmerTriangle` remains classic-Tcl-only
    (Python element registration → P1.5 with the ASD publish).
  - **Two new upstream quirks ledgered:** F2 — stage-0 Lysmer under implicit
    Newmark realizes only ~half the dashpot energy (C·v never enters the
    element residual; explicit integrators are consistent); F3 — the silent
    `eleLoad -type` accept. Consequence for gating: RES cannot close under
    Newmark+Lysmer for *solver* reasons, and the E_inject rebucket is
    RES-invariant by design → the P1 gate tests **IE truthfulness**, not RES.
  - **P1 Ricker gate ALL PASS** (`p1_ricker_closure.py`, worktree Tcl exe):
    G1 loader live; G2 legacy lie reproduced (`IE_end = −39.114 J` =
    `−E_ref`); G3 v2 IE truthful (`+1.2e−4`); G4 publisher↔recorder
    cross-check `E_inject = 39.114` vs `IE_v2 − IE_legacy = 39.114` (<0.01%);
    G5 radiation-only control clean (`E_inject = 0`, `IE ≥ 0`).
  - **Zone-A battery:** `test_energyBalanceRecorder.py` 4/4 (incl. the LNVD
    closure gate: published `E_lnvd` == legacy RES drift, v2 drift shrinks
    ≥5×); `_verify_explicit.py` 20/22 — the two failures (`dt_cr = −1`)
    reproduce EXACTLY on the pre-change main-checkout pyd → pre-existing
    environment issue, not this change.
- **2026-07-06 adversarial review (Opus 4.8, pre-push): NO CRITICAL/HIGH.**
  Verified independently: legacy byte-identity PASS (accumulation order,
  sizing, XML tags, echo, region offset algebra); RES algebra PASS (L cancels
  exactly, `RES_v2 = RES_v1 − D_lnvd`); producer lifecycles PASS (dt≤0 /
  revert / double-commit / recvSelf reseed; `updateCount` guard prevents the
  non-Linear-algorithm sub-step-1 over-count from ever publishing). Findings:
  MED-1 process-sticky column layout → user-facing doc added to the recorder
  header ("parse by column names, never by position"); LOW-1 pytest
  definition-order dependency (accepted, noted in the test file); LOW-2
  static-stage pseudo-time leak integration (moot — loader inactive at
  gravity, rate = 0); NIT-1 channels declared after first record correct IE/
  RES but emit no column (correct-but-quiet, normal flow unaffected).
- **2026-07-06 P1.5 BUILT + VERIFIED (same worktree, PR #488 merged first;
  synced onto `ladruno` before starting):**
  - **`LysmerTriangle` openseespy dispatch** (`OpenSeesElementCommands.cpp`)
    + the Python-side `eleLoad -lysmerVelocityLoader <dir>`
    (`OpenSeesPatternCommands.cpp`) — the element's parser already used
    interpreter-agnostic `OPS_Get*Input`, so only registration was needed;
    P0.5's Tcl-only gap is now closed in both languages.
  - **ASDAbsorbingBoundary2D/3D bottom compliant-base leak publisher.**
    Traced `addBaseActions()` end-to-end: a pure external source (time-series
    velocity × fixed coefficients, always additive into `R`, **no**
    stateful-member hazard unlike Lysmer's stage-3 mutation) — and confirmed
    from source that `addRff`/`addRffToSoil` (the LATERAL free-field
    mechanism) both `return` early whenever `m_boundary & BND_BOTTOM`, so the
    bottom-injection and lateral-injection source terms are **mutually
    exclusive**: publishing the bottom leak cannot double-count the lateral
    one. Scoped P1.5 to the bottom mechanism only — the lateral free-field
    injection remains an open, documented gap (P1.6+; recorder warning
    reworded to say so precisely rather than claim ASD is wholly
    uncovered). New `commitState()` override on both classes (neither had
    one) mirrors Lysmer's publisher exactly: recompute via `addBaseActions`
    into a scratch `Vector`, dot with `getVelocity()` (confirmed same index
    space via `m_dof_map`), trapezoid-integrate, publish to `ABSORB_LEAK`.
  - **ASD Ricker gate ALL PASS** (`p15_asd_ricker_closure.py`, openseespy —
    a 2D quad soil column on a ghost+real compliant-base pair, mirroring the
    element's own internal ghost/soil topology; ghost nodes additionally
    `fix()`-ed as belt-and-suspenders on top of the element's own internal
    penalty pin): G2 legacy lie reproduced (`IE_end = −82.3 J`); G3 v2 IE
    truthful (`+3.4e−5`); G4 cross-check `E_inject = 82.296` vs
    `IE_v2 − IE_legacy = 82.296` (6 sig figs); G5 zero-amplitude control
    clean.
  - **Zone-A gained a 5th test**
    (`test_energybalance_v2_asd_absorbing_injection_closure`) exercising the
    same gate as a proper pytest. **Caught its own version of the earlier
    column-location bug before it shipped**: with `LNVD_WORK` also declared
    by an earlier-sorting test file, `E_inject` is *not* the last channel
    column (`E_lnvd` always writes after it) — but since this test's element
    unconditionally declares `ABSORB_LEAK`, `E_inject` is reliably at the
    **fixed front-anchored offset** (col 7, right after `ULW`), which is what
    the fixed test asserts. General lesson reinforced: locate a channel by
    its own guaranteed-write-order position, not by an inherited assumption
    about what else is last.
  - **Full regression:** 5/5 (`_run_energy_tests.py` default set) and
    16/16 in the exact adr52-then-recorder ordering that caught the P1 bug;
    `_verify_explicit.py` 20/22 unchanged (same 2 pre-existing failures).
