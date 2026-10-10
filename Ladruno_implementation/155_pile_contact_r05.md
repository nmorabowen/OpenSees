---
title: "WP-155 — Pile-contact R0.5: mortar augmentation control, pairing guard, initial-gap shift"
project: Ladruno
type: ADR (amends the mortar lane of ADR-41; status row in the ADR-48 capstone)
status: "Merged (#897, d63f49750, 2026-09-30) — three opt-in -mortar flags, defaults byte-identical; review #897 (Fable) findings 1-5 addressed"
owner: nmora
related:
  - "[[48_ladruno_contact_capstone_adr]] (status-of-record row + command surface)"
  - "[[41_ladruno_mortar_alm_contact_adr]] (the mortar/ALM lane this amends)"
  - "[[39_ladruno_contact_domain_adr]] (NTS; why the flags are mortar-only)"
  - "[[78_ladruno_parallel_contact_adr]] (the definitions stream: v3 -> v4)"
  - "apeGmsh piles-validation ladder/R1_plumbing/VERDICT.md (N-1, N-2, G-9)"
tags: [adr, contact, mortar, alm, pile, interference, wp-155]
updated: 2026-09-29
---

# WP-155 — Pile-contact R0.5

> [!summary] The short version
> Three opt-in options on a `-mortar` contact. With none of them the build is byte-identical: the
> contact battery gives 340/340 on both the base and the new binary, and `setMortarContactOptions`
> is never called.
>
> | Flag | Fixes | What it does |
> |---|---|---|
> | `-augment commit\|request\|never` | **N-2**: a linear tie depended on the number of load steps | chooses when the commit-cycle Uzawa update (λ_N, λ_T, λ_tie) runs: every commit (the default), only inside the `analyze_augmented` bracket, or never (pure penalty) |
> | `-maxGap d` | **N-1**: one tie around a closed cylinder bound antipodal facets | refuses a (slave facet, master facet) pair whose slave centroid is farther than `d` from the master facet's plane, at `handle()` time |
> | `-gapOffset g0`, `-adjust [tol]` | **G-9**: no interference fit or initial-gap control | shifts the normal gap by `g0` (negative = interference); `-adjust` subtracts the as-meshed gap so a faceted interface starts exactly closed and stress-free |
>
> On the R1 pile deck (skin radius = hole radius, no geometry hack),
> `-adjust -gapOffset -0.001 -augment never` reproduces the R1 geometric 1 mm interference-fit
> result to within 0.02 %. It is independent of the step count; see §5.

## 1. Context — what R1 found

The apeGmsh pile ladder, rung R1, models a rigid cylindrical pile skin (quad facets, RBE2-coupled
to a beam) inside the cylindrical hole of a soil block. It surfaced three mortar-lane defects. The
verdict is `piles-validation/ladder/R1_plumbing/VERDICT.md`:

- **N-1.** `LadrunoContactHandler` pairs facets by brute force: every master facet against every
  slave facet. The kernel clip then accepts an anti-parallel facet, because its test is `|cos|`.
  On a closed surface, the far side of the cylinder projects onto each facet, so a single `-tie`
  bound each skin facet to its antipodal soil facets as well. The pile was **3.65× too stiff**.
  R1 worked around this by splitting the skin into four sectors.
- **N-2.** `LadrunoContactDomain::commit()` runs one Uzawa update per commit, unconditionally.
  A *linear* tied problem therefore depends on the step count: 0.7395 mm in 1 step against
  0.7131 mm in 4 steps at ε = 1e6. The cohesion-only C lane failed at step 2.
- **G-9.** There is no interference or adjust option. The prestress had to be built into the
  geometry (skin radius R + δ). A faceted skin in a faceted hole also starts with spurious gaps
  and penetrations.

## 2. Decisions

### D1 — `-augment commit|request|never` (N-2)

The Uzawa update is gated per contact inside `LadrunoContactDomain::commit(bool augmenting)`. The
four multiplier updates (λ_tie, λ_T, λ_N, and the ADR-57 E6 edge-edge λ_N of a `-edgeAlm` contact)
are skipped when:

- the contact is `never`, or
- the contact is `request` and the call is not inside the `ladrunoBeginAugment` bracket.

`Domain::commit()` passes its existing `contactAugmenting` flag, set by `ladrunoBeginAugment`. The
rest of the commit always runs: slip promotion (`gpT`), engagement-origin double-buffers, and the
edge-edge sign and friction promotion. That is path state, not augmentation.

The edge-edge E6 Uzawa was missed at 79daec57f and gated after review #897 (finding 1). The key
carries the contact tag, so it uses the same lookup. Regression test:
`test_augment_never_gates_the_edge_edge_alm`, which is the adr57 E7 real-brick ALM model. A mutant
that restores the ungated update fails it (`never`: the bracket status is 0 instead of 1).

- **Why three values and not two.** The request asked for "`never` = pure penalty, with
  `analyzeAugmented` still usable on request". Those are two different contracts:
  - `request` is pure penalty on physical steps and ALM inside the bracket. This is the
    capstone's contract #3, with the per-commit update removed.
  - `never` is pure penalty everywhere. The bracket is inert, so a driver cannot re-introduce
    augmentation by accident.
  - Keeping both costs one comparison.
- **Why the default stays `commit`.** Byte-identity. A deck that relied on "augments for free
  across load steps" keeps that behaviour.
- **Committed-only invariant (contract #2) is preserved.** The multipliers still change only in
  `commit()`. Skipping an update leaves them at their last committed value; with `never` that is 0.

### D2 — `-maxGap d`, a plane-distance guard at pairing time (N-1)

At `handle()`, in the 3D mortar loop, a pair is skipped when
`|n_m · (c_s − c_m)| > d`. Here `c_s` and `c_m` are the slave and master facet centroids and `n_m`
is the master facet's Newell unit normal, all in the **current** (committed) configuration.

Four alternatives were rejected:

- **A normal-opposition check.** Antipodal facets on a cylinder are exactly parallel, so `|cos|`
  cannot see them. A *signed* test needs a consistent facet winding, which apeGmsh does not
  guarantee: G-1/N-4 exist precisely because the winding and outward are unreliable. Proximity is
  winding-free.
- **A per-GP gap cutoff at runtime.** A tie gauss point would then switch on and off between
  iterates, which breaks contract #5 (the connectivity superset is frozen within an epoch) and
  makes the residual discontinuous. Pairing is an epoch decision, so the guard belongs where
  pairing happens.
- **A centroid-distance test.** A facet's legitimate partners in a coarse/fine pairing can lie a
  facet-length away tangentially. The plane distance measures only the normal separation, which
  is the quantity that is O(2R) for the antipodal case and O(h²/R) for a true neighbour.
- **Applying it to NTS.** NTS candidates come from the bucket sort, with 27-neighbour cells sized
  by the facet diagonal, so it is already proximity-bounded. The defect is specific to mortar's
  brute force.

The window is wide. On the test cylinder (R = 0.5, 16 vs 24 facets), proto T4 measures the
largest near-pair plane distance at 4.3e-3 and the smallest antipodal-overlap distance at 0.98. Any
`d` between them keeps exactly the true neighbours, so a value of a few element sizes is safe.

### D3 — `-gapOffset g0` and `-adjust [tol]` (G-9)

Both act in one place: `LadrunoContactFE::mortarActive()`, which every consumer reads (residual,
tangent active set, friction cone, viscous mask, SOFT=2, the λ_N accumulator, and
`ladrunoMortarPenetration`). The weighted gap is shifted in normalized form:

```
g~_I  <-  ( g~_I / a_I  +  s_I ) * a_I ,     s_I = gapOffset + adj_I ,   a_I = sum_J D_IJ
adj_I = -gbar_I(reference)   if (tol == 0 or |gbar_I(reference)| <= tol), else 0
```

- **The shift is on the nodal weighted gap, not per gauss point.** It is the quantity the pressure
  `p_I = min(0, λ_I + ε·ḡ_I)` reads. It needs no new state: the adapter computes `adj_I` lazily
  from the **reference coordinates** of its own two facets. That result is deterministic, so an
  adapter rebuilt at the next `handle()` recomputes the same bits. Contract #1 (adapter = stateless
  view) is kept.
- **Exact zero at the start.** At u = 0, `g/a + (−g0/a0)` subtracts a double from itself, so
  `ḡ_shift` is exactly 0.0 and not a round-off residue. Proto T3 checks this bitwise on 1e5 random
  samples, and the FE test measures `umax == 0.0` exactly.
- **Closed branch at p = 0.** A gap-shifted node that sits exactly at `p = 0` (every `-adjust` node
  at the start) takes the closed side of the kink in the *tangent*. The residual is 0 on both
  sides. The first Newton iterate therefore sees the interface stiffness instead of a floating
  pile. Unshifted contacts never reach this branch.
- **Why the reference and not "the first handle".** Capturing at the first handle needs Domain
  state keyed per (pair, node) that survives adapter rebuilds, and it changes meaning if a deck
  re-handles after a displaced stage. The reference is what Abaqus `ADJUST` and LS-DYNA
  `IGNORE` act on: the as-meshed geometry.
  - **Limitation, documented:** for a contact declared *after* a displaced stage, `-adjust`
    corrects to the undeformed mesh. Use `-gapOffset`, or declare the contact before the stage.
  - **Limitation, documented (review #897 finding 3):** the reference gap is captured **per facet
    pair**, and only for a pair whose facets **overlap at the reference configuration**. A pair
    that first overlaps later is unadjusted: `integratePair` at the reference returns empty, so
    `adj_I = 0`. This covers a pair formed by sliding and paired at a later `handle()`, or one
    that the reference clip misses. Its unadjusted as-meshed gap then enters the node's pressure
    next to the adjusted contribution of its neighbouring pair, so on a faceted curve a jump of
    order the chord sagitta can appear when such a pair engages. `-adjust` is therefore a
    small-sliding tool. For finite sliding, prefer the true geometry plus `-gapOffset`. A
    per-node reference that survives re-pairing is a possible follow-up.
- **The optional `tol`.** Only nodes with `|ḡ_ref| ≤ tol` are adjusted, as with Abaqus
  `ADJUST=value`. This keeps a genuinely open part of an interface open. `test_g9_adjust_tolerance_branch`
  exercises both sides. On the faceted slice every load-bearing node has |ḡ_ref| in
  (6e-5, 7e-5]. A tol of 5e-5 gives the raw response bit for bit, and 8e-5 gives exactly zero.
- **Refused on `-tie`.** A tie bonds the relative *displacement* `r = D u_s − M u_m`, never the
  gap. The as-meshed gap is already strain-free on a tie, so a shift has nothing to act on.
- **3D only.** The 2D mortar lane evaluates its own interval kernel (`mortarActive2D`). A 2D pair
  with `-maxGap`, `-gapOffset` or `-adjust` draws a named FATAL at `handle()`; it is never silently
  ignored. `-augment` works in 2D too, because it acts on the shared `MortarNormalState`.

### D4 — Serialization (ADR-78 P2)

The mortar record grows by five tail slots (augmentMode, maxGap, gapOffset, adjust, adjustTol),
from 38 to 43. `LCD_FMT_VERSION` is bumped from 3 to 4 per the protocol, so a v3 stream draws the
named version refusal.

The unpack re-applies the options through the same `setMortarContactOptions` choke point. The
harness `Ladruno_files/testbed/contact_p2/db_roundtrip_all_lanes.py` now declares the block on pair
D and adds four sensitivity cases. Its result:

- `resave_bitexact` true;
- `tips_ref == tips_rt`;
- all ten sensitivity cases true;
- `corrupt_restores_raised` 2.

## 3. Command surface

These are mortar options only. Each is refused without `-mortar`, with a named message.

```tcl
# Tcl (OpenSees.exe) -- the pile skin (master) against the soil hole (slave), per sector
contact 1 1 2 -mortar -epsN 1e7 -mu 1.0 -epsT 1e6 -outward 0.707 0.707 0 \
        -adjust -gapOffset -0.001 -augment never
# one tie around the whole closed skin, antipodal pairs refused
contact 5 9 10 -mortar -tie -epsTie 1e8 -maxGap 0.05 -augment never
# ALM only on request (pure penalty on physical steps):
contact 6 1 2 -mortar -epsN 1e7 -augment request      ;# then analyze_augmented(...)
```

```python
# openseespy (same tokens)
ops.contact(1, 1, 2, "-mortar", "-epsN", 1e7, "-mu", 1.0, "-epsT", 1e6,
            "-outward", 0.707, 0.707, 0.0, "-adjust", "-gapOffset", -0.001, "-augment", "never")
ops.contact(5, 9, 10, "-mortar", "-tie", "-epsTie", 1e8, "-maxGap", 0.05)
ops.contact(7, 1, 2, "-mortar", "-epsN", 1e7, "-adjust", 1e-3)   # adjust only |gap| <= 1 mm
```

| Option | Values | Default | Refused when |
|---|---|---|---|
| `-augment` | `commit` / `request` / `never` | `commit` | any other word; no `-mortar` |
| `-maxGap d` | `d > 0` | off | `d ≤ 0`; no `-mortar`; a 2D pair (at `handle()`) |
| `-gapOffset g0` | any real; `< 0` is interference | 0 | `-tie`; no `-mortar`; a 2D pair |
| `-adjust [tol]` | flag, optional `tol > 0` | off | `tol ≤ 0`; `-tie`; no `-mortar`; a 2D pair |

## 4. Gates

`tests/test_adr155_pile_contact_r05.py` has 27 cases. The oracle
`contact_prototypes/proto_adr155_r05.py` passes 11/11 checks. Run with the worktree's
`dist\bin` first on `sys.path`, as the conftest's `_testbed.ops` requires.

| Gate | Result |
|---|---|
| **Byte-identity** | Contact battery (58 files: `tests/*{contact,mortar,adr39,adr41,adr57,adr85,adr96,ladrunoTie}*`): base binary (fd87e396d) 340 passed, new binary 340 passed |
| **N-2**, linear tie, `never`, 1/2/5 steps | u_tip = 2.2e-3 for 1, 2 and 5 steps, bit-identical (max relative difference 0.0). The default gives 2.2e-3 / 2.1e-3 / 2.04e-3, exactly the oracle's Uzawa recursion (T2) |
| **N-2**, tie limit | relative error 1e-1, 1e-2, 1e-3, 1e-4 at ε = 1e5 … 1e8 (O(1/ε), ratio 10.000) |
| **N-2**, modes on the interference blocks | `request`: penalty σ = 9.900990099 on physical steps, then 10.0000000 (to 1e-10) after 5 held-load augmentations. `never`: the bracket is inert. `commit`: augments every step |
| **N-2**, cohesive bond | μ = 0, c = 100 under a prestress, shear ramp in 1 and 5 steps with `never`: same u to 1e-10 |
| **N-1**, closed-cylinder tie | tube-in-ring slice (16 / 24 facets): `-maxGap 0.1` single tie = 4-sector split to **−4.2e-9** relative. The unguarded single tie is 3.3 % stiffer on this slice (3.65× on the R1 pile, where the antipodal bond also locks rotation) |
| **G-9**, `-adjust` | the same faceted slice as a frictionless contact, zero load: raw umax = 3.70e-5 (spurious penetrations push), `-adjust` umax = **0.0**, penetration 0.0, 1 iteration |
| **G-9**, `-gapOffset` | two blocks, δ = 1e-3: σ = δ/(2L/E + 1/ε) = 9.900990099 (FE matches to 3e-14 relative); after augmentation δE/2L = 10 |
| **G-9**, usable under load | `-adjust` slice under lateral load converges; `-adjust -gapOffset −1e-4` gives a symmetric shrink fit (net ux < 1e-12, every tube node moves radially inward) |
| Serialization | `db_roundtrip_all_lanes.py`: bit-exact re-save, identical tips, 10/10 sensitivity |
| Command surface | bad values, non-mortar use and gap shift on `-tie` are all refused; every documented spelling parses |

## 5. The R1 pile deck, rerun

The deck is the R1 C deck regenerated at δ = 0: skin radius = hole radius, 4 sector contacts,
matched hole, μ = 1, ε_N = 1e7, ε_T = 1e6. It was run through **`OpenSees.exe` (the Tcl lane)** of
this build, with the flags appended to the four `contact` lines and no other edit. The flags are
`FIT = -adjust -gapOffset -0.001`. u_head is in mm; "net" is u_head minus the H = 0 prestress-only
run (−0.000845 mm).

| Variant (H = 100 kN) | Iterations per step | u_head | net |
|---|---|---|---|
| as meshed, no flags (no prestress; spurious facet gaps) | 21 | 1.1617 | — (+67 % over the B tie) |
| `FIT -augment never`, 1 step | 20 | 0.72030 | **0.72114** |
| `FIT -augment never`, 2 steps | 16, 20 | 0.72377 | 0.72461 |
| `FIT -augment never`, 5 steps | 16, 17, 15, 16, 18 | 0.72373 | 0.72457 |
| `FIT -augment never`, prestress step, then 1 step | 15 → 20 | 0.72365 | 0.72449 |
| `FIT -augment never`, prestress step, then 5 steps | 15 → 16, 15, 15, 16, 18 | 0.72385 | 0.72469 |
| `FIT` with the default `commit`, 5 steps | 16, 16, 27, 21, 25 | 0.72206 | — |
| R1, geometric fit (skin R + 1 mm), default, 1 step | 19 | 0.72037 | **0.7210** |
| R1, geometric fit, default, 2 / 5 steps (the N-2 drift) | 17, 22 / 16–25 | 0.7292 / 0.7347 | |

What the table shows:

- **The interference fit no longer needs the geometry hack.** In 1 step, `FIT` gives a net 0.72114
  mm against the geometric fit's 0.7210 mm, a 0.02 % difference.
- **The step-count drift is gone.** With `never`, 2 against 5 steps differ by 0.005 %, and 1 against
  5 steps with a prestress step differ by 0.03 %. The geometric fit under the default drifted
  +2 % (0.7204 → 0.7347 mm).
- **The remaining 1-step versus 2-step difference (0.48 %) is friction path dependence, not the
  Uzawa.** The gap offset is not ramped, so in one step the full prestress and the full head load
  build up together. With a prestress step first, 1 and 5 steps agree to 0.03 %.
- **Caveat: every μ = 1 number above runs on a FENCED friction lane (review #897 finding 2).** The
  skin and the hole do not match (quad skin against a triangulated hole), and the per-slave-node frictional state (`gT0`/`engaged`/`gpTtrial`/`lambdaTtrial`) is written once per facet pair, last-writer-wins (`LadrunoContactFE::addMortarFriction`, LEDGER_quirks MAJOR-1, fenced to matched meshes). `-augment never` removes only the **λ_T** half of that state. `gpT`,
  `gT0` and `engaged` are still per-node last-writer. The 0.02 % agreement with the geometric fit
  and the step-independence therefore show that the new flags reproduce the old recipe on the same
  lane. They do **not** validate the frictional mortar path on a non-matching mesh. No new test
  exercises the flags with μ > 0 on a non-matching mesh; the friction fix is its own slice (§6).

**The cohesion lane (μ = 0, c = 1000, ε_T = 1e7) is NOT fixed by this slice.** The evidence
contradicts R1's attribution of it to N-2:

- **Byte-identity holds on the real deck.** The R1 geometric deck `C_coh_d3_force` rerun through
  the new `OpenSees.exe` gives u_head = 0.70324876434071340 mm in 11 iterations, the R1 value to
  every printed digit.
- **Its 5-step run still fails at step 2, with `-augment never` too.** The default and `never` both
  fail after 8 and then 101 iterations, from the same step-1 state. The step-2 failure is therefore
  **not** the per-commit Uzawa.
- **The plausible cause** is the same fenced friction lane: the per-slave-node frictional state (`gT0`/`engaged`/`gpTtrial`/`lambdaTtrial`) is written once per facet pair, last-writer-wins (`LadrunoContactFE::addMortarFriction`, LEDGER_quirks MAJOR-1, fenced to matched meshes). A cohesive
  stick turns the inconsistent `gTeff` into a force immediately. It needs its own slice (see §6).
- **On the δ = 0 geometry the cohesion lane is worse.** The H = 0 prestress converges in 9
  iterations. The 1-step H = 100 run fails whether it uses `FIT`, `-gapOffset` alone, or a prestress
  step first. With `never` in 5 steps it reaches 80 kN, then fails in step 5. The default `commit`
  diverges in step 2.

## 6. Descoped, and why

- **Radial or per-facet orientation for one contact around a closed skin (G-1/N-4).** `-maxGap`
  makes a single *tie* correct. A single frictional *contact* still needs a per-facet allowed side,
  and a global `-outward` cannot describe a cylinder. That is an orientation feature
  (`-outward axis …`), not a pairing one, so the sector split stays the recipe for contact. It is
  a candidate follow-up, and apeGmsh's pile helper can emit the split.
- **`-maxGap` on NTS and on the edge-edge enumeration.** NTS is bucket-bounded (D2). The edge-edge
  lane has its own proximity band (`-edgeBand`).
- **NTS gap shift.** NTS never converged on the pile problem in R1 (N-3), so the value there is
  unproven. The H2 zero-gap fail-safe would also need a separate decision.
- **2D mortar gap shift and guard.** They are refused by name, not implemented.
- **The per-node friction-state defect (LEDGER_quirks MAJOR-1), deliberately NOT fixed here.**
  In 3D `addMortarFriction` (`LadrunoContactFE.cpp`, around lines 2170–2194), `gT0`/`engaged` are
  captured by whichever facet evaluates first, and `gpTtrial`/`lambdaTtrial` are overwritten by
  every later facet. It explains the cohesion-only multi-step failure plausibly (§5 shows the
  failure survives `-augment never`), and it underlies every μ > 0 number in §5. Recommended as
  the next fork slice: a per-(node, feTag) reconciliation or an area-weighted blend, with a
  shared-node friction regression and its own oracle.

## 7. Risks

- **`request` / `never` removes the across-step "free" ALM.** Penetration becomes O(1/ε) again
  unless the user calls `analyze_augmented`. This is the intended trade (step-independence over
  free accuracy), and it is documented at the flag.
- **`-adjust` hides a genuinely bad mesh fit.** It zeroes whatever gap exists. The `tol` form
  bounds it, and `ladrunoMortarPenetration` still reports the shifted gap.

## Implementation log

- 2026-09-29 — implemented. Files:
  - `SRC/domain/contact/LadrunoContactDomain.{h,cpp}`: `MortarContact` fields,
    `setMortarContactOptions`, `commit(bool)`, pack/unpack v4.
  - `SRC/domain/domain/Domain.cpp`: passes `contactAugmenting`.
  - `SRC/analysis/handler/LadrunoContactFE.{h,cpp}`: `setMortarGapShift`, the shift in
    `mortarActive`, the closed-branch tangent.
  - `SRC/analysis/handler/LadrunoContactHandler.cpp`: `-maxGap` guard, arming the gap shift, the 2D
    refusal.
  - `SRC/interpreter/OpenSeesOutputCommands.cpp`: parser.
  - Tests, oracle and harness as listed above.
- 2026-09-30 — the R1 pile deck rerun through `OpenSees.exe` (§5). Scripts are outside the repo
  (session scratch): the δ = 0 deck was regenerated with the R1 `r1_run.py`, and the flags were
  appended to its `contact` lines.
