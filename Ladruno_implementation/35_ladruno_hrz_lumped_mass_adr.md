---
title: Mass-conserving (HRZ) lumped mass for explicit dynamics
project: Ladruno
status: ready-to-implement
priority: high
owner: nmora
tags:
  - implementation
  - adr
  - integrator
  - explicit
  - mass
---

> [!warning] The `status:` above is STALE — this ADR has shipped
> Its frontmatter still carries the pre-implementation value. Trust
> [[LEDGER_implementations]] for *does it work / which PR*, and this ADR for *why*.
> Flagged 2026-08-23 by a ledger audit; see [[README]] §Conventions. Remove this
> banner when `status:` is corrected.

# Mass-conserving (HRZ) lumped mass for explicit dynamics

> **Design / ADR (pre-implementation).** Closes the follow-up explicitly flagged in
> both [[04_explicit_dynamics_and_energy_balance]] (Future Work §2, Risks) and
> [[05_robust_central_difference]] (C6 / inherited dt_cr caveats), and anticipated by
> the comment at `CriticalTimeStep.h:41` ("a true mass-conserving HRZ lump is a
> follow-up"). Adds the standard **Hinton–Rock–Zienkiewicz** diagonal lump as a third
> `CTSLumping` mode and a shared header-only utility, replacing the two existing
> lumps that are each defective on one axis. **Recommended for
> [[36_ladruno_selective_mass_scaling_adr]]** — mass scaling *functions* on the existing
> `Diagonal` lump (its per-node shares are already strictly positive), but HRZ is
> **required for robustness across the full element zoo** (row-sum can yield zero/negative
> shares on serendipity/beam-shell elements) and improves SMS fidelity by conserving total
> mass. Not a hard blocker, but it should land first. No new class, no new classTag — a
> lumping mode plus a util.

## What

A header-only `Ladruno::hrzLump(const Matrix& Mconsistent, const ID& dofTypes, ...)`
utility that takes a **consistent** element mass matrix and returns the
Hinton–Rock–Zienkiewicz diagonal lump: per spatial translational direction, scale the
consistent diagonal entries so their sum equals the element's total mass in that
direction. This is the lump that is simultaneously (a) **positive-definite**
(consistent diagonal entries are strictly positive), (b) **mass-conserving**
(Σ Mᵢᵢ = m_total per direction), and (c) **rotation-aware** (rotational DOFs get a
scaled, non-zero inertia rather than the zero that row-sum produces).

It is wired in two places:

1. **`CriticalTimeStep`** — a new `CTSLumping::HRZ` branch in the lumping block, so the
   per-element `K v = λ M v` pencil that drives the critical-time-step estimate uses an
   honest diagonal mass. Exposed as `-lump hrz` on the explicit integrators
   (`CentralDifferenceLadruno`, `ExplicitBathe`, `ExplicitBatheLNVD`).
2. **(optional, Phase 2)** elements that currently row-sum-lump in `formInertiaTerms`
   (e.g. `Brick` with `-lumped`) can call the same util so their `getMass()` itself
   returns a conserving lump.

**In scope:** the HRZ util; the `CTSLumping::HRZ` mode + `-lump hrz` parsing + the
2-bit `sendSelf/recvSelf` widening; element opt-in for `Brick`. **Out of scope:**
mass scaling (lives in [[36_ladruno_selective_mass_scaling_adr]]); changing any
element's *default* lump (defaults stay as-is for byte-compatibility); consistent
(non-diagonal) mass scaling.

## Why

The two lumps in the tree today are each wrong on one axis, and the code comments say
so:

- **`RowSum`** (`CriticalTimeStep.cpp:234-239`) — `Mᵢ = Σⱼ Mᵢⱼ`. Conserves total mass,
  but on higher-order / serendipity elements (8-node quad with midside nodes, Tet10,
  Bezier) shape functions go negative → **zero or negative diagonal masses**, and on
  beams/shells it produces **zero rotational inertia**. Both poison `M⁻¹` and the
  eigensolve.
- **`Diagonal`** ("diagonal-of-consistent", `CriticalTimeStep.cpp:230-231`) — `Mᵢ = Mᵢᵢ`.
  Always positive and dimensionally sane for rotations, but **does not conserve mass**
  (Σ Mᵢᵢ ≠ m_total), so it systematically mis-states inertia and therefore `dt_cr`.

HRZ is the textbook resolution (Hinton, Rock & Zienkiewicz 1976; Cook–Malkus–Plesha
§11.10; Felippa IFEM Ch. 31). It is also the *only* lump under which "distribute the
element's mass to its nodes" — the operation [[36_ladruno_selective_mass_scaling_adr]]
needs — is well-posed for the general element zoo.

## Where

- New code: `OpenSees/SRC/analysis/integrator/LadrunoMassLumping.h` — header-only,
  `namespace Ladruno`, sibling of `LadrunoFictitiousMass.h` / `LadrunoJ2Kernel.h`.
- Modify: `OpenSees/SRC/analysis/integrator/CriticalTimeStep.h` — widen
  `enum class CTSLumping { RowSum, Diagonal, HRZ }` (lines 34-42).
- Modify: `OpenSees/SRC/analysis/integrator/CriticalTimeStep.cpp` — add the `HRZ`
  branch in the lumping block (lines 227-240); needs each DOF's *type* (translation vs
  rotation) to group directions — obtain from the element's nodes / `getMass()`
  structure.
- Modify: `OpenSees/SRC/analysis/integrator/CentralDifferenceLadruno.cpp` (parser
  ~lines 116-123; `sendSelf/recvSelf` lumping encoding ~lines 559/581 — the bool field
  already exists, so this is a **3-way decode** of the existing `data(6)`, not a vector
  growth), `ExplicitBathe.cpp` / `ExplicitBatheLNVD.cpp` (same `-lump` surface, but these
  have **no** lumping field today — their `sendSelf` vector must **grow** to add one; see
  How). See finding B in the review log.
- Modify (Phase 2, optional): `OpenSees/SRC/element/brick/Brick.cpp`
  `formInertiaTerms` (lines 887-893) — route the `massType==1` branch through the util.
- Reference (pattern to copy): `LadrunoFictitiousMass.h` (header-only util consumed by
  integrators), `LadrunoHardening.h` (header-only shared numeric kernel).
- Build: header-only, no CMake/dependency change. New `.h` must be added to
  `Ladruno_scripts/stamp_headers.py` GLOBS and stamped (see
  [[feedback_always_stamp_header]]).

## How

### Formulation

Let **Mᶜ** be the consistent element mass (Voigt-assembled, n_dof × n_dof), and group
the DOFs by spatial translation direction *d* ∈ {x, y, z}. The element's total mass in
each translational direction is

- direct: m = ∫_Ωₑ ρ dΩ
- and equals the full sum of the consistent block in that direction: m = Σ_a Σ_b Mᶜ_{ab,d}
  (partition of unity, ΣN = 1).

**HRZ translational lump** (per direction *d*):

```
S_d   = Σ_a  Mᶜ_aa,d            (sum of consistent diagonal over nodes, direction d)
α_d   = m / S_d                  (the HRZ scale)
Mᴴᴿᶻ_aa,d = α_d · Mᶜ_aa,d        (lumped diagonal entry)
```

By construction Σ_a Mᴴᴿᶻ_aa,d = m (conserved) and Mᴴᴿᶻ_aa,d > 0 (since Mᶜ_aa,d > 0).

**Rotational DOFs** (beams/shells): α is computed *per translational direction* (α_x, α_y,
α_z), and on a beam these differ (e.g. `ElasticBeam3d` consistent mass: axial 140·m vs.
transverse 156·m ⇒ α_x ≠ α_y). A rotational DOF has no single owning direction, so "the α"
is ambiguous. **Convention (pinned):** scale each rotational DOF at node *a* by the
**mean of node a's translational α's**, `Mᴴᴿᶻ_rot,a = (⅓Σ_d α_d)·Mᶜ_rot,rot`. This gives a
positive, finite rotational inertia (vs. row-sum's zero) without claiming to conserve an
angular quantity the element doesn't track. If a future element makes even this mean
ill-defined, **fall back to `Diagonal` for its rotational DOFs and warn** rather than
forcing the convention.

**Reduction check:** row-sum gives `Mᴿᴬ_a = ρ∫N_a` while HRZ gives `Mᴴᴿᶻ_a ∝ ρ∫N_a²`;
these coincide *per node* only when the nodal contributions are equal — i.e. on
**geometrically regular** elements (square quad, parallelepiped hex, equal-length bar),
where all Mᶜ_aa,d are equal ⇒ both give m/n_nodes per node. So `-lump hrz` is byte-identical
to `-lump rowsum` on regular elements, and diverges on **distorted** elements (where the
two differ *per-node* but **both still conserve total mass**) and on higher-order/beam
elements (where row-sum is outright wrong: zero/negative or zero-rotational). HRZ is never
worse than row-sum, and strictly better wherever they differ.

### API

```cpp
// LadrunoMassLumping.h
namespace Ladruno {
  // Fill `mdiag` (length n) with the HRZ diagonal of the consistent mass `M`.
  // `dofDir[i]` = spatial direction (0/1/2) for a translational DOF, or -1 for a
  // rotational DOF (which then takes the mean of its node's translational α's).
  // Returns total conserved mass per direction for diagnostics.
  void hrzLump(const Matrix& M, const ID& dofDir, double* mdiag, int n);
}
```

**`dofDir` construction (pinned heuristic, not left open):** the caller builds `dofDir`
from the element's node layout — for each node `a = ele->getNodePtrs()[k]`, the first
`ndm` of its `getNumberDOF()` DOFs are translations (direction 0..ndm−1) and the
remainder are rotational (−1). For mixed-ndf or coupling elements (RBE2/RBE3) where this
is ill-defined, **fall back to `Diagonal` and warn** (per H-risk below).

`CriticalTimeStep.cpp` lumping block gains:

```cpp
case CTSLumping::HRZ:
    Ladruno::hrzLump(M, dofDir, mdiag, n);   // M already fetched + copied (D8 rule)
    break;
```

Integrator parser: `-lump rowsum|diagonal|hrz`. The serialization fix is **per-integrator**,
because the three integrators do not encode lumping the same way today:

- **`CentralDifferenceLadruno`** *already* serializes lumping as a bool in `data(6)`
  (`sendSelf` ~559, `recvSelf` ~581 decodes `(data(6)!=0)?Diagonal:RowSum`). Keep the
  `Vector(7)`; carry the int enum {RowSum=0, Diagonal=1, HRZ=2} in `data(6)` and make
  `recvSelf` a 3-way switch. There is no version slot, so the only unsafe pairing is an
  *old* binary decoding `HRZ=2` as `Diagonal` — document it.
- **`ExplicitBathe` / `ExplicitBatheLNVD`** serialize *no* lumping field today (`sendSelf`
  sends `Vector data(1)=p`). These must **grow** the vector to add a lumping element, and
  `recvSelf` must default-read a missing field to each integrator's current default
  (CDL default = `Diagonal`; Bathe default = `RowSum`). Bump each integrator's DB schema
  note in the ledger.

### Testing

Zone-A pytest (`tests/test_hrz_lumped_mass.py`, numpy oracle):

1. **Mass conservation** — for a distorted Tet10 / 8-node quad with midside nodes,
   assert Σ Mᴴᴿᶻ_ii (per direction) = ∫ρdV to ~1e-12; assert all entries > 0
   (row-sum fails this — include row-sum as a documented xfail/contrast).
2. **Rotational non-zero + direction-dependent α** — `ElasticBeam3d -cMass` consistent
   mass (axial 140·m vs transverse 156·m ⇒ α_x≠α_y): assert HRZ rotational diagonal is
   positive and equals the mean-α convention; row-sum gives zero (contrast).
3. **Reduction (regular only)** — **square** quad / **uniform** bar: `-lump hrz` ≡
   `-lump rowsum` bit-for-bit. Plus a **distorted bilinear quad** case: assert HRZ ≠
   row-sum per-node but Σ HRZ_ii = Σ rowsum_ii = ∫ρdV (both conserve).
4. **dt_cr effect** — `criticalTimeStep()` on a beam/shell model under `-lump rowsum`
   (rotational zero ⇒ inflated/garbage dt) vs `-lump diagonal` vs `-lump hrz`; assert
   HRZ gives the physically-correct, mass-conserving value and is stable.
5. **g++ standalone** — call `hrzLump` on a hand-built consistent matrix, compare to
   the numpy oracle (header-only, OpenSees-free — same gate style as `LadrunoJ2Kernel`).

## Decisions

| # | Decision | Rationale | Consequence / extension point |
|---|----------|-----------|-------------------------------|
| H1 | Add HRZ as a **third `CTSLumping` mode + shared util**, not a new class | It is a lumping *scheme*, not an integrator; the enum + comment at `CriticalTimeStep.h:41` already reserve the slot | Zero classTag/broker churn; all three explicit integrators get it via `-lump hrz` |
| H2 | Rotational DOFs scaled by the **mean of the node's translational α's** (α is direction-dependent, so no single α exists) | Gives positive finite rotational inertia (row-sum gives zero); well-defined even when α_x≠α_y; `Diagonal`+warn fallback if even the mean is ill-defined | Documented convention; not an angular-momentum-conserving claim |
| H3 | `-lump hrz` must **reduce to row-sum** on *geometrically regular* elements only (square quad, parallelepiped hex, equal-length bar) | Backward-compatibility gate; on distorted elements HRZ≠row-sum per-node but both conserve total mass | Safe to recommend as the new explicit default *later*, but not in this ADR |
| H4 | **Do not change any element's default lump** | Byte-compatibility with existing models; defaults are a separate, louder decision | Defaults stay; HRZ is opt-in via `-lump hrz` |
| H5 | Ship the util header-only in `namespace Ladruno` | Reuse by `CriticalTimeStep`, elements, and the mass-scaling integrator; g++-verifiable OpenSees-free | Future elements (Tet10, shells) call the same util in `getMass()` |

## Risks / open questions

> [!note]
> DOF-direction tagging (resolved, see How): `CriticalTimeStep` builds `dofDir` from
> `ele->getNodePtrs()[k]->getNumberDOF()` — first `ndm` DOFs per node are translations,
> the rest rotational. Mixed-ndf / coupling elements (RBE2/RBE3) where this is
> ill-defined fall back to `Diagonal` and warn. (Left as a `[!question]` in the first
> draft; pinned here after review finding `hrz-dofdir-api-presented-as-solved`.)

> [!question]
> Should HRZ eventually become the *default* `-lump` for the explicit integrators
> (replacing `diagonal`)? Defer to a follow-up once test 4 quantifies the dt_cr
> difference on representative beam/shell/solid models.

- **Numerical**: consistent mass diagonal `ρ∫N_a²` is strictly positive for properly
  integrated standard elements, so α is well-defined. Guard **both** `S_d > 0` *and*
  per-entry `Mᴴᴿᶻ_aa,d > 0` (an under-integrated/degenerate element can produce a
  non-positive consistent diagonal entry); fall back to `Diagonal` with a warning if
  either fails. Note the `Diagonal` fallback is itself non-conserving — it is a
  last-resort floor, not a silent substitute. The downstream eigensolve also self-guards
  (`CriticalTimeStep.cpp:89-91`).
- **Build/ABI**: header-only, no risk; remember the header stamp + `stamp_headers.py`
  GLOB ([[feedback_always_stamp_header]]).
- **Backwards compat**: new enum value + widened `sendSelf` field — old DBs/parallel
  streams predate HRZ; bump the schema note and default-read missing field as RowSum.
- **Ledger debt**: this PR updates `LEDGER_quirks.md` (the row-sum/diagonal entry at
  lines 249-269 gains the HRZ resolution) and, if Phase 2 lands, a
  `LEDGER_vanilla_files.md` row for `Brick.cpp` with a `// Ladruno` marker.

## Implementation log

- **2026-06-19 — adversarial review (45-agent workflow, 6 dimensions × refute-by-default verify).**
  Verdict: **sound-with-fixes**. Confirmed-correct (do not touch): the HRZ formula
  `Mᴴᴿᶻ_a = m·Mᶜ_aa/Σ Mᶜ_bb` (re-derived, positive + mass-conserving on distorted
  elements); the row-sum code reference; the CDL serialization seam; the two-ADR
  ordering. Folded in: (A) pinned the rotational-α convention (mean of node's
  translational α's; α is direction-dependent) — H2, Formulation, test 2; (B) corrected
  the serialization plan to be **per-integrator** (CDL = 3-way decode of the existing
  bool; ExplicitBathe/LNVD = *grow* the vector, no field today) — How, Where; (C)
  restricted the "reduces to row-sum" claim to **geometrically regular** elements + fixed
  test 3 — Formulation, H3, tests; pinned the `dofDir` heuristic; added a per-entry
  `Mᴴᴿᶻ>0` guard; softened the 35→36 dependency to "recommended, not hard". Refuted: the
  "Phase-2 Brick is broken" cluster (it's a no-op demonstrator on a regular hex, already
  marked optional).

*(filled in once executing; move to `Ladruno_internal/implemented_hrz_lumped_mass.md` when done)*
