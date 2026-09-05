---
title: "ADR 90 planning brief — Duvaut–Lions viscoplastic regularization wrapper (TIMs item 10)"
project: Ladruno
type: planning brief (pre-ADR)
status: "DRAFT 2026-09-04 — defines the planning phase; the ADR itself is not yet written"
owner: nmora
reserves: "ADR number 90 — NOT 88: PR #778 (H5DRM higher-order elements) cites ADR-88 throughout SRC/ and the ledgers, and 86_geomechanics_libraries_and_mpm_scoping_report §6.4 proposes ADR 89 for Track T. ND classTag 33022 (ND registry; the 33022 in the EigenSOE and PATTERN registries are per-registry, not collisions)"
related:
  - "[[86_ladruno_sanisand_adr]]"          # the material the wrapper serves first
  - "[[31_ladruno_concrete3d_adr]]"        # the shipped, material-specific Duvaut–Lions (-eta) and its gating template
  - "[[59_ladruno_gradient_concrete_adr]]" # the regularization ADR that was DESCOPED — the failure mode to avoid
  - "[[32_ladruno_dispbeamcolumn_regularization_adr]]"
  - "[[87_ladruno_depth_with_width_adr]]"  # warrant package a new material family owes
  - "[[testbed/00_canonical_testbed]]"
tags: [planning, adr, regularization, viscoplastic, duvaut-lions, sanisand, tims, localization]
updated: 2026-09-04
---

# ADR 90 planning brief — Duvaut–Lions viscoplastic regularization wrapper

> [!abstract] What this document is
> The **planning phase** for ADR 90, not the ADR. It records (1) what is already decided on both
> sides and must not be re-litigated, (2) what three source sweeps of the fork established on
> 2026-09-04, (3) the one contradiction the ADR has to resolve before anything else, (4) the
> acceptance case the fork owes TIMs a draft of, (5) the phases/gates, decisions and open
> questions the ADR will carry. Every claim below is cited to a file, a vault note, or a paper;
> nothing is asserted from memory. Numbers in `code` are `file:line` in this worktree.

External context: the TIMs / APE Gorini-macroelement calibration campaign (vault at
`C:\Users\nmb\Dropbox\UANDES EC\TIMs Workbench\vault\`, cited as **vault NN**). The fork's reply
of 2026-08-10 (vault 64 §6) offered a design view and asked TIMs to specify an acceptance case;
vault 65 D5 accepted the design and made the acceptance case its first work item. **As of
2026-09-04 that acceptance case has not been written on either side** (vault grep for
`imperfection|biaxial|acceptance case` hits only 64/65/73/Index; no fork commit, branch or
graveyard entry carries it). This brief exists to close that gap from the fork side.

---

## 1. Settled — do not re-open

| # | Settled item | Where |
|---|---|---|
| S1 | **Method:** a Duvaut–Lions viscoplastic **wrapper at the `NDMaterial` level**, not a modified `ManzariDafalias`. Wraps any inviscid material at the stress-update seam; serves DruckerPrager and future materials in the same stroke. | vault 64 §6 (fork's own words); vault 65 **D5** |
| S2 | **τ is a declared numerical parameter, not a soil property**; the deliverable is *"mesh-independent given a declared τ"*, characterised on a Deborah number. Uniqueness of the non-associated collapse load is **not** restored and is disclosed. | vault 10 line 75; vault 65 D3 + D5 |
| S3 | **Element frozen** `LadrunoBrick -formulation bbar`; **tetrahedra prohibited** on failure legs (0.286 of exact Prandtl–Reissner vs 1.084/0.994/0.951 for B-bar; turning B-bar off on the hex gives 0.400). | vault 65 D4; vault 71; fork R3 gate `tests/test_r3_prandtl_collapse_gate.py` (#722) |
| S4 | **Provenance is output, not documentation:** every regularized number ships with engine hash, element/formulation, τ, rate, three-mesh band, ultimate criterion. | vault 65 D3; fork `ops.ladrunoBuild()` idiom (#718) |
| S5 | **Characterisation protocol exists:** vary the relaxation parameter and the ramp duration independently; matched pairs must collapse onto the ratio (the vault did this for β_K/T to 0.0 / 0.7 / 4.2 %). | vault 14 |
| S6 | **Primary calibration instrument** is the displacement-controlled radial monotonic curve; probe procedure = radial probes, sequential swipe cross-check. | vault 65 D2, D7 |
| S7 | The fork's cost estimate stands: *the wrapper is the smaller half; the validation campaign proving a τ-controlled band width converging under refinement is the larger half.* | vault 64 §6 |
| S8 | The papers are in hand (Seafile `tim-macroelement-manuals\`): Perzyna 1966, de Borst–Sluys–Mühlhaus–Pamin 1993, Rudnicki & Rice 1975. **No paper acquisition phase is needed.** Targeted page reads happen while writing. | skill `tim-macroelement/references/library_map.md` |

---

## 2. What the 2026-09-04 source sweeps established

Three read-only sweeps (shipped `-eta` kernel; SANISAND attachment surface; the fork's ADR
acceptance bar). Findings that change the plan:

### F1 — The shipped Duvaut–Lions is material-specific. The wrapper is a **write, not a lift**.

`LadrunoConcrete3D -eta` is twelve lines inside the CDPM2 kernel
(`SRC/material/nD/LadrunoConcrete3DKernel.h:1492-1503`): a Simo–Hughes trial blend on the
**effective** stress, $\bar\sigma \leftarrow (1-\beta)\bar\sigma_{tr} + \beta\bar\sigma_{inv}$,
$\beta = \Delta t/(\eta+\Delta t)$, with the same blend on the hardening variable $\kappa_p$
(`:1497`) and on the effective tangent (`:1499-1501`) *before* damage is chained. It needs the
committed internal effective stress and $\kappa_p$ — nothing a generic `NDMaterial` seam exposes.
There is **no generic stress-update seam** in the fork today.

What *is* reusable wholesale is the **gating pattern** (ADR-31 §7 P3, `tests/_testbed/concrete3d_ref.py:3086-3210`):
PV1 τ=0 byte-exact; PV2 discrete→continuous order; **PV3 closed-form 1-D overstress**
$\sigma^*-\sigma_Y = E\dot\varepsilon\,\eta$ independent of Δt; PV4 blend identity to 1e-12;
PV5 tangent vs FD; PV6 overstress monotone in η; plus the C++ **non-tautology guard**
(viscous ≠ inviscid by > 1e-3, `tests/_testbed/concrete3d_kernel_check.cpp:386-390`).

### F2 — A generic wrapper is a *different* (two-track) model. That is a scoping decision, not a detail.

Textbook Duvaut–Lions relaxes the stress **and the internal variables** toward the inviscid
projection of the *viscoplastic* trial state. A wrapper that only sees `getStress()` cannot
re-seed the inner material's internal variables; it can only run the inner **inviscidly on the
total strain path** and let its own $\sigma_{vp}$ lag toward `inner.getStress()` ("two-track").
For perfect plasticity the two coincide (PV3 cannot tell them apart); for hardening / state-
dependent models (SANISAND: $\alpha$, fabric $z$, $e$, $\psi$-driven $M^b$) they do not.
ADR-31 refused the analogous shortcut for concrete because a nominal-stress blend relaxes
*damage* as well as plasticity (`31_ladruno_concrete3d_adr.md:291-292`). ⇒ **the ADR must
state which model it implements and scope the wrapper to plasticity-type inner materials**, and
P0 must measure whether two-track is an adequate regularizer for SANISAND (see §6 P0, OQ3).

### F3 — Δt in statics is controllable only under uniform pseudo-time.

`ops_Dt` = `Domain::dT` = current − committed pseudo-time (`Domain.cpp:2080-2082`, `:2125`,
`:2392`). Under `LoadControl(dλ)` it is a uniform, positive number and `-eta` is live — the
fork's own test relies on it (`tests/test_ladrunoConcrete3D_element.py:773-774`). Under
`DisplacementControl` the λ-increment is non-uniform (~1e5 swings, negative after `loadConst`)
— on the quirks ledger for IMPL-EX (`LEDGER_quirks.md:1440-1443`) and `-eta`
(`:1612`: `dt ≤ 0` ⇒ **inviscid**, deliberately not elastic). The ledger has **no** entry for
the converse risk — a rate artifact when τ is comparable to a static pseudo-step. ADR 90 adds
that quirk.

### F4 — The theory says quasi-static rate regularization is weak by construction.

de Borst et al. 1993 (Seafile `05_numerics_localization\deborst1993.pdf`, body pp. ~19–20):
rate dependence *"works equally well for both failure mechanisms … but its applicability is
obviously limited to transient loading conditions and the regularizing effect rapidly decreases
for slow loading rates or when we approach the rate-independent limit."* Their internal length
is $\ell = 2mc_e/E$ — a **wave-speed** quantity. Vault 65 D5's "$h \lesssim c\tau$" is that
dynamic statement. Under quasi-static pushovers there is no intrinsic length; the band width is
set by $\tau \times$ loading rate, i.e. by $\mathrm{De}$, and $\mathrm{De}\to 0$ recovers the
ill-posed problem. Also from the same pages: **the canonical acceptance deck** — plane-strain
biaxial specimen, smooth ends, symmetric half, an imperfect element (10 % cohesion reduction at
the left boundary near mid-height), Drucker–Prager, non-associated; their Fig. 18 ran it with
Duvaut–Lions.

### F5 — The SANISAND attachment surface, and its five traps.

- The base `ManzariDafalias` has **no** `setTrialStrain/getStress/getTangent`; only the
  dimensional wrappers do (`ManzariDafalias3D.cpp:78-149`). The strain-rate overload
  **discards** the rate (`:88-91`). The material reads no time at all.
- `getCopy(void)` is `exit(-1)`; only `getCopy("ThreeDimensional"|"3D"|"PlaneStrain"|"PlaneStrain2D")`
  exist (`ManzariDafalias.cpp:447-465`, `:534-555`). The exact bug `test_fspm_over_manzari_family.py` gates.
- `getStress()/getStrain()/getTangent()` return **static** class buffers (`ManzariDafalias3D.h:78-79`) — a wrapper must **copy, never alias**.
- `updateMaterialStage` resolves by string-tag match *inside the material* (`ManzariDafalias.cpp:791-828`); a wrapper must **forward `setParameter` to the inner on tag miss** (pattern: `FluidSolidPorousMaterial.cpp:320-333`, `InitStrainNDMaterial.cpp:436-473`) or the static `mElastFlag` never flips.
- `revertToLastCommit()` is an **empty stub** (`ManzariDafalias.cpp:513-517`) — the wrapper's viscoplastic history has no inner rollback to lean on (a pre-existing limitation, not a new one).
- The geotechnical sign flip happens inside the dimensional wrapper, so relaxation against `inner.getStress()` is already in element convention; ordering `{xx,yy,zz,xy,yz,zx}` matches `LadrunoBrick` (`LadrunoBrick.cpp:1869-1876`) — no permutation.
- **`LadrunoSANISANDPlaneStrain` (33021) exists**, order 3, so a 2D biaxial deck on the consumer's own material is viable — with two prerequisites: its `getCopy(void)` is uncovered (86 handoff §5) and `ManzariDafaliasPlaneStrain`'s null ctor sets the wrong classTag (`LEDGER_quirks.md:4826`).

### F6 — No width-objectivity gate exists anywhere in the fork.

Every objectivity gate the fork ships regularizes **dissipated energy** via `lch`
(`test_ladrunoRCConcrete_meshobj.py`, `test_ladrunoSolidShell_softening.py` G5,
`test_lemaitre_notched_bar.py`, `test_ladrunoBrick_asdconcrete*.py`); one of them says so in
its docstring: *"lch regularizes the dissipated ENERGY, not the localization WIDTH"*
(`tests/test_lemaitre_notched_bar.py:21-23`). Greps for `rudnicki|bifurcat|shear band|band width`
in `tests/` are empty. The **template** to copy is the G5 pairing: a positive gate at three mesh
densities **plus** a discriminating negative control (`test_g5_fixed_lch_not_objective`) that
proves the regularizer is what does the work.

### F7 — The acceptance bar, precisely.

- **ADR-59's kill list** (`59_ladruno_gradient_concrete_adr.md` §3, §8): an unverified "seam"
  claim (B1), an off-switch that is not byte-identical (B2), the wrong regularized variable (B3),
  "free" infrastructure reuse (M5), "mostly re-wiring" (M6), no named band-width consumer, no
  reconciliation with the fork's own shipped alternative (M1). Un-descope gate = named consumer +
  alternatives table + reconciliation. Every item has a direct analogue here (§9).
- **ADR-87 warrant package** for a new material family: verification manifest YAML
  (`Ladruno_implementation/ledger/<feature>.yml`, ≥1 non-self-referential oracle), a
  `-DLADRUNO_MUTATE_<FAMILY>` **mutation gate** that turns the suite red, a guide stub, both
  interpreters, Zone-A green, out-of-family verdict. Branch `wp/<n>-<slug>`, draft PR day one,
  merge commit, owner merges.
- **Docs-only scoping ADR** owes no manifest row (`ci/check_manifest.py` keys on `classTags.h`
  defines) — it owes a `LEDGER_implementations.md` row and a `README.md` index line.
- **Numbers:** next free ADR was **88** at sweep time (87 taken; 86 doubly used) — **superseded: 88 is cited by PR #778, 89 proposed for Track T ⇒ this ADR is 90**. Next free ND classTag =
  **33022** (`SRC/classTags.h:582-594`; other registries' 33022 are not collisions).
- **Slow tier / capacity rule** (`testbed/00_canonical_testbed.md` §1b–§1c): measured wall time
  in the docstring and ledger row; hard wall-clock stop reported as *inadmissible*; pinned
  controller budgets with calibration history; capacity = plateau + admissible termination +
  free advance.

### F8 — Stale documentation found in passing (fix in the ADR PR)

`LadrunoConcrete3D_guide.md:396,480,593` still lists `-eta` as "not yet (P3)"/"deferred";
`31_ladruno_concrete3d_adr.md:91` claims `-eta` works "under any tier" but the kernel gates it to
Tier-1 (`!mp.implex`, `LadrunoConcrete3DKernel.h:1492`) and it is inert in the BeamFiber view.

---

## 3. The problem the ADR must resolve first: **D5 contradicts D2/D7 under quasi-statics**

Vault 65 makes the displacement-controlled radial probe the primary instrument (D2, D7) and
Duvaut–Lions the regularizer (D5). On this engine those two are in tension three ways:

1. **F4** — the regularizing effect decays toward the rate-independent limit; a slow static push
   is that limit. The band width under Duvaut–Lions is a function of $\mathrm{De} = \tau/T$, not
   of the material.
2. **F3** — under `DisplacementControl` the pseudo-time increment that defines
   $\beta = \Delta t/(\tau+\Delta t)$ is erratic; De is not even well defined step to step.
3. **S2/S4** — every number must ship with τ *and rate*: the campaign already knows De is a
   coordinate of the answer, but has not confronted that De is set by the integrator, not the
   deck, on the static lane.

Candidate resolutions the ADR must weigh (an *alternatives table* is mandatory — F7):

| Option | What it means | Cost | Verdict to test in P0/P2 |
|---|---|---|---|
| **(a) Transient lane with τ replacing β_K** | Run the radial probes as slow damped transients — the reference model's *native* lane — with **Rayleigh β_K = 0** and the wrapper's τ as the one declared relaxation time. Δt is physical, De = τ/T is a deck quantity, and vault 14's machinery applies unchanged (its β_K artifact was *already* a Kelvin–Voigt relaxation of 23 ms — this replaces an accidental regularizer with a declared one). | Wall-clock (vault 14: ~17 min/run on the twin); P5's coupled lane still blocked by the contact `ndf` guard. | **Recommended primary.** |
| (b) Static lane under uniform `LoadControl` pseudo-time | Keep statics, but drive the push with a uniform pseudo-time so Δt = 1/nsteps is exact and De = τ (T = 1). Displacement targets via a displacement-controlled *load* series rather than the `DisplacementControl` integrator. | Loses arc-length-style path following past a limit point; vault 60's D16 stepping guards must be re-verified. | Acceptable **fallback**; must be gated identical-to-(a) at matched De on the P2 deck. |
| (c) Strain-driven internal clock in the wrapper | Replace Δt with an accumulated strain measure so the "rate" is path-length based, integrator-independent. | This is a different constitutive model (Perzyna-in-arclength), not Duvaut–Lions; no literature anchor; no closed-form gate. | **Rejected** unless (a)/(b) both fail. Listed for completeness. |
| (d) Do nothing to the lane; report the band | Run `DisplacementControl` as-is, report the three-mesh band without regularization (vault 10's disclosure stance). | Zero engine work; the residue stays as large as measured. | The **negative control** of the whole ADR, not an option. |

The ADR's headline claim therefore becomes: *"the regularized collapse load $q_u(\mathrm{De};h)$
converges in $h$ at fixed, declared De, and its De-dependence is measured and disclosed"* — not
"mesh-independent collapse load". This is consistent with D3 and with vault 10, and it is the
sentence that keeps ADR 90 out of ADR-59's failure mode.

---

## 4. Named consumer, and the claims the ADR may and may not make

**Named consumer (ADR-59 gate (a)):** vault 65 **P1** — the ultimate-surface loci fitted from
radial pushovers on SANISAND, whose post-peak/collapse values are mesh-sensitive by
ill-posedness. Band-**width** objectivity is required because the collapse load of a
non-associated softening strip depends on the mechanism's thickness, not only on the energy
it dissipates (F6 is exactly why the existing `lch` gates do not serve this consumer).

| Claim | Gate | Class |
|---|---|---|
| C1 τ = 0 is **byte-identical** to the inner material (off-switch) | PV1 clone at material point + element level; g++ byte check | correctness |
| C2 1-D steady overstress = $E\dot\varepsilon\tau$, Δt-independent | PV3 clone | correctness |
| C3 tangent = $(1-\beta)C_e + \beta C^{ep}_{inv}$, matches FD | PV4/PV5 clones | correctness |
| C4 **band width converges under h-refinement at fixed De** (positive) **and does not at τ = 0** (negative control) | new width gate, §5 | the deliverable |
| C5 collapse load at fixed De converges in h (three-mesh band collapses) | slow-tier gate, capacity rule | the deliverable |
| C6 De-dependence of width and $q_u$ is measured (τ × T sweep, matched pairs collapse) | vault-14-style sweep | disclosure |
| C7 the wrapper is transparent to `updateMaterialStage`, `getCopy(type)`, database round-trip, MP wire | Zone-A + roundtrip + `test_fspm_over_manzari_family` twin | wiring |

**Claims the ADR explicitly does NOT make:** restores uniqueness or the bound theorems (S2);
an intrinsic length in statics (F4); mesh-independence without a declared De; applicability to
damage models (F2); any change to the tetrahedron prohibition (S3).

---

## 5. The acceptance case — draft to hand TIMs as the OQ5 answer

The shape is fixed by de Borst 1993 Fig. 18 and vault 65 D5; the *numbers* are TIMs' to set
(OQ2). Proposed ladder, cheapest and most decisive first:

| Leg | Deck | Material | Element | Purpose |
|---|---|---|---|---|
| **A0** | 1-D softening bar, imperfect central element, N = 20/40/80/160 (de Borst Fig. 2 layout) | numpy oracle, two-track DL vs true DL | — | Settles F2/OQ3 before any C++: does two-track regularize width? |
| **A1** | Plane-strain biaxial, symmetric half, smooth ends, 10 % strength reduction in one boundary element at mid-height; 3 meshes h/h₂/h₄ + 2 orientations | `DruckerPragerPlaneStrain` (non-associated, cheap, no state dependence — clean bifurcation) | `quad PlaneStrain` | The canonical width gate; negative control τ = 0 |
| **A2** | Same deck | `LadrunoSANISANDPlaneStrain` (33021) | `quad PlaneStrain` | The consumer's material; needs the two 33021 prerequisites (F5) |
| **A3** | Same specimen extruded one element thick | `LadrunoSANISAND` (33019/20) | **`LadrunoBrick -formulation bbar`** | S3 compliance — the element the campaign actually uses |
| **A4** | R3 Prandtl–Reissner strip (existing gate) | SANISAND / DP | B-bar hex | Regression: at the De the campaign will declare, the associated-limit result must move only by the measured De-family, not by a defect |

**Measurements per leg:** band width (a threshold-free definition is OQ5 — candidate: second
moment of the deviatoric plastic-strain-rate profile across the band at mid-height), peak load,
post-peak curve, De = τ/T, `ladrunoBuild()`.
**Gates:** τ = 0 → width ∝ h (the negative control must *fail* objectivity — vault 65 D5's
own words); τ > 0 at fixed De → width(h₄)/width(h₂) within a band TIMs declares, and the
load–displacement curves converge; De × {½, 1, 2} → width scales monotonically and matched
(τ, T) pairs collapse (C6).
**Parameters TIMs must supply (OQ2):** target band width relative to B, ramp duration / strain
rate, the De they will run at, and the ultimate criterion (inherited OQ1, Prof. Gorini).

---

## 6. Phases and exit gates (ADR-31 P3 template: oracle → port → tier)

| Phase | Content | Exit gate | Warrant items |
|---|---|---|---|
| **P0 — ADR + numpy oracle** | Write ADR 90 with §3's alternatives table and §4's claim list. Numpy: two-track DL over a J2/DP point model, PV1–PV6 clones; **A0 bar** two-track vs true DL. | PV1–PV6 green; A0 shows width convergence for τ > 0 and ∝ h for τ = 0, and quantifies two-track vs true-DL difference on a hardening law | none (docs + `tests/_testbed`) |
| **P1 — C++ wrapper** | `LadrunoDuvautLions` (name OQ8), ND tag 33022, adopting ctor and forwarding modelled on `StagedStrainNDMaterial` (33014); `getCopy(type)` only; copy-not-alias; `setParameter` forward; `ops_Dt` with `dt ≤ 0 ⇒ inviscid`; tangent blend; τ as a `Parameter` (ramps for De sweeps — `-eta` cannot); response tokens `overstress`, `beta`, `dt`. Tcl + Py. | g++ byte check; Zone-A: byte, overstress, non-tautology, stage-forwarding over SANISAND, roundtrip, **mutation gate red** | manifest YAML, mutation build flag, guide stub, ledger rows, quirks entry (F3 converse), fix F8 |
| **P2 — acceptance case** | Legs A1→A3 in `tests/` (slow tier) with the G5-style positive + negative pairing; De sweep. Lane decision §3 (a)/(b) made here on measurement. | C4, C5, C6 green; A4 regression unchanged within the De family | slow-tier docstring wall times; capacity rule |
| **P3 — TIMs integration** | Provenance block fields (τ, De, rate) in the campaign's output schema; one SFIM/APE radial probe with declared τ on the chosen lane; out-of-family verdict. | Joint sign-off; vault note recording the acceptance case as *verified*, not claimed | `reviews/handoff_adr90.md` → verdict |
| **P4 — optional** | Explicit/transient tier (`CentralDifferenceLadruno`), `-implex` interplay, mesh-aware τ via the `ops_TheActiveElement->getCharacteristicLength()` latch (`LadrunoConcrete3D.cpp:347-350`). | only on a named consumer | — |

P0 and P1 are fork-only. P2's numbers and P3 are joint. Effort order-of-magnitude: P0 days,
P1 one WP, P2 the larger half (S7) dominated by run time, not code.

---

## 7. Decisions the ADR must make (draft D-list with recommendations)

| # | Decision | Recommendation | Why |
|---|---|---|---|
| D1 | Wrapper scope | **Plasticity-type inner materials only**; damage/plastic-damage models documented as out of scope (cannot be detected at the seam — document, do not guard) | F2; ADR-31 §4.4 precedent |
| D2 | Relaxation form | Two-track nominal-stress blend with trial $\sigma_n^{vp} + C_e\!:\!\Delta\varepsilon$; **P0/A0 decides** whether true DL is needed for SANISAND, and if so the fallback is a material-specific `-tau` inside `LadrunoSANISAND` (fork-owned subclass; base members are `protected`) | F2 |
| D3 | $C_e$ for the trial predictor | `inner.getInitialTangent()` read **before** `setTrialStrain` (committed-state $C_e$; SANISAND rewrites it every step, `ManzariDafalias3D.cpp:146-149`) | F5 |
| D4 | Δt source and fallback | `ops_Dt`; `dt ≤ 0 ⇒ β = 1` (inviscid), the fork convention (`LEDGER_quirks.md:1612`); plus a **non-uniform-Δt diagnostic** (warn once when Δt varies by > X% between commits) | F3 |
| D5 | Calibration lane | §3 (a) transient with β_K = 0 primary; (b) uniform-pseudo-time static fallback; gated equal at matched De | F3, F4, vault 14 |
| D6 | τ mutability | `Parameter` hook (`setParameter("tau")`) so De sweeps and staged activation do not need re-construction | S5 |
| D7 | Element for the final acceptance leg | B-bar hex (A3) — no result on another element is admissible for the campaign | S3 |
| D8 | Command shape | `nDMaterial LadrunoDuvautLions $tag $innerTag -tau $tau` (flags after positionals, ADR-86 parser rule) | 86 emitter guide |
| D9 | Off-switch | τ = 0 byte-identical **including instruction path** (ADR-31's "same instructions" rule) | ADR-59 B2 |

---

## 8. Open questions (owner in brackets)

- **OQ1** [Prof. Gorini] — the ultimate criterion (inherited vault OQ11/OQ1). Scales P1 by 3.3×; the acceptance case can run without it but its $q_u$ gate cannot be cited until it lands.
- **OQ2** [TIMs] — the acceptance-case numbers: target band width / B, ramp duration, De, tolerance bands. The fork drafts (§5); TIMs sets.
- **OQ3** [fork, P0] — is two-track DL an adequate width regularizer on a state-dependent hardening law, or does SANISAND need true DL (internal-variable relaxation)? A0 measures it.
- **OQ4** [both, P2] — lane (a) vs (b) on the actual radial probe; and whether vault 14's β_K family and the wrapper's τ family collapse onto one De when both are present (they should — same Kelvin–Voigt structure — but that is a claim to test, not assume).
- **OQ5** [fork] — a threshold-free band-width metric usable across DP and SANISAND and across element types.
- **OQ6** [fork] — interaction with SANISAND's own substepping / `-Pmin` whole-tensor resets: the inner is called once per trial and is unaware of the wrapper, so nominally none — verify on A2.
- **OQ7** [fork, prerequisite] — cover `LadrunoSANISANDPlaneStrain::getCopy(void)` and fix the null-ctor classTag quirk before A2.
- **OQ8** [fork] — class/command name: `LadrunoDuvautLions` (says what it is) vs `LadrunoViscoplastic` (leaves room for Perzyna). Recommendation: the former; a Perzyna variant would be a different class.

---

## 9. Risk register — mapped to the ADR-59 findings it must not repeat

| # | Sev | Risk | ADR-59 analogue | Mitigation |
|---|---|---|---|---|
| R1 | BLOCKING | Two-track ≠ Duvaut–Lions for hardening materials; the oracle passes on perfect plasticity and the claim silently over-reaches | B3 (wrong regularized variable) | A0 in P0 compares two-track vs true DL on a hardening law **before** C++ |
| R2 | BLOCKING | Regularization inert or erratic on the campaign's static lane (F3/F4); "mesh-converged" claimed from a lane where De is undefined | B2 (off-switch semantics) | §3 lane decision + Δt-uniformity diagnostic + De reported per run |
| R3 | MAJOR | Band width is De-dependent; a reviewer reads "converges in h" as "mesh-independent" | headline over-sell | Claim wording §3; C6 disclosure gate is mandatory, not optional |
| R4 | MAJOR | Width metric is threshold-dependent and reverses across meshes (note 82's checkerboard lesson) | M8 (boundary-layer artefact) | OQ5 settled in P0 on the A0 bar where the width is known analytically |
| R5 | MAJOR | 2D result on `quad` does not transfer to B-bar hex | M4 (cost/mesh economy) | A3 is the admissible leg; A1/A2 are ladders, not deliverables |
| R6 | MAJOR | Static-buffer aliasing / `getCopy(void)` / stage-flag forwarding — each a silent wrong answer | (new) | F5 traps become Zone-A tests in P1; `test_fspm_over_manzari_family` twin |
| R7 | MINOR | Slow-tier cost (R3 gate is 61 min; A3 at three meshes + De sweep will be hours) | — | wall times in docstrings; nightly Zone-B for the sweep, Zone-A slow for one leg |

---

## 10. Deliverables and ledger obligations

- `Ladruno_implementation/90_ladruno_viscoplastic_regularization_adr.md` (P0) — number 90 allocated on creation; this brief then becomes its planning appendix.
- `LEDGER_implementations.md` row (ADR 90, planning → P0 → …); `README.md` index line.
- `LEDGER_quirks.md`: (i) rate artifact when τ ~ static pseudo-step (F3 converse); (ii) two-track vs true DL semantics; (iii) the F8 stale-doc corrections.
- P1: `ledger/ladruno_duvaut_lions.yml`, `-DLADRUNO_MUTATE_DUVAUTLIONS`, guide stub, both interpreters, `SRC/classTags.h` 33022 with the `// N. Mora-Bowen (Ladruno) — …` comment, `banner_features.txt` line on ship.
- A vault-side note (TIMs) recording the acceptance case as **specified**, then as **verified** after P2 — the fork does not write vault notes, it hands the draft over.

---

## 11. Next actions, in order

1. **Hand §5 to TIMs** as the fork's draft of the OQ5 acceptance case, asking for OQ2's numbers. Unblocks vault 65's first work item.
2. **P0-A0 oracle** (two-track vs true DL, 1-D bar width convergence) — a day, decides D2.
3. **Write ADR 90** against the A0 result and §3's alternatives table; open the `wp/<n>-adr88-duvaut-lions` draft PR with the ADR + oracle only.
4. Fix OQ7 prerequisites on `LadrunoSANISANDPlaneStrain` as a small separate PR.
5. P1 only after the ADR passes an out-of-family read against §9.
