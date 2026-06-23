---
title: ADR-41 Track D2 design / handoff — viscous contact stabilization (p_visc = μ_c·v_rel)
project: Ladruno
status: D2 COMPLETE — D2.1 (NTS: rigid-plane + segment) #385 + D2.2 (mortar contact) #TBD shipped. Viscous stabilization on NTS + mortar; refused on -tie. Battery 83/83.
owner: nmora
tags:
  - implementation
  - contact
  - stabilization
  - viscous
  - handoff
---

# ADR-41 Track D2 — viscous contact stabilization

> Resolves the OPEN **Q-VISCOUS** residue ([[41_ladruno_mortar_alm_contact_adr]] §Q-VISCOUS, capstone
> [[48_ladruno_contact_capstone_adr]] Track D2 — "funded option, best ROI of the residue"). Picking D2
> as the next track FUNDS it: a gated, off-by-default `-visc <μ_c>` option, the flagship
> **pounding / rocking / uplift** chatter/snap-through damper. Abaqus TG §5.2.1.

## What it is

A velocity-proportional **normal** contact pressure that damps contact-status chatter and snap-through
near the active-set threshold (where penalty contact rings/flip-flops): `p_visc = μ_c · ġ`, where
`ġ = ḋ(gap)/dt` is the **normal relative velocity** (approach rate) at the contact point. It opposes the
approach/separation rate, bleeding the high-frequency energy that drives chatter, without changing the
converged static solution (in statics nodal velocity ≡ 0 ⇒ the term is identically inert).

It is NOT a material/structural damper (that is Rayleigh/modal on the elements); it is a **contact-local**
stabilizer active only on the contact interface, only while the pair is in contact.

## The mechanics — it reuses the normal-penalty operator, driven by velocity

The shipped normal penalty force is `f = −kₙ·g·B`, tangent `K = kₙ·B Bᵀ`, where `B = ∂g/∂u` is the gap
gradient (per mode: RIGID_PLANE `B = n`; SEGMENT `B = [n | −Nᵢ n]`; MORTAR `B̃ = [D,−M]⊗n` with `W=1/a`).
The gap **rate** is `ġ = B·v` (v = nodal velocities, `getTrialVel()`). So the viscous term is the SAME
operator with `kₙ → μ_c` and `u → v`:

  - **viscous force** (into the residual):  `f_visc = −μ_c · (B·v) · B`   (= `−μ_c·ġ·B`)
  - **viscous damping matrix** (into C):    `C_visc = μ_c · B Bᵀ`         (= `μ_c·∂ġ/∂v` ⊗ B)

Assembly seam (the adapter is a self-contained FE_Element view — it does NOT route through
`formEleResidual`/`getDamp`):
  1. **Residual** — add `f_visc` in `getResidual()` from `getTrialVel()`, beside the penalty/friction
     force. (`v ≡ 0` in statics ⇒ `f_visc ≡ 0` ⇒ static byte-identical.)
  2. **Tangent** — `addCtoTang(fact)` adds `fact·C_visc`. The transient integrators call
     `theEle->addCtoTang(c2)` in `formEleTangent` (Newmark/HHT/CentralDifference, verified), `c2` the
     damping coefficient. So implicit Newton gets the consistent damping tangent; the fork's mass-only
     `CentralDifferenceLadruno` never calls `addCtoTang` ⇒ **force-only under CDL** (exactly the P3
     friction explicit/implicit split).

**Active mask:** viscous fires ONLY where the normal force is active (RIGID_PLANE/SEGMENT `gap<0`; MORTAR
the per-facet `p_I<0` mask) — a separated pair has no contact and no damping (never damps a free-flying
node). **Not applied to `-tie`** (a permanent bond has no approach velocity to damp).

## Command surface

`contact … -visc <μ_c>` (a new trailing option on the existing `contact` parse), off by default (`μ_c=0`
⇒ byte-identical to no-viscous, the regression contract). Stored on `Contact`/`RigidPlane`, threaded into
the adapter ctor. `μ_c` units: pressure·time/length (a normal-dashpot coefficient); a physical default is
a small fraction of critical (`μ_c ≈ ζ·2√(kₙ·m)` order), but the MVP just takes the user value.
**Naturally inert in statics/arc-length** (velocity = 0); no explicit lane-disable needed (documented).

## Phasing

- **D2.1 (this session) — NTS: RIGID_PLANE + SEGMENT.** The flagship pounding/rocking/uplift regime
  (rocking/uplift on a foundation plane; deck/structure pounding via NTS facets). Viscous force in
  `getResidual` (from `getTrialVel`) + `addCtoTang` C matrix; `-visc` command + ctor threading. `-visc`
  on a MORTAR contact is REFUSED (clear message: "viscous stabilization is NTS-only in D2.1; mortar = D2.2")
  rather than silently ignored. Gates: explicit impact with `-visc` dissipates energy (e<1) vs `-visc 0`
  (e≈1) — chatter/restitution control; static solve byte-identical (vel=0 ⇒ inert); implicit Newmark
  converges with the consistent C tangent; `-visc 0` / `-visc`-absent byte-identical to the shipped path.
- **D2.2 (follow-on) — MORTAR.** Mirror the C2 `addMortarTang` penalty block with `epsN→μ_c`, `g→ġ`,
  the same active mask; pull only if a mortar chatter gate needs it (the mortar accuracy lane is
  implicit/quasi-static, far less chatter-prone than explicit NTS pounding).

## Scope discipline

- **Normal viscous only** (the chatter/snap-through driver). Tangential viscous is not Q-VISCOUS.
- **Off under arc-length** — naturally inert (vel=0), documented; if a future arc-length lane synthesises
  a pseudo-velocity, gate it off explicitly then.
- **Not a substitute** for ALM accuracy or friction — purely a high-frequency-chatter bleed.

## Gates (D2.1)

- explicit impact: `-visc μ_c>0` ⇒ measurable energy dissipation / restitution drop vs `-visc 0` (e≈1);
- static solve with `-visc μ_c>0` == without (velocity ≡ 0 ⇒ viscous identically zero — byte-identical);
- implicit Newmark with `-visc` CONVERGES (consistent `C_visc` via `addCtoTang`) and damps the transient;
- `-visc 0` / `-visc`-absent ⇒ byte-identical to the shipped contact (the regression contract);
- `-visc` on a mortar contact ⇒ clean refusal (D2.2 scope);
- full contact battery stays green.

## References

- ADR-41 §Q-VISCOUS, §Q-NORMAL (the chatter mitigation note); capstone [[48_ladruno_contact_capstone_adr]]
  Track D2 + the open question "D2 earns a committed phase". Abaqus Theory Guide §5.2.1 (viscous damping).
- Seams: `LadrunoContactFE.{h,cpp}` (RIGID_PLANE/SEGMENT getResidual + addCtoTang + the B/n operators
  already computed for the penalty force), `LadrunoContactDomain.{h,cpp}` (`Contact`/`RigidPlane` defs +
  the `-visc` field), `OpenSeesOutputCommands.cpp` (the `contact` parse). `Node::getTrialVel()`.
- Workflow: oracle `python3.12` `proto_d2_viscous.py`; build `Ladruno_scripts\build.bat OpenSeesPy`;
  watch Zone-A; PR on `ladruno`. See [[ladruno-adr41-mortar-c2]], [[ladruno-oracle-python]].
