---
wp: LEGACY
title: "Assumed-strain hourglass: the dev-projection vs reduced-shear interaction is nu-coupled"
legacy_seq: 17
---
### Assumed-strain hourglass: the dev-projection vs reduced-shear interaction is nu-coupled
- **Bites:** building the LadrunoBrick `physical` hourglass straight from Belytschko
  eq 8.7.26 (pointwise-*isochoric* assumed strain: 2/3,-1/3 dev-projection on the
  normal hourglass strains + mode-subset reduced shear) gave an element that was
  ~75% **too soft** in bending and got *worse* with nu (ratio 1.73 @nu=0 -> 2.59
  @nu=0.499 vs the analytic 1.0). Patch test + rank were exact, so the bug hid
  from the usual gates.
- **Why:** the dev-projection IS the proper B-bar mean-dilatation treatment for
  the hourglass normals (the algebra collapses to the same 2/3,-1/3). On its own
  (with full shear) it's fine; combined with the **reduced** assumed shear it
  removes too much energy -> over-soft, and the error scales with lambda(nu).
  Dropping the dev-projection (FULL compatible normal strains + reduced shear)
  gives a **correct shear-locking cure** — matches OpenSees `SSPbrick` to 3 digits
  and converges (0.94->1.005, nx=2..32) at nu=0 — but then VOLUMETRIC-locks at
  nu->0.5. There is **no single static projection** correct across nu with this
  shear field; a general all-nu element needs the coupled SSP/ASQBI operator
  (Belytschko sec 8.7.8 explicitly: 3D assumed-strain structure "not fully developed").
- **The validating oracle:** patch + rank CANNOT validate an assumed-strain
  element (gamma-orthogonality makes any variant pass). Use a **bending-convergence
  benchmark** and cross-check against `SSPbrick` (a proven OpenSees assumed-strain
  hex, ~1.0 across all nu). `tests/test_ladrunoBrick_bending.py`.
- **Status (2026-06-01):** shipped `physical` as the FULL-normals + reduced-shear
  **shear-locking cure** (verified vs SSPbrick); documented that near-incompressible
  needs `-formulation bbar`. A coupled general-nu operator is future work.
- **The definitive difference vs SSPbrick (read its source).** `SSPbrick.cpp` is an
  **EAS element = bbar + statically-condensed enhanced strain**. (1) Volumetric:
  its constant `Bnot` uses `dNmod` = mean-dilatation (B-bar) modified gradients
  (`SSPbrick.cpp:1254,1266`). (2) Shear/bending: 9 internal **enhanced-strain
  modes** `Fe`, condensed out — `interior = FCF − K_uα·K_αα⁻¹·K_αu`
  (`SSPbrick.cpp:1968`), then `Kstab = Mbenᵀ·interior·Mben`. The **static
  condensation** is why SSP works for ALL nu: the internal modes *adapt to C*. My
  `physical` is a single FIXED assumed-strain B (no internal DOFs/condensation) →
  can cure shear OR volumetric, never both across nu. **Upshot: a general-nu
  "physical" = our reserved `eas` formulation (v2), and SSPbrick is the production
  blueprint (bbar constant part + condensed EAS).** See `SSPbrick.cpp:1053`
  (`GetStab`), `:1243` (G/gamma), `:1647` (enhanced-strain block), `:1968`
  (condensation).
- **SHIPPED (2026-06-01): `LadrunoBrick -formulation eas` is the SSPbrick port.**
  Confirmed while porting: SSPbrick condenses the enhanced modes with the
  **initial** tangent (`GetStab` is called once in `setDomain`), so `Bnot`/`Kstab`
  are **constant** — there is **no per-step α internal state**, contrary to the
  general textbook EAS picture. That collapses the "heavy bit" (no
  `commitState`/`sendSelf` of α): the operators are deterministic from geometry +
  C(0), so the parallel receive side just rebuilds them in `setDomain` and
  `sendSelf` ships nothing extra. Validation gate: for a linear-elastic material
  the assembled `eas` stiffness is *identical* to `SSPbrick`, so the bending-
  benchmark tip matches SSPbrick to ~1e-6 across ν∈{0,0.3,0.45,0.499} (where
  `physical` vol-locks). One caveat: SSPbrick itself sends `Bnot`/`Kstab`/`J[]`
  over `sendSelf` (its null-ctor sets `mInitialize=false` → skips `GetStab` on
  recv); LadrunoBrick instead always rebuilds in `setDomain` — simpler, same
  result, smaller stream.
