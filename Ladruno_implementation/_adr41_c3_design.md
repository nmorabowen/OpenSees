---
title: ADR-41 C3 design / handoff — frictional mortar contact (Coulomb/Tresca on the λ_T form)
project: Ladruno
status: C3 COMPLETE — C3.0 design+oracle (#376) + C3.1 force (#377) + C3.2 symmetric tangent (#378) + C3.3 λ_T Uzawa & non-sym Csl (#379) shipped. Penalty + ALM friction, symmetric + consistent tangent, explicit + implicit. Next track: C4 mesh-tying.
owner: nmora
tags:
  - implementation
  - contact
  - mortar
  - friction
  - alm
  - handoff
---

# ADR-41 Track C3 — frictional mortar: wire the shipped friction kernel onto the mortar interface

> **START HERE (next session).** C2 shipped frictionless mortar ALM (#373/#374/#375): the `MORTAR`
> `LadrunoContactFE` mode assembles `F^s=−(D·p)n / F^m=+(Mᵀ·p)n` with `p_I = min(0, λ_I + epsN·ḡ_I)`,
> the per-global-slave-node normal multiplier `λ_N` lives on `LadrunoContactDomain::MortarNormalState`
> (committed-only, Uzawa once per `Domain::commit()`). **C3 adds the TANGENTIAL traction** — Coulomb/
> Tresca friction — by reusing the SHIPPED `LadrunoFrictionKernel` return map on a per-global-slave-node
> tangential slip + multiplier `λ_T`, exactly mirroring the `λ_N` home. C4 is mesh-tying. C3 adds **no new
> mechanics kernel** (the friction kernel shipped at A1 #367); it is enforcement + integration only.
>
> Read first: this doc; then the C2.2 handoff (`_adr41_c2_design.md` §"C2.2 handoff" + §"SHIPPED
> resolution") because C3 mirrors that exact per-node Domain-state pattern; the shipped
> `LadrunoFrictionKernel.h` (the return map + tangent block C3 calls); ADR-41 §"FrictionalLaw interface"
> + §"Coulomb vs Tresca = one cone"; the capstone [[48_ladruno_contact_capstone_adr]] C3 row.

## The decision that frames everything — per-node λ_T, NOT per-GP (mirror the shipped C2.2 λ_N)

ADR-41 §architecture sketched per-GP friction state on a `LadrunoMortarPair` object. **That object was
never built** — C2 keeps the narrow phase INLINE in the adapter and the path state on the Domain, keyed
by a rebuild-stable id (the shipped NTS + mortar pattern; there is no `LadrunoContactPair`/`LadrunoMortarPair`).
C3 follows the SHIPPED reality, not the ADR sketch:

- **Tangential state is per GLOBAL slave node `I`**, keyed `(contactTag, slaveNodeTag)` — the SAME key as
  `MortarNormalState`. A slave node gets a committed elastic-slip vector `gpT_I[3]` (in the 3D tangent
  plane), an ALM tangential multiplier `λ_T_I[3]`, and an engagement origin `gT0_I[3]`. This is the
  tangential analogue of `λ_N`, and it is variationally consistent at shared nodes for the same reason
  `λ_N` is (the multiplier is per global node; the per-facet `D/M` operators assemble it globally for free).
- **Per-GP slip is NOT stored.** The mortar-weighted tangential gap at node `I` is integrated through the
  same `D/M` operators that give the weighted normal gap — there is one slip vector per node, not per GP.
  (Per-GP would re-introduce the rebuild-instability the clip-polygon GPs have — they move every step.)

**Why this is the right MVP (and what it sacrifices).** Standard (not dual/biorthogonal) basis ⇒ `D` is
non-diagonal, so the per-node tangential traction is a weighted combination, not a clean nodal decoupling.
That is the SAME approximation C2 made for the normal multiplier and is fenced to ADR-47 (dual basis). It
passes the constant-friction patch test for the same `Σφ=1` reason the normal pressure does.

## The mechanics (per-node mortar friction, made concrete)

C1 gives per slave/master facet pair: `D_IJ`, `M_IK`, the weighted **normal** gap `g̃_I`, and the per-GP
normal `n`. C2 added the per-node normal pressure `t_N,I = −(λ_N,I + epsN·ḡ_I)` (a compression magnitude
≥ 0; the force is `−(D·p)n`). C3 adds the tangential traction:

1. **Weighted tangential gap (slip), per node — from DISPLACEMENTS, not positions.** The closest-point
   projection makes the weighted relative POSITION `∫N_I(x_s − x_m(ξ̄))` purely NORMAL (`n·r = g̃`), so
   positions carry **no** tangential slip information (the ADR-39 SEGMENT path hit this exact trap). The
   tangential slip is the tangential part of the weighted relative **displacement**:
   `ḡ_T,I = P_t · ( (Σ_J D_IJ u_s,J − Σ_K M_IK u_m,K) / a_I ) − gT0_I`, where `P_t = I − n⊗n`,
   `a_I = Σ_J D_IJ`, `u` are nodal DISPLACEMENTS (`getTrialDisp`), and `gT0_I` is the **engagement-config**
   weighted tangential displacement captured once at first contact (the C3 analogue of the NTS `gT0` — else
   a late-engaging node's pre-contact drift becomes a spurious stick traction; ADR-39 MAJOR-1, reused).
   **The slip fed to the return map is `gTeff_I = ḡ_T,I`.** (The ADR-39 SEGMENT form is `u_s − Σ N_i u_i`,
   tangential — the mortar form is the `D/M`-weighted generalisation.)
2. **Per-node return map (the SHIPPED kernel, verbatim).** Call
   `LadrunoFrictionKernel::frictionReturnMap(gTeff_I, gpT_I, N_I, epsT, mu, tFric_I, gpTtrial_I, cohesion, tauMax)`
   with the contact pressure `N_I = t_N,I = epsN·⟨−ḡ_I⟩ + (−λ_N,I)` (the C2 nodal normal pressure magnitude;
   the same one driving the normal force) and the tangential penalty `epsT` in the `kt` slot. The kernel
   returns the **already-negated** applied tangential force `tFric_I[3]` (opposes the slave motion — the
   A1 BLOCKER-1 sign lives in the kernel) and the trial elastic slip `gpTtrial_I` (a pure function of
   committed state ⇒ idempotent across re-evals). The unified cone `cap = min(μ·N_I + c, τmax)` is the
   kernel's; Tresca is `μ=0`, pure Coulomb is `τmax≤0`.
3. **Assemble like the normal force, but with the tangential traction vector** (not `±t·n`):
   `f^s_K += Σ_I D_KI · tFric_I` (a 3-vector per node), `f^m_L += −Σ_I M_IL · tFric_I`. Self-equilibrating
   by `Σφ=1` exactly as the normal force. (NOTE the tangential traction is a full 3-vector lying in the
   tangent plane, not a scalar×n — so it scatters as a vector, the only structural difference from C2.)
4. **Tangent** (C3.2): the friction block `K_ss = LadrunoFrictionKernel::frictionTangentBlock(...)` (3×3,
   `D_TT·P_t` + the optional non-symmetric `Csl = −∂cap/∂N·kn·n̂⊗n` Coulomb coupling) scattered through the
   mortar `B̃ = [D, −M]` operator: `K_c += B̃ᵀ (K_ss ⊗ per-node) B̃`. **Default SYMMETRIC** (drop `Csl` —
   any SOE is safe; the ADR-39 P3.5 Q2 BLOCKER: a symmetric solver silently drops the lower triangle);
   `-consistanttan` opts into the full non-symmetric Coulomb tangent (needs FullGeneral/UmfPack). The
   geometric `∂{D,M,n}/∂u` terms stay the deferred C1 stub (as in C2).
5. **λ_T outer update (ALM), once per `Domain::commit()`** — the tangential Uzawa, beside the `λ_N` update:
   `λ_T,I ← P_cap( λ_T,I + epsT·ḡ_T,I )`, where `P_cap` is the return-to-the-friction-cone projection
   (so the augmented tangential multiplier never exceeds `cap`). **MVP simplification (recommended):**
   start with **penalty-only friction** (`λ_T ≡ 0`, the `lambdaT=null` path the kernel already supports)
   — the normal ALM already gives the headline accuracy win, and frictional tangential drift is bounded
   by `epsT`; the tangential Uzawa is a refinement gated by a named test (Δt-independent tangential
   converged answer). Ship C3.1/C3.2 penalty-friction first (like NTS P3/P3.5), then add the `λ_T` Uzawa
   if a gate needs it. This keeps `commit()` change minimal: promote `gpT = gpTtrial` per node (the slip
   commit), exactly mirroring the shipped NTS `FrictionState` promotion.
6. **Active set / stick-slip.** A node is in NORMAL contact iff `λ_N,I + epsN·ḡ_I < 0` (the C2 mask,
   frozen within an augmentation sweep). Friction acts ONLY on normal-active nodes (`N_I > 0`). Stick vs
   slip is the kernel's per-node radial-return decision (local, re-evaluated each residual) — it does NOT
   change the SOE size or the normal active set, so the eqn-count-constant gate is preserved.

## The seams it plugs into (SHIPPED — grounded)

Mirror these; do not invent new architecture (the C2 lesson).

- **State on the Domain** (`LadrunoContactDomain`): add a `MortarFrictionState{gpT[3], gpTtrial[3], gT0[3],
  lambdaT[3], engaged}` map keyed `(contactTag, slaveNodeTag)` — the SAME `NodeKey` as `MortarNormalState`
  (consider folding the tangential fields INTO `MortarNormalState` to share the slot/GC, since they share
  the key and lifecycle). `commit()` promotes `gpT = gpTtrial` (+ the optional `λ_T` Uzawa); `revertToLast
  Commit()` drops `gpTtrial = gpT`. Reuse the `mortarNormalGC*` machinery verbatim (same key). **Committed-
  only `λ_T`** (mutated only in `commit()`), the EmbeddedRebar invariant — revert is then automatically safe.
- **Adapter `MORTAR` getResidual** (`LadrunoContactFE.cpp`): after the normal block, when `mu>0 && cd!=0`,
  for each slave node `I` build `gTeff_I` (the mortar-weighted tangential gap minus `gT0_I`, captured at
  first activation), call `frictionReturnMap` with `N_I` = the nodal normal pressure already computed,
  write `gpTtrial_I` to the slot, and add `Σ_I D_KI tFric_I` (slave) / `−Σ_I M_IL tFric_I` (master) to the
  residual. `mu≤0` short-circuits before any slot touch ⇒ byte-identical to frictionless C2 (the A1/P3
  byte-identity contract). The nodal tangential gap needs the per-node weighted position — accumulate it
  on the Domain like the normal gap (a `gtT[3]` accumulator beside `gtGlobal`), OR (simpler, since the
  return map is per-node and the slip is small) compute the facet-local tangential gap and feed the LOCAL
  gap (the C2.2 shipped resolution: local gap in the force, global only for the commit/query — the SAME
  order-dependence trap applies to a running global tangential gap, so use the LOCAL tangential gap in the
  force/return-map and reserve a global accumulation for the optional `λ_T` Uzawa + a query). **This is the
  C3 crux to pin in the oracle** (see oracle T-LOCAL below).
- **Adapter `addMortarTang`** (C3.2): scatter `K_ss` (per node) through `B̃=[D,−M]` ⊗ the tangent-plane,
  default symmetric. Reuse the C2 `addMortarTang` scatter skeleton (it already loops `B̃` over slave|master
  nodes); add the friction `K_ss` contribution beside the normal `epsN·diag(act/a)` block.
- **Handler**: no structural change — the same mortar pairing loop; thread the friction params (`mu`,
  `epsT`, `cohesion`, `tauMax`, `consistentTan`) from the `MortarContact` definition into the `MORTAR`
  ctor (extend the ctor arg list, as P3 did for the SEGMENT ctor). Mark the friction slot live with the
  existing `mortarNormalGCMark` (same key).
- **Command surface** (`OpenSeesOutputCommands.cpp`): the `-mortar` parse gains `-mu <v>`, `-epsT auto|<v>`,
  `-cohesion <v>`, `-tauMax <v>`, `-consistanttan` (mirror the NTS `contact` friction options). `-epsT auto`
  ← `γ_crit = 0.5%·L` sizing (ADR-41 A1 note) or simply `epsT = epsN` as a first default. Store on
  `MortarContact`. `mu=0` (default) ⇒ the frictionless C2 path, byte-identical.

## Phasing (mirror P3/P3.5 and C2)

- **C3.0 — design + oracle (THIS note + `proto_c3_mortar_friction.py`).** Pin the per-node mortar-friction
  mechanics in numpy: stick/slip return map driven through `D/M`, Coulomb incline `a=g(sinθ−μcosθ)`, Tresca
  cap, `Csl` vs FD (both sides of the `min`), self-equilibrium, the LOCAL-vs-global tangential-gap crux
  (T-LOCAL: local-gap friction force is deterministic + converges; running-global is order-dependent — the
  exact C2.2 lesson), `μ=0` ⇒ frictionless regression. Gate: oracle all-pass.
- **C3.1 — mortar friction FORCE (penalty, λ_T≡0).** The load-bearing first phase. `MortarFrictionState`
  on the Domain + the `getResidual` tangential block + the command-surface friction options + handler
  threading. Gates: incline `a=g(sinθ−μcosθ)` through a REAL solve (the sign gate); slip caps at `cap`;
  stick holds; `μ=0` byte-identical to C2; commit promotes slip; full battery green. Tangent UNCHANGED
  (explicit/CDL force-only, like P3) ⇒ implicit needs C3.2.
- **C3.2 — consistent friction TANGENT.** `addMortarTang` friction block (symmetric default + `-consistanttan`
  non-symmetric Coulomb). Gate: static frictional Newton CONVERGES (singular without it, the P3.5 gate);
  symmetric-solver-safe; `Csl` FD-checked on a slipping config; Newmark dynamic. **Two C3.1-gate TODOs that
  go live here (see [[LEDGER_quirks]]):** (1) `revertToLastCommit` must also revert `gT0`/`engaged`
  (double-buffer them) — a rejected implicit step otherwise latches a stale engagement origin; (2) the
  shared-node committed-slip order dependence (MAJOR-1) should get a non-matched friction regression +
  a per-node slip reconciliation before non-matched frictional meshes are trusted.
- **C3.3 (optional) — tangential Uzawa `λ_T`.** Only if a named gate needs the Δt-independent tangential
  converged answer beyond what penalty `epsT` gives. The `analyze_augmented` proc already drives it (the
  commit-cycle outer loop augments `λ_T` beside `λ_N`).

## Scope discipline (don't pull ADR-47 / C4 into C3)

- **Standard basis** (non-diagonal `D`); dual/biorthogonal, true-LM saddle-point, anisotropic friction,
  `mu(N,v)` rate/state friction wiring, self-contact = **ADR-47**.
- **No new DOFs / no custom algorithm** — `λ_T` is Domain-side per-node state, Uzawa on commit (the C2
  resolution, reused). Friction stick/slip is local (return map) ⇒ no SOE-size change.
- **Geometric tangent `∂{D,M,n}/∂u` deferred** (the C1 stub) — the material friction `K_ss` carries C3.2,
  exactly as C2 shipped without the geometric block.
- **Mesh-tying (`-tie`) = C4** (zero-gap, active set frozen ON).

## Gates (capstone C3 / ADR-41)

- incline `a = g(sinθ − μcosθ)` converges on a real mesh (the friction SIGN gate — friction OPPOSES motion);
- slip traction caps at `cap = min(μN+c, τmax)`; a sticking node holds (no drift beyond `cap/epsT`);
- Tresca cap (`μ=0`) and Coulomb+cohesion+cap (`min` both sides) both correct; `Csl` FD-checked (C3.2);
- `μ=0` ⇒ **byte-identical** to the frictionless C2 path (the regression contract);
- static frictional Newton CONVERGES (C3.2 — singular without the friction tangent);
- Δt-independent converged tangential answer (C3.3, if the λ_T Uzawa is built);
- SOE eqn count constant (stick/slip is local — no re-number);
- full contact battery stays green.

## References

- ADR-41 (`41_ladruno_mortar_alm_contact_adr.md`) §"FrictionalLaw interface", §"Coulomb vs Tresca = one
  cone", §"Mortar integration scheme", §"Integration points" (getTangent/c1, the unsymmetric `Csl`).
- Shipped seams: `LadrunoFrictionKernel.h` (`frictionReturnMap` + `frictionTangentBlock` + `frictionCap` +
  `tangentPart` — the verbatim kernel C3 calls), `LadrunoContactFE.{h,cpp}` (the `MORTAR` mode + the NTS
  `SEGMENT` friction wiring to mirror), `LadrunoContactDomain.{h,cpp}` (`MortarNormalState` + `mortarNormalGC*`
  + the commit Uzawa — extend in place), `LadrunoContactHandler.cpp` (mortar pairing loop).
- C2 precedent: `_adr41_c2_design.md` (the per-node Domain-state pattern + the LOCAL-vs-global-gap crux
  that C3 re-encounters for the tangential gap), [[LEDGER_quirks]] ("Per-facet mortar adapters reading a
  RUNNING global gap").
- Literature: Wriggers *Computational Contact Mechanics* (frictional augmented Lagrange); Popp/Gee/Wall
  (mortar friction). Skills: `fem-mechanics-expert`, `opensees-expert`.
- Workflow: oracle runs with **`python3.12`**; build `Ladruno_scripts\build.bat OpenSeesPy`; watch **Zone-A**
  before merging C++ (C3.1/C3.2 touch real TUs). See [[ladruno-adr41-mortar-c2]], [[ladruno-oracle-python]].

## Oracle-first checklist (do these in `proto_c3_mortar_friction.py`, in order)

1. **Return-map reuse**: the shipped `frictionReturnMap` logic (re-included standalone) gives stick
   (`‖tT*‖≤cap` ⇒ `tFric=−tT*`, slip preserved) and slip (`tFric=−cap·n̂`, `gpT += dλ·n̂`); FD-check the
   3×3 `K_ss` (stick `kt·P_t`; slip `D_TT·P_t` + the optional `Csl`) to 1e-6 — both sides of the `min` cone.
2. **Per-node assembly through D/M**: on a non-matched mesh, the tangential force `Σ_I D_KI tFric_I` is
   self-equilibrating (`ΣF^s + ΣF^m = 0`) and reduces to the consistent nodal traction for a constant slip.
3. **T-LOCAL crux**: the LOCAL-gap tangential force (per facet) is deterministic and converges; pin that a
   running-global tangential gap would be order-dependent (the C2.2 lesson, re-shown for `λ_T`).
4. **Incline**: a block on a μ-incline reaches `a = g(sinθ − μcosθ)` (slip) / stays (stick) — the sign gate.
5. **Tresca / Coulomb+cohesion+cap**: cap selection across the `min`; `μ=0` ⇒ frictionless (tFric=0).
6. **(If C3.3) λ_T Uzawa**: the tangential multiplier outer update drives the stick drift → 0 at finite epsT,
   Δt-independent.

Only after these pass in numpy, wire C3.1 (force) then C3.2 (tangent) into the C++ seams above, re-confirm
via the full contact battery + Zone-A.
