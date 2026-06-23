---
title: ADR-41 C4 design / handoff — mesh-tying (permanent mortar bond, zero-gap)
project: Ladruno
status: NOT STARTED — handoff for the next session (C1→C3 shipped; this is the last ADR-41 mortar track)
owner: nmora
tags:
  - implementation
  - contact
  - mortar
  - mesh-tying
  - handoff
---

# ADR-41 Track C4 — mesh-tying: a permanent mortar bond between non-matching meshes

> **START HERE (next session).** C1→C3 shipped the full mortar CONTACT lane (#369→#379): the clip→Gauss
> kernel (`LadrunoMortarKernel.h`, C1), frictionless commit-cycle ALM (per-node `λ_N`, epsN-independent
> penetration, C2.2 #375), and friction (Coulomb/Tresca force + symmetric & consistent tangent + λ_T Uzawa,
> C3.0–C3.3 #376→#379). **C4 is the last track: mesh-tying** — a PERMANENT bond that ties a slave surface to
> a non-matching master surface (no gap, no separation, no sliding). It is the **zero-gap limit of contact**:
> the active set is frozen ON for every slave node, and the FULL relative displacement (normal AND
> tangential) is driven to zero — there is no KKT inequality, no compression clamp, no friction cone.
> It reuses the SHIPPED mortar machinery almost entirely; C4 adds **no new kernel** and **no new DOFs**.
>
> Read first: this doc; the C2.2 handoff (`_adr41_c2_design.md` §"SHIPPED resolution") and the C3 design
> (`_adr41_c3_design.md`) — C4 is the same per-global-slave-node Domain-state + commit-cycle-Uzawa pattern,
> stripped of the inequality; the capstone [[48_ladruno_contact_capstone_adr]] mesh-tying row; the SHIPPED
> `LadrunoContactFE` MORTAR mode + `LadrunoContactDomain::MortarNormalState`.

## What mesh-tying IS (and how it differs from contact)

A tie welds the slave interface to the master interface so they deform as one continuous body across a
non-matching mesh transition (e.g. a fine-meshed region bonded to a coarse one, or dissimilar element
types). The weak (mortar) tie constraint is **zero weighted relative displacement** at every slave node:

  `r_I ≡ Σ_J D_IJ u_s,J − Σ_K M_IK u_m,K = 0`   (a 3-vector per slave node `I`)

Contrast with the shipped contact:

| | Contact (C2/C3) | Mesh-tying (C4) |
|---|---|---|
| Active set | KKT: active iff `λ_N+epsN·ḡ < 0` (compression) | **always ON** (every slave node, every step) |
| Constrained components | normal gap (`λ_N`) + tangential slip within the friction cone (`λ_T`) | **full 3-vector** relative displacement `r_I` (normal **and** tangential), no cone |
| Sign / clamp | `t_N = min(0, …)`; friction capped at `min(μN+c, τmax)` | **no clamp** — equality constraint, `t_I` unbounded |
| Gap reference | the NORMAL gap `g̃_I` (positions); slip from displacements | the FULL relative DISPLACEMENT `r_I` (displacements — same lesson as C3.1) |
| Tangent block | `epsN·(n⊗n)` (normal) + `K_ss` (friction) | `epsTie·(I₃)` — the **full identity**, symmetric, no active-set switching |

So C4 is *simpler* than C3 in the mechanics (no inequality, no cone, no stick/slip) but ties **all three
components**. It is essentially a weak (mortar) multi-point constraint `u_s = Π_mortar u_m`, enforced by
penalty + commit-cycle ALM rather than by the constraint handler's transformation.

## The mechanics (made concrete — reuse C2/C3, drop the inequality)

Per slave/master facet pair, C1 gives `D_IJ`, `M_IK`, `a_I = Σ_J D_IJ`, and the per-facet normal `n` (the
tie does NOT need `n` for the constraint — it ties all components — but the kernel still returns it; the tie
ignores it). The tie traction at slave node `I` is a full 3-vector:

  `t_I = λ_tie,I + epsTie · ( r_I / a_I )`     (no `min`, no projection — equality)

where `r_I = Σ_J D_IJ u_s,J − Σ_K M_IK u_m,K` is the weighted relative DISPLACEMENT (NOT position — the
C3.1 quirk: build it from `getTrialDisp()`; for a tie there is no engagement origin `gT0` because the bond
exists from step 0, so `gT0 ≡ 0` — the reference is the as-built configuration).

1. **Force** (assemble exactly like the friction force — a vector traction via D/−M):
   `f^s_K += Σ_I D_KI · t_I`,  `f^m_L += −Σ_I M_IL · t_I`  (self-equilibrating, `Σφ=1`).
2. **Tangent** (symmetric, no active-set, no `n⊗n` projection):
   `K[A][B] += epsTie · (b_IA b_IB / a_I) · I₃` summed over slave nodes `I`, `b` the `B̃=[D,−M]` row.
   This is the C3.2 scatter with `K_ss → epsTie·I₃` (the full identity, since all 3 components are tied).
   SPD, solver-safe, no `-consistanttan` needed (no pressure coupling — the tie has no friction cone).
3. **ALM (the accuracy win, mirrors C2.2/C3.3)**: one Uzawa step per `Domain::commit()`,
   `λ_tie,I ← λ_tie,I + epsTie · (r_I/a_I)` — **no clamp** (equality). Drives `r_I → 0` at FINITE epsTie
   (epsTie-INDEPENDENT bond — penalty alone leaves `r_I = O(t/epsTie)`). The held-load `analyze_augmented`
   proc augments it; ‖r‖ is the convergence measure (extend `ladrunoMortarPenetration` or add a `‖r‖` query).
4. **No active set, no stick/slip, no GC change** — every slave node of every tie pair is always live.

## The seams it plugs into (SHIPPED — reuse, do not re-architect)

- **State**: reuse `LadrunoContactDomain::MortarNormalState` (or add a parallel `MortarTieState{lambdaTie[3],
  lambdaTietrial[3]}` keyed the same `(contactTag, slaveNodeTag)`). λ_tie is a 3-vector, committed-only,
  Uzawa'd in `commit()` beside `λ_N`/`λ_T`. **Simplest: add `lambdaTie[3]/lambdaTietrial[3]` to
  `MortarNormalState`** (it already carries λ_N + the friction state on the same key) and a `bool isTie`.
- **Adapter**: a `MORTAR` tie MODE (or a flag on the existing MORTAR ctor). In `getResidual`, when `isTie`:
  skip the normal KKT/friction blocks entirely and assemble the tie force `Σ_I D_KI t_I` with
  `t_I = λ_tie,I + epsTie·(r_I/a_I)`, `r_I` from displacements. In `addMortarTang`, scatter `epsTie·I₃` via
  `b·b/a` (the C3.2 scatter with `K_ss=epsTie·I₃`). The geometric `∂{D,M}/∂u` terms stay deferred (as C2/C3).
- **Command surface**: `contact … -mortar -tie [-epsTie auto|<v>] [-augTol …]` (mesh-tie = the `-tie` flag on
  a mortar contact ⇒ `addMortarContact(..., isTie=true)`). `-tie` is mutually exclusive with friction
  (`-mu/-cohesion/-tauMax`) — a tie has no cone (refuse the combination with a clear error). epsTie auto ←
  the owning-solid stiffness (reuse `ladrunoResolveAutoKn`).
- **Handler / GC**: the same mortar pairing loop; a tied pair marks ALL its slave nodes live every handle()
  (mirror `mortarNormalGCMark`). No active-set epoch.
- **classTag**: none — rides the handler, like C2/C3.

## Phasing (mirror C2/C3)

- **C4.0 — oracle** (`proto_c4_mortar_tie.py`): the constant-stress PATCH TEST across a non-matching tied
  interface (the headline — a uniform stress field transmits EXACTLY; `D⁻¹M` recovers the master field, 1e-6);
  penalty tie `r = O(1/epsTie)`; ALM (no clamp) drives `r → 0` epsTie-INDEPENDENT; the assembled tie tangent
  `epsTie·B̃ᵀB̃⊗I₃` FD-checked symmetric. (Reuse the proto_c3 D/M assembly.)
- **C4.1 — penalty tie FORCE + tangent** (`λ_tie ≡ 0`): the `MORTAR -tie` mode, force + the symmetric
  `epsTie·I₃` tangent, the `-tie` command surface, handler threading. Gates: a split-bar / split-cantilever
  across a non-matching tied interface matches the MONOLITHIC solution (displacement + stress) to a penalty
  tol; static Newton converges (the tie tangent is SPD — no singular DOF); `-tie` refuses friction params.
- **C4.2 — ALM tie** (`λ_tie` Uzawa): drives the residual bond `‖r‖ → augTol` at finite epsTie (epsTie-
  INDEPENDENT), the held-load proc; the patch test passes to machine precision under augmentation.

## Scope discipline (don't pull ADR-47 in)

- **Standard basis** (non-diagonal `D`) — dual/biorthogonal tie is ADR-47. The penalty/ALM tie at finite
  epsTie is the mitigation (same as C2/C3).
- **No new DOFs / no transformation-handler tie** — λ_tie is Domain-side per-node state, Uzawa on commit
  (the C2/C3 pattern). A true-Lagrange tie via the constraint handler is NOT this track.
- **Tie is permanent** — no release, no re-pairing across steps (active set frozen ON). Bonded-then-debonding
  (cohesive-zone tie) is a separate future track, not C4.
- **Geometric tangent deferred** (the C1 stub), as in C2/C3.

## Gates (capstone mesh-tying row)

- **constant-stress patch test** across a non-matching tied interface — exact transmission (the headline);
- a split structure (bar/cantilever) tied across a non-matching mesh == the monolithic solution;
- penalty `‖r‖ = O(1/epsTie)`; **ALM drives `‖r‖ → augTol` at FINITE epsTie, epsTie-INDEPENDENT**;
- static Newton converges (SPD tie tangent); SOE eqn count constant (no new DOFs);
- `-tie` ⊥ friction (clean refusal); null/`-tie`-absent path byte-identical.

## References

- Shipped seams: `LadrunoMortarKernel.h` (C1 D,M,g̃), `LadrunoContactFE.{h,cpp}` (MORTAR mode + the friction
  force/tangent scatter to mirror), `LadrunoContactDomain.{h,cpp}` (`MortarNormalState` + commit-cycle Uzawa
  + GC), `LadrunoContactHandler.cpp` (mortar pairing), `analyze_augmented.py` (held-load Uzawa).
- C2/C3 design + the lessons C4 inherits: `_adr41_c2_design.md` (per-node Domain state, the LOCAL-gap
  order-dependence), `_adr41_c3_design.md` (the DISPLACEMENT-not-position rule for relative motion — C4's `r_I`
  is exactly the friction-slip construction, but for all 3 components and no gT0), [[LEDGER_quirks]].
- Capstone [[48_ladruno_contact_capstone_adr]] mesh-tying row (📋 planned → C4); ADR-41 §scope (C4 = mesh-tying).
- Literature: Puso & Laursen mortar mesh-tying; Popp/Wall. Skills: `fem-mechanics-expert`, `opensees-expert`.
- Workflow: oracle `python3.12`; build `Ladruno_scripts\build.bat OpenSeesPy`; watch Zone-A; PR on `ladruno`,
  merge `--admin --squash`. See [[ladruno-adr41-mortar-c2]], [[ladruno-oracle-python]].
