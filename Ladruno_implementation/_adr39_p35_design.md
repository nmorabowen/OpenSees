# ADR-39 P3.5 — implementation plan: implicit frictional Newton (consistent tangent)

> Makes the shipped explicit Coulomb friction (P3, #360) work under IMPLICIT integrators
> (static LoadControl, Newmark/HHT) by assembling the CONSISTENT friction tangent so
> Newton converges. P3 is FORCE-only: under an implicit solve a free tangential DOF is
> SINGULAR without a tangential stiffness (observed in P3 testing: `FullGeneral U(0,0)=0`).
> Grounded in two oracles: `proto_p3_friction.py` (per-traction blocks, 6/6) +
> `proto_p35_implicit_tangent.py` (the 3D ASSEMBLED tangent vs FD, 4/4). Parent ADR
> `39_..._adr.md`; driver `_adr39_loop_state.md` (iter 22). P3.5 = the friction TANGENT →
> critical junction (non-symmetry + active-set) → adversarial design gate before C++.

## Scope  — REVISED by the design gate (the default FLIPS to symmetric; see folds at end)

- **v1 P3.5 = the friction tangent in `addKtToTang`** (implicit path only). The explicit leg is
  UNCHANGED: CDL's `formEleTangent` calls only `addMtoTang` (no-op for contact) → the friction
  tangent is never assembled under explicit ⇒ **P3 explicit results are byte-identical** (gate-
  verified: `CentralDifferenceLadruno.cpp:222` calls addMtoTang only). The RESIDUAL is untouched.
- **DEFAULT = the SYMMETRIC tangent** `K_sym = Gᵀ D_TT P_t G` (the `d_TN⊗n` pressure-coupling
  column DROPPED). Symmetric ⇒ correct on ANY solver. `[GATE-Q2 BLOCKER]` the design's original
  "non-symmetric default" was a silent-wrong-answer foot-gun: symmetric SOEs
  (`ProfileSPDLinSOE.cpp:327`, `BandSPDLinSOE.cpp:269`, `SymSparse`) read only the upper triangle
  and SILENTLY drop the lower — so a non-symmetric tangent under the common `ProfileSPD` default
  factorizes a different system than the residual was linearized from (degraded convergence at
  best, wrong converged answer at worst), and the FE cannot detect the SOE type. So the default
  MUST be solver-safe. Converges superlinearly (the dropped column is the only inexactness).
- **OPT-IN `-consistanttan` = add the `d_TN⊗n` column** → the full CURRENT-N consistent tangent
  (true quadratic convergence). NON-SYMMETRIC ⇒ requires `FullGeneral`/`UmfPack`/`BandGen`; the
  PARSER emits an unconditional one-line WARNING stating that requirement (the one place that runs
  once and can warn; it cannot see the SOE either). Oracle-validated (`proto_p35_implicit_tangent.py`).
- **RESIDUAL ALWAYS uses CURRENT-N** (the physically correct equilibrium; `N=kn⟨−gN⟩₊`).
  `[GATE-Q4 resolved]` the tangent variants differ ONLY in which derivative terms they assemble;
  BOTH solve the SAME current-N root. Do NOT switch the residual's cap to committed-N (that would
  solve a physically wrong equilibrium — the cap would lag the true normal force). The symmetric
  default is thus a deliberately-approximate (modified-Newton-style) tangent for the exact
  residual: converges to the correct root, linearly/superlinearly not quadratically. Verified by
  reviewer-1 FD: same root, +1 iteration.
- **IMPL-EX (symmetric secant) = NOT in P3.5** — superseded; the symmetric consistent-minus-d_TN
  tangent is the robust default. Stays oracle-only.
- **DEFERRED to P3.5b** `[GATE-Q3]`: the per-step active-set FREEZE / chatter detector. The
  symmetric default already removes the chatter's main driver (the `d_TN` normal-coupling feedback
  that flips N within a step); the documented v1 mitigation for residual-level stick↔slip chatter
  is `algorithm NewtonLineSearch`. A flip-count freeze (per-pair, reset at commit) lands in P3.5b
  if a test demonstrates a limit cycle. Also deferred: the non-sym tangent's robustness under
  violent within-step N-swings (use `-symtan`/default there).

## The assembled tangent (oracle `proto_p35_implicit_tangent.py`, 4/4)

Per active pair, with `G = [I | −N_1 I | … | −N_nps I]` (3×ndof rel-disp operator), the
tangent-plane projector `P_t = I − n⊗n`, slip `gT = P_t·G·u`, normal gap `gN = n·G·u`:

```
r_fric = −Gᵀ tT                                  (P3 residual; slave −tT, master_i +N_i tT)
K_fric = Gᵀ [ D_TT·P_t + d_TN⊗n ] G              (K = −∂r/∂q ; ASSEMBLED ndof×ndof)
   STICK (‖tT*‖ ≤ μN):  D_TT = kt·I_T ,   d_TN = 0           → kt·Gᵀ P_t G   (SYMMETRIC)
   SLIP  (‖tT*‖ > μN):  D_TT = (μN·kt/‖tT*‖)(I_T − n̂⊗n̂),
                        d_TN = −μ·kn·n̂  (current-N)  |  0  (committed-N / `-symtan`)
```

`n̂` = the slip direction `tT*/‖tT*‖` (in the tangent plane). The SLAVE 3×3 block
`K_ss = D_TT·P_t + d_TN⊗n` carries the non-symmetry; the full `K_fric = Gᵀ K_ss G` scatters
it to the master DOFs (master block = `N_i N_j K_ss`, off-diagonal `−N_i K_ss`). Self-
equilibrium verified: `K_fric · u_rigid = 0` for a rigid tangential translation (ΣN_i=1).

Re-uses the EXISTING `n`, `N`, `N_i` (shape weights), `gTeff`, `gpT` already computed in the
friction residual block (one `segmentActive` + return map per evaluation). **Drops ∂n/∂u and
∂N_i/∂u** — the same small-slide main-term truncation as the shipped normal tangent (P2b-1
kn·BᵀB; the ∂n/∂u curvature block is the deferred P2b-2c). Residual stays exact; the dropped
terms are O(gₙ·κ) and the active-set main term carries the solve (P2b-2a proved this for normal).

## Where it lives (the FE)

- `LadrunoContactFE::addKtToTang(fact)` SEGMENT branch: after the existing `kn·BᵀB` normal
  block, when `mu>0` AND the pair is active (penetrating) AND the engine slot is reachable,
  add `fact · K_fric`. Recompute the return-map state (gTeff, tT*, slip/stick, n̂) exactly as
  the residual does — a PURE function of committed state + current trial disp (no new state).
- `addKiToTang` (initial-stiffness path): add the STICK friction tangent only (`kt·Gᵀ P_t G`,
  the initial tangential stiffness) — mirrors how the normal block adds `kn·BᵀB` there.
- Explicit (`addMtoTang`) UNCHANGED (no friction tangent). `getTangent` still routes through
  `formEleTangent` so the integrator picks the combination (CDL mass-only; Newmark c1·Kt; static Kt).
- **Flag plumbing (NOTE — the gate INVERTED the default, so the naming below flips):** the
  SHIPPED knob is **`-consistanttan`** → `bool consistentTan` on the `Contact`, threaded to the
  adapter. DEFAULT (`false`) = symmetric (drop `d_TN`); the flag opts INTO the non-symmetric
  consistent tangent (add `d_TN`) and the parser emits the FullGeneral/UmfPack warning. (The
  pre-flip text below uses the old `-symtan`/`symFrictionTan` names for the SAME mechanism with
  the opposite default — read it as `-consistanttan` with the default inverted.)

## Active-set within a Newton step (the chatter question)

- v1 P3.5 **re-projects the return map EVERY `getResidual`/`getTangent`** (consistent
  linearization — the tangent is the exact derivative of the residual it pairs with). The
  stick/slip decision is recomputed from the current trial disp each iteration. This is the
  standard consistent-tangent Newton and gives quadratic convergence near the solution.
- Chatter (a DOF flipping stick↔slip across iterations) is bounded by the penalty regularization
  (the cone is a smooth-ish radial return, not a hard Signorini). If a pathological case chatters,
  the documented mitigation is `algorithm NewtonLineSearch` (no code; standard OpenSees). An
  active-set FREEZE-per-step is NOT in v1 (it breaks quadratic convergence and the consistent-
  tangent contract); revisit only if a test demonstrates chatter.
- `N` for the cap: CURRENT penetration `N=kn⟨−gN⟩₊` for the default (matches the residual + the
  consistent tangent). `-symtan` uses committed N (frozen within the step) ⇒ symmetric.

## Tests — `tests/test_adr39_contact_p35_implicit.py`

1. **static frictional stick converges** — a slave on a fixed plane, normal load −P + a tangential
   load Q<μN, FREE tangential DOF, `system FullGeneral`, `Static`/Newton. WITHOUT the friction
   tangent this is singular (P3); WITH it, Newton converges to the elastic stick displacement
   `δ = Q/kt` (the gate — proves the tangent makes the implicit solve well-posed).
2. **quadratic convergence** — record the Newton residual norm per iteration on a frictional step;
   assert asymptotic quadratic rate (rₖ₊₁ ≲ C·rₖ²) for the consistent (current-N) tangent.
3. **static slip** — Q>μN with a tangential spring anchoring the slave (so equilibrium exists):
   the slave slides until the anchor balances μN; converges; friction force == μN.
4. **non-symmetry is real** — the slip step REQUIRES `FullGeneral`/`UmfPack`; document + a test
   that the consistent tangent solves the slip step (a symmetric-solver run is the `-symtan` leg).
5. **`-symtan` symmetric variant** — the same stick/slip problems converge with `system ProfileSPD`
   (symmetric) under `-symtan` (committed-N) — proves the symmetric knob works on a symmetric solver.
6. **Newmark implicit-dynamic friction** — a frictional transient under Newmark converges each step
   (the tangent is c1-scaled); cross-check the sliding against the explicit (CDL) P3 result.
7. **explicit byte-identity** — a CDL frictional run is byte-identical to P3 (the friction tangent
   is never assembled under explicit ⇒ no perturbation). Regression guard.

## Files (P3.5)
- `SRC/analysis/handler/LadrunoContactFE.{h,cpp}` — friction tangent in `addKtToTang`/`addKiToTang`
  (+ `symFrictionTan` member + ctor arg).
- `SRC/domain/contact/LadrunoContactDomain.{h,cpp}` — `Contact.symFrictionTan` + `addContact` arg.
- `SRC/interpreter/OpenSeesOutputCommands.cpp` — `contact … -symtan` flag (peek/un-read, like `-outward`).
- `SRC/analysis/handler/LadrunoContactHandler.cpp` — thread `ct.symFrictionTan` to the adapter.
- `tests/test_adr39_contact_p35_implicit.py`.
- Oracles (done): `contact_prototypes/proto_p3_friction.py` (6/6) + `proto_p35_implicit_tangent.py` (4/4).

## DESIGN GATE verdict + folds (2 source-grounded reviewers: tangent-mechanics + solver/active-set)

**tangent-mechanics → PASS** (no math wrong; doc-only). Independently re-derived `−∂r/∂q` and
FD-validated the FULL ndof×ndof assembly (not just the oracle's slave block) to 1.5e-9/kt using
the SHIPPED displacement-based slip. Confirmed: the `+Gᵀ[D_TT·P_t + d_TN⊗n]G` overall sign is
right (matches `+kn·BᵀB`); master scatter signs (`+N_iN_j`, `−N_i`) exact; `D_TT·P_t = D_TT` in
slip (idempotent — `n̂` in-plane, P_t·P_t=P_t — NOT a double-projection); c1 scaling correct
(rate-independent ⇒ all in `addKtToTang(c1)`, `Newmark.cpp:292`); `I_T ≡ P_t`.

**solver/active-set → SALVAGEABLE** (1 BLOCKER + 1 MAJOR, folded into Scope above):
- **Q1 sign/scatter → PASS.** **Q5 addKiToTang stick-only → PASS** (the stick tangent `kt·GᵀP_tG`
  is SPD and an upper bound on the slip stiffness ⇒ Modified/Initial-Newton contraction converges;
  adding the rank-deficient/non-sym slip block there would BREAK it; skipping friction → singular).
  **Explicit byte-identity → PASS** (no path assembles the friction tangent under CDL). **State
  reuse → PASS** (tangent must read committed `st.gpT`, NOT `gpTtrial`; add an FD-at-slip test).
- **Q2 → BLOCKER (folded):** flip the default to SYMMETRIC; non-sym is opt-in `-consistanttan` w/
  a parser warning. Grounded in the upper-triangle-only assemblers.
- **Q3 → MAJOR (folded/deferred):** per-step active-set freeze deferred to P3.5b; symmetric default
  removes the `d_TN` chatter driver; document `NewtonLineSearch`. (Reviewer wanted a freeze
  available in v1 — accepted risk: the symmetric default + penalty smoothing + line-search doc is
  the v1 mitigation; a flip-count freeze is the P3.5b fast-follow if a test shows a limit cycle.)
- **Q4 → resolved:** residual stays current-N (both reviewers); the symmetric default is an
  approximate tangent for the exact residual (same root, +1 iter).

**Status: design HARDENED + default flipped to solver-safe symmetric.** NEXT = code the symmetric
friction tangent in `addKtToTang` (+ `-consistanttan` opt-in adding `d_TN⊗n` + parser warning;
`addKiToTang` stick-only) → build → test (static stick converges = the gate; static slip; superlinear
default vs quadratic `-consistanttan` on FullGeneral; Newmark dynamic; explicit byte-identity; FD-at-
slip tangent==∂residual) → code gate → PR.
