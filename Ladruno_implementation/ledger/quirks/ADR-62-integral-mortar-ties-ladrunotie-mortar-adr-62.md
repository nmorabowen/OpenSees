---
wp: ADR-62
title: "Integral-mortar ties (LadrunoTie -mortar, ADR-62 P2): global D⁻¹ pre-inversion DODGES the \"needs handler chain support\" wall — but P is dense (one big handler…"
legacy_seq: 139
---
## Integral-mortar ties (LadrunoTie `-mortar`, ADR-62 P2): global D⁻¹ pre-inversion DODGES the "needs handler chain support" wall — but P is dense (one big handler group)

The ADR-62 plan assumed integral-mortar ties would "couple slave nodes ⇒ MP-chains ⇒ need the
projection handler's deferred chain support." **That premise is avoidable.** The mortar constraint
`D u_s = M u_m` (D = slave-interface consistent mass, `D_IJ=∫N_I^s N_J^s dΓ`; M = `∫N_I^s φ_K^m dΓ`)
is condensed **once at model-build** by pre-inverting D over the WHOLE interface: `u_s = P u_m`,
`P = D⁻¹M`. Each row of P then ties a slave to **master nodes only** (no slave appears as a retainer),
so it emits as an ordinary `EQ_Constraint` the SHIPPED `LadrunoProjectionHandler` already accepts
(verified: it allows dense rows + master-only retainers + many slaves sharing masters; only a DOF that
is both retained-master and constrained-slave, or constrained twice, is a "chain"/"double"). **⇒ P2
needed NO handler change and NO kernel change** — it reuses `LadrunoMortarKernel::integratePair` verbatim
and only adds a setup-time generator + one `Matrix::Solve` (DGESV). Confirmed in `proto_p2_mortar_tie.py`
(P·1=1, linear-completeness patch, master-only rows) + `tests/test_ladrunoTie_mortar.py` (genuinely
non-matching solid patch test, 6/6).

CONSEQUENCES / gotchas for the standard-basis condensation:
- **P is DENSE.** `D⁻¹` of a sparse SPD mass is full, so every slave couples to every master in the
  connected interface ⇒ the handler builds ONE large group (all interface DOFs) and factorizes
  `(LᵀML)` over all master DOFs once per `domainChanged`. Fine for typical tie interfaces; for a HUGE
  interface this is the cost a **dual/biorthogonal basis** would remove (diagonal D ⇒ sparse P ⇒ small
  local groups) — that's the deferred P2.1 optimization. (Row-sum LUMPING D is NOT a shortcut: it keeps
  partition-of-unity but BREAKS linear completeness ⇒ fails the constant-stress patch.)
- **DGESV only flags an EXACT zero pivot**, not near-singular/ill-conditioned D. So the generator must
  guard coverage BEFORE the solve: compute each slave node's FULL tributary area `fullCov[I]=∫N_I^s` over
  the whole slave surface via a **self-clip** (`integratePair(npsS,Xs,npsS,Xs,...)` — a facet clipped
  against itself = the full facet, reusing the kernel), then refuse if the master-overlapped `cover[I]` is
  `< (1−1e-3)·fullCov[I]` (the slave protrudes past the master — a partial/extrapolated bond). A SINGLE
  slave facet half-overlapping the master does NOT give any node `cover≈0` (its shape fn spans the
  overlap), so the cover≈0 test alone misses protrusion — the cover/fullCov RATIO is what catches it.
  Belt-and-braces: a **post-solve partition-of-unity check** (`|Σ_k P_Ik − 1| < 1e-6`) catches any
  ill-conditioned solve that slipped through (P·1=1 is algebraically exact, so drift ⇒ bad D).
- **Reference coords, not trial.** Feed `integratePair` the as-built `getCrds()` (NOT X+u): a mesh-tie
  freezes the bond at the reference config (same as P1).
- **refDir (mortar normal orientation)** defaults to the average MASTER facet normal; it only orients n
  (used for the gap g̃ — magnitude only — and the aux plane), so its sign is irrelevant to D/M/P. If the
  master normals cancel (folded/curved surface) the generator refuses and asks for `-outward ox oy oz`.
