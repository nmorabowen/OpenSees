---
title: ADR-97 P3 — mutation-gate record (the Hoek–Brown curvature term)
project: Ladruno
status: complete
owner: nmora
tags:
  - implementation
  - material
  - review
---

# ADR-97 P3 mutation gate

**Question the gate answers.** Do the P3 tests actually measure the
**curvature** of the Hoek–Brown surface, or would they pass for a map that
treated the surface as locally flat?

## The mutation

One line, in `ASDPlasticMaterial3D::hb_assemble`, on a scratch build:

```cpp
-            A.col(i) += (z[nshape + k] * dscale[k]) * contrib;
+            (void) contrib;  // ADR-97 P3 MUTATION: curvature term dropped
```

`A = I + Σ_k dλ_k · D3 · (d m̂_k/dy)` is the stress-row block of the Newton
Jacobian. `d m̂/dy` reduces to a single rank-one column update because `dm/dy`
has exactly one non-zero entry, at `(i,i)` — and **that entry is the curvature
of the meridian**: it is `−a(a−1)·mb_psi²/σci · arg^(a−2)`, which is identically
zero only for a *linear* surface. Deleting it leaves `A = I`, i.e. the Jacobian
of a map whose flow direction does not turn as the stress moves.

The mutation touches only the **Jacobian**. The residual is untouched, so the
converged point — and therefore the committed stress — must stay exact. That is
the shape of the P1 mutation, and it is what makes the survivors interpretable.

## Result — 11 killed, 21 survivors (`tests/test_adr97_p3_hoekbrown.py`)

### The tangent (what the gate is for)

| measurement | shipped | mutated |
|---|---|---|
| free-node FACE rig, axis aligned | **3.342e-09** | **the global Newton does not converge at all** (`analyze -> -3`) |
| free-node FACE rig, sheared | **3.869e-09** | **does not converge** |
| load-driven oedometric rig (degenerate edge) | **1.485e-10** (`rel_fro` 1.915e-10) | **8.309e-03** (`rel_fro` 9.985e-03) — **5.6e+07× worse** |

The two face rigs are the sharper kill: with a flat-surface Jacobian the
assembled tangent is wrong enough that the *element-level* Newton on a nearly
singular perfectly plastic state stops converging — the same outcome the P2
gate records for the shipped `Backward_Euler` tangents on that rig. The edge rig
still converges and gives a number, which is the one to quote: five decimal
orders of degradation.

### The local Newton (a second, independent kill)

An inexact Jacobian still converges — it just costs iterations, and the ADR's
`≤ 5` gate on the curved face is exactly what catches that:

| deck | shipped `cp_iterations` | mutated |
|---|---|---|
| face, no shear | 4 | **7** |
| face, sheared | 4 | **8** |
| edge `y1 == y2` | 3 | **7** |
| near-apex face | 5 | **7** |
| triaxial / shear / rotating paths | {0,3} / {4} / {0,3,4} | {0,5,6} / {8,9} / {0,4,5,6,7} |

### The committed stress — unchanged, as it must be

| deck | shipped rel. error vs the oracle | mutated |
|---|---|---|
| face, no shear | 1.448e-15 | **1.593e-15** |
| face, sheared | 1.735e-15 | **1.446e-15** |
| edge `y1 == y2` | 7.259e-16 | **5.142e-16** |
| edge `y2 == y3` | 1.268e-15 | **4.077e-16** |
| apex ×2 | 0.0 | **0.0** |
| near-apex face | 1.851e-13 | **1.877e-13** |

Worst committed `|f|` 1.5e-11. This is the control: it proves the mutation
landed on the **Jacobian only**, so every survivor below is a survivor for a
reason, not by luck.

### Survivors — 21, each provably out of reach

* the six region-return **value** assertions (above) — the residual is untouched;
* the tension plateau, the past-the-corner apex returns, the hydrostatic and
  near-degenerate NaN guards, the step-refinement sequence, and the potential-gap
  pin: all read only the **committed state**;
* the eight gate-6 refusals: decided at **parse time**, by compile-time family
  markers the mutation does not touch;
* `strict_convergence` byte-inertness and gate 3 (the back stress stays zero).

### Elsewhere in the tree — untouched

`test_adr97_p1_smooth.py`, `test_adr97_p2_principal.py`,
`test_adr97_p4_inertness.py`, `test_adr97_p6_failloud.py` and
`test_adr94_hlist_hb.py` on the **mutated** build: **86 passed, 0 failed.** The
mutation is confined to the Hoek–Brown tangent, which is where it was aimed, and
`Backward_Euler` stays byte-identical through it.

## Revert

The line was restored (`git checkout -- SRC/.../ASDPlasticMaterial3D.h`), the
tree rebuilt, and every number in [[reviews/adr97_p3_report]] re-measured on the
restored binary.
