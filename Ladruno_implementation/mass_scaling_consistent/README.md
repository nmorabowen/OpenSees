# Consistent / Olovsson selective mass scaling — Phase 0 oracle bundle

The v2 ("consistent") increment of the mass-scaling stack (lumped v1 = ADR 36,
`CentralDifferenceSMS`). This bundle is the **oracle-first de-risking** of the
formulation before any C++. See [[38_ladruno_consistent_mass_scaling_adr]].

## What the oracle proves (`oracle_olovsson_sms.py`, numpy-only)

Centroidal Olovsson scaling mass, per element:

```
M_bar_e = beta_e * [ diag(m_a) - m m^T / M_e ]        (per spatial direction)
beta_e  = (dtTarget / dt_e)^2 - 1     for sub-target elements (same factor as lumped v1)
```

Row sums of `M_bar_e` are zero ⇒ `M_bar_e · t_rigid = 0`: rigid-body translation gets
**no** added inertia, so the global low modes (f1) are preserved; the added inertia lands
only on the relative/deformation modes that govern the explicit step.

Run:

```
pythoncore-3.12-64\python.exe Ladruno_implementation\mass_scaling_consistent\oracle_olovsson_sms.py
```

### Results (2026-06-20)

| Case | lumped f1 drift | **Olovsson f1 drift** | lumped added | Olovsson added | PCG iters |
|---|---|---|---|---|---|
| A — single tiny interior elem (20× smaller), dtTarget 0.012 | −54.18% | **−0.135%** | 344.8% | 173.6% | 12 (resid 5e-16) |
| B — refined zone, 6 small elems, dtTarget 0.012 | −82.62% | **−0.212%** | 2830% | 1417% | 21 (resid 1e-11) |

- **C1 frequency preservation** — Olovsson f1 drift < 0.5% and > 5× better than lumped
  (here ~400×). This is the entire v2 selling point and substantiates ADR-36 M4
  ("consistent allows ~10× more aggressive scaling").
- **C2 less mass, same step** — Olovsson injects ~half the total mass of lumped for the
  same dtTarget and still raises the global stable step to dtTarget.
- **C3 cheap solve** — Jacobi(lumped-diagonal)-preconditioned CG converges in 12–21 iters
  to machine/1e-11 tolerance; a lumped-only "solve" mis-solves `M_tilde` by 50–1945%, so
  the cross-node coupling genuinely must be solved (the matrix-free PCG the C++ runs
  inside the leap-frog step — keeps `system Diagonal` as the preconditioner source).

The added-mass percentages look large because a 20×-smaller element needs ~400× mass to
reach the bulk step (inherent physics); the point is *where* that mass lands — lumped dumps
it on the diagonal (wrecks f1), Olovsson confines it to the relative mode (f1 intact).

## Why this matters for the C++

The lumped v1 path (`Node::setMass` diagonal + `system Diagonal`, trivial `a = M⁻¹r`)
**cannot** represent `M_bar_e` — it couples DOFs across nodes. The consistent path needs:

1. a per-element `M_bar_e` block sizing kernel (Phase 1), and
2. a matrix-free PCG `a = M̃⁻¹r` inside the integrator step using the lumped diagonal as
   preconditioner and the stored `M_bar_e` blocks (mapped via `FE_Element::getID`) for the
   coupling mat-vec (Phase 2).

This oracle is the test reference for both.
