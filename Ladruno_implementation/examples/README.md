# Explicit Noh–Bathe integrators — usage notes

`ExplicitBathe` (Noh–Bathe two-sub-step explicit) and `ExplicitBatheLNVD`
(same scheme + FLAC local non-viscous damping for dynamic relaxation).

See [`explicit_bathe_example.py`](explicit_bathe_example.py) for a runnable
OpenSeesPy demo (free vibration + LNVD relaxation).

## The explicit recipe (required)

Explicit integration only pays off — and is only correct — with:

| ingredient | why |
|---|---|
| **lumped / nodal mass** | the RHS is just `M`; a consistent mass on a diagonal solver is silently wrong, on a sparse solver it's pointless |
| `system Diagonal` | trivial `M⁻¹`; never a banded/sparse solver |
| `algorithm Linear` | the scheme performs **exactly 2 solves per step**; Newton/iterative will error |
| `test … 1` | single pass |
| `dt < dt_cr` | conditionally stable |

## Commands

```
integrator ExplicitBathe      <p=0.54> <-cfl> <-cflAbort> <-verbose> <-divergence f>
integrator ExplicitBatheLNVD  <p=0.54> <alpha=0.8> <-cfl> <-cflAbort> <-verbose> <-divergence f>
```

- **`p`** — sub-step parameter `0<p<1`; `0.54` gives good high-frequency dissipation.
- **`alpha`** (LNVD) — FLAC local damping `0≤alpha<1`; classic default `0.8`. Damps
  toward static equilibrium and vanishes as `v→0` (no bias). Watch `|unbalance|→0`
  (printed with `-verbose`) to judge convergence.
- **`-cfl`** — estimate the critical step (per-element eigenproblem) and print it.
  Two numbers are reported: the **conservative central-difference limit** `2/ω`
  and the **Noh–Bathe limit ≈2×** that (Noh & Bathe, *Comput. Struct.* 129, 2013).
- **`-cflAbort`** — stop the run if `dt` exceeds the Noh–Bathe limit.
- **`-divergence f`** — circuit breaker: stop if the kinetic-energy proxy grows by
  factor `f` in one step (spurious-energy blow-up). NaN/Inf accelerations always abort.
- **`-verbose`** — per-step reporting (off by default; explicit runs are millions of steps).

## Pair with the energy balance

Run a `recorder EnergyBalance -file energy.txt -time` alongside. For a stable run
the residual column `RES = ULW − (KE+IE+DW)` stays ≈0 and `ERR%` stays small;
a steadily growing `ERR%` is the instability signature of too-large `dt`.

## Parallel note

Under `openseesmp` each rank builds its own integrator locally; the broker
registration (tags 61/63) covers `openseessp`/database-restart reconstruction.
The `-cfl` element loop is per-partition (no cross-rank reduction) — the reported
`dt_cr` is a local-partition estimate.
