# ADR-43 L2 justification profile — is quadrature-parallel FEAST worth building?

**Question (ADR 43 §9 R0, umbrella ADR 45).** The shipped parallel FEAST path is
**L3-only** (#532): `eigen -feast … -rci` under `openseesmp` routes every contour
solve through a *distributed* `dmumps` (SYM=2, 2n block-real) across **all** ranks,
while the `dfeast_srci` outer loop runs replicated. **L2** (quadrature-parallel)
would instead split the world into `G` groups, run the L3 kernel in each group on a
*different* contour node concurrently, and `Allreduce` the projector. L2 is the
family's **largest novel-math surface** — it means abandoning `dfeast_srci` and
re-implementing FEAST's entire outer loop (Gauss–Legendre nodes/weights/Jacobian,
projector, Rayleigh–Ritz, refinement, trace convergence, m0 management, spurious
filtering). We build it **only if measured need justifies it** — never speculatively.

## Why the L3 strong-scaling curve settles it

FEAST's per-outer-iteration cost is `N_q` contour-node solves; conjugate symmetry
halves the distinct solves to `N_q_eff ≈ 4` (default `N_q = 8`). Each solve is a
factorize+solve of the 2n block system. Let `T_fac(P)` be the distributed
factorization time on `P` ranks, `S(P) = T_fac(1)/T_fac(P)` the speedup,
`E(P) = S(P)/P` the strong-scaling efficiency.

- **L3-only** runs the `N_q_eff` solves **sequentially**, each on all `P` ranks:
  `T_L3(P) ∝ N_q_eff · T_fac(1)/S(P)`.
- **L2** runs them **concurrently** in `G = N_q_eff` groups of `P/G` ranks:
  `T_L2(P) ∝ 1 · T_fac(1)/S(P/G)`.

So the achievable **L2-over-L3 speedup**, in the factorization-dominated limit, is

```
  speedup_L2/L3  ≈  N_q_eff · S(P/G) / S(P)
                 =  N_q_eff · (P/G)E(P/G) / (P·E(P))     with G = N_q_eff
                 =  E(P/G) / E(P).
```

**This is the whole decision, and it needs no L2 code.** `E(P/G)/E(P)` is read
straight off the L3 strong-scaling curve:

- **Perfect L3 scaling** (`E` flat) ⇒ ratio `= 1` ⇒ **L2 buys nothing** (concurrent
  small-group solves take exactly as long as sequential full-comm solves).
- **L3 saturated** at the budget (`E(P)` collapsed, `E(P/G)` still healthy) ⇒ ratio
  `≫ 1` ⇒ that wasted efficiency is exactly what L2 recovers, up to `≈ 1/E(P)` but
  bounded by `N_q_eff`.

Rule of thumb: `E(P/G)/E(P) ≳ 2×` on a **real** rank budget ⇒ L2 is worth the build;
`~1×` ⇒ it is not (spend ranks on L3, i.e. bigger per-node factorizations).

## Method

`solid_box.py` — a parameterized `ne³` `stdBrick` elastic box (fixed base, mass
density), a **solve-dominated 3D** generalized eigenproblem (3D fill-in makes the
distributed factorization the cost centre — the regime L3 parallelizes and L2 would
multiply). Replicated across ranks (the openseesmp assumption).

`l2_profile.py` — times the distributed `eigen -feast 0 fmax -rci` call at
`np = 1,2,4,8,…` on one fixed model (`np=1` = the serial PARDISO kernel, the L3=1
baseline), parses the `LADRUNO_FEAST_MPI` debug lines to confirm genuine
distribution + count factorizations, prints the strong-scaling table and the
`E(P/G)/E(P)` L2 ceiling.

```
# local (16 cores) — first estimate, np up to 8:
python Ladruno_implementation/feast_l2_profile/l2_profile.py 30 1,2,4,8 12
# large model on esmeralda (the real rank budget), reuse a coarse band:
L2_FMAX=<hz> mpi… python … l2_profile.py 70 1,2,4,8,16,32,64 12
```

Size guide (DOF ≈ `3·((ne+1)³−(ne+1)²)`): `ne=20 ≈ 26k`, `30 ≈ 86k`, `40 ≈ 206k`,
`55 ≈ 527k`, `70 ≈ 1.06M`. `L2_FMAX` skips the ARPACK band calibration (the lowest
modes are ~mesh-converged, so a coarse-mesh `fmax` is reusable and avoids a costly
ARPACK solve on a 1M-DOF mesh).

## Findings

_(appended after the runs)_
