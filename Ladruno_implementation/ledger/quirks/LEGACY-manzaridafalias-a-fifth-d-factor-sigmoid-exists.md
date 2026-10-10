---
wp: LEGACY
title: "ManzariDafalias: a FIFTH D_factor sigmoid exists, is dimensionally worse, is coupled to m_Pmin — and is dead code"
legacy_seq: 341
---
## ManzariDafalias: a FIFTH `D_factor` sigmoid exists, is dimensionally worse, is coupled to `m_Pmin` — and is dead code

**Found 2026-08-27, ADR-86 PR-3.** ADR 86 §7.2 enumerates four `D_factor` sites and PR-2
non-dimensionalised all four (`:3147`, `:3976`, `:4455` in the analytical Jacobian; `:4958` in
`GetStateDependent`). There is a fifth, at `ManzariDafalias.cpp:3556-3563` inside `NewtonSol2`,
and it is a *different* sigmoid:

```cpp
if (p < 0.001 * m_P_atm) {
    double be    = 207232.6584 * 2.0 * m_Pmin;
    double temp1 = exp(20.72326584 - be*p);
    double D_factor = MacauleyIndex(D) / (1+temp1);
```

- `be*p` carries **stress²**, so it is dimensionally worse than the four PR-2 repaired.
- Its steepness is **proportional to `m_Pmin`**. `207232.6584 = 20.72326584/1e-4` and
  `20.72326584 = -ln(1e-9)`, so `be` reduces to `2*20.723*P_atm` **only when `m_Pmin` has its
  vanilla value `1e-4*P_atm`**. `LadrunoSANISAND`'s `-Pmin` default would make it 10x steeper.
- **It is unreachable.** `NewtonSol2` is called only from `NewtonIter3`, and `NewtonIter3` has
  **no caller anywhere in `SRC/`** (verified by grep; `SAniSandMS` carries its own unrelated
  pair). `NewtonIter2_negP` is dead the same way — its only reference is a commented-out line
  at `:2244`.

So PR-2's "all four sites" is **correct about behaviour and incomplete about source**. Recorded,
not fixed: fixing dead code changes nothing and risks getting a derivative wrong. If
`NewtonIter3` is ever wired up, this becomes a live `-Pmin`-coupled shape defect.
