---
wp: LEGACY
title: "LadrunoBrick -b is a body force per unit VOLUME, not per unit mass — it is never multiplied by rho"
legacy_seq: 357
---
## `LadrunoBrick -b` is a body force **per unit VOLUME**, not per unit mass — it is never multiplied by rho

**Found 2026-09-05, ADR-90 WP-A2.**

`LadrunoBrick::addLoad` (`LadrunoBrick.cpp:723-745`) does
`appliedB[i] += loadFactor * data(i) * b[i]` for `LOAD_TAG_SelfWeight` and
`appliedB[i] += loadFactor * b[i]` for `LOAD_TAG_BrickSelfWeight`, and the residual then takes
`-= dvol * appliedB[p] * shp[3][j]`. **`rho` appears nowhere on that path** — it is read only by
the mass matrix and the inertia terms.

So gravity is

```python
ops.element('LadrunoBrick', e, *conn, mat, '-b', 0.0, 0.0, -gamma, ...)   # gamma = rho*g
ops.eleLoad('-ele', *tags, '-type', '-selfWeight', 0.0, 0.0, 1.0)
```

Passing `-b 0 0 -9.81` (the acceleration, expecting the element to multiply by rho) silently gives
`1/rho` of the intended weight — on a `rho = 2.0` deck, **half** the soil weight, with every
analysis converging and reporting success. The cheap catch is the resultant identity
`sum R_z(base) == gamma * V`, which reads 4.4e-16 when the convention is right.
