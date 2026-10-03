---
wp: LEGACY
title: "ManzariDafalias / LadrunoSANISAND TanType defaults to 0 = the ELASTIC tangent in the parser — so a deck that emits only the positional parameters runs algorith…"
legacy_seq: 355
---
## `ManzariDafalias` / `LadrunoSANISAND` `TanType` defaults to **0 = the ELASTIC tangent** in the parser — so a deck that emits only the positional parameters runs `algorithm Newton` as de-facto MODIFIED Newton

**Found 2026-09-05, ADR-90 WP-A2. Cost: a whole campaign relaunch, and ~7x of wall time on the
boundary-value legs.**

`ManzariDafalias3D::getTangent()` returns `mCe` (elastic) for `TanType == 0`, `mCep` for `1`, and
`mCep_Consistent` otherwise (`ManzariDafalias3D.cpp:134-142`; same shape in
`ManzariDafaliasPlaneStrain.cpp:137` and `ManzariDafalias3DRO.cpp:119`). The **parser default is
0** — `oData[1] = 0` at `LadrunoSANISAND.cpp:117`, deliberately copied from
`OPS_ManzariDafaliasMaterial` so that a renamed deck behaves identically — while the **null and
parallel constructors default to 2** (`ManzariDafalias.cpp:365`, `:426`). The two disagree, and
the one a deck actually hits is the parser's.

So `nDMaterial LadrunoSANISAND $tag <18 params>` with no positional optionals hands every
`algorithm Newton` an ELASTIC tangent: a modified-Newton iteration with linear convergence,
dressed as full Newton.

- **Where it hides:** on a **zero-free-DOF material-point deck it costs nothing** — there are no
  equations to iterate. Every existing fork SANISAND deck is of that shape
  (`tests/test_ladruno_sanisand.py` passes `*_PARAMS` plus flags and no positional optionals), so
  the defect has never had the chance to show itself. It appears the moment the material is put
  into a boundary-value problem.
- **Measured (ADR-90 WP-A2, strip footing on `LadrunoBrick -formulation bbar`):** at the parser
  default the h0 = 0.25 leg advanced 11 steps to s/B = 1.6e-4 in 350 s, spending ~40-65
  state-determination passes per step. With `TanType 2` and a reachable convergence test the same
  deck reaches s/B = 2.0e-3 in 36 s at h0 = 1.0 — about **7x**.
- **Emit the positional block explicitly:** `... $Rho $IntScheme $TanType $JacoType $TolF $TolR`
  **before** any `-flag` (the ADR-86 parser rejects a positional after a flag, by design).
- `mCep_Consistent` is **unsymmetric** under a non-associated flow rule: pair `TanType 2` with
  `system Pardiso -matrixType 0`, `UmfPack`, or another unsymmetric solver.
- The tangent changes the ITERATION PATH only — the substepped stress update is untouched — so
  this is free accuracy-wise. Measured on the WP-A2 deck, `q` at the matched `s/B = 0.002`
  checkpoint moved **0.52 %** between the two configurations.

> **FIXED for `LadrunoSANISAND` in WP-86b (ADR-86b, PR pending); NOT fixed in vanilla, on purpose.**
> `OPS_LadrunoSANISAND`'s `oData[1]` now defaults to **2** (the consistent tangent), and the
> construction echo NAMES the tangent it will run (`TanType = 2 (consistent mCep_Consistent
> (unsymmetric))`), so the change cannot be silent either.
> **`OPS_ManzariDafaliasMaterial` keeps its own default of `0`** (`ManzariDafalias.cpp:93`,
> re-verified at source during WP-86b) — every existing vanilla deck and every golden file produced
> by one depends on it, and moving it would be a silent answer-change in exactly the way this entry
> complains about. So **the two parsers now disagree on this one positional slot, deliberately.**
> An emitter or a deck that names all five positionals is immune to both defaults and to any future
> move; do that rather than relying on either.
> Note the null and parallel constructors already defaulted to 2 (`:365`, `:426`), so the fork
> parser is now the one that AGREES with them and vanilla's parser is the outlier.
