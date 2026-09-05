# ADR-92 P0 — the single-Gauss-point IMPL-EX oracle

**Question P0 answers, in one line:** *at the increment the TIMs campaign actually uses,
how wrong is an IMPL-EX SANISAND, does that wrongness converge at first order, and which of
the two extrapolation forms should the C++ implement?*

P0 is **not gated** on #792. It writes no C++ and touches no fork source. Its output is
`_adr92_p0_oracle_results.md` and a decision on **D1**.

---

## 1. What is being built

`sanisand_implex_oracle.py` — a numpy reimplementation of `ManzariDafalias`'s stress update
at ONE Gauss point, exact to the C++ (including the code's departures from Dafalias &
Manzari 2004), carrying three integrators:

| name | what it is | role |
|---|---|---|
| `reference` | the substepped return run at a tiny increment | the "exact" answer every error is measured against |
| `implicit` | the return at the test increment | the baseline the campaign runs today |
| `implex_A` | **D1 form** — extrapolate the plastic strain (`eps_p = eps − eps^e`, committed) | the ADR's recommendation |
| `implex_B` | textbook form — extrapolate `dGamma`, freeze `n`, `R`, `D` at state *n* | the alternative D1 must beat |

For both IMPL-EX forms the delivered tangent is `Ce(p_n)` — the elastic operator frozen at
the **committed** mean stress — and the companion return runs once per committed step.
**The stress update is INCREMENTAL**, `sigma~ = sigma_n + Ce(p_n):(d_eps - d_eps_p~)`, because
the code's elasticity is hypoelastic (`elastic_integrator` `:1008-1011`); a total-strain form
`Ce(p_n):(eps - eps_p~)` is wrong for this law and its error does not vanish as `dt -> 0`.

**G2's strain increment is a lower bound.** `1e-4 m / h` is the nominal element strain; the
Gauss point at the footing corner sees a strain concentration of one to two orders more.
Sweep the increment; do not report a single number as "the campaign's".

## 2. The material — the D-L cell's own constants

Retyped by the TIMs act from `oracle_sheet.yaml` and quoted in
`Session Notes/2026-09-05 — Relaxation and Memory.md` §2. **Provenance matters: these are
the constants the campaign's cluster legs ran, not a textbook sand.**

```
G0 = 264.32   nu = 0.3129   e_init = 0.6944   Mc = 1.3309   c = 0.71
lambda_c = 0.027   e0 = 0.83   xi = 0.45   P_atm = 101   m = 0.005
h0 = 1.3   ch = 0.968   nb = 3.5   A0 = 0.05   nd = 5.75
z_max = 12.5   cz = 1100   rho = 2.0
IntScheme 1   TanType 0   JacoType 1   TolF = TolR = 1e-10
-Pmin 0.0101   -Presidual 0.0   -honorTolR 0
```

A dense second density, `e_init = 0.60`, is GATE U's other leg and is a cheap second column.

## 3. Gates

**G0 — the oracle is not allowed to measure itself.** Before any IMPL-EX number is quoted,
the oracle's `implicit` path must reproduce the **fork binary's** `LadrunoSANISAND` on a
common strain path to `<= 1e-8` relative in every stress component and in `e`, `||alpha||`,
`||z||`. Drive the binary with a single-element or `nDMaterial`-level probe through the
installed build (`C:\Program Files\Ladruno\OpenSees\bin\opensees.pyd`, interpreter
`C:\Users\nmb\venv\opensees_env\Scripts\python.exe`), and open the probe with
`ops.ladrunoBuild()` so the build hash is on the record. **If G0 fails, P0 stops and reports
the discrepancy** — a mismatch here is a finding about the material, not a bug in the oracle.

**G1 — first order.** The IMPL-EX-minus-reference error must halve when the increment
halves. Report the measured exponent over at least four increments, not an assertion.

**G2 — the campaign's increment.** Report the error at the strain increment corresponding
to the act's `1e-4 m` platen step (`h = 1.0 m` coarse and `h = 0.5 m`), in stress, in `e`
and in the mobilised stress ratio.

**G3 — the corner.** A path driven toward `p' -> Pmin = 0.0101 kPa` — the free-surface ring
beside the footing edge, where the implicit solver loses its descent direction. Does
IMPL-EX stay bounded, and what is the error there?

**G4 — D1 vs D1-alternative.** `implex_A` and `implex_B` on the same three paths, same
increments. The decision goes to whichever is more accurate *and* better behaved at G3;
if they tie on accuracy, D1 stands on the zero-vanilla-footprint argument.

**G5 — the tangent identity.** The returned `Ce(p_n)` equals the numerical
`d sigma~ / d eps` to machine precision. Trivial by construction — measure it anyway, because
it is the property the C++ must not lose when it reads `G`, `K` off the trial stress.

## 4. Paths

| path | why |
|---|---|
| **T1** drained triaxial compression, `p0 = 100 kPa` | the act's Level-0 instrument; the calibrated regime |
| **T2** the same at `p0 = 5 kPa` and at `p0 = 1 kPa` | low confinement, where `Presidual` and `Pmin` bite |
| **T3** a path driven to the `Pmin` floor and held | G3, the corner |
| **T4** compression then reversal | the one-step lag of `alpha`, `z`, `alpha_in` under IMPL-EX |

## 5. The disclosure P0 also owes

G1's exponent **is** the regularization statement (ADR §2, scoping §C5): if the answer moves
with the step size at first order, an IMPL-EX leg is a regularized leg with `dt` as the
parameter. P0 writes the two or three sentences the act must print beside any WP1 verdict
taken off an IMPL-EX curve — in numbers from G1/G2, not in adjectives.

## 6. House rules

- `python3.12`, numpy. Not `python`.
- No fork source is edited by P0. Nothing is added to `tests/` at this phase.
- Every number in the results memo is re-derivable from a committed script and a seed.
