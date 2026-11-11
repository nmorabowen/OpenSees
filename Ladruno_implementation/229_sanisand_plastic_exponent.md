---
title: WP-229 — a separate pressure exponent for SANISAND's plastic modulus (-sasPlasticExp)
project: Ladruno
tags:
  - sanisand
  - sas-me
  - tims
---

# WP-229 — `-sasPlasticExp nh`: a separate pressure exponent for the plastic modulus

TIMs request F26 (`_tims_2d_model_requests_2026-11-11.md`). Opt-in, `LadrunoSANISAND` under SAS-ME
(`IntScheme 129`) only, default 0.5 = DM04, bit-identical at 0.5.

## 1. What it is

DM04 gives the elastic shear modulus and the plastic modulus the same pressure dependence:

- elastic: G = G0 p_at (2.97 − e)²/(1 + e) (p/p_at)^½;
- plastic: h = b0 / ((α − α_in):n), with b0 = G0 h0 (1 − ch e) (p/p_at)^−½, and Kp = (2/3) p h b:n.

So Kp ∝ p^½ and the elastic part of the loading denominator X = Q:C:R ∝ p^½: at a given state in
stress-ratio space (ψ, direction, (α − α_in):n) the ratio Kp/X does not depend on the pressure. The
stress-ratio/strain curve is then pressure-similar except through ψ.

This WP gives b0 its own exponent:

    b0 = G0 h0 (1 − ch e) (p/p_at)^−nh          nh = 0.5: DM04

b0 is unchanged at p = p_at. With nh < 0.5 it is smaller below p_at and larger above it:
Kp/X ∝ (p/p_at)^(0.5 − nh). That is the pressure trend in the strain at peak that F26's element
calibration found no bounded constant set removes.

The elastic exponent is NOT touched (it stays 0.5). Making it configurable belongs to WP-156 (#898),
which sets it from G_max data; this option is the plastic half, and #898 absorbs it.

## 2. Syntax

```tcl
nDMaterial LadrunoSANISAND $tag <17 constants> $Den  129 $TanType $JacoType $TolF $TolR  ... -sasPlasticExp $nh
```

- 0 ≤ nh ≤ 1. 0.5 is DM04 and takes the original `sqrt` line, bit for bit; the flag absent is 0.5.
- Refused for any `IntScheme` but 129 (the SAS-ME options rule). Under the CPPM schemes the
  Jacobians hard-code the ½, so the option could not be honest there.
- It is the 13th value of the `sasOptions` response (`sas_plasticExp`) and crosses `sendSelf`/`recvSelf`
  (wire tag 229.0).

## 3. Where it acts

One line, in vanilla `ManzariDafalias::GetStateDependent` (the only place b0 is computed), behind
`mLadrunoSas.opt.plasticExp != 0.5`. Everything downstream reads b0 through h: the stage, the drift
correction and the continuum tangent of SAS-ME need no change.

Oracle first: `Ladruno_scripts/sanisand_reference/model.py` (`Options.plastic_exp`, default 0.5, the
`sqrt` branch at 0.5) and the R1 testbed copy `Ladruno_files/testbed/sanisand_reseat_r1/sanisand_r1`.

## 4. How nh meets R1 (`-sasHFloor`, `-sasReseatHyst`, `-sasSoftCap`)

All three act on h downstream of b0, so nothing had to change in them. What changes is how often they
bind.

- **`-sasHFloor c_A`**: h = b0 / max((α − α_in):n, c_A √(2/3) m). The floor is on the DENOMINATOR, in
  yield-cone radii; whether a stage is floored does not involve b0. So the two compose
  multiplicatively: the largest h the floor allows is b0/(c_A √(2/3) m), and it carries the factor
  (p/p_at)^(0.5 − nh). With nh = 0.25 that bound is 0.56× DM04's at 10 kPa, 0.32× at 1 kPa and 0.18×
  at 0.1 kPa, and 2.3× at 3000 kPa (p_at 100). The set of floored stages is the same; the stiffness
  they get is not.
- **`-sasReseatHyst c_rev`**: a threshold on (α − α_in):n in cone radii. Independent of b0 and of nh.
- **`-sasSoftCap κ`**: where b:n < 0, h ≤ (1 − κ) X / ((2/3) p |b:n|). The cap is built from the
  ELASTIC part X of the loading denominator, which nh does not touch, so it is an absolute bound: where
  it binds, the response no longer depends on nh. It binds when b0/((α − α_in):n) exceeds it, so with
  nh < 0.5 it binds LESS often below p_at (b0 smaller) and MORE often above (b0 larger). Without the
  cap, the softening side of the model (b:n < 0, Kp < 0) reaches H = Kp + X ≤ 0, SAS-ME's refusal
  code 5, sooner above p_at and later below it, by the same factor.
- **Net effect at footing pressures** (p ≪ p_at): nh < 0.5 softens every plastic stage, moves the
  regime away from both the floor's bound and the cap, and so makes R1 bind less. Measured counts are
  in §5.

## 5. Gate: DM04's 17 triaxial simulations

Matrix: WP-150's gate 0 (DM04 Figs. 5–9, Verdugo & Ishihara 1996; `gate0_toyoura_oracle.py` on
`wp/150-sanisand-regularization-memo`), DM04 Table 1 constants, p_at 100 kPa, to 25 % axial strain in
500 equal increments. Gate: with nh in 0.2–0.3, q_peak and q at 25 % within ±10 % of the nh = 0.5
result.

**The gate FAILS as stated, on the oracle and in C++, in q_peak for three tests.** q at 25 % passes
for all 17.

Ratios to the nh = 0.5 run of the same code, q_peak / q at 25 % (oracle: `gate_nh_oracle.py`, `uw_model`; C++: `gate_nh_cxx.py`, SAS-ME, TolR 1e-4, R1 off):

| test | oracle nh 0.5: q_peak / q_25% (kPa) | oracle nh 0.2 | oracle nh 0.25 | oracle nh 0.3 | C++ nh 0.25 | C++ against oracle at nh 0.25 |
|---|---|---|---|---|---|---|
| `F5_u_e0.735_p100` | 3608.1 / 3589.0 | 0.992 / 0.997 | 0.992 / 0.997 | 0.992 / 0.998 | 0.992 / 0.997 | 8.2e-06 |
| `F5_u_e0.735_p1000` | 3594.8 / 3586.5 | 0.996 / 0.998 | 0.996 / 0.999 | 0.996 / 0.999 | 0.996 / 0.999 | 4.9e-06 |
| `F5_u_e0.735_p2000` | 3588.2 / 3585.2 | 0.998 / 0.999 | 0.998 / 0.999 | 0.998 / 0.999 | 0.998 / 0.999 | 3.9e-06 |
| `F5_u_e0.735_p3000` | 3585.3 / 3584.4 | 0.999 / 0.999 | 0.999 / 0.999 | 0.999 / 0.999 | 0.999 / 0.999 | 6.4e-06 |
| `F6_u_e0.833_p100` | 1353.9 / 1353.9 | 0.995 / 0.995 | 0.995 / 0.995 | 0.996 / 0.996 | 0.995 / 0.995 | 4.1e-06 |
| `F6_u_e0.833_p1000` | 1356.4 / 1356.4 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | 3.4e-06 |
| `F6_u_e0.833_p2000` | 1357.0 / 1357.0 | 1.099 / 1.003 | 1.051 / 1.002 | 1.002 / 1.002 | 1.051 / 1.002 | 7.5e-06 |
| `F6_u_e0.833_p3000` | 1408.6 / 1357.2 | **1.411 / 1.004** | **1.341 / 1.003** | **1.271 / 1.003** | **1.341 / 1.003** | 9.1e-06 |
| `F7_u_e0.907_p100` | 196.4 / 196.4 | 0.994 / 0.994 | 0.995 / 0.995 | 0.996 / 0.996 | 0.995 / 0.995 | 5.1e-06 |
| `F7_u_e0.907_p1000` | 425.9 / 209.1 | **1.300 / 1.027** | **1.246 / 1.023** | **1.194 / 1.019** | **1.246 / 1.023** | 1.2e-05 |
| `F7_u_e0.907_p2000` | 737.1 / 213.7 | **1.395 / 1.038** | **1.324 / 1.033** | **1.254 / 1.027** | **1.324 / 1.033** | 1.2e-05 |
| `F8_d_p500_e0.96` | 1047.4 / 1047.4 | 0.998 / 0.998 | 0.998 / 0.998 | 0.998 / 0.998 | 0.998 / 0.998 | 2.0e-04 |
| `F8_d_p500_e0.886` | 1069.3 / 1069.3 | 0.997 / 0.997 | 0.997 / 0.997 | 0.998 / 0.998 | 0.997 / 0.997 | 2.1e-04 |
| `F8_d_p500_e0.81` | 1146.8 / 1096.8 | 1.008 / 0.996 | 1.007 / 0.996 | 1.005 / 0.997 | 1.006 / 0.996 | 4.0e-04 |
| `F9_d_p100_e0.996` | 210.7 / 210.7 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | 3.2e-04 |
| `F9_d_p100_e0.917` | 215.7 / 214.8 | 0.999 / 1.000 | 0.999 / 1.000 | 0.999 / 1.000 | 0.999 / 1.000 | 2.0e-04 |
| `F9_d_p100_e0.831` | 243.1 / 220.6 | 1.003 / 0.999 | 1.003 / 0.999 | 1.002 / 1.000 | 1.003 / 0.999 | 4.1e-04 |

Bold = outside ±10 %. Oracle: 9 of 51 (test, nh) pairs fail, all in q_peak; the largest change of q at 25 % is 3.8 %.

C++ in full (`out_gate_nh_cxx.md`):

- 272 runs (17 tests × nh absent / 0.5 / 0.2 / 0.25 / 0.3 / 0.35 / 0.4 / 0.45 × R1 off / on), all completed, none refused.
- `-sasPlasticExp 0.5` against the flag absent: bit-identical, every q of every increment and the end state by `float.hex`, 34 runs.
- C++ against the oracle at the same nh, the whole curve, max |Δq|/q_max: 4.1e-04 worst over 119 runs (SAS-ME's own tolerance; nh = 0.5 alone is the same size).
- The same nine (test, nh) pairs fail in C++: 9 of 51 with R1 off, 9 of 51 with R1 on, with the oracle's ratios to the third digit (table).
- R1 on these paths: the floor acts on 2–7 stages per run (of 1368–57182) and the soft cap on 0, at nh 0.5 and at 0.25 alike; R1 on against off moves q by at most 6.9e-04 of q_max at either exponent. Monotonic triaxial loading does not exercise R1, so §4 is the statement and these counts only show that nh does not wake R1 up here.
- Drained runs: the lateral strain is solved per increment; the largest lateral-stress residual accepted is 7.7e-03 kPa (the bracket closed on a jump of the integrator's output, which is of the order TolR·σ).

What fails and why:

- The three failures are the loose, undrained tests that start at 10–30 p_at and contract to a
  quasi-steady state: e 0.833 from 3000 kPa, e 0.907 from 1000 and 2000 kPa. Their q_peak is the
  early undrained peak, which is set by Kp at p ≫ p_at, where nh < 0.5 makes b0 larger: at 3000 kPa
  by 30^0.25 = 2.3 with nh 0.25.
- The other 14, including every drained test and every test whose peak is the end of the curve, pass:
  13 of them within 0.8 % at every nh in 0.2–0.3, and `e 0.833, p0 2000` (which grows an early peak
  as nh falls) at +0.2 / +5.1 / +9.9 % for nh 0.3 / 0.25 / 0.2.
- With the elastic exponent at 0.5, all 17 pass only for nh ≥ about 0.42 (nh 0.45: worst +6.5 %;
  nh 0.40: two tests at +12–13 %).
- The intake's "all 17 within ±10 % with the constants at DM04" is its cell XC, which has the ELASTIC
  exponent free as well (nh 0.206, ng 0.756). Re-run here with a patched elastic law
  (`probe_elastic_exponent.py`; a probe, not an option), XC does pass: 0 of 17 outside, the three
  tests at +3.5 to +4.9 %. XB's exponents (0.251, 0.619) with DM04's constants fail the same three
  at +12 to +17 %. So the pass at XC rests on the elastic exponent, which the intake itself asks to
  keep out of the fit (its Asked 2): a larger ng stiffens G at p ≫ p_at and offsets the larger b0.
  With nh alone, the gate does not hold.

What this means for use:

- The exponent is anchored at p_at. The footing lives at p ≪ p_at, where nh < 0.5 only softens; the
  Verdugo–Ishihara tests live at 1–30 p_at, where it stiffens. One exponent cannot serve the
  plane-strain data at 5–400 kPa and leave the undrained peaks at 1000–3000 kPa alone.
- Two ways to keep both, neither built here (an owner / adviser choice):
  1. accept the three undrained peaks as out of the calibrated range (they are the tests DM04 itself
     fits least well in its Figs. 6–7), and gate on q at 25 % and on the 14;
  2. apply nh only below a reference pressure, b0 = b0_DM04 · min(1, (p/p_ref)^(0.5 − nh)), which is
     continuous at p_ref and leaves DM04 untouched above it. The 17 tests start at p ≥ p_at, but
     the loose undrained ones fall below it, so this would still need the gate re-run.

## 6. Verification

- `tests/test_ladruno_sanisand_plastic_exp.py`: bit-identity at 0.5 (R1 off and on); C++ against the
  oracle at nh 0.5 and 0.25 from 10, 100 and 1000 kPa; the direction of the effect; a one-stage
  identity on the tangent (the b0 ratio (p/p_at)^(0.5 − nh), exponents 0, 0.25, 1); parser refusals;
  the datastore round trip.
- `tests/test_ladruno_sanisand_reseat_r1.py`: the 643-replay byte-identity baseline, recorded on a
  build without WP-151, still passes with the flag absent.
- Testbed: `Ladruno_files/testbed/sanisand_plastic_exp/` (`gate_nh_oracle.py`, `gate_nh_cxx.py`,
  `probe_elastic_exponent.py`, outputs).
- Mutation: MUTATION_LINE

Not verified: any footing run (the Workbench runs it); plane-strain element paths (the option is in
the shared 3D kernel, and the tests drive the 3D material). `-implex` is refused with IntScheme 129, so
there is no IMPL-EX case to cover.
