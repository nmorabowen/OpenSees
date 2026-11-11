---
title: TIMs 2D model — the fork's answer to F25 and F26 (2026-11-11)
project: Ladruno
tags:
  - tims
  - sanisand
  - sas-me
---

# TIMs 2D model — answer to F25 and F26

Answer to `_tims_2d_model_requests_2026-11-11.md`. Both items are built, opt-in, on `LadrunoSANISAND`
under SAS-ME (`IntScheme 129`), each on its own draft PR against `ladruno`. Nothing is merged.

| | WP | PR | branch | plan |
|---|---|---|---|---|
| F26, plastic-modulus exponent | WP-229 | #1071 | `wp/229-sanisand-plastic-exp` | `229_sanisand_plastic_exponent.md` |
| F25, low-pressure clamp | WP-230 | #1072 | `wp/230-sanisand-lowp-clamp` | `230_sanisand_lowp_clamp.md` |

The two PRs are independent (both from `ladruno`). A deck that needs both options needs both merged;
the second to merge has a small mechanical rebase (the SAS option block and its wire layout).

## What the Workbench decks should use

The verified stack is unchanged. The two options are added to it:

```tcl
# F26 (WP-229): the stack, plus the plastic-modulus exponent
nDMaterial LadrunoSANISAND $tag <17 constants> $Den  129 $TanType $JacoType $TolF $TolR \
    -Presidual 0 -sasHFloor 1 -sasReseatHyst 1 -sasSoftCap 0.5 -sasTensionCutoff 0.1 1.0 \
    -sasPlasticExp 0.25

# F25 (WP-230), leg A: separation, the baseline (unchanged)
    ... -sasHFloor 1 -sasReseatHyst 1 -sasSoftCap 0.5 -sasTensionCutoff 0.1 1.0

# F25 (WP-230), leg B: clamp. p_floor, then the SAME p_sep and p_contact as leg A
    ... -sasHFloor 1 -sasReseatHyst 1 -sasSoftCap 0.5 -sasLowP clamp 0.1 0.1 1.0
#   optional, default leave:                           -sasLowPAlpha leave|project
```

- `-sasPlasticExp nh`: 0 ≤ nh ≤ 1; 0.5 is DM04 and is bit-identical to the flag absent.
- `-sasLowP clamp p_floor p_sep p_contact`: three values, no defaults (they are stresses). It REPLACES
  `-sasTensionCutoff`; giving both is refused. `-sasSepMaxP0` works with either.
- Read back: `sasOptions` (`sas_plasticExp`; `sas_lowPMode`, `sas_lowPFloor`, `sas_lowPAlpha`).
- Census for the A/B, per Gauss point, `sasStats`: `sas_sepEntriesTension`, `sas_sepEntriesLowP`,
  `sas_sepExits`, `sas_sepActive` (the same columns in both legs), and for leg B `sas_clampMaxFRatio`
  and `sas_clampMaxBound`.

## F26 — a separate pressure exponent for the plastic modulus

Asked 1 (the option), 5 (byte-identical at 0.5): **done.** b0 = G0 h0 (1 − ch e) (p/p_at)^−nh. The flag at
0.5 against the flag absent is bit-identical over 34 runs × 500 increments; the fork's 643-replay
baseline, recorded before WP-151, still passes.

Asked 2 (keep the elastic exponent out of it): **done by construction.** The option touches b0 only.

Asked 3 (the 17 triaxial tests within ±10 % for nh in 0.2–0.3): **the gate fails as stated.**

- q at 25 % strain: all 17 within 3.8 %.
- q_peak: 14 of 17 within 10 % (13 within 0.8 %). Three fail, at every nh in the range, on the oracle
  and in C++ with the same ratios:

  | test (undrained) | nh 0.30 | nh 0.25 | nh 0.20 |
  |---|---|---|---|
  | e 0.833, p0 3000 kPa | ×1.27 | ×1.34 | ×1.41 |
  | e 0.907, p0 1000 kPa | ×1.19 | ×1.25 | ×1.30 |
  | e 0.907, p0 2000 kPa | ×1.25 | ×1.32 | ×1.40 |

- Why: b0 is anchored at p_at. For nh < 0.5 it is smaller below p_at and larger above (×2.3 at 3000 kPa
  for nh 0.25). These three tests are loose, start at 10–30 p_at, and their q_peak is the early
  undrained peak, which follows Kp at high pressure.
- With the elastic exponent at 0.5 the 17 all pass only for nh ≥ about 0.42.
- The intake's "all 17 within ±10 %" is its cell XC, which also has the elastic exponent free (nh 0.206,
  ng 0.756). Re-run here as a probe it does pass (worst +4.9 %). XB's exponents with DM04's constants
  fail the same three (+12 to +17 %). So that pass rests on the elastic exponent, which Asked 2 keeps
  out. This checkout does not reproduce "17 of 17" for nh alone.
- C++ against the oracle at the same nh, over the whole curve, 119 runs: within 4.1e-4 of q_max. The
  C++ is the oracle's model; the failure is the model form, not the integrator.

What to decide (not built):

1. Gate on q at 25 % and on the 14, and declare the three undrained peaks at 1000–3000 kPa outside the
   calibrated range of an nh fitted at 5–400 kPa; or
2. apply nh below a reference pressure only, b0 = b0_DM04 · min(1, (p/p_ref)^(0.5 − nh)), continuous at
   p_ref. The footing lives at p ≪ p_at, so this keeps what F26 is for and leaves DM04 untouched above
   p_ref. The loose undrained tests fall below p_at late in the test, so the gate would be re-run.

Asked 4 (nh and R1): R1 acts on h downstream of b0 and needed no change.

- `-sasHFloor`: the floor is on the denominator (α − α_in):n, so which stages are floored does not
  depend on nh. The bound it sets, b0/(c_A √(2/3) m), carries (p/p_at)^(0.5 − nh): for nh 0.25, 0.56×
  DM04's at 10 kPa, 0.32× at 1 kPa, 0.18× at 0.1 kPa.
- `-sasReseatHyst`: a threshold on (α − α_in):n. Independent of nh.
- `-sasSoftCap`: an absolute bound built from the elastic part X of the loading denominator. Where it
  binds, the response no longer depends on nh. For nh < 0.5 it binds less often below p_at and more
  often above.
- Not measured on a footing. On the 17 triaxial tests R1 acts on 2–7 stages per run at either
  exponent, so those paths say nothing about the interaction at the wall states R1 was built for.

## F25 — a low-pressure clamp, for an A/B against separation

Asked 1 (the mode), 2 (the same triggers), 6 (byte-identical when off): **done.**

- One state machine, two modes. The clamp reuses WP-152's triggers E1 and E2, their gates, the closing
  law and the re-contact condition; only the state a point goes to differs.
- A clamped point sits at p = p_floor + K(p_contact)·max(g, 0) and keeps the deviator of its last
  committed state. The deviator is frozen while clamped.
- Off: WP-152's three element legs with `-sasTensionCutoff` (740 steps) and with R1 alone (25 steps)
  reproduce, bit for bit, a build without the clamp: stress, α, α_in, census and options of every step.

Two things in the request needed a decision. Both are stated here because they shape what the A/B can
show.

**"When a point's update would take p below the floor" against "the same triggers".** The triggers are
WP-152's refusals, as Asked 2 requires. A point that integrates cleanly to p = 0.05 kPa is not clamped
at 0.1; under the cutoff it is not separated either.

**Which deviator is kept.** The source rule keeps the deviator of a trial stress. SAS-ME has no trial
stress at a refusal, so the mode keeps the last committed deviator. What follows:

- After a tension entry (E1) the committed p is already below the floor. Clamping raises p and the
  stress ratio falls. On five element paths the clamped stress is at most 0.34 of the bounding stress
  ratio: inside the bounding surface.
- After an accuracy or cost entry (E2) at committed p0 > p_floor the clamp lowers p and multiplies the
  stress ratio by p0/p_floor. One Gauss point from p0 = 0.38 kPa with q/p = 0.64:

  | p_floor | q/p after the entry | stress ratio / bounding ratio |
  |---|---|---|
  | 0.02 | 12.2 | 5.73 |
  | 0.1 | 2.45 | 1.15 |
  | 0.5 | 0.49 | 0.23 |

- **With the pair above (p_sep = p_floor = 0.1) every E2 entry has p0 < p_floor**, so the clamp only
  ever raises p, and the clamped points stay inside the bounding surface. The A/B then measures a
  small kept shear, not an inadmissible state. To put clamped points outside the bounding surface the
  legs need p_sep > p_floor (for example `-sasTensionCutoff 0.5 1.0` against `-sasLowP clamp 0.1 0.5 1.0`,
  WP-152's own starting values). Which pair is the fair comparison is the act's call.
- In the exact model the question does not arise: the trajectory reaches p → 0 on the yield surface
  with a deviator that vanishes with p. A kept deviator at the floor belongs to the discrete rule.

Asked 3 (α and α_in): **`leave` is the default; `project` is built for comparison.** The stress change is
not added to α.

- `leave`: α and α_in stay as committed while clamped. No plastic flow happens, so no internal
  variable moves, and the state is reported as it is: outside the yield cone.
- `project`: α moves to the closest point of the yield cone around the clamped stress; α_in := α.
- At re-contact, under both, α is projected to the cone (SAS-ME needs a start with f ≤ 0) and
  α_in := α.
- The other's effect on one element path (`dev_open`: shear, opening, shear while clamped, closing):
  the stress path is bit-identical through re-contact, the work is the same (0.08395), and the end
  states agree to round-off. Only the reported inadmissibility differs: 172 cone radii under `leave`,
  0 under `project`. Against separation the same path does 0.9 % more work and ends with
  q = 0.45 kPa against 0.

Asked 4 (counters): entries, exits and the active flag are WP-152's columns, shared by both modes, so
the two legs are read with the same code. Added:

- `sas_clampMaxFRatio`: the largest f/(√(2/3)·m·p) a committed clamped state held.
- `sas_clampMaxBound`: the largest stress ratio over the bounding stress ratio. The first is 100–200
  whenever α is merely stale, because the yield cone is thin (m = 0.005). The second says whether the
  stress itself is admissible: above 1 it is not, whatever α is. Read them together.

Asked 5 (oracle first): a driver around the unmodified WP-134 oracle,
`Ladruno_files/testbed/sanisand_lowp_clamp/clamp_oracle.py`. Its separation mode reproduces WP-152's
driver record for record (740 records). The C++ follows it on five paths under both α rules: the same
entry and re-contact steps, the kept deviator to 1e-13 kPa, the work to four digits. E2 has no oracle
counterpart, as in WP-152, and is tested on the C++ alone.

For the A/B the act will run: the work done by clamped points is s_k : Δε_dev per step, with s_k
constant per point between its entry and its re-contact, so it can be summed from the recorded stress
and strain without anything new from the material.

## Not done, not verified

- No footing run of either option. The Workbench runs them.
- F26: the triaxial gate fails as stated (above). An owner / adviser decision.
- F25: no Newton test with free equations through a clamp entry and re-contact (WP-152 has one for
  separation); no plane-strain element path (the state machine is in the shared 3D kernel).
- Zone-A (Ubuntu) has not run on either PR: drafts skip the build. The Windows suites are in each PR.
