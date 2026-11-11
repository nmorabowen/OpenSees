---
title: WP-230 — a low-pressure clamp mode in SAS-ME, for an A/B against separation (-sasLowP clamp)
project: Ladruno
tags:
  - sanisand
  - sas-me
  - tims
---

# WP-230 — `-sasLowP clamp`: a low-pressure clamp, the comparison mode to WP-152's separation

TIMs request F25 (`_tims_2d_model_requests_2026-11-11.md`). Opt-in, `LadrunoSANISAND` under SAS-ME
(`IntScheme 129`) only, exclusive with `-sasTensionCutoff`, byte-identical when off.

## 1. What is compared

WP-152 (`-sasTensionCutoff p_sep p_contact`): a point that loses its confinement SEPARATES. It goes to
σ = p_min·I, carries no shear, and α = α_in = 0.

The other rule in use for bounding-surface sand models clamps the pressure at a floor and keeps the
deviator. This WP adds that as a second mode of the same state machine, so the two can be run on the same
footing legs with ONE thing changed: whether a point that has lost its confinement carries shear.

## 2. Syntax

```tcl
... 129 $TanType $JacoType $TolF $TolR -sasHFloor 1 -sasReseatHyst 1 -sasSoftCap 0.5 \
    -sasLowP clamp $pFloor $pSep $pContact  <-sasLowPAlpha leave|project>  <-sasSepMaxP0 $p0max>
```

- `p_floor`: the pressure a clamped point is held at (model stress units). `p_sep`, `p_contact` and
  `-sasSepMaxP0` are `-sasTensionCutoff`'s, with the same meaning and defaults.
- No default for `p_floor`: it is a stress, and a dimensional default would be right in one unit
  system only. The A/B pair for a kPa deck with the cutoff at `-sasTensionCutoff 0.1 1.0` is
  `-sasLowP clamp 0.1 0.1 1.0`.
- Requires p_floor > 0, p_floor ≥ p_min, p_floor < p_contact, 0 ≤ p_sep < p_contact.
- Refused together with `-sasTensionCutoff` (they are alternatives), and wherever the cutoff is
  refused: `-Presidual ≠ 0`, `-pRe ≠ 0`, without `-sasHFloor > 0`, with `-implex`, any scheme but 129.
- `-sasLowPAlpha` without the clamp is refused (a flag nothing reads).

## 3. The state machine

Triggers, gates, closing law and re-contact condition are WP-152's, line for line:

- **Entry.** E1: SAS-ME's tension refusal (code 6, or code 3 because the committed p ≤ 0) at committed
  p0 ≤ p0max. E2: an accuracy or cost refusal (code 4 or 9) at committed p0 < p_sep. Both only under a
  non-compressing increment; a qualifying refusal a gate holds back is counted and still refuses.
- **Clamped** (g = tr ε − tr ε_entry, compression positive):

      σ = p_c(g)·I + s_k,     p_c(g) = p_floor + K(p_contact)·max(g, 0),     s_k = dev(σ_committed at entry)

  The deviator is FROZEN. It does not follow the strain: a clamped point is not integrated, it holds
  the shear it had when it lost its confinement. So it is bounded by what the last admissible state
  carried, and no stress is invented. The committed stress itself carries s_k, so there is no extra
  state.
- **Re-contact** at g ≥ g_c = (p_contact − p_floor)/K(p_contact): SAS-ME resumes from
  σ = p_re·I + s_k, p_re = p_c(g) ≥ p_contact. Volumetric only, as WP-152.
- Fabric kept; e follows the strain; the tangent is the elastic stiffness at the state (WP-152's
  declared Newton regularisation; while closing, the consistent bulk K(p_contact)).

### What "the same triggers" implies

The request reads "when a point's update would take p below the floor". The triggers here are the
REFUSALS of WP-152, not a test p < p_floor, because the request also asks for the same triggers so
that the A/B changes one thing. A point that integrates cleanly down to p = 0.05 kPa is not clamped
at 0.1 kPa; under the cutoff it would not be separated either.

### Which deviator is kept

The source rule keeps the deviator of the TRIAL stress and resets the pressure. SAS-ME has no trial
stress at a refusal: the update did not complete. The deviator kept here is the last COMMITTED one.
Two consequences, both measured below:

- **E1 (tension).** The committed pressure is already below the floor (the trajectory was on its way
  to p = 0). Clamping RAISES p, so the stress ratio falls: the clamped stress is inside the bounding
  surface. On the five oracle paths the largest stress ratio / bounding ratio a clamped point holds is
  0.34.
- **E2 (accuracy or cost, committed p0 < p_sep).** If p0 > p_floor the clamp LOWERS p and the stress
  ratio grows by p0/p_floor. With p_sep > p_floor a clamped stress can be outside the bounding
  surface. This is the inadmissible state the request describes. With p_sep = p_floor (the pair
  above) it cannot happen: every E2 entry has p0 < p_floor.

In the exact model the question does not arise at all: the trajectory reaches p → 0 on the yield
surface, with a deviator that vanishes with p. A kept deviator at p_floor is a property of a discrete
rule, of how far the update overshoots, and of which state the rule reads.

## 4. α and α_in of a clamped point

The source rule adds the stress change to α. That adds a stress to a stress ratio and is not
reproduced. Two dimensionally consistent treatments were built; `leave` is the default.

| | `-sasLowPAlpha leave` (default) | `-sasLowPAlpha project` |
|---|---|---|
| while clamped | α, α_in as committed | α := closest point of the yield cone around the clamped stress, α = s_k/p − √(2/3)·m·n, n = (s_k − pα)/‖s_k − pα‖; α_in := α |
| f at the clamped state | > 0: the stress is outside the yield cone | 0 |
| what it measures | `clampMaxFRatio` = f/(√(2/3)·m·p), in cone radii | 0 by construction |
| at re-contact | the same projection, once, and α_in := α | already on the cone; α_in := α |

**Why `leave`.** α moves only with plastic flow, and a clamped point has none. Leaving it reports the
inadmissibility instead of hiding it, which is what the A/B has to measure. At re-contact SAS-ME needs
a start with f ≤ 0, so α has to be brought to the cone there under either rule; the closest-point
projection is the smallest change in stress-ratio space that does it, and α_in := α marks the new
loading process (WP-152 does the same with α = α_in = 0).

**What `project` changes, on one element path** (`dev_open`: triaxial shear, isotropic opening past the
entry, shear while clamped, closing past re-contact; oracle and C++):

| | `leave` | `project` | separation (`-sasTensionCutoff 0.5 1.0`) |
|---|---|---|---|
| entry / re-contact step | 34 / 206 | 34 / 206 | 34 / 206 |
| stress while clamped | p_floor·I + s_k | the same, bit for bit | p_min·I |
| net work, C++ (oracle) | 0.08395 (0.08395) | 0.08395 (0.08395) | 0.08318 |
| q, p at the end of the path, C++ (kPa) | 0.4531, 178.47 | 0.4531, 178.47 | 0.0000, 177.39 |
| α at the end of the path | — | differs from `leave` by < 1e-10 | — |
| `clampMaxFRatio`, C++ (oracle) | 172.4 (172.4) | 0.0 (0.0) | n/a |
| `clampMaxBound`, C++ (oracle) | 0.340 (0.340) | 0.340 (0.340) | n/a |

Under `leave` the stale α sits 172 cone radii from the clamped stress; `project` removes that at no
cost in stress, work or end state (the two stress paths are bit-identical through the re-contact step, and
the end states agree to round-off). The stress itself is at 0.34 of the bounding stress ratio under
both. Against separation the same path does 0.9 % more work and re-contacts carrying its
kept deviator (q = 0.45 kPa at the end, against 0).

So the choice does not touch the mechanics of a clamped point. It decides only whether the stale α is
reported (`leave`) or removed (`project`). The cone is thin (m = 0.005 in the campaign set), so a
stress one stress-ratio unit away from α is over a hundred cone radii out: `clampMaxFRatio` is large
whenever α is stale, and says nothing by itself about the stress being admissible. The second counter
does.

## 5. Counters (`sasStats`, appended; 44 → 46)

- Entries, exits and the active flag are WP-152's columns, shared by both modes: `sas_sepEntriesTension`
  (E1), `sas_sepEntriesLowP` (E2), `sas_sepExits`, `sas_sepActive`, with `sas_sepLastCode`,
  `sas_sepMaxP0`, `sas_sepHeldHighP`, `sas_sepHeldCompressing`. Counted once, when the transition
  commits.
- `sas_clampMaxFRatio` (column 44): the largest f/(√(2/3)·m·p) a COMMITTED clamped state held.
- `sas_clampMaxBound` (column 45): the largest stress ratio over the bounding stress ratio,
  √(3/2)‖s‖/p ÷ (g(θ_s)·M_c·exp(−n_b ψ)), a committed clamped state held. Above 1 the stress is
  outside the bounding surface whatever α is: the part of the price no treatment of α removes.
- `sasOptions` 12 → 15: `sas_lowPMode` (0 separation, 1 clamp), `sas_lowPFloor`, `sas_lowPAlpha`
  (0 leave, 1 project). Wire tag 230.0.

## 6. Verification

Oracle first: `Ladruno_files/testbed/sanisand_lowp_clamp/clamp_oracle.py`, a driver around the
unmodified WP-134 oracle (the R1 testbed copy), the twin of WP-152's `tc_oracle.py`. Its
`mode="separation"` reproduces `tc_oracle.drive()` record for record on WP-152's three paths (740
records). Fixture: `tests/data/wp230_oracle_paths.json`.

Element paths (one stdBrick, zero free equations, eight identical Gauss points; campaign constants; R1;
p_sep 0.5, p_contact 1.0, p_floor 0.1 kPa; from the post-flip state of WP-152, p = 2.0 kPa), `leave`:

| path | entries E1 / E2 / exits (C++) | net work, C++ (oracle) | worst stress difference to the oracle, normal phase | `clampMaxFRatio`, C++ (oracle) | `clampMaxBound`, C++ (oracle) |
|---|---|---|---|---|---|
| `iso` | 1 / 0 / 1 | 0.003815 (0.003815) | 9.7e-16 | 0.0 (0.0) | 0.000 (0.000) |
| `te` | 0 / 1 / 1 | 0.004858 (0.004858) | 3.1e-03 | 162.7 (160.0) | 0.033 (0.034) |
| `cyc` | 4 / 0 / 4 | 0.001603 (0.001603) | 1.4e-15 | 0.0 (0.0) | 0.000 (0.000) |
| `dev_open` | 1 / 0 / 1 | 0.083951 (0.083954) | 5.2e-05 | 172.4 (172.4) | 0.340 (0.340) |
| `te_dev` | 0 / 1 / 1 | 0.020472 (0.020470) | 3.1e-03 | 162.7 (160.0) | 0.033 (0.034) |

- Entry and re-contact fall on the oracle's steps on all five paths, under both α rules (ten runs). On `te`
  and `te_dev` the C++ enters through E2 in the step where the oracle reaches tension, as it does under
  the cutoff.
- While clamped the deviator is the one committed before the entry (to 1e-13 kPa) and p follows the
  closing law (1e-9 relative; exactly p_floor while open).
- The stress difference column uses WP-152's measure (relative above 4 kPa, absolute below) and stays
  under its 5e-3 bound.
- `iso` and `cyc` unload isotropically: there is no deviator to keep and the clamp holds p_floor·I.

E2, C++ only (a cost refusal, `-maxSubsteps 1`, under isochoric shear, from a K0-like state at
p0 = 0.382 kPa, q0 = 0.245 kPa, p_sep 0.5):

| p_floor (kPa) | stress after the entry: p, q (kPa) | q/p | `clampMaxBound` | `clampMaxFRatio` |
|---|---|---|---|---|
| 0.02 | 0.02, 0.245 | 12.2 | 5.73 | 2318 |
| 0.1 | 0.10, 0.245 | 2.45 | 1.15 | 360 |
| 0.5 | 0.50, 0.245 | 0.49 | 0.23 | 29 |

The deviator is kept in all three; what changes is how far the pressure drops. At p_floor 0.1 with p_sep
0.5 this one point already holds a stress 15 % outside the bounding surface. With p_floor = p_sep
(the last row) an E2 entry can only raise p.

Tests: `tests/test_ladruno_sanisand_lowp_clamp.py`, 44 tests. Among them, OFF is byte-identical: WP-152's three
element legs with `-sasTensionCutoff 0.5 1.0` (740 steps) and with R1 alone (25 steps, up to each leg's
first refusal) reproduce, bit for bit, what a build WITHOUT this WP gave: the stress, α, α_in, the 44
earlier `sasStats` columns and the 12 earlier `sasOptions` values of every step
(`tests/data/wp230_cutoff_byteid_baseline.json`, recorded by building the branch's base sources in this
worktree, then rebuilding with the change).

SANISAND suites on Windows (`ci_run_pytest.py`, 15 files + the new one): 250 passed, 3 skipped, 2 xfailed
before (4a468a246); 294 passed, 3 skipped, 2 xfailed after. Seven pinned widths in two existing test
files were extended (`sasOptions` 12 → 15, `sasStats` 44 → 46).

MUTATION_LINE

Not verified:

- No footing run. The Workbench runs the A/B.
- E2 has no oracle counterpart (the oracle integrates exactly and has no accuracy or cost failure),
  as in WP-152. It is tested on the C++ alone.
- A Newton column with free equations through a clamp entry and re-contact (WP-152 has that test for
  separation; the clamp's volumetric law is the same and its deviator is constant, but the test was
  not ported).
- Plane strain: the state machine is in the shared 3D kernel; the tests drive the 3D brick.
- A stage flip back to the elastic stage while clamped: the limitation WP-152 documents applies
  unchanged.
