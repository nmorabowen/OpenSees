---
title: WP-238 — a correlated void-ratio field for LadrunoSANISAND (per-element voidRatio + generator)
project: Ladruno
tags:
  - sanisand
  - tims
  - random-field
---

# WP-238 — a heterogeneous initial void ratio: per-element `voidRatio` and a field generator

Owner request, 2026-11-19: try a model with a Gaussian distribution on the material, as is done for
concrete, to get a better band on the TIMs footing.

## 1. What is perturbed, and why the void ratio

In concrete the noise goes on a strength, and a crack-band or nonlocal law fixes the width of what
then localises. For sand the quantity that varies in a pluviated bed is the density, so the field is put
on the **initial void ratio**, not on a constant. Through ψ = e − e_c(p) it moves the peak strength, the
dilatancy and the plastic stiffness together, as one sand at different densities, and the calibrated
constants stay as they are.

## 2. What it can and cannot do

- It is an **imperfection**. It lets a band start at the weakest path instead of from round-off, and
  breaks the symmetry a homogeneous half-model imposes. On the footing the peak is late (s/B 0.137
  against 0.091 measured) and the mechanism dives under the near-surface strip; this is the lever for
  both.
- It is **not a regularisation**. SANISAND here has no internal length, so a band is still one element
  wide and the post-peak branch still depends on the mesh. Read the result by the peak load and its
  settlement, not by the band's thickness.
- The answer is a **distribution**: several seeds per leg, reported as mean and scatter.
- The mesh must be the **full** footing: the field is not symmetric.

## 3. The field (`Ladruno_scripts/sanisand_void_ratio_field.py`)

    e(x) = mean + std · Z(x),     Z(x) = √(2/M) · Σ_j cos(k_j · x + φ_j)

Random Fourier modes: the wave vectors k_j are drawn once from the seed, from the spectral density of
the chosen correlation; the phases are uniform. Z has zero mean, unit variance and that correlation.

- **It is a function of position.** Two meshes sample the same realisation; refinement resolves it.
  White noise per element does not have this property: its variance over a fixed volume falls with
  refinement (tested as the negative control: block means of a 4× refined white-noise field keep less
  than half the standard deviation).
- **Correlation**, with the scale of fluctuation θ per axis:
  `gaussian` ρ = exp(−π Σ(r_i/θ_i)²) (smooth) or `markov` ρ = exp(−2 √Σ(r_i/θ_i)²) (rough).
  θ should span a few elements of the COARSEST mesh of the study, or that mesh cannot represent it.
- **Truncation**: values are clipped to [e_min, e_max] and the clipped fraction is reported.
- numpy only; 2D or 3D points; a CLI that reads `tag,x,y[,z]` and writes `tag,e`.

```python
import sanisand_void_ratio_field as vf
f = vf.field(centroids, mean=e0, std=s_e, theta=(theta_h, theta_v), seed=1, e_min=emin, e_max=emax)
vf.apply(ops, ele_tags, f.values)        # before the first analysis step
```

The values are sampled at the points given (element centroids), not averaged over the element. With
θ of a few elements the difference is the cell's curvature term (the test bounds it at 0.1 of a
standard deviation, rms, for θ = 3.2 h).

**The amplitude and θ are inputs.** This WP supplies no number for either: the standard deviation of
the void ratio and its horizontal and vertical scales in an air-pluviated bed have to come from the
specimen data of the test being modelled.

## 4. The material side: `voidRatio` per element

`setParameter -val $e -ele $eleTag voidRatio` sets the initial void ratio of every Gauss point of one
element, so a mesh keeps ONE material tag.

Vanilla `ManzariDafalias` already has a `voidRatio` parameter (id 8). What was wrong, and is fixed in
`LadrunoSANISAND` with no vanilla edit:

- **The untagged element form was a silent no-op.** The base registers the parameter only with the
  material tag in `argv[1]`. `setParameter -ele N voidRatio` hands the material one argument, the base
  returns −1, no parameter is created, and the command returns normally. A deck applying a field that
  way ran a homogeneous model. `LadrunoSANISAND::setParameter` now claims the untagged form, as it
  does for `stressCorrection` (ADR-92 P2-4, the same defect). The tagged forms still go to the base.
- **No validation.** NaN, 0 or a negative value was stored. Now refused.
- **Allowed at any time.** In the plastic stage a point with strain history carries α, fabric and
  stress integrated with its void ratio; changing e under them is a jump no path of the model takes.
  Now refused there. It stays allowed in the elastic stage (vanilla's meaning: the CURRENT void ratio
  becomes the value) and at a strain-free point in any stage.
- **The initial elastic operator was left stale.** G depends on e; at a strain-free point the operator
  is now rebuilt as `initialize()` builds it, and `mVoidRatio` is refreshed.

Result: a point given e before the first step is a point constructed with e, bit for bit.

## 5. Verification

Generator, `tests/test_sanisand_void_ratio_field.py`, 23 tests, numpy only (run on CPython 3.11 and 3.12):

- the value at a point does not depend on which other points are asked for, or in what order (1e-12);
- a 4× refined grid's block means reproduce the coarse values (rms difference under 0.1 of a standard
  deviation at θ = 3.2 h) and keep over 90 % of the standard deviation; per-element white noise with
  "the same seed" fails both (the negative control);
- over 4000 seeds: zero mean, unit variance, and the target correlation at seven lags within 0.06, for
  both kernels and for an anisotropic θ; the line integral of ρ is θ to 1e-3; skewness, kurtosis and
  the one-sigma fraction are Gaussian;
- seeds reproduce bit for bit; truncation bounds and clipped counts; bad input and an unreadable
  centroid file are refused; the CLI round trip.

Material, `tests/test_ladruno_sanisand_void_ratio_param.py`, 16 tests:

- a point set to e before the first step against a point CONSTRUCTED with e: stress and α of every
  step (5 elastic, 40 plastic with a reversal), bit for bit, by the untagged element form, the tagged
  element form and the `parameter` object; the same in plane strain (`quad`);
- three disconnected bricks under one material tag, one left alone and two set through `vf.apply`:
  each reproduces, bit for bit, the single-brick run constructed with its own void ratio;
- the initial elastic operator: a two-brick column with free equations, one Linear step, moves exactly
  as the column constructed with that void ratio;
- refused, with the run left bit-identical to the untouched one: 0, a negative value, NaN, inf, a
  point with strain history in the plastic stage (and the warning is printed), a tagged call naming
  another material;
- allowed: the plastic stage while the point is strain-free; the elastic stage with strain history,
  where `psi` moves by exactly the change of the current void ratio.

The gate fails without the change: on the unchanged build of `ladruno` (cfe478254), 7 of the first 15
of these tests fail, the untagged form and everything built on it among them.

SANISAND suites on Windows (`ci_run_pytest.py`, the 15 files, plus the two new ones after):
250 passed, 3 skipped, 2 xfailed before (cfe478254); 289 passed, 3 skipped, 2 xfailed after.

Not verified:

- No footing run and no boundary-value problem with a field. Whether the field moves the peak, and by
  how much, is the Workbench's to measure.
- Other elements: `stdBrick` and `quad` are tested. Elements forward `setParameter` arguments to their
  materials themselves; one that does not forward them leaves the call a no-op again, so check `psi` at
  one point of the deck's own element type.
- Parallel runs (`OpenSeesMP`): the void ratio crosses `sendSelf` in the base's block, but the per-element
  set on a partitioned domain was not run.
