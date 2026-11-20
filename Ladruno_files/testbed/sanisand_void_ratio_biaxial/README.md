# WP-238 demo: plane-strain biaxial test with a void-ratio field

A boundary-value check of the WP-238 field generator: does a correlated void-ratio field
give a shear band, and what does it do to the peak. Not a validation and not calibrated:
`std = 0.02` and `theta = 0.02 m` are illustrative inputs.

## Model
- 0.10 x 0.20 m, `LadrunoQuad -formulation bbar -type PlaneStrain`, `LadrunoSANISAND`
  SAS-ME (IntScheme 129), mean e = 0.6944, sigma_3 = 100 kPa.
- Smooth platens: bottom u_y fixed (u_x at the middle node), top u_y tied to one master.
- Stage 0 elastic confinement, then displacement control on the top platen.
- `biaxial.py` checks first that the element forwards `voidRatio` (psi moves by exactly
  the change of e). It does on `LadrunoQuad` bbar.

## Run
```
set BIAX_BIN=<dir with opensees.pyd or opensees.so>
python -S -u biaxial.py --out out/seed1 --seed 1 --tantype 0
python plot_biaxial.py out            # bands.png, curves.png, summary.md
python plot_stall.py out/seed1        # last state of one run
```
`esmeralda/` holds the build and the 1-CPU `srun` launcher used for the runs below.

## Result (2026-11-19, esmeralda, build c8ac8b36c, elastic tangent `--tantype 0`)
`results/bands.png`, `results/curves.png`, `results/summary.md`.

- **Every run forms a band, and every run stops on the minimum step soon after.** The
  20x40 runs end between 3.4 % and 4.8 % axial strain (target 14 %). The softening branch
  is not traced on that mesh. The peaks of seeds 3 and 4 are the last converged point,
  still rising: lower bounds.
- **The field selects the band.** Six seeds give six patterns (single bands of both
  signs, V and X pairs). The homogeneous specimen gives one band from a top corner.
- **Peak axial stress, 20x40:** homogeneous 716 kPa at 3.2 %; six seeds 598 to 690 kPa
  (mean 632). One white-noise realisation: 653 kPa, inside the range of the correlated
  ones; one case does not separate the two.
- **The peak depends on the mesh, with or without a field.** Homogeneous: 868 kPa at
  8.5 % on 1x2 (no band), 846 kPa at 6.0 % on 10x20, 716 kPa at 3.2 % on 20x40. Seed 1:
  742 kPa at 5.5 % on 10x20, 604 kPa at 3.1 % on 20x40. The band is one element wide.
  The field is an imperfection, not a regularisation.

## Solver finding
With the SAS-ME tangent (`--tantype 1`) Newton converges linearly (about 50 iterations
per step even on two elements) and the 20x40 runs stall near 1 % strain; a band is
already forming there (`results/stall_seed1_sasme_tangent.png`, shear max/mean 1.44).
Raising `-maxSubsteps`, dropping the R1 options or line search do not help. The elastic
tangent goes about three times further. Whether the stall is the tangent or the loss of
uniqueness of a non-associated model was not separated.

## Post-peak (2026-11-20, esmeralda, same build, `--tantype 0 --tol-relax 10 --tension-cutoff 0.5 1.0`)

Why the runs above stop: at the stall exactly ONE Gauss point of 3200 refuses the SAS-ME update
(`maxSubsteps`, code 9, or the minimum substep), at an in-plane mean stress below 0.001 kPa
(homogeneous: element 590 gp 4 at (0.0475, 0.1475) m; seed 2: element 211 gp 4 at
(0.0525, 0.0525) m). The stress is read after a failed attempt and excludes the out-of-plane
component. A point in the band has lost its confinement, so neither a smaller step, a tenfold
substep cap (`--maxsub 20000`) nor a looser unbalance tolerance (`--tol-relax 10`, alone) moves the
stall by more than 0.001 % strain. The deck logs this diagnostic at every such stop.

With the material's tension cutoff (WP-152, `-sasTensionCutoff p_sep p_contact` = 0.5 1.0 kPa) the
runs go through: `results/postpeak_bands.png`, `postpeak_curves.png`, `postpeak_summary.md`.

| case | mesh | peak (kPa) | at eps_a (%) | end (kPa) | end eps_a (%) | points separated at the end | steps at 3x / 10x tolerance |
|---|---|---|---|---|---|---|---|
| homogeneous | 20x40 | 715.8 | 3.23 | 567.5 | 13.91 (wall budget) | 23 | 132 / 757 |
| seed 1 | 20x40 | 604.3 | 3.13 | 495.3 | 14.00 | 59 | 0 / 0 |
| seed 2 | 20x40 | 620.2 | 4.86 | 512.9 | 14.00 | 43 | 0 / 0 |
| seed 3 | 20x40 | 602.7 | 4.18 | 495.4 | 14.00 | 52 | 0 / 0 |
| seed 4 | 20x40 | 675.2 | 4.85 | 551.7 | 13.21 (wall budget) | 57 | 143 / 116 |
| seed 5 | 20x40 | 656.0 | 4.50 | 538.9 | 14.00 | 53 | 0 / 0 |
| seed 6 | 20x40 | 689.8 | 3.91 | 565.7 | 14.00 | 77 | 0 / 0 |
| white noise | 20x40 | 687.6 | 5.68 | 591.5 | 12.40 (wall budget) | 43 | 67 / 0 |
| seed 1 | 10x20 | 742.1 | 5.47 | 671.2 | 14.00 | 4 | 0 / 0 |

(Of 3200 Gauss points on 20x40, 800 on 10x20.)

- Every case forms one through-going band and softens; the six seeds lose 17 to 18 % of their peak
  by 14 %. Six of nine runs needed no relaxed step at all.
- The homogeneous, seed 4 and white-noise runs are NOT clean: they hit the 8 h wall budget and
  oscillate by tens of kPa after the peak; the first two accepted steps at up to 10x the unbalance
  tolerance (0.1 kN), the white-noise run at 3x. Read them as a trend.
- The post-peak slope is mesh dependent, as expected of a local softening model: seed 1 ends at
  495 kPa on 20x40 and 671 kPa on 10x20. Nothing here regularises the band width.
- The separated points sit in the band. The cutoff was designed for footing edges (WP-152); inside a
  shear band it has not been validated against anything, so the softening branch is a traced
  numerical path, not a verified one.

## Dynamic relaxation: tried, rejected for this model

`--solver dr` (LadrunoDynamicRelaxation, continuous platen ramp with holds) is in the deck but gives a
wrong answer with SANISAND: the damped oscillations register as loading reversals and re-seat
`alpha_in`, which stiffens the response. Homogeneous 10x20 at 1 % strain, static 407.1 kPa
(0 % of commits move `alpha_in`):

| DR setting | sigma_a at 1 % (kPa) | commits that move alpha_in |
|---|---|---|
| kinetic, 1e-7 m/step | 607 | 6.0 % |
| kinetic, 2e-8 m/step | 581 | 0.156 % |
| viscous zeta = 1, 2e-8 m/step, held | 420.6 | 0.001 % |

The viscous run still reads 261.3 against 293.6 kPa at 0.5 % (lag). `-reversalTol` / `-reversalRel`
are refused under IntScheme 129. The deck logs the `commitStats` census on both paths.
