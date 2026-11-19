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
