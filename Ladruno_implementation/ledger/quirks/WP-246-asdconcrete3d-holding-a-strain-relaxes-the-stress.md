---
wp: WP-246
title: "ASDConcrete3D-spine return: holding a strain relaxes the stress under tension-compression (the committed effective stress sits outside its surface)"
date: 2026-11-22
---
### ASDConcrete3D-spine return: holding a strain relaxes the stress under tension–compression, and a ramp depends on its step count

**What bites.** `ASDConcrete3DMaterial::compute` (and the `LadrunoRCConcrete` kernel cloned from it
before WP-246) is not idempotent: re-evaluating the committed strain from the committed state
changes the stress whenever both cones are active and the tension backbone has softened with a
plastic share (`-Td` < 1). Each hold step advances `xt` again and the stress relaxes geometrically
toward a fixed point; a monotone ramp over-advances `xt` on every step, so 4 and 256 steps to the
same strain differ (`xt` by 13–31 % on a compression + shear path). In a fibre section it shows as
a fixed-strain force drift and a database round trip that disagrees with the live model, though
the wire is faithful.

**Why.** The return stores `σ̃ = (1 − dt_pl)·ST + (1 − dc_pl)·SC`, two cones scaled by different
factors. `dt_pl` is chosen so that the tensile measure of `(1 − dt_pl)·ST` alone would equal
`q(x)`, but the Lubliner measure of the reassembled mixed-sign stress (I1, √J2, σmax) is not
homogeneous in `ST`, so `f(σ̃_n)/E + x_pl(x_n) > x_n`: the committed state is outside its own
surface, and the next call reads that overshoot as new loading.

**Fix (Ladruno only, WP-246).** `returnMap3D` measures the committed overshoot
`(x_n − x_pl,n) − f(σ̃_n)/E` and, when negative, takes it off the trial measure, so the threshold
advances by the increment of the measure. Gated by `rc_beamfiber_gpp.cpp` B11/B12. Upstream
`ASDConcrete3D` is unchanged (it is a vanilla file); a model using it directly still drifts.
**How to catch it elsewhere:** any return-mapping material should be tested for
`re-solve(committed strain, committed state) == committed stress`, and for step-count
convergence of a ramp; neither is caught by one-pass comparisons to an oracle that mirrors the
same algorithm.
