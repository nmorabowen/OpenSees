---
wp: WP-232i
title: "MPM images on a fix face: A5's FLIP convention at a fixed face node (a = -v_n / dt) breaks the equivalence with the symmetric full domain; the face node's fixed components must project to 0"
date: 2026-11-16
---
### MPM images on a `fix` face: the face node's fixed components must project to 0, not take A5's `a = −vⁿ/Δt`

- **Bites:** any ghost-mirroring ("images") boundary treatment layered on the slice-1 grid conditions. Slaving the ghost to its mirror (`v_g = R v_image`) alone is not enough: the face-layer node lies ON the symmetry plane, and interface §6.5 defines a fixed component's FLIP increment from the one-sided projection, `a = (0 − p/m)/Δt`. In the body plus its reflection, which the images stand for, the reflection cancels that momentum exactly, so `a = 0`. With A5's convention the half model drifted from the full one by 2e-5 (positions) and 24 % (stresses) after 120 steps (WP-232i G8, g++ core check).
- **Fix (WP-232i):** on a B-spline grid, zero the fixed components of momentum, internal and external force on every `fix` face's face layer before the update (`LadrunoMPMCore.h`, `zeroFaceOdd`), and of the projected v̂, â in mode A (`LadrunoMPMImplicit.cpp`). G8 then holds to 1e-15. GIMP and linear keep A5 (their frozen values).
- **And its mass (critic A1, WP-232i):** the same node's mass is one-sided too (1/2^k of the body plus its reflections, k = the `fix` faces through it), so the active-node cutoff dropped in the half model nodes the full model keeps: 3–14 % stress difference in 2D, 9.6 % in a 3D octant at mode A's default cutoff 5e-2, 53 % in mode F at `massCutoff 0.05`. The cutoff now tests m × 2^k at all three sites (`imageCutoffMass`: mode F step 4, `activeNodes`, mode A's mass map); B-spline grids only.
- **Also:** an odd fold `S_image − S_ghost` evaluated from node coordinates `x0 + i·h` can come out −6e-16 at the face, because the two coordinates round differently. Measuring `r` from the nearer face and folding per axis (a tensor product of 1D folds) makes `S′ ≥ 0` exact in floating point (`axisSupport`, `foldedSupport`).
