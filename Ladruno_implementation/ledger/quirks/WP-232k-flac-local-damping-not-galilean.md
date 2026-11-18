---
wp: WP-232k
title: "FLAC local damping (mpm damping local) is not Galilean invariant: it brakes rigid motion"
date: 2026-11-18
---
### Local non-viscous damping acts on the absolute grid velocity

- **Bites:** any explicit MPM run that moves a body rigidly with local damping on. The damping
  force is -alpha |f| sign(v) per node component (interface §6.5): it uses the absolute
  velocity, not a velocity relative to the body. Under a uniform drift sign(v) never changes, so
  it cannot dissipate oscillations about the drift, and a body accelerated by a body force is
  braked to about (1 - alpha) of the intended acceleration.
- **Seen (WP-232k, K0 diagnosis):** a floating column driven by a uniform body acceleration with
  damping 0.8 reached a kinetic energy of 0.0045 against 0.115; kept on during the K0
  translation, damping stopped the column (drop +0.010 m against 2.5 m).
- **Practice:** switch local damping off for a rigid translation (`mpm damping local 0`) and back
  on to settle (amendment D21a).
