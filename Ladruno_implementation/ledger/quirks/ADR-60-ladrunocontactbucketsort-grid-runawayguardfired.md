---
wp: ADR-60
title: "LadrunoContactBucketSort::Grid::runawayGuardFired() is NOT a clean \"node ran away\" signal (ADR-60 R8)"
legacy_seq: 140
---
## `LadrunoContactBucketSort::Grid::runawayGuardFired()` is NOT a clean "node ran away" signal (ADR-60 R8)

The broad-phase grid's runaway guard clamps the centroid bbox to the `[clipPct, 100−clipPct]`
percentiles and sets `guardFired_` whenever that clamp **moves a bound**. With the shipped `clipPct=1.0`
that is the 1/99 percentile, so for ANY mesh with **>100 segment-centroids** the tails are clipped *by
design* and `guardFired_` is true — it does not mean a node diverged. So do **not** auto-surface
`runawayGuardFired()` as a warning (it would fire on every normal large model). ADR-60 R8 deliberately
leaves it a debug-only accessor. The real instability safety on the finite-sliding re-emit deformed feed
is `clipPct=0` (clip disabled so a genuinely-diverging node can't collapse the grid and silently drop
pairs) — a *behavior*, not a warning. Note also: with `clipPct=0` the guard never even computes
(`clip()` early-returns before the `guardFired_` test), so on the re-emit feed it is always false anyway.
