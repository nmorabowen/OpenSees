---
wp: PR-877
title: "ASDConcrete3D's damage response omits its plastic damage — a (1 − damage) secant scaler under-degrades it (upstream material, note only)"
pr: "#877"
date: 2026-10-10
---
### `ASDConcrete3D`'s `damage` response omits its plastic damage — a `(1 − damage)` secant scaler under-degrades it (upstream material, note only)
- **Bites:** a stiffness/stabilization scaler built on `setResponse("damage")` (the LadrunoQuad / LadrunoBrick Tier-A) sees 0.76 on an ASDConcrete3D point whose secant energy loss `1 − σ·ε/(ε·C0·ε)` is 0.98 (K&R coarse, the element under the load).
- **Why:** `ASDConcrete3DMaterial::getMaxDamage()/getAvgDamage()` return `HardeningLaw::evaluateAt(x).crackingDamage()` = `d = 1 − y/q`, the cracking (stiffness) part only; the plastic part of the backbone (`dt_plastic`, `dc_plastic`) is not in it.
- **Workaround:** the fork's Tier-A shadows also take the secant loss at points that report damage > 1 %. Upstream file not changed. 2026-09-26.
