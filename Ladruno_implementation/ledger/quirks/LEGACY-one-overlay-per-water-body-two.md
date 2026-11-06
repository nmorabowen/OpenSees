---
wp: LEGACY
title: "One overlay per water body: two LadrunoPorousOverlays must not share elements — and a shared water table needs ONE overlay, not one per soil layer"
legacy_seq: 183
---
### One overlay per water body: two `LadrunoPorousOverlay`s must not share elements — and a shared water table needs ONE overlay, not one per soil layer
- **Bites:** modeling a layered deposit as one overlay per layer disconnects the fluid: each overlay owns an independent p-field with its own drained set — no cross-layer flow, wrong consolidation. Conversely two overlays CLAIMING the same element double-count the fluid.
- **Why:** the p-field lives per-overlay (own CSR system); continuity exists only inside one region. Layered properties belong to `-layer` blocks INSIDE one overlay.
- **Workaround/status:** element overlap across overlays is a snapshot FATAL ⟨A-13⟩ (P1 battery-gated); the one-overlay-per-water-body rule is a modeling discipline the P4 guide owns. 2026-07-14.
