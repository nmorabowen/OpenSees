---
wp: WP-128
title: "The ModifiedEuler error norm already HAS a 1 kPa floor, and its substep count is stability-limited — an F18(a)-style -errFloor below 2‖σ‖ is inert (WP-128)"
legacy_seq: 499
---
### The ModifiedEuler error norm already HAS a 1 kPa floor, and its substep count is stability-limited — an F18(a)-style `-errFloor` below `2‖σ‖` is inert (WP-128)
- **Bites:** "absolute below `‖σ‖ = 0.5`, `/(2‖σ‖)` above" is exactly `‖dσ₂−dσ₁‖ / max(2‖σ‖, 1 kPa)`. It is continuous, not an abrupt switch; the port with `σ_ref = 1` is bit-identical on 640/640 increments. `‖σ‖` includes the hydrostatic part (`≥ √3·p`), so a floor `σ_ref` only acts where `p ≲ σ_ref/3.5`. The substeps do not scale as `TolE^−½`: at constant p = 2 kPa, η/M^b 1.00, TolE 1e-6 / 1e-5 / 1e-4 gives 60 / 54 / 52 substeps. The scheme sits on its explicit stability limit (`∝ 1/√p`), so a looser norm buys ~nothing until it is loose enough to stop seeing the instability.
- **Workaround/status:** if `-errFloor` ships, its byte-identical default is **1**, not 0. Cost at low p needs a stiffly-stable integrator, not a norm. `128_sanisand_ring_trace.md` §5.
