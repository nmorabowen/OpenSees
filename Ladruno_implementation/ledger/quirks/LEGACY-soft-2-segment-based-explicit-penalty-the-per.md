---
wp: LEGACY
title: "SOFT=2 (segment-based explicit penalty): the per-node Courant bound is necessary-NOT-sufficient — inter-node coupling in the assembled K_c raises ω_max·dt ~2×"
legacy_seq: 27
---
### SOFT=2 (segment-based explicit penalty): the per-node Courant bound is necessary-NOT-sufficient — inter-node coupling in the assembled K_c raises ω_max·dt ~2×
- **Bites:** sizing the ADR-39 B2 SOFT=2 segment-based penalty (`contact … -mortar -soft <SOFSCL>`).
  Each slave node I gets a per-node Courant penalty `k_soft,I = SOFSCL·4·m_eff,I/dt²` so that node's
  contact mode `ω_I·dt = 2√SOFSCL` (≤ 2 for SOFSCL ≤ 1, exactly like the B1 NTS rule). The instinct
  is therefore "SOFSCL up to ~1 is stable." It is NOT for the segment lane: the ASSEMBLED contact
  stiffness `K_c = Σ_I k_soft,I·B_Iᵀ B_I` couples the slave + master facet nodes that multiple `p_I`
  springs SHARE (via the dense D/M mortar matrices), so the max eigenvalue of `M⁻¹K_c` is larger than
  any single `ω_I²`. On the oracle's non-matched facet pair the coupled `ω_max·dt = 1.19` at SOFSCL=0.1
  (vs the per-node `2√0.1 = 0.63`) — a ~1.9× amplification — and a SOFSCL sweep crosses the
  central-difference limit `ω_max·dt = 2` near **SOFSCL ≈ 0.3**, not 1.0.
- **Why:** the B1 NTS gap operator `B = [n | −Nᵢ n]` couples one slave node to one segment's nodes; at
  a shared segment node the multiple pairs' springs add, but the coupling is mild. The B2 segment
  operator `B_I = [D,−M]/a_I` is DENSE (every slave node of a facet couples to every master node it
  overlaps), so the assembled Gram `Σ k_I B_Iᵀ B_I` concentrates much more on the shared modes. The
  per-node sizing bounds each node's OWN mode, not the coupled spectrum.
- **Workaround/status (2026-06-24):** keep the **default SOFSCL = 0.10** (comfortable margin — coupled
  `ω·dt ≈ 1.2`). The command surface warns when `-mortar -soft SOFSCL > 0.25` (the per-node `> 1`
  warning understates the segment lane); a denser / more-coupled mesh can push the amplification
  higher, so a user raising SOFSCL must verify a dt margin (or drop dt). Quantified by
  `proto_b2_soft2_segment.py` T4 (the assembled-`K_c` eigenbound). This is the segment-lane analogue
  of B1's per-pair-conservatism caveat ([[#402](https://github.com/nmorabowen/OpenSees/pull/402)]).
  Also inherited from B1 (not B2-specific): `ladrunoBuildNodalMass` reconstructs `m_eff` from the
  DIAGONAL-of-consistent element mass; under `CentralDifferenceLadruno -lump rowsum|hrz` the integrator
  inverts a different lumping, modestly eroding the SOFSCL margin — diagonal/HRZ is fine, rowsum drifts.
