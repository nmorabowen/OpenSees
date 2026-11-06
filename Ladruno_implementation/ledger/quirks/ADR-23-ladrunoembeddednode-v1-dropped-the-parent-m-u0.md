---
wp: ADR-23
title: "LadrunoEmbeddedNode v1 dropped the parent m_U0 offset-capture → absolute tie yanks on staged addition (FIXED, ADR 23 Phase 2c)"
legacy_seq: 69
---
### LadrunoEmbeddedNode v1 dropped the parent `m_U0` offset-capture → absolute tie yanks on staged addition (FIXED, ADR 23 Phase 2c)
- **The bug:** v1 computed every gap as a pure TRIAL-DISPLACEMENT difference
  (`g = u_c − Σ N_i u_host`, kernel `LadrunoEmbedded::computeGap`; likewise `g_p`, `g_r`), so
  the penalty enforced an ABSOLUTE tie `u_c = Σ N_i u_host`. An element added MID-ANALYSIS to a
  host that has already deformed (staged construction) activates with `g = −N·u_host ≠ 0` and the
  penalty **yanks the slave** by the full accumulated host displacement — a spurious force spike.
  The parent `ASDEmbeddedNodeElement` (and `equalDOF_Mixed`) already capture this offset
  (`m_U0` snapshot at `setDomain`, `getGlobalDisplacements()` returns `U − m_U0`); the fork's v1
  port silently dropped it.
- **The fix:** at `setDomain` capture each ACTIVE gap ONCE (`g0`/`gp0`/`gr0`, guarded by
  `g0Computed`) and drive ALL traction from the RELATIVE gap `(g − g0)`. Subtract the offset
  **inside** `computeGap`/`computeGapP`/`computeGapR` (NOT at each call site) so every consumer
  is covered in one place.
- **Force-free ≠ stress-free — the trap.** Zeroing only the penalty force is NOT enough in the
  D9 material mode: the gap also drives `matDir[d]->setTrialStrain(g·e_d)`. If the offset is an
  additive force correction, the material still sees the ABSOLUTE gap and is born PRE-STRAINED
  (a cohesive law partway up its backbone, a gap material already closed, bond pre-slipped) —
  force-corrected but NOT stress-free. Subtracting `g0` inside the gap (so the material's strain
  ORIGIN shifts) is what makes it genuinely stress-free. This is why "shift the canonical gap"
  beats "subtract at each consumer."
- **Default ON; no-op when undeformed.** Capture is ON by default (restores parent behavior);
  when the element is added at the undeformed state `g0 = 0` ⇒ byte-identical to the absolute
  tie, so the whole v1 battery is unaffected. `-absolute` (alias `-noInitGap`) opts out (legacy
  tie / a deliberate snap-to-host). `g0Computed` is serialized so `recvSelf` restores the
  captured offset instead of re-capturing. UR is linearized ⇒ `gr0` subtraction is exact for
  small inter-stage rotation, approximate for large. 2026-06-07.
