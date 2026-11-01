---
wp: PR-877
title: "A LAGGED (committed-only) damage probe on a stabilization scaler releases the stored hourglass force as a JUMP at the next step start — step halving cannot shr…"
pr: "#877"
date: 2026-10-10
---
### A LAGGED (committed-only) damage probe on a stabilization scaler releases the stored hourglass force as a JUMP at the next step start — step halving cannot shrink it
- **Bites:** `LadrunoBrick -formulation ssp` with the C3a Gauss-point shadows evaluated only at commit: the L&W beam (LadrunoConcrete3D, coarse) stopped at first cracking (0.6 mm, 26 kN) with "no convergence" through dλ0/64, while the centroid-sampled run (same floor) and the pre-C3a element went on to 51–56 kN.
- **Why:** the stabilization force is secant, `s·Kstab·u` on the TOTAL displacement. When the committed shadows crack, `s` drops (e.g. 1 → 0.1) between steps, so the first iteration of the next step carries an unbalanced force `Δs·Kstab·u` that does not depend on the step size; halving only shrinks the NEW increment, not the jump.
- **Fix (C3a):** evaluate the shadows at every trial `update()` once they have damaged (and while s is above the floor), so ω — and the force release — is resolved inside the Newton of the step that causes it. **LadrunoQuad's C1 shadows are still committed-only** (K&R converged with them; a candidate follow-up if a quad model stalls at cracking). 2026-09-26.
