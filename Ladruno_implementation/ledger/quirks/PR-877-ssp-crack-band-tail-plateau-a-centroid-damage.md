---
wp: PR-877
title: "ssp crack-band tail plateau: a CENTROID damage probe cannot see a bending/hourglass-mode crack — the undamaged element under the load becomes an elastic hourgl…"
pr: "#877"
date: 2026-10-10
---
### `ssp` crack-band tail plateau: a CENTROID damage probe cannot see a bending/hourglass-mode crack — the undamaged element under the load becomes an elastic hourglass hinge
- **Bites:** `LadrunoQuad -formulation ssp` on the Kormeling–Reinhardt notched beam (5–10 mm meshes) keeps a post-peak plateau at 0.28–0.36 kN where `std` and CDPM2 give ~0.1 kN, for LadrunoConcrete3D and ASDConcrete3D alike; it disappears on the 2.5 mm mesh. It looks like the 1 % Tier-A floor carrying force across the cracked band — lowering the floor to 1e-4 only moves the coarse tail from 0.32 to 0.27 kN.
- **Why:** measured with the `hourglass` response (`Σ d·f_stab / u_load`): 0.20 of the 0.32 kN is stabilization work, 95 % of it in the three top-row elements under the load point, which have ω = 0 at their centroid. The last uncracked ligament bends inside one element (the bilinear hourglass mode), the centroid sits on its neutral axis, so a centroid-sampled `(1 − ω)` never degrades `Kstab` there. The fully cracked band contributes < 0.01 kN.
- **Fix/status (2026-09-26, Ladruno C1):** sample ω also at 4 shadow Gauss-point material copies (committed, plus the secant loss where damage > 1 %), floor 1e-4. `-hgDamage centroid` / `-hgLegacy` reproduce the old behaviour. **LadrunoBrick's Tier-A still uses the centroid probe and floor 0.01** — expect the same hinge in brick crack bands.
- **Also:** ASDConcrete3D's `damage` response is `1 − y/q` (the cracking part only, its plastic damage is not in it): at the same state it reads 0.76 where the secant loss is 0.98 (K&R coarse, the element under the load), so a `(1 − damage)` stiffness scaler under-degrades ASD elements.
