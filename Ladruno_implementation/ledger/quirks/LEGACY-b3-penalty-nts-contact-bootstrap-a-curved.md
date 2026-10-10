---
wp: LEGACY
title: "B3: penalty NTS contact bootstrap — a curved indenter starts contact at ONE point"
legacy_seq: 133
---
### B3: penalty NTS contact bootstrap — a curved indenter starts contact at ONE point
- **Bites:** driving a curved (sphere/cylinder) indenter into a half-space by force. The first contact is a
  single point; the not-yet-contacting indenter material/columns are unsupported ⇒ free-fall under load ⇒
  Newton diverges at step 1 (seen across LoadControl / DisplacementControl / weak-spring-stabilized free
  indenters AND one-shot full pre-penetration of a fixed rigid sphere). The shipped `block_on_block` test
  converges only because its interface is FLAT (all slaves engage at once) at a tiny 1e-8 pre-penetration.
- **What works for a Hertz patch:** a FIXED rigid sphere (slaves pinned at the sphere surface ⇒ no free
  body) with a MODEST approach δ + moderate penalty kn + a looser convergence tol (1e-10 stalls on the
  stiff patch), and/or ramping the indentation gently. Robust quantitative 3D Hertz remains sensitive — it
  motivates pairing B3 with displacement control or D1 within-step augmentation. Found by the B3 gate 2.
