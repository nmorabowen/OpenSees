---
wp: ADR-66
title: "Solid-shell patch tests on a 1-element-thick mesh: the interior-node patch MUST use a traction-consistent field (ADR-66 P5.1)"
legacy_seq: 147
---
## Solid-shell patch tests on a 1-element-thick mesh: the interior-node patch MUST use a traction-consistent field (ADR-66 P5.1)

Every node of a one-element-thick patch lies ON the free top/bottom faces. A full affine gradient
carries `sigma·e_z != 0` there, so with no applied face tractions the TRUE solution legitimately
deviates from the affine field — a plain-displacement std brick "fails" this exactly like the ANS
element does (~40% at the interior node; replica-verified). This is an ILL-POSED TEST, not element
failure. **Fix:** choose the patch gradient with `eps_13 = eps_23 = 0` and
`eps_33 = −lam(eps_11+eps_22)/(lam+2mu)` (so `sigma·e_z = 0`); the interior node then lands on the
affine field to machine precision (1e-16 in the numpy replica; 1e-6 through the Penalty solve) for
ans and std alike, and the `E33` channel is still exercised (`eps_33 != 0`). Full-traction GP-level
exactness belongs to the FULLY-PRESCRIBED single-element patches. Corollary for reviewers: a
solid-shell "patch test failure" report must state the face-traction handling before it counts.
