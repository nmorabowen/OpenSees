---
wp: LEGACY
title: "A frictional verification model can be 95 % mobilised by its own gravity state, because the elastic K0 is what sets the initial stress ratio"
legacy_seq: 256
---
### A frictional verification model can be 95 % mobilised by its own gravity state, because the elastic K0 is what sets the initial stress ratio
- **Bites:** any limit-analysis / bearing-capacity check that picks a "nice
  easy friction angle" to validate the machinery against a closed-form answer.
  Under 1-D gravity or surcharge loading the initial stress ratio is the
  ELASTIC `K0 = nu/(1-nu)`, not `1 - sin(phi)`, so the mobilisation of the
  yield surface before anything is loaded is
  `m = (1-K0) / (sqrt(3) * alpha * (1+2*K0))` and depends only on nu and the
  surface — never on the load magnitude. At PDMY's moduli (`K = 1.5e5`,
  `G = 5.5e4` => nu = 0.3366, K0 = 0.507) that mobilises 19.1 deg. Validating
  against a phi_txc = 20 deg cone therefore starts with the WHOLE DOMAIN at
  **m = 0.950 of yield**, and the "verification" measures its own initial
  condition: 44.5 % of elements at m > 0.99 after 2 mm of footing settlement,
  the yielded zone already touching the roller boundary, and no convergence
  past s/B = 0.0019.
- **Tell:** a validation leg that dies almost immediately with a large fraction
  of the mesh yielded and the plastic zone at the boundary, while the SAME
  machinery runs fine on the (steeper, supposedly harder) real surface.
- **Rule:** compute `m` of the initial state and print it before the run; treat
  `m > 0.8` as void. The fix is free — a collapse load of an
  elastic-perfectly-plastic body does not depend on its elastic constants, so
  raise nu for the verification legs only (nu = 0.45 gives K0 = 0.818,
  m = 0.268) and MEASURE the independence with a nu-pair rather than asserting
  it. *2026-07-30 (ADR-79 collapse study).*
