---
wp: LEGACY
title: "Mass scaling multiplies the support-motion START SHOCK by sqrt(s) — a Linear sp ramp that is harmless unscaled can crack elements at the moving support under S…"
legacy_seq: 152
---
### Mass scaling multiplies the support-motion START SHOCK by sqrt(s) — a Linear sp ramp that is harmless unscaled can crack elements at the moving support under SMS
- **Bites:** quasi-static explicit runs driven by prescribed support motion (`sp` + timeSeries, the #333 recipe) under `CentralDifferenceSMS`. A `Linear` series applies a velocity STEP v at t=0; the wave it launches carries sigma ~ rho'*c'*v = sqrt(rho'*E)*v — and mass scaling inflates rho' by s = (dtTarget/dt_e)^2, so the shock stress grows by sqrt(s) = the dtTarget factor. Measured (ADR 66 G9c): a rate whose unscaled shock is a trivial 0.35 MPa hit 3.5 MPa > ft at the 10x target and cracked the pulled-face concrete element outright (omega_t -> 0.97) before any real loading happened.
- **Why:** impedance rho*c = sqrt(rho*E); uniform scaling multiplies rho by s and leaves E alone. The shock rides the SCALED impedance while the "quasi-static" rate was budgeted against the unscaled one.
- **Workaround/status (2026-07-06, ADR 66 G9):** drive support motion with a SMOOTHSTEP displacement protocol (lam_end * u^2(3-2u), zero start/end velocity — a ~60-point Path series), never a raw Linear ramp, whenever mass scaling is on; budget the ramp against the SCALED wave speed c' = c/sqrt(s). Same discipline applies to any velocity IC under SMS.
