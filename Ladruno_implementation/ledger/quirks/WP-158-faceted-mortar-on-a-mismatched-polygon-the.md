---
wp: WP-158
title: "Faceted mortar on a mismatched polygon: the geometric part of the pair force is non-smooth, so an \"exact\" FD tangent can be WORSE than the frozen-geometry one…"
legacy_seq: 542
---
### Faceted mortar on a mismatched polygon: the geometric part of the pair force is non-smooth, so an "exact" FD tangent can be WORSE than the frozen-geometry one (WP-158)
- **Bites:** at the meshed configuration a skin polygon (n) and a hole polygon (n+4) share vertices. The clipped
  overlap of a pair changes topology inside +-h, and the one-sided slopes of the pair force differ by about 100 %
  (|df/du| ~ 1e6 against epsN*a ~ 3e5 under a 1e4 kPa prestress). The central-FD tangent (the WP-158 FD oracle patch) then stalls
  Newton at 0.1-5 for every step size, while the analytic tangent (no geometric terms) plus `-consistanttan` converges.
  On a smooth crease the same FD tangent is quadratic.
- **Workaround/status:** for the pile use `-consistanttan` (Pardiso is mtype 11, non-symmetric). The FD pair tangent is
  not shipped; it is parked as `contact_prototypes/adr158_fd_pair_tangent_oracle.patch` (a diagnostic oracle, useful on
  smooth creases with cross-crease pairs). See [[158_mortar_tangent_diagnosis_consistanttan]] §3 D1, D3.
