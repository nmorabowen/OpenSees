---
wp: LEGACY
title: "Backward_Euler needs strict_convergence to fail loud on a starved deck; Closest_Point fails loud on the SAME starved deck by default -- do not reuse a Closest_…"
legacy_seq: 423
---
### `Backward_Euler` needs `strict_convergence` to fail loud on a starved deck; `Closest_Point` fails loud on the SAME starved deck by default -- do not reuse a `Closest_Point` starvation reproducer for `Backward_Euler` without adding the flag
- **Bites:** `test_adr97_p6_failloud.py::_starved` (`hiso=7000.0, niter=1` on the VM triaxial path) is a `Closest_Point` reproducer (`mat_vm`'s default `method`) -- it refuses without `strict_convergence` because `Closest_Point`'s Newton has no legacy silent-accept path. Reusing the exact same kwargs with `method="Backward_Euler"` (naively assuming the reproducer is method-agnostic) SILENTLY SUCCEEDS: `Backward_Euler` only fails loud on non-convergence when `strict_convergence` is explicitly turned on (ADR-84 P2a) -- by default it falls out of its Newton loop and commits the non-converged state as "success", which is the ORIGINAL upstream defect ADR-84 P2a's flag exists to opt out of.
- **Fix:** always pass `strict=1` (or the raw `strict_convergence 1` option) when writing a NEW `Backward_Euler` starvation test, even if copying an existing `Closest_Point` reproducer's `niter`/load values verbatim; `strict_convergence 1` is harmless (bit-identical) on a deck that already converges, so it is safe to add unconditionally to a starvation reproducer regardless of which integrator it targets.
