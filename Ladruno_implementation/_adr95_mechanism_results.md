# ADR-95 — collapse mechanism at the plateau vs the Prandtl–Reissner geometry (2026-09-08)

Repaired UW DruckerPrager, build c945f9a8b, h0 = 1.0 hex mesh / the tet mesh, no diagnostics.
Velocity field = Δu between the first converged steps past s/B 0.14 and 0.15; ε_q = incremental
deviatoric strain magnitude per element; analytical overlay for φ_ps = 27.470° (the harness's own
plane-strain angle from `mc_from_cone`): wedge depth 1.647 m, fan r0 = 1.927 m, passive outcrop
7.454 m = 3.73 B from each footing edge. Script `mechanism_snapshot.py`, dumps `mech_B?.npz`,
figures `adr95_mechanism.png` / `adr95_mechanism_full.png`.

| leg | element | q/q_exact at 0.15 | tail dq/d(s/B) kPa | vol/dev increment ratio | steps / failed | wall |
|---|---|---|---|---|---|---|
| B1 | LadrunoBrick -bbar | 1.0850 | 0.4 | 0.000 | 310 / 0 | 29 s |
| B2 | LadrunoBrick20 uri | 0.9618 | 3.0 | 0.095 | 1515 / 1034 | 522 s |
| B5 | BezierTet10 -bbar | 1.0403 | 0.9 | 0.001 | 1515 / 0 | 656 s |
| B4 | BezierTet10 std | 1.1824 | 5.2 | 0.029 | 1515 / 0 | 618 s |
| B3 | TenNodeTetrahedron | 1.1704 | 5.2 | 0.029 | 1515 / 0 | 1044 s |

Reading:
- Every leg is on a plateau (tail ≤ 5 kPa per unit s/B against an initial ~8000) and every velocity
  field is a shear mechanism whose extent matches the Prandtl passive outcrop (~7.5 m each side).
- **Nothing localises.** The ε_q ridge is 2–3 elements thick everywhere; no leg draws the log-spiral
  fan or the wedge faces — with 1 m elements against a 1.65 m wedge there are one to two elements
  across the whole fan. Confirms note 82 §6 with a velocity field instead of a total-strain one.
- **Isochoric check:** the two b-bar elements are exactly isochoric (0.000 / 0.001). The two
  standard-integration tets carry 3 % spurious volumetric increment and the H20 uri 9.5 % — the
  reduced-integration / plastic-rank signature (H4) showing up as a non-isochoric collapse field,
  and the same leg is the one with 1034 failed attempts and the lowest plateau.
- Coarse-mesh over-strength orders as expected: locked tets 1.17–1.18 > linear b-bar 1.085 >
  Bezier b-bar 1.040 > H20 uri 0.962 (under, the plastic-hourglass side).
Not claimed: convergence of any of these numbers; the fan cannot be resolved on this mesh.
