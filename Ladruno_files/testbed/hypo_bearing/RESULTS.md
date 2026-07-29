# ADR-79 bearing-backbone campaign — results

**Run 2026-07-28/29 on the merged `ladruno` state through PR #682
(`dist/bin` build `10db451c9`), graded 2816-element benchmark mesh, three legs
in parallel, 3.2–4.1 h each.**

## Verdict

**Climbing the geometry ladder does not produce a bearing limit point, and it
does not shrink the over-hardening — it makes it slightly worse.**

1. **No limit point on any rung.** Every leg is still hardening where it ends:
   the tangent `dq/ds` decays to 15 % (corot) / 25 % (hypo) of its initial
   value but never approaches zero, and `q` is monotone with its maximum at the
   last converged point. Genuine large strain does not bend this backbone over.
   `hypo` reaches **3.62 × Vesic at s/B = 7.6 %** and is still climbing.
2. **`hypo` is *stiffer* than `corot`, and the gap grows with penetration** —
   +2.1 % at s/B = 1 %, +6.2 % at 2.5 %, +14.8 % at 5 %. The sign matches
   ADR-79 P3's smoke (hypo stiffer; UL assembly on the *compacted*
   configuration plus the kinematic n(J) porosity), and it is the opposite of
   the "missing large deformation is why there is no limit point" hypothesis in
   ADR-78 §1. Large strain moves this model **away** from a limit point, not
   toward one.
3. **`-kozenyCarman` is inert here.** hypo and hypo+kc agree to ~1e-4 relative
   over the whole backbone. That is the physically expected answer, not a null
   result to explain away: the push is fully drained by construction, so
   permeability only sets how fast pore pressure dissipates, and at these rates
   every admissible k gives the same drained backbone. k(J) is a lever for
   *undrained or rapid* loading, and this campaign cannot see it.

So the ADR-78 §1 over-hardening is **not** a large-deformation-kinematics
artifact. Scoping finding 3 points at the much bigger lever: boundary
confinement. The same footing in a box with 0.5–1.5 B of clearance ran to
27 MPa (tens of × Vesic) purely as an oedometer, whereas this 4.5 B-clearance
benchmark is at 2.7–3.6 × Vesic at the same penetrations. We cannot attribute
a specific share of the reported 9.3 × to the SFIM model's boundaries without
that model — but the ladder is now excluded as the explanation.

*(A note on what "no limit point" does and does not mean here. A drained,
smooth, displacement-controlled footing on a hardening PDMY sand with
confinement-dependent moduli need not exhibit a load maximum at all: as the
footing sinks, the mechanism mobilizes deeper, better-confined soil whose
stiffness grows with p′. The campaign's finding is comparative and specific —
climbing the geometry ladder does not introduce one, and does not reduce the
over-hardening — not a claim that this model class must have one.)*

![backbones](bearing_backbone.png)

## Numbers

`q_Vesic = 637.5 kPa` = `q0·Nq·sq` (430.4) + `0.5·γ'·B·Nγ·sγ` (207.1),
φ = 33°, γ' = 9.81 kN/m³, B = 2 m, q0 = 10 kPa.

q/q_Vesic at settlement checkpoints (`--` = the leg never reached it):

| s/B | `-geom corot` | `-geom hypo` | `-geom hypo -kozenyCarman` | hypo/corot |
|---|---|---|---|---|
| 0.010 | 0.98 | 1.00 | 1.00 | 1.021 |
| 0.025 | 1.76 | 1.87 | 1.87 | 1.062 |
| 0.050 | 2.48 | 2.85 | 2.85 | 1.148 |
| 0.075 | -- | 3.59 | 3.59 | -- |
| 0.100 | -- | -- | -- | -- |
| 0.150 | -- | -- | -- | -- |

Per leg, over the span each actually covered:

| leg | reached s/B | steps | retries | q_end | q_end/q_Vesic | final dq/ds | ended by |
|---|---|---|---|---|---|---|---|
| `-geom corot` | 0.0602 | 677 | 111 | 1695 kPa | 2.66 | 14.8 % of initial | convergence floor |
| `-geom hypo` | 0.0762 | 927 | 166 | 2306 kPa | 3.62 | 25.2 % of initial | operator stop at 4 h |
| `-geom hypo -kozenyCarman` | 0.0760 | 923 | 165 | 2303 kPa | 3.61 | 25.2 % of initial | operator stop at 4 h |

`q` crosses the Vesic capacity at s/B ≈ 1 % on every rung and keeps climbing —
the over-hardening reproduces cleanly on a properly sized benchmark.

## Why the legs stop where they stop

**No leg reached the s/B = 0.15 target, and the two reasons are different — the
distinction matters, so it is recorded rather than smoothed over.**

- `corot` ended on its own at s/B = 0.0602: the adaptive increment fell below
  the 6.25 µm floor. That is a convergence failure on a heavily yielded,
  badly distorted near-footing mesh, *not* a limit load.
- `hypo` and `hypo -kozenyCarman` were **stopped by the operator** at
  s/B ≈ 0.076 after 4 h of wall clock, while both were still converging
  comfortably (0.1–0.2 mm increments). They did not fail; they were cut short.
  Their backbones are truncated by the campaign's time budget, not by the
  physics or the solver.

Neither ending is a limit point, and the two are distinguishable from one: at a
genuine limit point `dq/ds` → 0, and here it is still 15 % (corot) / 25 %
(hypo) of its initial value and clearly positive.

The honest caveat is the converse one: neither a convergence-limited nor a
time-limited truncation can prove there is no limit point *past* the last
converged step. What the campaign does establish is that no rung produces one
**within the span all three legs share** (s/B ≤ 0.060), and that over that span
the ladder's effect is to harden.

Notably `corot` hit its convergence floor at s/B = 0.060 while `hypo` was still
stepping at 0.076 — the rate-form UL lane is better conditioned than
rotate-only kinematics once near-footing elements are genuinely distorted,
which is a small independent point in `-geom hypo`'s favour.

## Validation

Five independent checks, all run against this exact model:

| check | result |
|---|---|
| global vertical equilibrium after gravity | ΣR_z = 66 784.00 kN vs ρ_sat·g·V + q0·A = 66 784.00 kN — **exact to 0.0000 %** |
| effective stress at depth | σ'_zz ≈ −80.3 kPa at the base element ≈ γ'·z + q0 (buoyant), K0 = 0.51 vs 1−sinφ = 0.455 — buoyancy is carried by the pore-pressure field, so normalizing Vesic by γ' is the right choice |
| mesh geometry | volume 3200.000 m³ vs the exact 20 × 20 × 8 box |
| reaction bookkeeping | R extrapolates to 62.5 kN at s → 0 = q0 × 6.25 m² (the footing-node tributary) — the closed-form surcharge correction is exact |
| displacement control | commanded `s_m` ≡ measured `s_meas_m` to machine precision at every one of ~2200 recorded steps |
| solver | Pardiso vs UmfPack agree to **5.2e-13** relative on 16 push steps (Pardiso 2.2× faster) — the solver is a wall-clock choice only |

The gravity state is also uniform to 1e-16 across the whole top face
(1.7638 mm), i.e. a clean 1-D response before the footing is engaged, and pore
pressure is hydrostatic (0 at the surface, 78.48 kPa at z = −8 m).

## Caveats

- **Mesh coarseness.** B is resolved by 4 elements (0.5 m hexes). Adequate for
  the cross-method comparison — all three legs share the mesh exactly — but the
  *absolute* q is not mesh-converged, and no refinement study was run. A
  coarse mesh generally over-predicts capacity for a punching mechanism.
- **Truncation.** See above; no leg reached s/B = 0.15. corot stopped on its
  convergence floor at 0.060; the two hypo legs were stopped by the operator at
  0.076 while still converging, so their backbones could be extended simply by
  spending more wall clock.
- **Smooth footing.** Only u_z is driven on the 25 footing nodes; u_x, u_y are
  free. A rough footing would raise q and change the mechanism.
- **Boundary distance.** 4.5 B of side clearance and 4 B of depth is much
  better than the scoping boxes but is still a finite domain with roller sides;
  some confinement contribution remains.
- **Drained only.** k = 1e-4 m/s with a slow push. The undrained branch — where
  `-kozenyCarman` would actually bite — is untested.
- **PDMY calibration** is the generic medium-loose set from the P3 smoke, not a
  site calibration; the absolute multiples of Vesic should not be read as a
  statement about any real footing.
- **q includes the footing's own surcharge share** (≤ 1 %: 5.6 kPa of the
  600–2200 kPa signal), because the regularizing surcharge must cover the
  footing patch (scoping finding 2).

## Reproducing

```bash
C:/Users/nmb/venv/opensees_env/Scripts/python.exe Ladruno_files/testbed/hypo_bearing/build_mesh.py
```

```bash
C:/Users/nmb/AppData/Local/Programs/Python/Python312/python.exe Ladruno_files/testbed/hypo_bearing/bearing_backbone.py all
```

```bash
C:/Users/nmb/venv/opensees_env/Scripts/python.exe Ladruno_files/testbed/hypo_bearing/analyze.py
```

The mesh is cached in `bearing_mesh.npz`, so the second command reproduces
these CSVs exactly. Run the legs separately (`corot` / `hypo` / `hypo_kc`) to
parallelize; each takes ~3.5 h.
