---
title: ADR-41 C1 design / handoff — the mortar kernel (clip → Gauss → D, M, g̃)
project: Ladruno
status: handoff (design, not yet implemented)
owner: nmora
tags:
  - implementation
  - contact
  - mortar
  - handoff
---

# ADR-41 Track C1 — the mortar kernel: `LadrunoMortarKernel.h`

> **START HERE (next session).** Tracks A1+A2 shipped (PRs #365/#367): the shared
> projection kernel (`LadrunoContactProjection.h`, with `projectFull`) and friction kernel
> (`LadrunoFrictionKernel.h`, unified cone) are merged on `ladruno`. C1 is the **headline
> accuracy feature** and the first *net-new mortar mechanics* (not an extraction). Build the
> **mortar integration kernel** + a numpy oracle; the hard part is the **overlap clip**. ALM
> enforcement is C2, friction is C3 — C1 is geometry + the patch test only.
>
> Read first: this doc, then ADR-41 §How "Mortar integration scheme (clipped Gauss-point
> mortar)", the capstone [[48_ladruno_contact_capstone_adr]] C1 row + contracts, and the
> shipped `SRC/domain/contact/LadrunoContactProjection.h` (the API C1 consumes).

## Deliverable (C1 scope only)

| File | Role | C1? |
|---|---|---|
| `SRC/domain/contact/LadrunoMortarKernel.h` | **header-only, OpenSees-free** mortar segment integration: slave/master **overlap clip** → sub-triangle Gauss → `D_IJ`, `M_IK`, weighted gap `g̃_I` (+ stubs for the linearization `dD,dM,dn,dξ`). | **✅ build now** |
| `Ladruno_implementation/contact_prototypes/proto_c1_mortar.py` | numpy oracle: clip-area, partition-of-unity `ΣM=ΣD`, **constant-pressure patch on a non-matched mesh**. | **✅ build now** |
| `LadrunoMortarPair.{h,cpp}`, `LadrunoMortarSegment.{h,cpp}` | per-GP state + narrow-phase driver / adapter wiring | ❌ **C2** (don't build in C1) |

C1 is **purely additive** (new header + oracle; no NTS-path touch) → null-mortar stays byte-identical;
no classTag, no manifest. Header-only, stamp via `Ladruno_scripts/stamp_headers.py` GLOBS (add the
new path, run it — see A1/A2 precedent).

## The mechanics (ADR-41 §How, made concrete)

For one **slave facet** (`nps_s` nodes, coords `Xs`) paired (via the ADR-39 broad phase) with a
**master facet** (`nps_m` nodes, `Xm`):

1. **Clip.** Project both facets onto a common **auxiliary plane**; Sutherland–Hodgman clip the
   slave∩master overlap polygon; **triangulate** it (fan from centroid). *(Details below — this is
   the crux.)*
2. **Gauss rule** per sub-triangle (default 3-pt; `-ngp` sets the slave-facet rule order). For each
   Gauss point with weight `w` and sub-triangle Jacobian `J` (physical-area scaling):
   - recover the slave parametric coord `ξ_s` of the GP (inverse-map within the slave facet) → slave
     shape fns `N_I^s(ξ_s)` and the physical slave point `x_s = Σ N_I^s X_I^s`;
   - project `x_s` onto the master facet via the **shipped** kernel:
     `LadrunoContactProjection::projectFull(nps_m, Xm, x_s, refDir, p)` → `{p.xi,p.eta (ξ_m), p.gap,
     p.n, p.t1,p.t2, p.g[2][2], p.phi (φ_K^m)}`.
3. **Accumulate** (standard, NOT dual, basis):
   - `D_IJ += w·J · N_I^s · N_J^s`            (slave–slave mortar matrix)
   - `M_IK += w·J · N_I^s · φ_K^m`            (slave–master mortar matrix)
   - **weighted gap, normal INSIDE the integral** (ADR-41 fixed D1's flat-facet factoring):
     `g̃_I += w·J · N_I^s · ( (x_s − x_m(ξ̄))·p.n )`   where `x_m(ξ̄) = Σ φ_K^m X_K^m`.
4. **B-operator** `B_gp` spans slave-facet nodes ∪ reachable master nodes (the ADR-39
   conservative-static-superset connectivity over a topology epoch).
5. Contact force / tangent (C2 consumes these): `F_c += w·J·B_gpᵀ·t`, `K_c += w·J·B_gpᵀ·D·B_gp`,
   with `t`,`D` from the normal law + (C3) `LadrunoFrictionKernel`. **In C1, just return `D`, `M`,
   `g̃` and the per-GP geometry** — no traction yet.

## The overlap clip (the crux — get this right first, in the oracle)

This is the only genuinely new algorithm. Standard segment-based mortar (Puso–Laursen 2004; Popp et
al. 2010; Wriggers, *Computational Contact Mechanics*):

1. **Auxiliary plane.** Build it from the master facet: a reference point `x0` (master centroid) and
   a unit normal `n0` (averaged master facet normal, or `projectFull`'s `p.n` at the slave centroid).
   Two in-plane axes `e1,e2` (orthonormal, e.g. Gram–Schmidt from `t1`).
2. **Project to 2D.** Map every slave node and master node onto `(e1,e2)` → 2D polygons `Ps`, `Pm`.
   Ensure consistent winding (CCW); reverse if the signed area is negative.
3. **Clip.** Sutherland–Hodgman `Ps ∩ Pm` (both convex for tri/quad ⇒ the intersection is a convex
   polygon, ≤ `nps_s+nps_m` vertices). Guard degeneracies: empty overlap (no contribution, skip),
   slivers (area < `tol·min(area_s,area_m)` ⇒ skip), near-duplicate vertices (merge within tol).
4. **Triangulate** the clip polygon: fan from its centroid → sub-triangles, each with a 3-pt (or
   higher) Gauss rule; `J` = 2·(sub-triangle 2D area)·(physical/parametric area ratio). Keep the
   area accounting honest — the patch test will catch a wrong `J`.
5. **Back-map each GP** to `ξ_s` (inverse bilinear/affine map within the slave facet — closed form
   for tri, Newton for quad) and to `ξ_m` (via `projectFull`). `projectFull` gives `ξ_m`, the gap,
   the metric, and `φ_K^m` directly — reuse it, do NOT re-derive projection.

**Do the clip + all of step 3–4 in `proto_c1_mortar.py` first** and pass the patch test in numpy
BEFORE transcribing to the header. That is the fork discipline and it will save you here.

## The interface C1 consumes (shipped, A2 — do not reimplement)

```cpp
// SRC/domain/contact/LadrunoContactProjection.h
struct LadrunoContactProjection::Projection {
    int status; double xi, eta, gap; double n[3], t1[3], t2[3], g[2][2], phi[4];
};
int LadrunoContactProjection::projectFull(int nps, const double X[4][3], const double xs[3],
                                          const double refDir[3], Projection &p,
                                          double tolR = 1e-12, int maxIt = 10);
// status: 0 conv in-bounds, 1 conv out-of-bounds, -1 invalid (degenerate / ambiguous normal).
```

## The gate (capstone C1 row + ADR-41 P1) — falsifiable, numeric

1. **Partition of unity:** row-sums `Σ_K M_IK == Σ_J D_IJ` to **1e-12** (the clip + back-map are
   consistent ⇔ this holds). This is the cheap always-on check.
2. **Constant-pressure patch test on a deliberately NON-MATCHED 2-block mesh** → recovered interface
   pressure constant to **≤ 1e-6 relative**. *This is THE headline falsifier* — a hard numeric
   pass/fail, not an oscillation-magnitude metric. Un-clipped "mortar-lite" provably fails it (that
   is why clipping is in the MVP, not deferred).
3. **`K_c` vs finite-difference on a ROTATED config to 1e-6** (once C2 assembles `K_c`; in C1 the
   oracle can FD-check `D,M,g̃` w.r.t. nodal coords).

## Scope discipline (don't pull C2/C3 work into C1)

- **Standard basis only** (D non-diagonal, per-facet `D`-solve). Dual/biorthogonal basis is **ADR-47**.
- **Frictionless geometry** — no traction, no `λ`, no ALM in C1. Enforcement (commit-cycle ALM) = C2;
  friction (adopt `LadrunoFrictionKernel` in the `λ_T` form) = C3.
- **Faceted per-GP normal** is fine for C1 (nodal-normal smoothing is Q-NORMAL → ADR-47).
- Keep `LadrunoMortarKernel.h` header-only + OpenSees-free (raw `double[]` + `<cmath>`), numpy-oracle
  testable — same discipline as the A1/A2 kernels.

## Risks / hard parts (flagged honestly)

- **Clip robustness** is the main risk: degenerate/empty overlaps, sliver polygons, coincident
  vertices, and the auxiliary-plane choice when slave and master are strongly skewed. Bound every
  loop, skip-with-tol on degeneracies, and make the oracle exercise a skewed non-matched mesh.
- **GP back-map for quads** needs a small bounded Newton (inverse bilinear map) — reuse the
  bounded-Newton pattern from `LadrunoContactProjection::project` (cap iters, degenerate guard).
- **Jacobian bookkeeping** (2D clip area → physical integral measure) is the silent-error spot; the
  partition-of-unity + patch test are the guards — trust them, not the algebra by eye.
- **D is not diagonal** (standard basis) → the per-node `λ` recovery is a per-facet local solve, not
  a nodal division. Fine for C1 (we only assemble D,M,g̃); relevant once C2 does the ALM update.

## References

- ADR-41 (`41_ladruno_mortar_alm_contact_adr.md`) §How "Mortar integration scheme"; capstone
  [[48_ladruno_contact_capstone_adr]] C1 row + the 8 contracts.
- Shipped: `SRC/domain/contact/LadrunoContactProjection.h` (`projectFull`),
  `LadrunoFrictionKernel.h` (C3 will consume), oracle style in `contact_prototypes/proto_a2_*.py` /
  `proto_a1_*.py`.
- Literature: Puso & Laursen (2004) mortar segmentation; Popp, Gee, Wall (2010) dual/std mortar;
  Wriggers, *Computational Contact Mechanics* (clip + segment integration). Skills:
  `fem-mechanics-expert`, `opensees-expert`.
- Workflow: `gh` is installed — `gh pr create/merge --admin`; oracles run with **`python3.12`**
  (numpy); watch **Zone-A** before merging C++ (auto-merge does not gate the build). See memories
  [[ladruno-pr-via-rest-api]], [[ladruno-oracle-python]] and `Ladruno_internal/WORKFLOW_GOTCHAS.md`.

## Oracle-first checklist (do these in `proto_c1_mortar.py`, in order)

1. Sutherland–Hodgman clip of two convex polygons → correct intersection area (vs a known case).
2. Tri+quad slave/master, matched mesh → `D`, `M`; `ΣM_IK == ΣD_IJ` to 1e-12.
3. **Non-matched mesh** (offset/rotated master) → still `ΣM=ΣD`; weighted gap `g̃` linear-exact.
4. **Constant-pressure patch test**: impose a uniform interface pressure on a 2-block non-matched
   mesh, recover it through `D⁻¹M` transfer → constant to ≤1e-6. (The headline gate.)
5. FD-check `∂{D,M,g̃}/∂(nodal coords)` (feeds the C2 `K_c` linearization).

Only after 1–4 pass in numpy, transcribe to `LadrunoMortarKernel.h`, then re-confirm via a C++
standalone `g++ -I SRC/domain/contact` compile + (in C2, when the adapter exists) Zone-A.
