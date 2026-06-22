# ADR-39 P2b — implementation plan: faceted-master NTS penalty contact (frictionless)

> Concrete build plan for P2b, grounded in the gated `_adr39_p2_design.md` (the 4-lens
> gate wp3cr60mf design decisions are LAW here) + the precedent map (SimpleContact3D,
> ZeroLengthContactASDimplex, ContactMaterial3D, LadrunoBrick). P2a shipped (#346 +
> gate-fix #350). Parent ADR `39_..._adr.md`; driver `_adr39_loop_state.md`.

## Sub-rung split (manage risk — mirrors the P2a/P2b split one level down)

P2b itself is large (projection kernel + faceted B + ∂n/∂u tangent + deformable master
+ auto-kn + Hertz). Split into shippable rungs, each design→code→build→test→gate:

- **P2b-1 — projection kernel + single slave vs single FIXED master segment.**
  Master nodes FIXED (no master compliance yet, but REAL bilinear/linear projection,
  derived outward normal, faceted B = [nᵀ | −φ_i nᵀ], penalty force, and the
  kₙBᵀB + ∂n/∂u tangent). Connectivity = {slave} ∪ {master-seg nodes}. `-kn $val`.
  ORACLE = hand-placed `ZeroLengthContactASDimplex` node-pair (rel 1e-6) when the
  slave projects onto a master node with the segment's flat normal. Isolates the
  projection + faceted assembly + tangent from deformable-master + auto-kn.
- **P2b-2 — deformable master + `-kn auto` + SOFT floor + Hertz + winding-flip gate.**
  Two deformable `LadrunoBrick` blocks; cache K/A/V at setDomain; FD-on-rotated tangent
  gate (flat PASS, curved O(gₙ·κ) shrinks); Hertz refinement benchmark.

Friction = P3. Bucket sort (multi-segment broad phase) = P2.5 — until then P2b-1 scans
all segments of the one master surface per slave (fine for the gate meshes).

## The kernel — `SRC/domain/contact/LadrunoContactKernel.h` (NEW, header-only, OpenSees-free)

Mirror `LadrunoJ2Kernel.h`: pure functions on raw `double[]`, no OpenSees types, so
both OPS_Analysis (the FE) and OPS_Domain (the surface) include it with no link
inversion. All math lives here; the FE is a thin adapter. Functions:

```
namespace LadrunoContactKernel {
  // bilinear (4-node) / linear (3-node) shape fns + tangents at (xi,eta)
  void shape(int nps, double xi, double eta, double N[4]);          // nps=3|4
  void interp(int nps, const double N[4], const double X[4][3], double xbar[3]);
  void tangents(int nps, double xi, double eta, const double X[4][3],
                double g1[3], double g2[3]);                        // ∂x̄/∂ξ, ∂x̄/∂η
  void normalRaw(const double g1[3], const double g2[3], double n[3], double &jac);

  // BOUNDED closest-point projection (GATE BLOCKER-2): cap 10 iters,
  //   reject if |detK| < eps*|g1||g2|, return convergence flag.
  //   tri-3: closed-form (de-risk); quad-4: bounded Newton on (x_s−x̄)·g_α=0.
  // returns: 0 = converged in-bounds, 1 = converged out-of-bounds, -1 = no valid proj
  int project(int nps, const double X[4][3], const double xs[3],
              double &xi, double &eta, double tolGap, int maxIt);

  // gap gₙ = n·(x_s − x̄); penetration ⇔ gₙ<0 AND in-bounds
  double gap(const double n[3], const double xs[3], const double xbar[3]);

  // outward-normal orientation (GATE BLOCKER-1): flip n so n·(x̄−c_elem)>0
  void orient(double n[3], const double xbar[3], const double cElem[3]);

  // penalty traction (Macaulay, GATE MINOR-15): tn = kn*<−gₙ>₊ ; 0 otherwise
  double traction(double kn, double g);

  // B-operator (1 x ndof), ndof = 3*(1+nps): [ nᵀ | −φ_1 nᵀ ... −φ_nps nᵀ ]
  void Bop(int nps, const double n[3], const double N[4], double B[/*3*(1+nps)*/]);

  // residual r = Bᵀ tn  (self-equilibrated: Σf=0, ΣM=0 at converged projection)
  // tangent  K_c = kn BᵀB + (∂n/∂u block), SYMMETRIC for frictionless,
  //          dropping only O(gₙ·κ) curvature.  See "Tangent" below.
  void residual(int nps, const double B[], double tn, double r[]);
  void tangent (int nps, const double X[4][3], const double n[3], const double N[4],
                double xi, double eta, double kn, double g, double K[][]);
}
```

### Tangent derivation (the crux — GATE BLOCKER-3 / Q-P2-tan)

Penalty contact virtual work (active, gₙ<0): δW = tₙ δgₙ with tₙ = kₙ⟨−gₙ⟩₊,
so δW = −kₙ gₙ δgₙ. Consistent linearization:
  Δ(δW) = −kₙ (Δgₙ)(δgₙ) − kₙ gₙ Δ(δgₙ)
        = kₙ (δgₙ)(Δgₙ)             [MAIN, since ∂tₙ/∂gₙ=−kₙ and δW=−kₙgₙδgₙ ⇒
                                      first term = kₙ δgₙ⊗Δgₙ = kₙ BᵀB]   ✓ [T1]
        + tₙ Δ(δgₙ)                  [GEOMETRIC: normal-variation + curvature]

δgₙ = nᵀ(δx_s − δx̄). At the converged projection (x_s−x̄)=gₙ n. The geometric term
tₙΔ(δgₙ) splits (Wriggers, *Computational Contact Mechanics* §9, NTS):
  - **normal-variation / metric part** (∂n/∂u): NONZERO even for a FLAT facet —
    this is what makes FD-on-rotated nonzero. Built from the tangent vectors g_α,
    the inverse metric m^{αβ}=（g_α·g_β)⁻¹, and the shape-fn tangential derivatives
    N_{,α}. KEEP it. Symmetric for frictionless.
  - **curvature part** O(gₙ·κ_αβ): for tri-3/quad-4 flat facets κ=0; for curved
    (warped quad) it is O(gₙ) → DROP (gate allows; FD-on-curved carries the residual).

Concrete (Wriggers NTS, single facet, frictionless), with T_α = [0 | −N_{,α,i} I]
(the ∂x̄/∂ξ_α displacement operator rows) and D_α = [nᵀ-deriv ...]:
  K_c = kₙ BᵀB  −  tₙ ( m^{αβ} (N_{,α}ᵀ n)(N_{,β}ᵀ n)ᵀ-type coupling )  + ...
IMPLEMENTATION STRATEGY (gate-sanctioned): code the MAIN term kₙBᵀB exactly +
the first-order normal-variation block from g_α, m^{αβ}, N_{,α}; then **the
FD-on-flat-ROTATED test is the oracle** — if it fails, refine the block. Do NOT
ship main-term-only (gate proved it fails FD-on-rotated by construction). For
P2b-1 the master segment is FIXED, so under a SLAVE-only perturbation the ∂n/∂u
(master-node) rows are inactive and the main term suffices for the slave block;
the full ∂n/∂u block is exercised/gated in P2b-2 (deformable master, rotated block).

> NOTE: getTangent returns BARE K_c (no c1) — Newmark scales via addKtToTang(c1),
> exactly as P2a. addKiToTang mirrors it. Routed through formEleTangent so explicit
> stays mass-only (the P2a lesson — keep it).

## Surface side — `LadrunoContactSurface` + `LadrunoContactDomain` (OPS_Domain)

- `LadrunoContactSurface` already stores MASTER_SEGMENTS flat connectivity + stride
  (3=tri,4=quad). ADD at setDomain (called from the FE/handler side, which has the
  Domain): resolve `Node*` per tag, cache reference coords, and the **derived outward
  normal per segment** (needs the master solid-element centroid — see below).
- `LadrunoContactDomain::addContact(tag, masterSurf, slaveSurf, kn, kt, mu)` already
  exists (P1b). P2b uses it (kt,mu ignored until P3). The rigid-plane path (P2a)
  stays; this is the segment path.
- **Outward normal orientation (GATE BLOCKER-1):** flip each segment normal so
  n·(x̄−c) > 0 where c = master SOLID element centroid. The master surface stores
  face nodes only → need the parent element. P2b-1 (fixed master): accept an explicit
  `-outward px py pz` reference point (or the centroid of ALL master-surface nodes as
  a proxy) to orient; P2b-2 wires the real master-element centroid. Winding-flip gate.

## FE side — `LadrunoContactFE` segment mode (3rd ctor)

Add `Mode SEGMENT`. Ctor: `(tag, Node* slave, ndm, masterNodes[], nps, kn, ContactDomain*)`.
- connectivity: myDOF_Groups = {slave DOF_Group, master node DOF_Groups}, ndof=3*(1+nps).
- getResidual: project slave onto the (one) master segment (P2b-1: the single segment;
  multi-seg later) via kernel; if no valid in-bounds penetrating projection → zero.
  Else tₙ=kₙ⟨−gₙ⟩₊, r = Bᵀtₙ.
- getTangent: route through integrator->formEleTangent (P2a pattern); addKtToTang /
  addKiToTang fill kernel `tangent(...)`.
- **Projection cadence (GATE MAJOR-9):** implicit RE-PROJECTS per Newton iter (in
  getResidual/getTangent, reading trial disp); explicit freezes per step (CDL only
  forms residual once). Keep (ξ̄,x̄,n,B) a consistent set per evaluation.
- **getID superset + immutability (GATE Q-P2-conn):** surface immutable after
  setDomain; assert projected segment nodes ⊆ getID in getResidual (opserr+skip).

## Parser

Reuse P1b `contactSurface tag -master n1 n2 n3 n4 ... -nps {3|4}` (verify the P1b
parser captures stride) + `contact tag masterSurf slaveSurf -kn {val|auto} [-outward px py pz]`.
P2b-1: `-kn $val` + `-outward` (or proxy). P2b-2: `-kn auto`.

## Tests — `tests/test_adr39_contact_p2b.py`

P2b-1 (this rung):
1. **ASDimplex oracle** (rel 1e-6): one slave node above a single FIXED flat master
   quad (4 fixed nodes), `-kn`. Place a hand `ZeroLengthContactASDimplex(masterNode,
   slaveNode, Kn, 0, 0, ndm=3, 0, nx,ny,nz)` with the SAME kn + normal where the slave
   projects onto a master node → static penetration force matches rel 1e-6.
2. **interior projection penetration** = P/kₙ (slave pressed onto segment interior,
   flat normal) static + explicit.
3. **self-equilibrium** |Σf| < 1e-10·|F_c| (slave + 4 master reactions sum to 0).
4. **oblique 30°/45° master** off-axis normal — static penetration along n.
5. **winding-flip → identical force** (GATE BLOCKER-1 gate): reverse the master node
   order → same contact force (normal re-derived, not trusted).
6. **out-of-bounds / no-projection** → zero force (slave beside the segment).
P2b-2: deformable blocks, FD-on-rotated tangent (flat PASS / curved shrinks), Hertz,
auto-kn, mesh-ratio sweep, convergence-rate.

## Files (P2b-1)
- NEW `SRC/domain/contact/LadrunoContactKernel.h` (header-only; stamp_headers + ledger).
- `LadrunoContactSurface.{h,cpp}` — setDomain coord/normal cache + segment iteration.
- `LadrunoContactDomain.{h,cpp}` — segment-contact accessor for the handler.
- `LadrunoContactFE.{h,cpp}` — SEGMENT ctor + project/residual/tangent via kernel.
- `LadrunoContactHandler.cpp` — build one segment adapter per slave node.
- parser (OpenSeesOutputCommands.cpp) — `-kn`/`-outward` on `contact`; verify `-nps`.
- `tests/test_adr39_contact_p2b.py`.

## Open items to resolve while coding P2b-1
- Q1: does the P1b `contactSurface -master` parser already store nodesPerSeg? (read it)
- Q2: tri-3 closed-form projection vs quad-4 bounded Newton — implement both; gate quad.
- Q3: orient proxy for fixed master (explicit `-outward` vs all-master-node centroid).
- Q4: exact ∂n/∂u block form — derive, then FD-gate in P2b-2 (master fixed in P2b-1).
