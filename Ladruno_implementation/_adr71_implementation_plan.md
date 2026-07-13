---
title: "ADR-71 implementation plan — LadrunoUP execution (P0 detailed, P1 concrete, P2–P4 sketch)"
status: working doc — execution plan, updated as phases land
---

# ADR-71 implementation plan

Companion working doc to [[71_ladruno_up_family_adr]] (LOCKED 2026-07-10, PR #547).
This is the *execution* plan: work packages (WP), owners, inputs, done-whens.
Delete or fold into the ADR's Implementation log when the runway closes.

## Execution model

- **Main session ("MAIN", Fable)** — interface pinning, cross-verification,
  discrepancy adjudication, builds, PRs, CI. Never delegated.
- **Opus 4.8 agents ("OPUS")** — all substantive production work packages
  (oracle, kernel, element, parser, test batteries, guide drafts) and the
  adversarial panels. Launched in background, parallel where files are disjoint.
- **Independence protocol (load-bearing for P0):** the numpy-oracle agent and
  the C++-kernel agent work from THIS DOC's interface contract + the ADR only.
  The oracle agent must NOT read `LadrunoUPKernel.h`; the kernel agent must NOT
  read `ladruno_up_reference.py`. Two independent implementations agreeing to
  1e-9 *is* the verification — shared code would fake it.
- **Branch/PR discipline** (fork traps): one PR per phase, one branch per PR,
  branch off `ladruno` AFTER #547 merges (if starting earlier, branch off
  `guppi/focused-gould-1fa488` — clean fast-forward). Verify PR state==OPEN
  before every push (auto-merge strands follow-up commits).
- **Windows gotchas** (from memory, non-negotiable): write files with the Write
  tool (Bash heredocs eat backslashes); build via `cmd /c Ladruno_scripts\build.bat …`
  from PowerShell; worktree Python tests run `python -S` + manual paths +
  `assert opensees.__file__` (boot `.pth` preloads a stale pyd); close VS Code
  before installer/DLL-touching steps; new source files added to
  `stamp_headers.py` GLOBS and stamped.

---

## P0 — kernel + providers + oracle (no OpenSees build)

### Work packages

| WP | Owner | Deliverable | Done when |
|---|---|---|---|
| 0.A | MAIN | Interface contract (§below) frozen; P0 branch created | this doc committed on the P0 branch |
| 0.B | OPUS-1 | `tests/ladruno_up_reference.py` — independent numpy oracle: all blocks per shape via high-order/exact quadrature, P0 gates (consolidation ODE, inf-sup smoke, FD dyn-seepage tangent, symmetry/pattern), self-check mode with hard asserts, `--emit tests/ladruno_up_cases.txt` | self-check green; cases file emitted; NO C++ read |
| 0.C | OPUS-2 | `SRC/element/ladrunoUP/LadrunoUPKernel.h` + `LadrunoUPShapes.h` (pure, per contract) + `tests/ladruno_up_kernel_check.cpp` (standalone g++, ingests cases.txt, ≤1e-9) + CMakeLists stub for the dir (headers only, no build impact) | compiles standalone with g++; harness runs; NO oracle read; CI wiring mirrors whatever `finitestrain2d_kernel_check.cpp` does (agent inspects and mirrors) |
| 0.D | MAIN | Cross-run: oracle self-check + emit → C++ check vs cases ≤1e-9; adjudicate every discrepancy (spec bug vs implementation bug), iterate via SendMessage to 0.B/0.C | all cases ≤1e-9; both sides' internal asserts green |
| 0.E | OPUS panel ×3 | Full adversarial gate (ADR §7 policy): (1) math-vs-ADR critic (blocks, α formula, DOF map, TH basis = barycentric L_i not {L_i²}); (2) oracle-correctness + independence critic (derivations vs papers; confirm zero cross-contamination); (3) degeneracy/robustness critic (zero-area, inverted, ndm mix-ups, ld bugs, h-definition on slivers) | findings adjudicated + fixed; refuted list recorded |
| 0.F | MAIN | Header stamps, `stamp_headers.py` GLOBS, ledger row amend (P0 status), PR off `ladruno`, CI watch | PR merged |

Launch order: 0.A → {0.B ∥ 0.C} → 0.D → 0.E → 0.F.

### P0 interface contract (FROZEN — agents implement exactly this)

**Files/namespaces.** `SRC/element/ladrunoUP/LadrunoUPKernel.h` and
`LadrunoUPShapes.h`, header-only, `namespace ladruno_up` (shapes in
`ladruno_up::shapes`), zero OpenSees includes, style of
`SRC/element/ladrunoPlane/LadrunoFiniteStrain2DKernel.h` (raw doubles, runtime
sizes, leading-dimension params, caller-owned buffers, contractual guards
documented in the header comment).

**Layouts.** Row-major. Per GP the caller supplies: `Nu[a]` (a = 0…nNu−1),
`dNu[a*ndm + i]` = ∂N_a/∂x_i (cartesian, current=reference config, v1 linear
geometry), `Np[b]`, `dNp[b*ndm + i]`, and `dv` = w·detJ (2D: ×thickness, folded
by the CALLER). Element matrices assemble per-node-interleaved `[u…, p?]` via
the DOF map below.

**Block API** (accumulate into caller-zeroed buffers; `ld*` = leading dim):

```cpp
// Q[(a*ndm+i), b] += dNu[a*ndm+i] * biotAlpha * Np[b] * dv        (coupling)
void addQ(const double* dNu, const double* Np, int nNu, int nNp, int ndm,
          double biotAlpha, double dv, double* Q, int ldQ);
// H[b,c] += sum_i dNp[b*ndm+i] * kbar[i] * dNp[c*ndm+i] * dv      (Darcy, diag k̄)
void addH(const double* dNp, int nNp, int ndm, const double* kbar,
          double dv, double* H, int ldH);
// S[b,c] += Np[b] * oneOverQbar * Np[c] * dv                       (storage)
void addS(const double* Np, int nNp, double oneOverQbar, double dv,
          double* S, int ldS);
// Htilde: same form as addH with scalar stabAlpha (isotropic)      (stabilization)
void addHtilde(const double* dNp, int nNp, int ndm, double stabAlpha,
               double dv, double* Ht, int ldHt);
// fseep[b] += sum_i dNp[b*ndm+i] * kbar[i] * rhoF * drive[i] * dv  (seepage source;
// drive = b (static/body) or b - ütrial (dynamic seepage), built by caller)
void addFseep(const double* dNp, int nNp, int ndm, const double* kbar,
              double rhoF, const double* drive, double dv, double* fseep);
// alpha = alpha0 * h*h / (Ks + 4.0*Gs/3.0)                          (McGann auto)
double stabAlphaAuto(double h, double Ks, double Gs, double alpha0);
// storage coefficient: 1/Qbar = n/Kf + (alpha-n)/Ks  (Ks<=0 => infinite grains)
double oneOverQbar(double n, double Kf, double biotAlpha, double Ks);
// Per-node-interleaved DOF map. pCarrier[n] in {0,1}. Returns total DOF count.
// uOff[n] = first u-DOF slot of node n; pOff[n] = p slot or -1.
int dofMap(int nNodes, int ndm, const int* pCarrier, int* uOff, int* pOff);
```

Everything else (K from the material, M, C_Ray, block placement into the
OpenSees APIs) is ELEMENT-side (P1) per ADR §3.2 — not kernel scope.

**Providers** (`LadrunoUPShapes.h`): one struct per shape, static-const data +
stateless eval:

```cpp
struct Provider {
  int nNodes, ndm, nGP;
  const double* gpCoord;   // nGP × (ndm or barycentric dim, per shape doc)
  const double* gpWeight;
  const int* pCarrierTH;   // Taylor–Hood carrier flags (all-ones for linear shapes)
  // Nu, dNu (cartesian), detJ at gp from nodal coords (row-major xy[ndm*nNodes])
  void eval(int gp, const double* xy, double* Nu, double* dNu, double* detJ);
  // pressure basis: equal-order => same as Nu; TH => barycentric-linear L_i on
  // the vertex carriers (NOT the Bernstein vertex subset — ADR §3.3 ⟨UP-8⟩)
  void evalP(int gp, const double* xy, bool taylorHood, double* Np, double* dNp);
  double elemSizeH(const double* xy);  // largest EDGE length (ADR §3.3 "largest
                                       // element dimension"; 0.E-F1 adjudication —
                                       // side not diagonal on Q4/H8; == vertex-pair
                                       // max on simplices)
};
```

Shapes v1: `T3` (1-pt), `Q4` (2×2), `H8` (2×2×2, via the same trilinear
formulas as shared `shp3d.h`), `BT6` (Bernstein quadratic, 3-pt interior rule,
vertices carriers), `BTET10` (Bernstein, 4-pt rule, vertex carriers). Bézier
providers assume straight sides (affine map — the element-level guard arrives
P3; the provider documents the precondition). Donor code to lift (read, adapt,
cite in comments): LadrunoQuad.cpp:1321–1377, LadrunoCST.cpp:655–688,
`shp3d.h`, BezierTri6.cpp:851–944, BezierTet10.cpp:1177–1236.

**Cases file** (`tests/ladruno_up_cases.txt`, oracle-emitted, whitespace/line
format, full `%.17g` precision):

```
CASE <name> SHAPE <T3|Q4|H8|BT6|BTET10> PORDER <equal|th> NDM <n> NNODES <n>
COORDS <nNodes*ndm floats>
PARAMS <biotAlpha> <oneOverQbar> <rhoF> <stabAlpha> <thick> <kbar_1..kbar_ndm>
DRIVE <ndm floats>
DOFMAP <totalDof> <uOff...> <pOff...>
Q <rows> <cols> <floats row-major>
H <rows> <cols> <...>
S <rows> <cols> <...>
HT <rows> <cols> <...>
FSEEP <len> <...>
END
```

**Canonical cases (≥12):** Q4 unit square + trapezoid-distorted; T3 unit +
skewed; H8 unit cube + distorted; BT6 straight (auto midpoints) equal-order +
TH; BTET10 straight equal-order + TH; parameter spreads: α=1.0 and α=0.79
(finite K_s), isotropic and anisotropic k̄, thick=1.0 and 2.5 (2D). Oracle
computes blocks with an INDEPENDENT quadrature (own shape functions, degree-
overkill Gauss or exact affine integration) — not by re-implementing the
kernel's loop structure.

**Oracle-only gates** (block equality ⇒ spectra equality, so these live in
numpy): (i) 1-element consolidation ODE — sealed-sides Q4, drained top,
assemble [S+H̃, H, Qᵀ] from oracle blocks, closed-form pressure decay vs
`scipy.integrate` of the semi-discrete system; (ii) inf-sup smoke exactly per
ADR §7 P0 (4×4 patch, Schur S_p = Qᵀ·K⁻¹·Q (+H̃), eigencount < 1e-10·λ_max:
equal-order-no-stab ≥ 2, TH/stab = 1) — K from a plain-elastic B-matrix the
oracle builds itself; (iii) FD quantification of the dropped dynamic-seepage
tangent (ADR §3.1); (iv) symmetry (H, S, H̃) + pattern (Q sign layout) asserts.

### 0.A addendum — emission & convention pins (MAIN, 2026-07-10)

Adjudicated before agent launch; both P0 agents implement exactly this.

1. **Cases = contractual-rule evaluation, not overkill.** Emitted blocks in
   `ladruno_up_cases.txt` are integrated with the shape's CONTRACTUAL GP rule
   (T3 1-pt, Q4 2×2, H8 2×2×2, BT6 3-pt, BTET10 4-pt). Rationale: several
   integrands are not exactly integrated by those rules (T3 S-block is
   quadratic vs the 1-pt rule; distorted-Q4 cartesian gradients are rational;
   equal-order BT6 S-block is quartic vs the degree-2 rule), so a true-integral
   oracle would legitimately differ from the kernel far above 1e-9.
   Independence lives in the implementation (independent numpy basis/Jacobian/
   assembly code), not in the GP points, which are contract data. The oracle
   MUST additionally self-check contractual-rule blocks against degree-overkill
   or exact-affine integration for case/block combos where the contractual rule
   is provably exact (e.g. T3 Q/H/FSEEP blocks, rectangular/parallelogram Q4,
   TH-linear blocks on straight Bézier shapes), assert ≤1e-12 there, and
   document known-inexact combos as expected quadrature deviation, not error.
2. **Pinned GP rules** (kernel agent verifies the donor code agrees and reports
   any discrepancy to MAIN rather than silently diverging):
   T3 = 1 pt at barycentric (⅓,⅓), w=½. Q4 = 2×2 Gauss (±1/√3, w=1).
   H8 = 2×2×2 Gauss. BT6 = 3-pt interior (ξ₁,ξ₂) ∈ {(⅙,⅙),(⅔,⅙),(⅙,⅔)},
   w=⅙ each (donor `BezierTri6.cpp` GP3). BTET10 = 4-pt, (L1,L2,L3) = perms of
   a=0.585410196624968 / b=0.138196601125011, w=1/24 each (donor
   `BezierTet10.cpp` GP4).
3. **Pinned basis order.** BT6 (ξ₃=1−ξ₁−ξ₂): N0=ξ₃², N1=ξ₁², N2=ξ₂²,
   N3=2ξ₁ξ₃, N4=2ξ₁ξ₂, N5=2ξ₂ξ₃ — vertices are nodes 0,1,2 (node0↔ξ₃,
   node1↔ξ₁, node2↔ξ₂). BTET10 (L4=1−L1−L2−L3): N0..N3 = L1²..L4²,
   N4=2L1L2, N5=2L2L3, N6=2L1L3, N7=2L1L4, N8=2L3L4, N9=2L2L4 — vertices are
   nodes 0–3. T3/Q4/H8: standard CCW orderings per the donor elements.
4. **TH pressure compaction.** `Np`/`dNp` are COMPACTED over carriers:
   nNp = #carriers, ordered by element node order restricted to carriers
   (BT6 TH: Np = [ξ₃, ξ₁, ξ₂]; BTET10 TH: Np = [L1, L2, L3, L4] per the
   vertex↔coordinate map above). Q/H/S/HT/FSEEP dims use this nNp; the case
   DOFMAP (pOff = −1 on non-carriers) carries assembly placement.
5. **dv folding.** Case harness computes dv = w·detJ·thick (2D) or w·detJ (3D)
   with `thick` from PARAMS; 3D cases set thick = 1.0.

### P0 exit = ADR §7 P0 gate row, verbatim.

---

## P1 — LadrunoUP element, Q4 lane end-to-end

Branch after P0 merges. Interface pinning (WP1.A, MAIN) happens in this doc
before agents launch; the ADR already fixes the hard parts (§3.2 contract table
+ obligations, §4.1 surface + legality matrix + loads/defaults).

| WP | Owner | Deliverable |
|---|---|---|
| 1.A | MAIN | Pin ctor signature, member layout, response IDs; refresh this doc |
| 1.B | OPUS | `LadrunoUP.{h,cpp}`: the six §3.2 APIs (incl. `getResistingForceIncInertia`, `getInitialStiff` [K₀,−Q;0,H] + cache-dirty, `addInertiaLoadToUnbalance` p-slots zeroed), solid-only Rayleigh (override `getDamp`/`getRayleighDampingForces`; never base helpers), `setParameter` material forwarding + self-co-registration (stage flip dirties α/K₀/H̃), `update()` (strain-only; p never enters material), commit/revert/revertToStart, zeroLoad/addLoad (SelfWeight replaces scaled solid+seepage values; others rejected), setResponse (`stresses`, `stressesTotal`, `porePressure`, `flux`), per-instance buffers (NO class statics) |
| 1.C | OPUS | `OPS_LadrunoUP.cpp` (legality matrix, FATAL unknown flags, `-permH/-gammaW` sugar, defaults, solver notice, `-stab auto <α0>` w/ elemSizeH + initial-tangent moduli) + registration: OpenSeesElementCommands + TclElementCommands + FEM_ObjectBrokerAllClasses + CMake + `classTags.h` #define (cross-registry comment; fix stale :943 note) + stamps + vanilla-ledger rows |
| 1.D | OPUS | `sendSelf/recvSelf` (EVERY ctor arg incl. -lumped/α₀/dynSeepage) + MP smoke test |
| 1.E | OPUS | `Ladruno_NodeResults.cpp` PressureSource contract-aware (disp-slot for LadrunoUP nodes, vel-slot upstream) ⟨FW-F4⟩ |
| 1.F | OPUS ×2 ∥ | Zone-B battery, split: **(i) analytic** — Terzaghi vs series, B2 ZS84 column (~1e-3), B3 Boone–Ingraffea hard numbers, static steady-seepage NEW-capability, all-impervious-static loud-singularity smoke, wrong-solver ProfileSPD-vs-UmfPack divergence doc + LEDGER_quirks row, Rayleigh contract test, FD consistency incl. IncInertia; **(ii) equivalence + stabilization** — two-leg quadUP equivalence (γ=½/β=¼ tight ≤1e-6 w/ `-dynSeepage off` + unbalance-based tests; γ=0.6 Δt-halving mutual convergence), SSPquadUP pressure-block cross-check, B4 CB-index + α₀ sweep, bbar behavior gate, sendSelf round-trip |
| 1.G | MAIN | Build (`cmd /c build.bat OpenSees`), run battery (`python -S`), integrate, ledger/quirks updates, PR, CI watch |

No separate adversarial panel at P1 (ADR §7 policy — the battery carries it);
MAIN does a focused self-review of 1.B's Rayleigh/IncInertia code against the
§3.2 obligations before the PR.

### WP1.A pins (FROZEN 2026-07-10 — agents implement exactly this)

**Scope guard.** P1 = the unified element with the Q4 lane gated end-to-end.
The element is shape-generic over the P0 providers, but the PARSER at P1
rejects nNodes ∈ {6,10} with a named error ("Bézier Taylor–Hood lanes land at
P3 — heterogeneous-ndf plumbing not yet built"); T3/H8 parse and run (their
gates arrive P2). All nodes must have ndf = ndm+1 at P1 (equal-order only);
`setDomain` errors loudly otherwise.

**Files/ownership.** 1.B owns `SRC/element/ladrunoUP/LadrunoUP.{h,cpp}` ONLY.
1.C owns `SRC/element/ladrunoUP/OPS_LadrunoUP.cpp` + the registration edits
(`classTags.h`, `OpenSeesElementCommands.cpp`, `TclElementCommands.cpp`,
`FEM_ObjectBrokerAllClasses.cpp`, `ladrunoUP/CMakeLists.txt`) + vanilla-ledger
rows. 1.D (wave 2) owns the `sendSelf/recvSelf` bodies inside LadrunoUP.cpp +
`tests/test_ladruno_up_mp_smoke.py`. 1.E owns `Ladruno_NodeResults.cpp` only.
1.F owns `tests/test_ladruno_up_element_analytic.py` (i) and
`tests/test_ladruno_up_element_equiv.py` (ii).

**Ctor signature (C++, pinned):**

```cpp
LadrunoUP(int tag, int ndm, const ID& nodeTags,        // 3|4|8 nodes at P1
          NDMaterial& theMaterial,                     // copied per GP via getCopy(ndm==2?"PlaneStrain":"ThreeDimensional")
          double thickness,                            // 2D only; 1.0 in 3D
          double Kf, double poro, double rhoF,
          const Vector& perm,                          // ndm entries, k_hydraulic/gammaW
          double biotAlpha, double Ks,                 // Ks<=0 => infinite grains
          const Vector& body, const Vector& fluidBody, // ndm entries each
          int formulation,                             // 0=std, 1=bbar
          int pOrder,                                  // 0=equal (P1); 1=TH reserved P3
          bool lumped,
          int stabMode, double stabValue,              // 0=off,1=auto(stabValue=alpha0),2=manual(stabValue=alpha)
          bool dynSeepage);                            // default true
```

Default ctor `LadrunoUP()` for the broker (tag 0, null members) as family idiom.

**Member layout (pinned, per-instance — NO class statics, ADR-40):** provider
selected by (ndm,nNodes) into a small tagged union/switch (no virtual provider
class — the P0 structs are stateless); `NDMaterial** theMaterials[nGP]`;
per-instance `Matrix K_, damp_, mass_` sized totalDof×totalDof at setDomain via
kernel `dofMap` + persistent `Vector resid_`; committed-solid-K copy
(`Matrix KcSolid_`) maintained at `commitState` for βKc Rayleigh; caches:
`Matrix K0_` + `bool k0Dirty_`, `double stabAlpha_` + `bool stabDirty_`
(both dirtied by `updateMaterialStage`/perm `setParameter`); GP scratch arrays
(Nu/dNu/Np/dNp/dv) as per-instance std::vector<double>, sized once in
setDomain. uOff/pOff int vectors from kernel dofMap.

**DOF order:** per-node-interleaved [ux uy (uz) p] — identical to the P0
kernel dofMap. p slot = index ndm within each node block.

**Response IDs (setResponse, pinned):** 1 = `stresses` (effective, per-GP,
material stress vector); 2 = `stressesTotal` (σ_total = σ′ − α·p on the normal
components per ADR §3.1 — tension-positive σ′, compression-positive p; pin
wording corrected at 1.G, agent flag adjudicated); 3 = `porePressure` (per-GP scalar, Np·p_nodal);
4 = `flux` (per-GP Darcy flux vector, −k̄·(∇p − ρ_f(b_f − ü)) — ndm comps);
plus `material $gpNum …` forwarding (family idiom). IDs 1–4 chosen to avoid
the plane-family's 21 (`stressZZ`) namespace; anything else → null response.

**The six §3.2 APIs (contract table verbatim from the ADR):** getTangentStiff
[K,−Q;0,H]; getInitialStiff [K₀,−Q;0,H] cached+dirty; getDamp [C_Ray,0;Qᵀ,S+H̃]
solid-only Rayleigh (own committed copy; never base helpers); getMass [M,0;0,0]
(consistent or HRZ row-sum if `-lumped`, solid block only);
getResistingForce [∫BᵀσdV − Q·p − body; H·p − f_seep] (NO rate terms);
getResistingForceIncInertia = resisting + M·ü + getDamp()·v̇-terms (rates from
trial vel/accel); addInertiaLoadToUnbalance −M·(R·üg), p-slots zeroed;
getRayleighDampingForces override (solid-only). `update()` = strains→
setTrialStrain only, p never enters the material. zeroLoad/addLoad:
LOAD_TAG_SelfWeight REPLACES scaled body+fluidBody; all other elemental loads
rejected with named error; zeroLoad restores ctor values.

**Parser (1.C) pins:** surface per ADR §4.1 verbatim (incl. `-permH/-gammaW`
sugar, `-stab auto <alpha0>` default-on for equal-order with one-line notice,
unknown-flag-FATAL, solver notice at creation, `-geom linear` only-value,
`-dynSeepage on|off`). classTags.h: `#define ELE_TAG_LadrunoUP 33017` with
cross-registry comment; FIX the stale line-943 "33016–33019 free" note.
`-stab` on TH → fatal (moot at P1, TH rejected anyway — keep the check for P3).
stabAlphaAuto consumes h = provider elemSizeH (largest edge) + skeleton
moduli from the MATERIAL INITIAL TANGENT (probe D0: Gs from D0 shear diag,
Ks from λ+2G/3 — document extraction in code).

**1.E pin:** `Ladruno_NodeResults.cpp` PressureSource: if the node carries a
LadrunoUP element (query via node's element list class tags == 33017) read
disp-slot ndm; else legacy vel-slot ndm. Flag-free, automatic. ⟨FW-F4⟩.

---

## P2–P4 sketch (each = own branch/PR; agents pinned at phase start)

### WP3.A pins (FROZEN 2026-07-11 — P3 agents implement exactly this)

P3 = Taylor-Hood on Bezier BT6/BTET10 (vertex-p, heterogeneous ndf). Branch
`feature/adr71-p3`. Full adversarial panel #2 runs before the PR (ADR policy).

**Ownership.** 3.B owns SRC/element/ladrunoUP/LadrunoUP.{h,cpp} +
OPS_LadrunoUP.cpp (all three, one agent — the changes interlock). 3.C-B1 owns
tests/test_ladruno_up_element_th_b1.py. 3.C-PL owns
tests/test_ladruno_up_element_th.py. Waves: 3.B -> MAIN build -> {3.C-B1 ||
3.C-PL} -> panel x3 -> PR.

**3.B pins:**
1. Parser flip: (2,6)/(3,10) now ACCEPTED and REQUIRE `-pOrder linear`
   (omitted or `equal` on a quadratic shape = FATAL naming the reserved
   equal-order-quadratic axis). `-stab` on TH stays fatal. Linear shapes
   unchanged.
2. Heterogeneous-ndf validation (setDomain, STRICT): vertex/carrier nodes
   ndf == ndm+1 exactly; mid-edge nodes ndf == ndm exactly; anything else =
   loud error naming node + expectation. (Rationale: a uniform-ndf model
   would leave floating mid-edge p-DOFs -> silent singularity; the fork
   philosophy is loud. The equalDOF mixed-ndf example in 3.C-PL shows the
   modeling pattern.)
3. TH assembly: dofMap fed the provider pCarrierTH (pOrder==1); nP_ =
   #carriers; evalP(taylorHood=true); Htilde NEVER assembled on TH.
4. Straight-side guard (setDomain, BT6+BTET10, both pOrders): every mid-edge
   node must satisfy ||x_m - (x_a+x_b)/2|| <= 1e-6 * ||x_b-x_a|| for its
   edge (a,b); violation = loud error naming element/node/distance.
   (Provider affine-map precondition, ADR section 3.3.)
5. **BTET10 winding ADJUDICATION (the P0 trap, decided):** the pinned node<->L
   map makes conventionally RIGHT-handed tets evaluate detJ < 0. The element
   accepts BOTH windings: at setDomain, if all GPs have detJ < 0 the element
   folds dv = w*|detJ| (gradients from the signed J are orientation-correct;
   the measure must be positive) and prints a ONE-TIME informational note;
   all-positive proceeds as-is; mixed sign or |detJ| ~ 0 = loud error.
   T3/Q4/H8/BT6 keep the P1 detJ>0-per-GP rejection (their standard
   CCW/right-handed winding IS map-positive).
6. recvSelf/configureSizing must reproduce the TH sizing (pCarrier path) —
   extend the existing helper, no hand-copies. sendSelf already ships pOrder.
7. Responses on TH: porePressure/flux use the carrier-compacted Np (nP_
   values); stresses/stressesTotal unchanged (u-side quadratic).

**3.C-B1 pins (B1 ZCB80 phasor gate, ADR section 7.1 verbatim):** BT6-TH
column 30 m x 1 m (crossed or strip mesh of BT6 — agent picks the cleanest
straight-sided triangulation and documents it), q(t)=100*sin(3.379t) Pa top,
drained top, rigid impermeable base, E=7.492e8 Pa, nu=0.2, n=0.333,
k'=1e-7 m/s (kbar = k'/(rho_w*g)), rho_s=2000, rho_w=1000; TWO legs:
Qbar=1e4 MPa (compressible; 1/Qbar via equivalent Kf = n*Qbar) and Qbar=1e9
MPa (undrained-limit checkerboard stressor). FE realization per ADR: ring-up
until per-cycle amplitude drift < 0.5%, LS-fit A*sin+B*cos per node over the
last cycle, compare |p(z)|, arg p(z), |u(z)| vs the ZCB80 closed form
(eqs 12-30 — implement independently in numpy from the Book/paper formulas;
document the equations used). Gates: Q1e4 leg L2-rel <= 2%; Q1e9 leg:
no-checkerboard (CB metric ~ 0 on the p field) + pressure L2 mesh-refinement
convergence (2 refinements) + quadratic-u rate preserved (u L2 error order
~3 or at least >2 on refinement); one dt-halving check pins dt <= T/100
adequacy. gamma=0.6/beta=0.3025.

**3.C-PL pins (TH plumbing battery):** (i) sendSelf/recvSelf DB round-trip
on a BT6-TH AND a BTET10-TH element (non-default args; mirror the 1.D serial
gate; exact restore + re-solve); (ii) recorders: node disp-slot p on
carriers, eleResponse porePressure/flux lengths = nP_/carrier-consistent,
PressureSource disp-slot on TH carriers; (iii) numberers: same tiny TH model
solved under Plain vs RCM numberer -> identical results (mixed-ndf
renumbering proof); (iv) equalDOF mixed-ndf interface example: dry region
(ndf=2 solid elements) + saturated BT6-TH region sharing an interface,
explicit-DOF-list equalDOF on shared u-DOFs (ADR FW-F10 — bare form
mis-sizes), static+transient smoke, THE example the P4 guide will cite;
(v) BTET10 winding acceptance: same one-element problem right-handed vs
left-handed input -> identical p/u (1e-12) + the one-time note observed;
(vi) straight-side guard: curved mid-edge node -> loud error (assert message,
not crash); (vii) BT6-TH static steady seepage exactness (P1 gate 4 analog);
(viii) BTET10-TH consolidation smoke vs series (loose 5%).

**Run mechanics:** hermetic bootstrap donors as P2; UmfPack; batteries
target <= 120 s each. MAIN rebuilds after 3.B before batteries run.

### WP2.A pins (FROZEN 2026-07-11 — P2 agents implement exactly this)

P2 = the H8 and T3 lanes GATED (both already parse+run since P1 — this phase
is battery work; element-source edits are NOT in agent scope: bugs get
REPORTED to MAIN for adjudication). Branch `feature/adr71-p2`.

**Ownership.** OPUS-H8 owns `tests/test_ladruno_up_element_h8.py` ONLY.
OPUS-T3 owns `tests/test_ladruno_up_element_t3.py` ONLY. Style/bootstrap donor
for both: `tests/test_ladruno_up_element_analytic.py` + `_equiv.py` (hermetic
py -3.12 bootstrap, printA order='F', zone_b markers, tolerances rationale).

**OPUS-H8 gates (ADR §7 P2 row):**
1. brickUP two-leg equivalence, P1 methodology: leg 1 tight (γ=½/β=¼,
   `-dynSeepage off -stab off`, consistent mass, NormUnbalance, ≤1e-6 —
   expect machine-level) + leg 2 production (γ=0.6/β=0.3025, Δt-halving
   mutual first-order convergence). Upstream = `brickUP` (bulk = pre-combined
   Q̄ ≈ Kf/n like quadUP — agent verifies the arg mapping in upstream source
   and pins it in the docstring). PLUS the `-lumped` leg both ways (inspect
   how upstream brickUP lumps; if it has no lumped option, run OUR -lumped
   vs consistent as a Δt-convergence pair instead and document).
2. 3D Terzaghi: H8 1×1×10 column vs the series (Tv sweep as P1, 3D BCs:
   sealed sides via no-p-fix, drained top).
3. B4-3D (McGann 2015 cube-footing analog): 30 m cube quarter-model, 10³
   mesh h=3 m, strip→square footing load ramped+held, drained top only,
   CB_lap 3D (face-interior Laplacian roughness) + α₀ ∈ {off, 0.25}:
   monotone CB suppression + auto-α twin-run identity (P1 trick), α value
   = 0.25·9/(Ks+4Gs/3) pinned ±5%. Keep runtime sane (≤10³ mesh, ≤60 s).
**OPUS-T3 gates (ADR §7 P2 row ⟨scope-F12⟩):**
1. B4-T3: the P1 B4 footing re-meshed with CROSSED-diagonal T3 split (each
   quad → 4 triangles, center node), CB gate + α₀-sweep {off, 0.05, 0.25,
   0.5} — pins that the §3.3 α transfer to simplices holds with h = largest
   triangle edge; auto-vs-manual twin-run identity.
2. T3 honest-baseline documentation gate: measured undrained locking —
   drained vs undrained footing settlement T3-mesh vs Q4-mesh reference;
   ASSERT the direction (T3 stiffer, worse under undrained constraint,
   inherited CST pathology) and PIN the measured ratios in the docstring
   (documentation assert, generous bands).
3. T3 Terzaghi consolidation smoke on the crossed mesh (loose gates — the
   point is no-blowup + right decay class, not accuracy records).

**Run mechanics (both):** existing `dist/bin` pyd is current for this branch
(tip = #557 squash; NO rebuild needed unless MAIN says so). Same bootstrap
as the P1 batteries. Battery target runtime ≤90 s per file.


| Phase | Notable WPs (owner) |
|---|---|
| **P2** H8+T3 | providers already exist → element lanes + gates (OPUS); B4-3D + B4-T3 crossed-diagonal CB gates (OPUS); brickUP two-leg equivalence w/ `-lumped` (OPUS); 3D Terzaghi (OPUS) |
| **P3** TH on Bézier | heterogeneous-ndf DOF plumbing in element + setDomain vertex/mid-edge validation + straight-side guard (OPUS); B1 ZCB80 phasor gate (ring-up + LS fit per §7.1) (OPUS); recorders/numberers proving (OPUS); equalDOF mixed-ndf example for the guide (OPUS); **full OPUS adversarial panel #2** (ADR policy) |
| **P4** dynamics + ecosystem | **B5 Simon oracle** (OPUS — §7.1 spec + the FOUR errata; û(0,τ) measured-first protocol); PDMY liquefaction column vs quadUP (OPUS); both-Newmark-sets oscillation pinning (OPUS); gravity/hydrostatic init recipe + `LadrunoUP_guide.md` draft (OPUS, MAIN polish); PressureSource/Monitor gates; banner line + `patch_banner.py` + ledger status→shipped (MAIN) |
| P5–P7 | demand/reserved/research per ADR — not planned here |

## Standing notes

- Cost posture: OPUS for production + panels (user directive 2026-07-10);
  MAIN reserves itself for contracts, adjudication, build/PR mechanics.
- Every OPUS production agent gets: the ADR, this doc's relevant contract
  section, the named donor files, and the explicit list of files it OWNS
  (disjoint ownership — no two agents touch one file in the same wave).
- Agents report structured summaries; MAIN never merges unverified agent output
  — 0.D-style cross-checks or the phase battery must pass first.

## Log

- 2026-07-10 — plan created (P0 contract frozen; ADR locked via PR #547).
- 2026-07-10 — **P0 SHIPPED** (PR #551 merged, squash b2ac3ba04). WP1.A pins frozen; P1 branch `feature/adr71-p1`.
- 2026-07-11 — **P1 EXECUTED** (waves 1+2, five Opus agents). Battery: analytic
  10 pass + 1 documented xfail; equivalence 9/9; MP smoke + serial DB
  round-trip green. Highlights: quadUP two-leg equivalence machine-identical
  at γ=½/β=¼ (3.6e-15); B4 checkerboard monotone suppression, auto-α within
  1.7% of the ADR pin; Rayleigh contract bit-exact; FD tangent 1e-10 class.
  **MAIN adjudications:**
  - `stressesTotal` sign = ADR §3.1 (σ_total = σ′ − α·p on normals); pin
    wording was loose, corrected.
  - **ADR "loud singularity" claim REFUTED** (1.F-i finding 5): all-impervious
    static returns rc=0 with a consistent-Neumann arbitrary p level
    (solver-dependent); pinned as strict xfail; guide/P4 must say "silently
    arbitrary p level", not "loud failure". σ_min/σ_max < 1e-12 confirmed.
  - **`-dynSeepage on` diverges under Δt-refinement in QUASI-static
    consolidation** (noise-fed ü term): default stays `on` per LOCKED ADR;
    quirks row + guide rule "consolidation runs use -dynSeepage off"; P4 B5
    revisits the default.
  - **Transformation + staged nonzero p-`sp` = silently wrong steady state**
    (Penalty/Lagrange correct): quirks row; init recipe (P4) must specify
    Penalty.
  - `printA -ret` is COLUMN-major: quirks row (false transposed-Q alarms).
  - α₀ out-of-range warning made once-per-process (was per-element spam).
  - 1.D `configureSizing()` ctor/recvSelf shared helper: APPROVED (anti-drift).
  - Rayleigh shadow of non-virtual getRayleighDampingForces: verified safe —
    zero external callers in SRC/domain + SRC/analysis (element-internal only).
- 2026-07-11 — **P3 EXECUTED** (3.B element TH enablement + 2 battery agents +
  MAIN integration fix: builder-NDF gate deferred past shape detection — the
  TH modeling dance legally leaves the builder at ndf=ndm). Suite 50+1xfail.
  B1 ZCB80 gate: independent u-p closed form (per-mode exponential
  referencing), graded BT6-TH column — Q̄=1e4 leg |p̂| L2 1.47% / |û| 0.06%;
  Q̄=1e9 leg no-checkerboard + monotone convergence; Bernstein edge-thirds
  consistent loading; uniform meshes CANNOT pass B1 (5 cm boundary layer —
  graded mesh is load-bearing).
- 2026-07-11 — **P3 adversarial panel #2 (3 Opus critics): ALL PASS, zero
  defects.** Highlights: C1 re-derived the dispersion quadratic + independent
  FD solver confirming the B1 oracle (<0.5%); winding-fold reasoning airtight
  (only the measure folds; signed-J gradients are true spatial gradients);
  C2 — the ADR-70 βK-Rayleigh P-clobber class is STRUCTURALLY ABSENT
  (distinct resid_/rayForce_, experiment 3.2e-30) and the |detJ| fold is
  measure-exact (ΣRz = ρgV at 0 relerr on a folded RH tet); C3 — three-route
  oracle cross-check table all ✓, element-level TH transient blocks correct
  (honest-p transpose Cpu=(−Kup)ᵀ at 6.2e-8). **Adjudications:**
  - "quadratic-u rate > 2" pin REFUTED for B1 (u inherits the pressure
    boundary-layer rate ~1.0–1.2 through coupling); replaced by the STRONGER
    parabola-exactness proof (1.5e-14; exact P2 reproduction) — C1 concurs;
    ADR/guide rewording lands with the P4 guide.
  - porePressure response = per-GP (matches ADR §4.4 "porePressure (GP)");
    WP1.A pin wording corrected; note for future shapes where nGP ≠ nP.
  - B1 density literalism (ρ=2000 as mixture, ADR lists grain+water):
    self-consistent, gated quantities ρ-insensitive; one-line note for the
    P4 guide when citing B1.
  - Panel-recommended hardening APPLIED (3.C-PL amendment): TH transient-
    block FD gate, -dynSeepage-on default-config smoke, porePressure value
    gate on the linear-field patch.
  - Guide notes for P4: -stab off is fatal on TH (arg lists can't be shared
    across mixed meshes); re-recv different-shape material-leak wrinkle
    (unreachable via broker path); bbar-on-TH legal but ungated (formulation
    axis is P1/P2 scope).
- 2026-07-11 — **P1 SHIPPED** (PR #557 merged, squash 94dcf9b8d; manifest-row
  CI gate learned: every new ELE_TAG #define needs a testbed/manifest.yaml
  row). P2 branch `feature/adr71-p2`; WP2.A pins frozen.
- 2026-07-11 — **P2 EXECUTED** (2 Opus agents, H8 + T3 lanes; 14 new tests,
  full suite 36 pass + 1 xfail). H8: brickUP tight-leg equivalence 1.1e-15,
  3D Terzaghi in the P1 ladder, B4-3D CB ratio 3.20 + auto-α twin bit-exact
  (6.686e-5 = ADR pin ±2%). T3: crossed B4 CB_cell suppression 8.23 monotone
  + α-transfer-to-simplices pinned (h = largest edge = 3 m), Terzaghi smoke
  <0.1%. **MAIN adjudications:**
  - **brickUP `-lumped` also row-sum-lumps S** (BrickUP.cpp:916); ours lumps
    solid mass only (ADR §3.2 "solid block only") ⇒ lumped legs gated as a
    Δt-convergence pair (1.2e-6 finest), not operator identity. Contract
    difference, not a bug — documented in the H8 test docstring.
  - **Crossed-T3 spurious mode is CENTER-NODE-localized**, not the corner
    checkerboard: corner-lattice CB stays smooth (off/0.25 = 1.34,
    non-monotone) while the cell-center metric CB_cell shows the instability
    (0.80 → 0.10, ratio 8.23). WP2.A's corner-lattice gate literalism was
    wrong for this topology; CB_cell APPROVED as the primary suppression
    gate, corner metric computed + documented for Q4 comparability.
  - **ADR ⟨scope-F12⟩ locking direction REFUTED vs std-Q4**: crossed-T3
    locks LESS than full-integration Q4 (T3x/Q4std settlement ratio 0.998 →
    1.262 as ν→0.49; the union-jack center bubble relieves dilatation;
    Q4std/Q4bbar itself drops to 0.751). The pin's spirit holds only vs the
    locking-free bbar reference (0.960 → 0.948). ADR/guide wording at P4:
    "T3 crossed locks, but LESS than full-integration Q4"; honest-baseline
    framing survives, the vs-Q4 direction does not.

- 2026-07-10 — 0.A addendum committed (GP/basis pins, contractual-rule emission).
- 2026-07-10 — 0.B/0.C landed (twin Opus agents); 0.D cross-run green first try:
  16 cases, 96/96 blocks + 16/16 dofMaps ≤1e-16 (gate 1e-9).
- 2026-07-10 — 0.E panel (3 Opus critics, all PASS-WITH-FIXES; independence
  verdict INDEPENDENT with content-level evidence). **Adjudications:**
  - **h-definition (C1-F1, MAJOR)**: plan contract's "largest vertex-pair
    distance" mis-cited ADR §3.3 "largest element dimension" (B4 calibration
    ⇒ side, not diagonal). FIXED both sides: h = largest EDGE length.
    Correlated-error class the twin-agent gate can't see — caught by panel.
  - **BTET10 case orientation (C2-A1, MAJOR)**: canonical tets negatively
    oriented (detJ<0, inverted blocks emitted). FIXED + permanent detJ>0 and
    S/H-PSD asserts added. Signed-detJ fold (addendum §5) STANDS for P0;
    orientation legality is element-side (P1/P3).
  - **TH inf-sup leg missing (C1-F2 = C2-A5, MAJOR)**: ADR §7 P0 demands
    per-pair eigencounts incl. TH = exactly 1. FIXED: generic patch assembler,
    all equal-order pairs (no-stab ≥2 / stab =1) + BT6-TH/BTET10-TH (=1).
  - **Harness NaN-blindness (C3-F1, MAJOR)**: comparator accepted all-NaN
    blocks as "ok". FIXED + per-block relative tolerance (C1-F3) + parser
    hardening (C3-F2/F3) + case-count assert in pytest.
  - **Coverage gaps (C2-B3, MAJOR/MINOR)**: ρ_f≡1 and vertical-only drive in
    every case; aniso k̄ never on distorted geometry. FIXED: new cases.
  - Doc/nit batch: stabAlphaAuto precondition + KsDrained/GsDrained rename
    (grain-vs-skeleton Ks collision, C1-F5/C3-F4), oneOverQbar domain,
    zero-gradient wording, C++17 pin, consolidation-gate honesty (C2-A3:
    structural sanity check, not a block discriminator), overkill-gate
    missed-exact combos, hydrostatic H·p=f_seep sign assert.
  - **P1/P3 TRAP (record, no P0 action)**: under the pinned BTET10 node↔L map,
    detJ = det[v0−v3, v1−v3, v2−v3] = −det[v1−v0, v2−v0, v3−v0] — positive
    detJ requires the conventionally LEFT-handed vertex winding. A P1/P3
    element-side detJ>0 guard would reject standard right-handed tets; the
    element must either normalize winding, use |detJ| for the measure (donor
    O11 convention), or define its guard against the map-consistent sign.
    Decide at P3 (BTET10 element lane). Cross-check unaffected (both sides
    share the map; canonical cases re-wound positive).
  - **Refuted / no-action**: inf-sup count 8 on Q4 patch = mathematically
    expected (rank argument verified, C2-A4), not an assembly bug; D-matrix
    choice immaterial to null count; sliver over-stabilization is the ADR's
    deliberate pinned definition (C3 attack 5); inverted-element gradient
    behavior better-defined than the BezierTri6 donor's (C3 attack 2);
    collapsed-quad-to-triangle is legitimate usage, not a defect; ld-semantics
    and dofMap survived hostile-input experiments (C3 attacks 3/4).
