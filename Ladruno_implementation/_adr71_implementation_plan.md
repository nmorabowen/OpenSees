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

---

## P2–P4 sketch (each = own branch/PR; agents pinned at phase start)

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
