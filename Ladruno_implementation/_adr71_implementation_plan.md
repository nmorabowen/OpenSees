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
  double elemSizeH(const double* xy);  // largest vertex-pair distance (ADR §3.3)
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
