# 21 — RC-3D Pre-Build Validation Gates (Gate 1 + Gate 2)

**Status:** RUN — both gates executed 2026-06-03 (build `f6700775`, `stdBrick` +
`ASDConcrete3D`, openseespy 3.12). Decisive results below.
**Why:** the adversarial review of [[20_ladruno_embedded_reinforcement_adr]] (workflow
`wf_32a898cd`) demanded two empirical gates *before* any embedded-element C++ — one to
decide whether the element is even needed (conformal escape hatch), one to settle the
disputed "emergent confinement" blocker. Scripts: `rc3d_gates/gate1_conformal_column.py`,
`rc3d_gates/gate2_concrete_confinement.py`. Run with the 3.12 build:
`pythoncore-3.12-64/python.exe gateN_*.py` (auto-bootstraps `dist/bin`).

---

## Gate 1 — is conformal meshing viable (so the embedded element is only for hard cases)?

**Method.** Square column (b×b×H), n×n×nz `stdBrick` concrete (ASDConcrete3D), vertical
`Steel02` truss rebar along every perimeter vertical **node-line** → rebar shares concrete
nodes ⇒ perfect bond with **NO embedding element and NO new C++**. Constant axial load
(N/fcAg = 0.1) then lateral displacement-control pushover. The structural point: with
shared nodes a bar can only sit on a mesh node-line, so **perimeter bars = 4·n** — bar
count is locked to the cross-section mesh.

**Result** (total steel held CONSTANT at `ρ·B² = 0.00160 m²` for every n — `a_bar =
ρB²/nbar` — so this isolates mesh/placement, it does **not** vary steel):

| n | bars | elems | ndof | converged | drift | V_peak [kN] |
|---|------|-------|------|-----------|-------|-------------|
| 1 | 4 | 25 | 72 | yes | 0.020 | 128.1 |
| 2 | 8 | 60 | 162 | yes | 0.020 | 136.8 |
| 3 | 12 | 105 | 288 | yes | 0.020 | 129.2 |

**Verdict — conformal is VIABLE for regular grid-aligned cages, zero new code.** All three
meshes converge (implicit, KrylovNewton + step-halving) through a sane axial+pushover.
`V_peak` is ~flat (≈128–137 kN); at constant steel the ~7% spread is mesh refinement +
corner→perimeter bar redistribution, **not** more steel (an earlier draft wrongly said
"capacity rises with steel" — corrected).

**What is MEASURED vs INFERRED:**
- *Measured:* conformal converges for 4, 8, **and 12** longitudinal bars on a regular grid.
  Pure bar count is **not** the wall — 12 grid-bars work fine.
- *Inferred (not yet run):* conformal breaks on (a) bars at **non-grid positions** /
  non-uniform spacing (they must land on mesh node-lines), (b) **ties/hoops** forcing
  combinatorial horizontal meshing, and (c) the **O(N²·nz) mesh-cost** blow-up
  (4→25, 8→60, 12→105 concrete elems). These are reasoning, not measurement.

**So:** conformal = the v1 path for regular grid-aligned columns/beams/walls; the embedded
element (ADR 20) earns its keep for **non-grid cages, hoops, and beam-column joints** —
which still wants a hoop/non-grid demonstration before it's a measured boundary. Explicit
conformal run is also a documented follow-up (mass + `dt_cr`), not yet exercised.

## Gate 2 — does ASDConcrete3D produce EMERGENT confinement (no fc inflation)?

**Method.** Single `stdBrick` + ASDConcrete3D under **constant lateral confining pressure
p** (held via `loadConst`) + axial displacement control; symmetry rollers; +x/+y faces kept
planar by `equalDOF`. Peak axial stress = confined strength `fcc(p)`, vs Mander (calibrated
to the same unconfined backbone) and Richart `fcc = fc + 4.1p`. Backbone calibrated so the
**unconfined** peak == fc (first compression point = elastic limit at 0.4fc, d=0; peak −fc
at −0.002).

**Result.**

| p/fc | fcc_sim | Mander | sim/Mander | ε_peak |
|------|---------|--------|------------|--------|
| 0.00 | 30.00 | 30.00 | **1.00** | 0.0020 |
| 0.05 | 37.24 | 39.30 | 0.95 | 0.0025 |
| 0.10 | 44.50 | 46.95 | 0.95 | 0.0029 |
| 0.15 | 51.71 | 53.47 | 0.97 | 0.0034 |
| 0.20 | 58.96 | 59.16 | **1.00** | 0.0039 |

**Verdict — PASS.** Unconfined recovers fc exactly; confined `fcc` tracks **Mander within
5%** across p/fc∈[0,0.20]; peak strain grows 0.0020→0.0039 (the **ductility** gain too).
**This refutes the adversarial review's #1 blocker** ("ASDConcrete3D is scalar-damage, can't
confine"). Confinement *does* self-emerge via the Lubliner Kc triaxial surface. The narrower
true caveat stands: there is **no tunable dilation-angle input**, so the *amount* of
confinement is governed by Kc + the compression backbone, not a dilation knob — but it does
not need pre-inflated fc, and it is quantitatively close to Mander here. **NB:** plain Newton
+ fixed steps **diverged** at high p; convergence needed **KrylovNewton + recursive
step-halving** — direct evidence for the blessed adaptive-stepping wrapper ([[19_rc3d_modeling_recipe]] §5).

---

## Consequences for ADR 20

1. **Confinement blocker → resolved (Gate 2 PASS).** R8 downgraded: emergent confinement is
   real to ~5% vs Mander; keep the "no dilation knob, validate Kc+backbone" caveat, drop the
   "must calibrate fc per member like Mander" framing. The §4 wording is corrected accordingly.
2. **Scope sharpened (Gate 1).** Conformal covers regular grid-aligned cages with no new code
   (*measured* to 12 bars), so the embedded element's justification narrows to **non-grid bar
   positions, ties/hoops, and the O(N²·nz) mesh cost** (*inferred* — still wants a hoop/non-grid
   demo). Build order: ship the conformal recipe + a bar-layout helper first; build
   `LadrunoEmbeddedRebar` for the cases conformal can't mesh.
3. **Solver gate confirmed.** Both gates needed KrylovNewton + step-halving; the adaptive
   wrapper is a real prerequisite, not a nicety.

## Reproduce
```
cd Ladruno_implementation/rc3d_gates
pythoncore-3.12-64/python.exe gate2_concrete_confinement.py   # PASS table above
pythoncore-3.12-64/python.exe gate1_conformal_column.py       # conformal column table
```
(Scripts auto-bootstrap `dist/bin`; set `LADRUNO_DIST` to override.)
