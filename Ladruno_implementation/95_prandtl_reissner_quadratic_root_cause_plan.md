# ADR-95 — Prandtl–Reissner, exact campaign: root cause of the quadratic (Bezier) wall

**Status:** PLAN, pre-registered 2026-09-06. Branch `wp/95-prandtl-bezier-root-cause`.
**Owner decision points:** D1 (go on P0), D2 (fix lane after P1), D3 (P5 wall-clock spend).
**Predecessors:** notes 81 (#721), 82 (#725), 83 (#727), gate #722, TIMs T0–T4, ADR-90 report.

## 0. The question, restated so it can be answered

"Why do the Bezier elements fail the Prandtl–Reissner test" is three stacked deficits,
and only the first is unexplained:

| deficit | size | status |
|---|---|---|
| A. The **quadratic-class wall**: every quadratic leg (H20 Lagrange std/uri, BezierTet10 std/bbar, TenNodeTet) dies on the step floor while still hardening, at s/B ≈ 0.001–0.011 | linear bbar hex 0.998 vs best quadratic 0.70–0.77 | **OPEN — this campaign** |
| B. The **tet penalty** on top of A (Bezier -bbar 0.49–0.58 vs H20 uri 0.70–0.77, both walled) | ~1.4× | measured only at unequal allowances; quantify at matched settlement (P3) |
| C. **Controller allowance** (budget 24→48 moved Bezier +17%) | ±5–15% | understood (note 82 §7.1); rule: no walled number is a capacity |

Bezier itself is exonerated as a *basis* (the pure-Lagrange H20 walls identically). The
campaign therefore identifies the event on the cheap in-tree H20 deck, then **transfers**
the identity to BezierTet10 (P3). A result that names the event in H20 but does not
reproduce in BezierTet10 is a failure of the campaign, not a success.

## 1. What is already established (do not re-derive)

Eliminated, with the note that killed each: constraint ratio r (82 §6), isochoric span
(82 §6.1), volumetric constraint twice (82 §6 + TIMs flow-rule pair), Bernstein basis
(81), tet geometry as the *cause* (81, ~2× only), the driver (linear plateaus in the same
harness), the limit point (82 §7.3: σ_min flat then 1000× drop across 0.06 mm),
checkerboard (82 §6.3, weak, reverses with h), nondeterminism (marginal-decision
mechanism, TIMs n=3), a 27-node lane (82 §9), spurious uri modes (81).

Established positives: the wall is **in the Newton path** — DR walks through it to
settled equilibria (83 §5.2); the tangent event is **abrupt** (σ_min 2.16e-4 → 2.29e-7
between s/B 0.01115 and 0.01118 on `h20uri h1.0`); DR's own ramp residual grows 10×
past the same settlement (83 §5.3) ⇒ something in the **material/discretisation state**
changes there, not only the operator.

## 2. Pre-registered hypotheses and their predictions

The deck: UW `DruckerPrager`, φ_txc = 20°, ν = 0.45, ρ̄ = 0, **SY = 0.2 kPa** (apex
regulariser). UW DP is a **two-surface** model: f1 = cone, f2 = tension cutoff at
I1 = T, with `mTo = √(2/3)·σ_y/ρ` — at SY = 0.2 the cutoff sits at I1 ≈ 0, so any
Gauss point whose mean stress reaches zero (heave zone beside the footing edge) enters
the f2 or the f1∧f2 **corner** branch, whose consistent tangent is a different operator.
The `count > 3` forced-accept bailout prints `Jact =` and appears in **no** campaign
log, so that path is excluded as the trigger.

| # | hypothesis | P1 prediction (must be seen at the event, absent before) | P2 knob that must MOVE the wall |
|---|---|---|---|
| **H1** corner/tension-cutoff branch switch | step change in #GPs on f2 or f1∧f2 at s/B ≈ 0.01118; linear control shows no such step (or shows it and survives — which kills H1) | SY 0.2 → 2 → 20 kPa; disabling the cutoff (`mTo = 1e10`) pushes the wall OUT |
| **H2** loss of ellipticity (Rudnicki–Rice: non-associated + ν = 0.45) | min over GPs of det(acoustic tensor n·D_ep·n) crosses 0 at the event; branch counts smooth | ν 0.45 → 0.30 delays the wall; ρ̄ = ρ/2 delays it (TIMs ρ̄ = ρ walled EARLIER — H2 already constrained, so this is the falsifier) |
| **H3** Newton algebra only (no GP event) | branch counts and det(A) both smooth; per-iteration residual localises on the footing-edge nodes and diverges from iteration 1 | line search / `NormUnbalance` / `-initial` tangent do not move it (82 §7.2 says they do not) → H3 dies and H1/H2 must be re-examined |

A prediction not met is reported as such. Numbers from walled legs are **allowances**;
only branch/ellipticity **states** at matched settlement are compared across elements.

## 3. Phases, agents, gates

Orchestrator = this session. Agents write `Ladruno_implementation/_adr95_<phase>_results.md`
with a ≤20-line VERDICT header; the orchestrator reads headers only. One build per phase.

### P0 — Instrumentation (Opus, ~1 build) → D1
1. `DruckerPrager::setResponse` (UWmaterials, **vanilla file** → ledger row + `// Ladruno`
   comments): new response `ladrunoBranch` = `[branch(0 el,1 f1,2 f2,3 corner), gamma0,
   gamma1, f1_trial, f2_trial, forcedAccept, I1, detAmin]`, where `detAmin` is the
   minimum over ~200 unit directions of det(n·D_ep·n) using the *consistent* tangent
   from the last update. Computed on request only (sampling cost, not per-step cost).
2. Unit test: single-element DP driven onto each branch, asserting the flag and that
   `detAmin > 0` in the elastic range. Zone-A pytest.
3. Python-side Newton forensics in `quad_path_diag.py`: per-iteration `||R||` and the
   ten largest `nodeUnbalance` DOFs at the first failed step (`--forensics`).
Build `OpenSees OpenSeesPy` in this worktree (WMI launch + `until -f` watcher), verify
`opensees.pyd` mtime, `ladrunoBuild()` hash in every probe.

### P1 — The trajectory (Sonnet, ~2 h wall) → decides H1/H2/H3
`quad_path_diag.py --elem h20uri --h0 1.0 --branch` sampled at note 82 §7.3's stations
plus every 5e-6 of s/B in [0.0108, 0.0113]; **linear control** `h8bbar --h0 1.0
--branch` at the same stations (mandatory: a branch step that the linear element also
takes and survives is not the cause). Output: branch histogram, detAmin, σ_min vs s/B
on one CSV. Decision per §2 table.

### P2 — Confirmatory knob (Sonnet, ~4 leg-hours)
Only the knob of the surviving hypothesis, three levels, `h20uri` **and** `h8bbar` each
(the linear leg must keep plateauing at 0.99–1.00 or the knob changed the physics).
Wall position (s/B at floor) must move monotonically in the predicted direction; report
as reach, never as capacity.

### P3 — Bezier transfer (Sonnet, `opensees_env`, ~3 h)
`build_mesh_tet10.py` (apeGmsh, Bernstein-consistent loads per ADR 0091) → same
`--branch` sampler on `BezierTet10 -bbar` and `TenNodeTetrahedron` at matched DOF.
Deliverables: (a) same event identity at the tet wall — yes/no; (b) deficit B measured
honestly: q at matched s/B = 0.008 for h8bbar / h20uri / bezier-bbar / tet10 (all still
on-path there), which separates the tet penalty from the wall for the first time.

### P4 — Fix lane (Opus) → D2, conditional
- H1 → apex/corner smoothing on the fork side: hyperbolic DP (Abbo–Sloan) or a
  `-noTensionCutoff` option; gate = wall gone on `h20uri` AND BezierTet10 at unchanged
  linear answer; the standing gate #722 must stay bit-identical.
- H2 → physics, not a bug: the quadratic elements resolve a bifurcation the linear
  bbar hex cannot represent. Fix = path regularisation (ADR-90 wrapper already exists as
  a tangent regulariser even though it does not give a converged band width); report
  the linear 0.998 as the *smeared* answer, not the exact one.
- H3 → Newton variant lane; decide after forensics.

### P5 — Close note 83 §7 item 1 (background, ~6 leg-hours) → D3
Continue the DR `dtfac = 0.25` quadratic leg to s/B = 0.05 with a zero-increment hold
every 0.00125. Launched only after P1, so the DR sampler also records the branch state.
Answers "does a quadratic element plateau at all" independently of the fix.

### Review
Red/blue adversarial review of `_adr95_results` before the PR flips to ready (the
three prior retractions on this problem are the reason).

## 4. Rules in force
- Termination mode on every leg (TARGET / BUDGET / FLOOR / WALL / CLOCK); capacity only
  if TARGET with flat tail and free advance; allowance legs quoted with the allowance.
- Controls run FIRST and are quoted in the same table.
- Bit-identical repeat of any leg whose number enters a conclusion (n = 2 minimum).
- `ladrunoBuild()` in every probe; `import opensees` after `sys.path.insert(dist/bin)`.
- No Prandtl number from a walled leg is compared across elements; only states are.
- Build in this worktree only; never on the shared checkout.

## 5. Deliverables and ledgers
- `_adr95_p0..p5_results.md`; final `95_prandtl_reissner_quadratic_root_cause.md` note.
- `LEDGER_vanilla_files.md`: `DruckerPrager.cpp/.h` (branch response). `LEDGER_quirks.md`:
  UW DP two-surface corner at small SY; `Jact =` print = forced accept. `LEDGER_implementations.md`
  only if P4 ships an option (+ banner line).
- Cost estimate: P0 1 build + 1 agent; P1–P3 ~9 leg-hours across two Sonnet agents;
  P4/P5 on owner decision.
