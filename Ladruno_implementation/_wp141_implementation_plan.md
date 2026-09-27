# WP-141 implementation plan — tasks, models, effort, oracles

Companion to [[141_rc_layered_shell_program]]. Drafted 2026-09-27. Quality first, then tokens.

## 0. Rules for every phase

- **Kernel-oracle doctrine** (`LadrunoConcrete3D_guide.md` §10): numpy oracle first → fixture → header
  kernel → g++ byte-check → Zone-A pytest → one Windows build → structural rerun. No OpenSees build
  until the g++ check is green (a build costs ~25–30 min; the oracle loop costs seconds).
- **Oracle tiers.** Every gate names its tier:
  - **T1** a closed form or a published rule, derived by hand;
  - **T2** an independent numpy oracle + g++ byte-check of the C++ kernel against its fixture;
  - **T3** in-OpenSees: flags-off bit identity, same-binary limit cases, FD tangents, structural
    energy/force checks, experiments.
  A T3 test that passes on a T2 fixture it generated itself proves nothing; say which tier carries
  the claim.
- **Every new flag defaults off and is byte-identical to the previous binary** on the full RC and
  Concrete3D suites.
- **Mutation gate** (`ci/mutation_gate.py`, floor 0.60) on every code phase. At least one named
  mutant per oracle below must die. A surviving mutant means the test observes the wrong path
  (WP-124 lesson).
- **Warrant package before ready** (ADR-87): verification manifest row, mutation score, guide update,
  ledgers (vanilla / implementations / quirks), Zone-A green. The owner merges.
- **Token rules:**
  - one fresh session per phase, handoff through this file and the ledgers, not transcripts;
  - subagents return conclusions of 700 words or fewer and read files in ranges;
  - pytest one module per call, with long runs in the background (the watchdog stalls seen on
    2026-09-27 were silent long runs);
  - batch all C++ edits of a phase into one build;
  - reuse the validation study's harness for structural runs.

## 1. Model and effort policy

| Work type | Model | Effort | Why |
|---|---|---|---|
| New math: derivation, oracle authoring, tangent derivation (P3b, P4) | Opus 5.5 | **xhigh** | errors are silent and survive tests that reuse the same derivation |
| C++ port against a frozen oracle; vanilla seam; structural rigs | Opus 5.5 | high | exacting but checkable against T1/T2 |
| Adversarial formulation review (P3b, P4 only) | **Fable 5.1** | high | a different lens; caught four wrong-reason passes in ADR-66 G8 |
| Code review of vanilla edits (P1b, P3a) | Opus 5.5 via `/code-review high` | high | core/vanilla surface |
| Mechanical code (response forwarding, parser plumbing), rig variants | Sonnet 5 | medium | an oracle catches mistakes |
| Docs, ledgers, banner, manifest rows | Sonnet 5 | low | template work, then review in the PR |
| Lookups, grep audits | Haiku 4.5 (Explore) | low | read-only |

Effort is a **session** setting; subagents inherit it and only the model can be chosen per agent. So
each phase runs as its own session at the effort of its hardest task, and cheaper work inside it goes
to Sonnet/Haiku subagents. No task uses `max`.

## 2. Dependencies (external)

| Needs | Blocks | Owner |
|---|---|---|
| PR #877 merged (RC schema v6, `-crackedNu/-betaC`, TS `c=200`; CDPM2 B1/B2 via #867; element-side return codes C3). Zone-A was red on 3 tests on 2026-09-27; a test-only fix sits unpushed on the study's local branch `wp/concrete3d-integration-zonea` (894d1b138) — **owner action** | P1a, P1d, P2 baseline runs, P3, P4 | owner |
| apeGmsh #1185 merged (`RCLayeredShell`, `PlateFiber`, guards, `no_eas`) | P2 rig authoring convenience (not correctness) | owner |

**Agreed with the concrete study session (2026-09-27):**
- WP-141 owns the PlateFiber part of study B4 (P1a). The study will not touch `LadrunoConcrete3D.cpp` until P1a lands.
- After #877 the study plans no edits to `LadrunoRCKernel.h` / `LadrunoRCConcrete.cpp`. It may still edit `LadrunoConcrete3DKernel.h` (CDPM2 compression calibration), under the oracle doctrine.
- No one else plans to touch `Element.h/.cpp`, `ASDShellQ4`, `PlateRebarMaterial`, or the hoop/BeamFiber code.
- The structural gates (P2) live in the study repo as `benchmarks/03_layered_shell/<case>/`, following its conventions:
  - `model.py` with an argparse CLI, a README, and `plot.py` writing `results.md`;
  - commits go to a try branch `work/<gupi>/<topic>`, never `main`, with commits ending in a `Gupi: <name>` line + Co-Authored-By;
  - prescribed displacements run through `ops.integrator.LadrunoLoadControl(dlam, tangent_predictor=True)` (the owner's standing rule).
- The kernel-level oracles stay as fork tests.
- Baseline numbers already measured on the #877 build, PV20 τ_max:

  | Model | τ_max (MPa) |
  |---|---|
  | `LadrunoRCConcrete` + MCFT + `-crackedNu` | 4.53–4.62 (= MCFT hand solution) |
  | `LadrunoConcrete3D` PlateFiber | 4.67, still rising |
  | `ASDConcrete3D` | 4.41 |
  | test | 4.26 |

## 3. Tasks

| ID | Task | Model / effort | Oracles | Est. sessions |
|---|---|---|---|---|
| P1a | `LadrunoConcrete3D` PlateFiber σ33 condensation: damped Newton, tolerance relative to `ft`, honest nonzero return (**study B4, handed to WP-141**); first-ever PlateFiber tests | Opus 5.5 high | O7 | 0.75 |
| P1b | `PlateRebarMaterial` forwards unknown `setResponse` keys to its uniaxial (vanilla, additive, `// Ladruno`) | Sonnet 5 medium + `/code-review high` | O6 | 0.5 |
| P1c | Verify the `ASDShellQ4` default in-plane frame in fork vs upstream 3.7.1 sources; RC guide + quirk rows (`-local` with rebar; `ASDSteel1D` in shells) | Opus 5.5 medium | O8 | 0.25 |
| P1d | Doc drift G8 (lch "= edge", midpoint 1/n², in-kernel steel claim, ADR-19 banner) — after #877 | Sonnet 5 low | grep-clean | 0.25 |
| P2a | G-A inclined-crack rig: periodic unit cell (study repo) | Opus 5.5 high | O2c | 1 |
| P2b | G-B reinforced tension tie rig | Sonnet 5 high | O3b | 0.5 |
| P2c | G-C lightly reinforced panel rig | Sonnet 5 high | O4b | 0.25 |
| P2d | **Baseline run + decision gate D-P2** (below) | Opus 5.5 medium | — | 0.25 |
| P3a | `Element::getCharacteristicLength(const Vector& n)` + `ASDShellQ4` override (vanilla) | Opus 5.5 high + `/code-review high` | O1 | 0.75 |
| P3b | Directional latch at crack onset (`LadrunoRCConcrete`) — design+oracle, then C++ | Opus 5.5 **xhigh** → high; Fable 5.1 review | O2 | 1.5 |
| P3c | Crack-spacing band in reinforced layers | Opus 5.5 high | O3 | 0.75 |
| P3d | Tension-stiffening cutoff at steel yield | Opus 5.5 high | O4 | 0.5 |
| P3e | (optional) cyclic TS: ε1,max envelope + secant unload | Opus 5.5 high | O4c | 0.75 |
| P4a | Confined plate view — derivation + numpy oracle + fixture | Opus 5.5 **xhigh**; Fable 5.1 review | O5a–d, O5g | 1 |
| P4b | C++ kernel + wrapper + g++ check + Zone-A | Opus 5.5 high | O5b, O5c, O5e | 1 |
| P4c | Physics validation (shell tied-column vs brick+ties vs test) | Opus 5.5 high | O5f | 0.75 |

Total ≈ 10 sessions (P3e optional).

**Decision gate D-P2** (before any P3 code). If the baseline G-A error on realistic meshes is
below 10 % at every angle, or the error is dominated by mesh-bias *locking* (the band snaps to mesh
lines, which h(n) cannot fix), then re-rank P3a/P3b below P3c/P3d and record why in the ADR.

## 4. Oracle catalogue

**O1 — directional length seam (P3a)**
- O1a (T1) For a parallelogram/rectangle with edges a, b, the projection formula
  `h(n) = 2 / Σ_a |∇N_a(ξ=0)·n|` equals the centroidal chord `min(a/|n·x̂|, b/|n·ŷ|)`. On a square:
  a, 1.1547·a, √2·a at 0°, 30°, 45°. Under EAS the upstream factor ½ applies.
- O1b (T2) An independent numpy `h(n)` on 1000 random convex quads; C++ within 1e-12. Probe it through
  a material-side response (`lchUsed`) so no second vanilla touch is needed.
- O1c (T3) Rigid rotation of the element in 3-D with `-local` rotated alike: h unchanged (1e-12).
  Every other element keeps the scalar default, and the **whole Zone-A battery is bit-identical**.
- O1d (T3) Under `-corotational`: in-plane tension along local x′ of a rotated element gives h = a,
  not a global projection. This checks the frame of n (ADR open question 1).
- Mutants: drop the EAS factor; use `n` in the global frame; return the projected width
  a(|c|+|s|) instead of the chord.

**O2 — directional latch (P3b)**
- O2a (T2) `rc_shell_ref.py`: the latch fires in the step where ε1 first reaches ε_cr and uses that
  step's principal direction. The rescale happens before any softening (elastic pre-crack response
  is lch-free). Fixture → `rc_reg_gpp.cpp` byte-check 1e-14.
- O2b (T1, plumbing only) One element, homogeneous uniaxial tension at θ: dissipated energy per
  volume × h(n) = Gf (1e-6) for θ ∈ {0, 14.04, 26.57, 45}°.
- O2c (T3, **the physics claim**) G-A periodic cell:
  - square L×L, structured n×n `ASDShellQ4` (EAS on and `-noeas`), plain `LadrunoRCConcrete -autoRegularization`, membrane only;
  - periodic `equationConstraint` pairs to a macro node; macro uniaxial stress along θ with tan θ ∈ {0, ¼, ½, 1} (angles compatible with square periodicity);
  - one 5 %-thinner element seeds the crack.
  - **Oracle:** `W_diss / (t · L / cos θ) = Gf`. Gate after P3: ±10 % at all angles. The baseline is recorded, not gated.
  - The damaged-element map must follow θ; if the band locks to mesh lines, the test reports locking instead of passing.
  - Smoke first: θ = 0 must reproduce the existing Bažant-bar result (`test_ladrunoRCConcrete_meshobj.py`).
- O2d (T3) Wrong-instant guard: a non-proportional pre-crack path rotates the principal axis by 30°
  between the first call and crack onset. Latching at the first call gives a different h, and the
  test must see it.
- O2e Flags-off byte identity on the full RC suite (90/90 after #877).
- Mutants: latch at first call; latch every step; skip the rescale under `-implex`.

**O3 — crack-spacing band (P3c)**
- O3a (T1, RTD 1016 §2.4.3.1) `h_eff = min(h(n), s_θ)`, `s_θ = 1/(|cosθ|/s_rx + |sinθ|/s_ry)`, with θ
  between the crack normal and x. s_rx, s_ry are mean crack spacings, user-given or computed by the
  apeGmsh builder from EC2 7.11. Limits θ = 0 → s_rx and θ = 90° → s_ry are unit tests. The rule
  applies only to layers flagged reinforced.
- O3b (T3) G-B tie:
  - a membrane strip with `PlateRebar` (ρ ≈ 1 %), mesh h/s_rm ∈ {1, 2, 4} — the coarse-shell regime this rule targets; h < s_rm is plain crack band, already covered by `_meshobj`;
  - **Oracle:** the EC2 eq. 7.9 closed form for mean strain, `ε_m = max(σ_s/E_s − k_t·f_ct·(1+α_e·ρ)/(ρ·E_s), 0.6·σ_s/E_s)`, with k_t = 0.6 and σ_s = N/A_s;
  - gate: within ±15 % in stabilized cracking, and mesh spread ≤ 10 %;
  - run with TS off (the cyclic regime) and with TS on + cutoff.
- Mutants: `max` instead of `min`; spacing applied to plain layers.

**O4 — tension-stiffening cutoff (P3d/P3e)**
- O4a (T1) σ_ts = 0 for ε1 ≥ ε_y, continuous, with no jump: a fade over a band chosen in the P3 design
  session. Below the fade start the old curve holds to 1e-15. g++ one-sided FD tangent at both
  fade boundaries.
- O4b (T3) G-C: a lightly reinforced (ρ 0.3–0.5 %) panel in uniaxial tension with elastic-perfectly-plastic steel.
  - Gate with the cutoff: the plateau `N/(A_s·f_y) = 1.00 ± 0.01`.
  - Baseline without the cutoff: `A_s·f_y + A_c·σ_ts(ε_y)`, recorded as the documented overestimate.
- O4c (T1, P3e) Unloading from (ε1,max, σ_ts) runs along the secant to the origin and reloads along
  it. No re-inflation on unload (monotone check).

**O5 — confined plate view (P4)**
- O5a (T1) Elastic closed form. With in-plane strains prescribed and the tie active (ε33 > 0):
  `ε33 = −λ(ε11+ε22)/(λ+2μ+K)`, then σ11, σ22, τ12 in closed form; 1e-12. The tie is tension-only
  (slack for ε33 ≤ ε_p).
- O5b (T3, limits, nonlinear, same binary) `hoopK = 0` is byte-identical to the post-B4 plate view.
  `hoopK = 1e6·E` matches the 3-D view driven with ε33 = 0, retained components within 1e-6,
  **damage included**. This is the strongest single check of the whole chain.
- O5c (T2) `concrete3d_ref.py::confined_plate_step`, mirroring `confined_step` (:1501) with one
  confined direction, **nominal** balance and tie plastic memory. Fixture → g++ byte-check →
  condensed-tangent FD 1e-6.
- O5d (T1) Tie memory: load past yield, then unload. Slope K, residual strain `ε_max − f_y/K`, zero
  force while ε33 ≤ ε_p.
- O5e (T3, pre-peak cross-view) The plate view with ε22 set to the BeamFiber view's converged
  lateral strain reproduces its σ22 = σ33 = −σ_hoop and σ11 (1e-8). This holds only where ω ≈ 0,
  because the BeamFiber view balances the effective stress.
- O5f (T3, physics) A tied column as a shell: plane = axis × one lateral, thickness = the other
  lateral.
  - Through-thickness tie legs → `-hoop K fy`; in-plane legs → `PlateRebar` at 90°.
  - Compare against study E3 (Sheikh & Uzumeri 1980 2A1-1 / 4B3-19 / 4D6-24, and the `LadrunoBrick` + embedded-tie model).
  - Gates: peak within 10 % of test; post-peak ductility ordered by tie spacing; tie stress rising at the peak.
  - Also PV20 with `LadrunoConcrete3D` layers and `-tcTemper proj`, `hoopK = 0`: must not regress
    (the study saw the strut collapse without that flag).
- O5g Mutants: effective instead of nominal balance (killed by a damaged state with ω > 0.3); tie
  memory dropped; tie stiffness left out of `condenseTangent()`.

**O6 — PlateRebar forwarding (P1b)** (T3 identity) A standalone uniaxial driven with
`ε11c² + ε22s² + γ12cs` from the same history gives the same forwarded response (bit-exact). The
plate `stress` stays 5 components. Mutant: forward to a stale copy.

**O7 — plate-view honest failure (P1a)**
- (T1) Every converged point has |σ33| ≤ tol·ft (a self-certifying residual).
- (T2) A numpy plate-fiber driver in `concrete3d_ref.py` (nested ε33 Newton over
  `damaged_step_tensor`). Fixture block → g++ byte-check → condensed-tangent FD. The PlateFiber view
  has had no test at all until now.
- (T3, same binary) Drive the 3-D view with the plate's converged ε33: σ33 = 0 and the retained
  components are identical.
- (T3) A forced snap-back state returns −1, and `ASDShellQ4` cuts the step. PV20 with
  `LadrunoConcrete3D` layers and `-tcTemper proj` shows 0 silent failures (`returnFailures` response
  from #877).
- Mutants: return 0 on non-convergence; absolute tolerance; undamped step.

**O8 — rebar frame (P1c)** (T3) One `ASDShellQ4` with `PlateRebar` at 0° and `-local` along a
rotated direction: membrane N/ε along `-local` equals `E_c(h−Σt) + E_s·t` (1e-6), and without
`-local` the frame follows the fork's documented default. The apeGmsh live test already does this;
mirror it in fork Zone-A.

## 5. Sequencing

```
now:            P1b, P1c, P2a–c (rigs, study repo try branch)   [independent of #877]
#877 merged:    P1a, P1d, P2d baseline + D-P2  ->  P3a -> P3b -> P3c -> P3d (-> P3e)
P1a landed:     P4a -> P4b -> P4c              [B1/B2 arrive with #877]
```

Branches: each code phase is its own `wp/<n>-<slug>`, numbered when opened, cut from fresh
`ladruno`, with a draft PR on day one. #882 stays the design record (ADR + this plan).
