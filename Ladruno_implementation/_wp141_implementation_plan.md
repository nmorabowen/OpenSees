# WP-141 implementation plan — tasks, models, effort, oracles

Companion to [[141_rc_layered_shell_program]]. Drafted 2026-09-27; revised 2026-10-31 after the #882 review (H1–H3). Quality first, then tokens.

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
| P2f | Re-baseline on the #877 build (`-noeas`, lch = h, `LadrunoLoadControl`, `-local`): G-B both arms (O3b), G-C (O4b), G-PV PV19/PV20/PV27 (O4d) | Sonnet 5 high | O3b, O4b, O4d | 0.5 |
| P2d | **Baseline run + decision gate D-P2** (below) | Opus 5.5 medium | — | 0.25 |
| P3a | `Element::getCharacteristicLength(const Vector& n)` + `ASDShellQ4` override returning the full h(n) and reporting EAS (vanilla); after P3c/P3d and R2 | Opus 5.5 high + `/code-review high` | O1 | 0.75 |
| P3b | **Parked.** Directional latch at crack onset (`LadrunoRCConcrete`) — design+oracle, then C++ | Opus 5.5 **xhigh** → high; Fable 5.1 review | O2 | (1.5) |
| P3c | Crack-spacing band in reinforced layers (blocked on #877; Q5/Q6 first) | Opus 5.5 high | O3 | 0.75 |
| P3d | Tension-stiffening crack check (MCFT reserve over the bar sets; blocked on #877) | Opus 5.5 high | O4 | 0.75 |
| P3e | (optional) cyclic TS: ε1,max envelope + secant unload | Opus 5.5 high | O4c | 0.75 |
| P4a | Confined plate view — derivation + numpy oracle + fixture | Opus 5.5 **xhigh**; Fable 5.1 review | O5a–d, O5g | 1 |
| P4b | C++ kernel + wrapper + g++ check + Zone-A | Opus 5.5 high | O5b, O5c, O5e | 1 |
| P4c | Physics validation (shell tied-column vs brick+ties vs test) | Opus 5.5 high | O5f | 0.75 |

Total ≈ 10 sessions (P3e optional).

**Decision gate D-P2** (before any P3 code). If the baseline G-A error on realistic meshes is
below 10 % at every angle, or the error is dominated by mesh-bias *locking* (the band snaps to mesh
lines, which h(n) cannot fix), then re-rank P3a/P3b below P3c/P3d and record why in the ADR.

**D-P2 outcome (2026-09-28, from the partial G-A baseline): TRIGGERED.** The inclined-crack error is
dominated by effects h(n) cannot fix:
- the EAS lch factor, 2× energy (ADR G10);
- drilling-spring bridging (ADR G11);
- mesh locking: tanθ = ½ followed mesh lines, and 45° formed a 3-element-wide staircase.

Re-ranking:
- **New P2e (first):** establish when, if ever, an EAS `ASDShellQ4` localizes into one Gauss-point
  row, using strain gradients, asymmetric seeds and aspect ratios, with `LadrunoRCConcrete` and
  `ASDConcrete3D`. Quantify drilling contamination on PV20 and a diagonal-cracking wall. The outcome
  decides the correct length factor.
- **P3a re-scoped:** the new directional method carries the *corrected* length (EAS factor first,
  direction second), leaving the vanilla scalar untouched. (P2e then settled the factor: the full
  h(n), paired with `-noeas`; see the amendment below.)
- **P3b deferred** until a seeded-strip gate shows a single-element-wide band along θ exists.
- **G-A gate replaced** by a seeded straight strip or band at θ, `-noeas` or with a corrected factor,
  with drilling energy measured. The periodic cell stays as the locking detector.
- **P3c/P3d move ahead of P3b.** G-B and G-C already show the size of their errors (G-C: +117 % at ε_y).

**P2e outcome (2026-09-28, study try `work/claude/wp141-p2e-eas-drilling` @ 352dc4d, 71 runs).**

*EAS length factor.* EAS localizes into one Gauss-point column only through a bifurcation. That needs
the two node lines bounding the band to be free to stretch differently, as near a free edge (strips,
wide strips, notched beam). There h/2 is right (1.01–1.14). In constrained kinematics — the fork's 1-D
Bažant chain, the periodic cell, tied node lines, ν = 0 — both columns soften and h is right; h/2 gives
1.85–1.97.
- One constraint (`--tie-y`) moves the same mesh from 1.01 to 1.97.
- Strain gradients do not select a column; drilling plays no part.
- `LadrunoRCConcrete` and `ASDConcrete3D` agree within 0.01.
- ⇒ **Under EAS no fixed factor is right; `-noeas` + lch = h is deterministic (0.99–1.04).**

*Drilling.*
- PV20 4×4 with `-noeas`: the field stays uniform and irrotational, drilling energy ≡ 0, and all three
  drilling variants are identical (τ_max 4.276 vs test 4.26).
- Under EAS the drilling setting changes only the post-peak response (τ at γ = 0.02: 2.64 default /
  2.00 stab / 1.61 NL).
- On the 45° cell, `-drillingNL` releases the crack. It is a cumulative per-GP drilling damage, clipped
  at 0.999 and applied one step late (`ASDShellQ4.cpp:2046-2073`), which leaves θz nearly free — a
  conditioning risk.
- `-drillingStab 1.0` still stores 0.75·Gf·t·ℓ_c.
- G-A's "3.25 diffuse" did not reproduce.

*Re-ranking (proposed, owner to confirm):*
- **R1 (cheap, now):**
  - quirk rows + guide: crack-band RC on `ASDShellQ4` ⇒ `-noeas`, and `-drillingNL` where localized
    inclined cracks form;
  - `tests/test_ladrunoRCConcrete_meshobj.py` gates the ABSOLUTE energy under `-noeas`; its u_y = 0
    chain is exactly the 2× configuration;
  - apeGmsh warns when a crack-band layer section meets EAS.
- **R2 (gate before R1 becomes guidance):** the in-plane bending cost of `-noeas` — a cantilever wall
  with n = 1/2/4/8 elements across, EAS vs `-noeas` vs Timoshenko closed form. apeGmsh measured
  0.74× tip displacement for one element.
- ~~**P3a/P3b parked.**~~ **Amended 2026-10-31 (#882 review H1): P3b parked, P3a not.** The factor
  problem is solved by `-noeas` with no seam, and directional h(n) fixes neither mesh locking nor the
  3-element 45° staircase on the periodic cell, so P3b (the inclined-crack latch) waits. But P2e
  finding 3 measured a separate error that `-noeas` does not touch: on 2:1 elements the band is h_x
  wide while lch is the min node distance h_y, 2.06 under EAS and 1.58 under `-noeas` (partial; the
  length error alone predicts 2.0) — ADR G12. That is P3a's job: the seam returns the full h(n),
  paired with `-noeas`, and under EAS it reports EAS so the material warns or refuses (ADR Q2). P3a
  runs after P3c/P3d and R2; its consumer is the P3c onset rescale.
- **P3d, then P3c**, after #877 (both edit the paths #877 changes; their gates are re-baselined on
  its build in P2f). P3d first: it is the smaller change, and G-C/G-PV fully specify it.
- **New P3f (research):** a drilling law that follows the section's current secant stiffness instead
  of NL damage. This is vanilla `ASDShellQ4` — owner decision; an upstream suggestion is the
  alternative.
- **Upstream report** to the `ASDShellQ4` authors with the P2e evidence (evidence, not a patch).
  Owner approval before posting.

## 4. Oracle catalogue

**O1 — directional length seam (P3a)**
- O1a (T1) For a parallelogram/rectangle with edges a, b, the projection formula
  `h(n) = 2 / Σ_a |∇N_a(ξ=0)·n|` equals the centroidal chord `min(a/|n·x̂|, b/|n·ŷ|)`. On a square:
  a, 1.1547·a, √2·a at 0°, 30°, 45°. The seam returns the full h(n) with EAS on or off (no factor ½;
  P2e); under EAS it reports EAS.
- O1b (T2) An independent numpy `h(n)` on 1000 random convex quads; C++ within 1e-12. Probe it through
  a material-side response (`lchUsed`) so no second vanilla touch is needed.
- O1c (T3) Rigid rotation of the element in 3-D with `-local` rotated alike: h unchanged (1e-12).
  Every other element keeps the scalar default, and the **whole Zone-A battery is bit-identical**.
- O1d (T3) Under `-corotational`: in-plane tension along local x′ of a rotated element gives h = a,
  not a global projection. This checks the frame of n (ADR open question 1).
- O1e (T3, **the physics claim**, replaces O2c as P3a's gate) Seeded strip, `-noeas`, plain
  `LadrunoRCConcrete -autoRegularization` consuming the seam, elements 2:1 and 1:2, crack normal along
  each side: W/(Gf·t·ℓ_c) = 1.00 ± 0.05. Baseline with the scalar (min node distance) on 2:1 with the
  normal along the long side: ≈ 2.0 (P2e beam rows 1.58, partial); square elements: 1.00 (control).
- O1f (T3) Under EAS, a material with `-autoRegularization` behaves as ADR Q2 decides (one warning,
  or the sentinel); `-noeas` is silent.
- Mutants: apply the upstream factor ½ under EAS (killed by O1a/O1f); use `n` in the global frame;
  return the projected width a(|c|+|s|) instead of the chord; use min(h_x, h_y) (killed by O1e).

**O2 — directional latch (P3b, parked; O2c is now the recorded locking detector, not a gate)**
- O2a (T2) `rc_shell_ref.py`: the latch fires in the step where ε1 first reaches ε_cr and uses that
  step's principal direction. The rescale happens before any softening (elastic pre-crack response
  is lch-free). Fixture → `rc_reg_gpp.cpp` byte-check 1e-14.
- O2b (T1, plumbing only) One element, homogeneous uniaxial tension at θ: dissipated energy per
  volume × h(n) = Gf (1e-6) for θ ∈ {0, 14.04, 26.57, 45}°.
- O2c (T3, **the physics claim**) G-A periodic cell:
  - square L×L, structured n×n `ASDShellQ4` (EAS on and `-noeas`), plain `LadrunoRCConcrete -autoRegularization`, membrane only;
  - periodic `equationConstraint` pairs to a macro node; macro uniaxial stress along θ with tan θ ∈ {0, ¼, ½, 1} (angles compatible with square periodicity);
  - one 5 %-thinner element seeds the crack.
  - **Oracle:** `W_diss / (t · ℓ_c) = Gf`, with ℓ_c = L·√(p²+q²) for tanθ = p/q in lowest terms.
    The periodic crack is a closed line of direction (−p, q) on the torus: L at 0°, L√2 at 45°, L√5
    at tanθ = ½, L√17 at tanθ = ¼. It equals L/cosθ only at 0° and 45° (corrected 2026-09-27).
  - Recorded, not gated (P3b parked; amended 2026-10-31). The former "±10 % at all angles after P3" was unattainable: 45° `-noeas` reads 2.13 from locking and drilling, which h(n) does not fix.
  - The damaged-element map must follow θ; if the band locks to mesh lines, the test reports locking instead of passing.
  - Smoke first: θ = 0 must reproduce the existing Bažant-bar result (`test_ladrunoRCConcrete_meshobj.py`).
- O2d (T3) Wrong-instant guard: a non-proportional pre-crack path rotates the principal axis by 30°
  between the first call and crack onset. Latching at the first call gives a different h, and the
  test must see it.
- O2e Flags-off byte identity on the full RC suite (90/90 after #877).
- Mutants: latch at first call; latch every step; skip the rescale under `-implex`.

**O3 — crack-spacing band (P3c)**
- O3a (T1, RTD 1016 §2.4.3.1) `h_eff = min(h, s_θ)`, `s_θ = 1/(|cosθ|/s_rx + |sinθ|/s_ry)`, with θ
  between the crack normal at onset and x; h is the scalar today, h(n) after P3a. s_rx, s_ry are
  **mean** crack spacings, user-given or computed by the apeGmsh builder as s_r,max(EC2 7.11)/1.7
  (ADR Q5 — owner to confirm). Limits θ = 0 → s_rx and θ = 90° → s_ry are unit tests. The rule
  applies only to layers flagged reinforced, to tension only.
- O3b (T3) G-B tie, re-run on the #877 build (P2f) with `-noeas`, lch = h, `LadrunoLoadControl`:
  - a membrane strip with `PlateRebar` (ρ ≈ 1 %), **elastic** steel, length L = 12·s_rm, mesh
    h/s_rm ∈ {½, 1, 2, 4}, one element 1 % thinner as the seed;
  - **TS-off arm (the P3c gate), energy oracle (T1):** with perfect bond every element cracks once,
    the tie forms L/s_rm cracks. Run until every concrete GP has d_t > 0.99; `W_conc` = Σ_gp ∫σ:dε
    over the concrete layers minus their remaining elastic energy. Gate
    `W_conc / (G_F·A_c·L/s_rm) = 1.00 ± 0.05` at h/s_rm ∈ {1, 2, 4}. Baseline expected
    1.00 / 0.50 / 0.25 (= s_rm/h: each element dissipates G_F·A_c once); h/s_rm = ½ recorded only,
    expected 2.0 before and after P3c (ADR Q6);
  - **TS-on arm (no-regress for P3c and P3d):** the EC2 eq. 7.9 closed form for mean strain,
    `ε_m = max(σ_s/E_s − k_t·f_ct·(1+α_e·ρ)/(ρ·E_s), 0.6·σ_s/E_s)`, k_t = 0.6, σ_s = N/A_s, within
    ±15 % in stabilized cracking (σ_s ∈ [1.3·σ_sr, 0.8·f_y]) and mesh spread ≤ 10 pp. The
    2026-09-28 run passed it at baseline (+7.0–7.1 %, 0.1 pp);
  - EC2 7.9 is **not** an oracle for the TS-off arm: it contains tension stiffening, so TS off sits
    +55–62 % above it by construction (the 2026-09-28 numbers), whatever P3c does.
- O3c (T3, optional physics) a tension-tie test with counted cracks, to check that the chosen s_rm
  is physical — O3b only checks the rule against the s_rm it is given.
- Mutants: `max` instead of `min` (O3b TS-off fails at h/s_rm = 2, 4); spacing applied to plain layers
  (a plain-layer arm keeps 1/h dissipation); latch kept after a cut step (a forced cut-and-retry at
  crack onset must reproduce the uncut result).

**O4 — tension-stiffening crack check (P3d) and cyclic TS (P3e)**
P3d is the MCFT crack check (ADR §4), not a cutoff on ε1:
`σ_ts ≤ σ_cap = Σ_i ρ_i·max(0, f_y,i − f_s,i)·cos²θ_i`, `f_s,i = clip(E_s,i·n_iᵀεn_i, ±f_y,i)`.
- O4a (T1, material point) σ_cap in closed form for one bar set at θ ∈ {0°, 30°, 60°} to the crack
  normal and for two orthogonal sets with one yielded; 1e-14. Where σ_ts < σ_cap the stress and
  tangent equal the uncapped path bit for bit. σ_cap → 0 continuously as the last bar yields; g++
  one-sided FD tangent on both sides of the activation point.
- O4b (T3) G-C: a lightly reinforced (ρ 0.3–0.5 %) panel in uniaxial tension along the bars with
  elastic-perfectly-plastic steel.
  - Gate with the crack check: `N/(A_s·f_y) ≤ 1.01` at every strain after cracking, and the plateau
    `1.00 ± 0.01` on [ε_y, 10·ε_y]. In a tie along the bars the check is exactly "steel at the
    crack ≤ f_y", so the cap may bind before ε_y — that is the physics, not a defect.
  - Baseline without the rule: `A_s·f_y + A_c·σ_ts(ε_y)`, +117 / +70 % at ε_y pre-#877 (c = 500);
    re-measured in P2f (c = 200).
  - Uniaxial only: it cannot tell a crack check from a cutoff on ε1, because here ε1 is the bar strain.
- O4d (T3) G-PV: PV19/PV20/PV27 (study `pv_panels`, MCFT flag set, `-local` pinned, EAS and `-noeas`).
  - Gate: τ_max within ±7 % of test, and the pre-peak τ–γ curve bit-identical to the pre-P3d binary.
  - Expected: the cap is inactive up to the peak. From the study CSVs, σ_cap/σ_ts at τ_max is
    1.93 / 2.13 / 3.67 (PV20 / PV19 / PV27, ν = 0.2) and 1.23 / 1.35 / 2.95 (ν = 0). Re-measured in P2f.
  - Why it discriminates: ε1 passes the y-bar ε_y = 1.5e-3 at γ ≈ 1.7–2.0e-3 and the y bars yield at
    γ = 3.4e-3 / 4.75e-3, all before γ_peak = 8.1e-3 / 10.6e-3. A cutoff on ε1 or on "any bar
    yielded" removes tension stiffening there and gives the `-betaC`-alone curve, 3.61 MPa
    (−15 %), which fails.
- Mutants: cutoff on ε1 (killed by O4d); cutoff when any bar yields (O4d); cos² dropped (O4a at 30°,
  60°); cap left out (O4b).
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
#877 merged:    P1a, P1d, P2f re-baseline  ->  P3d -> P3c (-> P3e)
R2 + P3c done:  P3a (seam, -noeas, elongated elements)
parked:         P3b (until a seeded strip shows a one-element band along theta)
P1a landed:     P4a -> P4b -> P4c              [B1/B2 arrive with #877]
```

Branches: each code phase is its own `wp/<n>-<slug>`, numbered when opened, cut from fresh
`ladruno`, with a draft PR on day one. #882 stays the design record (ADR + this plan).
