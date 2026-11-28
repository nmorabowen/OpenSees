# Dual-number program — dispatch briefs

Copy one brief into the Agent prompt, then fill `<PR>` and `<NEXT>` (the last handoff's Next, or
"start"). Every brief begins with the shared header below.

---

## Shared header (prepend to every brief)

You are a worker in the Ladruno dual-number program (forward-mode automatic differentiation in the
OpenSees fork `github.com/nmorabowen/OpenSees`, default branch `ladruno`). Read, in order:
`AGENTS.md`, `.claude/skills/ladruno-dual/SKILL.md` (the worker rules there are binding),
`Ladruno_implementation/_wp262_dual_orchestration_plan.md` (your row in §2), and the WP-262
recipe: `SRC/material/nD/LadrunoDual.h`, `SRC/material/nD/LadrunoJ2Kernel.h` (`returnMapT`,
`plasticResidual`, the IFT step), `SRC/material/uniaxial/LadrunoUniaxialJ2Kernel.h` +
`LadrunoUniaxialJ2.cpp` (kernel extraction + DDM by `Dual<1>`), `tests/ladrunoj2_dual_check.cpp`
and `tests/ladruno_uniaxialj2_dual_check.cpp` (oracle + bit-identity + break-the-gate patterns),
and the two WP-262 quirk fragments in `Ladruno_implementation/ledger/quirks/`.

Your worktree was created from WP-262's head. First act: `git switch -c wp/<n>-<slug>` there.
Open your draft PR yourself: `gh pr create --repo nmorabowen/OpenSees --base ladruno --draft
--head wp/<n>-<slug>` (until #1172 merges its diff includes WP-262; say so in the body). Do not
run `build.bat` (rule 9): when you need a binary run, finish your code, write the test, and say in
your report exactly which tests the orchestrator must run after building. Report at the end: what
changed (files), the gates you ran with their exact numbers, what the orchestrator must build/run,
open questions, and any `LadrunoDual.h` primitive you needed (do not edit that header).

---

## WP-263 · D-RC — `LadrunoRCConcrete` dual tangent + the first convergence measurement

**Goal.** Give `LadrunoRCConcrete` / `LadrunoRCFiniteStrain` an exact algorithmic tangent by pushing
`Dual<6>` through `returnMap3D`, behind a new option, and produce the harness that measures whether
it cuts Newton iterations (plan decision D2). Model · effort: opus · high.

**What is known (fork map, `ladruno` @ 757c6845d).** `LadrunoRCKernel.h:656 returnMap3D` is a pure
function of total strain + committed thresholds (no local Newton). `tangentMode` 0 algorithmic /
1 secant / 2 numerical (`LadrunoRCKernel.h:522`; parser flags `LadrunoRCConcrete.cpp:133-134`,
`LadrunoRCFiniteStrain.cpp:126-127`; option table `SRC/material/LadrunoOptSpec.h:168`). Mode 0 is the
secant `W_B·C0` (`:960-972`) + three hand cross-terms (β-compression `:976-987`, tension-stiffening
pin `:988-1010`, interlock cap `≈:1012+`) and **omits** `dP/dε` projector rotation (`:643-653`, `:996`),
`d(dt_bar)/dε` and `dv_ci,max/dε_n` (`≈:1020-1025`). Mode 2 is a forward FD (`:948-958`), refused
under IMPL-EX (`:943-947`). Hazard: `eig3x3sym` (`:313`) at equal principal strains — the existing
`degen` branch. Existing gates: `tests/test_ladrunoRCConcrete_*_cpp.py` with `tests/_testbed/rc_*_gpp.cpp`
(header-only g++ checks — your bit-identity base), `test_ladrunoRCConcrete_implex.py`,
`test_ladrunoRCConcrete_wall.py`, `tests/_testbed/rc_wall_harness.py`, quirk
`LEGACY-cyclic-softening-rc-shell-wall-asdshellq4.md` (the wall stalls past first crack for a
STRUCTURAL reason — do not promise to fix that).

**Deliver.**
1. `returnMap3D` (and the helpers it calls) templated on `T`; the double path bit-identical (all
   `rc_*_gpp.cpp` outputs byte-identical vs pre-WP headers).
2. `tangentMode 3` = dual: the tangent is `∂σ/∂ε` from one `Dual<6>` pass, valid under IMPL-EX
   (differentiates the extrapolated map actually executed). Option `-dualTangent` in both parsers
   and `LadrunoOptSpec.h`. Default unchanged (mode 0).
3. A header-only oracle `tests/_testbed/rc_dual_gpp.cpp` + `tests/test_ladrunoRCConcrete_dual_cpp.py`:
   dual tangent vs central FD of `returnMap3D` (step plateau) on cracking, crushing, interlock-capped,
   tension-stiffened and degenerate-principal states; report the size of each term mode 0 omits
   (‖C_dual − C_mode0‖ / ‖C_dual‖ per state class); break-the-gate: mode 0 must differ from the dual
   by more than the oracle tolerance on at least the projector-rotation states.
4. **M1 harness** `tests/_testbed/rc_newton_count.py`: runs an existing RC deck (panel under shear,
   and the wall harness's first cycles) twice — default tangent vs `-dualTangent` — and prints
   Newton iterations per step, total iterations, failed steps and wall time, with `ladrunoBuild()`.
   The orchestrator runs it after building.
**Done when.** Bit-identity green; oracle green with the omitted-term table; the harness written and
dry-run-able; ledger fragments; draft PR with the numbers. **Not in scope:** making mode 3 the default
(that is decision D2, after the numbers).

---

## WP-264 · D-C3D — remove the production finite differences: Concrete3D drivers, RebarBuckling factor

**Goal.** Replace two finite differences shipped inside production tangents with dual numbers.
Model · effort: opus · high.

**What is known.** `LadrunoConcrete3DKernel.h:2227 dscalarDsig` central-differences each Voigt slot
of `scalarDriver(which, σ, mp)` (12 driver calls per gradient, each with an eigendecomposition),
step `DRIVER_FD_REL·S` (`:2246-2264`: the long unit-system history, #2621 / #2639); on the hydrostatic
axis `dscalarDsigAxis` (`:2273`) handles the cone kink; both reached via `driverGrad` (`:2555`). PR-877
measured these FDs as part of why Concrete3D runs 3–4× ASDConcrete3D. `LadrunoRebarBuckling.cpp:401-410`:
`dfactor/dε` by central FD (h = 1e-8) of a closed-form scalar factor. Byte-exact lanes:
`tests/_testbed/concrete3d_replay.cpp` with `concrete3d_replay_records*.txt`,
`concrete3d_kernel_check.cpp`, `test_ladrunoConcrete3D_*.py`. Quirks: `ADR-31-ladrunoconcrete3d-f2-port-damaged-tangent.md`,
`PR-877-*`, `LEGACY-spectral-eigenprojection-gradients-need-the.md`.

**Deliver.**
1. The scalar drivers templated on `T`. **Hazard:** if they go through `eig3sym`, repeated principal
   stresses (uniaxial, biaxial-equal) make eigenvector derivatives undefined. Find whether each driver
   depends on σ only through symmetric functions of the eigenvalues; if so, evaluate it (or its
   derivative) from invariants (I1, J2, θ) so the dual is well defined off the axis; keep
   `dscalarDsigAxis`'s dedicated treatment on the axis (a kink: one-sided by construction).
2. `driverGrad` uses the dual gradient off the axis. Keep the FD path selectable (a compile-time or
   option switch) for one WP so the replay comparison can be made.
3. Oracle: dual gradient vs the shipped FD and vs a step-plateau central FD, in MPa, kPa and Pa, over
   the replay record states; the damaged tangent with dual drivers vs with FD drivers (report the
   difference distribution); per-call time of `driverGrad` dual vs FD (g++ `-O2`, no `-march=native`).
4. **Byte-exact replay:** the stress path must not change. If the replay records include tangents, a
   tangent change is expected — re-record ONLY after the orchestrator and owner agree, and say so.
5. `LadrunoRebarBuckling` `dfactor/dε` → `Dual<1>` through the closed-form factor; a check vs the old FD.
**Done when.** No production FD left in those two places (except the documented axis branch);
oracle tables; replay stress byte-exact; timing table; ledger + a quirk if the eigen hazard bit.

---

## WP-265 · D-IMPLEX — tangent gates that bite

**Goal.** Today no test has to fail when a tangent is wrong: all 8 families of the mutation gate set
`"IDENT": 0.0` (`Ladruno_scripts/mutation_gate.py:64,79,88,109,116,…`; `SRC/Ladruno_mutation.h:27-73`
gates `getTangentStiff` / `getInitialStiff` / `getResistingForce*` in LadrunoBrick, Brick20, CPML,
ShellModifierSection and MPM). Give IDENT a real floor, and add the IMPL-EX tangent-identity gate.
Model · effort: opus · high.

**Deliver.**
1. Read the mutation framework end to end (`mutation_gate.py`, `Ladruno_mutation.h`,
   `tests/test_mutation_framework.py`, `tests/test_mutation_membrane.py`). Report why IDENT is 0:
   which tests run per family and why none sees a replaced tangent (expect: residual-only checks,
   or Newton that converges anyway with an identity tangent on linear problems).
2. Per family, a tangent-sensitive test that a correct tangent passes and an IDENT mutant fails:
   an element-level consistency check (assembled `K` vs a central FD of the resisting force at a
   deformed, nonlinear state — use the WP-158 quirk for probing from Python), or a Newton-iteration
   bound on a nonlinear deck. Set each family's IDENT floor from measurement (the orchestrator runs
   the mutation builds; give exact commands).
3. IMPL-EX identity gate: for every IMPL-EX material whose kernel is templated today (`LadrunoJ2`,
   `LadrunoUniaxialJ2`), assert at kernel level that the shipped IMPL-EX tangent equals the dual
   derivative of the extrapolated map with the clock latched (`LadrunoImplexClock.h` is read-only
   for you). For the non-templated IMPL-EX materials (`LadrunoCohesiveHinge*`, `ModIMKPeakOriented`,
   `Bilin`, RC, Concrete3D, HS, SANISAND) list how each could join once templated — do not template
   them (other WPs own those files).
**Done when.** Every family's IDENT floor > 0 by measurement; the IMPL-EX gate green for J2 and
UniaxialJ2 with a planted-mutant leg; ledger. **Files you own:** `mutation_gate.py`,
`SRC/Ladruno_mutation.h`, new `tests/test_*tangent_gate*`. Touch no material or element source.

---

## WP-266 · D-UP1 — the Jacobian blocks `LadrunoUP -geom hypo` drops

**Goal.** Measure, then make available, the full element Jacobian of the H8 u-p hypo lane by pushing
`Dual<32>` (24 u + 8 p DOFs) through the per-Gauss-point residual. Model · effort: opus · high.

**What is known.** `LadrunoUP.cpp:921-1030 updateHypo`: per GP the residual depends on u through the
porosity-evolved storage `1/Q̄(n(J))` (`:1001-1003`), the Kozeny–Carman permeability `kc(n(J))` floored
at 1e-6 (`:1004-1012`), the rotated permeability `R k̄ Rᵀ` (`:1013-1022`), the current-config Q/H/S
and spatial gradients through `F⁻¹` (`:976-998`). The tangent policy is explicit — "residual exact,
tangent dominant-term", rate-convection omitted (`:1102-1107`). The storage-coupling residual is the
INCREMENTAL form `Qᵀ(u − u_committed)/Δt` (ADR-78; do not differentiate a velocity form). Kernels:
`LadrunoUPKernel.h` (blocks, double, header-only, numpy oracle `tests/ladruno_up_reference.py`),
`SRC/element/solidTransformation/LadrunoHypoKernel.h`. Tests: `test_ladruno_up_element_hypo.py`,
`test_ladruno_up_kernel_cpp.py`. Research: `gitAPE/dual-numbers-opensees/research/04_porous_biot.md`.

**Deliver.**
1. Extract the hypo per-GP residual into a header-only templated function (new kernel header in
   `SRC/element/ladrunoUP/`), with the class calling its double instantiation; bit-identity of the
   double path (a verbatim-transcription check, as WP-262 did for UniaxialJ2).
   `LadrunoHypoKernel.h` may need templating for this — it is shared with LadrunoBrick: if you
   template it, its double path must stay bit-identical and you list every caller.
2. `tests/ladruno_up_hypo_dual_check.cpp` + pytest: the dual element Jacobian vs a central FD of the
   residual (step plateau) at large-strain, compacted (near the Kozeny–Carman floor) and rotated
   states; the block-by-block size of what the shipped tangent drops (`∂(H p)/∂u`, `∂(S ṗ)/∂u`,
   `∂(Q·)/∂u`, `∂(R k̄ Rᵀ)/∂u`, rate-convection) relative to the full Jacobian.
3. A binary deck for the orchestrator: large-strain consolidation (H8 hypo, Kozeny–Carman on) whose
   Newton iteration count can later be compared with a full-Jacobian option (the option itself is
   wave-2 P4 — do not add it now unless trivial; report).
**Done when.** Bit-identity green; dual Jacobian = FD at round-off; the dropped-block table; ledger.

---

## WP-267 · D-SENS2 — DDM sensitivities for `LadrunoCohesiveHinge` and `LadrunoBondSlip`

**Goal.** Close two more silent-zero DDM materials by the WP-262 recipe (quirk
`WP-262-uniaxial-ddm-stubs-silent-zero.md`). Model · effort: opus · medium.

**What is known.** `SRC/material/uniaxial/LadrunoCohesiveHinge.cpp` (678 lines): analytic envelope
tangent on loading (`:379`, may be negative), secant-to-origin on unload (`:381`), IMPL-EX `Ksec_tilde`
(`:343-364`, `LadrunoImplexClock.h`), parameters at `:653-660`. `LadrunoBondSlip.cpp` (317 lines):
envelope tangent with a `−1e-12·k0` floor (`:182-183`), elastic unload/reload (`:192`); WP bondslip
reversal slider (PR #1115). Tests: `test_ladrunoCohesiveHinge_material.py`,
`test_ladruno_cohesive_hinge_implex_bk.py`, `test_ladrunoBondSlip_material.py`.

**Deliver, per material (CohesiveHinge first).**
1. Extract the trial-state update into a templated kernel header (`<Material>Kernel.h`, `T` state, `P`
   constants); the class calls its double instantiation; a verbatim-transcription bit-identity check.
2. `getStressSensitivity`, `getTangentSensitivity`, `getInitialTangentSensitivity`, `commitSensitivity`
   by one `Dual<1>` pass seeded on the active parameter and on stored history sensitivities (`SHV`,
   copied in `getCopy`, cleared in `revertToStart`). Under IMPL-EX the extrapolation uses committed
   history — seed its sensitivity too.
3. Kernel check: DDM through history vs FD of a whole cyclic path per parameter; a NOHIST mutant fails.
4. A binary end-to-end test the orchestrator runs: a `zeroLength` (BondSlip) or a hinge-carrying beam
   (CohesiveHinge in a `LadrunoForceBeamColumn` or `zeroLength` rotational spring) under cyclic load
   control, `sensNodeDisp` vs FD, with a Steel01 control leg in the same harness
   (`tests/test_ladrunoUniaxialJ2_ddm_fiber.py` is the template).
**Not in scope:** `LadrunoRebarBuckling` (WP-264 owns its file this wave), `LadrunoCohesiveHingeBiaxial`
(nD), the damage-free extras. **Done when.** Both materials pass all four legs; ledger; the quirk's
"still stubbed" list updated (edit that fragment).
