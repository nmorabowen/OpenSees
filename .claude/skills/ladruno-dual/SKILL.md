---
name: ladruno-dual
description: >
  Run and implement the Ladruno dual-number program (forward-mode automatic differentiation,
  WP-262 onward): LadrunoDual.h, kernels templated on the scalar, dual tangent oracles, dual
  production tangents replacing finite differences / secants / stubs, and DDM sensitivities by
  parameter-seeded duals. Use when the user runs "/ladruno-dual <wp>" (262, 263, 264, 265, 266,
  267, …) or asks to start, continue, orchestrate, review or unblock a dual-number work package;
  and before templating any fork kernel on a scalar type, writing a dual tangent oracle, replacing
  a numerical tangent with a dual one, or implementing getStressSensitivity / commitSensitivity
  with duals. Holds the work-package map, the roster, the protocol, the worker rules and the
  dispatch briefs.
---

# `/ladruno-dual <wp>` — run one dual-number work package

**Authority:** the plan `Ladruno_implementation/_wp262_dual_orchestration_plan.md` and the
program issue ("Program: dual numbers (WP-262…)", label `program`). This skill is the operating
layer. When they disagree, the plan wins; fix the skill in the same PR.

**Two modes.**
- **Orchestrator:** the session that runs `/ladruno-dual <wp>`. It plans, claims, dispatches,
  launches builds, verifies and records. It does not edit a worker's files while that worker runs.
- **Worker:** an agent dispatched with a brief from [references/wp-briefs.md](references/wp-briefs.md).
  It implements one WP in its own worktree under the worker rules below.

## The map (plan §1)

```
WAVE 0  WP-262  LadrunoDual.h + J2/Lemaitre kernels + UniaxialJ2 DDM + this skill   (PR #1172)
WAVE 1  WP-263 D-RC      RC returnMap3D Dual<6> tangent + Newton-iteration measurement  -> decision D2
        WP-264 D-C3D     Concrete3D FD driver gradients -> dual; RebarBuckling dfactor -> Dual<1>
        WP-265 D-IMPLEX  tangent gates that bite: IMPL-EX identity + non-zero mutation IDENT floor
        WP-266 D-UP1     LadrunoUP hypo GP residual templated; Dual<32> element-Jacobian oracle
        WP-267 D-SENS2   DDM for LadrunoCohesiveHinge + LadrunoBondSlip (WP-262 recipe)
WAVE 2  D-HS, D-SAN (after D2), D-COROT, D-MORTAR, D-UPDDM, apeGmsh G1     (numbers claimed later)
WAVE 3  decided at the wave-2 close
```

- Wave-1 branches are cut from **WP-262's head** (they need `LadrunoDual.h`). Each PR is opened
  with `--base ladruno`, never `--base wp/262-…`. Until #1172 merges, the PR diff includes WP-262's
  commits. After it merges: `git merge origin/ladruno`.
- **File ownership in wave 1:** 263 owns `LadrunoRCKernel.h`, `LadrunoRCConcrete.*`,
  `LadrunoRCFiniteStrain.*`. 264 owns `LadrunoConcrete3DKernel.h`, `LadrunoConcrete3D.*`,
  `LadrunoRebarBuckling.*`. 265 owns `Ladruno_scripts/mutation_gate.py`, `SRC/Ladruno_mutation.h`
  and new `tests/test_*tangent_gate*`. 266 owns `SRC/element/ladrunoUP/*`. 267 owns
  `LadrunoCohesiveHinge.*` (uniaxial) and `LadrunoBondSlip.*`.
  **`LadrunoDual.h` belongs to nobody in wave 1:** a worker that needs a new primitive reports it
  and the orchestrator adds it in a separate small commit.
- `LadrunoOptSpec.h` is shared. A worker that adds an option there touches only its own lines and
  says so in its report.

## Roster (plan §4) — set `model` and `effort` on every Agent call

| Role | `model` · `effort` | WPs |
|---|---|---|
| Orchestrator | this session (Opus) or Fable · high; max at D2 / D3 | all |
| Kernel templating + oracle authors | opus · high | 263, 264, 265, 266 |
| Recipe ports | opus · medium | 267 |
| Harnesses | sonnet · medium | M0, M1 |
| Critic at ready-flip | the other family: fable · high for Opus authors | every WP |
| Read-only sweeps | `Explore`, "very thorough" | any |

## Orchestrator loop

1. **Boot (≤ 8 tool calls):**
   - `ListAgents`, and `gh issue list --label in-flight`, for duplicates;
   - the program issue with its comments;
   - the WP's draft PR and its last `### Handoff`;
   - `git fetch origin`.
2. **Claim:**
   - post the WP number on the program issue;
   - branch `wp/<n>-<slug>` from the base named in the brief;
   - open a draft PR, `--base ladruno`, on day one.
3. **Dispatch:**
   - copy the brief, then fill in the PR number and the last handoff's **Next**;
   - `subagent_type: general-purpose`, the roster's model and effort, `isolation: "worktree"`;
   - independent workers go in one message; never poll them.
4. **Verify — never accept a worker's word:**
   - re-run the worker's g++ checks and the static gates yourself;
   - **you launch every MSVC build**, in the worker's worktree:
     `cmd /c "call Ladruno_scripts\setup_env.bat && Ladruno_scripts\build.bat" > <log>` in the background;
   - run binary tests with `py -3.12 -S Ladruno_scripts/ci_run_pytest.py -- <test files>`, from `tests/`;
   - check `ladrunoBuild` = the worker's head;
   - every reported number names the build hash, the deck, the output path and the date.
5. **Review:** one out-of-family critic at ready-flip, briefed with the diff and the plan row; then `/code-review high` for C++.
6. **Land:**
   - Zone-A on demand with `gh workflow run ladruno.yml --ref <branch>`;
   - read CI with `ccd_pr` and Auto-fix;
   - **merges are the owner's** (plan §8 item 5, 2026-11-28).
7. **Record:**
   - a handoff on the PR;
   - the board row on #937 and the program issue refreshed;
   - measurements go into the research bench (`gitAPE/dual-numbers-opensees/research/`).

## Worker rules (verbatim in every brief; plan §6)

1. **Template, don't fork.** Write the kernel once as `template <class T>`, plus a `P` for the
   constants when parameters are seeded. Keep the `double` entry point with its old name and
   signature, and keep every caller unchanged.
2. **Bit-identity is the first gate.** Before any derivative work, prove the `double` path
   unchanged, in one of two ways:
   - diff the existing C++ checks' `%.17g` output against the pre-WP headers
     (`git show origin/ladruno:<path>` copies on a first `-I`);
   - or compare bit for bit against a verbatim transcription of the old code inside your check
     (the WP-262 `tests/ladruno_uniaxialj2_dual_check.cpp` pattern).

   Byte-exact replay lanes, Concrete3D's above all, must stay byte-exact.
3. **Branches and tolerances read `value(x)`.** There is no implicit `Dual → double`; never add one.
   Mixed `double`/`Dual` operators exist in `LadrunoDual.h`.
4. **A local Newton gets the IFT step:** after convergence, one Newton step with its value
   zeroed: `T c = R/dR; ladruno_dual::zero_value(c); x -= c;` inside
   `if constexpr (ladruno_dual::is_dual<T>::value)`. A vector local solve uses the converged
   Jacobian. Quirk: `WP-262-dual-derivative-lags-local-newton.md`.
5. **Every oracle has a break-the-gate leg:** a mutant it must flag, such as a symmetrized
   tangent, dropped history sensitivities or a dropped term.
6. **Kinks keep their dedicated handling:** apex, the hydrostatic axis, repeated eigenvalues,
   `sqrt(0)`, `acos(±1)`. Never push a dual through an eigen-solver at a degenerate state. Prefer
   invariant forms or the closed-form isotropic-function derivative
   (`LogStrainKernel.h::isoFunction` handles degeneracy).
7. **Production replaces only an FD, a secant or a stub.** An analytic tangent that passes the
   oracle stays. A new production tangent ships behind an option first, and the WP reports its
   numbers.
8. **Zero vanilla files**, unless your brief says otherwise.
9. **Do not run `build.bat`.** The orchestrator builds. Compile your header-only checks with g++
   (`-O2 -std=c++17`, `*cxx_runtime_flags()` in the pytest wrapper; **no `-march=native`**).
   Write the binary-level tests; the orchestrator runs them.
10. **Ledger:**
    - `Ladruno_implementation/ledger/implementations/WP-<n>-<slug>.md`, plus a quirk fragment
      for every trap;
    - never edit `LEDGER_*.md`;
    - run `python ci/check_ledger_fragments.py`, `python ci/check_quirk_patterns.py`,
      `python ci/check_duplicate_defs.py` and `python Ladruno_scripts/stamp_headers.py --check`
      (a new fork file goes into its `GLOBS`).
11. **Commit on your WP branch only**, staged by path, with the trailer lines the session gives
    you. Push, and keep the draft PR body current. Finish with a `### Handoff` comment on the PR
    (template: `Ladruno_internal/ISSUES.md` §3).
12. Read code from `origin/ladruno` (or your base) when in doubt. Quote the code you make claims
    about, with `file:line`.
