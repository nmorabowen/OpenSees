# `ci/` — Ladruno test-bed gates

Dependency-light gates run by `.github/workflows/ladruno.yml`. All are runnable
locally from the repo root. Source of truth:
[`Ladruno_implementation/testbed/00_canonical_testbed.md`](../Ladruno_implementation/testbed/00_canonical_testbed.md).

| Script | Gate | Deps | Exit |
|---|---|---|---|
| `check_classtags.py` | classTag value collisions (Axis 1, G2) + cross-header drift + ladruno-band policy | none | 1 on a ladruno-involved collision |
| `check_manifest.py` | every ladruno classTag has a manifest row + a test (or WAIVED) | PyYAML | 1 on an unaccounted/active-but-untested tag |
| `check_tcl_results.py` | turns a `FAILED` line in `results.out` into a nonzero exit (G1) | none | 1 if any FAILED |
| `check_quirk_patterns.py` | `LEDGER_quirks` entries with a greppable pattern, enforced on fork-stamped sources: L1 Rayleigh snapshot (#562), L2 singleton reset on `wipe`, L3 task-guide pointers resolve (WP-115), L4 an Element subclass's `commitState()` chains to `Element::commitState()` so `betaKc`'s Kc is refreshed (WP-118), L5 no element subtracts its load vector in both `getResistingForce()` and the `getResistingForceIncInertia()` that calls it — the ground-motion load counted twice (WP-119; scans ALL element files, vanilla included), L6 the ground-motion inertia load reaches the residual as +M·R·a_g — accumulation sign × application sign must be +1 (WP-117; all element files; unreadable signs are skipped, never guessed). L9 no declaration as the whole unbraced body of an `if`/`else`/`for`/`while` that shadows an outer variable — `ManzariDafalias::ForwardEuler`'s `Vector r` (WP-158; all SRC files). L10 a fork-stamped integrator with its own `Vector*` march state that inherits `IncrementalIntegrator`'s no-op `revertToLastStep()` — a retried step marches from uncommitted state (WP-153, LadrunoDynamicRelaxation #899). Self-test: `test_check_quirk_patterns.py` | none | 1 on any unwaived finding |
| `check_viewer_ledger.py` | V1: a change (`base...head`) that brings a new source file into `Ladruno_tools/<tool>/` (added, copied, or moved in from elsewhere) must add a `LEDGER_implementations.md` line naming `Ladruno_tools/<tool>` — that tool's row; another row, a whitespace-only edit or a deleted ledger does not count (WP-121; #35, #53, #485, #487 did not). Runs on every event: the PR range on a PR, `origin/ladruno...HEAD` otherwise. Modification-only changes are not checked. Self-test: `test_check_viewer_ledger.py` | git | 1 on a finding, 2 if the range cannot be evaluated, 3 if git cannot run |
| `check_mkl_compat.py` | The MKL-gated sources build against the OLDEST MKL a fork build uses, esmeralda's oneMKL 2024.2.2 (WP-148; #864 used a 2025.0 macro, Windows built, Linux did not, #886). Needs a build dir configured with `-DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DLADRUNO_MKL_PARDISO_LINUX=ON -DLADRUNO_MKL_FEAST_LINUX=ON -DMKL_RT_HINT=<mkl>/lib`; builds nothing else. M1: every unit that touches MKL (includes an MKL header or `ProfilerRunMeta.h`, or declares an MKL routine) compiles with its CMake command. M2: every MKL symbol those objects reference is defined by the linked MKL layers (`nm`), which is the only way to see a too-new routine in the FEAST sources, since they declare MKL routines by hand. Refuses to pass on nothing: wrong MKL version, or a required unit missing or compiled without its macro. `--self-test` proves M1 (on the #864 incident) and M2 still fire. CI job `mkl-compat` gets the headers and libraries from the PyPI `mkl`/`mkl-include` wheels, unzipped rather than pip-installed. | a configured build dir, gcc, `nm` | 1 on a finding, 2 if it cannot evaluate |

```bash
python ci/check_classtags.py            # default: actionable only
python ci/check_quirk_patterns.py      # L1+L2+L3; --only L1 / --root DIR / --list-waivers
python ci/check_classtags.py --verbose  # also list inherited-upstream collisions
python ci/check_classtags.py --strict   # warnings become errors
python ci/check_manifest.py
python ci/check_viewer_ledger.py        # this branch vs origin/ladruno; --base/--head for any range
# Linux, after configuring build/mkl as in the mkl-compat job of ladruno.yml:
python ci/check_mkl_compat.py --build-dir build/mkl --self-test
# after running the Tcl suite into EXAMPLES/verification/results.out:
python ci/check_tcl_results.py
```

**Why these exist:** the three defect classes they catch are invisible to every
runtime test and shipped past human review in this very repo — a 205-line
`classTags.h` drift, value-collision hacks, and a Tcl "test" that can't fail CI.
Discipline alone demonstrably failed here; these make the rules machine-checked.
