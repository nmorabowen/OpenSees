# `ci/` — Ladruno test-bed gates

Dependency-light gates run by `.github/workflows/ladruno.yml`. All are runnable
locally from the repo root. Source of truth:
[`Ladruno_implementation/testbed/00_canonical_testbed.md`](../Ladruno_implementation/testbed/00_canonical_testbed.md).

| Script | Gate | Deps | Exit |
|---|---|---|---|
| `check_classtags.py` | classTag value collisions (Axis 1, G2) + cross-header drift + ladruno-band policy | none | 1 on a ladruno-involved collision |
| `check_manifest.py` | every ladruno classTag has a manifest row + a test (or WAIVED) | PyYAML | 1 on an unaccounted/active-but-untested tag |
| `check_tcl_results.py` | turns a `FAILED` line in `results.out` into a nonzero exit (G1) | none | 1 if any FAILED |
| `check_quirk_patterns.py` | `LEDGER_quirks` entries with a greppable pattern, enforced as named rules (slugs; the old L-numbers are deprecated `--only` aliases for one release, WP-162). The rules are the generated table under [Quirk-lint rules](#quirk-lint-rules). Self-tests: one `test_quirk_<slug>.py` per rule + `test_quirk_rules.py` (registry, CLI, this table) | none | 1 on any unwaived finding, 2 on an unknown rule |
| `check_viewer_ledger.py` | V1: a change (`base...head`) that brings a new source file into `Ladruno_tools/<tool>/` (added, copied, or moved in from elsewhere) must add a `LEDGER_implementations.md` line naming `Ladruno_tools/<tool>` — that tool's row; another row, a whitespace-only edit or a deleted ledger does not count (WP-121; #35, #53, #485, #487 did not). Runs on every event: the PR range on a PR, `origin/ladruno...HEAD` otherwise. Modification-only changes are not checked. Self-test: `test_check_viewer_ledger.py` | git | 1 on a finding, 2 if the range cannot be evaluated, 3 if git cannot run |
| `check_mkl_compat.py` | The MKL-gated sources build against the OLDEST MKL a fork build uses, esmeralda's oneMKL 2024.2.2 (WP-148; #864 used a 2025.0 macro, Windows built, Linux did not, #886). Needs a build dir configured with `-DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DLADRUNO_MKL_PARDISO_LINUX=ON -DLADRUNO_MKL_FEAST_LINUX=ON -DMKL_RT_HINT=<mkl>/lib`; builds nothing else. M1: every unit that touches MKL (includes an MKL header or `ProfilerRunMeta.h`, or declares an MKL routine) compiles with its CMake command. M2: every MKL symbol those objects reference is defined by the linked MKL layers (`nm`), which is the only way to see a too-new routine in the FEAST sources, since they declare MKL routines by hand. Refuses to pass on nothing: wrong MKL version, or a required unit missing or compiled without its macro. `--self-test` proves M1 (on the #864 incident) and M2 still fire. CI job `mkl-compat` gets the headers and libraries from the PyPI `mkl`/`mkl-include` wheels, unzipped rather than pip-installed. | a configured build dir, gcc, `nm` | 1 on a finding, 2 if it cannot evaluate |
| `check_ledger_fragments.py` | The build-control ledgers are one fragment per entry under `Ladruno_implementation/ledger/<kind>/` (WP-161), so no two PRs write the same file. F1 schema (front matter, body shape per kind), F2 file name `<wp>-<slug>.md` matches `wp`, F3 no duplicate name/body/quirk heading, F4 `class_tags` vs `SRC/classTags.h` (via `check_classtags.parse`), F5 `banner` is a `banner_features.txt` line, F6 templates + `ledger.py build` renders, F7 the committed `LEDGER_*.md` stubs are untouched. Self-test (incl. the round-trip proof on the real ledgers): `test_ledger.py` | none | 1 on an error |
| `ledger.py` | Not a gate: `build` renders `LEDGER_*.md` into `Ladruno_implementation/ledger/_build/` (gitignored; the static-gates job uploads it as the `ledgers` artifact), `split` is the one-time migration, `migrate --from-diff <base>` turns a pre-WP-161 branch's ledger edits into fragments, `roundtrip` proves build(split(ledgers)) == ledgers | none (git for `migrate`/`roundtrip`) | 2 on an error |

```bash
python ci/check_classtags.py            # default: actionable only
python ci/check_quirk_patterns.py      # all rules; --only dead-decl,rayleigh / --root DIR / --list-waivers
ruff check --select F811 ci tests Ladruno_scripts   # duplicate definitions (pip install ruff==0.16.10)
python ci/check_duplicate_defs.py      # D1: a def/class defined twice in one scope (what F811 misses)
python ci/check_classtags.py --verbose  # also list inherited-upstream collisions
python ci/check_classtags.py --strict   # warnings become errors
python ci/check_manifest.py
python ci/check_viewer_ledger.py        # this branch vs origin/ladruno; --base/--head for any range
# Linux, after configuring build/mkl as in the mkl-compat job of ladruno.yml:
python ci/check_mkl_compat.py --build-dir build/mkl --self-test
# after running the Tcl suite into EXAMPLES/verification/results.out:
python ci/check_tcl_results.py
python ci/check_ledger_fragments.py    # ledger fragments (WP-161)
python ci/ledger.py build              # -> Ladruno_implementation/ledger/_build/LEDGER_*.md
```

**Why these exist:** the three defect classes they catch are invisible to every
runtime test and shipped past human review in this very repo — a 205-line
`classTags.h` drift, value-collision hacks, and a Tcl "test" that can't fail CI.
Discipline alone demonstrably failed here; these make the rules machine-checked.

## Duplicate definitions (WP-162)

`ci/check_duplicate_defs.py` (D1) fails on a `def`/`class` defined twice as direct statements of
one scope in `ci/`, `tests/`, `Ladruno_scripts/` — Python keeps the last silently. The static-gates
step also runs `ruff check --select F811` (pinned `ruff==0.16.10`), which catches a duplicated test
name or re-import but NOT the incident: on the real #899/#901 merge (two `def _l9`, the first one
referenced by tests in between) F811 passes. Self-test: `test_duplicate_defs.py`. Waive with
`# ladruno-lint: redef-ok <reason>`. Exit 1 on a finding.

## Quirk-lint rules

Named by slug (WP-162): a new rule takes a new slug in `RULES` (`ci/check_quirk_patterns.py`) and its own `ci/test_quirk_<slug>.py`, never "the next free" L-number — two PRs both took L9 once, and their test helpers merged silently, the later `def _l9` shadowing the earlier. `ruff check --select F811` (a static gate) now fails on that shape. This table is generated; `test_quirk_rules.py` fails if it drifts.

<!-- quirk-rules:start (generated: python ci/check_quirk_patterns.py --rules-table) -->
| Rule | Alias (deprecated) | Waiver | Scope | Catches |
|---|---|---|---|---|
| `rayleigh` | L1 | `// ladruno-lint: rayleigh-ok <reason>` | fork-stamped | Rayleigh forces accumulated into a buffer that is not a function-local vector seeded before the first `getRayleighDampingForces()` call (#562) |
| `wipe` | L2 | `// ladruno-lint: wipe-ok <reason>` | fork-stamped | a process-wide singleton whose state is not reset on `wipe` |
| `pointers` | L3 | — | task guides | a `Quirks: "..."` pointer in `.claude/skills/*/SKILL.md` that no longer matches the quirks ledger (WP-115) |
| `commit` | L4 | `// ladruno-lint: commit-ok <reason>` | fork-stamped | an Element subclass's `commitState()` that does not chain to `Element::commitState()`, so betaKc's Kc is never refreshed (WP-118) |
| `double-load` | L5 | `// ladruno-lint: double-ok <reason>` | all element files | an element that subtracts its load vector in both `getResistingForce()` and the `getResistingForceIncInertia()` that calls it (WP-119) |
| `ground-sign` | L6 | `// ladruno-lint: sign-ok <reason>` | all element files | the ground-motion inertia load reaching the residual with the wrong sign; unreadable signs are skipped (WP-117) |
| `sequence` | L7 | `// ladruno-lint: sequence-ok <reason>` | all SRC files | arithmetic combining two element-accessor calls in one statement: unspecified call order over shared storage (WP-124 C15, WP-140) |
| `ci-coverage` | L8 | `# ci-coverage: <kind> <reason>` | zone_a tests | a `zone_a` test that branches on the platform without `# ci-coverage: <kind> <reason>` (WP-143) |
| `dead-decl` | L9 | `// ladruno-lint: decl-ok <reason>` | all SRC files | a declaration as the whole unbraced body of an if/else/for/while that shadows an outer variable, `ManzariDafalias::ForwardEuler`'s `Vector r` (WP-158) |
| `revert` | L10 | `// ladruno-lint: revert-ok <reason>` | fork-stamped | a fork integrator with its own `Vector*` march state that inherits the no-op `revertToLastStep()` (WP-153, #899) |
| `unknown-token` | — | `// ladruno-lint: unknown-ok <reason>` | fork-stamped | a fork parser that lets an unknown token through: a "unknown tokens are ignored" comment, or an option ladder ending its token loop with no final `else` (WP-167) |
<!-- quirk-rules:end -->
