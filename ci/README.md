# `ci/` — Ladruno test-bed gates

Dependency-light gates run by `.github/workflows/ladruno.yml`. All are runnable
locally from the repo root. Source of truth:
[`Ladruno_implementation/testbed/00_canonical_testbed.md`](../Ladruno_implementation/testbed/00_canonical_testbed.md).

| Script | Gate | Deps | Exit |
|---|---|---|---|
| `check_classtags.py` | classTag value collisions (Axis 1, G2) + cross-header drift + ladruno-band policy | none | 1 on a ladruno-involved collision |
| `check_manifest.py` | every ladruno classTag has a manifest row + a test (or WAIVED) | PyYAML | 1 on an unaccounted/active-but-untested tag |
| `check_tcl_results.py` | turns a `FAILED` line in `results.out` into a nonzero exit (G1) | none | 1 if any FAILED |
| `check_quirk_patterns.py` | `LEDGER_quirks` entries with a greppable pattern, enforced on fork-stamped sources: L1 Rayleigh snapshot (#562), L2 singleton reset on `wipe`, L3 task-guide pointers resolve (WP-115), L4 an Element subclass's `commitState()` chains to `Element::commitState()` so `betaKc`'s Kc is refreshed (WP-118), L5 no element subtracts its load vector in both `getResistingForce()` and the `getResistingForceIncInertia()` that calls it — the ground-motion load counted twice (WP-119; scans ALL element files, vanilla included), L6 the ground-motion inertia load reaches the residual as +M·R·a_g — accumulation sign × application sign must be +1 (WP-117; all element files; unreadable signs are skipped, never guessed). Self-test: `test_check_quirk_patterns.py` | none | 1 on any unwaived finding |
| `check_viewer_ledger.py` | V1: a change (PR diff, `base...head`) that adds a file under `Ladruno_tools/` must also edit `LEDGER_implementations.md` (WP-121; #35, #53, #485, #487 did not). Runs on `pull_request` only; modification-only changes are not checked. Self-test: `test_check_viewer_ledger.py` | git | 1 on a finding, 2 if the range does not resolve |

```bash
python ci/check_classtags.py            # default: actionable only
python ci/check_quirk_patterns.py      # L1+L2+L3; --only L1 / --root DIR / --list-waivers
python ci/check_classtags.py --verbose  # also list inherited-upstream collisions
python ci/check_classtags.py --strict   # warnings become errors
python ci/check_manifest.py
python ci/check_viewer_ledger.py        # this branch vs origin/ladruno; --base/--head for any range
# after running the Tcl suite into EXAMPLES/verification/results.out:
python ci/check_tcl_results.py
```

**Why these exist:** the three defect classes they catch are invisible to every
runtime test and shipped past human review in this very repo — a 205-line
`classTags.h` drift, value-collision hacks, and a Tcl "test" that can't fail CI.
Discipline alone demonstrably failed here; these make the rules machine-checked.
