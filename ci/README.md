# `ci/` — Ladruno test-bed gates

Dependency-light gates run by `.github/workflows/ladruno.yml`. All are runnable
locally from the repo root. Source of truth:
[`Ladruno_implementation/testbed/00_canonical_testbed.md`](../Ladruno_implementation/testbed/00_canonical_testbed.md).

| Script | Gate | Deps | Exit |
|---|---|---|---|
| `check_classtags.py` | classTag value collisions (Axis 1, G2) + cross-header drift + ladruno-band policy | none | 1 on a ladruno-involved collision |
| `check_manifest.py` | every ladruno classTag has a manifest row + a test (or WAIVED) | PyYAML | 1 on an unaccounted/active-but-untested tag |
| `check_tcl_results.py` | turns a `FAILED` line in `results.out` into a nonzero exit (G1) | none | 1 if any FAILED |

```bash
python ci/check_classtags.py            # default: actionable only
python ci/check_classtags.py --verbose  # also list inherited-upstream collisions
python ci/check_classtags.py --strict   # warnings become errors
python ci/check_manifest.py
# after running the Tcl suite into EXAMPLES/verification/results.out:
python ci/check_tcl_results.py
```

**Why these exist:** the three defect classes they catch are invisible to every
runtime test and shipped past human review in this very repo — a 205-line
`classTags.h` drift, value-collision hacks, and a Tcl "test" that can't fail CI.
Discipline alone demonstrably failed here; these make the rules machine-checked.
