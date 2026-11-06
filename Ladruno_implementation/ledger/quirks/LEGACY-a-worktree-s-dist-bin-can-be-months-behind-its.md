---
wp: LEGACY
title: "A worktree's dist/bin can be MONTHS behind its branch, and an \"engine guard\" that checks the module's PATH will happily wave it through"
legacy_seq: 254
---
### A worktree's `dist/bin` can be MONTHS behind its branch, and an "engine guard" that checks the module's PATH will happily wave it through
- **Bites:** any harness that pins itself to a worktree build. The `hypo_bearing`
  runner already carries a guard (scoping finding 5) against the installed
  Ladruno's site `.pth` pre-importing `opensees` — but that guard asserts
  `os.path.dirname(ops.__file__) == dist/bin`, i.e. *where* the module loaded
  from. It says nothing about *what is in it*. A fresh worktree checked out at
  a branch with ADR-78/79 merged had a `dist/bin/opensees.pyd` from an earlier
  build, which passed the location guard and then refused the feature under
  test: `-geom 'corot' not supported: only 'linear' is accepted (the axis is
  reserved ... ADR 71 §2.4)` — an ADR-71-era binary answering for an
  ADR-79-era branch.
- **Tell:** a capability error naming an ADR *older* than your branch, on a
  worktree you never built in. `git log` looks right; the `.pyd` mtime predates
  the feature commits.
- **Rule:** a location guard is necessary but not sufficient — assert
  CAPABILITY. Cheapest form is to construct the thing under test (one element
  with the flags the campaign needs) before committing hours to a run, which is
  also what catches a stale build in a *shared* checkout that another agent
  rebuilt on a different branch. Rebuild the worktree (`build.bat OpenSeesPy`)
  and verify `dist/bin/opensees.pyd` mtime moved. Cross-ref
  [[ladruno-build-in-worktree-not-shared-checkout]]. *2026-07-30 (ADR-79 locking leg).*
