---
wp: WP-225
title: "`ubuntu-latest` moved to 26.04 / GCC 15.2 under an unchanged workflow: 16 Zone A tests failed that no PR touched"
date: 2026-11-06
---
### `ubuntu-latest` moved to 26.04 / GCC 15.2 under an unchanged workflow: 16 Zone A tests failed that no PR touched

- **Bites:** Zone A on PR #1063 (WP-224, a header-only change to the refusal
  counter) failed 16 tests in code it does not reach: plateRebar bit-exactness,
  Brick20 eigenvalues, LadrunoRigidBody NaN, SANISAND CPPM baselines, CPML at
  1.43e-10, ADR-97 inertness. The last green `ladruno` run (42041136405) ran the
  same tests on ubuntu-24.04; the failing run's log names image ubuntu-26.04 and
  GCC 15.2.0. In the same run the Concrete3D replay job, pinned to
  `ubuntu-24.04`, passed. Every PR and push after the move would have failed
  the same way.
- **Why:** `runs-on: ubuntu-latest` is a label GitHub moves without a commit in
  the repository. Bit-exact and tight-baseline tests are compiler-sensitive (a
  new GCC reorders or contracts floating-point operations), so the label move
  reads as a regression in whatever PR runs first.
- **Fix/status:** WP-225 pins every fork workflow job to `ubuntu-24.04`
  (`ladruno.yml` static gates, MKL compat and Zone A; `ladruno_mutation.yml`;
  `ladruno_release.yml`, whose Linux tarball would otherwise need the newer
  glibc to load). Quirk-lint rule `floating-runner` fails a `*-latest` image in
  any `.github/workflows/ladruno*.yml`. Moving to 26.04 is its own PR that
  re-baselines the 16 tests and checks the two NaN / domain-error failures.
