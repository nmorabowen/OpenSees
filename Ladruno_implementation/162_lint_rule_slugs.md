# WP-162 — lint rules by slug, one test file per rule, a duplicate-definition gate

Status: DRAFT PR #916 (owner merges after #894). Sibling: WP-161 (#915, per-WP ledger
fragments); the two touch disjoint lines and merge cleanly in either order.

## The problem

- WP-153 and WP-158 both took quirk rule **L9** — "the next free number" —
  and both added `def _l9(...)` to `ci/test_check_quirk_patterns.py`. The #899
  merge-up of #901 merged that file with **no conflict and two `def _l9`**: every
  WP-153 test then called WP-158's helper. It was caught by eye and renumbered
  L10 by hand.
- One shared self-test file and one numbering sequence mean every PR that adds
  a rule edits the same lines. No Python linter ran in CI.

## Change

- **Slugs.** `ci/check_quirk_patterns.py` has ONE registry, `_RULE_LIST` ->
  `RULES` (`rayleigh`, `wipe`, `pointers`, `commit`, `double-load`,
  `ground-sign`, `sequence`, `ci-coverage`, `dead-decl`, `revert`), with
  import-time asserts: unique slug, unique alias, unique check function,
  unique waiver token, and the `WAIVER` regex agrees with the registry.
  Findings print with the slug. `--only` takes slugs; the L-numbers are
  accepted as deprecated aliases for one release (a stderr note).
  Waiver syntax is unchanged (`// ladruno-lint: <token> <reason>`).
- **One test file per rule**: `ci/test_quirk_<slug>.py` (104 cases, the same
  104 as before), shared helpers in `ci/_quirk_testkit.py` (constants and pure
  functions only). `ci/test_quirk_rules.py` checks the registry, the CLI and
  that the `ci/README.md` rule table equals `--rules-table` (generated).
- **Duplicate-definition gate**, its own static-gates step before the classTag
  step: `ci/check_duplicate_defs.py` (D1) + `ruff check --select F811`
  (`ruff==0.16.10`) over `ci tests Ladruno_scripts`. Both clean on `ladruno`
  9fcb6cfaa (D1: 525 files).

## Where the design changed, and why

**F811 alone does not catch the incident.** Pyflakes F811 is "redefinition of
UNUSED name"; in the merged file the first `_l9` is referenced by the WP-153
tests between the two definitions, so it counts as used. Reproduced on the real
three-way merge (`git merge-file` of e2560dd2b / merge-base / 5300da720, the
#899 merge-up): two `def _l9`, and ruff 0.16.10 `--select F811` reports "All
checks passed". Ruff has no function-redefined rule (pylint E0102). So D1, a
40-line `ast` check: a `def`/`class` defined twice as direct statements of one
scope fails, unless the definitions are alternatives (inside `if`/`try`),
`@overload`, a property accessor or `@x.register`. F811 stays in the step for
what it does catch (a duplicated test name, a re-import); a test pins that it
still misses the incident, so a future ruff that catches it shows up.

## Verification

See the PR body.
