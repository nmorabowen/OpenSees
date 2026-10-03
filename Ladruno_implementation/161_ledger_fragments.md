# WP-161 — one ledger fragment per WP; the ledgers are generated

Status: DRAFT PR (owner merges after #894). Sibling: WP-162 (quirk-lint rule
slugs, per-rule tests, ruff F811), which touches none of the same lines.

## The problem (measured over the last 60 PR merges into `ladruno`)

- 77 branch merge-ups; 64 conflicted, **55 of them on the ledgers only**.
- Touch rates: `LEDGER_quirks.md` 50/60 PRs, `LEDGER_implementations.md` 45/60,
  `LEDGER_vanilla_files.md` 35/60. Implementations rows go in at the TOP of the
  table, quirks and vanilla rows at the END: every PR writes at the same anchor.
- `merge=union` (in `.gitattributes` since 114f609b5) does not help where it
  matters: GitHub's server-side merge ignores merge drivers, so the PR reads
  CONFLICTING, no workflow runs, and the author does a local merge-up, an
  hour-long Windows rebuild and Zone-A again (WORKFLOW_GOTCHAS §9). Locally,
  union turns in-place edits into silent duplicates (§2a).
- The union history shows in the files: implementations rows past the
  build-history bullets (outside any table), blank lines inside the vanilla
  table, ~360 vanilla rows appended under the "Upstreamable bugfixes" table.

Invariant this WP establishes: **no two PRs write the same ledger file.**

## Design

- `Ladruno_implementation/ledger/{implementations,quirks,vanilla}/WP-<nnn>-<slug>.md`.
  Migrated content is keyed `ADR-<n>`, `PR-<n>` or `LEGACY` (no id in the entry).
- Front matter: `wp`, `title`, `pr`, `date`, `status`, `files`, `class_tags`,
  `banner`, plus `section` (implementations: `table`/`history`), `table`
  (vanilla: `main`/`upstreamable`) and `legacy_seq` (migrated only). A strict
  YAML subset (JSON values), parsed without PyYAML.
- Body = today's row or entry, verbatim. A vanilla fragment is one PR's rows;
  the build groups rows by file, so two WPs touching one file are two
  fragments, never an edit.
- `<kind>/_template.md` holds each ledger's intro and conventions, with
  placeholder lines where the entries go.
- Follow-ups edit the OLDER fragment (`status: "fixed-by WP-<n>"`, a
  `#### Follow-up (WP-<n>)` section): a small-file 3-way merge, rarely concurrent.
- `ci/ledger.py`: `build` (newest WP first by `date`, then the migrated
  entries in their old order; vanilla grouped by file), `split` (the one-time
  migration), `migrate --from-diff` (an open branch's old-style edits ->
  fragments), `roundtrip` (the proof).
- Generated `LEDGER_*.md` go to `ledger/_build/` (gitignored); CI builds them
  and uploads the `ledgers` artifact. The committed `LEDGER_*.md` are ~15-line
  stubs, so `[[LEDGER_quirks]]` links still resolve; they record the split source.
- Gate `ci/check_ledger_fragments.py` (F1 schema, F2 name <-> wp, F3 unique
  name/body/heading, F4 class tags vs `SRC/classTags.h` via
  `check_classtags.parse`, F5 banner line exists, F6 templates + build, F7 the
  stub is untouched). Retargeted: quirk lint L3 (task-guide pointers) reads the
  quirk fragments; viewer gate V1 counts an added line under
  `ledger/implementations/`.
- `.gitattributes`: the three ledger `merge=union` lines are gone;
  `banner_features.txt` keeps union (append-only data, not generated).

## Where the design changed, and why

- **Generated files live in `ledger/_build/`, not at `LEDGER_*.md`.** One path
  cannot be both a gitignored build output and a committed stub.
- **`banner_features.txt` is not generated.** It is hand-written (the banner's
  source of truth, `patch_banner.py` reads it); the ledger only says each
  shipped row should have a line. A fragment may name its line in `banner`,
  and F5 checks it exists.
- **Migrated entries keep their old order** (`legacy_seq`); only new fragments
  sort newest-first. The old order is not chronological (union appends), but
  re-sorting would detach the `###` sub-entries that follow some `##` quirk
  entries from their parent.
- **Follow-up heading is `####`, not `##`.** A `##` inside a `###` quirk breaks
  the heading hierarchy of the generated ledger and would re-split as a new
  entry.
- **The "Upstreamable bugfixes" table keeps the rows its position gave it.**
  Most of its ~360 rows are ordinary vanilla rows appended at the end of the
  file; nothing in a row says which table it meant. A fragment's `table:` is
  now one field to fix, per PR.
- **classTag check:** `check_classtags.py` reads only the headers; the fragment
  gate imports its parser rather than changing it.

## Merge-time procedure (owner)

The split commit is LAST on the branch and must be regenerated against the
then-current `ladruno`, so PRs merged in between are captured:

```bash
git fetch origin
git switch wp/161-ledger-fragments
git log -1 --format=%s                    # must read "ledger(WP-161): split ..."
git reset --hard HEAD~1                   # drop that split commit (it is the tip)
git merge origin/ladruno                  # tooling commits only: no ledger conflict
python ci/ledger.py split                 # refuses if the round trip fails
python ci/check_ledger_fragments.py && python -m pytest -q ci/test_ledger.py
python ci/ledger.py build
git add -A Ladruno_implementation && git commit -m "ledger(WP-161): split the ledgers into fragments at ladruno <sha>"
git push --force-with-lease
```

Then every open PR that edited a `LEDGER_*.md` merges `ladruno` once more and
runs `ci/ledger.py migrate` (`ledger/README.md`, "A branch that still edits the
old LEDGER_*.md").

## Verification

See the PR body (CI runs, the merge-tree matrix, fragment counts).
