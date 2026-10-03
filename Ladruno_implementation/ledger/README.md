---
title: Ledger fragments — how to write a ledger entry
project: Ladruno
tags:
  - ledger
---

# Ledger fragments (WP-161)

The three build-control ledgers are **one file per entry**, so no two PRs ever
write the same file. `LEDGER_implementations.md`, `LEDGER_quirks.md` and
`LEDGER_vanilla_files.md` one folder up are stubs; the full ledgers are
generated.

```
ledger/
  implementations/   one feature row (or one build-history bullet) per file
  quirks/            one quirk per file, `### ` heading first
  vanilla/           one PR's vanilla-file rows per file, one row per file touched
  <kind>/_template.md   the intro + conventions; edit only to change those
  _build/            generated LEDGER_*.md (gitignored)
```

- **Read:** `rg <pattern> Ladruno_implementation/ledger/`, or
  `python ci/ledger.py build` and open `ledger/_build/LEDGER_*.md`.
- **Write:** add `ledger/<kind>/WP-<nnn>-<slug>.md` in the same PR as the
  change. Never edit `LEDGER_*.md` (the gate fails on an edited stub).
- **Follow-up to an older entry** (a fix, a status flip, a retraction): edit
  THAT entry's fragment — set `status: "fixed-by WP-<nnn>"` and add a
  `#### Follow-up (WP-<nnn>)` section (quirks: `####`, so it stays inside the
  `###` entry) or update the row in place (tables).
  Never copy the entry into a new fragment; the gate rejects duplicate bodies.
- **Check:** `python ci/check_ledger_fragments.py` (a static gate in `ladruno.yml`).

## Format

YAML front matter, then the body exactly as the old ledger held it. One
`key: value` per line; quote strings as JSON (`"..."`), lists as JSON arrays.

```markdown
---
wp: WP-161
title: "Per-WP ledger fragments"
pr: "#915"
date: 2026-10-02
status: "draft"
class_tags: ["ELE_TAG_LadrunoBrick=33002"]
banner: "the exact line in Ladruno_scripts/banner_features.txt"
---
| **WP-161 — per-WP ledger fragments** | tooling | — | `ci/ledger.py` | draft | #915 |
```

| key | needed | meaning |
|---|---|---|
| `wp` | always | `WP-<n>`; migrated content uses `ADR-<n>`, `PR-<n>` or `LEGACY`. Must match the file name prefix. |
| `title` | always | short plain-text title |
| `date` | new fragments | `YYYY-MM-DD`; the build puts the newest first |
| `status` | implementations rows | as in the row's Status cell |
| `pr` | when known | `"#915"`; fill it after the PR opens |
| `files` | vanilla | the file cell of each row, in order |
| `class_tags` | if the entry owns tags | `SYMBOL` or `SYMBOL=value`, checked against `SRC/classTags.h` (say RESERVED/PLANNED for a tag not yet in the header) |
| `banner` | optional | the matching `banner_features.txt` line (checked) |
| `section` | implementations | `table` (default) or `history` |
| `table` | vanilla | `main` (default) or `upstreamable` |
| `legacy_seq` | migrated only | position in the old single-file ledger; do not add it to new fragments |

Bodies: **implementations** — one `| Feature | Kind | Class tag | Files |
Status | PR(s) |` row. **quirks** — a `### <symptom you would search for>`
heading and the entry (no second `##`/`###` inside; use `####`).
**vanilla** — `| file | why | PR |` rows, one per vanilla file this PR touched;
two WPs touching one file are two fragments, the build groups the rows by file.

## A branch that still edits the old LEDGER_*.md

Merging `ladruno` into a branch cut before WP-161 conflicts on the three
ledgers. Turn the branch's ledger edits into fragments, then take the stubs:

```bash
git fetch origin
git merge origin/ladruno                          # conflicts in LEDGER_*.md
python ci/ledger.py migrate --from-diff "$(git merge-base HEAD MERGE_HEAD)" --wp WP-<nnn>
git checkout MERGE_HEAD -- Ladruno_implementation/LEDGER_implementations.md \
    Ladruno_implementation/LEDGER_quirks.md Ladruno_implementation/LEDGER_vanilla_files.md
python ci/check_ledger_fragments.py && git add -A Ladruno_implementation && git commit
```

`migrate` writes added entries as new fragments (`--wp`, today's date) and
applies in-place edits of old entries to their fragments; it prints `MANUAL`
for what it cannot place (a deleted entry, a changed intro).
