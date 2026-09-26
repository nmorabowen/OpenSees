# WP-120 code-health survey tooling

Read-only, dependency-free (Python 3.10+, git). Survey scripts, **not** CI gates. See
[`../120_code_health_survey.md`](../120_code_health_survey.md) for the findings.

| Script | What it does | Runtime |
|---|---|---|
| `common.py` | shared helpers; imports `clean()` / `functions()` from `ci/check_quirk_patterns.py` (one C++ scanner in the repo) | — |
| `inventory.py` | fork-added files per git vs the stamp vs `stamp_headers.py` GLOBS; writes `unstamped_fork_files.txt` | ~2 s |
| `clones.py` | CPD-style token-window clone finder (`--mode exact|norm`, `--scope fork|vanilla|cross`, `--min-tokens`, `--no-intra`, `--add-fork`, `--json`) | fork ~12 s; vanilla / cross ~10 min |
| `history.py` | clone families × first-parent PR history × `LEDGER_quirks` sections | ~7 s |
| `deadcode.py` | unbuilt files, never-referenced fork functions, dead guards, never-defined `#ifdef`s, commented-out statements | ~10 min (cleans all of SRC) |

```
python inventory.py
python clones.py --json results/exact.json
python clones.py --mode norm --json results/norm.json
python clones.py --scope cross --add-fork unstamped_fork_files.txt --json results/cross_plus.json
python history.py results/exact.json results/norm.json results/cross_plus.json --out results/history.json
python deadcode.py --json results/dead.json
```

Survey a historical tree (acceptance / trend runs): extract it with
`git archive <sha> SRC | tar -x -C <dir>` and set `WP120_ROOT=<dir>`.

`results/` holds the 2026-09-25 outputs on `ladruno` @ `bc5c33453` (the vanilla-scope JSONs, 0.6 and 1.5 MB,
are not committed; their one-line summaries are in `results/vanilla_summary.txt`).
