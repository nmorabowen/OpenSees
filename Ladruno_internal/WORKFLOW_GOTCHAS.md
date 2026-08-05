---
title: Ladruno PR / CI / git workflow gotchas
project: Ladruno
tags:
  - internal
  - workflow
  - ci
---

# Workflow gotchas — PR, CI, and git on the `ladruno` fork

Hard-won process lessons for working on this fork. None of these are about the
*code* — they are about the **PR/CI/git workflow**, which behaves differently
here than on a normal repo because `ladruno` **auto-merges open PRs very fast**
(an automated babysitter merges `guppi/*` worktree branches within minutes,
usually squash, several per hour). That single fact causes most of the traps
below.

> These were previously only in an agent's machine-local memory. Ported here so
> the knowledge survives a fresh clone on another machine. See also the global
> `~/.claude/CLAUDE.md` stacked-PR `--base` note and the project `CLAUDE.md`
> ledger rules.

---

## 1. Stranded commits — don't pile increments onto a merged PR branch

**The trap:** you open a PR from a `guppi/*` branch; it auto-merges (often
squash) within minutes. You then push *more* commits to the **same branch**
assuming they'll append to the open PR. They don't — a merged PR is terminal,
pushes neither reopen it nor create a new PR, so those commits sit **stranded**
on an orphaned branch, absent from `ladruno`, while `ladruno` advances without
them.

This has bitten repeatedly (PR #65 LadrunoBrick, #82 LadrunoJ2, #151 EAS, #163
InitDefGrad) — even follow-ups pushed *minutes* after opening race the
auto-merge.

**The tells:**
- `gh pr view <n> --json state,commits` shows `state: MERGED` with `nCommits: 1`
  (squash) while you have more local commits.
- `git branch -r --contains <your-code-sha> | grep ladruno` prints **nothing**.

**Rules:**
- **One logical PR per branch.** For *any* follow-up increment, branch fresh
  from current `ladruno` — never reuse a branch whose PR you already opened.
- Re-checking `state==OPEN` before pushing is **not enough** — it can merge
  between the check and the push. Just always branch fresh.
- After the final push, **prove the code landed**:
  `git branch -r --contains <code-sha> | grep ladruno` must print a line.

**Recovery:** open a NEW PR from the leaf branch → `ladruno`. First
`git merge origin/ladruno` locally to pull current and resolve conflicts in the
shared vanilla files (`SRC/classTags.h`, `FEM_ObjectBrokerAllClasses.cpp`,
`OpenSeesElementCommands.cpp`, the ledgers — the material team edits these too),
rebuild, re-run the battery, then push and PR. Stranded commits cherry-pick
cleanly onto the fresh branch in practice.

**Two fresh-branch hazards that bit during recovery (ladruno advancing mid-session):**
1. `git stash -u` + `stash pop` can sweep in *stale* files from a parallel
   session and silently revert ladruno's newer copies. **Always**
   `git diff --cached origin/ladruno --stat` before committing and confirm the
   staged set is EXACTLY your files.
2. If your branch base is an *ancestor* of current `ladruno`, the diff attributes
   others' work to your PR. Rebase onto `origin/ladruno` tip:
   `git diff HEAD -- <my src> > p.patch; git reset --hard origin/ladruno; git apply p.patch`,
   then redo doc/ledger edits.
3. `git reset --hard` **deletes staged-but-uncommitted new files** (untracked
   survive, but `git add`ed-then-uncommitted do not). Copy new files aside or
   `git stash` them first. (Recovered a lost doc once via
   `git fsck --lost-found` → `git cat-file -p <blob>`.)

---

## 2. "Ledger CI issue" = a stale, conflicting PR (CI never ran)

Because `ladruno` merges several PRs per hour, a feature branch goes stale within
hours and its PR re-conflicts on the shared ledgers.

**Why it looks like a CI problem:** an un-mergeable PR (content conflict vs
current `ladruno`) makes GitHub **refuse to run** the `pull_request` workflow on
the merge ref. `gh pr checks` then shows **"no checks reported"** — which looks
like CI is just slow, but actually means the PR is CONFLICTING and CI is blocked.
That presents to a human as a vague "ledger CI issue."

**Recovery recipe (proven on PR #108):**
1. `git fetch origin ladruno`; check distance with
   `git rev-list --count HEAD..origin/ladruno`.
2. `git merge --no-commit --no-ff origin/ladruno` and resolve:
   - **`LEDGER_implementations.md`** almost always conflicts (every sibling PR
     edits it). Cleanest: `git checkout --theirs` the ledger, then **re-apply
     only YOUR row deltas** with targeted edits — don't hand-merge the giant rows.
   - **`FEM_ObjectBrokerAllClasses.cpp`** conflicts when a sibling added a broker
     include+case next to yours → **keep-both** (additive).
   - Other ledgers usually auto-merge.
3. **Re-run the gates on the MERGED tree** (`python ci/check_classtags.py` +
   `ci/check_manifest.py`) — CI runs them on the merge, so this is the real check.
4. **Inherited manifest debt** (see §3): backfill any missing rows the gate
   surfaces, not just your own.
5. Commit the merge, push. Verify `gh pr view <n> --json mergeable` → MERGEABLE
   and CI starts. (`mergeable=UNKNOWN` right after a push is just GitHub
   recomputing; `statusCheckRollup` conclusions are authoritative.)
6. Churn can re-bite between fetch and push — just merge again; the second
   conflict is usually tiny.

NB cross-namespace classTags do NOT collide in `check_classtags` (family = text
before the final `_<Name>`), so e.g. `ELE_TAG_BezierTri6`=33000 and
`MAT_TAG_LadrunoUniaxialJ2`=33000 coexist fine.

---

## 3. CI gates are fast-check-only — `ladruno` HEAD can be broken

**Auto-merge gates ONLY on the classTag+manifest fast check — NOT the Zone-A
(Ubuntu) job (neither the build nor the pytest).** Consequences:

- A PR that **fails to compile** can still merge (PR #175 did). A broken build on
  `ladruno` HEAD then makes EVERY later PR's Zone-A red — and since the build
  dies, the pytest phase never runs, so test bugs stay masked until the build is
  fixed.
- `check_manifest.py` (gate **G9**) is **NON-BLOCKING** — PRs merge red — so a
  classTag shipped without its `manifest.yaml` row leaves `ladruno` HEAD already
  failing G9, and every later PR *inherits* the failure.

**Rules:**
- After pushing C++ to a fork PR, **watch the Zone-A job** and fix-forward. Green
  fast-gate ≠ green build. A fast (~1–2 min) Zone-A fail = compile error; a slow
  (~5–6 min) fail = test failure. If `ladruno` HEAD is broken, your fix unblocks everyone.
- **`gh pr checks <n> --watch` exits PREMATURELY — do not trust its exit 0 as "Zone-A
  passed".** It watches only the check-runs registered *at the moment it polls*. Zone-A
  registers a few seconds after push, so if classTag (fast) has already passed and Zone-B
  is `skipping`, the watch sees "all known checks done" and exits 0 **before Zone-A even
  appears** (observed on #384/#385: watch exited 0, then `gh pr checks` showed Zone-A
  *pending*). For a code PR, confirm the ACTUAL Zone-A run instead:
  `gh pr checks <n>` must list Zone-A with a terminal `pass`, OR watch the run id directly
  — `gh run watch <run-id> && gh run view <run-id> --json conclusion`. Merging on a
  premature watch-exit risks landing a Linux-only break on `ladruno` HEAD (the C4/D2 merges
  got lucky — Zone-A passed after the fact, but verify, don't assume).
- **Put the `classTags.h` `#define`, the `manifest.yaml` row, AND the ledger row
  in the SAME commit.** Never a follow-up commit — a failing gate is exactly what
  tempts you to "just push the fix" into the auto-merge race (→ §1 stranding).
- When recovering, run `python ci/check_manifest.py` on the fresh branch FIRST —
  it lists ALL missing rows (yours + inherited); backfill them all in one PR.

---

## 4. `gh pr` multi-line body — use `--body-file -`, never `--body @-`

To give `gh pr create`/`gh pr edit` a multi-line body from a heredoc, use
**`--body-file -`** (reads stdin). **NEVER `--body @-`** — that is curl syntax;
`gh` takes `@-` as the LITERAL body string and ignores stdin, so the PR ships
with the 2-char body `@-` and your whole description is silently dropped (no
error). This bit across PRs #175–#181.

```sh
gh pr create --base ladruno --title "..." --body-file - <<'EOF'
...body...
EOF
```
(closing `EOF` at column 0). Sanity-check once:
`gh pr view <n> --json body --jq '.body|length'` — a length of ~2 means it failed.
`gh pr edit <n> --body-file -` backfills an already-created/merged PR.

---

## 5. Mandatory: stamp the LADRUNO header on every new authored file

Every **new fork-authored** class/source file (any element / material /
integrator / recorder / utility we author — the `LEDGER_implementations.md` set)
MUST carry the stamped LADRUNO header (ASCII art + the four-author credit
"Nicolas Mora Bowen · Patricio Palacios · José Abell · Guppi"). This is
**non-optional**, part of "done" for the implementation, not a later cleanup
(standing user instruction).

**How:** add the new file's path/glob to the `GLOBS` list in
`Ladruno_scripts/stamp_headers.py`, then run `python Ladruno_scripts/stamp_headers.py`
(idempotent; inserts the block below the OpenSees/PEER header). `--check` exits 1
if any authored file is unstamped — good as a pre-commit / CI gate. Do **not**
stamp vanilla upstream files (→ §6); those keep their original header + inline
`// Ladruno` markers.

---

## 6. Vanilla-file footprint — keep it minimal and additive

Different change budgets by surface, all in service of long-term upstreamability:

- **Frozen `MPCORecorder.cpp` (MPCO):** only touch what pertains to *our own
  elements*. Deliberately frozen (byte-identical to upstream). Do not add features
  or refactor it.
- **Ladruno recorder (`MPCORecorderLadruno.cpp` + `MPCOL_*`):** ours — change freely.
- **Vanilla OpenSees files** (element/material/transf/etc.): keep edits to the
  **minimum needed**, strictly additive, never alter existing behavior. Mark every
  touch with a `// Ladruno` comment and add a `LEDGER_vanilla_files.md` row (so the
  set is reconstructable via `grep -rn "Ladruno" SRC/`).
- **EXCEPTION — genuine upstream bugs:** if something in vanilla OpenSees is an
  actual error (not just inconvenient), surface it, **ask, and fix only with the
  user's approval** — then still mark `// Ladruno` + ledger it. Real bugfixes are
  the most upstreamable kind of change; the minimal-footprint rule is about not
  bolting fork-only *features* into vanilla files, not leaving known bugs in place.

**Litmus test before editing a vanilla `SRC/` file:** "is this the smallest
additive hook that makes our feature work?" If a change belongs in our own code,
put it there, not in the vanilla file.

### 6a. `.gitattributes` marks ~5,300 files `-text`, so line-ending churn inflates the recorded footprint ~170x

**Root cause, and it is not anybody's editor being wrong.** The root
`.gitattributes` opens with `* text=auto !eol` and then spends **5,323 of its
5,345 lines** listing individual files as `-text` — a CVS/SVN→git conversion
artifact inherited from upstream. `-text` tells git to store bytes **verbatim**
and normalize nothing, so those files have **no canonical line ending**: whatever
the last committer's editor emitted is what git keeps. Any contributor whose
editor disagrees rewrites the whole file, and the §6 footprint audit then reads a
routine registration hook as a core rewrite.

Measured case: PR #700 records **10,543** changed lines across
`SRC/interpreter/OpenSeesCommands.{cpp,h}`; the real change is **63 lines**. The
parent blob was ordinary CRLF (5,240 LF / 5,240 CR), the child is pure LF
(5,303 LF / 0 CR) — a plain CRLF→LF flip, stored verbatim because `-text` was in
force on both files.

- **Detect:** `diff --strip-trailing-cr <(git show REV^:path) <(git show REV:path)`.
  Tiny result + huge `git show --numstat` = ending churn. To see it directly,
  `git cat-file blob $(git rev-parse REV:path) | tr -cd '\r' | wc -c` — and use
  `git ls-files --eol` for the index-vs-worktree view across many files.
  Beware `grep -c $'\r$'` inside nested command substitution: if the `$'...'`
  quoting doesn't survive, the pattern degenerates and matches *every* line, which
  reads exactly like "all lines are CRLF". If a CR count equals the total line
  count on both sides of a comparison, suspect the measurement, not the file.
- **Audit with it:** any footprint review, upstream-port sizing, or
  `LEDGER_vanilla_files` row for a big-numstat interpreter file must quote the
  strip-trailing-cr count, not the raw numstat, or the fork looks far more
  divergent from upstream than it is.
- **Fixed for the registration surface (PR #704).** The eight files every fork
  feature touches to register a class — `OpenSeesCommands.{cpp,h}`,
  `{Python,Tcl}Wrapper.cpp`, `OpenSees{Element,Pattern}Commands.cpp`,
  `tcl/commands.cpp`, `classTags.h` — now carry an explicit `text` override at the
  END of `.gitattributes` (last matching line wins). All eight were already LF in
  the blob, so it was a **zero-diff** change that simply makes it permanent.
  Verified A/B: CRLF-ify a covered file and `git add` it → **empty** staged diff;
  do the same to a still-`-text` file (`OpenSeesCommandsTcl.cpp`) → **188/188**
  whole-file churn.
- **Still exposed:** every other `-text` file. If you find yourself fighting
  ending churn on one, add it to that block rather than re-normalizing by hand —
  but check `git ls-files --eol` first: adding `text` to a file whose blob is
  CRLF *does* produce a one-time whole-file diff, and for a vanilla file that
  becomes permanent divergence from upstream. Zero-diff only holds for blobs that
  are already `i/lf`. A repo-wide `git add --renormalize .` would rewrite ~447
  blobs and is **not** a casual cleanup on a fork that tracks upstream.

## `gh pr merge --auto` on this repo merges IMMEDIATELY (auto-merge is disabled), and ledger-placeholder fills race the squash

Repo settings have GitHub auto-merge DISABLED (`enablePullRequestAutoMerge`
errors). When you run `gh pr merge --auto --squash` on a PR whose base checks
are already satisfiable, gh falls back to a DIRECT merge — the PR lands
instantly, **without waiting for Zone-A**. Two consequences (both bit us on
2026-07-22, PRs #595/#598):

1. **Do not rely on `--auto` as a CI gate.** If Zone-A must pass first, watch
   the checks yourself and merge after.
2. **The ledger-placeholder race:** the "open PR → get number → fill the
   `| ADR-xx PR |` placeholder → push" flow loses when the PR merges on step 1.
   The fill commit lands on a branch whose base predates the squash ⇒ its own
   PR conflicts. Recovery: fresh branch off updated `ladruno`, redo the
   one-line fill (sed), separate docs PR (#600 pattern). Better: write the row
   with the placeholder, merge, then fill in the NEXT unit's PR.

Also: GNU `sed -i` on a CRLF `SRC/` file rewrites EVERY line ending — a
one-word edit became a 9,900-line diff (OpenSeesCommands.cpp). Use the
byte-preserving Edit path for CRLF sources; check `git diff --stat` before
committing after any sed.

**This bites from BOTH directions, and the warning above did not stop it.**
Hit again 2026-07-26 (ADR-77 C0-5), inverted: a Python patch script doing
`Path(p).write_text("\n".join(lines))` on **LF** sources rewrote every line as
**CRLF**, because `write_text` defaults to `newline=None` -> `os.linesep`. Five
inserted profiler scopes became **1119 changed lines** in `BFGS.cpp` (~1000 more
in `Broyden.cpp`/`KrylovNewton.cpp`).

Why LF sources are the norm here at all: **`.gitattributes` marks most vanilla
`SRC/` files `-text`** (253 KB of such entries), so git does *no* EOL
normalisation in either direction -- whatever byte the tool writes is a real
content change. **The repo's own config, not git's defaults, decides what counts
as a change.**

The cost is not cosmetic on this fork: it defeats `LEDGER_vanilla_files.md`'s
entire purpose (keep the divergence from upstream small and *greppable*) and
makes the PR unreviewable.

- **Avoid:** patch existing sources through bytes, not text --
  `p.write_bytes(p.read_bytes().decode().replace(...).encode())`, or pass
  `newline=""` to `open`. The Edit tool is already byte-preserving; prefer it.
- **Detect:** `git show --stat HEAD` immediately after committing. If a file's
  changed-line count is near its total length, this is why. *Do this even when
  the edit "obviously" touched five lines* -- that assumption is exactly what
  suppressed the check both times.
- **Repair:** `p.write_bytes(p.read_bytes().replace(b"\r\n", b"\n"))` (or the
  reverse) then `git commit --amend`, which is safe while the PR is unmerged --
  confirm with `gh pr view <n> --json state,mergedAt` before force-pushing.

**Related trap, same family (hit 2026-07-27, ADR-77 G2 ext):** the Bash tool's
heredoc layer EATS ONE LEVEL OF BACKSLASH even with a quoted delimiter, so an
inline python script whose anchor contains the two characters backslash-n
receives a real newline instead -- the anchor silently never matches (or worse,
matches wrongly). Write patch scripts to a file (Write tool is verbatim) and
construct backslashes with `chr(92)`.

## 7. The G1 Tcl gate has been green on a DEAD suite — "OK (0 checks passed)"

Found 2026-07-26 by the first marker-gated Tcl step (the ADR-76 LAPACK
regression, #643). The Zone-A `build/Release/OpenSees` binary aborts Tcl
initialization on the CI runner — `application-specific initialization failed:
Can't find a usable init.tcl` — because `TCL_LIBRARY` is not set and the conan
Tcl runtime is not on the binary's search list. After that failure NO OpenSees
command is registered (`invalid command name "wipe"`), so every deck dies on
its first command.

The "Tcl verification suite + ladruno tcl (gated, G1)" step has been hitting
this on EVERY run — visible as two `init.tcl` complaints in even fully green
logs (e.g. run 30173461993) — and still passing, because `runVerificationSuite`
leaves an empty `results.out` and `ci/check_tcl_results.py` reports
`OK (0 checks passed)`. A gate that counts failures must also refuse an empty
result set; "no failures because nothing ran" is the silent-truncation trap.

- **Honest pattern** (the LAPACK step): `export TCL_LIBRARY="$(dirname "$(find
  ~/.conan2 -name init.tcl -path '*/tcl8.6/*' | head -1)")"` before invoking the
  binary, and gate on a positive terminal marker, never the exit code.
- **RESURRECTED (#653, 2026-07-26, stacked on #643):** the G1 step got the same
  `TCL_LIBRARY` export and `check_tcl_results.py` now fails on `total == 0`
  (either fix alone is insufficient — without the export the suite dies empty;
  without the zero-check the death is invisible). First genuine execution:
  **`OK (19 checks passed)`** — all nine verification decks + the ladruno
  cantilever clean, including under the ADR-76 LAPACK singular fix. The gate is
  live; a future `0 checks ran` failure means the Tcl runtime went missing
  again, not that the decks regressed.
