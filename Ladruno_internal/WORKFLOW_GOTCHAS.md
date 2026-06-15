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
- After pushing C++ to a fork PR, **watch the Zone-A job**
  (`gh pr checks <n> --watch`) and fix-forward. Green fast-gate ≠ green build.
  A fast (~1–2 min) Zone-A fail = compile error; a slow (~5–6 min) fail = test
  failure. If `ladruno` HEAD is broken, your fix unblocks everyone.
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
