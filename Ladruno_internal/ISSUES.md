# GitHub issues are the shared memory (WP-174)

**Board:** #937 (pinned; `gh issue list --label board`).

Agents on this fork run in parallel sessions, on different machines, and lose
their context at every compaction. Whatever one session knows and the next one
needs must be somewhere every session reads. That place is GitHub: issues,
draft-PR comments and labels. The model is the same as apeGmsh's
(`internal_docs/program/PROGRAM.md` there).

This file is the **static** part: where things live, the labels, the templates.
It changes by PR. **Live state never goes in this file, in the repo, or in a
session's private memory.**

## 1. What lives where

| Kind of knowledge | Lives in | Lifetime |
|---|---|---|
| A live WP: its state, gates, next action | its **draft PR** (body + handoff comments) | until merge |
| A multi-WP program (e.g. TIMs F18–F23) | one **`program` issue**; each WP links it | until the program closes |
| What is open across the fork right now | the **board** (pinned `board` issue) | permanent; rows refreshed |
| A bug, gap, or finding not being fixed yet | one **issue** per finding (§4 template) | until fixed or `wontfix` |
| Long-running jobs (esmeralda SLURM, campaigns) | a handoff comment on the owning PR/issue: job ids, roots, logs, pass gates | until the jobs are read |
| A decision only the owner can make | `human-gate` label + one comment stating the question | until answered |
| A failure class seen twice | a **`lesson`** issue → a lint rule or quirk entry within 7 days | until the rule lands |
| What shipped, what diverges from upstream, a gotcha learned | the **ledger fragments** (AGENTS.md) | permanent, reviewed |
| User preferences, pointers ("program X is #n") | private memory | permanent |

Two consequences:

- **Private memory holds pointers, not state.** "TIMs program: #n" is fine;
  open-PR lists, job ids, SHAs or "resume here" notes are not. If a memory file
  starts accumulating state, move the state to the issue and leave the pointer.
- **Issues are mutable, ledgers are reviewed.** Something learned in an issue
  that future agents must not rediscover graduates to a `quirks/` fragment or a
  `ci/check_quirk_patterns.py` rule, in a PR. Close the issue with a link to it.

## 2. Session protocol

1. **Boot.** Read the board, then the issue/PR you were pointed at with
   `gh issue view <n> --comments` (or `gh pr view <n> --comments`). The last
   handoff comment is where you resume. Check `in-flight` before starting work
   someone else may be on.
2. **Claim.** Add `in-flight` to the issue/PR you are working. Remove it when
   you stop, even if unfinished.
3. **Findings as you go.** A bug or gap outside your WP's scope is a new issue
   (§4), not a paragraph in your report and not a memory file.
4. **Hand off.** Before the session ends (or nears its usage limit), post the
   §3 handoff comment and refresh the board row. A session that ends without a
   handoff has lost its work for everyone else.

## 3. Handoff template

Post as a comment on the WP's draft PR (or the `program` issue):

```markdown
### Handoff — WP-<nnn> — <YYYY-MM-DD> — <model> @ <effort>
- **Head:** <sha> — <one line on what it contains>
- **Verified:** <commands run + result; Windows/Linux logs cited>
- **Running:** <host, job ids, roots, log paths — or "nothing">
- **Blocked / human gates:** <what, waiting on whom>
- **Decisions needed from the owner:** <question, options>
- **Lessons** (second occurrence → open a `lesson` issue): <...>
- **Next:** <the exact next action, runnable without this session's context>
```

## 4. Finding template

The fork's existing bug issues (#904–#913) set the shape; the issue form in
`.github/ISSUE_TEMPLATE/finding.md` carries it:
**Symptom · Where · Repro · Expected · Evidence · Workaround in use.**
Cite code as `path::symbol` where you can; line numbers drift.

## 5. Labels

| Label | Meaning |
|---|---|
| `board` | The one pinned board issue. |
| `program` | A multi-WP program issue. |
| `in-flight` | A session is working it now. |
| `blocked` | Waiting on a gate, a run, or another WP. |
| `human-gate` | Needs the owner: decide, merge, configure. |
| `lesson` | A failure class seen twice; must become a rule within 7 days. |

Plus GitHub's defaults (`bug`, `enhancement`, `wontfix`, ...). Create a
missing label with `gh label create <name> --force`.

## 6. Board format

The board issue body is one table, one row per live WP or program:

```markdown
| WP / program | PR / issue | State | Next action | Gate |
|---|---|---|---|---|
| WP-173 unknown-system gate | #929 | draft, test-only | wait for WP-172 | WP-172 merge |
```

Edit the row (`gh issue edit <board> --body-file -`) when you post a handoff;
delete it when the WP merges or is dropped.
