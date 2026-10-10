---
wp: WP-107
title: "A whole-file CRLF→LF rewrite of a vanilla file is invisible in review and permanent in git blame (WP-107, red-team S1)"
legacy_seq: 458
---
### A whole-file CRLF→LF rewrite of a vanilla file is invisible in review and permanent in `git blame` (WP-107, red-team S1)

An editor helpfully normalised `SRC/material/nD/UWmaterials/ManzariDafalias.{cpp,h}`
while making a 78-line change. Both files are pinned `-text` in `.gitattributes`, so
git stores the bytes verbatim and nothing normalised them back: the PR diff read
**5635/5584 and 431/410**, ~11,000 of its 11,474 additions were line-ending noise,
the review surface was inflated ~50x, `git blame` pointed the entire file at the WP,
and a concurrent PR touching the same file conflicted **wholesale** (`git merge-tree`
confirmed both before and after the fix).

**Check before every PR on this fork:** `git diff origin/ladruno..HEAD --numstat` —
any file whose additions ≈ deletions ≈ its own line count is a conversion, not a
change. `git diff -w --ignore-cr-at-eol` shows what really changed. The fix is to
rewrite the file with its original endings and re-commit; the content is unaffected.
