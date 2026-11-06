---
wp: LEGACY
title: "Build gotcha: cp (git-bash) / Copy-Item to the main checkout can silently skip the rebuild"
legacy_seq: 79
---
### Build gotcha: `cp` (git-bash) / `Copy-Item` to the main checkout can silently skip the rebuild
- **Bites:** editing fork source in the WORKTREE, copying to the main checkout, running `build.bat`, and testing — but the binary still shows OLD behavior (a new `-flag` is silently ignored, a new response returns empty). Two distinct traps stack: (1) `cp "src" "C:\…\dst"` from the Bash tool can mis-resolve the Windows backslash destination and write nowhere useful (the file you think you copied is unchanged in the main checkout — `grep -c <newsymbol>` there returns 0); (2) even after a correct `Copy-Item`, ninja compares mtimes and **`Copy-Item` PRESERVES the source's (older) mtime**, so if the worktree file was edited before the last build's `.obj`, ninja sees the object as newer and SKIPS recompiling — `opensees.pyd`'s timestamp never advances.
- **Tells:** the built `.pyd` LastWriteTime does not change after a "successful" build; `grep -c "<your new symbol>" <main-checkout-source>` returns 0; a parser silently ignores your new token (the RC parser treats unknown tokens as no-ops — forward-compat — so a non-compiled flag fails OPEN, not loud).
- **Fix (proven):** copy via `Copy-Item` (reliable on Windows), then explicitly bump the destination mtime `(Get-Item $dst).LastWriteTime = Get-Date` (and the `.cpp` that `#include`s a changed header) BEFORE `build.bat`, and confirm the build log shows `Building CXX object …<file>.cpp.obj` + `copying OpenSeesPy.dll -> opensees.pyd` and that the `.pyd` timestamp advanced. Run `build.bat` via the PowerShell tool, not the Bash heredoc (the latter captured only the banner here). Learned 2026-06-17.
