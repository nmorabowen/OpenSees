---
wp: LEGACY
title: "ParallelPlain numbers vertex 0 LAST (the tag-0 quirk) and ParallelNumberer::numberDOF(ID&) never worked at all"
legacy_seq: 205
---
### `ParallelPlain` numbers vertex 0 LAST (the tag-0 quirk) and `ParallelNumberer::numberDOF(ID&)` never worked at all
- **Bites:** (a) stock `ParallelPlain` checks "already ordered" against a zero-filled ID, so merged tag 0 always reads as present and gets pushed to the end — a valid but surprising permutation (bandwidth outlier on the first node). (b) The `numberDOF(ID& lastDOFs)` variant (SP/DomainDecomposition lane) has its numbering call commented out upstream and a mismatched recv layout — any caller gets garbage start-DOFs.
- **Workaround/status:** `LadrunoParallelPlain` fixes (a) (G1b-gated: valid bijection, differs from stock by design); the Ladruno override hard-errors on (b). *2026-07-22.*
