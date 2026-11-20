# WP-239 session handoff: LadrunoMembrane holds its frames by tag

Written in the ISSUES.md §3 shape so it can be pasted, unchanged, as the first
handoff comment on the WP-239 draft PR once the owner allows a push. No PR
exists yet; nothing is pushed.

### Handoff — WP-239 — 2026-11-20 — Claude Opus 5.5 @ default

- **Head:** branch `wp/239-membrane-frame-tags` (cut from `ladruno` cfe478254).
  - `0fbf08ffa`: the fix, the new tests and the ledger fragments.
  - `d335e4fa3`: the restore test compares committed displacements raw against raw, bit for bit.
  - This handoff commit.
- **The bug:** `LadrunoMembrane` resolved its supporting frames once, in `update()` from `Domain::addElement`, and cached
  a raw `Element*`. `remove element` deletes the beam, so a beam removed and re-added with the same tag left the
  membrane writing its slab share through a freed pointer. The NLRHA hinge step does exactly that when it moves each
  beam onto hinge-twin nodes. Linux (glibc) happened to give the freed block back to the same-tag beam; on Windows the
  same model crashes every time. A model restored from a database re-resolved against the twins instead. The twins are
  not membrane nodes, so each frame lost its link line and printed one "lies on the panel but is not one of its nodes"
  warning per end.
- **The fix** (`SRC/element/ladrunoMembrane/LadrunoMembrane.{h,cpp}`; fork-authored, so no vanilla row):
  1. `CachedFrameInfo` keeps `frameTag` plus the end node tags and coordinates; there is no `Element*`. `addLoad` calls
     `Domain::getElement(frameTag)` on every call. A missing frame warns once ("supporting frame <t> is no longer in the
     domain; its slab share is not applied") and is skipped. This was chosen over a Domain generation counter: it needs
     no upstream Domain edit, costs one map lookup per emitted segment, and keeps the geometry resolved at the start.
  2. A frame end that is not a membrane node but lies within 1e-6·L of one (a hinge twin) ties to that membrane node,
     the original. It keeps its WP-208 link line, its WP-209 end loads and its overhang-torque recipient, with no
     warning. Frame ends that are on the panel but coincide with no membrane node still warn.
  3. There is no new state, so `sendSelf`/`recvSelf` are unchanged; the frame tags were already sent.
- **Ledger:** `ledger/implementations/WP-239-membrane-frame-tags.md` (status draft) and
  `ledger/quirks/WP-239-element-pointer-after-remove.md`. `python ci/check_ledger_fragments.py`: 0 errors.
- **Verified** on a Windows full `Ladruno_scripts\build.bat` of this branch, under `py -3.12`:
  - **New tests:** `tests/test_wp239_membrane_frame_replace.py` passes 8/8. The model is
    `tests/_testbed/membrane_twin_model.py`: a 2×1 bay floor on columns, two membranes and seven beams, re-created on
    twins after the membranes exist.
    - The `equalDOF` and rigid `zeroLength` twins, with and without decoy allocations between the remove and the re-add,
      match the no-replay model in beam end shears and slab displacements to 1e-6.
    - The total vertical reaction equals the applied load, 495.0.
    - Membranes declared after the replay give the same answer with no warning.
    - A removed frame warns once.
    - Save in one process and restore in a fresh one, with 5 Plain patterns and 21 `-beamUniform` loads: the committed
      displacements are bit-identical, forces match to 1e-8, and there is no "not one of its nodes" or EleLoad error.
  - **Existing membrane suites** (`test_ladruno_membrane`, `_kernel_cpp`, `test_mutation_membrane`, `test_wp204`,
    `test_wp208`, `test_wp222`, `test_membrane_load_oracle` and `test_partial_trapezoidal_beam_loads`): all pass.
  - **Pre-fix binary** (main checkout `dist`, Nov 15):
    - Every replay case dies with "Windows fatal exception: access violation".
    - The restore case sends the slab load of a removed beam to `ZeroLength 5020` ("load type unknown"), the object
      allocated next into the freed block.
  - **Not run:** ASan and valgrind. There is no ASan build on Windows, and nothing was built on Esmeralda.
- **"LoadPattern::recvSelf - channel failed to recv the EleLoad ID":** not reproduced with a single writer.
  - A probe with 13 292 beams and 4 Plain patterns (53 168 `-beamUniform` loads) saved and restored in a fresh process
    cleanly, with identical reactions and displacements.
  - FileDatastore names its files `<db>.IDs.<size>.<tag>`, so two jobs on the same database path overwrite each other's
    records; that is the likely cause.
  - Workaround: give each job its own database path. No fix was made.
- **Running:** nothing.
- **YOO pin:**
  - YOO runs keep the current pin (`6050d78e1`) until a WP-239 pin is built and accepted. Its live results are not
    affected, because it resolved before the replay against the original nodes. Do not restore it from a database.
  - **Pin branch** `pin/yoo-6050d78e1-wp239` (local only): `6050d78e1` plus the same two commits (`bcece7a2a`,
    `5e56f06e1`). The membrane code is unchanged between `6050d78e1` and `ladruno`, so the cherry-pick is clean.
  - **Build kit:** `build_yoo_6050_wp239.sh` and `wp239-membrane-frame-tags.patch`, in the side-session notes folder
    for this session group. The script applies the patch with `git am` on a depth-1 clone of `6050d78e1`, which is
    already on GitHub, so no push is needed. It checks tree `02918011cfac7a64abf131c1970ae27fc4f97355` and writes to a
    new pin, `~/pins/yoo-ladruno-6050d78e1-wp239-py311-debian13`. The old pin and `env.sh` are untouched.
  - **Acceptance:** the G1 preload sums its reactions to 1 411 401.97 kN, per-beam end shears equal the old pin's, and a
    save/restore smoke test prints no "not one of its nodes".
- **Blocked / human gates:** the owner decides on the push, the draft PR, and the Esmeralda pin build. The session
  that did this work was told not to do any of the three without the owner's OK.
- **Decisions needed from the owner:**
  - Push `wp/239-membrane-frame-tags` and open the draft PR (base `ladruno`)?
  - Build the WP-239 pin on Esmeralda and switch YOO runs to it?
- **Lessons:** an element must never cache another element's `Element*` across commands. `remove element` frees it,
  and the allocator decides whether the bug shows. This is recorded as quirk WP-239; a lint rule is a candidate if it
  recurs.
- **Next:**
  1. With the owner's OK: `git push -u origin wp/239-membrane-frame-tags`.
  2. `gh pr create --draft --base ladruno`, and paste this file as the first handoff comment.
  3. Add a board (#937) row.
  4. Separately, the pin build above.
  5. Clean-up: the leftover `.claude/worktrees/wp239-pin` directory, if it is still on disk (`git worktree prune`).
