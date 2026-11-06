---
wp: LEGACY
title: "A fully-prescribed model (zero free DOFs) under constraints Transformation FATALLY exits the process — FullGenLinSOE::getX - vectX == 0"
legacy_seq: 164
---
### A fully-prescribed model (zero free DOFs) under `constraints Transformation` FATALLY exits the process — `FullGenLinSOE::getX - vectX == 0`

- **Bites:** any prescribed-displacement rig that pins/`sp()`s EVERY DOF of a small
  patch model (single-element material-response probes are the classic case: all x
  fixed for uniaxial strain, bottom y fixed, top y driven by `ops.sp` in a pattern).
  Under `constraints Transformation` the sp-handled DOFs are condensed OUT, the
  equation count hits 0, and `FullGenLinSOE::getX()` hits a raw
  `opserr << "FATAL ..."; exit()` — killing the whole Python kernel/pytest run with
  no traceback (surfaced 2026-07-07 while building `tests/test_planestrain_sigma_zz.py`,
  PR #525). Other SOEs have sibling zero-size exits; this is not FullGeneral-specific.
- **Why:** Transformation removes constrained DOFs from the numbered system; a model
  where every DOF is fixed or sp-prescribed leaves size-0 vectors that the SOE treats
  as an allocation failure, and OpenSees's error path is `exit`, not a recoverable
  analysis error.
- **Workaround/status (2026-08-13): FIXED — the workarounds below are no longer
  required.** Root cause was never the solvers (every LAPACK solver already had its
  `if (n == 0) return 0;` quick return) but `setSize(Graph&)` in **six** SOEs:
  `FullGenLinSOE`, `BandGenLinSOE`, `BandSPDLinSOE`, `ProfileSPDLinSOE`,
  `SProfileSPDLinSOE` build their `vectX`/`vectB` (FullGen also `matA`) wrappers only
  under `if (size != oldSize)`, and `DiagonalSOE` excludes `size == 0` from that block
  *explicitly*. With zero equations `size == oldSize == 0`, so the wrappers keep the
  null value the default constructor gave them and the first `getX()`/`getB()` takes
  the `FATAL ... exit(-1)` branch. `exit(-1)` is a clean process exit, **not** a signal
  -- hence no traceback and nothing from `faulthandler` -- and the FATAL text goes to
  `opserr`, which the Python module redirects, so nothing prints at all. Fixed with six
  one-line `|| vectX == 0` guards, plus a `size > 0` guard on the two ProfileSPD
  variants' `profileSize = iDiagLoc[size-1]` (which read `iDiagLoc[-1]` when an existing
  SOE was resized *down* to zero). Provably inert for any model with free equations: the
  new branch needs `vectX == 0 && size == oldSize`, reachable only at the first
  `setSize` with `size == 0`. A zero-equation solve is trivially successful and now
  returns rc=0 with the `sp` values enforced, matching `UmfPack` (which always worked).
  Gate: `tests/test_soe_zero_free_equations.py` (zone_a, subprocess-isolated; 6 of its
  13 cases fail pre-fix). See the six `SRC/system_of_eqn/linearSOE/**` rows in
  [[LEDGER_vanilla_files]].
- **NOTE -- this quirk was rediscovered THREE times before anyone fixed it** (2026-06-03
  PR #155 while building LadrunoRCConcrete; 2026-07-07 PR #525 building
  `tests/test_planestrain_sigma_zz.py`; 2026-08-12 on an `sp`-driven MC/MCTC brick
  probe), each time landing a *workaround* in a different section of this file. If you
  find yourself writing a fourth, fix the code instead.
- **Superseded workaround** (kept for context): for fully-prescribed rigs use `constraints("Penalty", 1e15, 1e15)` -- the
  prescribed DOFs then STAY in the system (size > 0) and the penalty violation at
  1e15 vs typical stiffness is ~1e-10 relative, invisible to material-response
  checks. Alternatively leave at least one genuinely free DOF in the model.
