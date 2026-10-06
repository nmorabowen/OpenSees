---
wp: PR-735
title: "An INSTALLED Ladruno hijacks import opensees in every venv it has wired — sys.path.insert cannot win *(ROOT-CAUSED and FIXED 2026-08-11, #735 — see the last bu…"
date: 2026-08-11
legacy_seq: 250
---
### An INSTALLED Ladruno hijacks `import opensees` in every venv it has wired — `sys.path.insert` cannot win  *(ROOT-CAUSED and FIXED 2026-08-11, #735 — see the last bullet; the workarounds below are kept because they still apply to any venv wired by an older installer)*
- **Bites:** any script that bootstraps a *worktree* build with the standard
  `sys.path.insert(0, "<worktree>/dist/bin"); import opensees`. The Ladruno
  installer writes `ladruno_opensees.pth` into the venv's `site-packages`,
  which imports `_ladruno_opensees_boot` — and that module prepends
  `C:\Program Files\Ladruno\OpenSees\bin` to `sys.path` **and runs
  `import opensees`** (to alias `openseespy` onto it) at INTERPRETER STARTUP.
  By the time line 1 of your script executes, `sys.modules["opensees"]` is
  already the installed build; the `sys.path.insert` is a no-op.
- **Tell:** the feature you just built is "missing". Measured this session:
  `ERROR LadrunoUP 1 -- -geom 'corot' not supported: only 'linear' is accepted
  (ADR 71 §2.4)` on a worktree whose own `dist/bin` build accepts it fine. The
  banner is the giveaway — it prints the *installed* build's commit, not the
  worktree's.
- **Workaround (still applies for a venv you cannot touch):** run campaign/
  testbed scripts with an interpreter that has no Ladruno `.pth` (the base
  `C:\Users\<u>\AppData\Local\Programs\Python\Python312\python.exe`), and
  **assert which engine loaded** — compare `os.path.dirname(opensees.__file__)`
  against the intended `dist/bin` and `raise SystemExit` on mismatch.
  `Ladruno_files/testbed/hypo_bearing/bearing_backbone.py` carries that guard
  as the reference pattern. Note the venvs that DO have the `.pth` (e.g.
  `opensees_env`) are still the right ones for apeGmsh mesh work — just not
  for running a worktree build. *2026-07-28 (ADR-79 bearing campaign).*
- **Fix (2026-08-10):** `wire_venv_pth.py`'s generated boot module now checks
  `LADRUNO_OPENSEES_BIN` / `LADRUNO_OPENSEESMP_BIN` FIRST, before its baked-in
  install dirs — a runtime escape hatch that needs no re-run of the wirer and
  does not disturb any OTHER session sharing the venv:
  `set LADRUNO_OPENSEES_BIN=<worktree>\dist\bin` before `import opensees` (or
  `import openseespy.opensees`, since that alias chains through the same
  eager import) binds the worktree build instead of the install. Verified
  end-to-end in `opensees_env`: without the override, `ladrunoBuild()` reads
  the installed hash; with it set, the SAME process reads the worktree's.
  Regenerate an already-wired venv's boot script with the fixed
  `wire_venv_pth.py <bin-dir> [<mp-dir>]` to pick up the fix without a full
  installer re-run. Gate: `tests/test_wire_venv_pth_override.py` (no built
  engine needed — renders `BOOT_TEMPLATE` and asserts which dir wins).
  apeGmsh's live-backend resolver (`opensees/emitter/live.py`) was
  independently hardened the same day: it no longer trusts a pre-bound
  `opensees` module just because it is *some* fork build (`criticalTimeStep`
  present) — it now also checks the module's `__file__` directory matches
  `APEGMSH_OPENSEES_BIN` before reusing it, closing the gap for code paths
  that reach `sys.modules['opensees']` before apeGmsh's own resolver runs.

- **ROOT-CAUSE FIX (2026-08-11, #735): the boot module no longer imports anything at startup.** The
  2026-08-10 entry above treats `LADRUNO_OPENSEES_BIN` as *the* escape hatch, which conceded the premise
  — that `import opensees` must happen at interpreter startup. It did not. The eager import existed only
  to alias `openseespy`/`openseespy.opensees` onto the sequential build; that alias is now resolved by a
  lazy `sys.meta_path` finder which imports the engine on the FIRST request for the name and not before.
  The rest of the boot module (`sys.path`, `add_dll_directory`, process-local `PATH`) only REGISTERS
  search locations — it loads nothing — so the module is now passive in the sense BUILD_GOTCHAS §5 asks
  for. **Both symptoms go at once:** `sys.path.insert(0, <worktree>/dist/bin)` wins again unaided
  (verified: the same venv that used to force the install now resolves to the worktree), and a bare venv
  interpreter stops pinning the install's DLLs, which is what made installer UPGRADES fail with
  `DeleteFile failed; code 5`. `LADRUNO_OPENSEES_BIN` survives as a deliberate override, demoted from
  crutch. **Re-run `wire_venv_pth.py <bin-dir> [<mp-dir>]` once per venv** — an already-wired venv keeps
  the old eager boot script until you do. Gated by `tests/test_wire_venv_pth_override.py` (5 tests): the
  two override tests above, plus an AST check that no `import opensees` sits outside a function
  (verified non-vacuous against three regression shapes — module-level, inside an `if`, and `from`-import),
  a behavioural check that nothing is aliased at exec but both names resolve afterwards, and one that the
  alias is skipped under `PMI_RANK`. **Measuring this needs care: two obvious probes lie.**
  `Get-Process($pid).Modules` reported 7 modules and 0 held for a Python whose own stdout proved it had
  imported the installed `.pyd`; and `tasklist /m X /fi "PID eq N"` misses because the venv launcher
  re-spawns under a different PID than `Start-Process` returns. Use unfiltered `tasklist /m opensees.pyd`
  as a set difference against a live baseline, and always run the eager-import CONTROL — a probe that
  reports "no holders" for both cases is measuring nothing.
- **Update (F2-c):** the `openseespy` alias is now OPT-IN. Default wiring no longer aliases it; wire with `wire_venv_pth.py --alias-openseespy` or set `LADRUNO_OPENSEESPY_ALIAS=1` at run time. Re-wire already-wired venvs only if you want the new default (an old boot module keeps the alias).
