---
wp: LEGACY
title: "gmsh.initialize() REPLACES the process's native Win32 PATH — later child processes fail DLL/exe resolution battery-wide"
legacy_seq: 160
---
### `gmsh.initialize()` REPLACES the process's native Win32 PATH — later child processes fail DLL/exe resolution battery-wide

- **Bites:** any pytest battery (or long-lived python session) that runs gmsh
  and LATER spawns a child by bare name or runs a freshly-built MinGW exe.
  Symptoms observed (2026-07-07): `test_hrz_standalone_kernel` exe dies at load
  with `0xC0000135` (STATUS_DLL_NOT_FOUND — libstdc++ not on the inherited
  PATH) and `test_cpp_kernel_matches_oracle_dump` raises
  `FileNotFoundError: WinError 2` on bare `"g++"`. Both PASS in isolation and
  fail in the full battery — the classic order-dependence smell.
- **Why (proven with a ctypes probe):** `gmsh.initialize()` calls
  `SetEnvironmentVariable("PATH", <system dirs + C:\ + sys.prefix>)` — a ~315
  char stub replacing the real ~2 kB PATH — and `finalize()` does NOT restore
  it. Python's `os.environ` is a STARTUP SNAPSHOT, so python-side code (and
  `shutil.which`) still sees the good PATH; only native consumers break:
  `CreateProcess` bare-name lookup and the DLL loader of child processes both
  read the LIVE native block. Import alone is harmless — the nuke fires at
  `initialize()`. The worst variant was a zone_b module running a module-level
  gmsh mesh: pytest COLLECTION imports every module (deselection happens
  after), so the corruption hit even `-m zone_a` batteries that never run a
  single gmsh test.
- **Fixes shipped:** (1) never call gmsh at module level in a test file (the
  zoneb mesh is now lazy via `_mesh_once()`); (2) `tests/conftest.py` autouse
  `_native_path_resync` re-syncs the native PATH from `os.environ` before
  every test on Windows; (3) the two g++ checker tests are hermetic anyway —
  absolute `shutil.which("g++")` path, `env=os.environ.copy()` on every
  `subprocess.run`, exe outputs under `tmp_path`.
- **General rule:** in any long-lived process that touched gmsh, pass
  `env=os.environ.copy()` and absolute executable paths to `subprocess`.
