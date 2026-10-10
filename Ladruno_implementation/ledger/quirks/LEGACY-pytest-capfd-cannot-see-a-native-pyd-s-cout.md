---
wp: LEGACY
title: "pytest capfd cannot see a native .pyd's cout/cerr on this Windows build"
legacy_seq: 393
---
### pytest `capfd` cannot see a native `.pyd`'s `cout`/`cerr` on this Windows build
- **Bites:** a test using `capfd` to assert on `opensees.pyd` stdout/stderr sees nothing, though the same code prints under a piped shell or `subprocess`. The `.pyd`'s own linked CRT writes through a stream the mid-process `dup2` swap does not reach.
- **Rule:** run the model in a child process and capture OS-level stdout/stderr; helper `_run_child()` in `tests/test_adr94_hlist_mechanical.py`. Pass `stdin=subprocess.DEVNULL` — without it `subprocess.run` under pytest intermittently raises `OSError: [WinError 6] The handle is invalid` from `_winapi.DuplicateHandle` (~1 in 3 runs observed).
