#!/usr/bin/env python3
"""Run pytest against THIS worktree's dist/bin/opensees.pyd, never another build.

A boot `.pth` in the CPython 3.12 site-packages can import `opensees` from a different
worktree at interpreter start (LEDGER_quirks, "An INSTALLED Ladruno hijacks
`import opensees`"), so launch this with `python -S`: site processing is skipped,
the paths below are wired by hand, and the loaded module is asserted.

    <py3.12> -S Ladruno_implementation/wp123_undamped/run_pytest.py tests/test_x.py [pytest args]
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
DIST = os.path.join(ROOT, "dist", "bin")
TESTS = os.path.join(ROOT, "tests")
SITE = os.path.join(os.path.dirname(sys.executable), "Lib", "site-packages")

assert sys.flags.no_site, "run with python -S"
assert "opensees" not in sys.modules
os.add_dll_directory(DIST)
sys.path.insert(0, DIST)
sys.path.insert(0, TESTS)
sys.path.append(SITE)                      # pytest/numpy; appending does not run its .pth files
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")

import opensees                            # noqa: E402
assert os.path.normcase(opensees.__file__) == os.path.normcase(os.path.join(DIST, "opensees.pyd")), \
    opensees.__file__
print("opensees from", opensees.__file__, flush=True)

import pytest                              # noqa: E402

sys.exit(pytest.main(["-p", "no:cacheprovider", "-q", *sys.argv[1:]]))
