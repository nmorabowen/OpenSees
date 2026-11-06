"""WP-176 -- fresh-interpreter children run the PARENT's engine, or fail loudly.

WP-175 found that the ADR-97 gate-4 children (`[sys.executable, dump_hist.py]`,
no `-S`) resolved `opensees` through whatever `site` did first: a boot `.pth`
in site-packages (`wire_venv_pth`'s `_ladruno_opensees_boot`) put a worktree's
`dist\\bin` at `sys.path[0]` and imported it before `PYTHONPATH` counted. A
byte-identity gate could then compare a different build against the baseline
and say nothing.

`_testbed.subprocess_run.pinned_child()` starts the child with `-S`, passes the
parent's `sys.path` as `PYTHONPATH`, and pins the parent's `opensees.__file__`
in `LADRUNO_EXPECT_OPENSEES`. `_testbed._ops` raises `ImportError` on any other
engine. These tests pin all three halves, including the hijack itself: a decoy
`opensees.py` placed ahead of the real engine must be REFUSED, not used.
"""
import os
import subprocess
import sys

import pytest

from _testbed import ops
from _testbed.subprocess_run import EXPECT_ENV, pinned_child

pytestmark = [
    pytest.mark.zone_a,
    # the pip-openseespy fallback has no module named `opensees` to pin
    pytest.mark.skipif("opensees" not in sys.modules,
                       reason="no local `opensees` engine loaded"),
]

HERE = os.path.dirname(os.path.abspath(__file__))

PROBE = ("import sys; from _testbed import ops; "
         "print('NO_SITE', sys.flags.no_site); print('ENGINE', ops.__file__)")


def _run(argv, env, cwd=HERE):
    return subprocess.run(argv, cwd=cwd, env=env, capture_output=True,
                          text=True, stdin=subprocess.DEVNULL,
                          encoding="utf-8", errors="replace", timeout=120)


def test_pinned_child_runs_the_parents_engine_without_site():
    argv, env = pinned_child(["-c", PROBE])
    assert "-S" in argv
    assert env[EXPECT_ENV] == os.path.abspath(ops.__file__)
    p = _run(argv, env)
    assert p.returncode == 0, p.stdout + p.stderr
    assert "NO_SITE 1" in p.stdout
    got = [ln.split(" ", 1)[1] for ln in p.stdout.splitlines()
           if ln.startswith("ENGINE ")]
    assert got and os.path.normcase(os.path.abspath(got[0])) \
        == os.path.normcase(os.path.abspath(ops.__file__)), p.stdout


def test_a_decoy_engine_ahead_on_the_path_is_refused(tmp_path):
    """The WP-175 hijack in miniature: something ELSE named `opensees` wins the
    import. The child must die with the pin message, not run the decoy."""
    (tmp_path / "opensees.py").write_text("ladrunoBuild = lambda: 'decoy'\n")
    argv, env = pinned_child(["-c", PROBE])
    env["PYTHONPATH"] = str(tmp_path) + os.pathsep + env["PYTHONPATH"]
    # cwd too: `-c` puts the cwd at sys.path[0], and CI copies the real
    # opensees.so into tests/ -- run from tests/ and the decoy never loads.
    p = _run(argv, env, cwd=str(tmp_path))
    out = p.stdout + p.stderr
    assert p.returncode != 0, out
    assert "ENGINE" not in p.stdout, out
    assert "parent pinned" in out, out


def test_no_pin_means_no_guard():
    """Without the env var `_testbed._ops` is unchanged (every non-pinned
    caller, and the parent itself)."""
    argv, env = pinned_child(["-c", PROBE])
    env.pop(EXPECT_ENV)
    p = _run(argv, env)
    assert p.returncode == 0, p.stdout + p.stderr
