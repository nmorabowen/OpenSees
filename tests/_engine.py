"""Worktree-engine binder shared by the zone_b batteries.

Implements the boot-.pth stale-pyd eviction WITHOUT the same-DLL re-import
that crashes the interpreter at shutdown (LEDGER_quirks entry
"sys.modules.pop(\"opensees\") + re-import ... crashes Python AT SHUTDOWN"):
re-importing the SAME opensees.pyd re-runs PyInit_opensees (moduledef.m_size
> 0), which re-registers Py_AtExit(cleanupFunc); the doubled cleanup then
wipes a dangling PythonModule* at Py_FinalizeEx (0xC0000005/0xC0000409,
heap-state dependent — small batteries pass by allocation luck). So: reuse
the cached module when it is already this worktree's pyd; only
evict-and-reimport a genuinely foreign one (a DIFFERENT DLL — one atexit
each, harmless).

Deliberately a TOP-LEVEL module (not part of the _testbed package):
importing _testbed runs its __init__, which imports opensees eagerly —
that would defeat the path setup this helper performs first, and would
break the overlay batteries' pytest-free standalone (`python -S`) mode.
This module imports nothing beyond os/sys.
"""
import os
import sys


def bind_worktree_engine(dist):
    """Bind and return THIS worktree's `opensees` module.

    `dist` is the worktree's dist/bin directory holding opensees.pyd. The
    caller is responsible for skipping (or exiting) when the engine is not
    built — by the time this runs, the pyd must exist.
    """
    dist = str(dist)
    os.environ["PATH"] = dist + os.pathsep + os.environ.get("PATH", "")
    try:
        os.add_dll_directory(dist)
    except (FileNotFoundError, OSError, AttributeError):
        pass
    if dist not in sys.path:
        sys.path.insert(0, dist)

    cached = sys.modules.get("opensees")
    if cached is not None and os.path.normcase(
            os.path.dirname(getattr(cached, "__file__", "") or "")
    ) == os.path.normcase(dist):
        return cached  # already this worktree's pyd — NEVER re-import it

    for m in ("opensees", "openseespy", "openseespy.opensees"):
        sys.modules.pop(m, None)
    import opensees as ops

    assert os.path.normcase(os.path.dirname(ops.__file__)) == os.path.normcase(dist), (
        f"wrong opensees.pyd imported: {ops.__file__} (want {dist})"
    )
    return ops
