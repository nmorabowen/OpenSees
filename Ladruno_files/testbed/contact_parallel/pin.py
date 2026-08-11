"""Verify the OpenSees Python module came from the intended build.

The venv's `ladruno_opensees.pth` eagerly imports `opensees` at interpreter
startup from a baked-in directory -- which on this machine points at a copy left
in ANOTHER session's scratchpad. sys.path juggling cannot undo that (the module
is already in sys.modules), so the caller must set LADRUNO_OPENSEES_BIN /
LADRUNO_OPENSEESMP_BIN before Python starts. This just fails loud if they didn't.
"""
import os

LADRUNO = r"C:\Program Files\Ladruno\OpenSees"
DIRS = {"opensees": LADRUNO + r"\bin", "openseesmp": LADRUNO + r"\openseesmp"}


def load(name):
    mod = __import__(name)
    want, got = os.path.normcase(DIRS[name]), os.path.normcase(mod.__file__)
    if not got.startswith(want):
        raise RuntimeError(
            f"{name} resolved to {mod.__file__}\n"
            f"  expected under {DIRS[name]}\n"
            f"  set LADRUNO_OPENSEES_BIN / LADRUNO_OPENSEESMP_BIN before launching."
        )
    return mod
