"""pytest setup for the robust-solve torture battery (ADR-31 acceptance gate).

Puts Ladruno_scripts/ on sys.path (for ladruno_solve / robust_drive and the
sibling fixture modules) and silences the splash banner before opensees loads.
The fixtures resolve the built `opensees.pyd` from LADRUNO_DIST (defaulting to
the main worktree's dist/bin), so this gate runs against the live build.
"""
import os
import sys

os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")

_HERE = os.path.dirname(os.path.abspath(__file__))
_SCRIPTS = os.path.abspath(os.path.join(_HERE, ".."))
for _p in (_SCRIPTS, _HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)
