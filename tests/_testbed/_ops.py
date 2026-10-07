"""Single source of the dual-import idiom upstream uses everywhere.

Local build first (CI copies opensees.so into tests/), pip openseespy otherwise.
"""
import os

try:  # local source build on CI / dev box
    import opensees as ops
except ModuleNotFoundError:  # installed wheel
    import openseespy.opensees as ops

# Ladruno (WP-176): a child started by `subprocess_run.pinned_child()` must run
# the parent's engine.  Fail at import, not as a silent cross-build comparison.
_expect = os.environ.get("LADRUNO_EXPECT_OPENSEES")
if _expect and os.path.normcase(os.path.abspath(getattr(ops, "__file__", ""))) \
        != os.path.normcase(_expect):
    raise ImportError("child imported opensees from %r, parent pinned %r"
                      % (getattr(ops, "__file__", None), _expect))

__all__ = ["ops"]
