"""Single source of the dual-import idiom upstream uses everywhere.

Local build first (CI copies opensees.so into tests/), pip openseespy otherwise.
"""
try:  # local source build on CI / dev box
    import opensees as ops
except ModuleNotFoundError:  # installed wheel
    import openseespy.opensees as ops

__all__ = ["ops"]
