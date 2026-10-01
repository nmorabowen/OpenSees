"""O1 -- continuum rate oracle for LadrunoNORSAND (WP-144, equation sheet 144a sections 12-14).

Shared interface (identical names in O2):
  Params, Params.validate(), State, initial_state, run_path, tangent, acoustic_min_det,
  triaxial, k2_path.
"""
from .params import Params
from .integrator import State, integrate_increment
from .api import initial_state, run_path, tangent, triaxial
from .localization import acoustic_min_det, k2_path

__all__ = ["Params", "State", "initial_state", "run_path", "tangent", "acoustic_min_det",
           "triaxial", "k2_path", "integrate_increment"]
