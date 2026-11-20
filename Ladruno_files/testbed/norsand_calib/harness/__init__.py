"""WP-144 P3 calibration harness for LadrunoNORSAND (Zone B: numpy/scipy, runs on Esmeralda).

Modules:
  sand       material constants of the sand that are NOT fitted (CSL, M, elastic targets) + Toyoura placeholders
  energy     the ENERGY plug: maps the DM04 elastic targets onto the oracle's energy parameters (BA06, HAR) and onto deck flags
  model      ParamVector -> O2/O1 Params (one place where the plug, the sand and the fitted vector meet)
  drivers    element drivers on O2 (drained PS / TC / TE, undrained PS / TC) and their O1 counterparts
  data       CSV schema loader/writer ("eps_a_pct, sr, eps_v_pct" + metadata) and synthetic curves
  objective  residuals (stress ratio, eps_v, peak phi, eps_peak) with the weights stated in one dict
  fit        constrained reparametrisation, least_squares multi-start, identifiability (Jacobian, profiles)
Scripts:
  validate_drivers.py   O2 drivers -> O1 at first order (one path per driver; --energy BA06|HAR)
  check_plugs.py        seconds-long wiring checks (energy plugs, unified pi_i0, global policy)
  smoke_recovery.py     the harness gate: recover known parameters from an O2-generated synthetic set (clean, and noisy at another step)
"""
