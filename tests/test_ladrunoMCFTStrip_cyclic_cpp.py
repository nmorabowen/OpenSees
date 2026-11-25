"""Ladruno WP-256: LadrunoMCFTStrip -cyclicDegradation in the g++ gate, without OpenSees.

Core and cover strips with the option on, on the WP-254 paths (tests/_testbed/mcft_strip_gpp.cpp, mode
"stripc"), against the numpy reference (tests/_testbed/mcft_strip_ref.py, c_cyc): bit-identical. The gate's
self-checks C1-C3 (a monotonic path unchanged, unchanged until the first re-opening, repeated cycles lose
strength) run in test_ladrunoMCFTStrip_cpp.py::test_self_checks. Skipped without a C++ compiler.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_ladrunoMCFTStrip_cpp as G  # the WP-254 g++ helpers (exe fixture, paths, strip parameters)

pytestmark = [pytest.mark.zone_a, pytest.mark.t1]

CCYC = 300.0
exe = G.exe          # the module-scoped g++ fixture


# ------------------------------------------------------------------------------------------ g++ paths
@pytest.mark.parametrize("m, drv", [(0.0, 0), (1.0, 0), (1.0, 1)], ids=["m0", "m1", "m1both"])
@pytest.mark.parametrize("which", ["core", "cover"])
@pytest.mark.parametrize("path", list(G._paths()))
def test_strip_path_with_the_option_vs_reference(exe, which, path, m, drv):
    cc = CCYC if m == 0.0 else CCYC * 100.0          # dP ~ re-opening x opening: a ~100x larger rate
    i = G.CORE if which == "core" else G.COVER
    rows = G._paths()[path]
    text = ("STRIP " + " ".join(G._f(v) for v in G._strip_params(i))
            + f" {G._f(cc)} {G._f(m)} {drv}\n{len(rows)}\n")
    text += "".join(f"{G._f(ex)} {G._f(g)} {c}\n" for ex, g, c in rows)
    got = G._run(exe, "stripc", text)
    st_ = G._ref_strip(i)
    st_.c_cyc = cc
    st_.m_cyc = m
    st_.d_cyc = drv
    st = st_.new_state()
    exp = []
    for ex, g, c in rows:
        o = st_.stresses(np.array([ex]), np.array([g]), st)
        exp.append([o["sx"][0], o["txy"][0], o["ey"][0]])
        if c:
            st_.commit(np.array([ex]), np.array([g]), st)
    exp = np.array(exp)
    ds = np.abs(got[:, :2] - exp[:, :2]).max() / G.SEC.conc.fc
    dy = np.abs(got[:, 2] - exp[:, 2]).max()
    assert ds < G.TOL_STRESS and dy < G.TOL_EY, f"{which} {path}: max |ds|/fc = {ds:.3e}, max |d eps_y| = {dy:.3e}"
    if path == "cyclic_shear":
        # the path re-opens cracks, so the option acts (P is a strain for m = 0, a strain^2 for m = 1)
        assert st.reopen.max() > (1e-3 if m == 0.0 else 1e-4)
