"""WP-134 -- reference vs the fork's C++ (ladrunoSANISANDReplay, WP-127).

SKIPS unless the WP-127 opensees.pyd and the CPython 3.12 runner are present
(Ladruno_scripts/sanisand_reference/cxx.py; override with
SANISAND_REF_OPENSEES_BIN / SANISAND_REF_PY312 / SANISAND_REF_SITE312).
Gate: on monotonic plastic increments from benign states (class A of the doc,
section 5) the ModifiedEuler-faithful continuum `Options.uw_me()` equals the C++
ModifiedEuler at TolR 1e-8 (-honorTolR 1) to 1e-5 relative in the stress
increment; the paper model differs by the named UW additions (U4, U9)."""
import os
import sys

import numpy as np
import pytest

pytest.importorskip("scipy")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                                "Ladruno_scripts"))

from sanisand_reference import CAMPAIGN, Control, Options, integrate  # noqa: E402
from sanisand_reference import cxx  # noqa: E402
from sanisand_reference.crosscheck import benign_states, job_of  # noqa: E402
from sanisand_reference.model import t2v, v2t  # noqa: E402

pytestmark = pytest.mark.skipif(not cxx.available(),
                                reason="WP-127 opensees.pyd / CPython 3.12 runner not found")

CASES = [("p20_TC", [1e-4, 0, 0, 0, 0, 0]), ("p100_TC", [1e-3, 0, 0, 0, 0, 0]),
         ("p50_TE", [-1e-4, 0, 0, 0, 0, 0]), ("p50_TCshear", [0, 0, 0, 1e-3, 0, 0])]


def test_uw_me_continuum_equals_cpp_modified_euler_on_monotonic_increments():
    states = dict(benign_states())
    jobs = [job_of(states[s], de, "ME8") for s, de in CASES]
    out = cxx.run_jobs(jobs, CAMPAIGN.as_opensees())
    for (s, de), c in zip(CASES, out):
        st = states[s]
        r = integrate(st, Control.strain(de), CAMPAIGN, Options.uw_me(), record=False)
        assert r.status == "ok" and c["rc"] == 0
        assert [m["mode"] for m in r.segments][0] == "plastic"   # class A
        d_ref = t2v(r.state.sigma) - t2v(st.sigma)
        d_cpp = np.array(c["sigma"]) - t2v(st.sigma)
        assert np.linalg.norm(d_ref - d_cpp) < 1e-5 * np.linalg.norm(d_cpp), s
        a_ref = t2v(r.state.alpha) - t2v(st.alpha)
        a_cpp = np.array(c["alpha"]) - t2v(st.alpha)
        assert np.linalg.norm(a_ref - a_cpp) < 1e-4 * np.linalg.norm(a_cpp), s
        # the PAPER model is measurably different (U4 G(e_init), U9 frozen K,G)
        rp = integrate(st, Control.strain(de), CAMPAIGN, Options(), record=False)
        d_p = t2v(rp.state.sigma) - t2v(st.sigma)
        assert np.linalg.norm(d_p - d_cpp) > 1e-3 * np.linalg.norm(d_cpp), s
