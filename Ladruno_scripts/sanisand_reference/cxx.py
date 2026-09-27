"""Cross-check against the fork's C++ through WP-127's `ladrunoSANISANDReplay`.

The fork's opensees.pyd is a CPython-3.12 extension and the 3.12 test runner has
no scipy, so the C++ side runs in a SUBPROCESS (`_cxx_runner.py`, 3.12 with -S
and manual paths, asserting opensees.__file__) while the reference runs here.
The binary is used READ-ONLY.

Environment overrides:
  SANISAND_REF_OPENSEES_BIN   dist\\bin holding opensees.pyd with the replay command
  SANISAND_REF_PY312          the CPython 3.12 interpreter
  SANISAND_REF_SITE312        its site-packages (numpy)
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile

DEFAULT_BIN = (r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\.claude\worktrees"
               r"\tims-implementation-review-3733c6\dist\bin")
DEFAULT_PY = r"C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\python.exe"
DEFAULT_SITE = r"C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\Lib\site-packages"

# C++ prototypes: name -> (IntScheme, TolR, honorTolR, maxSubsteps)
PROTOS = {
    "ME8": (1, 1.0e-8, 1, 200000),    # ModifiedEuler, TolE = 1e-8 honoured: the TIGHT one
    "ME": (1, 1.0e-7, 0, 20000),      # the campaign setting (TolE hard-coded 1e-4)
    "RK45": (45, 1.0e-10, 1, 20000),  # RungeKutta45 (dT_min hard-coded 1e-3 -- WP-128)
}


def paths():
    return (os.environ.get("SANISAND_REF_OPENSEES_BIN", DEFAULT_BIN),
            os.environ.get("SANISAND_REF_PY312", DEFAULT_PY),
            os.environ.get("SANISAND_REF_SITE312", DEFAULT_SITE))


def available():
    b, py, _ = paths()
    return os.path.isfile(os.path.join(b, "opensees.pyd")) and os.path.isfile(py)


def run_jobs(jobs, params, protos=None, pmin=0.0101, presidual=0.0, timeout=3600):
    """jobs: list of dicts (proto, sigma, alpha, alpha_in, z, e, deps), all
    compression positive, Voigt order xx yy zz xy yz zx.  Returns the list of
    replay results (dicts) in the same order."""
    b, py, site = paths()
    runner = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_cxx_runner.py")
    spec = dict(params=list(params), protos=protos or PROTOS, pmin=pmin,
                presidual=presidual, jobs=jobs)
    with tempfile.TemporaryDirectory() as td:
        jp = os.path.join(td, "jobs.json")
        op = os.path.join(td, "out.json")
        json.dump(spec, open(jp, "w"))
        env = dict(os.environ)
        env.pop("PYTHONPATH", None)
        cp = subprocess.run([py, "-S", runner, b, site, jp, op], capture_output=True,
                            text=True, timeout=timeout, env=env)
        if cp.returncode != 0 or not os.path.isfile(op):
            raise RuntimeError("C++ replay subprocess failed:\n" + cp.stdout[-2000:]
                               + cp.stderr[-4000:])
        return json.load(open(op))
