"""ADR-97 P6 -- Zone-A smoke test for ``Ladruno_implementation/adr97_oracle/
measure_p6.py``, the Cerro-Lindo-scale measurement driver behind the D1
default-flip decision (``97_ladruno_asdp_closest_point_adr.md`` D1, Phases
table P6 row).

This does NOT run the real measurement (that is a 2-scale x 2-family x
5-configuration sweep, each its own subprocess, minutes per Cerro-Lindo-scale
run -- see ``Ladruno_implementation/reviews/adr97_p6_measurement.md``).  It
runs the driver's tiny ``smoke`` mesh (a 3x3x2-element box, 144 DOF) for 3
push steps, one FRESH ``python3.12`` child process per configuration --
required, not a style choice: ``ASDPlasticMaterial3D``'s per-tag
``integration_method``/``tangent_type`` option maps are process-global
statics keyed by material tag (ADR-94 lesson, reused by ADR-97 P1/P2), so
reusing one process across configurations would silently read back a
previous configuration's option map.

Asserts:
  * every one of the four configurations (``Backward_Euler``/Secant,
    ``Backward_Euler``/Continuum, ``Closest_Point``/Continuum,
    ``Closest_Point``/Algorithmic) converges (gravity + all 3 push steps,
    rc == 0 throughout);
  * ``Closest_Point`` + ``Algorithmic``'s total Newton-iteration count is
    <= the shipped default's (``Backward_Euler`` + ``Secant``) -- the
    ORDERING ADR-97 P1/P2 measured on smaller rigs (7.1x and comparable
    fewer iterations), pinned here at mesh scale without pinning exact
    counts (those are measured, not test-worthy, in the real P6 run).

Not itself an input to D1 -- see ``reviews/adr97_p6_measurement.md`` for that.
"""
import json
import os
import subprocess
import sys
import tempfile

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
_DIST_DIR = os.path.dirname(os.path.abspath(ops.__file__))
_DRIVER = os.path.join(_TESTS_DIR, os.pardir, "Ladruno_implementation",
                        "adr97_oracle", "measure_p6.py")

CONFIGS = ["BE_Secant", "BE_Continuum", "CP_Continuum", "CP_Algorithmic"]


def _run_config(config, out_dir):
    out_path = os.path.join(out_dir, "smoke_%s.json" % config)
    env = dict(os.environ)
    env["PYTHONPATH"] = _DIST_DIR + os.pathsep + env.get("PYTHONPATH", "")
    env["PATH"] = _DIST_DIR + os.pathsep + env.get("PATH", "")
    env["LADRUNO_OPENSEES_QUIET"] = "1"
    cmd = [sys.executable, _DRIVER,
           "--family", "MC", "--scale", "smoke", "--config", config,
           "--push-steps", "3", "--grav-steps", "2", "--out", out_path]
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env,
                           timeout=60, cwd=_TESTS_DIR, stdin=subprocess.DEVNULL)
    assert proc.returncode == 0, (
        "measure_p6.py child for %s exited %d\nstdout tail:\n%s\nstderr tail:\n%s"
        % (config, proc.returncode, proc.stdout[-2000:], proc.stderr[-2000:]))
    with open(out_path) as fh:
        return json.load(fh)


def test_smoke_all_configs_converge_and_cp_algorithmic_is_cheaper():
    with tempfile.TemporaryDirectory() as out_dir:
        results = {c: _run_config(c, out_dir) for c in CONFIGS}

    for config, r in results.items():
        assert r.get("error") is None, "%s: %r" % (config, r.get("error"))
        assert r["grav"]["first_fail_step"] is None, (
            "%s: gravity stage failed to converge at step %s"
            % (config, r["grav"]["first_fail_step"]))
        assert r["push"] is not None and r["push"]["first_fail_step"] is None, (
            "%s: push stage failed to converge at step %s"
            % (config, r["push"]["first_fail_step"] if r["push"] else "?"))
        assert len(r["push"]["iters"]) == 3, (
            "%s: expected 3 converged push steps, got %d"
            % (config, len(r["push"]["iters"])))

    be_secant_total = results["BE_Secant"]["total_iters"]
    cp_algo_total = results["CP_Algorithmic"]["total_iters"]
    assert cp_algo_total <= be_secant_total, (
        "Closest_Point+Algorithmic (%d total Newton iterations) should not cost "
        "more than the shipped default Backward_Euler+Secant (%d) on this rig "
        "-- ADR-97 P1/P2 measured the opposite ordering at smaller scale"
        % (cp_algo_total, be_secant_total))
