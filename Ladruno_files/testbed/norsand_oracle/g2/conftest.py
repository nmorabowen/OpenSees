"""G2 gate (WP-144, Zone B): pytest plumbing for Ladruno_files/testbed/norsand_oracle/g2/.

The fork's `--runslow` option and the zone / slow markers live in tests/conftest.py, which does not apply to this
folder; this file provides the same contract (opt-in slow tier: `--runslow` or LADRUNO_RUN_SLOW=1) so that

    python -m pytest <file> -q -rA --tb=short -p no:cacheprovider --runslow

works from anywhere.  It also registers the markers used here and keeps the G2 modules importable (g2_common.py
puts norsand_oracle/, kernel_parity/ and g2/ on sys.path).
"""
import os

import pytest

import g2_common  # noqa: F401  (side effect: sys.path, DLL directory, the opensees import)


def pytest_addoption(parser):
    try:
        parser.addoption("--runslow", action="store_true", default=False,
                         help="run @pytest.mark.slow cases (opt-in tier; LADRUNO_RUN_SLOW=1 does the same)")
    except ValueError:      # already registered (the whole norsand_oracle tree run with another conftest)
        pass


def pytest_configure(config):
    for m in ("zone_a: dep-free test (OpenSees + numpy only)",
              "zone_b: needs the oracles (scipy / sympy / fork-local Python); fork Zone B",
              "no_gmsh: zone_b case that never meshes: exempt from tests/conftest.py's gmsh-missing skip (the fork's "
              "zone_b gate is session-wide, so this folder needs the marker when collected from the repo root)",
              "slow: opt-in slow tier (--runslow); the docstring states the measured wall time"):
        config.addinivalue_line("markers", m)


def pytest_collection_modifyitems(config, items):
    run_slow = config.getoption("--runslow", default=False) or os.environ.get("LADRUNO_RUN_SLOW")
    skip_slow = pytest.mark.skip(reason="slow tier: opt in with --runslow or LADRUNO_RUN_SLOW=1")
    for item in items:
        if not run_slow and item.get_closest_marker("slow"):
            item.add_marker(skip_slow)
