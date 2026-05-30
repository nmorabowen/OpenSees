"""Ladruno test-bed pytest configuration.

Enforces the two-zone contract WITHOUT breaking upstream's unmarked tests:
  * a test may not carry both @zone_a and @zone_b,
  * @zone_b tests are auto-skipped when gmsh / apeGmsh are unavailable and run
    with APEGMSH_SKIP_VIEWER=1 so a headless box never tries to open a viewer,
  * unmarked tests (the existing upstream ones) are left alone.

See Ladruno_implementation/testbed/00_canonical_testbed.md.
"""
import os
import importlib.util

import pytest


def _has(mod: str) -> bool:
    return importlib.util.find_spec(mod) is not None


def pytest_collection_modifyitems(config, items):
    for item in items:
        zones = {m.name for m in item.iter_markers()} & {"zone_a", "zone_b"}
        if zones == {"zone_a", "zone_b"}:
            item.add_marker(
                pytest.mark.skip(reason="test marked both zone_a and zone_b — pick one")
            )
            continue
        if "zone_b" in zones and not (_has("gmsh") and _has("apeGmsh")):
            item.add_marker(
                pytest.mark.skip(reason="zone_b deps (gmsh/apeGmsh) not installed")
            )


@pytest.fixture(autouse=True)
def _zone_b_headless(request):
    """zone_b cases must never pop a viewer (kernel-crash on headless/CI)."""
    if request.node.get_closest_marker("zone_b"):
        prev = os.environ.get("APEGMSH_SKIP_VIEWER")
        os.environ["APEGMSH_SKIP_VIEWER"] = "1"
        yield
        if prev is None:
            os.environ.pop("APEGMSH_SKIP_VIEWER", None)
        else:
            os.environ["APEGMSH_SKIP_VIEWER"] = prev
    else:
        yield
