"""Ladruno test-bed pytest configuration.

Enforces the two-zone contract WITHOUT breaking upstream's unmarked tests:
  * a test may not carry both @zone_a and @zone_b,
  * @zone_b tests are auto-skipped when gmsh / apeGmsh are unavailable and run
    with APEGMSH_SKIP_VIEWER=1 so a headless box never tries to open a viewer,
  * unmarked tests (the existing upstream ones) are left alone.

SLOW TIER (`@pytest.mark.slow`).  Orthogonal to the zones: a case may be
`zone_a` (dep-free, upstreamable) and still cost minutes because the PHYSICS is
expensive — a collapse load pushed to a plateau, a convergence sequence over
three meshes.  Such cases must not be paid for on every push, but must also not
be exiled to zone_b, which would make them invisible to a box that simply lacks
gmsh.  So `slow` is opt-in and auto-skipped by default:

    pytest tests/                       # slow cases skipped
    pytest --runslow tests/             # slow cases run
    LADRUNO_RUN_SLOW=1 pytest tests/    # same, for CI/env-driven runners
    pytest --runslow -m slow tests/     # only the slow tier

A `slow` case MUST state its measured wall time in its module docstring and in
its ledger row — a runtime nobody has measured is a runtime nobody budgeted.

See Ladruno_implementation/testbed/00_canonical_testbed.md.
"""
import os
import importlib.util

import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--runslow", action="store_true", default=False,
        help="run @pytest.mark.slow cases (minutes-scale physics; see the "
             "conftest docstring). LADRUNO_RUN_SLOW=1 does the same.",
    )


def _has(mod: str) -> bool:
    return importlib.util.find_spec(mod) is not None


def pytest_collection_modifyitems(config, items):
    run_slow = config.getoption("--runslow") or os.environ.get("LADRUNO_RUN_SLOW")
    skip_slow = pytest.mark.skip(
        reason="slow tier — opt in with --runslow or LADRUNO_RUN_SLOW=1"
    )
    for item in items:
        if not run_slow and item.get_closest_marker("slow"):
            item.add_marker(skip_slow)
            continue
        zones = {m.name for m in item.iter_markers()} & {"zone_a", "zone_b"}
        if zones == {"zone_a", "zone_b"}:
            item.add_marker(
                pytest.mark.skip(reason="test marked both zone_a and zone_b — pick one")
            )
            continue
        # zone_b needs meshing (gmsh). A case that additionally needs apeGmsh /
        # LS-DYNA / etc. guards itself with pytest.importorskip — we only gate on
        # gmsh here so a gmsh-only case (e.g. the notched-bend crack-band study)
        # runs on a box that has gmsh but not the apeGmsh wrapper.
        if "zone_b" in zones and not _has("gmsh"):
            item.add_marker(
                pytest.mark.skip(reason="zone_b meshing dep (gmsh) not installed")
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


# ---- Windows native-PATH self-heal --------------------------------------
# gmsh.initialize() REPLACES the process's native Win32 PATH with a short
# stub (system dirs + sys.prefix) and never restores it. Python's os.environ
# is a startup snapshot, so python code never notices — but every child
# process inherits the nuked native block: bare-name CreateProcess lookups
# fail (WinError 2: "g++" not found) and freshly-built MinGW exes die at
# load with 0xC0000135 (libstdc++ not on the inherited PATH). Bit us via a
# module-level gmsh mesh in a zone_b file (import-time = pytest collection),
# breaking the g++ checker tests battery-wide; zone_b tests legitimately
# running gmsh mid-battery would do the same. Re-sync the native PATH from
# os.environ before every test. Cost: one GetEnvironmentVariableW per test.
if os.name == "nt":
    import ctypes as _ct

    @pytest.fixture(autouse=True)
    def _native_path_resync():
        want = os.environ.get("PATH", "")
        buf = _ct.create_unicode_buffer(65536)
        n = _ct.windll.kernel32.GetEnvironmentVariableW("PATH", buf, 65536)
        if (buf.value if n else "") != want:
            _ct.windll.kernel32.SetEnvironmentVariableW("PATH", want)
        yield
