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
