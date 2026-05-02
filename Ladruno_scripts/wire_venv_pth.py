"""Drop a .pth file into a venv's site-packages so `import opensees` works.

Invoked by installer.iss at install time:
    <venv>\\Scripts\\python.exe wire_venv_pth.py <bin-dir-to-add-to-syspath>

Kept tiny on purpose - the Inno Pascal Script can't conveniently do
sysconfig lookups, so it shells out to the venv's own python and asks it.
"""
import sys
import sysconfig
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: wire_venv_pth.py <bin-dir>", file=sys.stderr)
        return 2

    bin_dir = sys.argv[1]
    site_packages = Path(sysconfig.get_paths()["purelib"])
    if not site_packages.is_dir():
        print(f"site-packages not found: {site_packages}", file=sys.stderr)
        return 1

    # Two-line .pth:
    #   line 1: path entry — adds bin_dir to sys.path so `import opensees` finds opensees.pyd
    #   line 2: executable hook — site.py runs any .pth line starting with `import` at
    #           interpreter startup. We register our `opensees` module under the
    #           `openseespy` and `openseespy.opensees` aliases so packages that hardcode
    #           `import openseespy.opensees` (aprGmsh, openseespy-style helpers) get our
    #           build (with MPCO/HDF5) instead of vanilla pip openseespy.
    #
    # The `opensees.opensees = opensees` self-attribute is required because
    # `import openseespy.opensees as X` compiles to IMPORT_NAME + IMPORT_FROM,
    # and IMPORT_FROM calls getattr(sys.modules['openseespy'], 'opensees')
    # rather than reading sys.modules['openseespy.opensees']. Without the
    # attribute, the binding step raises ImportError.
    pth = site_packages / "ladruno_opensees.pth"
    alias_hook = (
        "import sys, opensees; "
        "sys.modules['openseespy'] = opensees; "
        "sys.modules['openseespy.opensees'] = opensees; "
        "opensees.opensees = opensees"
    )
    pth.write_text(bin_dir + "\n" + alias_hook + "\n", encoding="ascii")
    print(str(pth))

    if sys.version_info[:2] != (3, 12):
        v = f"{sys.version_info.major}.{sys.version_info.minor}"
        print(f"WARNING: this venv is Python {v}; opensees.pyd is built for 3.12", file=sys.stderr)
        return 3

    return 0


if __name__ == "__main__":
    sys.exit(main())
