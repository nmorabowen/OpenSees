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

    pth = site_packages / "ladruno_opensees.pth"
    pth.write_text(bin_dir + "\n", encoding="ascii")
    print(str(pth))

    if sys.version_info[:2] != (3, 12):
        v = f"{sys.version_info.major}.{sys.version_info.minor}"
        print(f"WARNING: this venv is Python {v}; opensees.pyd is built for 3.12", file=sys.stderr)
        return 3

    return 0


if __name__ == "__main__":
    sys.exit(main())
