"""Ladruno (WP-177, #960): the key a gate-4 baseline is filed under.

The same commit built with the same toolchain versions gives different bits on
another machine (WP-175), so a byte-identity baseline belongs to a HOST.
`hosts.json` maps `host_key()` to that host's baseline file; the test compares
`==` only against this host's own file.

The key is the platform plus the CPU brand string. A toolchain upgrade on the
same host does not change the key; if it moves the bits, the `==` leg fails
and the host's baseline is re-dumped deliberately (README "Adding a host").
"""
import platform
import re
import sys


def cpu_brand():
    if sys.platform == "win32":
        try:
            import winreg
            k = winreg.OpenKey(
                winreg.HKEY_LOCAL_MACHINE,
                r"HARDWARE\DESCRIPTION\System\CentralProcessor\0")
            return winreg.QueryValueEx(k, "ProcessorNameString")[0]
        except OSError:                                 # pragma: no cover
            pass
    return platform.processor() or platform.machine() or "unknown"


def host_key():
    return "%s | %s" % (sys.platform, re.sub(r"\s+", " ", cpu_brand()).strip())


if __name__ == "__main__":
    print(host_key())
