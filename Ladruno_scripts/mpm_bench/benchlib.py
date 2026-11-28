"""WP-253: what the benchmark scripts share on top of the M0b harness.

The M0b harness (`Ladruno_scripts/mpm_m0b/common.py`, WP-232g) is REUSED, not forked: fresh child
process per run (`--child cfg.json`), affinity mask, thread pins before import, the profiler HDF5
reader (`safe_read_profile`), `median_spread`, result-JSON layout. It is loaded here under the
module name `m0b_common` (its sibling `models.py` would otherwise collide with ours).

Added here: the provenance every row must carry -- `ladrunoBuild` of OpenSees.exe AND of the
.pyd/.so (refuse on a mismatch), machine (CPU model, cores, RAM, OS), MKL version, thread count.
"""
from __future__ import annotations

import ctypes
import importlib.util
import os
import platform
import re
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
M0B = HERE.parent / "mpm_m0b"
IS_WIN = os.name == "nt"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


C = _load("m0b_common", M0B / "common.py")           # the M0b harness, as is
M0B_MODELS = _load("m0b_models", M0B / "models.py")    # its constants (block / column parameters)

DEFAULT_PYD = REPO / "dist" / "bin" / ("opensees.pyd" if IS_WIN else "opensees.so")
DEFAULT_EXE = REPO / "dist" / "bin" / ("OpenSees.exe" if IS_WIN else "OpenSees")


# --------------------------------------------------------------------------- #
#  ladrunoBuild of the executable                                             #
# --------------------------------------------------------------------------- #
def exe_build(exe):
    """`[ladrunoBuild]` printed by OpenSees.exe on a one-line Tcl script; None if it cannot be
    run or printed nothing recognisable (then the caller refuses)."""
    exe = Path(exe).resolve()
    if not exe.exists():
        return None, f"no executable at {exe}"
    with tempfile.TemporaryDirectory() as td:
        tcl = Path(td) / "build.tcl"
        tcl.write_text('puts "LADRUNOBUILD=[ladrunoBuild]"\nexit\n', encoding="utf-8")
        env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1")
        try:
            p = subprocess.run([str(exe), str(tcl)], capture_output=True, text=True, timeout=120,
                               cwd=td, env=env, stdin=subprocess.DEVNULL)
        except Exception as e:                       # noqa: BLE001
            return None, repr(e)
    m = re.search(r"LADRUNOBUILD=(\S+)", p.stdout + p.stderr)
    return (m.group(1) if m else None), (None if m else "no LADRUNOBUILD= in output")


# --------------------------------------------------------------------------- #
#  machine                                                                    #
# --------------------------------------------------------------------------- #
def cpu_model():
    if IS_WIN:
        try:
            import winreg
            with winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE,
                                r"HARDWARE\DESCRIPTION\System\CentralProcessor\0") as k:
                return winreg.QueryValueEx(k, "ProcessorNameString")[0].strip()
        except OSError:
            return platform.processor()
    try:
        for ln in Path("/proc/cpuinfo").read_text().splitlines():
            if ln.lower().startswith("model name"):
                return ln.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor()


def ram_bytes():
    if IS_WIN:
        return C.total_ram_bytes()
    try:
        for ln in Path("/proc/meminfo").read_text().splitlines():
            if ln.startswith("MemTotal"):
                return int(ln.split()[1]) * 1024
    except OSError:
        pass
    return None


def machine_info():
    return {"cpu_model": cpu_model(), "logical_cores": os.cpu_count(),
            "ram_bytes": ram_bytes(), "os": platform.platform(), "hostname": platform.node()}


# --------------------------------------------------------------------------- #
#  MKL version                                                                #
# --------------------------------------------------------------------------- #
def _file_version(path):
    """Windows file version resource 'a.b.c.d' of a DLL, or None."""
    if not IS_WIN:
        return None
    try:
        ver = ctypes.WinDLL("version")
        size = ver.GetFileVersionInfoSizeW(str(path), None)
        if not size:
            return None
        buf = ctypes.create_string_buffer(size)
        if not ver.GetFileVersionInfoW(str(path), 0, size, buf):
            return None
        ptr, ln = ctypes.c_void_p(), ctypes.c_uint()
        if not ver.VerQueryValueW(buf, "\\", ctypes.byref(ptr), ctypes.byref(ln)) or not ptr.value:
            return None
        v = (ctypes.c_uint32 * 8).from_address(ptr.value)   # VS_FIXEDFILEINFO: [sig, struc, msV, lsV...]
        ms, ls = v[2], v[3]                                 # dwFileVersionMS / LS
        return f"{ms >> 16}.{ms & 0xFFFF}.{ls >> 16}.{ls & 0xFFFF}"
    except Exception:                                       # noqa: BLE001
        return None


def _loaded_modules():
    """Paths of the modules mapped into THIS process."""
    if not IS_WIN:
        try:
            return sorted({ln.split()[-1] for ln in Path("/proc/self/maps").read_text().splitlines()
                           if ln.split()[-1].startswith("/")})
        except OSError:
            return []
    try:
        psapi = ctypes.WinDLL("psapi")
        k32 = ctypes.WinDLL("kernel32")
        k32.GetCurrentProcess.restype = ctypes.c_void_p
        h = k32.GetCurrentProcess()
        mods = (ctypes.c_void_p * 2048)()
        need = ctypes.c_ulong()
        psapi.EnumProcessModules.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_ulong,
                                             ctypes.POINTER(ctypes.c_ulong)]
        psapi.GetModuleFileNameExW.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_wchar_p,
                                               ctypes.c_ulong]
        if not psapi.EnumProcessModules(h, mods, ctypes.sizeof(mods), ctypes.byref(need)):
            return []
        out = []
        for i in range(need.value // ctypes.sizeof(ctypes.c_void_p)):
            b = ctypes.create_unicode_buffer(1024)
            if psapi.GetModuleFileNameExW(h, mods[i], b, 1024):
                out.append(b.value)
        return out
    except Exception:                                       # noqa: BLE001
        return []


def mkl_info(search_dir=None):
    """The MKL this process actually loaded (call it AFTER a linear solver has run, so a lazily
    loaded MKL is in the process), with `mkl_get_version_string` when mkl_rt is mapped; plus the
    file versions of every mkl_*.dll next to the binary (`search_dir`) as a fallback."""
    mods = [m for m in _loaded_modules() if "mkl" in os.path.basename(m).lower()]
    info = {"loaded": {os.path.basename(m): _file_version(m) for m in mods}, "version_string": None}
    for m in mods:
        if "mkl_rt" in os.path.basename(m).lower():
            try:
                lib = ctypes.CDLL(m)
                buf = ctypes.create_string_buffer(512)
                lib.mkl_get_version_string(buf, 512)
                info["version_string"] = buf.value.decode(errors="replace").strip()
            except Exception as e:                          # noqa: BLE001
                info["version_string_error"] = repr(e)
            break
    if search_dir and Path(search_dir).is_dir():
        pats = ("mkl_*.dll", "mkl_*.so*", "libmkl_*.so*")
        files = sorted({f for p in pats for f in Path(search_dir).glob(p)})
        info["binary_dir"] = {f.name: _file_version(f) for f in files[:40]}
    return info


def mkl_summary(info):
    """One short string for the table: the version string if there is one, else the file
    version of the first mkl_core / mkl_rt seen."""
    if not info:
        return None
    if info.get("version_string"):
        return info["version_string"]
    for src in (info.get("loaded") or {}, info.get("binary_dir") or {}):
        for name, v in src.items():
            if v and ("core" in name.lower() or "_rt" in name.lower()):
                return f"{name} {v}"
    return None


# --------------------------------------------------------------------------- #
#  result helpers                                                             #
# --------------------------------------------------------------------------- #
def first(v):
    """ops.mpm returns a scalar or a one-element list, depending on the command."""
    return float(v[0] if isinstance(v, (list, tuple)) else v)


def row_key(c):
    """The identity of a table row: everything but the repeat and the profile kind."""
    return (c["mode"], c["dim"], c["basis"], c["mat"], c["N"], c["threads"])


def key_str(k):
    mode, dim, basis, mat, N, th = k
    return f"{mode}/{dim}D/{basis}/{mat}/N{N:g}/t{th}"
