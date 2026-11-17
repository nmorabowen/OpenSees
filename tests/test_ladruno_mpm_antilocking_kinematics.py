"""WP-232j: anti-locking kinematics and surfaces (gate rows L1, L6 spot rows, L8).

Written from tests/mpm_antilocking_interface.md (J1-J13, frozen 2026-11-17) and
tests/mpm_antilocking_cases.txt; never from tests/mpm_reference.py or
tests/mpm_implicit_reference.py.

- L1 a/b/c header-level: tests/mpm_antilocking_core_check.cpp (g++, every basis, 2D and 3D).
- L1 a/b through the engines: every `kind affine` case, mode F and mode A, `mpm debug affineCheck`.
- L8: `mpm info antiLocking` round trip (amendment 2); classic Tcl = Python on one L1 case and on a
  mode-F run with anti-locking on; loud refusals (keyword, mid-epoch change, a symmetric-storage
  system in mode A); the mode-A restart pack (version 3 round trip; a version-2 pack restores off).
- L6 spot: off leaves the counts at zero and the per-particle path at 0.
Cost on Zone-A: about 30 s (the g++ check about 10 s, 16 small affine cases, two short Tcl decks).
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_mpm_antilocking_harness as A  # noqa: E402
import ladruno_mpm_implicit_harness as H  # noqa: E402
from _cxx import cxx_runtime_flags  # noqa: E402

pytestmark = pytest.mark.zone_a

REPO = Path(__file__).resolve().parents[1]
_HERE = os.path.dirname(os.path.abspath(__file__))
_CORE = str(REPO / "SRC" / "domain" / "pattern" / "ladrunoMPM")
_KERNEL = str(REPO / "SRC" / "element" / "solidTransformation")


def _tcl_exe():
    exe = os.environ.get("LADRUNO_TCL_EXE")
    if exe:
        return exe if Path(exe).exists() else None
    cand = REPO / "dist" / "bin" / ("OpenSees.exe" if os.name == "nt" else "OpenSees")
    return str(cand) if cand.exists() else None


def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    if not hasattr(ops, "mpm"):
        pytest.skip("this opensees build has no `mpm` command")
    return ops


AFFINE = [c["name"] for c in A.parse_cases() if c["kind"] == "affine"]


def test_core_check_cpp():
    """L1 a, b, c and J1, J3 on the header (LadrunoMPMCore.h): affine fields give dFbar = dF to
    1e-12, a 30 deg rotation gives deps <= 1e-12, plane strain embeds with F_zz = 1, the corners
    lie in every basis's support."""
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "alc.exe")
        cmd = [gpp, "-O2", "-std=c++17", "-Wall", "-Wextra", "-I", _CORE, "-I", _KERNEL,
               os.path.join(_HERE, "mpm_antilocking_core_check.cpp"), "-o", exe,
               *cxx_runtime_flags()]
        cc = subprocess.run(cmd, stdin=subprocess.DEVNULL, capture_output=True, text=True,
                            encoding="utf-8", errors="replace")
        assert cc.returncode == 0, "compile failed: " + cc.stderr
        rr = subprocess.run([exe], stdin=subprocess.DEVNULL, capture_output=True, text=True,
                            encoding="utf-8", errors="replace", timeout=600)
    out = rr.stdout + rr.stderr
    assert rr.returncode == 0 and "RESULT: PASS" in out, out
    assert "[checks] " in out and int(out.split("[checks] ")[1].split()[0]) >= 600, out


@pytest.mark.parametrize("name", AFFINE)
def test_affine_engine(name):
    """L1 a, b through opensees.pyd (amendment 1: untruncated supports, `masscut 1e-12`)."""
    ops = _ops()
    c = A.case_by_name(name)
    rows = A.evaluate_affine(c, A.run_affine(c, ops))
    for lab, v, det in rows:
        print(f"{v.upper():6s} {lab}: {det}")
    assert rows and all(v == "pass" for _, v, _ in rows), rows


def _deck(cmds, tail):
    lines = [name + " " + " ".join(H._tcl_arg(a) for a in args) for name, args in cmds]
    return "\n".join(lines + tail + ['puts "DECK-DONE"']) + "\n"


def _run_tcl(text, name):
    exe = _tcl_exe()
    if exe is None:
        pytest.skip("classic OpenSees exe not found (build it, or set LADRUNO_TCL_EXE)")
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as td:
        return H._run_deck(text, exe, td, name, 600)


def _tcl_list(out, tag):
    for line in out.splitlines():
        t = line.split()
        if t and t[0] == tag:
            return [float(v) for v in t[1:]]
    raise AssertionError(f"no '{tag}' line in the Tcl output:\n{out[-3000:]}")


def test_affine_tcl_equals_python():
    """L8: classic Tcl and Python give the same L1 answer (affine_A_gimp_2d, every field)."""
    ops = _ops()
    c = A.case_by_name("affine_A_gimp_2d")
    py = A.run_affine(c, ops)
    nd = c["ndm"]
    tail = []
    for line in c["affine"]:
        ia, ib = line.index("A"), line.index("b")
        vals = line[ia + 1:ia + 1 + nd * nd] + line[ib + 1:ib + 1 + nd]
        tail.append(f'puts "AFF_{line[0]} [mpm debug affineCheck {" ".join(vals)}]"')
    out = _run_tcl(_deck(A.model_commands(c), tail), "affine_tcl")
    for line in c["affine"]:
        tc = _tcl_list(out, "AFF_" + line[0])
        pv = py[line[0]]["out"]
        assert tc[2:] == pv[2:], (line[0], tc, pv)            # the counts, exactly
        assert tc[0] <= 1e-12 and pv[0] <= 1e-12, (tc, pv)      # both within L1a
        assert abs(tc[0] - pv[0]) <= 1e-15 and abs(tc[1] - pv[1]) <= 1e-15, (tc, pv)


def _modef_cmds(antilocking):
    return [("wipe", []), ("model", ["basic", "-ndm", 2, "-ndf", 2]),
            ("nDMaterial", ["ElasticIsotropic", 1, 1.0e7, 0.45]),
            ("mpm", ["grid", H.HOST, "-origin", 0.0, 0.0, "-h", 0.5, "-cells", 4, 6]),
            ("mpm", ["body", 1, "box", 0.0, 0.0, 2.0, 2.5, "-ppc", 2, "-mat", 1, "-rho", 1000.0]),
            ("mpm", ["fix", "ymin", "xy"]), ("mpm", ["fix", "xmin", "x"]),
            ("mpm", ["gravity", 0.0, -1000.0]),
            ("mpm", ["antiLocking", antilocking])]


def test_modef_run_tcl_equals_python():
    """L8: a mode-F run with anti-locking on, 5 steps, is the same in classic Tcl and Python
    (positions, F, stress); the averaged / fallback counts agree; on differs from off."""
    ops = _ops()
    for name, args in _modef_cmds("fbar"):
        getattr(ops, name)(*H._py_args(args))
    ops.mpm("analyze", 5, 2e-4)
    py = {f: A._as_list(ops.mpm("particleData", f)) for f in ("x", "F", "sig")}
    py_counts = A._as_list(ops.mpm("info", "antiLockingCounts"))
    tail = ["mpm analyze 5 2e-4"] + [f'puts "P_{f} [mpm particleData {f}]"' for f in ("x", "F", "sig")]
    tail.append('puts "COUNTS [mpm info antiLockingCounts]"')
    out = _run_tcl(_deck(_modef_cmds("fbar"), tail), "modef_tcl")
    for f in ("x", "F", "sig"):
        tc = _tcl_list(out, "P_" + f)
        scale = max(abs(v) for v in py[f]) or 1.0
        assert len(tc) == len(py[f]) and max(abs(a - b) for a, b in zip(tc, py[f])) <= 1e-14 * scale, f
    assert _tcl_list(out, "COUNTS") == py_counts, (out[-500:], py_counts)
    assert py_counts[0] > 0, py_counts
    for name, args in _modef_cmds("off"):
        getattr(ops, name)(*H._py_args(args))
    ops.mpm("analyze", 5, 2e-4)
    off_sig = A._as_list(ops.mpm("particleData", "sig"))
    assert A._as_list(ops.mpm("info", "antiLockingCounts")) == [0.0, 0.0]          # L6 spot
    assert set(A._as_list(ops.mpm("particleData", "antiLocking"))) == {0.0}
    assert max(abs(a - b) for a, b in zip(off_sig, py["sig"])) > 0.0                # not the identity


def test_setting_round_trip_and_refusals():
    """L8 / K3: default off; `mpm info antiLocking` round trip; an unknown keyword and a change
    after the first step of an epoch are refused loudly; `mpm reset` re-opens the setting."""
    ops = _ops()
    cmds = [x for x in _modef_cmds("off") if not (x[0] == "mpm" and x[1][0] == "antiLocking")]
    for name, args in cmds:
        getattr(ops, name)(*H._py_args(args))
    assert int(ops.mpm("info", "antiLocking")) == 0
    ops.mpm("antiLocking", "fbar")
    assert int(ops.mpm("info", "antiLocking")) == 1
    with pytest.raises(Exception):
        ops.mpm("antiLocking", "on")
    ops.mpm("analyze", 1, 2e-4)
    with pytest.raises(Exception):
        ops.mpm("antiLocking", "off")
    assert int(ops.mpm("info", "antiLocking")) == 1
    ops.mpm("reset")
    ops.mpm("antiLocking", "off")
    assert int(ops.mpm("info", "antiLocking")) == 0


def _modea_cmds(system="UmfPack"):
    return _modef_cmds("fbar")[:-2] + [
        ("mpm", ["gravity", 0.0, -10.0]), ("mpm", ["mode", "quasiStatic"]),
        ("mpm", ["antiLocking", "fbar"]),
        ("constraints", ["Plain"]), ("numberer", ["RCM"]), ("system", [system]),
        ("test", ["NormDispIncr", 1e-12, 30]), ("algorithm", ["Newton"]),
        ("integrator", ["LoadControl", 0.25]), ("analysis", ["Static"])]


def test_symmetric_system_refused():
    """K3 / amendment 5: mode A with `fbar` refuses a symmetric-storage system loudly (the F-bar
    tangent is unsymmetric) and leaves the particles unmoved; UmfPack runs."""
    ops = _ops()
    for name, args in _modea_cmds("ProfileSPD"):
        getattr(ops, name)(*H._py_args(args))
    x0 = A._as_list(ops.mpm("particleData", "x"))
    with pytest.raises(Exception):
        ops.mpm("analyze", 1, 0.25)
    assert A._as_list(ops.mpm("particleData", "x")) == x0
    for name, args in _modea_cmds("UmfPack"):
        getattr(ops, name)(*H._py_args(args))
    ops.mpm("analyze", 1, 0.25)
    assert A._as_list(ops.mpm("info", "antiLockingCounts"))[0] > 0


def _restart(ops, wd, pack_version=None):
    for name, args in _modea_cmds():
        getattr(ops, name)(*H._py_args(args))
    for _ in range(2):
        ops.mpm("analyze", 1, 0.25)
    if pack_version is not None:
        ops.mpm("debug", "packVersion", pack_version)
    ops.database("File", str(wd / "aldb"))
    ops.save(2)
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.database("File", str(wd / "aldb"))
    ops.restore(2)
    for name, args in _modea_cmds()[-7:]:
        getattr(ops, name)(*H._py_args(args))


def test_restart_pack_v3_and_v2():
    """J13 / L8: a version-3 database restores with the setting and continues as the
    uninterrupted run; a version-2 pack (the 232i layout) restores with anti-locking off."""
    ops = _ops()
    for name, args in _modea_cmds():
        getattr(ops, name)(*H._py_args(args))
    for _ in range(4):
        ops.mpm("analyze", 1, 0.25)
    ref = A._as_list(ops.mpm("particleData", "sig"))
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as td:
        wd = Path(td)
        _restart(ops, wd)
        assert int(ops.mpm("info", "antiLocking")) == 1
        for _ in range(2):
            ops.mpm("analyze", 1, 0.25)
        got = A._as_list(ops.mpm("particleData", "sig"))
        scale = max(abs(v) for v in ref)
        assert max(abs(a - b) for a, b in zip(got, ref)) <= 1e-10 * scale
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as td:
        _restart(ops, Path(td), pack_version=2)
        assert int(ops.mpm("info", "antiLocking")) == 0
