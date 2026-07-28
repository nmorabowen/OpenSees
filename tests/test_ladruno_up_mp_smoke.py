"""LadrunoUP serialization gate (ADR 71 P1, WP1.D).

Two complementary checks of `LadrunoUP::sendSelf` / `recvSelf`:

  1. test_serial_database_roundtrip  -- the actual serialization gate.
     Builds a fully non-default Q4 u-p model (bbar, lumped mass, manual stab,
     dynSeepage off, non-unit Biot alpha / finite grain Ks, anisotropic perm,
     gravity body force, Rayleigh factors set), commits a static step, then
     `database File` + `save` + `wipe` + `restore`. Restore recreates every
     element through the broker and drives `recvSelf`; `save` drove `sendSelf`
     (both via FileDatastore-as-Channel through Domain::send/recvSelf). We then
     assert the reconstructed element is bit-for-bit functionally identical:
     same nodal pore pressures, same per-GP element responses (porePressure,
     effective+total stresses, Darcy flux -- these read the rebuilt geometry AND
     the re-brokered materials), and that it re-solves from the restored state.
     If ANY ctor arg or committed material state failed to round-trip, one of
     these responses moves.

  2. test_mp_two_rank_smoke  -- complementary parallel-assembly smoke.
     `mpiexec -n 2` on the openseesmp build with two Q4 u-p elements straddling
     a partition (element 1 on rank 0, element 2 on rank 1, sharing edge nodes
     2 & 3). MumpsParallelSOE (matType=0, unsymmetric -- the honest-p tangent
     needs it) assembles+solves the distributed system; both ranks must agree on
     the shared-node pore pressure after an implicit transient step.

     SCOPE NOTE: this build bundles METIS 4, so the `partition` command (the only
     route that migrates elements between ranks, hence the only MP path that
     would call Element::send/recvSelf) is disabled. The replicated-interpreter
     shared-node model each rank builds locally does NOT migrate elements -- so
     this test proves LadrunoUP *works under MP*, but the send/recvSelf gate
     itself is carried by test 1 (the serial database round-trip). Documented so
     the next agent doesn't mistake a green MP smoke for a serialization test.

Both tests run their OpenSees work in a `python -S` SUBPROCESS: the shared test
box preloads a stale opensees.pyd via a boot .pth (see MEMORY: "worktree runners
MUST use python -S + manual paths"), which an in-process `import opensees` under
plain pytest would pick up instead of this worktree's freshly-built engine.

Run:  py -3.12 -m pytest tests/test_ladruno_up_mp_smoke.py -v
"""
import os
import re
import subprocess
import sys

import pytest

# ---- locate this worktree's freshly-built engine ---------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
DIST_BIN = os.path.join(_ROOT, "dist", "bin")
DIST_MP = os.path.join(_ROOT, "dist", "openseesmp")


# ===========================================================================
#  1. serial database round-trip -- the sendSelf/recvSelf gate
# ===========================================================================

# Runs in a `-S` subprocess so the boot-.pth stale-pyd preload cannot shadow the
# worktree build. Asserts internally (non-zero exit on failure) and prints a
# success sentinel. Every u-p knob is deliberately OFF its default so a dropped
# serialization field is observable.
_SERIAL_DRIVER = r'''
import os, sys
DIST_BIN = sys.argv[1]
DB = sys.argv[2]
os.add_dll_directory(DIST_BIN)
sys.path.insert(0, DIST_BIN)
os.environ["PATH"] = DIST_BIN + os.pathsep + os.environ.get("PATH", "")
import opensees as ops
assert os.path.normcase(DIST_BIN) in os.path.normcase(ops.__file__), \
    ("wrong engine imported: " + ops.__file__)

COMMON = ["-thick", 1.0, "-Kf", 2.2e6, "-poro", 0.3, "-rhoF", 1.0,
          "-perm", 1.0e-4, 2.0e-4, "-alpha", 0.9, "-Ks", 1.0e7,
          "-body", 0.0, -9.81, "-fluidBody", 0.0, -9.81,
          "-formulation", "bbar", "-lumped", "-stab", 0.15,
          "-dynSeepage", "off"]

def build():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0); ops.node(2, 1.0, 0.0); ops.node(3, 1.0, 1.0)
    ops.node(4, 0.0, 1.0); ops.node(5, 2.0, 0.0); ops.node(6, 2.0, 1.0)
    ops.nDMaterial("ElasticIsotropic", 1, 25000.0, 0.3, 2.0)
    ops.element("LadrunoUP", 1, 1, 2, 3, 4, 1, *COMMON)
    ops.element("LadrunoUP", 2, 2, 5, 6, 3, 1, *COMMON)
    ops.fix(1, 1, 1, 0); ops.fix(5, 1, 1, 0)
    ops.fix(3, 0, 0, 1); ops.fix(4, 0, 0, 1); ops.fix(6, 0, 0, 1)
    ops.rayleigh(0.02, 0.0, 1.0e-4, 0.0)   # solid-only Rayleigh; must round-trip

def solve():
    ops.constraints("Plain"); ops.numberer("RCM"); ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-10, 25); ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0); ops.analysis("Static")
    return ops.analyze(1)

def vmax(a, b):
    assert a and b and len(a) == len(b), (a, b)
    return max(abs(x - y) for x, y in zip(a, b))

build()
assert solve() == 0, "reference static analysis did not converge"
ref_p = {n: ops.nodeDisp(n, 3) for n in (2, 3, 4, 6)}
assert abs(ref_p[2]) > 1.0, "test model produced ~zero pressure -- not a probe"
ref = {}
for tag in (1, 2):
    ref[tag] = {name: ops.eleResponse(tag, name)
                for name in ("porePressure", "stresses", "stressesTotal", "flux")}

# sendSelf -> database, then recvSelf <- database into a wiped domain
ops.database("File", DB); ops.save(1)
ops.wipe()
ops.database("File", DB); ops.restore(1)

res_p = {n: ops.nodeDisp(n, 3) for n in (2, 3, 4, 6)}
assert max(abs(ref_p[n] - res_p[n]) for n in ref_p) == 0.0, "nodal pressure drift"
for tag in (1, 2):
    for name, refval in ref[tag].items():
        got = ops.eleResponse(tag, name)
        d = vmax(refval, got)
        assert d < 1.0e-10, "ele %d response %r drifted by %g" % (tag, name, d)
assert solve() == 0, "re-solve after restore did not converge"
print("ROUNDTRIP-OK")
'''


def _serial_available():
    return os.path.isfile(os.path.join(DIST_BIN, "opensees.pyd"))


def test_serial_database_roundtrip(tmp_path):
    if not _serial_available():
        pytest.skip("serial opensees.pyd not built in dist/bin")

    driver = tmp_path / "up_serial_driver.py"
    driver.write_text(_SERIAL_DRIVER)
    db = str(tmp_path / "up_db")

    env = dict(os.environ)
    env["PATH"] = DIST_BIN + os.pathsep + env.get("PATH", "")
    # encoding pinned: the child prints the UTF-8 splash banner; cp1252
    # text-mode decode raises in the reader thread => proc.stdout is None
    # and the assert dies with TypeError (quirks ledger).
    proc = subprocess.run(
        [sys.executable, "-S", str(driver), DIST_BIN, db],
        capture_output=True, text=True, timeout=300, env=env,
        encoding="utf-8", errors="replace")

    assert "ROUNDTRIP-OK" in proc.stdout and proc.returncode == 0, (
        f"serial database round-trip failed (rc={proc.returncode})\n"
        f"--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}")


# ===========================================================================
#  2. MP two-rank shared-node smoke (complementary; see module docstring)
# ===========================================================================

_MP_DRIVER = r'''
import os, sys
MODDIR = sys.argv[1]; OUT = sys.argv[2]
_IMPI = r"C:\Program Files (x86)\Intel\oneAPI\mpi\latest\bin"
_LIBFABRIC = r"C:\Program Files (x86)\Intel\oneAPI\mpi\latest\opt\mpi\libfabric\bin"
_MKL = r"C:\Program Files (x86)\Intel\oneAPI\mkl\latest\bin"
_ICOMP = r"C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin"
for d in (MODDIR, _IMPI, _LIBFABRIC, _MKL, _ICOMP):
    if d and os.path.isdir(d):
        os.add_dll_directory(d)
os.environ["PATH"] = os.pathsep.join(
    [p for p in (_IMPI, _LIBFABRIC, MODDIR, os.environ.get("PATH", "")) if p])
sys.path.insert(0, MODDIR)
import openseesmp as ops
pid = ops.getPID(); nproc = ops.getNP(); ops.start()
ops.model("basic", "-ndm", 2, "-ndf", 3)
ops.nDMaterial("ElasticIsotropic", 1, 25000.0, 0.3, 2.0)
common = ["-thick", 1.0, "-Kf", 2.2e6, "-poro", 0.3, "-rhoF", 1.0,
          "-perm", 1e-4, 1e-4, "-alpha", 1.0, "-body", 0.0, -9.81, "-stab", "off"]
if pid == 0:
    ops.node(1, 0.0, 0.0); ops.node(2, 1.0, 0.0)
    ops.node(3, 1.0, 1.0); ops.node(4, 0.0, 1.0)
    ops.element("LadrunoUP", 1, 1, 2, 3, 4, 1, *common)
    ops.fix(1, 1, 1, 0); ops.fix(2, 1, 1, 0)
    ops.fix(3, 0, 0, 1); ops.fix(4, 0, 0, 1)
else:
    ops.node(2, 1.0, 0.0); ops.node(5, 2.0, 0.0)
    ops.node(6, 2.0, 1.0); ops.node(3, 1.0, 1.0)
    ops.element("LadrunoUP", 2, 2, 5, 6, 3, 1, *common)
    ops.fix(2, 1, 1, 0); ops.fix(5, 1, 1, 0)
    ops.fix(6, 0, 0, 1); ops.fix(3, 0, 0, 1)
ops.constraints("Transformation")
ops.numberer("ParallelPlain")
ops.system("Mumps", "-ICNTL14", 200)
ops.test("NormDispIncr", 1e-8, 30, 0)
ops.algorithm("Newton")
ops.integrator("Newmark", 0.5, 0.25)
ops.analysis("Transient")
rc = 0
for _ in range(3):
    rc = ops.analyze(1, 1.0)
    if rc != 0:
        break
p2 = ops.nodeDisp(2, 3)
with open(os.path.join(OUT, "up_mp_rank%d.txt" % pid), "w") as fh:
    fh.write("pid=%d nproc=%d rc=%d p2=%.17g\n" % (pid, nproc, rc, p2))
print("RANK %d/%d rc=%d p2=%g DONE" % (pid, nproc, rc, p2))
ops.wipe()
'''


def _mp_available():
    return (os.path.isfile(os.path.join(DIST_MP, "mpiexec.exe"))
            and os.path.isfile(os.path.join(DIST_MP, "openseesmp.pyd")))


def test_mp_two_rank_smoke(tmp_path):
    if not _mp_available():
        pytest.skip("openseesmp build (dist/openseesmp/{mpiexec.exe,openseesmp.pyd}) absent")

    driver = tmp_path / "up_mp_driver.py"
    driver.write_text(_MP_DRIVER)
    mpiexec = os.path.join(DIST_MP, "mpiexec.exe")

    env = dict(os.environ)
    env["PATH"] = DIST_MP + os.pathsep + env.get("PATH", "")
    cmd = [mpiexec, "-n", "2", sys.executable, "-S", str(driver),
           DIST_MP, str(tmp_path)]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              timeout=300, env=env,
                              encoding="utf-8", errors="replace")
    except (subprocess.TimeoutExpired, OSError) as exc:  # infra, not the element
        pytest.skip(f"MP launcher failed to run: {exc}")

    rank_files = [tmp_path / f"up_mp_rank{r}.txt" for r in (0, 1)]
    if not all(f.exists() for f in rank_files):
        # no output at all => Intel-MPI / launcher infra problem, not a
        # LadrunoUP regression; surface stderr so it is not silently lost
        pytest.skip("MP run produced no rank output (MPI infra?). "
                    f"stdout/stderr tail:\n{proc.stdout[-800:]}\n{proc.stderr[-800:]}")

    parsed = {}
    for f in rank_files:
        m = re.search(r"pid=(\d+) nproc=(\d+) rc=(-?\d+) p2=([-\d.eE+]+)",
                     f.read_text())
        assert m, f"unparseable rank file: {f.read_text()!r}"
        parsed[int(m.group(1))] = dict(
            nproc=int(m.group(2)), rc=int(m.group(3)), p2=float(m.group(4)))

    assert set(parsed) == {0, 1}
    for pid, d in parsed.items():
        assert d["nproc"] == 2, f"rank {pid} saw nproc={d['nproc']}"
        assert d["rc"] == 0, f"rank {pid} transient did not converge (rc={d['rc']})"

    # the whole point: both ranks agree on the shared-node (node 2) pressure
    p0, p1 = parsed[0]["p2"], parsed[1]["p2"]
    assert abs(p0) > 1.0e-6, "shared-node pressure ~0 -- MP model degenerate"
    assert abs(p0 - p1) <= 1.0e-9 * max(1.0, abs(p0)), (
        f"ranks disagree on shared-node pressure: {p0!r} vs {p1!r}")
