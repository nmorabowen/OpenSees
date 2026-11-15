"""WP-232c phase c2: the two MPM file contracts through the shipped binaries
(Ladruno_implementation/232a_mpm_contracts.md; ADR 232 §8 S1 rows "contract round
trips", "particle-id stability across outputs", "seeding guard", "recorder
synchrony").

* `mpm particles -file <path>` reads 232a's Contract-1 fixture
  (tests/data/mpm/particle_set_v1.h5) all or nothing, opening the path as given,
  and refuses every malformed variant loudly, naming the rule;
* `mpm output` writes the .ladruno particle group (FORMAT_VERSION 1,
  FORMAT_MINOR_VERSION 2): every file passes ladruno_format.validate() and
  tests/mpm_contracts.py::validate_particle_results;
* the round trip Contract 1 -> run -> Contract 2 keeps SET_ID, bodies, units,
  the registry and output 0 (232a §3);
* positions and stresses in an output are the committed state of the same step;
* particle ids are identical across outputs; `mpm reset` opens a new stage with
  PARTICLE_EPOCH + 1.

Needs h5py (Zone-A gets it with #1077).
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "Ladruno_scripts" / "ladruno_recorder_tests"))
import ladruno_format as lf  # noqa: E402
import mpm_contracts as mc  # noqa: E402

pytestmark = pytest.mark.zone_a

REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "tests" / "data" / "mpm" / "particle_set_v1.h5"
TAG = 232
FIELDS = ("POSITION", "VELOCITY", "STRESS", "STRAIN", "VOLUME", "DEFORMATION_GRADIENT",
          "GIMP_HALF_LENGTH")


def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    if not hasattr(ops, "mpm"):
        pytest.skip("this opensees build has no `mpm` command")
    return ops


def _tcl_exe():
    exe = os.environ.get("LADRUNO_TCL_EXE")
    if exe:
        return exe if Path(exe).exists() else None
    cand = REPO / "dist" / "bin" / ("OpenSees.exe" if os.name == "nt" else "OpenSees")
    return str(cand) if cand.exists() else None


def _deck(ops, mat2=("ElasticIsotropic", 2, 2.0e4, 0.3), grid=None):
    """The fixture's model: plane strain, 4 x 4 cells of 0.5 m, uGIMP; tags 1 and 2."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 2.0e5, 0.3)
    ops.nDMaterial(*mat2)
    g = grid or ("-origin", 0.0, 0.0, "-h", 0.5, "-cells", 4, 4)
    ops.mpm("grid", TAG, *g, "-basis", "gimp")


def _set():
    with h5py.File(FIXTURE, "r") as f:
        P = f["PARTICLES"]
        return {
            "SET_ID": f.attrs["SET_ID"].decode() if isinstance(f.attrs["SET_ID"], bytes) else str(f.attrs["SET_ID"]),
            "ID": P["ID"][...], "POSITION": P["POSITION"][...], "VOLUME": P["VOLUME"][...],
            "MASS": P["DENSITY"][...] * P["VOLUME"][...], "STRESS": P["STRESS"][...],
            "BODY_ID": P["BODY_ID"][...], "MATERIAL_TAG": P["MATERIAL_TAG"][...],
            "VELOCITY": P["VELOCITY"][...] if "VELOCITY" in P else np.zeros_like(P["POSITION"][...]),
            "BODIES": [(int(i), n.decode() if isinstance(n, bytes) else str(n))
                       for i, n in zip(f["BODIES/ID"][...], f["BODIES/NAME"][...])],
        }


def _run_and_write(ops, out, steps=5, every=2, dt=1.0e-4):
    ops.mpm("output", str(out), "-every", every, "-fields", *FIELDS)
    ops.mpm("analyze", steps, dt)
    ops.mpm("output", "-close")


# --------------------------------------------------------------------------- #
#  Contract 1 in                                                                #
# --------------------------------------------------------------------------- #
def test_loads_the_contract1_fixture_exactly():
    ops = _ops()
    _deck(ops)
    n = ops.mpm("particles", "-file", str(FIXTURE))
    S = _set()
    assert int(n) == len(S["ID"]) == 12
    ids = np.array(ops.mpm("particleData", "id"))
    order = np.argsort(ids)
    assert np.array_equal(ids[order].astype(np.int64), S["ID"])
    x = np.array(ops.mpm("particleData", "x")).reshape(-1, 2)[order]
    assert np.array_equal(x, S["POSITION"])
    assert np.array_equal(np.array(ops.mpm("particleData", "V0"))[order], S["VOLUME"])
    assert np.array_equal(np.array(ops.mpm("particleData", "m"))[order], S["MASS"])
    assert np.array_equal(np.array(ops.mpm("particleData", "sig0")).reshape(-1, 3)[order], S["STRESS"])


def test_opens_a_relative_path_as_given(tmp_path):
    ops = _ops()
    shutil.copy(FIXTURE, tmp_path / "case1.particles.h5")
    _deck(ops)
    here = os.getcwd()
    try:
        os.chdir(tmp_path)
        assert int(ops.mpm("particles", "-file", "case1.particles.h5")) == 12
    finally:
        os.chdir(here)


def _variant(tmp_path, name, edit):
    p = tmp_path / f"{name}.particles.h5"
    shutil.copy(FIXTURE, p)
    with h5py.File(p, "r+") as f:
        edit(f)
    return p


def _replace(f, path, data):
    del f[path]
    f[path] = data


REFUSALS = {
    "format": ("C1-FMT", lambda f: f.attrs.__setitem__("FORMAT", "something-else")),
    "version": ("C1-VER", lambda f: f.attrs.__setitem__("FORMAT_VERSION", np.int32(2))),
    "setid": ("C1-SETID", lambda f: f.attrs.__setitem__("SET_ID", "not-a-uuid")),
    "volume": ("C1-VOL", lambda f: _replace(f, "PARTICLES/VOLUME", -f["PARTICLES/VOLUME"][...])),
    "unknown": ("C1-UNKNOWN", lambda f: f.create_dataset("PARTICLES/COLOR", data=np.zeros(12))),
    "flags": ("C1-FLAGS", lambda f: f.create_dataset("PARTICLES/FLAGS", data=np.ones(12, dtype="<u4"))),
    "id_dtype": ("C1-DTYPE", lambda f: _replace(f, "PARTICLES/ID", f["PARTICLES/ID"][...].astype("<i4"))),
    "id_repeat": ("C1-ID", lambda f: _replace(f, "PARTICLES/ID", np.r_[0, np.arange(11)].astype("<i8"))),
    "n": ("C1-N", lambda f: f.attrs.__setitem__("NUM_PARTICLES", np.int64(13))),
    "mat": ("C1-MAT", lambda f: _replace(f, "PARTICLES/MATERIAL_TAG", np.full(12, 77, dtype="<i4"))),
    "unequal_h": ("CELL_SIZE components differ",
                  lambda f: (f["GRID"].attrs.__setitem__("CELL_SIZE", np.array([0.5, 0.25])),
                             f["GRID"].attrs.__setitem__("NUM_CELLS", np.array([4, 8], dtype="<i4")))),
}


@pytest.mark.parametrize("case", sorted(REFUSALS))
def test_refuses_a_malformed_set_and_loads_nothing(case, tmp_path, capfd):
    ops = _ops()
    code, edit = REFUSALS[case]
    p = _variant(tmp_path, case, edit)
    _deck(ops)
    with pytest.raises(Exception):
        ops.mpm("particles", "-file", str(p))
    assert int(ops.mpm("info", "numParticles")) == 0
    assert code in capfd.readouterr().err


def test_refuses_a_grid_that_differs_from_the_deck(capfd):
    ops = _ops()
    _deck(ops, grid=("-origin", 0.0, 0.0, "-h", 0.5, "-cells", 5, 4))
    with pytest.raises(Exception):
        ops.mpm("particles", "-file", str(FIXTURE))
    assert int(ops.mpm("info", "numParticles")) == 0
    assert "differs from the deck" in capfd.readouterr().err


def test_refuses_stress33_on_a_non_elastic_material(capfd):
    ops = _ops()
    _deck(ops, mat2=("DruckerPrager", 2, 1.0e4, 5.0e3, 1.0e3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
    with pytest.raises(Exception):
        ops.mpm("particles", "-file", str(FIXTURE))
    assert int(ops.mpm("info", "numParticles")) == 0
    assert "non-elastic" in capfd.readouterr().err


def test_refuses_a_second_set_whose_ids_meet_the_epoch(capfd):
    ops = _ops()
    _deck(ops)
    ops.mpm("particles", "-file", str(FIXTURE))
    with pytest.raises(Exception):
        ops.mpm("particles", "-file", str(FIXTURE))
    assert int(ops.mpm("info", "numParticles")) == 12
    assert "already used" in capfd.readouterr().err


# --------------------------------------------------------------------------- #
#  Contract 2 out                                                               #
# --------------------------------------------------------------------------- #
def test_round_trip_contract1_run_contract2(tmp_path):
    ops = _ops()
    _deck(ops)
    ops.mpm("particles", "-file", str(FIXTURE))
    out = tmp_path / "run.ladruno"
    _run_and_write(ops, out, steps=5, every=2)
    assert lf.validate(str(out)) == []
    assert mc.validate_particle_results(str(out)) == []
    S = _set()
    with h5py.File(out, "r") as f:
        info = f["INFO"].attrs
        assert int(np.atleast_1d(info["FORMAT_VERSION"])[0]) == 1
        assert int(np.atleast_1d(info["FORMAT_MINOR_VERSION"])[0]) == 2
        st = f["MODEL_STAGE[0]"]
        assert lf._attr(st, "KIND") == "transient"
        assert int(lf._attr(st, "EMPTY_PARTITION")) == 1
        assert int(lf._attr(st, "PARTICLE_EPOCH")) == 0
        mp = st["MODEL/PARTICLES"]
        assert lf._attr(mp, "SOURCE_SET_ID") == S["SET_ID"]
        assert lf._attr(mp, "MODE") == "explicit"
        assert int(lf._attr(mp, "PATTERN_TAG")) == TAG
        assert (lf._attr(mp, "UNIT_LENGTH"), lf._attr(mp, "UNIT_FORCE"), lf._attr(mp, "UNIT_TIME")) == ("m", "kN", "s")
        assert float(lf._attr(mp, "PLANE_STRAIN_THICKNESS")) == 1.0
        names = [(int(i), n) for i, n in zip(mp["BODIES/ID"][...], lf.str_dataset(mp["BODIES/NAME"]))]
        assert sorted(names) == sorted(S["BODIES"])
        for k in ("ID", "BODY_ID", "MATERIAL_TAG", "MASS"):
            assert np.array_equal(mp[k][...], S[k]), k
        op = st["RESULTS/ON_PARTICLES"]
        assert list(op["STEP"][...]) == sorted(op["STEP"][...]) and len(op["STEP"]) == 4   # 0, 2, 4, 5
        a, b = op["OFFSETS"][0], op["OFFSETS"][1]
        assert np.array_equal(op["ID"][a:b], S["ID"])
        assert np.array_equal(op["POSITION/DATA"][a:b], S["POSITION"])
        assert np.array_equal(op["VOLUME/DATA"][a:b, 0], S["VOLUME"])
        assert np.array_equal(op["VELOCITY/DATA"][a:b], S["VELOCITY"])
        assert np.array_equal(op["STRESS/DATA"][a:b], S["STRESS"])


def test_outputs_are_synchronous_and_ids_stable(tmp_path):
    """Output k holds the committed positions AND stresses of the same step."""
    ops = _ops()
    _deck(ops)
    ops.mpm("particles", "-file", str(FIXTURE))
    out = tmp_path / "sync.ladruno"
    ops.mpm("output", str(out), "-every", 1, "-fields", "POSITION", "STRESS")
    ids = np.array(ops.mpm("particleData", "id"))
    order = np.argsort(ids)
    seen = []
    for _ in range(4):
        ops.mpm("analyze", 1, 2.0e-4)
        seen.append((np.array(ops.mpm("particleData", "x")).reshape(-1, 2)[order],
                     np.array(ops.mpm("particleData", "sig")).reshape(-1, 3)[order]))
    ops.mpm("output", "-close")
    assert lf.validate(str(out)) == []
    with h5py.File(out, "r") as f:
        op = f["MODEL_STAGE[0]/RESULTS/ON_PARTICLES"]
        off, allids = op["OFFSETS"][...], op["ID"][...]
        assert len(off) == 6          # output 0 (initial) + 4 steps
        for k in range(5):
            assert np.array_equal(allids[off[k]:off[k + 1]], ids[order].astype(np.int64))
        for k, (x, s) in enumerate(seen, start=1):
            a, b = off[k], off[k + 1]
            assert np.array_equal(op["POSITION/DATA"][a:b], x), k
            assert np.array_equal(op["STRESS/DATA"][a:b], s), k
        assert np.any(op["STRESS/DATA"][off[4]:off[5]] != op["STRESS/DATA"][off[0]:off[1]])


def test_mpm_reset_opens_a_new_epoch_stage(tmp_path):
    ops = _ops()
    _deck(ops)
    ops.mpm("particles", "-file", str(FIXTURE))
    out = tmp_path / "epochs.ladruno"
    ops.mpm("output", str(out), "-every", 2)
    ops.mpm("analyze", 3, 1.0e-4)
    ops.mpm("reset")
    ops.mpm("analyze", 2, 1.0e-4)
    ops.mpm("output", "-close")
    assert lf.validate(str(out)) == []
    with h5py.File(out, "r") as f:
        stages = sorted(k for k in f if k.startswith("MODEL_STAGE"))
        assert stages == ["MODEL_STAGE[0]", "MODEL_STAGE[1]"]
        assert [int(lf._attr(f[s], "PARTICLE_EPOCH")) for s in stages] == [0, 1]


def test_lattice_bodies_write_a_valid_group(tmp_path):
    """No Contract-1 file: units from `mpm units`, bodies named body<id>."""
    ops = _ops()
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e5, 0.3)
    ops.mpm("grid", TAG, "-origin", 0.0, 0.0, 0.0, "-h", 0.5, "-cells", 4, 4, 4)
    ops.mpm("body", 7, "box", 0.5, 0.5, 0.5, 1.5, 1.5, 1.5, "-ppc", 2, "-mat", 1, "-rho", 1000.0)
    ops.mpm("velocity", 7, "uniform", 0.1, 0.0, -0.1)
    ops.mpm("units", "m", "N", "s")
    out = tmp_path / "lattice.ladruno"
    ops.mpm("output", str(out), "-dt", 3.0e-4, "-fields", *FIELDS)
    ops.mpm("analyze", 10, 1.0e-4)
    ops.mpm("output", "-close")
    assert mc.validate_particle_results(str(out)) == []
    with h5py.File(out, "r") as f:
        mp = f["MODEL_STAGE[0]/MODEL/PARTICLES"]
        assert "SOURCE_SET_ID" not in mp.attrs
        assert lf.str_dataset(mp["BODIES/NAME"]) == ["body7"]
        steps = list(f["MODEL_STAGE[0]/RESULTS/ON_PARTICLES/STEP"][...])
        assert len(steps) == 5            # t = 0, 3, 6, 9 x 1e-4 and the last step


def test_classic_tcl_round_trip(tmp_path):
    exe = _tcl_exe()
    if exe is None:
        pytest.skip("classic OpenSees exe not found (build it, or set LADRUNO_TCL_EXE)")
    out = tmp_path / "tcl.ladruno"
    deck = tmp_path / "c2.tcl"
    deck.write_text(
        "model basic -ndm 2 -ndf 2\n"
        "nDMaterial ElasticIsotropic 1 2.0e5 0.3\n"
        "nDMaterial ElasticIsotropic 2 2.0e4 0.3\n"
        f"mpm grid {TAG} -origin 0 0 -h 0.5 -cells 4 4 -basis gimp\n"
        f"puts \"LOADED [mpm particles -file {{{FIXTURE.as_posix()}}}]\"\n"
        f"mpm output {{{out.as_posix()}}} -every 2 -fields {' '.join(FIELDS)}\n"
        "mpm analyze 5 1.0e-4\n"
        "mpm output -close\n"
        "puts DECK-DONE\n", encoding="utf-8")
    proc = subprocess.run([exe, str(deck)], capture_output=True, text=True, timeout=300,
                          env=dict(os.environ, LADRUNO_OPENSEES_QUIET="1"), cwd=str(tmp_path),
                          stdin=subprocess.DEVNULL)
    log = proc.stdout + proc.stderr
    assert "DECK-DONE" in log and "LOADED 12" in log, log
    assert mc.validate_particle_results(str(out)) == []
    with h5py.File(out, "r") as f:
        assert lf._attr(f["MODEL_STAGE[0]/MODEL/PARTICLES"], "SOURCE_SET_ID") == _set()["SET_ID"]


# --------------------------------------------------------------------------- #
#  The c2 review (Fable, ba7ab310e..742b9981c): writer lifecycle               #
# --------------------------------------------------------------------------- #
def _counts(path, stage="MODEL_STAGE[0]"):
    with h5py.File(path, "r") as f:
        op = f[f"{stage}/RESULTS/ON_PARTICLES"]
        return np.diff(op["OFFSETS"][...]).tolist(), op


def test_particles_added_after_a_refused_analyze_reach_the_file(tmp_path):
    """Review 1: the stage opens only once a step passes its refusals."""
    ops = _ops()
    _deck(ops)
    ops.mpm("body", 1, "box", 0.5, 0.5, 1.0, 1.0, "-ppc", 2, "-mat", 1, "-rho", 2.0)
    out = tmp_path / "late.ladruno"
    ops.mpm("output", str(out), "-every", 1)
    dtc = float(ops.mpm("criticalTimeStep"))
    with pytest.raises(Exception):
        ops.mpm("analyze", 1, 10 * dtc)
    ops.mpm("body", 2, "box", 1.0, 1.0, 1.5, 1.5, "-ppc", 2, "-mat", 2, "-rho", 2.0)
    ops.mpm("analyze", 3, 0.2 * float(ops.mpm("criticalTimeStep")))
    n = int(ops.mpm("info", "numParticles"))
    ops.mpm("output", "-close")
    assert mc.validate_particle_results(str(out)) == []
    counts, _ = _counts(out)
    assert counts == [n] * 4 and n == 8


def test_output0_holds_initial_conditions_set_after_a_zero_step_analyze(tmp_path):
    """Review 1: `analyze 0` opens no stage, so a later velocity is in output 0."""
    ops = _ops()
    _deck(ops)
    ops.mpm("body", 1, "box", 0.5, 0.5, 1.0, 1.0, "-ppc", 2, "-mat", 1, "-rho", 2.0)
    out = tmp_path / "zero.ladruno"
    ops.mpm("output", str(out), "-every", 1)
    ops.mpm("analyze", 0, 1e-5)
    ops.mpm("velocity", 1, "uniform", 0.3, 0.0)
    ops.mpm("analyze", 2, 1e-5)
    ops.mpm("output", "-close")
    with h5py.File(out, "r") as f:
        op = f["MODEL_STAGE[0]/RESULTS/ON_PARTICLES"]
        a, b = op["OFFSETS"][0], op["OFFSETS"][1]
        assert np.allclose(op["VELOCITY/DATA"][a:b], [0.3, 0.0])


def test_a_failed_analyze_keeps_its_accepted_steps(tmp_path):
    """Review 7: a block leaving the grid is refused mid-run; the steps it took are
    written (the last one too), and the file validates."""
    ops = _ops()
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e6, 0.25)
    ops.mpm("grid", TAG, "-origin", 0.0, 0.0, "-h", 0.1, "-cells", 10, 10)
    ops.mpm("body", 1, "box", 0.1, 0.4, 0.3, 0.6, "-ppc", 2, "-mat", 1, "-rho", 2000.0)
    ops.mpm("velocity", 1, "uniform", -20.0, 0.0)
    out = tmp_path / "fail.ladruno"
    ops.mpm("output", str(out), "-every", 4)
    dt = 0.3 * float(ops.mpm("criticalTimeStep"))
    with pytest.raises(Exception):
        ops.mpm("analyze", 400, dt)
    t_end = float(ops.mpm("info", "time"))
    ops.mpm("output", "-close")
    assert mc.validate_particle_results(str(out)) == []
    with h5py.File(out, "r") as f:
        times = f["MODEL_STAGE[0]/RESULTS/ON_PARTICLES/TIME"][...]
        assert times[-1] == t_end and len(times) >= 2


def test_no_file_without_a_stage(tmp_path):
    """Nit: `mpm output` then `-close` with no analyze writes no INFO-only file."""
    ops = _ops()
    _deck(ops)
    ops.mpm("particles", "-file", str(FIXTURE))
    out = tmp_path / "never.ladruno"
    ops.mpm("output", str(out))
    ops.mpm("output", "-close")
    assert not out.exists()


def test_dt_cadence_rearms_when_time_goes_back(tmp_path):
    """Nit: after the Domain's `reset` (time back to 0, same epoch) the -dt cadence
    keeps writing."""
    ops = _ops()
    _deck(ops)
    ops.mpm("particles", "-file", str(FIXTURE))
    out = tmp_path / "rearm.ladruno"
    ops.mpm("output", str(out), "-dt", 2.0e-4)
    ops.mpm("analyze", 6, 1.0e-4)
    ops.reset()
    ops.mpm("analyze", 6, 1.0e-4)
    ops.mpm("output", "-close")
    assert lf.validate(str(out)) == []
    with h5py.File(out, "r") as f:
        times = f["MODEL_STAGE[0]/RESULTS/ON_PARTICLES/TIME"][...]
    after = [t for k, t in enumerate(times) if k > 0 and times[k] < times[k - 1]]
    assert after, times                      # the time went back inside the stage
    restart = int(np.argmax(np.diff(times) < 0)) + 1
    assert len(times) - restart >= 3, times   # 2e-4, 4e-4, 6e-4 after the restart


_KILL_CHILD = r"""
import os, sys
sys.path.insert(0, os.environ["LADRUNO_MPM_TESTS"])
from _testbed import ops
out, fixture = sys.argv[1], sys.argv[2]
ops.wipe(); ops.model("basic", "-ndm", 2, "-ndf", 2)
ops.nDMaterial("ElasticIsotropic", 1, 2.0e5, 0.3)
ops.nDMaterial("ElasticIsotropic", 2, 2.0e4, 0.3)
ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 0.5, "-cells", 4, 4, "-basis", "gimp")
ops.mpm("particles", "-file", fixture)
ops.mpm("output", out, "-every", 1, "-fields", "POSITION", "VELOCITY", "STRESS", "VOLUME")
print("CHILD START", flush=True)
ops.mpm("analyze", 50000000, 1e-6)
"""


def _killed_run(tmp_path, wait):
    """Run the kill child, kill it `wait` s after its analyze starts, return the file."""
    import time
    _ops()
    out = tmp_path / "killed.ladruno"
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1",
               LADRUNO_MPM_TESTS=os.path.dirname(os.path.abspath(__file__)))
    env.pop("HDF5_USE_FILE_LOCKING", None)          # the reader rule needs no lock override
    p = subprocess.Popen([sys.executable, "-c", _KILL_CHILD, str(out), str(FIXTURE)],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
                         stdin=subprocess.DEVNULL, env=env)
    try:
        assert "CHILD START" in p.stdout.readline()
        time.sleep(wait)
    finally:
        p.kill()
        p.wait()
    return out


@pytest.mark.parametrize("wait", [0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0])
def test_a_killed_run_leaves_a_readable_file(tmp_path, wait):
    """Review 2 (and its follow-up ruling): under the DEFAULT settings (a flush after
    every output) a run killed mid-analyze at any point leaves a file that opens by the
    reader rule (`open_particle_file`: SWMR-read after a kill), and the contract's crash
    rule (232a §2.5) recovers every completed output, data included. STRICT: every
    output the crash rule counts is exact. The writer's barrier (a flush after the rows
    and before the OFFSETS entry) makes a counted entry always have its rows; SWMR
    makes the chunk-index updates crash-safe (an index split no longer tears)."""
    out = _killed_run(tmp_path, wait)
    S = _set()
    with mc.open_particle_file(out) as f:                # must open
        T = mc.committed_outputs(f["MODEL_STAGE[0]"])     # the crash rule
        op = f["MODEL_STAGE[0]/RESULTS/ON_PARTICLES"]
        steps, offs = op["STEP"][:T], op["OFFSETS"][:T + 1]
        assert T >= 1, (T, op["STEP"].shape, op["OFFSETS"].shape)
        ids = op["ID"][:offs[T]]
        for k in ("POSITION", "VELOCITY", "STRESS", "VOLUME"):
            data = op[k]["DATA"][:offs[T]]                 # readable, every completed output
            assert np.all(np.isfinite(data)), k
        assert np.all(np.diff(offs) >= 0)
        vol = op["VOLUME"]["DATA"][:offs[T], 0]
        for k in range(T):                                 # every COUNTED output, the last too
            assert np.array_equal(ids[offs[k]:offs[k + 1]], S["ID"]), (k, T)
            assert np.all(vol[offs[k]:offs[k + 1]] > 0.0), (k, T)   # fill would be 0
        assert np.all(np.diff(steps) > 0)


def test_a_killed_file_needs_the_swmr_read_open(tmp_path):
    """The reader rule (232a §2.5): a killed run leaves the SWMR consistency flags set,
    so a plain read-only open is refused ("already open for write"); the rule's
    SWMR-read open reads it, with no HDF5_USE_FILE_LOCKING override."""
    out = _killed_run(tmp_path, 1.0)
    with pytest.raises(OSError, match="already open for write"):
        h5py.File(out, "r")
    with mc.open_particle_file(out) as f:
        assert f.swmr_mode
        assert f["MODEL_STAGE[0]/RESULTS/ON_PARTICLES/OFFSETS"].shape[0] >= 2
    assert mc.committed_outputs(out) >= 1


def test_a_closed_file_opens_plainly_in_the_latest_format(tmp_path):
    """A clean close clears the SWMR flags: a plain read-only open works (no SWMR
    mode). The file uses the latest HDF5 format (superblock 3), as SWMR requires."""
    ops = _ops()
    _deck(ops)
    ops.mpm("particles", "-file", str(FIXTURE))
    out = tmp_path / "closed.ladruno"
    _run_and_write(ops, out, steps=5, every=2)
    with h5py.File(out, "r") as f:
        assert not f.swmr_mode
        assert f.id.get_create_plist().get_version()[0] >= 3
    with mc.open_particle_file(out) as f:
        assert not f.swmr_mode


# --------------------------------------------------------------------------- #
#  WP-232h (#1085): the presized root and the guard lock                       #
# --------------------------------------------------------------------------- #
_RESET_CHILD = r"""
import os, sys
sys.path.insert(0, os.environ["LADRUNO_MPM_TESTS"])
from _testbed import ops
out, fixture, stages, steps = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
ops.wipe(); ops.model("basic", "-ndm", 2, "-ndf", 2)
ops.nDMaterial("ElasticIsotropic", 1, 2.0e5, 0.3)
ops.nDMaterial("ElasticIsotropic", 2, 2.0e4, 0.3)
ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 0.5, "-cells", 4, 4, "-basis", "gimp")
ops.mpm("particles", "-file", fixture)
ops.mpm("output", out, "-every", 1, "-fields", "POSITION", "VELOCITY", "STRESS", "VOLUME")
for s in range(stages):
    ops.mpm("analyze", steps, 1e-6)
    if s == 1:
        print("STAGE 1", flush=True)
    if s < stages - 1:
        ops.mpm("reset")
ops.mpm("output", "-close")
print("DONE", flush=True)
"""


def _reset_child(out, stages, steps, stderr=subprocess.PIPE):
    _ops()
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1",
               LADRUNO_MPM_TESTS=os.path.dirname(os.path.abspath(__file__)))
    env.pop("HDF5_USE_FILE_LOCKING", None)
    return subprocess.Popen([sys.executable, "-c", _RESET_CHILD, str(out), str(FIXTURE),
                             str(stages), str(steps)],
                            stdout=subprocess.PIPE, stderr=stderr, text=True,
                            stdin=subprocess.DEVNULL, env=env)


def _h5lib():
    import ctypes
    lib = ctypes.CDLL(os.path.join(os.path.dirname(h5py.__file__), "hdf5.dll")) \
        if os.name == "nt" else ctypes.CDLL(h5py.h5.__file__)
    return ctypes, lib


def test_the_root_is_presized_and_stays_compact(tmp_path):
    """#1085: the root gains a MODEL_STAGE link at every `mpm reset`, in a layout written
    without SWMR. Its header is presized (``plist::presizeLinks``: links compact for good,
    sized for 1024) so adding one is ONE in-place write of the header. Unsized, it goes
    dense after 8 links (here 13: INFO + 12 stages) or grows a continuation chunk, and a
    kill in that flush could make the root, and every earlier stage, unreadable
    (``hdf5_crash_replay.py --resets 10 --no-presize``)."""
    out = tmp_path / "presized.ladruno"
    o, e = _reset_child(out, 12, 2).communicate(timeout=600)
    assert "DONE" in o, e[-2000:]
    assert lf.validate(str(out)) == []
    ctypes, lib = _h5lib()
    with h5py.File(out, "r") as f:
        assert sum(1 for k in f if k.startswith("MODEL_STAGE")) == 12
        gcpl = f["/"].id.get_create_plist()     # the root group's (stored) link info
        mx, mn = ctypes.c_uint(), ctypes.c_uint()
        lib.H5Pget_link_phase_change.argtypes = (ctypes.c_int64, ctypes.c_void_p, ctypes.c_void_p)
        assert lib.H5Pget_link_phase_change(gcpl.id, ctypes.byref(mx), ctypes.byref(mn)) >= 0
        assert (mx.value, mn.value) == (65535, 65534)

        class GInfo(ctypes.Structure):    # H5G_info_t
            _fields_ = [("storage_type", ctypes.c_int), ("nlinks", ctypes.c_uint64),
                        ("max_corder", ctypes.c_int64), ("mounted", ctypes.c_bool)]
        gi = GInfo()
        lib.H5Gget_info.argtypes = (ctypes.c_int64, ctypes.c_void_p)
        assert lib.H5Gget_info(f.id.id, ctypes.byref(gi)) >= 0
        assert gi.nlinks == 13
        assert gi.storage_type == 1          # H5G_STORAGE_TYPE_COMPACT (2 = dense)
        # sized up front (H5Pset_est_link_info): the 13 links sit in ONE header chunk; an
        # unsized compact root would have grown continuation chunks instead
        assert h5py.h5o.get_info(f.id).hdr.nchunks == 1


@pytest.mark.parametrize("locking", [True, False])
def test_a_live_reader_does_not_disturb_the_writer(tmp_path, locking):
    """#1085: readers poll a live particle file through 15 stages (an `mpm reset` every
    200 steps: ~1 s of writing, so the poll overlaps it on a slow runner too). locking=True (h5py's default): HDF5's own locking would let the reader's
    shared lock break the SWMR writer (Linux: H5Fstart_swmr_write refused, EAGAIN;
    Windows: the writes themselves fail). The writer turns HDF5 locking off and holds the
    guard lock for the file's life, so the reader is refused cleanly
    (``WriterBusyError``). locking=False (the live monitor's way): no lock, the reader
    reads through SWMR. Either way the writer reports no error and the finished file
    equals a run nobody read."""
    import time
    out = tmp_path / "live.ladruno"
    err = tmp_path / "live.err"
    with open(err, "w") as fe:
        p = _reset_child(out, 15, 200, stderr=fe)
        assert "STAGE 1" in p.stdout.readline()
        outcomes = {"open": 0, "busy": 0, "other": 0}
        deadline = time.monotonic() + 600.0
        while p.poll() is None:
            if time.monotonic() > deadline:
                p.kill()
                pytest.fail(f"the writer did not finish in 600 s; reader outcomes {outcomes}")
            try:
                f = (lf.open_ladruno(str(out)) if locking
                     else h5py.File(out, "r", swmr=True, locking=False))
            except lf.WriterBusyError:
                outcomes["busy"] += 1
                continue
            except OSError:
                outcomes["other"] += 1        # a stage being laid out: retry
                continue
            with f:
                outcomes["open"] += 1
                time.sleep(0.002)             # hold it a moment, as a viewer would
        o, _ = p.communicate()
    e = err.read_text()
    assert "DONE" in o, e[-3000:]
    assert "ERROR" not in e and "Permission denied" not in e, e[-3000:]
    if locking:
        assert outcomes["busy"] > 0, outcomes     # refused, never let in to break the writer
    else:
        assert outcomes["open"] > 0, outcomes
    ref = tmp_path / "ref.ladruno"
    o, e2 = _reset_child(ref, 15, 200).communicate(timeout=600)
    assert "DONE" in o, e2[-2000:]
    assert lf.validate(str(out)) == []
    with h5py.File(out, "r") as f, h5py.File(ref, "r") as r:
        names = sorted(k for k in r if k.startswith("MODEL_STAGE"))
        assert names == sorted(k for k in f if k.startswith("MODEL_STAGE"))
        assert len(names) == 15
        for n in names:
            op, rp = f[f"{n}/RESULTS/ON_PARTICLES"], r[f"{n}/RESULTS/ON_PARTICLES"]
            for ds in ("STEP", "TIME", "OFFSETS", "ID", "POSITION/DATA", "VELOCITY/DATA",
                       "STRESS/DATA", "VOLUME/DATA"):
                assert np.array_equal(op[ds][...], rp[ds][...]), (n, ds)
