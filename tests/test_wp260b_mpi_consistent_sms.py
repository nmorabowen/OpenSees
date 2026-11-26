"""WP-260b Finding 3 — OpenSeesMP `MPIDiagonalSOE` consistent-SMS gate (`np=3`).

Exercises the distributed-PCG paths touched by WP-260b under real MPI:
  - packed `localNeedW + 2^24 * localBad` decode in `LadrunoConsistentRefine.h`;
  - lockstep `gNeedW > 0` multiplicity-weight build on `MPIDiagonalSOE`;
  - `LinearSOE::globalReduceSum2` / `MPIDiagonalSOE::globalReduceSum2` 2-double
    `MPI_Allreduce`;
  - preallocated `sharedSumRequests` / `sharedSumStatuses` across 3 ranks (rank 1
    communicates with BOTH rank 0 and rank 2);
  - reordered parallel PCG (`consistentParPCG`) with scaled elements straddling
    both partition boundaries (shared nodes 4 and 7).

Known upstream blocker (#1161):
  In-place mid-run `remove element` without `wipeAnalysis` under `system MPIDiagonal`
  crashes (`0xc0000005`) or zeroes `A` even on baseline `origin/ladruno` (`7e91aaf1c`)
  with plain `CentralDifferenceLadruno`, because upstream `MPIDiagonalSOE::setSize()`
  never resets `isAfactored`, `MPIDiagonalSolver::setSize()` never resets `notSet`,
  and `MPIDiagonalSolver::solve()` caches `static int size`/`numShared`.
  Accordingly, this module tests:
    (1) steady 35-step `np=3` vs serial (`CentralDifferenceSMSConsistent` and
        `ExplicitBathe -sms -consistent`, scaled elements on both partition boundaries);
    (2) mid-run `remove element` across `wipeAnalysis` + fresh `MPIDiagonal` (same
        `neq`, comparing pre- and post-removal trajectories against serial).

Where it runs:
  - Runs automatically on the 5-target Windows build (`dist/bin/OpenSeesMP.exe` +
    `dist/openseesmp/mpiexec.exe`) and any Linux/esmeralda MP build via
    `LADRUNO_MP_EXE` / `LADRUNO_MPIEXEC`.
  - Zone-A Ubuntu CI builds only serial `OpenSees` + `OpenSeesPy` (no MPI target),
    so `@needs_mp` skips cleanly there while keeping the test committed in the repo.

Run:  py -3.12 -m pytest tests/test_wp260b_mpi_consistent_sms.py -v
"""
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]


def _find_binaries():
    mp_exe = os.environ.get("LADRUNO_MP_EXE")
    mpiexec = os.environ.get("LADRUNO_MPIEXEC")
    ser_exe = os.environ.get("LADRUNO_SERIAL_EXE")
    if not mp_exe:
        cand = REPO / "dist" / "bin" / ("OpenSeesMP.exe" if os.name == "nt" else "OpenSeesMP")
        mp_exe = str(cand) if cand.exists() else None
    if not mpiexec:
        cand = REPO / "dist" / "openseesmp" / ("mpiexec.exe" if os.name == "nt" else "mpiexec")
        mpiexec = str(cand) if cand.exists() else shutil.which("mpiexec")
    if not ser_exe:
        cand = REPO / "dist" / "bin" / ("OpenSees.exe" if os.name == "nt" else "OpenSees")
        ser_exe = str(cand) if cand.exists() else None
    return mp_exe, mpiexec, ser_exe


MP_EXE, MPIEXEC, SER_EXE = _find_binaries()
needs_mp = pytest.mark.skipif(
    MP_EXE is None or MPIEXEC is None or SER_EXE is None,
    reason="needs OpenSeesMP, mpiexec, and OpenSees (5-target build or LADRUNO_*_EXE)",
)

_TCL_DECK = r"""
set is_mp %(is_mp)d
set do_midrun_remove %(do_midrun_remove)d
if {$is_mp} {
    set pid [getPID]
    set np  [getNP]
    if {$np != 3} {
        puts stderr "ERROR expected np=3, got $np"
        exit 2
    }
} else {
    set pid -1
    set np  1
}

model BasicBuilder -ndm 2 -ndf 2

set lengths {1.0 1.0 0.05 0.05 1.0 0.05 0.05 1.0 1.0}
set xs {0.0}
set acc 0.0
foreach L $lengths {
    set acc [expr {$acc + $L}]
    lappend xs $acc
}

set E 1.0e4
set rho 1.0
set A 1.0
uniaxialMaterial Elastic 1 $E
uniaxialMaterial Elastic 2 4.0e4

if {$is_mp} {
    set e_start [expr {$pid * 3 + 1}]
    set e_end   [expr {($pid + 1) * 3}]
    set n_start $e_start
    set n_end   [expr {$e_end + 1}]
} else {
    set e_start 1
    set e_end   9
    set n_start 1
    set n_end   10
}

for {set n $n_start} {$n <= $n_end} {incr n} {
    set x [lindex $xs [expr {$n - 1}]]
    node $n $x 0.0
    if {$n == 1} {
        fix $n 1 1
    } else {
        fix $n 0 1
    }
}

for {set e $e_start} {$e <= $e_end} {incr e} {
    element truss $e $e [expr {$e + 1}] $A 1 -rho $rho
}

if {$do_midrun_remove} {
    if {!$is_mp || $pid == 0} { element truss 101 2 3 0.5 2 -rho $rho }
    if {!$is_mp || $pid == 1} { element truss 102 5 6 0.5 2 -rho $rho }
    if {!$is_mp || $pid == 2} { element truss 103 8 9 0.5 2 -rho $rho }
}

if {!$is_mp || $pid == 2} {
    timeSeries Constant 1
    pattern Plain 1 1 {
        load 10 50.0 0.0
    }
}

proc setup_analysis {is_mp} {
    constraints Transformation
    if {$is_mp} {
        numberer ParallelRCM
        system MPIDiagonal
    } else {
        numberer Plain
        system Diagonal
    }
    test NormDispIncr 1.0e-12 1
    algorithm Linear
    %(integrator_cmd)s
    analysis Transient
}

setup_analysis $is_mp

set dt 0.012
set out_path "%(out_prefix)s"
if {$is_mp} {
    set out_file "${out_path}_rank${pid}.txt"
} else {
    set out_file "${out_path}_serial.txt"
}
set fd [open $out_file w]

for {set step 1} {$step <= 35} {incr step} {
    if {$do_midrun_remove && $step == 16} {
        wipeAnalysis
        if {!$is_mp || $pid == 0} { remove element 101 }
        if {!$is_mp || $pid == 1} { remove element 102 }
        if {!$is_mp || $pid == 2} { remove element 103 }
        setup_analysis $is_mp
    }
    set ok [analyze 1 $dt]
    if {$ok != 0} {
        puts stderr "ERROR analyze failed at step $step (ok=$ok, pid=$pid)"
        close $fd
        exit 3
    }
    set line "$step"
    for {set n $n_start} {$n <= $n_end} {incr n} {
        append line [format " %%d:%%.17g" $n [nodeDisp $n 1]]
    }
    puts $fd $line
}
close $fd
puts "WP260B-MPI-DONE pid=$pid"
"""


def _parse_trace(path: Path):
    out = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.strip().split()
        if not parts:
            continue
        step = int(parts[0])
        vals = {}
        for tok in parts[1:]:
            nt, u = tok.split(":")
            vals[int(nt)] = float(u)
        out[step] = vals
    return out


def _env_with_dist():
    env = dict(os.environ)
    env["LADRUNO_OPENSEES_QUIET"] = "1"
    paths = [str(REPO / "dist" / "openseesmp"), str(REPO / "dist" / "bin")]
    if env.get("PATH"):
        paths.append(env["PATH"])
    env["PATH"] = os.pathsep.join(paths)
    return env


def _run_and_compare(tmp_path, integrator_kind, integrator_cmd, do_midrun_remove):
    tag = f"{integrator_kind}_{'midrem' if do_midrun_remove else 'steady'}"
    out_prefix = (tmp_path / f"sms_{tag}").as_posix()
    ser_tcl = tmp_path / f"ser_{tag}.tcl"
    mp_tcl = tmp_path / f"mp_{tag}.tcl"

    ser_tcl.write_text(
        _TCL_DECK % {
            "is_mp": 0,
            "do_midrun_remove": int(do_midrun_remove),
            "integrator_cmd": integrator_cmd,
            "out_prefix": out_prefix,
        },
        encoding="utf-8",
    )
    mp_tcl.write_text(
        _TCL_DECK % {
            "is_mp": 1,
            "do_midrun_remove": int(do_midrun_remove),
            "integrator_cmd": integrator_cmd,
            "out_prefix": out_prefix,
        },
        encoding="utf-8",
    )

    env = _env_with_dist()

    p_ser = subprocess.run(
        [SER_EXE, str(ser_tcl)],
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        timeout=120,
        env=env,
        encoding="utf-8",
        errors="replace",
    )
    ser_out = p_ser.stdout + p_ser.stderr
    assert p_ser.returncode == 0 and "WP260B-MPI-DONE" in ser_out, (
        f"serial run failed (rc={p_ser.returncode})\nSTDOUT:\n{p_ser.stdout}\nSTDERR:\n{p_ser.stderr}"
    )

    try:
        p_mp = subprocess.run(
            [MPIEXEC, "-n", "3", MP_EXE, str(mp_tcl)],
            stdin=subprocess.DEVNULL,
            capture_output=True,
            text=True,
            timeout=120,
            env=env,
            encoding="utf-8",
            errors="replace",
        )
    except subprocess.TimeoutExpired as exc:
        pytest.fail(f"OpenSeesMP np=3 hung (collective/lockstep deadlock?): {exc}")

    mp_out = p_mp.stdout + p_mp.stderr
    assert p_mp.returncode == 0 and mp_out.count("WP260B-MPI-DONE") == 3, (
        f"MP np=3 run failed (rc={p_mp.returncode})\nSTDOUT:\n{p_mp.stdout}\nSTDERR:\n{p_mp.stderr}"
    )

    ser = _parse_trace(Path(f"{out_prefix}_serial.txt"))
    ranks = [_parse_trace(Path(f"{out_prefix}_rank{r}.txt")) for r in range(3)]
    assert len(ser) == 35
    for r in range(3):
        assert len(ranks[r]) == 35

    max_rel_pre = 0.0
    max_rel_post = 0.0
    max_shared_diff = 0.0
    for step in range(1, 36):
        ser_vec = np.array([ser[step][n] for n in range(1, 11)])
        scale = max(float(np.linalg.norm(ser_vec)), 1e-30)

        max_shared_diff = max(
            max_shared_diff,
            abs(ranks[0][step][4] - ranks[1][step][4]) / scale,
            abs(ranks[1][step][7] - ranks[2][step][7]) / scale,
        )

        mp_nodes = {}
        for r in range(3):
            mp_nodes.update(ranks[r][step])
        mp_vec = np.array([mp_nodes[n] for n in range(1, 11)])
        rel = float(np.linalg.norm(mp_vec - ser_vec) / scale)
        if step <= 15:
            max_rel_pre = max(max_rel_pre, rel)
        else:
            max_rel_post = max(max_rel_post, rel)

    assert max_shared_diff < 1e-13, f"shared boundary node mismatch: {max_shared_diff:.3e}"
    assert max_rel_pre < 1e-11, f"pre-step16 np=3 vs serial rel err: {max_rel_pre:.3e}"
    assert max_rel_post < 1e-11, f"post-step16 np=3 vs serial rel err: {max_rel_post:.3e}"
    tip_final = abs(ser[35][10])
    assert 1e-4 < tip_final < 1.0, f"unexpected tip disp {tip_final}"


@needs_mp
@pytest.mark.parametrize(
    "integrator_kind,integrator_cmd",
    [
        (
            "cdl_sms_cons",
            "integrator CentralDifferenceSMSConsistent 0.012 -maxAddedMass 50.0 -pcgTol 1.0e-13 -pcgMaxIt 500",
        ),
        (
            "eb_sms_cons",
            "integrator ExplicitBathe 0.54 -sms 0.012 -consistent -maxAddedMass 50.0 -pcgTol 1.0e-13 -pcgMaxIt 500",
        ),
    ],
)
def test_mpidiagonal_consistent_sms_np3_matches_serial(tmp_path, integrator_kind, integrator_cmd):
    _run_and_compare(tmp_path, integrator_kind, integrator_cmd, do_midrun_remove=False)


@needs_mp
@pytest.mark.parametrize(
    "integrator_kind,integrator_cmd",
    [
        (
            "cdl_sms_cons",
            "integrator CentralDifferenceSMSConsistent 0.012 -maxAddedMass 50.0 -pcgTol 1.0e-13 -pcgMaxIt 500",
        ),
        (
            "eb_sms_cons",
            "integrator ExplicitBathe 0.54 -sms 0.012 -consistent -maxAddedMass 50.0 -pcgTol 1.0e-13 -pcgMaxIt 500",
        ),
    ],
)
def test_mpidiagonal_consistent_sms_np3_midrun_remove_after_wipe_analysis(
    tmp_path, integrator_kind, integrator_cmd
):
    _run_and_compare(tmp_path, integrator_kind, integrator_cmd, do_midrun_remove=True)
