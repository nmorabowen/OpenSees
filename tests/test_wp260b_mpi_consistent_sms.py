"""WP-260b Finding 3 — OpenSeesMP `MPIDiagonalSOE` consistent-SMS gate (`np=3`).

Exercises every distributed-PCG path touched by WP-260b under real MPI:
  - packed `localNeedW + 2^24 * localBad` decode in `LadrunoConsistentRefine.h`;
  - lockstep `gNeedW > 0` multiplicity-weight rebuild both at step 1 AND after a
    mid-run topology change (`remove element` on all ranks + `remove node` on rank 2
    so `neq` shrinks and `ConsistentPCGScratch` reallocates);
  - `LinearSOE::globalReduceSum2` / `MPIDiagonalSOE::globalReduceSum2` 2-double
    `MPI_Allreduce`;
  - preallocated `sharedSumRequests` / `sharedSumStatuses` across 3 ranks (rank 1
    communicates with BOTH rank 0 and rank 2);
  - reordered parallel PCG (`consistentParPCG`) with scaled elements straddling
    both partition boundaries (shared nodes 4 and 7).

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

# 9-element 2D axial truss chain (nodes 1..10, y fixed, node 1 clamped) partitioned
# across 3 ranks (3 elements per rank):
#   rank 0: elems 1(L=1.0), 2(L=1.0), 3(L=0.05) -> nodes 1..4  + redundant scaled elem 101 (2-3)
#   rank 1: elems 4(L=0.05), 5(L=1.0), 6(L=0.05) -> nodes 4..7 + redundant scaled elem 102 (5-6)
#   rank 2: elems 7(L=0.05), 8(L=1.0), 9(L=1.0) -> nodes 7..10 + redundant scaled elem 103 (9-100) + free node 100
# Short elements 3,4 straddle the rank 0|1 boundary (node 4); short elements 6,7
# straddle the rank 1|2 boundary (node 7). Every element has dt_e <= 0.01 < dtTarget=0.012.
# After step 15, each rank removes its redundant element (and rank 2 removes free node 100),
# forcing domainChanged -> sc.invalidate() -> gNeedW > 0 rebuild AFTER step 1!
_TCL_DECK = r"""
set is_mp %(is_mp)d
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

# Per-rank redundant scaled elements removed at step 15 (triggers mid-run gNeedW > 0)
if {!$is_mp || $pid == 0} {
    element truss 101 2 3 0.5 2 -rho $rho
}
if {!$is_mp || $pid == 1} {
    element truss 102 5 6 0.5 2 -rho $rho
}
if {!$is_mp || $pid == 2} {
    set x9 [lindex $xs 8]
    node 100 [expr {$x9 + 0.05}] 0.0
    fix 100 0 1
    element truss 103 9 100 0.5 2 -rho $rho
}

# External load ONLY on tip node 10 (owned exclusively by rank 2 — no double load)
if {!$is_mp || $pid == 2} {
    timeSeries Constant 1
    pattern Plain 1 1 {
        load 10 50.0 0.0
    }
}

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

set dt 0.012
set out_path "%(out_prefix)s"
if {$is_mp} {
    set out_file "${out_path}_rank${pid}.txt"
} else {
    set out_file "${out_path}_serial.txt"
}
set fd [open $out_file w]

for {set step 1} {$step <= 35} {incr step} {
    if {$step == 16} {
        if {!$is_mp || $pid == 0} { remove element 101 }
        if {!$is_mp || $pid == 1} { remove element 102 }
        if {!$is_mp || $pid == 2} {
            remove element 103
            remove node 100
        }
    }
    set ok [analyze 1 $dt]
    if {$ok != 0} {
        puts stderr "ERROR analyze failed at step $step (ok=$ok, pid=$pid)"
        close $fd
        exit 3
    }
    set line "$step"
    for {set n $n_start} {$n <= $n_end} {incr n} {
        append line [format " %d:%.17g" $n [nodeDisp $n 1]]
    }
    puts $fd $line
}
close $fd
puts "WP260B-MPI-DONE pid=$pid"
"""


def _parse_trace(path: Path):
    """Return dict[step -> dict[nodeTag -> ux]]."""
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
def test_mpidiagonal_consistent_sms_np3_matches_serial_with_midrun_topology_change(
    tmp_path, integrator_kind, integrator_cmd
):
    out_prefix = (tmp_path / f"sms_{integrator_kind}").as_posix()
    ser_tcl = tmp_path / f"ser_{integrator_kind}.tcl"
    mp_tcl = tmp_path / f"mp_{integrator_kind}.tcl"

    ser_tcl.write_text(
        _TCL_DECK % {"is_mp": 0, "integrator_cmd": integrator_cmd, "out_prefix": out_prefix},
        encoding="utf-8",
    )
    mp_tcl.write_text(
        _TCL_DECK % {"is_mp": 1, "integrator_cmd": integrator_cmd, "out_prefix": out_prefix},
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
    assert p_ser.returncode == 0 and "WP260B-MPI-DONE" in p_ser.stdout, (
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

    assert p_mp.returncode == 0 and p_mp.stdout.count("WP260B-MPI-DONE") == 3, (
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

        # Shared partition-boundary nodes 4 (ranks 0|1) and 7 (ranks 1|2) must agree across ranks
        max_shared_diff = max(
            max_shared_diff,
            abs(ranks[0][step][4] - ranks[1][step][4]) / scale,
            abs(ranks[1][step][7] - ranks[2][step][7]) / scale,
        )

        # Assemble global MP vector and compare against serial DiagonalSOE
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
    assert max_rel_pre < 1e-11, f"pre-removal np=3 vs serial rel err: {max_rel_pre:.3e}"
    assert max_rel_post < 1e-11, f"post-removal (gNeedW>0) np=3 vs serial rel err: {max_rel_post:.3e}"
    # Non-vacuity: displacement is nonzero and bounded (plain CD diverges at dt=0.012 > 0.0005)
    tip_final = abs(ser[35][10])
    assert 1e-4 < tip_final < 1.0, f"unexpected tip disp {tip_final}"
