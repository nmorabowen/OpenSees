"""Rank-env probe gate -- model runner.

The recorder's partitioned-output detection (LadrunoRecorder::initialize)
indexes interpreter-per-rank runs by the launcher's per-rank environment.
Launchers disagree on the variable names (PMI_* for Intel/MS-MPI, OMPI_* for
OpenMPI, SLURM_* for srun under any plugin -- including pmix, where no PMI_*
exist; that gap silently killed a 16-rank cluster run). This gate fakes each
launcher's environment around an otherwise-sequential run and asserts the
output filename reacts correctly (the checker then validates the manifest).

The pyd links a static CRT that snapshots the environment at DLL load, so a
mutated os.environ is invisible to std::getenv after import -- each case
therefore runs in its own child process with the env set at launch, exactly
like a real MPI launcher does.

Run (matches run_regression.bat :gate1):
    python rank_env_model.py <dist_bin> <out_dir>
"""

from __future__ import annotations

import os
import subprocess
import sys

DIST = sys.argv[1]
OUT = sys.argv[2]

ENV_VARS = (
    "PMI_SIZE", "PMI_RANK",
    "OMPI_COMM_WORLD_SIZE", "OMPI_COMM_WORLD_RANK",
    "SLURM_NTASKS", "SLURM_PROCID",
)

# (case name, env to set, expected output filename)
CASES = [
    ("plain", {}, "rank_env_plain.ladruno"),
    ("pmi", {"PMI_SIZE": "4", "PMI_RANK": "2"},
     "rank_env_pmi.part-2.ladruno"),
    ("ompi", {"OMPI_COMM_WORLD_SIZE": "8", "OMPI_COMM_WORLD_RANK": "5"},
     "rank_env_ompi.part-5.ladruno"),
    ("slurm", {"SLURM_NTASKS": "16", "SLURM_PROCID": "7"},
     "rank_env_slurm.part-7.ladruno"),
    # a real MPI launcher's own pair outranks sbatch's ambient SLURM_NTASKS
    ("precedence",
     {"PMI_SIZE": "4", "PMI_RANK": "2",
      "SLURM_NTASKS": "16", "SLURM_PROCID": "7"},
     "rank_env_precedence.part-2.ladruno"),
]


def child(case: str) -> None:
    """Run one tiny model with the recorder; env was set by the parent."""
    os.add_dll_directory(DIST)
    sys.path.insert(0, DIST)
    import opensees as ops

    declared = os.path.join(OUT, f"rank_env_{case}.ladruno")
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 1.0, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.uniaxialMaterial("Elastic", 1, 1000.0)
    ops.element("truss", 1, 1, 2, 1.0, 1)
    ops.recorder("ladruno", declared)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0, 0.0)
    ops.system("BandSPD")
    ops.numberer("RCM")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    ops.analyze(1)
    ops.wipe()  # destroys the recorder -> flush + close


if len(sys.argv) > 4 and sys.argv[3] == "--child":
    child(sys.argv[4])
    sys.exit(0)

ok = True
for case, env, expected in CASES:
    for fname in (f"rank_env_{case}.ladruno", expected):
        stale = os.path.join(OUT, fname)
        if os.path.exists(stale):
            os.remove(stale)
    child_env = {k: v for k, v in os.environ.items() if k not in ENV_VARS}
    child_env.update(env)
    child_env["LADRUNO_OPENSEES_QUIET"] = "1"
    proc = subprocess.run(
        [sys.executable, os.path.abspath(__file__), DIST, OUT, "--child", case],
        env=child_env, capture_output=True, text=True)
    if proc.returncode != 0:
        print(f" FAIL {case}: child rc={proc.returncode}\n{proc.stderr[-2000:]}")
        ok = False
        continue
    if not os.path.exists(os.path.join(OUT, expected)):
        print(f" FAIL {case}: expected {expected}, dir has "
              f"{sorted(f for f in os.listdir(OUT) if f.startswith('rank_env_'))}")
        ok = False
        continue
    print(f"  ok  {case}: wrote {expected}")

if not ok:
    sys.exit(1)
print("RANK_ENV_MODEL OK")
