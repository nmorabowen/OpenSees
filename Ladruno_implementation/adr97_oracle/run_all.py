"""ADR-97 P0 -- run every oracle and write ``reference_output.txt``.

    python3.12 Ladruno_implementation/adr97_oracle/run_all.py

The first four oracles are pure numpy.  ``fd_tangent_driver.py`` needs the built
OpenSeesPy module; give it a PYTHONPATH, e.g.::

    LADRUNO_OPENSEES_QUIET=1 \\
    PYTHONPATH=<worktree>/dist/bin \\
    python3.12 Ladruno_implementation/adr97_oracle/run_all.py

Stderr is DISCARDED: the C++ side writes a page of diagnostics per
``ASDPlasticMaterial3D`` construction on ``opserr``, which is noise here (the
Python-side stdout of every script is kept verbatim).  Any non-zero exit status
is recorded and re-raised at the end -- every oracle asserts its own gates, so a
clean run is itself the pass/fail signal.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ["cppm_vm.py", "cppm_dp.py", "cppm_mc.py", "path_independence.py",
           "fd_tangent_driver.py"]
OUT = os.path.join(HERE, "reference_output.txt")


def _ops_available():
    r = subprocess.run([sys.executable, "-c",
                        "import opensees" if True else ""],
                       capture_output=True)
    if r.returncode == 0:
        return True
    r = subprocess.run([sys.executable, "-c", "import openseespy.opensees"],
                       capture_output=True)
    return r.returncode == 0


def main():
    env = dict(os.environ)
    env.setdefault("LADRUNO_OPENSEES_QUIET", "1")
    chunks, failures = [], []
    for name in SCRIPTS:
        path = os.path.join(HERE, name)
        if name == "fd_tangent_driver.py" and not _ops_available():
            chunks.append(
                f"\n\n##### {name} #####\nSKIPPED -- no OpenSeesPy module on "
                f"PYTHONPATH.\nRe-run with e.g.\n  LADRUNO_OPENSEES_QUIET=1 "
                f"PYTHONPATH=<worktree>/dist/bin python3.12 "
                f"{os.path.join('Ladruno_implementation', 'adr97_oracle', name)}\n")
            print(f"[skip] {name} (no opensees module)")
            continue
        print(f"[run ] {name}")
        r = subprocess.run([sys.executable, path], capture_output=True,
                           text=True, env=env, cwd=HERE)
        chunks.append(f"\n\n##### {name}  (exit {r.returncode}) #####\n"
                      + r.stdout)
        if r.returncode != 0:
            failures.append(name)
            chunks.append("\n--- STDERR TAIL ---\n"
                          + "\n".join(r.stderr.strip().split("\n")[-25:]))
            print(f"[FAIL] {name}")
    header = ("ADR-97 P0 oracle reference output\n"
              "Regenerate with:  python3.12 "
              "Ladruno_implementation/adr97_oracle/run_all.py\n"
              "Every gate in these scripts is an assert; a clean exit IS the "
              "pass.\n" + "=" * 78)
    with open(OUT, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(header + "".join(chunks) + "\n")
    print(f"\nwrote {OUT}")
    if failures:
        raise SystemExit("FAILED: " + ", ".join(failures))


if __name__ == "__main__":
    main()
