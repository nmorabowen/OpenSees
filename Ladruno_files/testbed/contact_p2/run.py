"""ADR-78 P2 gate driver.

    python run.py <build-dir>

<build-dir> holds OpenSees.exe / OpenSeesMP.exe / opensees.pyd -- normally
`build\\build\\Release` under the tree you built (named explicitly: a bare
`OpenSees` on PATH has resolved to another session's build on this machine
more than once).

Cases (decks in this dir + the preserved P0 decks in ../contact_parallel):

  serial-auto      serial.tcl        -kn auto reference       -> runs, values
  mp-auto          mp.tcl (np2)      -kn auto partitioned     -> SAME penalty as
                                     serial, measured through the response
                                     (w5/w11/w15/R parity to ~1e-12): the P2
                                     "-kn auto needs no change" gate.
  serial-soft      serial_soft.tcl   -soft at np=1            -> RUNS (no refusal)
  mp1-soft         serial_soft.tcl   -soft, OpenSeesMP -n 1   -> RUNS (np-keyed,
                                     not build-keyed)
  mp2-soft         mp_soft.tcl (np2) -soft under partitioning -> REFUSED with the
                                     named FATAL; the deck's tail line must not
                                     print. The P2 "refusal fires" gate.
  db-roundtrip     db_roundtrip.py   definitions survive Domain::sendSelf /
                                     recvSelf through `database File` save ->
                                     wipe -> restore (P2b).
"""

import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
P0 = os.path.normpath(os.path.join(HERE, "..", "contact_parallel"))


def oneapi_dirs():
    # See ../auto_penalty_mpi/run.py for why each entry exists (impi.dll,
    # libfabric under opt\mpi, compiler runtime). Missing libfabric dies in
    # MPIDI_OFI_mpi_init_hook with an error that names nothing.
    root = os.environ.get("I_MPI_ROOT")
    oneapi = os.environ.get("ONEAPI_ROOT", r"C:\Program Files (x86)\Intel\oneAPI")
    mpiroot = root or os.path.join(oneapi, "mpi", "latest")
    cands = [
        os.path.join(mpiroot, "bin"),
        os.path.join(mpiroot, "opt", "mpi", "libfabric", "bin"),
        os.path.join(oneapi, "compiler", "latest", "bin"),
    ]
    return [c for c in cands if os.path.isdir(c)]


def find_mpiexec(dirs):
    for d in dirs:
        c = os.path.join(d, "mpiexec.exe")
        if os.path.exists(c):
            return c
    return "mpiexec"


def run(cmd, env, cwd=HERE):
    p = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True,
                       errors="replace", timeout=300)
    return p.returncode, (p.stdout or "") + (p.stderr or "")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    bd = os.path.abspath(sys.argv[1])
    ops = os.path.join(bd, "OpenSees.exe")
    opsmp = os.path.join(bd, "OpenSeesMP.exe")
    for exe in (ops, opsmp):
        if not os.path.exists(exe):
            print("MISSING %s" % exe)
            return 2

    env = dict(os.environ)
    tcllib, d = None, bd
    for _ in range(6):
        cand = os.path.join(d, "dist", "lib", "tcl8.6")
        if os.path.isdir(cand):
            tcllib = cand
            break
        nd = os.path.dirname(d)
        if nd == d:
            break
        d = nd
    if tcllib:
        env["TCL_LIBRARY"] = tcllib
    dirs = oneapi_dirs()
    env["PATH"] = os.pathsep.join([bd] + dirs + [env.get("PATH", "")])
    mpiexec = find_mpiexec(dirs)

    failures = []

    def verdict(name, ok, detail=""):
        print("%-16s %-6s %s" % (name, "PASS" if ok else "FAIL", detail))
        if not ok:
            failures.append(name)

    # --- serial -kn auto reference (P0 deck) ---
    rc, out = run([ops, os.path.join(P0, "serial.tcl")], env)
    m = re.search(r"P0TCL analyze=(-?\d+) w15=(\S+) w11=(\S+) w5=(\S+) R=(\S+)", out)
    ok = rc == 0 and m and int(m.group(1)) == 0
    verdict("serial-auto", ok, m.group(0) if m else out[-300:])
    ref = {k: float(m.group(i)) for i, k in ((2, "w15"), (3, "w11"), (4, "w5"), (5, "R"))} if m else {}

    # --- 2-rank -kn auto: penalty parity through the response ---
    rc, out = run([mpiexec, "-n", "2", opsmp, os.path.join(P0, "mp.tcl")], env)
    m0 = re.search(r"P0TCLMP pid=0 analyze=(-?\d+) w5=(\S+) w11ghost=(\S+) R=(\S+)", out)
    m1 = re.search(r"P0TCLMP pid=1 analyze=(-?\d+) w15=(\S+) w11=(\S+)", out)
    ok = m0 and m1 and int(m0.group(1)) == 0 and int(m1.group(1)) == 0 and ref
    detail = ""
    if ok:
        pairs = [("w5", float(m0.group(2))), ("R", float(m0.group(4))),
                 ("w15", float(m1.group(2))), ("w11", float(m1.group(3)))]
        rel = max(abs(v - ref[k]) / max(abs(ref[k]), 1e-30) for k, v in pairs)
        ok = rel < 1e-12
        detail = "max rel delta vs serial = %.3e" % rel
    else:
        detail = out[-300:]
    verdict("mp-auto", ok, detail)

    # --- -soft serial: must run ---
    rc, out = run([ops, "serial_soft.tcl"], env)
    ok = rc == 0 and "P2SOFTSERIAL analyze=0" in out and "partitioned host" not in out
    verdict("serial-soft", ok, "" if ok else out[-300:])

    # --- -soft single-rank MP: np=1, must run (np-keyed refusal) ---
    rc, out = run([mpiexec, "-n", "1", opsmp, "serial_soft.tcl"], env)
    ok = "P2SOFTSERIAL analyze=0" in out and "partitioned host" not in out
    verdict("mp1-soft", ok, "" if ok else out[-300:])

    # --- -soft 2-rank: must be REFUSED, loudly, and tear the job down ---
    rc, out = run([mpiexec, "-n", "2", opsmp, "mp_soft.tcl"], env)
    named = ("not supported under a partitioned host" in out) and ("-soft" in out)
    ok = named and ("SHOULD-NOT-PRINT" not in out)
    verdict("mp2-soft", ok, "named refusal%s, tail line %s" %
            ("" if named else " MISSING",
             "suppressed" if "SHOULD-NOT-PRINT" not in out else "PRINTED"))

    # --- -soft 2-rank, MORTAR lane: the lane P1 never guarded. Pre-P2 this deck
    #     RUNS (tail line prints); post-P2 it is refused with the same named FATAL.
    rc, out = run([mpiexec, "-n", "2", opsmp, "mp_soft_mortar.tcl"], env)
    named = ("not supported under a partitioned host" in out)
    ok = named and ("SHOULD-NOT-PRINT" not in out)
    verdict("mp2-soft-mortar", ok, "named refusal%s, tail line %s" %
            ("" if named else " MISSING",
             "suppressed" if "SHOULD-NOT-PRINT" not in out else "PRINTED"))

    # --- P2b: definitions survive save -> wipe -> restore ---
    rc, out = run([sys.executable, "db_roundtrip.py", bd], env)
    m = re.search(r"P2DB (\{.*\})", out)
    if not m:
        verdict("db-roundtrip", False, out[-400:])
    else:
        r = json.loads(m.group(1))
        checks = []
        checks.append(("contact analyses converged",
                       all(r[k] == 0 for k in ("ref_ok", "rt_ok", "lv_ok"))))
        # the discriminator must discriminate: WITHOUT contact the top block is
        # floating, so the static solve legitimately fails to converge (nc_ok
        # != 0) and/or lands far from the reference. Either signature means a
        # silently-dropped definitions block cannot masquerade as `ref`.
        checks.append(("discriminator live",
                       r["nc_ok"] != 0 or
                       abs(r["w15_nocontact"] - r["w15_ref"]) > 1e-6 * abs(r["w15_ref"])))
        checks.append(("roundtrip == ref",
                       abs(r["w15_roundtrip"] - r["w15_ref"]) <= 1e-14 * abs(r["w15_ref"])))
        checks.append(("liveverify == ref",
                       abs(r["w15_liveverify"] - r["w15_ref"]) <= 1e-14 * abs(r["w15_ref"])))
        checks.append(("liveverify silent", "differ from the live ones" not in out))
        checks.append(("vanilla db roundtrip",
                       r["van_ref_ok"] == 0 and r["van_rt_ok"] == 0 and
                       abs(r["w5_van_rt"] - r["w5_van_ref"]) <= 1e-14 * abs(r["w5_van_ref"])))
        bad = [n for n, okc in checks if not okc]
        verdict("db-roundtrip", not bad,
                ("w15 ref=%.9e rt=%.9e nc=%.9e" % (r["w15_ref"], r["w15_roundtrip"],
                                                   r["w15_nocontact"]))
                + ((" FAILED: " + ", ".join(bad)) if bad else ""))

    # --- P2b review follow-up: ALL-LANES round trip + field sensitivity ---
    rc, out = run([sys.executable, "db_roundtrip_all_lanes.py", bd], env)
    m = re.search(r"^P2DBALL (\{.*\})", out, re.M)
    if not m:
        verdict("db-all-lanes", False, out[-400:])
    else:
        r = json.loads(m.group(1))
        mark = out.find("P2DBALL-MARK pre-mutation-restore-done")
        warn = out.find("differ from the live ones")
        checks = [
            ("analyses converged", r["ref_ok"] == 0 and r["rt_ok"] == 0),
            ("all four lanes exact", r["tips_ref"] == r["tips_rt"]),
            ("lanes distinct + tie in tension",
             r["tips_ref"]["B"] > 0 > r["tips_ref"]["A"]
             and len(set(r["tips_ref"].values())) == 4),
            ("verify silent until mutated, then warns",
             mark != -1 and warn > mark),
            ("pack deterministic", r["determinism_same"]),
            ("every probed field packed", all(r["sensitivity"].values())),
        ]
        bad = [n for n, okc in checks if not okc]
        verdict("db-all-lanes", not bad,
                ("tips=" + ",".join("%s=%.6e" % kv
                                    for kv in sorted(r["tips_ref"].items())))
                + ((" FAILED: " + ", ".join(bad)) if bad else ""))

    print()
    if failures:
        print("FAILED: " + ", ".join(failures))
        return 1
    print("ALL PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
