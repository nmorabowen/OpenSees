"""ADR-78 P4 validation-battery driver.

    python run.py <build-dir>

<build-dir> holds OpenSees.exe / OpenSeesMP.exe -- normally `dist\\bin` under
the tree you built (named explicitly: a bare `OpenSees` on PATH has resolved
to another session's build on this machine more than once).

Cases (decks in this dir; outputs land in ./_out, wiped per run):

  pound-serial     pound.tcl, OpenSees.exe        explicit pounding reference
  pound-np2        pound.tcl, mpiexec -n 2        2-rank pounding vs serial:
                                                  disp/vel histories + peak base
                                                  (impact) reaction + ghost==native
  pound-np4        pound.tcl, mpiexec -n 4        4-rank (both bodies CUT) vs serial
  energy-serial    closure of the serial energy balance (contact-inactive samples)
  energy-np2       per-rank EnergyBalance files summed == serial channels + closure
  energy-np4       same at np=4
  tie-serial       tie.tcl, OpenSees.exe          mortar mesh-tie + held-load ALM
  tie-np2          tie.tcl, mpiexec -n 2          partitioned tie vs serial + analytic
  scale-np124      pound.tcl NB=8/NSTEP=40000     np=1/2/4 parity + coarse wall-time
                                                  table (DESKTOP SANITY, not a study)
  ctl-soft         pound.tcl + P4_SOFT, np2       P2 refusal must fire (named FATAL)
  ctl-noghost      pound.tcl + P4_NOGHOST, np2    declaration-time refusal must fire
                                                  (named, loud) AND tear the MPI job
                                                  down on its own -- a tree-kill by this
                                                  driver is now a FAILURE
  ctl-ghostmass    pound.tcl + P4_GHOSTMASS, np2  INV-6 mutation: runs, but the
                                                  serial-comparison instrument must
                                                  FLAG it (silent-wrong-answer check)
"""

import glob
import os
import re
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "_out")


def oneapi_dirs():
    # See ../auto_penalty_mpi/run.py for why each entry exists (impi.dll,
    # libfabric under opt\mpi -- NOT mpi\latest\libfabric, compiler runtime).
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


class Runner:
    def __init__(self, bd):
        self.ops = os.path.join(bd, "OpenSees.exe")
        self.opsmp = os.path.join(bd, "OpenSeesMP.exe")
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
        self.env = env
        self.mpiexec = find_mpiexec(dirs)

    def run(self, deck, n=None, extra_env=None, timeout=600):
        """Returns (rc, combined-output, wall-seconds); rc == -9 means the job
        TIMED OUT and its whole process tree was killed.

        Output goes to a FILE, never a pipe: an MPI job whose ranks outlive
        mpiexec (measured here -- the P0 noghost mutation hangs post-#737)
        leaves grandchildren holding an inherited pipe handle, and
        subprocess.run(capture_output=True) then blocks FOREVER, even after
        its own TimeoutExpired killed mpiexec. taskkill /T reaps the tree
        (hydra_pmi_proxy + ranks are descendants of mpiexec)."""
        env = dict(self.env)
        if extra_env:
            env.update(extra_env)
        if n is None:
            cmd = [self.ops, os.path.join(HERE, deck)]
        else:
            cmd = [self.mpiexec, "-n", str(n), self.opsmp, os.path.join(HERE, deck)]
        logf = os.path.join(OUT, "_log_%s.txt" % (extra_env or {}).get("P4_TAG", "run"))
        t0 = time.perf_counter()
        with open(logf, "w") as fh:
            p = subprocess.Popen(cmd, cwd=OUT, env=env, stdout=fh,
                                 stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
            try:
                rc = p.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                subprocess.run(["taskkill", "/F", "/T", "/PID", str(p.pid)],
                               capture_output=True)
                rc = -9
        dt = time.perf_counter() - t0
        with open(logf, "r", errors="replace") as fh:
            return rc, fh.read(), dt


failures = []


def verdict(name, ok, detail=""):
    print("%-14s %-6s %s" % (name, "PASS" if ok else "FAIL", detail))
    if not ok:
        failures.append(name)


def load_hist(tag, kind, name):
    """(time, value) columns of a per-node recorder file."""
    d = np.loadtxt(os.path.join(OUT, "%s_%s_%s.txt" % (kind, tag, name)))
    return np.atleast_2d(d)


def hist_delta(tag_a, tag_b, kind, name):
    a = load_hist(tag_a, kind, name)
    b = load_hist(tag_b, kind, name)
    if a.shape != b.shape:
        return None
    return float(np.max(np.abs(a[:, 1] - b[:, 1])))


def load_energy(tag, ranks):
    """Sum the v1 EnergyBalance columns (KE IE DW ULW RES ERR) over ranks.
    Returns (time, KE, IE, DW, ULW) -- RES/ERR are recomputed globally."""
    tot = None
    t = None
    for r in ranks:
        # under MPI the file handler appends a per-process ".part-<rank>" suffix
        cands = sorted(glob.glob(os.path.join(OUT, "energy_%s_r%d*.txt" % (tag, r))))
        if not cands:
            raise RuntimeError("no energy file for %s rank %d" % (tag, r))
        d = np.atleast_2d(np.loadtxt(cands[0]))
        if tot is None:
            t = d[:, 0]
            tot = d[:, 1:5].copy()
        else:
            if d.shape[0] != tot.shape[0]:
                raise RuntimeError("energy row mismatch for %s rank %d" % (tag, r))
            tot += d[:, 1:5]
    return t, tot  # cols: KE IE DW ULW


def closure(tot):
    """RES(t) = ULW - (KE + IE + DW); E_peak = max(|KE|+|IE|+|DW|+|ULW|)."""
    ke, ie, dw, ulw = tot[:, 0], tot[:, 1], tot[:, 2], tot[:, 3]
    res = ulw - (ke + ie + dw)
    epeak = float(np.max(np.abs(ke) + np.abs(ie) + np.abs(dw) + np.abs(ulw)))
    return res, epeak


def separated_mask(tag, gap):
    """True where the contact is OPEN: gap(t) = w_s1 - w_m1 + GAP0 > 0.
    While the penalty spring is loaded its stored energy is a PHYSICAL RES
    excursion (contact forces are assembled at the SOE, invisible to the
    element/node energy sweep), so closure is gated on separated samples."""
    s1 = load_hist(tag, "disp", "s1")[:, 1]
    m1 = load_hist(tag, "disp", "m1")[:, 1]
    return (s1 - m1 + gap) > 0.0


GAP0 = 5.0e-3
NODES = ("m1", "s1", "t1")


def parse_vals(out, tag):
    vals = {}
    for m in re.finditer(r"P4VAL %s pid=\d+ (\w+)=(\S+)(?: (\w+)=(\S+))?" % re.escape(tag), out):
        vals[m.group(1)] = float(m.group(2))
        if m.group(3):
            vals[m.group(3)] = float(m.group(4))
    return vals


def analyze_ok(out, tag, np_):
    stat = re.findall(r"P4POUND %s pid=(\d+) np=%d analyze=(-?\d+)" % (re.escape(tag), np_), out)
    return len(stat) == np_ and all(int(s[1]) == 0 for s in stat)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    bd = os.path.abspath(sys.argv[1])
    for exe in ("OpenSees.exe", "OpenSeesMP.exe"):
        if not os.path.exists(os.path.join(bd, exe)):
            print("MISSING %s" % os.path.join(bd, exe))
            return 2
    if os.path.isdir(OUT):
        shutil.rmtree(OUT)
    os.makedirs(OUT)
    R = Runner(bd)

    # ---------------- pounding: serial reference -----------------------------
    rc, out, _ = R.run("pound.tcl", extra_env={"P4_TAG": "ser"})
    ok = rc == 0 and analyze_ok(out, "ser", 1)
    verdict("pound-serial", ok, "" if ok else out[-300:])
    if not ok:
        print("\nFAILED: " + ", ".join(failures))
        return 1

    # ---------------- pounding: 2 ranks --------------------------------------
    rc, out2, _ = R.run("pound.tcl", n=2, extra_env={"P4_TAG": "np2"})
    ok = analyze_ok(out2, "np2", 2)
    detail = out2[-300:]
    if ok:
        deltas = {}
        for name in NODES:
            for kind in ("disp", "vel"):
                d = hist_delta("ser", "np2", kind, name)
                deltas["%s.%s" % (kind, name)] = d
        # peak impact force transmitted to the base (sum of the 4 base reactions)
        rs = np.atleast_2d(np.loadtxt(os.path.join(OUT, "reac_ser.txt")))
        r2 = np.atleast_2d(np.loadtxt(os.path.join(OUT, "reac_np2.txt")))
        pk_s = float(np.max(np.abs(rs[:, 1:].sum(axis=1))))
        pk_2 = float(np.max(np.abs(r2[:, 1:].sum(axis=1))))
        # ghost copy of the slave interface node == its native value
        gs = load_hist("np2", "disp", "s1g")[:, 1]
        sn = load_hist("np2", "disp", "s1")[:, 1]
        ghost_d = float(np.max(np.abs(gs - sn)))
        bad = [k for k, v in deltas.items() if v is None or v > 1e-12]
        okc = not bad and abs(pk_2 - pk_s) <= 1e-9 * max(pk_s, 1.0) and ghost_d == 0.0
        detail = ("max|hist delta|=%.3e  peakR ser=%.6e np2=%.6e  ghost-native=%.3e"
                  % (max(v for v in deltas.values() if v is not None), pk_s, pk_2, ghost_d))
        ok = okc
    verdict("pound-np2", ok, detail)

    # ---------------- pounding: 4 ranks (both bodies cut) --------------------
    rc, out4, _ = R.run("pound.tcl", n=4, extra_env={"P4_TAG": "np4"})
    ok = analyze_ok(out4, "np4", 4)
    detail = out4[-300:]
    if ok:
        dmax = 0.0
        for name in NODES:
            for kind in ("disp", "vel"):
                d = hist_delta("ser", "np4", kind, name)
                dmax = max(dmax, 1.0 if d is None else d)
        ok = dmax <= 1e-12
        detail = "max|hist delta| vs serial = %.3e" % dmax
    verdict("pound-np4", ok, detail)

    # ---------------- energy closure + serial/parallel parity ----------------
    CLOSURE_FRAC = 0.01   # gate: |RES| at contact-open samples < 1% of E_peak
    PARITY_FRAC = 1e-9    # gate: per-channel serial-vs-summed-ranks mismatch
    try:
        t_s, tot_s = load_energy("ser", [0])
        res_s, ep_s = closure(tot_s)
        m_s = separated_mask("ser", GAP0)
        cl_s = float(np.abs(res_s[m_s][-1])) if m_s.any() else float("inf")
        ok = cl_s < CLOSURE_FRAC * ep_s and m_s.any()
        verdict("energy-serial", ok,
                "closure@open=%.4e (%.3f%% of E_peak=%.4e), max|RES| in-contact=%.4e"
                % (cl_s, 100.0 * cl_s / ep_s, ep_s, float(np.max(np.abs(res_s)))))
    except Exception as e:  # noqa: BLE001
        verdict("energy-serial", False, repr(e))

    for tag, ranks in (("np2", [0, 1]), ("np4", [0, 1, 2, 3])):
        try:
            t_p, tot_p = load_energy(tag, ranks)
            res_p, ep_p = closure(tot_p)
            m_p = separated_mask(tag, GAP0)
            cl_p = float(np.abs(res_p[m_p][-1])) if m_p.any() else float("inf")
            par = float(np.max(np.abs(tot_p - tot_s))) / max(ep_s, 1e-30)
            ok = cl_p < CLOSURE_FRAC * ep_p and par < PARITY_FRAC
            verdict("energy-%s" % tag, ok,
                    "closure@open=%.4e (%.3f%%), channel parity vs serial=%.3e of E_peak"
                    % (cl_p, 100.0 * cl_p / ep_p, par))
        except Exception as e:  # noqa: BLE001
            verdict("energy-%s" % tag, False, repr(e))

    # ---------------- mortar mesh-tie + ALM ----------------------------------
    UTIP = 2.0e-3  # analytic tip disp (tension PL/EA per unit block)
    rc, outt, _ = R.run("tie.tcl", extra_env={"P4_TAG": "ser"})
    mb = re.search(r"P4TIEB ser pid=0 ok=(-?\d+) okA=(-?\d+) res=(\S+) u5=(\S+) u9g=(\S+)", outt)
    mt = re.search(r"P4TIET ser pid=0 ok=(-?\d+) okA=(-?\d+) u9=(\S+) u15=(\S+)", outt)
    ok = (mb and mt and int(mb.group(1)) == 0 and int(mb.group(2)) == 0)
    tser = {}
    if ok:
        tser = dict(res=float(mb.group(3)), u5=float(mb.group(4)),
                    u9=float(mt.group(3)), u15=float(mt.group(4)))
        ok = (tser["res"] < 1e-9
              and abs(tser["u9"] - tser["u5"]) < 1e-9          # the bond holds
              and abs(tser["u15"] - UTIP) < 1e-6 * UTIP)       # analytic bar solution
    verdict("tie-serial", ok,
            ("res=%.3e u5=%.10e u9=%.10e u15=%.10e (analytic %.1e)"
             % (tser["res"], tser["u5"], tser["u9"], tser["u15"], UTIP)) if tser
            else outt[-300:])

    rc, outm, _ = R.run("tie.tcl", n=2, extra_env={"P4_TAG": "np2"})
    mb = re.search(r"P4TIEB np2 pid=0 ok=(-?\d+) okA=(-?\d+) res=(\S+) u5=(\S+) u9g=(\S+)", outm)
    mt = re.search(r"P4TIET np2 pid=1 ok=(-?\d+) okA=(-?\d+) u9=(\S+) u15=(\S+)", outm)
    ok = (mb and mt and int(mb.group(1)) == 0 and int(mb.group(2)) == 0
          and int(mt.group(1)) == 0 and int(mt.group(2)) == 0 and tser)
    if ok:
        tmp = dict(res=float(mb.group(3)), u5=float(mb.group(4)), u9g=float(mb.group(5)),
                   u9=float(mt.group(3)), u15=float(mt.group(4)))
        rel = max(abs(tmp[k] - tser[k]) / max(abs(tser[k]), 1e-30) for k in ("u5", "u9", "u15"))
        ok = tmp["res"] < 1e-9 and rel < 1e-12 and tmp["u9g"] == tmp["u9"]
        verdict("tie-np2", ok,
                "res=%.3e  max rel delta vs serial=%.3e  ghost-native=%.3e"
                % (tmp["res"], rel, abs(tmp["u9g"] - tmp["u9"])))
    else:
        verdict("tie-np2", False, outm[-300:])

    # ---------------- strong-scaling sanity (DESKTOP, coarse) ----------------
    # NB=8 (16 bricks), 40k steps. This is a correctness-first sanity check on
    # a 4-brick-deep partition; the timing column is COARSE wall time on a
    # desktop (startup + tiny model => not a scaling study, and says so).
    senv = {"P4_NB": "8", "P4_NSTEP": "40000"}
    times, finals = {}, {}
    for n in (1, 2, 4):
        tag = "b8np%d" % n
        rc, o, wall = R.run("pound.tcl", n=n, extra_env=dict(senv, P4_TAG=tag), timeout=900)
        times[n] = wall
        finals[n] = parse_vals(o, tag) if analyze_ok(o, tag, n) else None
    ok = all(finals[n] for n in (1, 2, 4))
    detail = ""
    if ok:
        rel = 0.0
        for n in (2, 4):
            for k in ("t1", "vt1", "s1", "m1"):
                rel = max(rel, abs(finals[n][k] - finals[1][k]) / max(abs(finals[1][k]), 1e-30))
        ok = rel < 1e-12
        detail = ("np1=%.1fs np2=%.1fs np4=%.1fs  max rel delta vs np1=%.3e"
                  % (times[1], times[2], times[4], rel))
    else:
        detail = "a run failed: " + repr({n: bool(finals[n]) for n in times})
    verdict("scale-np124", ok, detail)

    # ---------------- controls (every PASS needs an instrument that can fail)
    # (a) P2 refusal: partitioned -soft pounding deck is REFUSED, loudly.
    rc, o, _ = R.run("pound.tcl", n=2, extra_env={"P4_TAG": "soft", "P4_SOFT": "1"},
                     timeout=180)
    named = "not supported under a partitioned host" in o
    tail = re.search(r"P4POUND soft .* analyze=", o) is not None
    verdict("ctl-soft", named and not tail,
            "named refusal%s, tail %s" % ("" if named else " MISSING",
                                          "suppressed" if not tail else "PRINTED"))

    # (b) fail-loud: ghost declarations stripped. Post-#737 (removal lane) the
    #     missing nodes are refused AT DECLARATION with a named error -- loud,
    #     never a silent pass. The P4 battery ORIGINALLY found that refusal to be
    #     rank-local (a parser `return -1`, no MPI teardown): the other rank
    #     blocked in its next collective and the job had to be reaped by this
    #     driver's tree-kill, regressing P1's 1.1 s FATAL teardown. Fixed by
    #     routing the declaration verbs through ladrunoContactFatal(), so the
    #     gate now requires all three: named refusal, no results, AND the job
    #     ending itself. rc == -9 (tree-killed) is a FAILURE, not an expectation
    #     -- that is what makes this case a standing guard on the teardown.
    #     WALL_MAX is generous against P1's measured 1.1 s: it must separate a
    #     teardown from a hang, not time the teardown.
    WALL_MAX = 20.0
    rc, o, wall = R.run("pound.tcl", n=2, extra_env={"P4_TAG": "nogh", "P4_NOGHOST": "1"},
                        timeout=45)
    named = "does not exist" in o
    tail = re.search(r"P4POUND nogh .* analyze=", o) is not None
    torn_down = rc != -9 and wall < WALL_MAX
    verdict("ctl-noghost", named and not tail and torn_down,
            "declaration refusal%s, tail %s, %s (%.1fs)"
            % ("" if named else " MISSING",
               "suppressed" if not tail else "PRINTED",
               "job HUNG -> tree-killed (TEARDOWN REGRESSED)" if rc == -9
               else ("job exited rc=%d" % rc if wall < WALL_MAX
                     else "job exited rc=%d but took >%.0fs (teardown too slow)"
                          % (rc, WALL_MAX)),
               wall))

    # (c) INV-6 mutation: mass on ghosts runs to completion (silently wrong);
    #     the serial-comparison instrument itself must FLAG the divergence.
    rc, o, _ = R.run("pound.tcl", n=2, extra_env={"P4_TAG": "gmas", "P4_GHOSTMASS": "1"},
                     timeout=180)
    ok = analyze_ok(o, "gmas", 2)
    detail = o[-300:]
    if ok:
        d = hist_delta("ser", "gmas", "disp", "t1")
        ref = float(np.max(np.abs(load_hist("ser", "disp", "t1")[:, 1])))
        ok = d is not None and d > 1e-3 * ref
        detail = ("mutation moved t1 by %.3e (%.2f%% of |t1|max) -> instrument %s"
                  % (d, 100.0 * d / ref, "CATCHES it" if ok else "MISSES IT"))
    verdict("ctl-ghostmass", ok, detail)

    print()
    if failures:
        print("FAILED: " + ", ".join(failures))
        return 1
    print("ALL PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
